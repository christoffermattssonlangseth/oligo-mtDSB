#!/usr/bin/env python
# coding: utf-8

# # spatial neighbourhoods

# In[1]:


import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.neighbors import radius_neighbors_graph
import scanpy as sc

def spatial_neighbourhoods_fast(
    ad,
    cluster_label="leiden_0.5",
    max_distance=200.0,
    x_col="x_centroid",
    y_col="y_centroid",
):
    """
    Returns an AnnData where X is (cells x clusters) with counts of neighbors
    within `max_distance` per cluster. Uses sparse ops (fast + memory-light).
    """
    # --- coordinates ---
    coords = ad.obs[[x_col, y_col]].to_numpy()

    # Sparse adjacency (N x N) with 1 if within radius, excluding self
    A = radius_neighbors_graph(
        coords, radius=max_distance, mode="connectivity",
        include_self=False, n_jobs=-1
    ).tocsr()

    # --- cluster one-hot (N x K) ---
    cats = ad.obs[cluster_label].astype("category")
    cat_codes = cats.cat.codes.to_numpy()
    N = ad.n_obs
    K = len(cats.cat.categories)
    rows = np.arange(N, dtype=int)
    C = sp.csr_matrix((np.ones(N, dtype=np.int32), (rows, cat_codes)), shape=(N, K))

    # --- neighbor counts per cluster: F = A @ C  (N x K) ---
    F = A @ C  # integer counts

    # Build output AnnData
    out = sc.AnnData(
        X=F,
        obs=ad.obs.copy(),
        var=pd.DataFrame(index=[str(c) for c in cats.cat.categories])
    )
    # Keep the adjacency if you want it later (in obsp, not layers)
    out.obsp["neighbors_within_radius"] = A
    out.uns["spatial_neighbourhoods"] = {
        "cluster_label": cluster_label,
        "max_distance": float(max_distance),
        "x_col": x_col,
        "y_col": y_col,
    }
    return out



# In[2]:


adata = sc.read_h5ad('../data/mtDNA_DSB_5k_clustered_manual_annotation.h5ad')


# # run spatial neigbourshoods

# In[3]:


sc.pl.umap(adata, color = 'cell_class')


# In[4]:


list(adata.obs.sample_id.unique())


# In[5]:


get_ipython().run_cell_magic('time', '', 'sample_ids = [\'RB4282\',\n \'RB4350\',\n \'RB4401\',\n \'RB4403\',\n \'RB4405\',\n \'RB4498\',\n \'RB4620\',\n \'RB4627\',\n \'RB4630\',\n \'RB4653\',\n \'RB4658\',\n \'RB4676\']\n# or: adata.obs["sample_id"].unique().tolist()\nad_sn_list = [\n    spatial_neighbourhoods_fast(\n        adata[adata.obs["sample_id"] == sid],\n        cluster_label="cell_class",\n        max_distance=100,\n        x_col="x_centroid",\n        y_col="y_centroid",\n    )\n    for sid in sample_ids\n]\n')


# In[6]:


ad_sn_list


# In[7]:


get_ipython().run_cell_magic('time', '', "ad_sn = sc.concat(ad_sn_list, join = 'outer', fill_value = 0)\n")


# In[8]:


get_ipython().run_cell_magic('time', '', "spatial = np.array(ad_sn.obs[['x_centroid','y_centroid']])\nad_sn.obsm['spatial'] = spatial\n")


# In[14]:


ad_sn.X = ad_sn.X.astype('float32')


# In[15]:


get_ipython().run_cell_magic('time', '', 'sc.pp.pca(ad_sn, n_comps=20, svd_solver="randomized")\n')


# In[18]:


get_ipython().run_cell_magic('time', '', 'sc.pp.neighbors(ad_sn, n_neighbors=10, use_rep="X_pca", n_pcs=20)\n')


# In[23]:


sc.tl.umap(
    ad_sn,
    init_pos="random",
    min_dist=0.2,
    negative_sample_rate=2,    # default 5; lower = faster
    maxiter=200                # additional cap on optimization steps
)


# In[ ]:


get_ipython().run_cell_magic('time', '', "for i in [0.01,  0.04, 0.06,0.1,0.15,0.2,0.3,0.5]: \n    key = 'local_neighborhood_'+str(i)\n    if key in ad_sn.obs.columns: \n        sc.pl.umap(ad_sn,color=[key],  s = 10)#,save='UMAP_10X_colors.svg')\n    else: \n        sc.tl.leiden(ad_sn,resolution=i, key_added = key)\n        sc.pl.umap(ad_sn,color=[key],  s = 10)#,save='UMAP_10X_colors.svg')\n")


# In[15]:


import seaborn as sns

key = "local_neighborhood_0.15"
cats = ad_sn.obs[key].astype("category").cat.categories
palette = sns.color_palette("tab20", n_colors=len(cats)).as_hex()

ad_sn.uns[f"{key}_colors"] = palette
sc.pl.umap(ad_sn, color=key, frameon=False, size=6)


# In[69]:


import matplotlib.pyplot as plt
import scanpy as sc

target_cluster = "27"   # <-- put your cluster name here
samples_ = ad_sn[ad_sn.obs['local_neighborhood_0.1'] == target_cluster].obs.sample_id.value_counts().head(5)

for run in list(samples_.index):
    print(run)
    ad_int = ad_sn[ad_sn.obs['sample_id'] == run].copy()
    # build a palette: target cluster → red, others → grey
    clusters = ad_int.obs['local_neighborhood_0.1'].unique().tolist()
    palette = {
        c: ("red" if c == target_cluster else "lightgrey")
        for c in clusters
    }

    with plt.rc_context({'figure.figsize': (20, 10)}):
        sc.pl.spatial(
            ad_int,
            color="local_neighborhood_0.1",
            spot_size=15,
            palette=palette
        )
    plt.show()


# In[25]:


import matplotlib.pyplot as plt


# In[18]:


# 5) Plot a spatial map for your favorite resolution:
for sample_id in ad_sn.obs.sample_id.unique():
    print(sample_id)
    with plt.rc_context({"figure.figsize": (10, 10)}):
        ad_sn_ = ad_sn[ad_sn.obs['sample_id'] == sample_id]
        sc.pl.spatial(
            ad_sn_,
            color='local_neighborhood_0.1',
            spot_size=15,
            title="Neighbours-based domains (resolution=0.1)"
        )


# In[45]:


import scanpy as sc
import pandas as pd
import scipy
import matplotlib.pyplot as plt

def post_merge(df, labels, post_merge_cutoff, linkage_method='single', 
               linkage_metric='correlation', fcluster_criterion='distance', name='', save=False):
    Z = scipy.cluster.hierarchy.linkage(df.T, method=linkage_method, metric=linkage_metric)
    merged_labels_short = scipy.cluster.hierarchy.fcluster(Z, post_merge_cutoff, criterion=fcluster_criterion)

    # map cluster name to merged label
    label_conversion = dict(zip(df.columns, merged_labels_short))
    new_labels = [label_conversion[i] for i in labels]

    # Plot
    fig, ax = plt.subplots(figsize=(20, 10))
    scipy.cluster.hierarchy.dendrogram(Z, labels=df.columns, color_threshold=post_merge_cutoff, ax=ax)
    ax.hlines(post_merge_cutoff, 0, ax.get_xlim()[1])
    ax.set_title('Merged clusters')
    ax.set_ylabel(linkage_metric, fontsize=20)
    ax.set_xlabel('Pre-merge cluster labels', fontsize=20)
    ax.tick_params(axis="x", labelrotation=90, labelsize=10)  # <-- rotate tick labels
    plt.show()
    if save:
        fig.savefig(f"{name}.svg", dpi=500)

    return new_labels


# In[47]:


# Compute cluster-level mean expression
cluster_key = 'local_neighborhood_0.1'
expr = pd.DataFrame(ad_sn.X, index=ad_sn.obs_names, columns=ad_sn.var_names)
cluster_means = expr.groupby(ad_sn.obs[cluster_key]).mean().T  # genes x clusters

# Run post-merge
merged_labels = post_merge(cluster_means, labels=ad_sn.obs[cluster_key].tolist(),
                           post_merge_cutoff=0.01,
                           linkage_metric='correlation',
                           linkage_method='average',
                           name='SupFig3Dend',
                           save=False)

# Add new merged label to AnnData object
ad_sn.obs['merged_cluster'] = merged_labels
ad_sn.obs['merged_cluster'] = ad_sn.obs['merged_cluster'].astype(str)


# In[48]:


import warnings
warnings.filterwarnings("ignore")


# In[49]:


sc.tl.rank_genes_groups(ad_sn, groupby='local_neighborhood_0.1', method='wilcoxon')
# See top 5 marker genes per cluster
sc.pl.rank_genes_groups(ad_sn, n_genes=25, sharey=False)

marker_genes = pd.DataFrame({
    group: ad_sn.uns['rank_genes_groups']['names'][group][:10]
    for group in ad_sn.uns['rank_genes_groups']['names'].dtype.names
})
marker_genes.head()
#marker_genes.to_csv('../data/broad_markers_leiden0-5.csv')
for col in marker_genes.columns: 
    print(col)
    genes = marker_genes[col].tolist()
    print(" ".join(genes))
    print(' ')


# In[50]:


dict(zip(ad_sn.obs["local_neighborhood_0.1"].cat.categories, ad_sn.obs["local_neighborhood_0.1"].cat.categories))


# In[70]:


anno = {'0': 'Isocortex',
 '1': 'Olfactory areas',
 '2': 'Thalamus',
 '3': 'Hypothalamus',
 '4': 'Palladium + corpus callosum',
 '5': 'Striatum',
 '6': 'Striatum like amygdalar nuclei',
 '7': 'Medial habenula',
 '8': 'Meninges',
 '9': 'Lateral ventricle I',
 '10': 'Third ventricle',
 '11': 'Dentate gyrus',
 '12': 'Hindbrain', # not too sure about this one
 '13': 'Uknown',
 '14': 'Cortical subplate',
 '15': 'Lateral venricle II',
 '16': 'Perireunensis nucleus',
 '17': 'Isocortex',
 '18': 'Meninges',
 '19': 'Unknown',
 '20': 'Isocortex',
 '21': 'Unknown',
 '22': 'Unknown',
 '23': 'Unknown',
 '24': 'Unknown',
 '25': 'Isocortex',
 '26': 'Unknown',
 '27': 'Unknown'}


# In[71]:


ad_sn.obs['compartment'] = ad_sn.obs["local_neighborhood_0.1"].map(anno)


# In[79]:


import seaborn as sns

key = "compartment"
cats = ad_sn.obs[key].astype("category").cat.categories
palette = sns.color_palette("Paired", n_colors=len(cats)).as_hex()

ad_sn.uns[f"{key}_colors"] = palette
sc.pl.umap(ad_sn, color=key, frameon=False, size=6)


# In[80]:


for run in ad_sn.obs['sample_id'].unique():
    print(run)
    ad_int = ad_sn[ad_sn.obs['sample_id'] == run]

    with plt.rc_context({'figure.figsize': (20, 10)}):
        sc.pl.spatial(ad_int, spot_size=15, color = 'compartment')
    plt.show()



# In[81]:


adata.obs['compartment'] = adata.obs.index.map(dict(zip(ad_sn.obs.index,ad_sn.obs.compartment)))


# In[83]:


adata.obs['compartment']


# In[86]:


adata.write('../data/mtDNA_DSB_5k_clustered_manual_annotation_compartment.h5ad')


# In[2]:


adata = sc.read('../data/mtDNA_DSB_5k_clustered_manual_annotation_compartment.h5ad')


# In[4]:


import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec

def plot_spatial_compact(
    ad,
    color="leiden_2",
    groupby="sample_id",
    spot_size=18,
    cols=3,
    height=8,
    legend_col_width=1.2,
    palette=None,          # <- NEW: dict {category: "#hex"} or list in desired order
):
    # --- categories from the FULL object (freeze order) ---
    cats = ad.obs[color].astype("category").cat.categories

    # --- build a shared color list in that order ---
    if isinstance(palette, dict):
        cols_list = [palette[c] for c in cats]
    elif isinstance(palette, (list, tuple)):
        if len(palette) < len(cats):
            raise ValueError("Palette list shorter than number of categories.")
        cols_list = list(palette)[:len(cats)]
    elif f"{color}_colors" in ad.uns:
        # reuse existing palette if present on the full object
        cols_list = list(ad.uns[f"{color}_colors"])
        if len(cols_list) != len(cats):
            raise ValueError(f"{color}_colors length does not match categories.")
    else:
        base = sc.plotting.palettes.default_64
        reps = int(np.ceil(len(cats) / len(base)))
        cols_list = (base * reps)[:len(cats)]

    # store on the full object too (handy if you plot elsewhere later)
    ad.uns[f"{color}_colors"] = cols_list

    # --- layout ---
    sids = list(ad.obs[groupby].unique())
    n = len(sids)
    rows = int(np.ceil(n / cols))

    panel_w = height * cols * 0.6 / rows
    fig_w = panel_w + legend_col_width
    fig = plt.figure(figsize=(fig_w, height), constrained_layout=False)

    gs = GridSpec(
        rows, cols + 1, figure=fig,
        width_ratios=[1]*cols + [legend_col_width / (fig_w - legend_col_width)],
        wspace=0.02, hspace=0.02
    )

    # --- panels ---
    for i, sid in enumerate(sids):
        r, c = divmod(i, cols)
        ax = fig.add_subplot(gs[r, c])

        ad_sub = ad[ad.obs[groupby] == sid].copy()

        # CRUCIAL: force identical categories & palette on each subset
        ad_sub.obs[color] = ad_sub.obs[color].astype("category")
        ad_sub.obs[color] = ad_sub.obs[color].cat.set_categories(cats)
        ad_sub.uns[f"{color}_colors"] = cols_list

        sc.pl.spatial(
            ad_sub, color=color, spot_size=spot_size,
            show=False, ax=ax, legend_loc=None, frameon=False, title=str(sid)
        )
        ax.set_xlabel(""); ax.set_ylabel("")
        ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)

    # blank unused
    for j in range(n, rows*cols):
        r, c = divmod(j, cols)
        fig.add_subplot(gs[r, c]).axis("off")

    # --- legend ---
    ax_leg = fig.add_subplot(gs[:, -1])
    ax_leg.axis("off")
    handles = [
        Line2D([0], [0], marker='o', color='w',
               markerfacecolor=cols_list[k], markersize=7, label=str(cat))
        for k, cat in enumerate(cats)
    ]
    ax_leg.legend(handles=handles, title=color, frameon=False, loc="center left")

    fig.subplots_adjust(left=0.01, right=0.98, top=0.98, bottom=0.02, wspace=0.02, hspace=0.02)
    plt.show()


# In[21]:


adata.obs["compartment"] = (
    adata.obs["compartment"]
    .replace({"Olfactory areas": "Olfactory areas + Hippocampal formation"})
    .astype("category")
)


# In[23]:


adata.obs["compartment"]


# In[27]:


plot_spatial_compact(adata, color="cell_class", groupby="sample_id",spot_size=20,
    cols=6,
    height=8,
    legend_col_width=1.0,)


# In[25]:


compartment_palette = {
    "Cortical subplate": "#1f77b4",
    "Dentate gyrus": "#ff7f0e",
    "Hindbrain": "#2ca02c",
    "Hypothalamus": "#d62728",
    "Isocortex": "#b963fc",
    "Lateral venricle II": "#8c564b",
    "Lateral ventricle I": "#e377c2",
    "Meninges": "#29c0d2",
    "Medial habenula": "#bcbd22",
    "Olfactory areas + Hippocampal formation": "#b4c9e7",
    "Palladium + corpus callosum": "#aec7e8",
    "Perireunensis nucleus": "#ffbb78",
    "Striatum": "#98df8a",
    "Striatum like amygdalar nuclei": "#ff9896",
    "Thalamus": "#c5b0d5",
    "Third ventricle": "#c49c94",
    "Unknown": "#c7c7c7",
    "Uknown": "#c7c7c7"
}
plot_spatial_compact(adata, color="compartment", groupby="sample_id", palette=compartment_palette,spot_size=20,
    cols=6,
    height=8,
    legend_col_width=1.0,)


# In[31]:


adata[adata.obs.compartment == 'Hindbrain'].obs.cell_class.value_counts()


# In[93]:


del ad_sn.uns['compartment_colors']


# In[94]:


ad_sn.write('../data/mtDNA_DSB_5k_neigh.h5ad')


# In[13]:


ad_sn = sc.read('../data/mtDNA_DSB_5k_neigh.h5ad')


# In[95]:


sc.tl.rank_genes_groups(adata, groupby='compartment', method='t-test')
# See top 5 marker genes per cluster
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)

marker_genes = pd.DataFrame({
    group: adata.uns['rank_genes_groups']['names'][group][:10]
    for group in adata.uns['rank_genes_groups']['names'].dtype.names
})
marker_genes.head()
#marker_genes.to_csv('../data/broad_markers_leiden0-5.csv')
for col in marker_genes.columns: 
    print(col)
    genes = marker_genes[col].tolist()
    print(" ".join(genes))
    print(' ')


# In[96]:


mtDSB_genes = [
    "Hspa5", "Hspa9", "Hsph1",
    "Atf5", "Trib3", "Zbtb16", "Ddit3",
    "Cdkn1a", "Bcl2l1",
    "Sgk1", "Nmu", "Plin4", "Aldoa",
    "Serpina3n", "Mt2", "Gstp1",'Ldha'
]


# In[103]:


sc.pl.dotplot(
        adata,
        var_names=mtDSB_genes,
        groupby="compartment",
        standard_scale="var",
        #dot_max=0.5,
        #dot_min=0.05,
        color_map="Reds",

        dendrogram=True,
        figsize=(15, 4)
    )


# In[102]:


sc.pl.dotplot(
        adata[adata.obs.condition == 'mtDSB'],
        var_names=mtDSB_genes,
        groupby="compartment",
        standard_scale="var",
        #dot_max=0.5,
        #dot_min=0.05,
        color_map="Reds",

        dendrogram=True,
        figsize=(15, 4)
    )


# In[104]:





# In[105]:


adata_mTDSB = adata[adata.obs.condition == 'mtDSB']
sc.tl.rank_genes_groups(adata_mTDSB, groupby='compartment', method='t-test')
# See top 5 marker genes per cluster
sc.pl.rank_genes_groups(adata_mTDSB, n_genes=25, sharey=False)

marker_genes = pd.DataFrame({
    group: adata_mTDSB.uns['rank_genes_groups']['names'][group][:10]
    for group in adata_mTDSB.uns['rank_genes_groups']['names'].dtype.names
})
marker_genes.head()
#marker_genes.to_csv('../data/broad_markers_leiden0-5.csv')
for col in marker_genes.columns: 
    print(col)
    genes = marker_genes[col].tolist()
    print(" ".join(genes))
    print(' ')


# In[ ]:




