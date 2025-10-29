#!/usr/bin/env python
# coding: utf-8

# # spatial neighbourhoods

# In[2]:


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



# In[3]:


adata = sc.read_h5ad('../data/mtDNA_DSB_5k_clustered_manual_annotation.h5ad')


# # run spatial neigbourshoods

# In[4]:


sc.pl.umap(adata, color = 'cell_class')


# In[5]:


adata.obs.sample_id.unique()


# In[54]:


get_ipython().run_cell_magic('time', '', 'sample_ids = ["RB4405"]  # or: adata.obs["sample_id"].unique().tolist()\nad_sn_list = [\n    spatial_neighbourhoods_fast(\n        adata[adata.obs["sample_id"] == sid],\n        cluster_label="cell_class",\n        max_distance=100,\n        x_col="x_centroid",\n        y_col="y_centroid",\n    )\n    for sid in sample_ids\n]\n')


# In[55]:


ad_sn_list


# In[56]:


get_ipython().run_cell_magic('time', '', "ad_sn = sc.concat(ad_sn_list, join = 'outer', fill_value = 0)\n")


# In[57]:


get_ipython().run_cell_magic('time', '', "spatial = np.array(ad_sn.obs[['x_centroid','y_centroid']])\nad_sn.obsm['spatial'] = spatial\n")


# In[58]:


get_ipython().run_cell_magic('time', '', 'sc.pp.neighbors(ad_sn, n_neighbors=10)\nsc.tl.umap(ad_sn,min_dist=0.2)\n')


# In[88]:


get_ipython().run_cell_magic('time', '', "for i in [0.01,  0.04, 0.06,0.1,0.15,0.2,0.3,0.5]: \n    key = 'local_neighborhood_'+str(i)\n    if key in ad_sn.obs.columns: \n        sc.pl.umap(ad_sn,color=[key],  s = 10)#,save='UMAP_10X_colors.svg')\n    else: \n        sc.tl.leiden(ad_sn,resolution=i, key_added = key)\n        sc.pl.umap(ad_sn,color=[key],  s = 10)#,save='UMAP_10X_colors.svg')\n")


# In[60]:


import matplotlib.pyplot as plt


# In[91]:


# 5) Plot a spatial map for your favorite resolution:
with plt.rc_context({"figure.figsize": (10, 10)}):
    sc.pl.spatial(
        ad_sn,
        color='local_neighborhood_0.1',
        spot_size=15,
        title="Neighbours-based domains (resolution=0.1)"
    )


# In[92]:


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


# In[95]:


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


# In[96]:


for run in ad_sn.obs['sample_id'].unique():
    print(run)
    ad_int = ad_sn[ad_sn.obs['sample_id'] == run]

    with plt.rc_context({'figure.figsize': (20, 10)}):
        sc.pl.spatial(ad_int, spot_size=15, color = 'local_neighborhood_0.1')
    plt.show()



# In[99]:


sc.tl.rank_genes_groups(ad_sn, groupby='merged_cluster', method='wilcoxon')
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


# In[101]:


dict(zip(ad_sn.obs["local_neighborhood_0.1"].cat.categories, ad_sn.obs["local_neighborhood_0.1"].cat.categories))


# In[106]:


anno = {'0': 'Isocortex',
 '1': 'Hypothalamus',
 '2': 'Olfactory areas',
 '3': 'Thalamus',
 '4': 'Striatum',
 '5': 'Striatum-like amygdalar nculei',
 '6': 'Palladium',
 '7': 'Corpus callosum',
 '8': 'Fornix system',
 '9': 'Meninges',
 '10': 'Medial habenula',
 '11': 'Piriform area',
 '12': 'Isocortex',
 '13': 'Isocortex',
 '14': 'Fornix system'}


# In[107]:


ad_sn.obs['compartment'] = ad_sn.obs["local_neighborhood_0.1"].map(anno)


# In[109]:


for run in ad_sn.obs['sample_id'].unique():
    print(run)
    ad_int = ad_sn[ad_sn.obs['sample_id'] == run]

    with plt.rc_context({'figure.figsize': (20, 10)}):
        sc.pl.spatial(ad_int, spot_size=15, color = 'compartment')
    plt.show()



# In[115]:


ad_sub = adata[adata.obs["sample_id"] == 'RB4405']


# In[116]:


ad_sub.obs['compartment'] = ad_sub.obs.index.map(dict(zip(ad_sn.obs.index,ad_sn.obs.compartment)))


# In[117]:


ad_sub.obs['compartment']


# In[118]:


sc.tl.rank_genes_groups(ad_sub, groupby='compartment', method='wilcoxon')
# See top 5 marker genes per cluster
sc.pl.rank_genes_groups(ad_sub, n_genes=25, sharey=False)

marker_genes = pd.DataFrame({
    group: ad_sub.uns['rank_genes_groups']['names'][group][:10]
    for group in ad_sub.uns['rank_genes_groups']['names'].dtype.names
})
marker_genes.head()
#marker_genes.to_csv('../data/broad_markers_leiden0-5.csv')
for col in marker_genes.columns: 
    print(col)
    genes = marker_genes[col].tolist()
    print(" ".join(genes))
    print(' ')


# In[122]:


mtDSB_genes = [
    "Hspa5", "Hspa9", "Hsph1",
    "Atf5", "Trib3", "Zbtb16", "Ddit3",
    "Cdkn1a", "Bcl2l1",
    "Sgk1", "Nmu", "Plin4", "Aldoa",
    "Serpina3n", "Mt2", "Gstp1",'Ldha'
]


# In[127]:


sc.pl.dotplot(
        ad_sub,
        var_names=mtDSB_genes,
        groupby="compartment",
        standard_scale="var",
        #dot_max=0.5,
        #dot_min=0.05,
        color_map="Reds",

        dendrogram=True,
        figsize=(7, 3)
    )


# In[ ]:




