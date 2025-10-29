#!/usr/bin/env python
# coding: utf-8

# In[1]:


import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import warnings
warnings.filterwarnings("ignore")


# In[2]:


adata = sc.read_h5ad('../data/mtDNA_DSB_5k_clustered_higher_res.h5ad')


# In[3]:


# --- Find marker genes per cluster ---
sc.tl.rank_genes_groups(adata, groupby="leiden_3", method="t-test")
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)

# get a tidy table of top markers
markers = sc.get.rank_genes_groups_df(adata, group=None)
markers.head()


# In[4]:


marker_genes = pd.DataFrame({
    group: adata.uns['rank_genes_groups']['names'][group][:20]
    for group in adata.uns['rank_genes_groups']['names'].dtype.names
})


# In[5]:


for col in marker_genes.columns: 
    print(col)
    genes = marker_genes[col].tolist()
    print(" ".join(genes))
    print(' ')


# In[6]:


from mllmcelltype import annotate_clusters, setup_logging


# In[7]:


from dotenv import load_dotenv
import os

load_dotenv()  # will read .env into environment
api_key = os.getenv("../OPENAI_API_KEY")


# In[8]:


marker_genes = marker_genes.T


# In[9]:


marker_genes = marker_genes.reset_index().rename(columns={'index':'cluster'})


# In[10]:


marker_genes = marker_genes.set_index("cluster").T.to_dict("list")


# In[11]:


import os

# Annotate clusters with a single model
annotations = annotate_clusters(
    marker_genes=marker_genes,  # DataFrame or dictionary of marker genes
    species='mouse',               # Organism species
    provider='openai',            # LLM provider
    model='gpt-4o-mini',               # Specific model
    tissue='brain'                #Tissue context (optional but recommended)
)

# Print annotations
for cluster, annotation in annotations.items():
    print(f"Cluster {cluster}: {annotation}")


# In[12]:


adata.obs['cell_class'] = adata.obs['leiden_3'].map(annotations)


# In[13]:


sc.pl.umap(adata, color='cell_class', frameon=False)


# In[14]:


spatial = np.array(adata.obs[['x_centroid','y_centroid']])
adata.obsm['spatial'] = spatial


# In[36]:


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
        base = sc.pl.palettes.default_64 if hasattr(sc.pl.palettes, "default_64") else sc.pl.palettes.default_102
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

    # RETURN instead of showing
    return fig


# In[16]:


list(adata.obs['run'].unique())


# In[17]:


adata.obs["sample_id"] = adata.obs["run"].str.extract(r"__([A-Za-z0-9]+)__\d{8}__")


# In[38]:


import matplotlib.pyplot as plt
import os

out_dir = "../results/figures"
os.makedirs(out_dir, exist_ok=True)

for domain in ["cell_class"]:
    if domain in adata.obs.columns:
        print(f"Plotting and saving: {domain}")
        fig = plot_spatial_compact(
            adata,
            color=domain,
            groupby="sample_id",
            spot_size=20,
            cols=6,
            height=8,
            legend_col_width=1.0,
        )
        save_base = os.path.join(out_dir, f"spatial_{domain}")
        for ext in ["png"]:
            fig.savefig(f"{save_base}.{ext}", bbox_inches="tight", dpi=300)
        plt.close(fig)
    else:
        print(f"Skipping {domain}: not in adata.obs")


# In[39]:


for domain in ['cell_class']:
    if domain in adata.obs.columns: 
        plot_spatial_compact(adata, color=domain, groupby="sample_id",spot_size=20,
            cols=6,
            height=8,
            legend_col_width=1.0,)
    else:
        continue


# In[19]:


# --- Find marker genes per cluster ---
sc.tl.rank_genes_groups(adata, groupby="cell_class", method="t-test")
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)

# get a tidy table of top markers
markers = sc.get.rank_genes_groups_df(adata, group=None)
markers.head()


# In[20]:


marker_genes = pd.DataFrame({
    group: adata.uns['rank_genes_groups']['names'][group][:4]
    for group in adata.uns['rank_genes_groups']['names'].dtype.names
})


# In[21]:


cellclass_markers_top4 = {}
for col in marker_genes.columns: 
    print(col)
    genes = marker_genes[col].tolist()
    print(" ".join(genes))
    print(' ')
    cellclass_markers_top4[col] = genes


# In[25]:


adata.X.max()


# In[30]:


sc.pl.dotplot(
    adata,
    var_names=cellclass_markers_top4,
    groupby="cell_class",
    color_map="coolwarm",
    vmax = 0.8,
    vmin = 0, 
    dendrogram=True
)


# In[29]:


adata.write('../data/mtDNA_DSB_5k_clustered_higher_res_LLM_anno.h5ad')


# In[ ]:




