#!/usr/bin/env python
# coding: utf-8

# In[13]:


import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import warnings
warnings.filterwarnings("ignore")


# In[2]:


adata = sc.read_h5ad('../data/mtDNA_DSB_5k_clustered_LLM_anno.h5ad')


# In[3]:


adata_OL = adata[adata.obs['cell_class'].str.contains('ligo')]


# In[6]:


# dimensionality reduction
sc.pp.pca(adata_OL, n_comps=30)
sc.pp.neighbors(adata_OL, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata_OL)


# In[51]:


resolutions = [0.3, 0.5,0.8, 1]

for resolution in resolutions:
    key = f'OL_leiden_{resolution}'

    if key in adata_OL.obs.columns:
        print(f"Skipping {resolution}: {key} already exists.")
    else:
        print(f"Clustering at resolution {resolution}...")
        sc.tl.leiden(adata_OL, resolution=resolution, key_added=key)
        print("Done.")

    # plot UMAP
    sc.pl.umap(adata_OL, color=key, legend_loc='on data', frameon=False)


# In[52]:


# --- Find marker genes per cluster ---
sc.tl.rank_genes_groups(adata_OL, groupby="OL_leiden_1", method="t-test")
sc.pl.rank_genes_groups(adata_OL, n_genes=25, sharey=False)

# get a tidy table of top markers
markers = sc.get.rank_genes_groups_df(adata_OL, group=None)
markers.head()


# In[53]:


marker_genes = pd.DataFrame({
    group: adata_OL.uns['rank_genes_groups']['names'][group][:20]
    for group in adata_OL.uns['rank_genes_groups']['names'].dtype.names
})


# In[54]:


for col in marker_genes.columns: 
    print(col)
    genes = marker_genes[col].tolist()
    print(" ".join(genes))
    print(' ')


# In[62]:


sc.pl.dotplot(
       adata_OL,
       var_names=['Mbp','B2m', 'Epas1', 'H2-D1', 'H2-K1', 'Cd74','Serpina3n','Gfap','C4b','Slc16a3','Slc16a1','Ldha','Ldhb'],
       groupby="OL_leiden_1",
       standard_scale="var",
       #dot_max=0.5,
       #dot_min=0.05,
       color_map="Reds",
       dendrogram=False,
       figsize=(7, 5),

   )


# In[59]:


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# pick your cluster key (change if needed)
cluster_key = "OL_leiden_1"  

# count cells per replicate × cluster
ct = pd.crosstab(adata_OL.obs[cluster_key], adata_OL.obs["run"])

# plot counts
plt.figure(figsize=(8,5))
ct.T.plot(kind="bar", stacked=True, figsize=(8,5))
plt.ylabel("Number of cells")
plt.title("Oligodendrocyte subclusters per replicate")
plt.legend(title="Subcluster", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.show()

# plot fractions
fractions = ct.div(ct.sum(axis=0), axis=1)

plt.figure(figsize=(8,5))
fractions.T.plot(kind="bar", stacked=True, figsize=(8,5))
plt.ylabel("Fraction of cells")
plt.title("Oligodendrocyte subclusters per replicate (fractions)")
plt.legend(title="Subcluster", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.show()


# In[61]:


sc.pl.dotplot(
       adata,
       var_names=['B2m', 'Epas1', 'H2-D1', 'H2-K1', 'Cd74','Serpina3n','Gfap','C4b','Slc16a3','Slc16a1','Ldha','Ldhb'],
       groupby="run",
       standard_scale="var",
       #dot_max=0.5,
       #dot_min=0.05,
       color_map="Reds",
       dendrogram=False,
       figsize=(5, 5),

   )


# In[64]:


for run in adata_OL.obs['run'].unique():
    print(run)
    ad_int = adata_OL[adata_OL.obs['run'] == run]
    with plt.rc_context({'figure.figsize': (20, 10)}):
        sc.pl.spatial(ad_int, spot_size=25, color = 'OL_leiden_1')
    plt.show()



# In[67]:


to_remove = [
    "Acadm", "Acadvl", "Atp5f1a", "Cox4i1",
    "Fis1", "Idh3a", "Mdh2", "Prdx3",
    "Suclg1", "Uqcrc1"
]

mito_genes_mouse = [g for g in mito_genes_mouse if g not in to_remove]

print(mito_genes_mouse)


# # mitochondrial expression

# In[70]:


sc.pl.dotplot(
       adata_OL,
       var_names=mito_genes_mouse,
       groupby="OL_leiden_1",
       standard_scale="var",
       #dot_max=0.5,
       #dot_min=0.05,
       color_map="Reds",
       dendrogram=False,
       figsize=(7, 5),

   )


# In[71]:


sc.pl.dotplot(
       adata,
       var_names=mito_genes_mouse,
       groupby="run",
       standard_scale="var",
       #dot_max=0.5,
       #dot_min=0.05,
       color_map="Reds",
       dendrogram=False,
       figsize=(7, 5),

   )


# In[ ]:




