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


adata = sc.read_h5ad('../data/mtDNA_DSB_5k_clustered_LLM_anno.h5ad')


# In[3]:


adata_MG = adata[adata.obs['cell_class'].str.contains('Micro')]


# In[4]:


# dimensionality reduction
sc.pp.pca(adata_MG, n_comps=30)
sc.pp.neighbors(adata_MG, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata_MG)


# In[5]:


resolutions = [0.3, 0.5,0.8, 1]

for resolution in resolutions:
    key = f'MG_leiden_{resolution}'

    if key in adata_MG.obs.columns:
        print(f"Skipping {resolution}: {key} already exists.")
    else:
        print(f"Clustering at resolution {resolution}...")
        sc.tl.leiden(adata_MG, resolution=resolution, key_added=key)
        print("Done.")

    # plot UMAP
    sc.pl.umap(adata_MG, color=key, legend_loc='on data', frameon=False)


# In[8]:


# --- Find marker genes per cluster ---
sc.tl.rank_genes_groups(adata_MG, groupby="MG_leiden_0.8", method="t-test")
sc.pl.rank_genes_groups(adata_MG, n_genes=25, sharey=False)

# get a tidy table of top markers
markers = sc.get.rank_genes_groups_df(adata_MG, group=None)
markers.head()


# In[9]:


marker_genes = pd.DataFrame({
    group: adata_MG.uns['rank_genes_groups']['names'][group][:20]
    for group in adata_MG.uns['rank_genes_groups']['names'].dtype.names
})


# In[10]:


for col in marker_genes.columns: 
    print(col)
    genes = marker_genes[col].tolist()
    print(" ".join(genes))
    print(' ')


# In[13]:


sc.pl.dotplot(
       adata_MG,
       var_names=['B2m', 'Epas1', 'H2-D1', 'H2-K1', 'Cd74','Serpina3n','Gfap','C4b','Slc16a3','Slc16a1','Ldha','Ldhb'],
       groupby='MG_leiden_0.8',
       standard_scale="var",
       #dot_max=0.5,
       #dot_min=0.05,
       color_map="Reds",
       dendrogram=False,
       figsize=(7, 3),

   )


# In[15]:


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# pick your cluster key (change if needed)
cluster_key = "MG_leiden_1"  

# count cells per replicate × cluster
ct = pd.crosstab(adata_MG.obs[cluster_key], adata_MG.obs["run"])

# plot counts
plt.figure(figsize=(8,5))
ct.T.plot(kind="bar", stacked=True, figsize=(8,5))
plt.ylabel("Number of cells")
plt.title("Microglia subclusters per replicate")
plt.legend(title="Subcluster", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.show()

# plot fractions
fractions = ct.div(ct.sum(axis=0), axis=1)

plt.figure(figsize=(8,5))
fractions.T.plot(kind="bar", stacked=True, figsize=(8,5))
plt.ylabel("Fraction of cells")
plt.title("Microglia subclusters per replicate (fractions)")
plt.legend(title="Subcluster", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.show()


# In[17]:


for run in adata_MG.obs['run'].unique():
    print(run)
    ad_int = adata_MG[adata_MG.obs['run'] == run]
    with plt.rc_context({'figure.figsize': (20, 10)}):
        sc.pl.spatial(ad_int, spot_size=25, color = 'MG_leiden_1')
    plt.show()



# In[ ]:




