#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import warnings
warnings.filterwarnings("ignore")


# In[2]:


adata = sc.read_h5ad('../data/mtDNA_DSB_5k_clustered.h5ad')


# In[3]:


# --- Find marker genes per cluster ---
sc.tl.rank_genes_groups(adata, groupby="leiden_2", method="t-test")
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


adata.obs['cell_class'] = adata.obs['leiden_2'].map(annotations)


# In[13]:


sc.pl.umap(adata, color='cell_class', legend_loc='on data', frameon=False)


# In[14]:


spatial = np.array(adata.obs[['x_centroid','y_centroid']])
adata.obsm['spatial'] = spatial


# In[22]:


for run in adata.obs['run'].unique():
    print(run)
    ad_int = adata[adata.obs['run'] == run]
    with plt.rc_context({'figure.figsize': (20, 10)}):
        sc.pl.spatial(ad_int, spot_size=15, color = 'leiden_2')
    plt.show()



# In[30]:


adata.write('../data/mtDNA_DSB_5k_clustered_LLM_anno.h5ad')


# In[ ]:




