#!/usr/bin/env python
# coding: utf-8

# In[1]:


import warnings
warnings.filterwarnings('ignore')
import os
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.cluster import KMeans
import numpy as np
import scanpy as sc
from scipy import sparse
from sklearn.decomposition import TruncatedSVD
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import NearestNeighbors
import re


# In[2]:


adata_RRmap = sc.read_h5ad('../data/RREAE_5k_raw.h5ad')
adata_EAE5K = sc.read_h5ad('../data/eae_5k_raw.h5ad')
adata_mtDSB = sc.read_h5ad('../data/mtDNA_DSB_5k_raw.h5ad')


# In[3]:


adata_EAE5K.obs["sample_extracted"] = adata_EAE5K.obs["run"].str.extract(r"__(G\d+_L\d+)__")
adata_RRmap.obs["sample_extracted"] = adata_RRmap.obs["sample"].str.extract(r"__(S\d+-[A-Za-z]+\d+)__")
adata_mtDSB.obs["sample_extracted"] = adata_mtDSB.obs["run"].str.extract(r"__(RB\d+)__")


# In[4]:


adata_mtDSB.obs_names_make_unique()
adata_RRmap.obs_names_make_unique()
adata_EAE5K.obs_names_make_unique()


# In[5]:


adata_mtDSB.obs['project'] = 'mtDSB'
adata_RRmap.obs['project'] = 'RREAE'
adata_EAE5K.obs['project'] = 'Chronic EAE'
ad = sc.concat(
    {
        "RRmap": adata_RRmap,
        "EAE5K": adata_EAE5K,
        "mtDSB": adata_mtDSB,
    },
    label="dataset",       # new column in .obs
    keys=None,             # use dict keys above
    index_unique="-",      # make indices unique: barcode-dataset
    join="outer",          # keep all genes (union)
    merge="same"           # only merge .obs/.var columns that are identical
)


# In[6]:


ad.write('../data/integrated_2.h5ad')


# In[7]:


# 2) (Optional) light filtering
sc.pp.filter_cells(ad, min_counts=40)      # tune if needed
sc.pp.filter_genes(ad, min_cells=5)

# 3) Normalize & log1p (sparse-friendly in modern Scanpy)
sc.pp.normalize_total(ad, target_sum=1e4)
sc.pp.log1p(ad)


# In[8]:


plt.rcdefaults()
sc.tl.pca(ad)
sc.pl.pca_variance_ratio(ad, n_pcs=50, log=True)
sc.pp.neighbors(ad, n_neighbors=15, n_pcs=30)


# In[9]:


sc.tl.umap(ad, min_dist=0.1)


# In[21]:


resolutions = [0.5, 1,1.5, 2]

for resolution in resolutions:
    key = f'leiden_{resolution}'

    if key in ad.obs.columns:
        print(f"Skipping {resolution}: {key} already exists.")
    else:
        print(f"Clustering at resolution {resolution}...")
        sc.tl.leiden(ad, resolution=resolution, key_added=key)
        print("Done.")

    # plot UMAP
    sc.pl.umap(ad, color=key, legend_loc='on data', frameon=False)


# In[34]:


ad


# In[36]:


sc.pl.umap(ad, color = ['leiden_1.5','Serpina3n', 'project'], legend_loc='on data', frameon=False)


# In[32]:


ad.write('../data/integrated_clustered_full.h5ad')


# In[23]:


spatial = np.array(ad.obs[['x_centroid','y_centroid']])
ad.obsm['spatial'] = spatial


# In[25]:


# --- Find marker genes per cluster ---
sc.tl.rank_genes_groups(ad, groupby="leiden_1.5", method="t-test")
sc.pl.rank_genes_groups(ad, n_genes=25, sharey=False)

# get a tidy table of top markers
markers = sc.get.rank_genes_groups_df(ad, group=None)
markers.head()


# In[26]:


marker_genes = pd.DataFrame({
    group: ad.uns['rank_genes_groups']['names'][group][:20]
    for group in ad.uns['rank_genes_groups']['names'].dtype.names
})


# In[27]:


for col in marker_genes.columns: 
    print(col)
    genes = marker_genes[col].tolist()
    print(" ".join(genes))
    print(' ')


# In[28]:


ad.obs


# In[29]:


ad_mtDSB = ad[ad.obs.project == 'mtDSB']
for run in ad_mtDSB.obs['sample_extracted'].unique():
    print(run)
    ad_int = ad_mtDSB[ad_mtDSB.obs['sample_extracted'] == run]
    with plt.rc_context({'figure.figsize': (20, 10)}):
        sc.pl.spatial(ad_int, spot_size=15, color = 'leiden_1.5')
    plt.show()



# In[30]:


ad_RR = ad[ad.obs.project == 'RREAE']
for run in ad_RR.obs['sample_extracted'].unique():
    print(run)
    ad_int = ad_RR[ad_RR.obs['sample_extracted'] == run]
    with plt.rc_context({'figure.figsize': (20, 10)}):
        sc.pl.spatial(ad_int, spot_size=15, color = 'leiden_1.5')
    plt.show()



# In[31]:


ad_chronic = ad[ad.obs.project == 'Chronic EAE']
for run in ad_chronic.obs['sample_extracted'].unique():
    print(run)
    ad_int = ad_chronic[ad_chronic.obs['sample_extracted'] == run]
    with plt.rc_context({'figure.figsize': (20, 10)}):
        sc.pl.spatial(ad_int, spot_size=15, color = 'leiden_1.5')
    plt.show()



# In[ ]:




