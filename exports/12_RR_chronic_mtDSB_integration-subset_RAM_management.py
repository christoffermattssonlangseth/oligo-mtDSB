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

# 4) No scaling (keeps sparsity). Compute PCs via TruncatedSVD (no centering; OK for huge sparse)
n_comps = 50
svd = TruncatedSVD(n_components=n_comps, random_state=0)
X_pca = svd.fit_transform(ad.X)  # returns dense (cells × PCs), but only 50 dims
ad.obsm["X_pca"] = X_pca.astype(np.float32)
ad.uns["pca"] = {
    "variance_ratio": svd.explained_variance_ratio_.astype(np.float32),
    "variance": (svd.explained_variance_ratio_ * svd.explained_variance_ratio_.sum()).astype(np.float32)  # placeholder
}

# 5) Subsample for graph building
rng = np.random.default_rng(0)
sub_n = min(200_000, ad.n_obs)  # adjust if you have more/less RAM
sub_idx = rng.choice(ad.n_obs, size=sub_n, replace=False)
sub_idx.sort()
mask_sub = np.zeros(ad.n_obs, dtype=bool)
mask_sub[sub_idx] = True

# 6) Neighbors + Leiden on the subsample only (fast & RAM-light)
ad_sub = ad[mask_sub].copy()
ad_sub.obsm["X_pca"] = ad.obsm["X_pca"][mask_sub]
sc.pp.neighbors(ad_sub, use_rep="X_pca", n_neighbors=30, method="umap")


# In[8]:


sc.tl.leiden(ad_sub, resolution=2.5, key_added="leiden_sub")

# 7) Map all cells to the subsample: label transfer via kNN in PCA space
#    Fit index on subsample PCs and query all cells (streamable, RAM-friendly)
nbrs = NearestNeighbors(n_neighbors=5, algorithm="auto", n_jobs=-1)
nbrs.fit(ad_sub.obsm["X_pca"])

# For very tight memory, query in chunks
batch = 200_000
labels_sub = ad_sub.obs["leiden_sub"].to_numpy()
all_labels = np.empty(ad.n_obs, dtype=labels_sub.dtype)

for start in range(0, ad.n_obs, batch):
    stop = min(start + batch, ad.n_obs)
    dists, inds = nbrs.kneighbors(ad.obsm["X_pca"][start:stop])
    # majority vote from nearest subsample neighbors
    voted = np.array([
        np.argmax(np.bincount(labels_sub[inds[i]].astype(int)))
        for i in range(inds.shape[0])
    ], dtype=labels_sub.dtype)
    all_labels[start:stop] = voted

ad.obs["leiden"] = all_labels  # cluster for ALL 1.7M cells

# 8) (Optional) UMAP on subsample only, then embed the rest by kNN-averaging
sc.tl.umap(ad_sub, min_dist=0.3, spread=1.0)
ad_sub.obsm["X_umap"] = ad_sub.obsm["X_umap"].astype(np.float32)

# project full set to UMAP by averaging neighbors' coordinates
umap_full = np.empty((ad.n_obs, 2), dtype=np.float32)
for start in range(0, ad.n_obs, batch):
    stop = min(start + batch, ad.n_obs)
    _, inds = nbrs.kneighbors(ad.obsm["X_pca"][start:stop])
    umap_full[start:stop] = ad_sub.obsm["X_umap"][inds].mean(axis=1)
ad.obsm["X_umap"] = umap_full


# In[9]:


import leidenalg


# In[10]:


ad.obs.leiden = ad.obs.leiden.astype(str)


# In[11]:


ad.write('../data/integrated_clustered.h5ad')


# In[12]:


spatial = np.array(ad.obs[['x_centroid','y_centroid']])
ad.obsm['spatial'] = spatial


# In[13]:


sc.pl.umap(ad, color = 'leiden')


# In[14]:


# --- Find marker genes per cluster ---
sc.tl.rank_genes_groups(ad, groupby="leiden", method="t-test")
sc.pl.rank_genes_groups(ad, n_genes=25, sharey=False)

# get a tidy table of top markers
markers = sc.get.rank_genes_groups_df(ad, group=None)
markers.head()


# In[15]:


marker_genes = pd.DataFrame({
    group: ad.uns['rank_genes_groups']['names'][group][:20]
    for group in ad.uns['rank_genes_groups']['names'].dtype.names
})


# In[16]:


for col in marker_genes.columns: 
    print(col)
    genes = marker_genes[col].tolist()
    print(" ".join(genes))
    print(' ')


# In[17]:


ad.obs


# In[18]:


ad_mtDSB = ad[ad.obs.project == 'mtDSB']
for run in ad_mtDSB.obs['sample_extracted'].unique():
    print(run)
    ad_int = ad_mtDSB[ad_mtDSB.obs['sample_extracted'] == run]
    with plt.rc_context({'figure.figsize': (20, 10)}):
        sc.pl.spatial(ad_int, spot_size=15, color = 'leiden')
    plt.show()



# In[19]:


ad_RR = ad[ad.obs.project == 'RREAE']
for run in ad_RR.obs['sample_extracted'].unique():
    print(run)
    ad_int = ad_RR[ad_RR.obs['sample_extracted'] == run]
    with plt.rc_context({'figure.figsize': (20, 10)}):
        sc.pl.spatial(ad_int, spot_size=15, color = 'leiden')
    plt.show()



# In[20]:


ad_chronic = ad[ad.obs.project == 'Chronic EAE']
for run in ad_chronic.obs['sample_extracted'].unique():
    print(run)
    ad_int = ad_chronic[ad_chronic.obs['sample_extracted'] == run]
    with plt.rc_context({'figure.figsize': (20, 10)}):
        sc.pl.spatial(ad_int, spot_size=15, color = 'leiden')
    plt.show()



# In[ ]:




