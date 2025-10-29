#!/usr/bin/env python
# coding: utf-8

# In[11]:


import os
import numpy as np
import scanpy as sc
import squidpy as sq
from anndata import AnnData
import scipy.sparse as sp
import warnings
warnings.filterwarnings('ignore')
import os
import pandas as pd
import matplotlib.pyplot as plt
def format_data_neighs_colapse(
    adata: AnnData,
    spatial_key: str = "spatial",
    neighs: int = 10
) -> AnnData:
    """Pseudobin each cell’s expression by summing over its spatial neighbors."""
    if spatial_key not in adata.obsm:
        coords = adata.obs[['X','Y']].values.astype(float)
        adata.obsm[spatial_key] = coords

    ad = adata.copy()
    print(f"    → Computing spatial neighbors (k={neighs})…", end="", flush=True)
    sq.gr.spatial_neighbors(ad, spatial_key=spatial_key, n_neighs=neighs)
    print("done")

    C = ad.obsp['spatial_connectivities'].astype(bool).astype(int)
    C = C + sp.eye(C.shape[0], format='csr')
    X = ad.X.copy() if not isinstance(ad.X, np.ndarray) else ad.X
    Xc = C.dot(X)

    print("    → Pseudobinning complete")
    return AnnData(
        X    = Xc,
        obs  = ad.obs.copy(),
        var  = ad.var.copy(),
        obsm = ad.obsm.copy(),
        uns  = ad.uns.copy()
    )


def domains_by_rbd(
    adata: AnnData,
    hyperparams: dict,
    sample_key: str = 'sample',
    cell_id_key: str = None,
    checkpoint_dir: str = './checkpoints',
    force_rerun: bool = False
) -> (AnnData, AnnData):
    """
    Read-based domains with checkpointed saves, force-rerun option, and printouts.
    """
    os.makedirs(checkpoint_dir, exist_ok=True)

    # A) ensure unique IDs
    if cell_id_key is None:
        print("1) Ensuring unique cell IDs…", end="", flush=True)
        adata.obs_names_make_unique()
        cell_id_key = 'cell_id'
        adata.obs[cell_id_key] = adata.obs_names
        print("done\n")

    def load_or_compute(path, compute_fn, step_name):
        if (not force_rerun) and os.path.exists(path):
            print(f"{step_name}: Loading from checkpoint…", end="", flush=True)
            ad = sc.read_h5ad(path)
            print("loaded\n")
            return ad, False
        else:
            print(f"{step_name}: Computing…", flush=True)
            ad = compute_fn()
            ad.write_h5ad(path)
            print(f"{step_name}: Saved to {path}\n")
            return ad, True

    # 1) PSEUDOBINNING
    pseudobin_path = os.path.join(checkpoint_dir, 'pseudobinned.h5ad')
    def compute_pseudobin():
        parts = []
        for samp in adata.obs[sample_key].unique():
            print(f"  Sample '{samp}':")
            sub = adata[adata.obs[sample_key] == samp].copy()
            parts.append(format_data_neighs_colapse(
                sub, spatial_key='spatial', neighs=hyperparams['neighbors']
            ))
        return sc.concat(parts, merge='same')

    pseudobinned, _ = load_or_compute(
        pseudobin_path, compute_pseudobin, "Step 1) Pseudobinning"
    )

    # 2) FILTERING
    filtered_path = os.path.join(checkpoint_dir, 'filtered.h5ad')
    def compute_filtered():
        ad_neigh = pseudobinned.copy()
        ad_neigh.X = np.nan_to_num(ad_neigh.X)
        sc.pp.filter_cells(ad_neigh, min_counts=4)
        ad_neigh.raw = ad_neigh
        return ad_neigh

    adata_neigh, _ = load_or_compute(
        filtered_path, compute_filtered, "Step 2) Filtering"
    )

    # 3) UMAP
    umap_path = os.path.join(checkpoint_dir, 'umap.h5ad')
    def compute_umap():
        nn  = hyperparams.get('n_neighbors', 20)
        pcs = hyperparams.get('n_pcs', None)
        if pcs is None:
            sc.pp.neighbors(adata_neigh, n_neighbors=nn, use_rep='X')
        else:
            sc.tl.pca(adata_neigh, n_comps=pcs, svd_solver='arpack')
            sc.pp.neighbors(adata_neigh, n_neighbors=nn, n_pcs=pcs)
        sc.tl.umap(adata_neigh, min_dist=hyperparams.get('min_dist',0.1))
        return adata_neigh

    adata_neigh, _ = load_or_compute(
        umap_path, compute_umap, "Step 3) UMAP"
    )

    # 4) CLUSTERING at one or many resolutions
    algo = hyperparams['clustering_algorithm'].lower()
    res_list = hyperparams.get('resolutions',
               [hyperparams.get('resolution')])

    for res in res_list:
        key = f"rbd_domain_{res}"
        chk = os.path.join(checkpoint_dir, f"clustered_{res}.h5ad")

        def compute_cluster():
            if algo == 'leiden':
                sc.tl.leiden(adata_neigh, resolution=res, key_added=key)
            else:
                sc.tl.louvain(adata_neigh, resolution=res, key_added=key)
            return adata_neigh

        adata_neigh, _ = load_or_compute(
            chk, compute_cluster,
            f"Step 4) Clustering (res={res})"
        )

        # map back to original
        print(f"  → Mapping rbd_domain_{res} back to original…", end="", flush=True)
        mapping = dict(zip(
            adata_neigh.obs[cell_id_key],
            adata_neigh.obs[key]
        ))
        adata.obs[key] = adata.obs[cell_id_key].map(mapping).astype(str)
        print("done\n")

    print("All steps complete! Returning results.")
    return adata, adata_neigh

hyperparams = {
    'neighbors': 20,
    'clustering_algorithm': 'leiden',
    'resolutions': [0.1, 0.2, 0.5, 1.0],
    'n_neighbors': 20,
    'n_pcs': 30,
    'min_dist': 0.5
}


# In[2]:


adata = sc.read_h5ad('../data/mtDNA_DSB_5k_clustered_manual_annotation.h5ad')


# In[3]:


adata_sub = adata[adata.obs['sample_id'] == 'RB4405']


# In[4]:


adata_sub.obs['cell_id'] = adata_sub.obs.index


# In[12]:


# 3) Call domains_by_rbd, forcing a fresh run
adata_annotated, adata_pseudobinned = domains_by_rbd(
    adata_sub,
    hyperparams,
    sample_key='sample_id',
    cell_id_key='cell_id',               # or whatever your unique-ID column is
    checkpoint_dir='../results/rbd_checkpoints',  # where to store (or load) .h5ad files
    force_rerun=True                     # overwrite any existing checkpoints
)


# In[13]:


# 5) Plot a spatial map for your favorite resolution:
with plt.rc_context({"figure.figsize": (10,10)}):
    sc.pl.spatial(
        adata_pseudobinned,
        color='rbd_domain_0.1',
        spot_size=15,
        title="Read-based domains (resolution=0.1)"
    )


# In[14]:


# 5) Plot a spatial map for your favorite resolution:
with plt.rc_context({"figure.figsize": (10,10)}):
    sc.pl.spatial(
        adata_pseudobinned,
        color='rbd_domain_0.2',
        spot_size=15,
        title="Read-based domains (resolution=0.2)"
    )


# In[15]:


# 5) Plot a spatial map for your favorite resolution:
with plt.rc_context({"figure.figsize": (10,10)}):
    sc.pl.spatial(
        adata_pseudobinned,
        color='rbd_domain_0.2',
        spot_size=15,
        title="Read-based domains (resolution=0.5)"
    )


# In[16]:


# 5) Plot a spatial map for your favorite resolution:
with plt.rc_context({"figure.figsize": (10,10)}):
    sc.pl.spatial(
        adata_pseudobinned,
        color='rbd_domain_0.5',
        spot_size=15,
        title="Read-based domains (resolution=0.5)"
    )


# In[17]:


# 5) Plot a spatial map for your favorite resolution:
with plt.rc_context({"figure.figsize": (10,10)}):
    sc.pl.spatial(
        adata_pseudobinned,
        color='rbd_domain_1.0',
        spot_size=15,
        title="Read-based domains (resolution=1.0)"
    )


# In[ ]:




