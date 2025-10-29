#!/usr/bin/env python
# coding: utf-8

# In[1]:


import scanpy as sc


# In[53]:


adata = sc.read_h5ad('../data/mtDNA_DSB_5k_clustered_LLM_anno.h5ad')


# In[54]:


adata_raw = sc.read_h5ad('../data/mtDNA_DSB_5k_raw.h5ad')


# In[55]:


adata_raw = adata_raw[adata_raw.obs.index.isin(adata.obs.index)]
adata = adata[adata.obs.index.isin(adata_raw.obs.index)]


# In[56]:


adata_raw


# In[59]:


import numpy as np
import scipy.sparse as sp

# 1) Make all names unique (both objects)
adata.obs_names_make_unique();  adata.var_names_make_unique()
adata_raw.obs_names_make_unique();  adata_raw.var_names_make_unique()

# Quick sanity: check duplicates (should be empty after the calls above)
dup_cells_ad    = adata.obs_names[adata.obs_names.duplicated()]
dup_genes_ad    = adata.var_names[adata.var_names.duplicated()]
dup_cells_raw   = adata_raw.obs_names[adata_raw.obs_names.duplicated()]
dup_genes_raw   = adata_raw.var_names[adata_raw.var_names.duplicated()]
print("dups adata cells:", len(dup_cells_ad), "| genes:", len(dup_genes_ad))
print("dups raw   cells:", len(dup_cells_raw), "| genes:", len(dup_genes_raw))

# 2) Align to the intersection (shared cells & genes)
shared_cells = adata.obs_names.intersection(adata_raw.obs_names)
shared_genes = adata.var_names.intersection(adata_raw.var_names)
print(f"Shared cells: {len(shared_cells)}  Shared genes: {len(shared_genes)}")

# Subset 'adata' to the shared indices (preserve 'adata' order by reindex_like)
adata = adata[shared_cells, shared_genes].copy()

# 3) Build position indexers into adata_raw for those same names
row_idx = adata_raw.obs_names.get_indexer(adata.obs_names)
col_idx = adata_raw.var_names.get_indexer(adata.var_names)

# Ensure all found
if (row_idx < 0).any() or (col_idx < 0).any():
    missing_cells = list(adata.obs_names[np.where(row_idx < 0)[0]])
    missing_genes = list(adata.var_names[np.where(col_idx < 0)[0]])
    raise ValueError(f"Missing in adata_raw — cells: {len(missing_cells)}, genes: {len(missing_genes)}")

# 4) Slice the raw matrix and assign to a layer (shape will now match exactly)
rawX = adata_raw.X
counts_aligned = rawX[row_idx][:, col_idx] if sp.issparse(rawX) else rawX[row_idx][:, col_idx]
adata.layers["counts"] = counts_aligned.copy()  # sparse or dense OK

print("Assigned adata.layers['counts'] with shape:", adata.layers["counts"].shape)


# In[62]:


adata.write('../data/mtDNA_DSB_5k_clustered_LLM_anno_with_counts.h5ad')

