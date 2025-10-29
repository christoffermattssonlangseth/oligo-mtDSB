#!/usr/bin/env python
# coding: utf-8

# In[58]:


import warnings

import decoupler as dc
import pertpy as pt
import scanpy as sc

warnings.filterwarnings("ignore")


# In[59]:


adata = sc.read_h5ad('../data/mtDNA_DSB_5k_clustered_manual_annotation_with_raw.h5ad')


# In[60]:


adata


# In[4]:


pdata = dc.pp.pseudobulk(adata, sample_col="sample_id", groups_col="cell_type", layer="counts", mode="sum")
dc.pp.filter_samples(pdata, inplace=True)
pdata


# In[5]:


dc.pl.filter_samples(pdata, groupby=["sample_id", "cell_type"], figsize=(12, 4))


# In[6]:


pdata.layers["counts"] = pdata.X.copy()

sc.pp.normalize_total(pdata, target_sum=1e4)
sc.pp.log1p(pdata)
sc.pp.scale(pdata, max_value=10)
sc.pp.pca(pdata)


# In[7]:


# Return raw counts to X
dc.pp.swap_layer(pdata, "counts", inplace=True)


# In[8]:


sc.pl.pca(pdata, color=["sample_id", "age", "sex", "condition"], ncols=1, size=300)
sc.pl.pca_variance_ratio(pdata)


# In[9]:


pdata_subset = pdata.copy()
pdata_subset.obs = pdata.obs[["sample_id", "age", "sex", "condition"]]

dc.tl.rankby_obsm(
    pdata_subset,
    key="X_pca",
    uns_key="pca_anova",
)


# In[10]:


dc.pl.obsm(
    pdata_subset,
    key="pca_anova",
    names=["sample_id", "age", "sex", "condition"],
    titles=["Principle component scores", "Adjusted p-values from ANOVA"],
    cmap_obs={},
)


# In[18]:


pdata.obs


# In[46]:


pds2 = pt.tl.PyDESeq2(adata=pdata, design="~condition")


# In[47]:


pds2.fit()


# In[48]:


res_df = pds2.test_contrasts(pds2.contrast(column="condition", baseline="control", group_to_compare="mtDSB"))


# In[49]:


res_df.head(10)


# In[50]:


pds2.plot_volcano(res_df, log2fc_thresh=0)


# In[51]:


pds2.plot_fold_change(res_df, n_top_vars=15)


# In[52]:


res_df = pds2.compare_groups(pdata, column="condition", baseline="control", groups_to_compare=["mtDSB"])
edgr.plot_multicomparison_fc(res_df, figsize=(12, 1.5))


# ## with age as well

# In[53]:


pds = pt.tl.PyDESeq2(adata=pdata, design="~ age + condition")


# In[54]:


pds2.fit()


# In[56]:


res_df = pds2.test_contrasts(pds2.contrast(column="condition", baseline="control", group_to_compare="mtDSB"))


# In[57]:


pds2.plot_volcano(res_df, log2fc_thresh=0)


# In[ ]:




