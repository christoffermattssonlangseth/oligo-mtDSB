#!/usr/bin/env python
# coding: utf-8

# In[1]:


import warnings

import decoupler as dc
import pertpy as pt
import scanpy as sc

warnings.filterwarnings("ignore")


# In[2]:


adata = sc.read_h5ad('../data/mtDNA_DSB_5k_clustered_manual_annotation_with_raw.h5ad')


# In[6]:


adata = adata[adata.obs.cell_type.str.contains('ligo')]


# In[8]:


pdata = dc.pp.pseudobulk(adata, sample_col="sample_id", groups_col="cell_type", layer="counts", mode="sum")
dc.pp.filter_samples(pdata, inplace=True)
pdata


# In[9]:


dc.pl.filter_samples(pdata, groupby=["sample_id", "cell_type"], figsize=(12, 4))


# In[10]:


pdata.layers["counts"] = pdata.X.copy()

sc.pp.normalize_total(pdata, target_sum=1e4)
sc.pp.log1p(pdata)
sc.pp.scale(pdata, max_value=10)
sc.pp.pca(pdata)


# In[11]:


# Return raw counts to X
dc.pp.swap_layer(pdata, "counts", inplace=True)


# In[12]:


sc.pl.pca(pdata, color=["sample_id", "age", "sex", "condition"], ncols=1, size=300)
sc.pl.pca_variance_ratio(pdata)


# In[13]:


pdata_subset = pdata.copy()
pdata_subset.obs = pdata.obs[["sample_id", "age", "sex", "condition"]]

dc.tl.rankby_obsm(
    pdata_subset,
    key="X_pca",
    uns_key="pca_anova",
)


# In[14]:


dc.pl.obsm(
    pdata_subset,
    key="pca_anova",
    names=["sample_id", "age", "sex", "condition"],
    titles=["Principle component scores", "Adjusted p-values from ANOVA"],
    cmap_obs={},
)


# In[15]:


pdata.obs


# In[16]:


pds2 = pt.tl.PyDESeq2(adata=pdata, design="~condition")


# In[17]:


pds2.fit()


# In[18]:


res_df = pds2.test_contrasts(pds2.contrast(column="condition", baseline="control", group_to_compare="mtDSB"))


# In[19]:


res_df.head(10)


# In[27]:


pds2.plot_volcano(res_df, log2fc_thresh=0,to_label = 10 )


# In[21]:


pds2.plot_fold_change(res_df, n_top_vars=15)


# In[23]:


res_df = pds2.compare_groups(pdata, column="condition", baseline="control", groups_to_compare=["mtDSB"])


# In[32]:


res_df = res_df.set_index('variable')


# In[33]:


data = res_df[["stat"]].T.rename(index={"stat": "disease.vs.normal"})
data


# In[35]:


collectri = dc.op.collectri(organism="mouse")
collectri


# In[36]:


# Run
tf_acts, tf_padj = dc.mt.ulm(data=data, net=collectri)

# Filter by sign padj
msk = (tf_padj.T < 0.05).iloc[:, 0]
tf_acts = tf_acts.loc[:, msk]

tf_acts


# In[38]:


dc.pl.barplot(data=tf_acts, name="disease.vs.normal", figsize=(12, 8), top = 50)


# In[47]:





# In[53]:


dc.pl.network(
    net=collectri,
    data=data,
    score=tf_acts,
    sources=list(tf_acts.T.sort_values(by = 'disease.vs.normal', ascending = False).head(5).index),
    targets=5,
    figsize=(6, 6),
    vcenter=True,
    by_abs=True,
    size_node=15,
)


# In[60]:


dc.pl.volcano(
    data=res_df,
    x="log_fc",
    y="p_value",
    net=collectri,
    name="Stat1",
    top=20,
    figsize=(6, 6),
)


# In[61]:


progeny = dc.op.progeny(organism="mouse")
progeny


# In[62]:


# Run
pw_acts, pw_padj = dc.mt.ulm(data=data, net=progeny)

# Filter by sign padj
msk = (pw_padj.T < 0.05).iloc[:, 0]
pw_acts = pw_acts.loc[:, msk]

pw_acts


# In[63]:


dc.pl.barplot(data=pw_acts, name="disease.vs.normal", figsize=(3, 3))


# In[65]:


import numpy as np


# In[66]:


# Transform to df
df = pw_acts.melt(value_name="score").merge(
    pw_padj.melt(value_name="pvalue")
    .assign(logpval=lambda x: x["pvalue"].clip(2.22e-4, 1))
    .assign(logpval=lambda x: -np.log10(x["logpval"]))
)
dc.pl.dotplot(df=df, x="score", y="variable", s="logpval", c="score", scale=1, figsize=(4, 4))


# In[67]:


dc.pl.source_targets(data=res_df, x="weight", y="stat", net=progeny, name="MAPK", top=15, figsize=(8, 8))


# In[68]:


hallmark = dc.op.hallmark(organism="mouse")
hallmark


# In[69]:


# Run
hm_acts, hm_padj = dc.mt.ulm(data=data, net=hallmark)

# Filter by sign padj
msk = (hm_padj.T < 0.05).iloc[:, 0]
hm_acts = hm_acts.loc[:, msk]

hm_acts


# In[70]:


dc.pl.barplot(data=hm_acts, name="disease.vs.normal", figsize=(6, 8))


# In[73]:


# Tranform to df
df = hm_acts.melt(value_name="score").merge(
    hm_padj.melt(value_name="pvalue")
    .assign(padj=lambda x: x["pvalue"].clip(2.22e-16, 1))
    .assign(padj=lambda x: np.log10(x["pvalue"]))
)
dc.pl.dotplot(df=df, x="score", y="variable", s="padj", c="score", scale=1, figsize=(6, 6))


# In[ ]:




