#!/usr/bin/env python
# coding: utf-8

# In[40]:


import scanpy as sc
import decoupler as dc
import pertpy as pt


sc.set_figure_params(figsize=(3, 3), frameon=False)


# In[41]:


adata = sc.read_h5ad('../data/mtDNA_DSB_5k_clustered_manual_annotation_with_raw.h5ad')


# In[42]:


adata.X = adata.layers['counts']
adata.X = adata.X.astype(int)


# In[43]:


adata = adata[adata.obs.cell_class.str.contains('ligo')]
adata = adata[~adata.obs.cell_class.str.contains('mic')]


# In[ ]:


adata.obs_names_make_unique()


# In[44]:


sc.pl.umap(adata, color=["cell_class", "condition"], ncols=1)


# ## Pseudobulking
# The pseudo-bulk approach involves the following steps:
# 
# Subsetting the cell type of interest
# Extracting their raw integer counts
# Summing their counts per gene into a single profile if they pass quality control
# Then, DEA can be performed if there are at least two biological replicates per condition (more replicates are recommended).
# 
# Pseudobulking can easily be performed using the function {func}decoupler.pp.pseudobulk. In this example, the counts are just summed, though other modes such as the mean or any custom aggregation function are available. For more information, refer to the mode argument.

# In[45]:


pdata = dc.pp.pseudobulk(
    adata=adata,
    sample_col="sample_id",
    groups_col="cell_class",
    mode="sum",
)


# In[46]:


pdata.obs["age"]


# In[47]:


# make sure these are in pdata.obs (one row = one biological sample / pseudobulk)
pdata.obs["condition"] = pdata.obs["condition"].astype("category")   # e.g. ["control","mtDSB"]
pdata.obs["condition"] = pdata.obs["condition"].cat.reorder_categories(["control","mtDSB"], ordered=True)

# age as categorical (p21/p60)
pdata.obs["age"] = pdata.obs["age"].astype("category")       # e.g. ["p21","p60"]
pdata.obs["age"] = pdata.obs["age"].cat.reorder_categories(["21","60"], ordered=True)

# OR age as continuous (if you have numeric ages)
if "age" in pdata.obs:
    pdata.obs["age"] = pdata.obs["age"].astype(float)
    pdata.obs["age_z"] = (pdata.obs["age"] - pdata.obs["age"].mean()) / pdata.obs["age"].std()


# In[48]:


# Store raw counts in layers
pdata.layers["counts"] = pdata.X.copy()

# Normalize, scale and compute pca
sc.pp.normalize_total(pdata, target_sum=1e4)
sc.pp.log1p(pdata)
sc.pp.scale(pdata, max_value=10)
sc.tl.pca(pdata)

# Return raw counts to X
dc.pp.swap_layer(adata=pdata, key="counts", inplace=True)


# In[49]:


import warnings
warnings.filterwarnings("ignore")


# In[50]:


pdata.obs["age"] = (
    pdata.obs["age"]
    .astype(str)
    .replace({"21.0": "p21", "60.0": "p60"})
    .astype("category")
    .cat.reorder_categories(["p21","p60"], ordered=True)
)


# In[51]:


pdata.obs["age"] = pdata.obs["age"].astype(str)          # convert numbers → strings
pdata.obs["age"] = pdata.obs["age"].astype("category")   # convert → categorical
pdata.obs["age"] = pdata.obs["age"].cat.reorder_categories(["p21","p60"], ordered=True)


# In[52]:


cols = pds.formulaic_contrasts.design_matrix.model_spec.column_names
print(cols)


# In[53]:


import pertpy as pt

pds = pt.tl.PyDESeq2(adata=pdata, design="~ age + condition", layer="counts")  # drop layer=... if X has counts
pds.fit()

# condition effect (mtDSB vs control), adjusted for age
res_df = pds.test_contrasts(pds.contrast(column="condition", baseline="control", group_to_compare="mtDSB"))


# In[54]:


# 1️⃣ condition effect (mtDSB vs control)
res_condition = pds.test_contrasts(
    pds.contrast(column="condition", baseline="control", group_to_compare="mtDSB")
)

# 2️⃣ age effect (60 vs 21)
# NOTE: use the exact strings that exist in your pdata.obs['age'] column (likely "21.0" / "60.0")
res_age = pds.test_contrasts(
    pds.contrast(column="age", baseline="p21", group_to_compare="p60")
)

# 3️⃣ interaction (does mtDSB effect differ by age)
# model with interaction
pds = pt.tl.PyDESeq2(adata=pdata, design="~ age * condition", layer="counts")
pds.fit()

res_df = pds.test_contrasts(
    pds.contrast(column="condition", baseline="control", group_to_compare="mtDSB")
)
import numpy as np

# the exact column names in your model
cols = pds.formulaic_contrasts.design_matrix.model_spec.column_names
print(cols)
# ('Intercept', 'age[T.p60]', 'condition[T.mtDSB]', 'age[T.p60]:condition[T.mtDSB]')

# build the contrast vector for the interaction term
contrast_vector = np.zeros(len(cols))
contrast_vector[cols.index('age[T.p60]:condition[T.mtDSB]')] = 1
print(contrast_vector)  # should now show [0., 0., 0., 1.]

# run the test
res_interaction = pds.test_contrasts(contrast_vector)
res_interaction.head()


# In[74]:


res_interaction.T['Ldhb']


# In[56]:


res_interaction = res_interaction.set_index('variable')


# In[57]:


import matplotlib.pyplot as plt
import decoupler as dc

import matplotlib.pyplot as plt
plt.rcParams["figure.figsize"] = (10, 6)    
dc.pl.volcano(
        res_interaction,
        x="log_fc",
        y="p_value",
        top=30,
    figsize=(12, 8)
    )


# # 🧬 Biological narrative
# 
# “In oligodendrocytes, mtDNA double-strand breaks elicit an age-dependent transcriptional remodeling.
# At p21, the mtDSB response engages developmental regulators and neuropeptide-like programs, whereas at p60 it shifts toward inflammatory and dedifferentiation-associated genes (Irf8, Il1a, Hoxb3, Fcgr2b).
# This suggests that aging transforms the oligodendrocyte mtDNA-damage response from a regenerative to a stress-reactive state.”
