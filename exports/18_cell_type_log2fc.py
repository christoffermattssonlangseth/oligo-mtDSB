#!/usr/bin/env python
# coding: utf-8

# In[7]:


import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import warnings
warnings.filterwarnings("ignore")


# In[8]:


ad_raw = sc.read_h5ad('../data/mtDNA_DSB_5k_raw.h5ad')


# In[10]:


ad = sc.read_h5ad('../data/mtDNA_DSB_5k_clustered_manual_annotation.h5ad')


# In[11]:


import pandas as pd

# contingency table
ct_counts = pd.crosstab(ad.obs["cell_class"], [ad.obs["condition"], ad.obs["age"]])

# normalize within condition/age
ct_frac = ct_counts.div(ct_counts.sum(axis=0), axis=1)

# difference mtDSB - control at each age
diff = ct_frac.xs("mtDSB", level="condition", axis=1) - ct_frac.xs("control", level="condition", axis=1)

print(diff.sort_values(by="21"))  # or "60"


# In[12]:


import scanpy as sc
import numpy as np

effect_size = {}

for ct in ad.obs["cell_type"].unique():
    sub = ad[ad.obs["cell_type"] == ct]
    if sub.n_obs < 50:  # skip small groups
        continue

    sc.tl.pca(sub, svd_solver="arpack")
    # condition separation on PC1
    from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
    lda = LDA(n_components=1)
    X = sub.obsm["X_pca"][:, :10]
    y = sub.obs["condition"].values
    score = lda.fit(X, y).score(X, y)
    effect_size[ct] = 1 - score  # higher = more affected


# In[ ]:


effect_size


# In[13]:


from collections import Counter
import scanpy as sc

de_counts = {}

for ct in ad.obs["cell_class"].unique():
    print(ct)
    sub = ad[ad.obs["cell_class"] == ct]
    if sub.n_obs < 50:
        continue
    sc.tl.rank_genes_groups(sub, groupby="condition", method="t-test")
    sc.pl.rank_genes_groups(sub)
    res = sc.get.rank_genes_groups_df(sub, group="mtDSB")
    de_counts[ct] = (res["pvals_adj"] < 0.05).sum()

# cell types with most DE genes
sorted(de_counts.items(), key=lambda x: x[1], reverse=True)


# In[22]:


sub = ad[ad.obs["cell_class"].isin(['Mature oligodendrocytes','Oligodendrocytes precursor cells'])]


# In[23]:


sub.obs['condition-age'] = sub.obs['condition'].astype(str) + '-'+ sub.obs['age'].astype(str)


# In[24]:


sc.tl.rank_genes_groups(sub, groupby="condition-age", method="t-test")
sc.pl.rank_genes_groups(sub)


# In[17]:


mtDSB_genes = [
    "Hspa5", "Hspa9", "Hsph1",
    "Atf5", "Trib3", "Zbtb16", "Ddit3",
    "Cdkn1a", "Bcl2l1",
    "Sgk1", "Nmu", "Plin4", "Aldoa",
    "Serpina3n", "Mt2", "Gstp1",'Ldha'
]


# In[25]:


sc.pl.dotplot(
        sub,
        var_names=mtDSB_genes,
        groupby="condition-age",
        standard_scale="var",
        #dot_max=0.5,
        #dot_min=0.05,
        color_map="Reds",
        dendrogram=True,
        figsize=(8, 2)
    )


# In[44]:


sc.pl.dotplot(
        sub,
        var_names=mtDSB_genes,
        groupby="sample_id",
        standard_scale="var",
        #dot_max=0.5,
        #dot_min=0.05,
        color_map="Reds",
        dendrogram=True,
        figsize=(8, 2)
    )


# In[19]:


sub.obs["condition-age"]


# In[26]:


import numpy as np
import pandas as pd

genes = mtDSB_genes  # your list

# subset to oligos + opcs

logfc = {}
for g in genes:
    if g not in sub.var_names:
        continue

    # get raw expression (use layer="counts" if you have raw counts stored)
    expr = sub[:, g].X
    if hasattr(expr, "toarray"):
        expr = expr.toarray().ravel()
    else:
        expr = np.asarray(expr).ravel()

    cond = sub.obs["condition-age"].values
    mean_ctrl = expr[cond == "control-21"].mean()  # add small offset
    mean_ds   = expr[cond == "mtDSB-21"].mean()

    logfc[g] = np.log2(mean_ds / mean_ctrl)

df_logfc = pd.Series(logfc).sort_values(ascending=False)
print(df_logfc)


# In[27]:


import numpy as np
import pandas as pd

genes = mtDSB_genes  # your list

# subset to oligos + opcs

logfc = {}
for g in genes:
    if g not in sub.var_names:
        continue

    # get raw expression (use layer="counts" if you have raw counts stored)
    expr = sub[:, g].X
    if hasattr(expr, "toarray"):
        expr = expr.toarray().ravel()
    else:
        expr = np.asarray(expr).ravel()

    cond = sub.obs["condition-age"].values
    mean_ctrl = expr[cond == "control-60"].mean()  # add small offset
    mean_ds   = expr[cond == "mtDSB-60"].mean()

    logfc[g] = np.log2(mean_ds / mean_ctrl)

df_logfc = pd.Series(logfc).sort_values(ascending=False)
print(df_logfc)


# In[28]:


for i in sub.obs.age.unique():
    sub_s = sub[sub.obs.age == i]
    sc.tl.rank_genes_groups(sub_s, groupby="condition-age", method="t-test")
    sc.pl.rank_genes_groups(sub_s)


# In[30]:


results = []
for age in sub.obs["age"].unique():
    for cond in sub.obs["condition"].unique():
        sub_ = sub[(sub.obs["age"] == age) & (sub.obs["condition"] == cond)]
        for g in genes:
            if g not in sub_.var_names:
                continue
            expr = sub_[:, g].X
            if hasattr(expr, "toarray"):
                expr = expr.toarray().ravel()
            else:
                expr = np.asarray(expr).ravel()
            results.append({
                "gene": g,
                "age": age,
                "condition": cond,
                "mean_expr": expr.mean()
            })

df_means = pd.DataFrame(results)

# pivot to get logFC for each age
df_pivot = df_means.pivot_table(index="gene", columns=["age","condition"], values="mean_expr")
df_pivot["logFC_mtDSB_vs_ctrl_21"] = np.log2((df_pivot['21',"mtDSB"]+1e-6)/(df_pivot['21',"control"]+1e-6))
df_pivot["logFC_mtDSB_vs_ctrl_60"] = np.log2((df_pivot['60',"mtDSB"]+1e-6)/(df_pivot['60',"control"]+1e-6))


# In[31]:


df_pivot


# In[43]:


for run in ad.obs['sample_id'].unique():
    print(run)
    ad_int = ad[ad.obs['sample_id'] == run]

    with plt.rc_context({'figure.figsize': (20, 10)}):
        sc.pl.spatial(ad_int, spot_size=15, color = 'Ldha', vmax = 1)
    plt.show()



# In[ ]:




