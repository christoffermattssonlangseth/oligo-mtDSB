#!/usr/bin/env python
# coding: utf-8

# In[1]:


mitokine_like = [
    "Gdf15", "Adm", "Cst7", "Igfbp3", "Serpina3n",
    "Cdkn1a", "Maff", "Eif4ebp1", "Aldh1l2"
]


# In[ ]:


import pandas as pd

records = []
for cell, pack in results.items():
    se = pack.get("simple_effects", {})
    for age in ("21", "60"):
        if age in se and isinstance(se[age], pd.DataFrame):
            res = se[age]
            for gene in mitokine_like:
                if gene in res.index:
                    r = res.loc[gene]
                    records.append({
                        "cell_class": cell,
                        "age": age,
                        "gene": gene,
                        "log2FC": r.get("log2FC_mtDSB_vs_control", None),
                        "padj": r.get("pvalue_mtDSB_vs_control", None),
                        "mean_cpm": r.get("mean_cpm_all", None)
                    })
df = pd.DataFrame(records)

