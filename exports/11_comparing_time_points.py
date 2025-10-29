#!/usr/bin/env python
# coding: utf-8

# In[1]:


import scanpy as sc


# In[2]:


adata = sc.read_h5ad('../data/mtDNA_DSB_5k_clustered_LLM_anno_with_counts.h5ad')


# In[3]:


adata = adata[adata.obs.cell_class.str.contains('ligo')]


# In[7]:


sc.pl.dotplot(
        adata,
        var_names=['Serpina3n','Ldha','Mog','Mal','Mbp', 'Aldoa','Pkm'],
        groupby="condition",
        standard_scale="var",
        #dot_max=0.5,
        #dot_min=0.05,
        color_map="Reds",
        #categories_order = order,
        dendrogram=False,
        figsize=(7, 5)
    )


# In[8]:


# --- Find marker genes per cluster ---
sc.tl.rank_genes_groups(adata, groupby="condition", method="t-test")
sc.pl.rank_genes_groups(adata, n_genes=30, sharey=False, fontsize = 10)

# get a tidy table of top markers
markers = sc.get.rank_genes_groups_df(adata, group=None)
markers.head()


# In[9]:


adata.obs['cond_age'] = adata.obs['condition'].astype(str) + '_' +adata.obs['age'].astype(str)


# In[10]:


sc.pl.dotplot(
        adata,
        var_names=['Serpina3n','Ldha','Mog','Mal','Mbp', 'Aldoa','Pkm'],
        groupby="cond_age",
        standard_scale="var",
        #dot_max=0.5,
        #dot_min=0.05,
        color_map="Reds",
        #categories_order = order,
        dendrogram=False,
        figsize=(4, 4)
    )


# In[11]:


# --- Find marker genes per cluster ---
sc.tl.rank_genes_groups(adata, groupby="cond_age", method="t-test")
sc.pl.rank_genes_groups(adata, n_genes=30, sharey=False, fontsize = 10)

# get a tidy table of top markers
markers = sc.get.rank_genes_groups_df(adata, group=None)
markers.head()


# In[12]:


markers


# In[13]:


# --- Find marker genes per cluster ---
sc.tl.rank_genes_groups(adata, groupby="condition", method="t-test")
sc.pl.rank_genes_groups(adata, n_genes=30, sharey=False, fontsize = 10)

# get a tidy table of top markers
markers = sc.get.rank_genes_groups_df(adata, group=None)
markers.head()


# In[14]:


markers[markers['names'] == 'Ldha']


# In[15]:


for age in adata.obs['age'].unique():
    print(age)
    adata_sub = adata[adata.obs['age'] == age]
    # --- Find marker genes per cluster ---
    sc.tl.rank_genes_groups(adata_sub, groupby="condition", method="t-test")
    sc.pl.rank_genes_groups(adata_sub, n_genes=30, sharey=False, fontsize = 10)

    # get a tidy table of top markers
    markers = sc.get.rank_genes_groups_df(adata_sub, group=None)
    markers.head()


# In[18]:


import scanpy as sc
import pandas as pd

ad = adata.copy()

# make sure metadata columns are clean
ad.obs['age'] = ad.obs['age'].astype(int)
ad.obs['condition'] = ad.obs['condition'].astype(str).str.strip().str.lower()

# build a new AnnData with counts as the "raw" slot
ad_raw = sc.AnnData(X=ad.layers["counts"], obs=ad.obs.copy(), var=ad.var.copy())
ad.raw = ad_raw

def de_condition_within_age(ad, age_val):
    """mtDSB vs control at a given age."""
    ad_sub = ad[ad.obs['age'] == int(age_val)].copy()
    ad_sub = ad_sub[ad_sub.obs['condition'].isin(['control','mtdsb'])].copy()
    ad_sub.obs['condition'] = pd.Categorical(
        ad_sub.obs['condition'], categories=['control','mtdsb']
    )
    sc.tl.rank_genes_groups(
        ad_sub,
        groupby='condition',
        reference='control',
        method='wilcoxon',
        use_raw=True
    )
    df = sc.get.rank_genes_groups_df(ad_sub, group=None)
    df['contrast'] = f"mtDSB_vs_control@{age_val}"
    return df

# Run DE
de_21 = de_condition_within_age(ad, 21)
de_60 = de_condition_within_age(ad, 60)

print("Top DE genes at 21d:")


# In[23]:


# Example metabolic modules
glycolysis = ["Gapdh","Eno1","Pkm","Ldha","Pgk1","Pfkp","Aldoa","Tpi1","Gpi1","Pdk1"]
tca_cycle = ["Idh1","Idh2","Idh3a","Idh3b","Idh3g","Sdha","Sdhb","Sdhc","Ogdh","Mdh2"]
oxphos = ["Ndufs1","Ndufs4","Ndufv2","Uqcrc1","Uqcrfs1","Cox4i1","Cox7a1","Atp5f1","Atp5o"]
lactate_shuttle = ["Slc16a1","Slc16a3","Ldha","Ldhb"]
lipid_metabolism = ["Abca1","Hmgcr","Apoe","Soat1","Cyp46a1"]


# In[24]:


def filter_metabolic(de_df, gene_lists):
    metabolic_genes = set().union(*gene_lists.values())
    return de_df[de_df["names"].isin(metabolic_genes)].copy()

gene_sets = {
    "Glycolysis": glycolysis,
    "TCA cycle": tca_cycle,
    "OXPHOS": oxphos,
    "Lactate shuttle": lactate_shuttle,
    "Lipid metabolism": lipid_metabolism
}

# Example: check DE in 21d mtDSB vs control
met_21 = filter_metabolic(de_21, gene_sets)
met_60 = filter_metabolic(de_60, gene_sets)

print("21d metabolic DE genes:\n", met_21[["names","logfoldchanges","pvals_adj"]])
print("60d metabolic DE genes:\n", met_60[["names","logfoldchanges","pvals_adj"]])


# In[27]:


print(met_60.sort_values(by = 'logfoldchanges', ascending = False))


# In[ ]:


ad_60.obs.cond_age.unique()


# In[ ]:




