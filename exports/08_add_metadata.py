#!/usr/bin/env python
# coding: utf-8

# In[46]:


import pandas as pd
import scanpy as sc


# In[47]:


anno = pd.read_csv('/Users/christoffer/Downloads/mtDSB_anno.csv')


# In[48]:


anno


# In[49]:


adata = sc.read_h5ad('../data/mtDNA_DSB_5k_clustered_LLM_anno.h5ad')


# In[89]:


adata.write('../data/mtDNA_DSB_5k_clustered_LLM_anno.h5ad')


# In[50]:


anno.sample_id = 'RB'+anno.sample_id.astype(str)


# In[51]:


for meta in ['age', 'sex','genotype','condition']:
    mapping_dict = dict(zip(anno['sample_id'], anno[meta]))
    adata.obs[meta] = adata.obs['sample_id'].map(mapping_dict)


# In[55]:


adata.obs.age = adata.obs.age.astype(str)


# In[56]:


adata[adata.obs.age.isna()].obs.sample_id.unique()


# In[ ]:





# In[59]:


sc.pl.dotplot(
        adata,
        var_names=['Serpina3n'],
        groupby="age",
        standard_scale="var",
        #dot_max=0.5,
        #dot_min=0.05,
        color_map="Reds",
        #categories_order = order,
        dendrogram=False,
        figsize=(7, 5)
    )


# In[88]:


sc.pl.dotplot(
        adata,
        var_names=['Ldha','Serpina3n','Gfap','Pgk1','Mbp','Mog'],
        groupby="condition",
        standard_scale="var",
        #dot_max=0.5,
        #dot_min=0.05,
        color_map="Reds",
        #categories_order = order,
        dendrogram=False,
        figsize=(4, 4)
    )


# In[71]:


# --- Find marker genes per cluster ---
sc.tl.rank_genes_groups(adata, groupby="condition", method="t-test")
sc.pl.rank_genes_groups(adata, n_genes=30, sharey=False, fontsize = 10)

# get a tidy table of top markers
markers = sc.get.rank_genes_groups_df(adata, group=None)
markers.head()


# In[74]:


markers[markers.names == 'Pgk1']


# In[75]:


adata_60 = adata[adata.obs.age == '60']


# In[76]:


# --- Find marker genes per cluster ---
sc.tl.rank_genes_groups(adata_60, groupby="condition", method="t-test")
sc.pl.rank_genes_groups(adata_60, n_genes=30, sharey=False, fontsize = 10)

# get a tidy table of top markers
markers = sc.get.rank_genes_groups_df(adata_60, group=None)
markers.head()


# In[82]:


marker_genes = pd.DataFrame({
    group: adata.uns['rank_genes_groups']['names'][group][:50]
    for group in adata.uns['rank_genes_groups']['names'].dtype.names
})
marker_genes.head()
#marker_genes.to_csv('../data/broad_markers_leiden0-5.csv')


# In[83]:


combined_list = []
for col in marker_genes.columns: 
    print(col)
    genes = marker_genes[col].tolist()
    print(" ".join(genes))
    print(' ')
    combined_list.append(genes)
combined_list = [item for sublist in combined_list for item in sublist]
print()


# In[ ]:




