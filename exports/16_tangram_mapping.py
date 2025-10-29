#!/usr/bin/env python
# coding: utf-8

# In[1]:


import scanpy as sc
import tangram as tg


# In[2]:


ad_sp = sc.read_h5ad('../data/mtDNA_DSB_5k_clustered_LLM_anno.h5ad')
ad_sc = sc.read_h5ad('/Users/christoffer/work/karolinska/development/data/abc_atlas/combined_scRNAseq_for_tangram_subset_genes.h5ad')
tg.pp_adatas(ad_sc, ad_sp, genes=None)


# In[3]:


ad_sp.X.max()


# In[ ]:


ad_map = tg.map_cells_to_space(
                   ad_sc, 
                   ad_sp,         
                   mode='clusters',
                   cluster_label='class')


# In[ ]:




