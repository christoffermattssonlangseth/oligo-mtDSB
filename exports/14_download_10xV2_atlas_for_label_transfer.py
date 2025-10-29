#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
from pathlib import Path
import numpy as np
import anndata
import time

from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache


# In[2]:


download_base = Path('../../data/abc_atlas')
abc_cache = AbcProjectCache.from_cache_dir(download_base)

abc_cache.current_manifest


# In[ ]:


cell = abc_cache.get_metadata_dataframe(directory='WMB-10Xv2', file_name='cell_metadata')
cell.set_index('cell_label', inplace=True)


# In[ ]:


matrices = cell.groupby(['dataset_label', 'feature_matrix_label'])[['library_label']].count()
matrices.columns  = ['cell_count']
matrices


# In[4]:


files = abc_cache.list_data_files('WMB-10Xv2')
print(files)


# In[6]:


files = [
 'WMB-10Xv2-OLF/log2',
 'WMB-10Xv2-OLF/raw',
 'WMB-10Xv2-TH/log2',
 'WMB-10Xv2-TH/raw']


# In[7]:


for file in files: 
    if 'raw' in file:
        print('get this file')
        print(file)
        abc_cache.get_data_path(directory='WMB-10Xv2', file_name=file)


# In[ ]:


cell = abc_cache.get_metadata_dataframe(directory='WMB-10Xv2', file_name='cell_metadata')
cell.set_index('cell_label', inplace=True)


# In[ ]:


ad = anndata.read_h5ad(file,backed='r')
gene = ad.var

