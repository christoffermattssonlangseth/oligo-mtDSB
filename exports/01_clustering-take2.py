#!/usr/bin/env python
# coding: utf-8

# In[1]:


import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import warnings
warnings.filterwarnings("ignore")


# In[2]:


ad = sc.read('../data/mtDNA_DSB_5k_clustered.h5ad')


# In[3]:


ad


# In[4]:


resolutions = [0.5, 1,1.5, 2, 2.1,2.3,2.5,3]

for resolution in resolutions:
    key = f'leiden_{resolution}'

    if key in ad.obs.columns:
        print(f"Skipping {resolution}: {key} already exists.")
    else:
        print(f"Clustering at resolution {resolution}...")
        sc.tl.leiden(ad, resolution=resolution, key_added=key)
        print("Done.")

    # plot UMAP
    sc.pl.umap(ad, color=key, legend_loc='on data', frameon=False)


# In[6]:


ad.write('../data/mtDNA_DSB_5k_clustered_higher_res.h5ad')


# In[ ]:




