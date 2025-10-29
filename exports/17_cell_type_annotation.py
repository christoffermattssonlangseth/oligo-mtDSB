#!/usr/bin/env python
# coding: utf-8

# In[6]:


import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import warnings
warnings.filterwarnings("ignore")


# In[10]:


ad = sc.read_h5ad('../data/mtDNA_DSB_5k_clustered_LLM_anno.h5ad')


# In[13]:


ad.X.max()


# In[12]:


# --- Find marker genes per cluster ---
sc.tl.rank_genes_groups(ad, groupby="leiden_2", method="t-test")
sc.pl.rank_genes_groups(ad, n_genes=25, sharey=False)

# get a tidy table of top markers
markers = sc.get.rank_genes_groups_df(ad, group=None)
markers.head()


# In[15]:


marker_genes = pd.DataFrame({
    group: ad.uns['rank_genes_groups']['names'][group][:20]
    for group in ad.uns['rank_genes_groups']['names'].dtype.names
})
for col in marker_genes.columns: 
    print(col)
    genes = marker_genes[col].tolist()
    print(" ".join(genes))
    print(' ')


# In[20]:


annotation = {
    '0': 'Telencephalon astrocytes I',
    '1': 'Pericytes I',
    '2': 'Mature oligodendrocytes I',
    '3': 'Striatal neurons I',
    '4': 'Neurons I',
    '5': 'Olfactory astrocytes I',
    '6': 'Excitatory neurons I',
    '7': 'Excitatory neurons II',
    '8': 'Excitatory neurons (cortex) I',
    '9': 'Inhibitory neurons I',
    '10': 'Mature oligodendrocytes II',
    '11': 'Vascular leptomeningeal cells I',
    '12': 'Inhibitory neurons II',
    '13': 'Microglia I',
    '14': 'Mature oligodendrocytes III',
    '15': 'Excitatory neurons III',
    '16': 'Vascular endothelial cells I',
    '17': 'Neurons II',
    '18': 'Excitatory neurons (thalamus) I',
    '19': 'Excitatory neurons (cortex) II',
    '20': 'Oligodendrocytes precursor cells I',
    '21': 'Excitatory neurons (thalamus) II',
    '22': 'Telencephalon astrocytes II',
    '23': 'Mature oligodendrocytes IV',
    '24': 'Excitatory neurons IV',
    '25': 'Excitatory neurons V',
    '26': 'Mature oligodendrocytes V',
    '27': 'Chorid plexus epithelial cells I',
    '28': 'Ependymal cells I',
    '29': 'Excitatory neurons VI',
    '30': 'Cholinergic neurons I',
    '31': 'Excitatory neurons VII',
    '32': 'Excitatory neurons (cortex) III',
    '33': 'Mature oligodendrocytes VI',
    '34': 'Olfactory astrocytes II',
    '35': 'Vascular endothelial cells II',
    '36': 'D1 medium spiny neurons (striatum) I',
    '37': 'Mature oligodendrocytes VII',
    '38': 'Mature oligodendrocytes VIII',
    '39': 'Endothelial cells I',
    '40': 'Unknown',
    '41': 'Inhibitory neurons III',
    '42': 'Neural progenitors I',
    '43': 'Mature oligodendrocytes (mic) I',
    '44': 'Unknown',
    '45': 'Unknown'
}
# Strip roman numerals automatically
import re
annotation_base = {k: re.sub(r"\s+[IVXLCDM]+$", "", v) for k, v in annotation.items()}


# In[22]:


ad.obs['cell_type'] = ad.obs.leiden_2.map(annotation)
ad.obs['cell_class'] = ad.obs.leiden_2.map(annotation_base)


# In[32]:


import matplotlib.pyplot as plt

with plt.rc_context({"figure.figsize": (20, 15)}):
    sc.pl.umap(
        ad,
        color="cell_type",
        legend_loc="on data",
        frameon=False,
        legend_fontoutline = 5,
        size=10
    )


# In[33]:


import matplotlib.pyplot as plt

with plt.rc_context({"figure.figsize": (20, 15)}):
    sc.pl.umap(
        ad,
        color="cell_class",
        legend_loc="on data",
        frameon=False,
        legend_fontoutline = 5,
        size=10
    )


# In[38]:


for run in ad.obs['run'].unique():
    print(run)
    ad_int = ad[ad.obs['run'] == run]
    with plt.rc_context({'figure.figsize': (20, 10)}):
        sc.pl.spatial(ad_int, spot_size=15, color = 'cell_class')
    plt.show()



# In[44]:


ad.write('../data/mtDNA_DSB_5k_clustered_manual_annotation.h5ad')


# In[ ]:




