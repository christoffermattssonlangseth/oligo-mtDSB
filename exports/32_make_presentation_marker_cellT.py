#!/usr/bin/env python
# coding: utf-8

# In[2]:


panel.iloc[:,0]


# In[18]:


ct_counts


# In[20]:


genes


# In[22]:


hits


# In[24]:


hits


# In[28]:


ct_counts


# In[37]:


from wordcloud import WordCloud
import matplotlib.pyplot as plt

# Assume you already have ct_counts (a pandas Series)
# e.g. ct_counts = ad.obs["cell_types"].value_counts()

# Limit to top 50
freqs = ct_counts.to_dict()

# Create word cloud
wc = WordCloud(
    width=1200,
    height=600,
    background_color="white",
    colormap="Set2",      # 🔥 try: "coolwarm", "plasma", "tab20", "cividis", "magma", "cubehelix", "Set2", "Pastel1"
    prefer_horizontal=1.0,
).generate_from_frequencies(freqs)

plt.figure(figsize=(20,20))
plt.imshow(wc, interpolation="bilinear")
plt.axis("off")
plt.title("Cell type abundance (font size ∝ abundance)", fontsize=14)
plt.show()


# In[23]:


import pandas as pd
import numpy as np

# 1) Your panel
panel = pd.read_csv("markers_mouse_xenium.csv")
genes = panel.iloc[:,1].astype(str).str.upper().unique()

# 2) PanglaoDB (downloaded CSV)
pang = pd.read_csv("/Users/christoffer/Downloads/PanglaoDB_markers_27_Mar_2020.tsv", sep = '\t')  # path to your local file
# Harmonize columns (check exact headers in the CSV you downloaded)
cols = [c.lower() for c in pang.columns]
pang.columns = cols
# Typical columns: 'species', 'official gene symbol', 'cell type', 'tissue'
pang['official gene symbol'] = pang['official gene symbol'].str.upper()

# Filter to mouse & nervous system tissues (broader than just CNS)
pang_mouse = pang[pang['species'].str.contains("Mm", case=False, na=False)]
pang_brainish = pang_mouse[pang_mouse['organ'].str.contains("brain|cns|cortex|hippocampus|midbrain|cerebellum|spinal", case=False, na=True)]

# Intersect
hits = pang_mouse[pang_mouse['official gene symbol'].isin(genes)].copy()

# Summaries
ct_counts = hits['cell type'].value_counts().sort_values(ascending=False)
broad_map = {
    "neuron":"Neuronal","astrocyte":"Glial","oligodendrocyte":"Glial","opc":"Glial",
    "microglia":"Immune","t cell":"Immune","b cell":"Immune","macrophage":"Immune",
    "endothelial cell":"Vascular","pericyte":"Vascular","smooth muscle cell":"Vascular",
    "fibroblast":"Stromal","ependymal cell":"Glial"
}
hits['broad'] = hits['cell type'].str.lower().map(broad_map).fillna("Other")
broad_counts = hits['broad'].value_counts().sort_values(ascending=False)

ct_counts.head(20), broad_counts


# In[ ]:




