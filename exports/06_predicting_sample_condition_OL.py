#!/usr/bin/env python
# coding: utf-8

# In[1]:


import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import warnings
warnings.filterwarnings("ignore")


# In[15]:


adata = sc.read_h5ad('../data/mtDNA_DSB_5k_clustered_LLM_anno.h5ad')


# In[16]:


adata = adata[adata.obs.cell_class.str.contains('ligo')]


# In[79]:


adata.obs.cell_class.unique()


# In[17]:


# --- Find marker genes per cluster ---
sc.tl.rank_genes_groups(adata, groupby="sample_id", method="t-test")
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, fontsize = 15)

# get a tidy table of top markers
markers = sc.get.rank_genes_groups_df(adata, group=None)
markers.head()


# In[18]:


markers = sc.get.rank_genes_groups_df(adata, group=None)
markers.head()


# In[21]:


marker_genes = pd.DataFrame({
    group: adata.uns['rank_genes_groups']['names'][group][:10]
    for group in adata.uns['rank_genes_groups']['names'].dtype.names
})
marker_genes.head()
#marker_genes.to_csv('../data/broad_markers_leiden0-5.csv')


# In[22]:


combined_list = []
for col in marker_genes.columns: 
    print(col)
    genes = marker_genes[col].tolist()
    print(" ".join(genes))
    print(' ')
    combined_list.append(genes)
combined_list = [item for sublist in combined_list for item in sublist]
print()


# In[24]:


sc.pl.dotplot(
        adata,
        var_names=combined_list,
        groupby="sample_id",
        standard_scale="var",
        #dot_max=0.5,
        #dot_min=0.05,
        color_map="Reds",

        dendrogram=True,
        figsize=(20, 4)
    )


# In[26]:


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# pseudobulk per sample_id
pb = adata.to_df().groupby(adata.obs["sample_id"]).mean()

# keep only marker genes
all_markers =combined_list#[g for genes in marker_modules.values() for g in genes]
pb_markers = pb.loc[:, pb.columns.intersection(all_markers)]

# scale by gene
pb_scaled = (pb_markers - pb_markers.mean(axis=0)) / pb_markers.std(axis=0)

# clustered heatmap
sns.clustermap(
    pb_scaled.T,
    col_cluster=True, row_cluster=True,
    cmap="vlag", center=0,
    figsize=(12, 10)
)
plt.show()


# In[27]:


sample_annotations_balanced_OL = {
    "RB4282": "DSB",
    "RB4350": "DSB",
    "RB4403": "DSB",
    "RB4405": "DSB",
    "RB4630": "DSB",
    "RB4401": "DSB",        # swing sample, assigned to DSB for balance
    "RB4498": "Control",
    "RB4620": "Control",
    "RB4627": "Control",
    "RB4653": "Control",
    "RB4658": "Control",
    "RB4676": "Control"
}


# In[29]:


adata.obs["condition_predicted"] = adata.obs["sample_id"].map(sample_annotations_balanced_OL)


# In[32]:


sc.tl.rank_genes_groups(adata, groupby="condition_predicted", method="wilcoxon")
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)


# In[33]:


markers = sc.get.rank_genes_groups_df(adata, group=None)
markers.head()


# In[36]:


marker_genes = pd.DataFrame({
    group: adata.uns['rank_genes_groups']['names'][group][:15]
    for group in adata.uns['rank_genes_groups']['names'].dtype.names
})
marker_genes.head()
#marker_genes.to_csv('../data/broad_markers_leiden0-5.csv')


# In[37]:


combined_list = []
for col in marker_genes.columns: 
    print(col)
    genes = marker_genes[col].tolist()
    print(" ".join(genes))
    print(' ')
    combined_list.append(genes)
combined_list = [item for sublist in combined_list for item in sublist]
print()


# In[84]:


sc.pl.dotplot(
        adata,
        var_names=['Ldha','Serpina3n'],
        groupby="condition_predicted",
        standard_scale="var",
        #dot_max=0.5,
        #dot_min=0.05,
        color_map="Reds",

        dendrogram=True,
        figsize=(4, 2)
    )


# In[39]:


sc.pl.dotplot(
        adata,
        var_names=combined_list,
        groupby="condition_predicted",
        standard_scale="var",
        #dot_max=0.5,
        #dot_min=0.05,
        color_map="Reds",

        dendrogram=True,
        figsize=(20, 4)
    )


# In[40]:


sc.pl.dotplot(
        adata,
        var_names=combined_list,
        groupby="sample_id",
        standard_scale="var",
        #dot_max=0.5,
        #dot_min=0.05,
        color_map="Reds",

        dendrogram=True,
        figsize=(20, 4)
    )


# In[ ]:


sample_annotations_balanced_OL = {
    "RB4282": "DSB",
    "RB4350": "DSB",
    "RB4403": "DSB",
    "RB4405": "DSB",
    "RB4630": "DSB",
    "RB4401": "DSB",        # swing sample, assigned to DSB for balance
    "RB4498": "Control",
    "RB4620": "Control",
    "RB4627": "Control",
    "RB4653": "Control",
    "RB4658": "Control",
    "RB4676": "Control"
}


# In[31]:


sc.pl.pca(adata, color="condition_predicted")


# In[50]:


sc.pl.pca(adata, color=["Mbp"])


# In[44]:


import pandas as pd

# PC loadings (genes x PCs)
loadings = pd.DataFrame(
    adata.varm["PCs"], 
    index=adata.var_names,
    columns=[f"PC{i+1}" for i in range(adata.varm["PCs"].shape[1])]
)

# Top drivers of PC1
top_pc1 = loadings["PC1"].abs().sort_values(ascending=False).head(20)
print("Top PC1 drivers:", top_pc1.index.tolist())

# Top drivers of PC2
top_pc2 = loadings["PC2"].abs().sort_values(ascending=False).head(20)
print("Top PC2 drivers:", top_pc2.index.tolist())


# In[68]:


import scanpy as sc
import numpy as np
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt

# choose the embedding you want to split on
XY = adata.obsm['X_umap']          # or: adata.obsm['X_pca'][:, :2]

# k-means into 2 blobs
labels = KMeans(n_clusters=2, random_state=42).fit_predict(XY)
adata.obs['blob2'] = labels.astype(str)

# quick check
sc.pl.pca(adata, color=["blob2"])
sc.pl.embedding(adata, basis='umap', color='blob2', frameon=False)



# In[69]:


sc.pl.dotplot(adata, ["Mbp","Cnp","Mog","Cldn11","Mal","Ermn",
                      "Slc1a3","Gja1","Aldoc","Mt2"], groupby="blob2")


# In[76]:


# --- Find marker genes per cluster ---
sc.tl.rank_genes_groups(adata, groupby="blob2", method="t-test")
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, fontsize = 15)

# get a tidy table of top markers
markers = sc.get.rank_genes_groups_df(adata, group=None)
markers.head()


# In[77]:


marker_genes = pd.DataFrame({
    group: adata.uns['rank_genes_groups']['names'][group][:20]
    for group in adata.uns['rank_genes_groups']['names'].dtype.names
})
marker_genes.head()
#marker_genes.to_csv('../data/broad_markers_leiden0-5.csv')


# In[78]:


combined_list = []
for col in marker_genes.columns: 
    print(col)
    genes = marker_genes[col].tolist()
    print(" ".join(genes))
    print(' ')
    combined_list.append(genes)
combined_list = [item for sublist in combined_list for item in sublist]
print()


# In[85]:


adata


# In[86]:


sc.pl.dotplot(
        adata,
        var_names='Ldha',
        groupby="sample_id",
        standard_scale="var",
        #dot_max=0.5,
        #dot_min=0.05,
        color_map="Reds",

        dendrogram=True,
        figsize=(7, 4)
    )


# In[87]:


adata


# In[88]:


sc.pl.dotplot(
        adata,
        var_names='Ldha',
        groupby="condition_predicted",
        standard_scale="var",
        #dot_max=0.5,
        #dot_min=0.05,
        color_map="Reds",

        dendrogram=True,
        figsize=(7, 4)
    )


# In[ ]:


['RB4676', 'RB4653', 'RB4627', 'RB4620', 'RB4658', 'RB4350', 'RB4630', 'RB4282', 'RB4405', 'RB4401', 'RB4498', 'RB4403']
['Control', 'Control', 'Control', 'Control', 'Control', 'mtDSB', 'Control', 'mtDSB', 'mtDSB', 'mtDSB', 'mtDSB', 'mtDSB']

