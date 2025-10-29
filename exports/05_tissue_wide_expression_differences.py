#!/usr/bin/env python
# coding: utf-8

# In[34]:


import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import warnings
warnings.filterwarnings("ignore")


# In[35]:


adata = sc.read_h5ad('../data/mtDNA_DSB_5k_clustered_LLM_anno.h5ad')


# In[36]:


# --- Find marker genes per cluster ---
sc.tl.rank_genes_groups(adata, groupby="sample_id", method="t-test")
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, fontsize = 15)

# get a tidy table of top markers
markers = sc.get.rank_genes_groups_df(adata, group=None)
markers.head()


# In[37]:


markers = sc.get.rank_genes_groups_df(adata, group=None)
markers.head()


# In[38]:


marker_genes = pd.DataFrame({
    group: adata.uns['rank_genes_groups']['names'][group][:10]
    for group in adata.uns['rank_genes_groups']['names'].dtype.names
})
marker_genes.head()
#marker_genes.to_csv('../data/broad_markers_leiden0-5.csv')


# In[39]:


combined_list = []
for col in marker_genes.columns: 
    print(col)
    genes = marker_genes[col].tolist()
    print(" ".join(genes))
    print(' ')
    combined_list.append(genes)
combined_list = [item for sublist in combined_list for item in sublist]
print()


# In[8]:


sc.pl.dotplot(
        adata,
        var_names=combined_list,
        groupby="sample_id",
        standard_scale="var",
        #dot_max=0.5,
        #dot_min=0.05,
        color_map="Reds",

        dendrogram=True,
        figsize=(15, 4)
    )


# In[9]:


marker_modules = {
    "Stress_ImmediateEarly": [
        "Sgk1", "Ddit4", "Tsc22d3", "Nr4a1", "Egr1", "Fos", "Arc"
    ],
    "Mitochondrial_UPR": [
        "Hspa5", "Plin4", "Pdia3", "Pink1"
    ],
    "Cytoskeletal_Myelin": [
        "Kif5a", "Kif5b", "Kif5c", "Map2", "Tppp", "Stmn4"
    ],
    "Metabolic_Remodeling": [
        "Ldha", "Eno2", "Dbp", "Nr1d1"
    ],
    "Glial_Immune": [
        "Gfap", "Serpina3n", "Il33"
    ],
    "Neuronal_Synaptic": [
        "Gnb1", "Cyfip2", "Syn1", "Gpm6a", "Rtn3", "App", "Aplp1", "Adgrg1", "Sez6l2"
    ],
    "Hormone_Endocrine": [
        "Gh", "Prl", "Ppp1r1b", "Gpr88", "Kdr"
    ],
    "Sex_Linked": [
        "Xist"
    ]
}


# In[10]:


sc.pl.dotplot(
        adata,
        var_names=marker_modules,
        groupby="sample_id",
        standard_scale="var",
        #dot_max=0.5,
        #dot_min=0.05,
        color_map="Reds",

        dendrogram=True,
        figsize=(15, 4)
    )


# In[11]:


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# pseudobulk per sample_id
pb = adata.to_df().groupby(adata.obs["sample_id"]).mean()

# keep only marker genes
all_markers =[g for genes in marker_modules.values() for g in genes]
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


# In[12]:


sample_annotations = {
    "RB4282": "DSB model",
    "RB4350": "DSB model",
    "RB4401": "Control",
    "RB4403": "DSB model",
    "RB4405": "DSB model",
    "RB4498": "Uncertain (maybe DSB)",
    "RB4620": "Control/uncertain",
    "RB4627": "Control/uncertain",
    "RB4630": "Uncertain",
    "RB4653": "Control",
    "RB4658": "Control",
    "RB4676": "Control",
}


# In[13]:


sample_annotations_balanced = {
    "RB4282": "DSB",
    "RB4403": "DSB",
    "RB4405": "DSB",
    "RB4350": "DSB",
    "RB4498": "DSB",
    "RB4630": "DSB",
    "RB4620": "Control",
    "RB4627": "Control",
    "RB4676": "Control",
    "RB4653": "Control",
    "RB4658": "Control",
    "RB4401": "Control"
}


# In[14]:


adata.obs["condition_predicted"] = adata.obs["sample_id"].map(sample_annotations_balanced)


# In[16]:


sc.pl.pca(adata, color="condition_predicted")


# In[18]:


sc.tl.rank_genes_groups(adata, groupby="condition_predicted", method="wilcoxon")
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)


# In[68]:


import pandas as pd

def get_DEGs(adata, n_genes=100):
    res = adata.uns['rank_genes_groups']
    groups = res['names'].dtype.names
    dfs = []
    for g in groups:
        df = pd.DataFrame({
            "gene": res['names'][g],
            "logfoldchange": res['logfoldchanges'][g],
            "pval_adj": res['pvals_adj'][g],
            "group": g
        })
        dfs.append(df.head(n_genes))
    return pd.concat(dfs)

deg_table = get_DEGs(adata, n_genes=10)


# In[ ]:





# In[69]:


deg_table


# In[70]:


import numpy as np

# compute variance across runs for each gene
mean_per_run = adata.to_df().groupby(adata.obs["sample_id"]).mean()
var_across_runs = mean_per_run.var(axis=0)

# top varying genes
top_var_genes = var_across_runs.sort_values(ascending=False).head(50)


# In[71]:


top_var_genes


# In[49]:


import re

# extract "RB####" from the run string and put into new column
adata.obs["sample_id"] = adata.obs["run"].str.extract(r"(RB\d+)")


# In[51]:


sc.pl.umap(adata, color="sample_id")


# In[54]:


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

# Average expression per sample_id
pb = adata.to_df().groupby(adata.obs["sample_id"]).mean()

# Run PCA
pca = PCA(n_components=2)
coords = pca.fit_transform(pb)

df_pca = pd.DataFrame(coords, index=pb.index, columns=["PC1","PC2"])

# Plot
plt.figure(figsize=(6,6))
sns.scatterplot(x="PC1", y="PC2", data=df_pca)
for i, txt in enumerate(df_pca.index):
    plt.annotate(txt, (df_pca.PC1[i], df_pca.PC2[i]))
plt.title("Pseudobulk PCA by sample_id")
plt.show()


# In[56]:


import pandas as pd
loadings = pd.Series(pca.components_[0], index=pb.columns)  # PC1 loadings
print(loadings.sort_values(ascending=False).head(20))


# In[78]:


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# --- Get top PC1 loadings (absolute value to get both positive & negative drivers) ---
pc1_loadings = pd.Series(pca.components_[0], index=pb.columns)
top_genes = pc1_loadings.abs().sort_values(ascending=False).head(30).index

# --- Subset pseudobulk matrix ---
pb_pc1 = pb[top_genes]

# --- Z-score per gene ---
pb_pc1_scaled = (pb_pc1 - pb_pc1.mean(axis=0)) / pb_pc1.std(axis=0)

# --- Clustered heatmap ---
sns.clustermap(
    pb_pc1_scaled.T,
    col_cluster=True, row_cluster=True,
    cmap="vlag", center=0,
    figsize=(12, 10),
    xticklabels=True, yticklabels=True
)
plt.show()


# In[77]:


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# --- Step 1. Pseudobulk: average per sample_id ---
pb = adata.to_df().groupby(adata.obs["sample_id"]).mean()

# --- Step 2. Select top variable genes ---
top_var = pb.var(axis=0).sort_values(ascending=False).head(30).index
pb_top = pb[top_var]

# --- Step 3. Scale by gene (z-score) ---
pb_scaled = (pb_top - pb_top.mean(axis=0)) / pb_top.std(axis=0)

# --- Step 4. Clustered heatmap ---
sns.clustermap(
    pb_scaled.T,  # genes as rows
    col_cluster=True, row_cluster=True,
    cmap="vlag", center=0,
    figsize=(12, 10),
    xticklabels=True, yticklabels=True
)
plt.show()


# In[22]:


import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

# ---------- 0) inputs ----------
# assumes you already have: adata.obs["sample_id"]
# (optionally subset to OLs beforehand)

# define marker modules (edit if you like)
stress_genes = [
    "Sgk1","Ddit4","Tsc22d3","Pink1","Mt2","Hspa1a","Hspa1b","Hsph1","Dnajb1","Atf4","Atf5"
]
myelin_genes = [
    "Mbp","Mog","Mag","Cnp","Cldn11","Mal","Ermn","Pllp","Aspa","Ugt8a","Ptgds","Car2","Sept4"
]

# ---------- 1) pseudobulk by sample ----------
# average expression per sample (use mean; you can use sum for counts if desired)
pb = adata.to_df().groupby(adata.obs["sample_id"]).mean()
# store list of sample_ids
sample_ids = pb.index.tolist()

# ---------- 2) build a samples×genes AnnData ----------
adata_pb = sc.AnnData(pb.values)
adata_pb.obs.index = sample_ids
adata_pb.var.index = pb.columns

# ---------- 3) dimensionality reduction & pseudotime ----------
# scale (z-score genes) for stability at sample level
sc.pp.scale(adata_pb, max_value=10)
sc.tl.pca(adata_pb, n_comps=10)
sc.pp.neighbors(adata_pb, n_neighbors=min(6, adata_pb.n_obs-1), n_pcs=10)  # small graph; n_obs is small (12)
sc.tl.diffmap(adata_pb)
# pick a root sample automatically: lowest stress / highest myelin, or just the first
# we'll choose the sample with minimal stress-minus-myelin score as root (control-like)
# compute quick module summaries on *unscaled* pb
def safe_mean(mat, genes):
    g = [g for g in genes if g in mat.columns]
    return mat[g].mean(axis=1) if len(g) else pd.Series(0, index=mat.index)

stress_score_pb = safe_mean(pb, stress_genes)
myelin_score_pb = safe_mean(pb, myelin_genes)
# after neighbors + diffmap on adata_pb
# choose a root sample (e.g., lowest stress − myelin)
root_sample = (stress_score_pb - myelin_score_pb).idxmin()
root_index  = list(adata_pb.obs_names).index(root_sample)

# set the root index for DPT in .uns['iroot']
adata_pb.uns['iroot'] = root_index

# now run DPT (no root args here)
sc.tl.dpt(adata_pb, n_dcs=5)

# pseudotime is stored here:
pt = adata_pb.obs['dpt_pseudotime']
print("Root:", root_sample)
print("Ordered samples:", list(pt.sort_values().index))
# ---------- 4) plots: PCA & Diffusion colored by pseudotime ----------
sc.pl.pca(adata_pb, color="dpt_pseudotime", title="Sample PCA colored by pseudotime")



# In[23]:


# diffusion components scatter
import matplotlib.pyplot as plt
X = adata_pb.obsm["X_diffmap"][:, :2]
plt.figure(figsize=(5,5))
sc = plt.scatter(X[:,0], X[:,1], c=adata_pb.obs["dpt_pseudotime"], s=80)
for i, sid in enumerate(adata_pb.obs_names):
    plt.text(X[i,0]+0.01, X[i,1]+0.01, sid, fontsize=9)
plt.xlabel("DC1"); plt.ylabel("DC2"); plt.title("Samples (Diffusion map)")
plt.colorbar(sc, label="pseudotime")
plt.tight_layout(); plt.show()

# ---------- 5) module trends along pseudotime ----------
# align to pseudotime order
order = adata_pb.obs["dpt_pseudotime"].sort_values().index
stress_aligned = stress_score_pb.loc[order]
myelin_aligned = myelin_score_pb.loc[order]
ptime = adata_pb.obs.loc[order, "dpt_pseudotime"]

plt.figure(figsize=(6,4))
plt.plot(ptime.values, (stress_aligned - stress_aligned.mean())/stress_aligned.std(), marker='o', label="Stress module (z)")
plt.plot(ptime.values, (myelin_aligned - myelin_aligned.mean())/myelin_aligned.std(), marker='o', label="Myelin module (z)")
plt.xlabel("Pseudotime (DPT)"); plt.ylabel("Z-scored module")
plt.title("Module trajectories across inferred disease time")
plt.legend(); plt.tight_layout(); plt.show()

# ---------- 6) handy outputs ----------
print("\nRoot sample chosen (control-like):", root_sample)
print("\nSamples ordered by pseudotime (early → late):")
print(list(order))

# also store results for downstream
adata_pb.obs["stress_score"] = stress_score_pb
adata_pb.obs["myelin_score"] = myelin_score_pb
# save if you want
# adata_pb.write_h5ad("samples_pseudotime.h5ad", compression="gzip")


# In[24]:


# correlations between pseudotime and modules
pt = adata_pb.obs['dpt_pseudotime']
from scipy.stats import spearmanr, kendalltau

for name, vec in {"Stress": adata_pb.obs["stress_score"],
                  "Myelin": adata_pb.obs["myelin_score"]}.items():
    print(name, "Spearman r, p=", spearmanr(pt, vec))
    print(name, "Kendall tau, p=", kendalltau(pt, vec))


# In[25]:


# top loadings on DC1 (or PC1 if you used PCA)
dc1 = pd.Series(adata_pb.obsm['X_diffmap'][:,0], index=adata_pb.obs_names)
# already computed stress/myelin lists earlier; you can also get per-gene correlations:
mean_expr = pb  # samples x genes (unscaled)
cors = mean_expr.apply(lambda g: spearmanr(pt, g)[0], axis=0).sort_values()
print("Most decreasing along time:", cors.head(15).index.tolist())
print("Most increasing along time:", cors.tail(15).index.tolist())


# In[26]:


order = pt.sort_values().index
stages = {sid: ("early" if i < 4 else "mid" if i < 8 else "late")
          for i, sid in enumerate(order)}
adata_pb.obs["stage_blind"] = adata_pb.obs_names.map(stages)


# In[27]:


adata_pb.obs["stage_blind"]


# In[ ]:


sample_annotations_balanced = {
    "RB4282": "DSB",
    "RB4403": "DSB",
    "RB4405": "DSB",
    "RB4350": "DSB",
    "RB4498": "DSB",
    "RB4630": "DSB",
    "RB4620": "Control",
    "RB4627": "Control",
    "RB4676": "Control",
    "RB4653": "Control",
    "RB4658": "Control",
    "RB4401": "Control"
}


# In[29]:


import pandas as pd
import scanpy as sc

# pseudobulk (samples x genes)
pb = adata.to_df().groupby(adata.obs["sample_id"]).mean()

stress = ["Sgk1","Ddit4","Tsc22d3","Pink1","Mt2","Hspa1a","Hspa1b","Hsph1","Dnajb1","Atf4","Atf5", 'Serpina3n','Ldha']
myelin = ["Mbp","Mog","Mag","Cnp","Cldn11","Mal","Ermn","Pllp","Aspa","Ugt8a","Ptgds","Car2","Sept4"]

def score(mat, genes):
    g = [x for x in genes if x in mat.columns]
    return mat[g].mean(axis=1) if g else pd.Series(0, index=mat.index)

s_stress = score(pb, stress)
s_myelin = score(pb, myelin)

anchor_control = (s_stress - s_myelin).idxmin()
print("Anchor control:", anchor_control)


# In[ ]:


import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

# ---------- 0) inputs ----------
# assumes you already have: adata.obs["sample_id"]
# (optionally subset to OLs beforehand)

# define marker modules (edit if you like)
stress_genes = [
    "Sgk1","Ddit4","Tsc22d3","Pink1","Mt2","Hspa1a","Hspa1b","Hsph1","Dnajb1","Atf4","Atf5"
]
myelin_genes = [
    "Mbp","Mog","Mag","Cnp","Cldn11","Mal","Ermn","Pllp","Aspa","Ugt8a","Ptgds","Car2","Sept4"
]

# ---------- 1) pseudobulk by sample ----------
# average expression per sample (use mean; you can use sum for counts if desired)
pb = adata.to_df().groupby(adata.obs["sample_id"]).mean()
# store list of sample_ids
sample_ids = pb.index.tolist()

# ---------- 2) build a samples×genes AnnData ----------
adata_pb = sc.AnnData(pb.values)
adata_pb.obs.index = sample_ids
adata_pb.var.index = pb.columns

# ---------- 3) dimensionality reduction & pseudotime ----------
# scale (z-score genes) for stability at sample level
sc.pp.scale(adata_pb, max_value=10)
sc.tl.pca(adata_pb, n_comps=10)
sc.pp.neighbors(adata_pb, n_neighbors=min(6, adata_pb.n_obs-1), n_pcs=10)  # small graph; n_obs is small (12)
sc.tl.diffmap(adata_pb)
# pick a root sample automatically: lowest stress / highest myelin, or just the first
# we'll choose the sample with minimal stress-minus-myelin score as root (control-like)
# compute quick module summaries on *unscaled* pb
def safe_mean(mat, genes):
    g = [g for g in genes if g in mat.columns]
    return mat[g].mean(axis=1) if len(g) else pd.Series(0, index=mat.index)

stress_score_pb = safe_mean(pb, stress_genes)
myelin_score_pb = safe_mean(pb, myelin_genes)
# after neighbors + diffmap on adata_pb
# choose a root sample (e.g., lowest stress − myelin)
root_sample = (stress_score_pb - myelin_score_pb).idxmin()
root_index  = list(adata_pb.obs_names).index(root_sample)

# set the root index for DPT in .uns['iroot']
adata_pb.uns['iroot'] = root_index

# now run DPT (no root args here)
sc.tl.dpt(adata_pb, n_dcs=5)

# pseudotime is stored here:
pt = adata_pb.obs['dpt_pseudotime']
print("Root:", root_sample)
print("Ordered samples:", list(pt.sort_values().index))
# ---------- 4) plots: PCA & Diffusion colored by pseudotime ----------
sc.pl.pca(adata_pb, color="dpt_pseudotime", title="Sample PCA colored by pseudotime")


