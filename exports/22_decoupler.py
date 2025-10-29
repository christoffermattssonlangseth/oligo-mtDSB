#!/usr/bin/env python
# coding: utf-8

# # Pseudobulk Enrichment Analysis
# When cell identity clusters are well-defined, it can be advantageous to perform analyses at the pseudobulk level rather than at the single-cell level. Pseudobulking involves aggregating counts across cells of the same type within each sample, effectively creating sample-level gene expression profiles per cell type. This approach helps mitigate the effects of technical noise and dropouts common in single-cell data, enabling the detection of lowly expressed genes that might otherwise be missed.
# 
# Moreover, conducting differential expression analysis (DEA) at the pseudobulk level, treating each biological sample as the unit of observation, is statistically more robust. Unlike single-cell DEA, which assumes cells are independent (an assumption that is violated when cells originate from the same individual), sample-level pseudobulk analysis avoids inflation of p-values by reducing the number of observations and by correctly modeling biological replication {cite:p}psbulk.
# 
# The resulting gene-level statistics from pseudobulk DEA can then be used as input for downstream enrichment analyses.
# 
# In this notebook, we demonstrate how to use decoupler to infer transcription factor (TF) and pathway enrichment scores from a multi-sample scRNA-seq human dataset.
# 
# 

# In[1]:


import scanpy as sc
import decoupler as dc

sc.set_figure_params(figsize=(3, 3), frameon=False)


# In[ ]:


adata = sc.read_h5ad('../data/mtDNA_DSB_5k_clustered_manual_annotation_with_raw.h5ad')


# In[33]:


adata.X = adata.layers['counts']


# In[34]:


adata.X = adata.X.astype(int)


# In[14]:


adata.obs_names_make_unique()


# In[99]:


sc.pl.umap(adata, color=["cell_class", "condition"], ncols=1)


# ## Pseudobulking
# The pseudo-bulk approach involves the following steps:
# 
# Subsetting the cell type of interest
# Extracting their raw integer counts
# Summing their counts per gene into a single profile if they pass quality control
# Then, DEA can be performed if there are at least two biological replicates per condition (more replicates are recommended).
# 
# Pseudobulking can easily be performed using the function {func}decoupler.pp.pseudobulk. In this example, the counts are just summed, though other modes such as the mean or any custom aggregation function are available. For more information, refer to the mode argument.

# In[35]:


pdata = dc.pp.pseudobulk(
    adata=adata,
    sample_col="sample_id",
    groups_col="cell_class",
    mode="sum",
)


# ## Variability Exploration
# With pseudobulk profiles generated for each cell type and sample, variability across them can now be explored.
# 
# This involves some basic preprocessing followed by principal component analysis (PCA).

# In[36]:


# Store raw counts in layers
pdata.layers["counts"] = pdata.X.copy()

# Normalize, scale and compute pca
sc.pp.normalize_total(pdata, target_sum=1e4)
sc.pp.log1p(pdata)
sc.pp.scale(pdata, max_value=10)
sc.tl.pca(pdata)

# Return raw counts to X
dc.pp.swap_layer(adata=pdata, key="counts", inplace=True)


# In[37]:


dc.tl.rankby_obsm(pdata, key="X_pca")


# In[38]:


pdata.obs


# In[39]:


sc.pl.pca_variance_ratio(pdata)
dc.pl.obsm(adata=pdata, return_fig=True, nvar=5, titles=["PC scores", "Adjusted p-values"], figsize=(10, 5))


# In[40]:


pdata.obs


# In[41]:


sc.pl.pca(
    pdata,
    color=[ "condition", 'age', 'sex'],
    ncols=1,
    size=300,
    frameon=True,
)


# In[42]:


OL = pdata[pdata.obs['cell_class'] == 'Mature oligodendrocytes'].copy()


# In[43]:


dc.pl.filter_by_expr(
    adata=OL,
    group="condition",
    min_count=10,
    min_total_count=15,
    large_n=10,
    min_prop=0.7,
)
dc.pl.filter_by_prop(
    adata=OL,
    min_prop=0.1,
    min_smpls=2,
)


# In[44]:


# Import DESeq2
from pydeseq2.dds import DeseqDataSet, DefaultInference
from pydeseq2.ds import DeseqStats

# Build DESeq2 object
inference = DefaultInference(n_cpus=8)
dds = DeseqDataSet(
    adata=OL,
    design_factors=["condition"],
    refit_cooks=True,
    inference=inference,
)

# Compute LFCs
dds.deseq2()

# Extract contrast between conditions
stat_res = DeseqStats(dds, contrast=["condition", "mtDSB", "control"], inference=inference)

# Compute Wald test
stat_res.summary()


# In[45]:


# Extract results
results_df = stat_res.results_df
results_df


# In[3]:


# Build DESeq2 dataset
inference = DefaultInference(n_cpus=8)
dds = DeseqDataSet(
    adata=OL,
    metadata="condition",               # column in adata.obs
    ref_level=["condition", "control"], # reference level
    refit_cooks=True,
    inference=inference,
)


# In[52]:


import matplotlib.pyplot as plt


# In[63]:


import matplotlib.pyplot as plt
import decoupler as dc

import matplotlib.pyplot as plt
plt.rcParams["figure.figsize"] = (10, 6)    
dc.pl.volcano(
        results_df,
        x="log2FoldChange",
        y="pvalue",
        top=20,
    figsize=(12, 8)
    )


# In[101]:


data = results_df[["stat"]].T.rename(index={"stat": "disease.vs.normal"})
data


# ## Enrichment analysis
# Enrichment analysis tests whether a specific set of omics features is "overrepresented" or "coordinated" in the measured data compared to a background distribution. These sets are predefined based on existing biological knowledge and may vary depending on the omics technology used.
# 
# Enrichment analysis requires the use of an enrichment method, and several options are available. In the original manuscript of decoupler {cite:p}decoupler, we benchmarked multiple methods and found that the univariate linear model (ulm) outperformed the others. Therefore, we will use ulm in this vignette.
# 
# The scores from {func}decoupler.mt.ulm should be interpreted such that larger magnitudes indicate greater significance, while the sign reflects whether the features in the set are overrepresented (positive) or underrepresented (negative) compared to the background.
# 
# Transcription factor scoring from gene regulatory networks
# Transcription factors (TFs) are genes that, once translated into proteins, bind to DNA and regulate the expression of other genes by either promoting or inhibiting their transcription. Gene Regulatory Networks (GRNs) capture these TF-gene interactions and can be constructed from prior knowledge or inferred from omics data. The fundamental unit of a GRN is a TF and its associated target genes, collectively known as a regulon. Each regulon functions as a gene set in enrichment analysis.
# 
# Although TFs are measured in transcriptomic data, their transcript levels often do not reflect their actual activity in a given cell. Instead, scoring TFs through enrichment analysis based on the expression of their target genes provides a more accurate representation of their regulatory activity {cite:p}grn_review.
# 
# CollecTRI network
# CollecTRI is a comprehensive resource containing a curated collection of TFs and their transcriptional targets compiled from 12 different resources {cite:p}collectri. This collection provides an increased coverage of transcription factors and a superior performance in identifying perturbed TFs compared to other literature based GRNs such as DoRothEA {cite:p}dorothea. Similar to DoRothEA, interactions are weighted by their mode of regulation (activation or inhibition).
# 
# In this tutorial we will use the human version but other organisms are available. We can use decoupler to retrieve it from the OmniPath server {cite:p}omnipath.
# 
# 

# In[102]:


collectri = dc.op.collectri(organism="mouse")
collectri


# ## Scoring
# Pathway scores can be readily computed by running the ulm method.

# In[103]:


# Run
tf_acts, tf_padj = dc.mt.ulm(data=data, net=collectri)

# Filter by sign padj
msk = (tf_padj.T < 0.05).iloc[:, 0]
tf_acts = tf_acts.loc[:, msk]

tf_acts


# In[104]:


dc.pl.barplot(data=tf_acts, name="disease.vs.normal", figsize=(12, 8), top = 50)


# In[108]:


dc.pl.network(
    net=collectri,
    data=data,
    score=tf_acts,
    sources=["Stat1", "Elk1", "Stat5b", "Myc",'Id2','Stat3'],
    targets=5,
    figsize=(5, 5),
    vcenter=True,
    by_abs=True,
    size_node=15,
)


# In[109]:


dc.pl.volcano(
    data=results_df,
    x="log2FoldChange",
    y="pvalue",
    net=collectri,
    name="Stat1",
    top=10,
    figsize=(6, 6),
)


# ## Pathway Scoring
# The same approach used for TF scoring can also be applied to pathways. Numerous databases provide curated pathway gene sets, with one of the most well-known being MSigDB, which includes several collections {cite:p}msigdb. These and many other resources can be accessed using the function {func}decoupler.op.resource. To view the list of available databases, use {func}decoupler.op.show_resources.
# 
# ## PROGENy Pathway Genes
# PROGENy is a comprehensive resource that provides a curated collection of pathways and their target genes, along with weights for each interaction {cite:p}progeny.
# 
# Below is a brief description of each pathway:
# 
# Androgen: involved in the growth and development of the male reproductive organs
# EGFR: regulates growth, survival, migration, apoptosis, proliferation, and differentiation in mammalian cells
# Estrogen: promotes the growth and development of the female reproductive organs
# Hypoxia: promotes angiogenesis and metabolic reprogramming when O2 levels are low
# JAK-STAT: involved in immunity, cell division, cell death, and tumor formation
# MAPK: integrates external signals and promotes cell growth and proliferation
# NFkB: regulates immune response, cytokine production and cell survival
# p53: regulates cell cycle, apoptosis, DNA repair and tumor suppression
# PI3K: promotes growth and proliferation
# TGFb: involved in development, homeostasis, and repair of most tissues
# TNFa: mediates haematopoiesis, immune surveillance, tumour regression and protection from infection
# Trail: induces apoptosis
# VEGF: mediates angiogenesis, vascular permeability, and cell migration
# WNT: regulates organ morphogenesis during development and tissue repair

# In[110]:


progeny = dc.op.progeny(organism="mouse")
progeny


# In[111]:


# Run
pw_acts, pw_padj = dc.mt.ulm(data=data, net=progeny)

# Filter by sign padj
msk = (pw_padj.T < 0.05).iloc[:, 0]
pw_acts = pw_acts.loc[:, msk]

pw_acts


# In[112]:


dc.pl.barplot(data=pw_acts, name="disease.vs.normal", figsize=(3, 3))


# In[113]:


import numpy as np


# In[114]:


# Transform to df
df = pw_acts.melt(value_name="score").merge(
    pw_padj.melt(value_name="pvalue")
    .assign(logpval=lambda x: x["pvalue"].clip(2.22e-4, 1))
    .assign(logpval=lambda x: -np.log10(x["logpval"]))
)
dc.pl.dotplot(df=df, x="score", y="variable", s="logpval", c="score", scale=1, figsize=(4, 4))


# In[115]:


dc.pl.source_targets(data=results_df, x="weight", y="stat", net=progeny, name="MAPK", top=15, figsize=(8, 8))


# In[116]:


_, pos_le = dc.pl.leading_edge(
    results_df,
    stat="stat",
    net=progeny[progeny["weight"] > 0],
    name="MAPK",
)
print("(+) leading edge:", pos_le[:5])
_, neg_le = dc.pl.leading_edge(
    results_df,
    stat="stat",
    net=progeny[progeny["weight"] < 0],
    name="MAPK",
)
print("(-) leading edge:", neg_le[:5])


# ## Hallmark gene sets
# Hallmark gene sets are curated collections of genes that represent specific, well-defined biological states or processes. They are part of MSigDB and were developed to reduce redundancy and improve interpretability compared to older, more overlapping gene set collections {cite:p}msigdb.
# 
# A total of 50 gene sets are provided, designed to be non-redundant, concise, and biologically coherent.
# 
# This is how to access them.

# In[117]:


hallmark = dc.op.hallmark(organism="mouse")
hallmark


# In[118]:


# Run
hm_acts, hm_padj = dc.mt.ulm(data=data, net=hallmark)

# Filter by sign padj
msk = (hm_padj.T < 0.05).iloc[:, 0]
hm_acts = hm_acts.loc[:, msk]

hm_acts


# In[119]:


dc.pl.barplot(data=hm_acts, name="disease.vs.normal", figsize=(6, 8))


# In[120]:


# Tranform to df
df = hm_acts.melt(value_name="score").merge(
    hm_padj.melt(value_name="pvalue")
    .assign(padj=lambda x: x["pvalue"].clip(2.22e-16, 1))
    .assign(padj=lambda x: np.log10(x["pvalue"]))
)
dc.pl.dotplot(df=df, x="score", y="variable", s="padj", c="score", scale=0.25, figsize=(6, 6))


# In[121]:


_, le = dc.pl.leading_edge(
    results_df,
    stat="stat",
    net=hallmark,
    name="HYPOXIA",
)
print("leading edge:", le[:5])


# In[ ]:




