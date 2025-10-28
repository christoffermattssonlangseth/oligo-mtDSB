## 4. Spatial domain discovery and compartment calling

To define spatial compartments within the mtDSB tissues, we utilized a previously established **read-based domain discovery workflow**, implemented in Python. This pipeline operates directly on the Xenium-derived per-cell transcriptomes while preserving spatial context.  

In brief, for each sample, we generated a *pseudobinned* expression matrix by aggregating each cell’s transcript counts with those of its spatial neighbors (typically the **20 nearest cells**, determined using **Squidpy**). This local neighborhood summation smooths cell-to-cell variability into spatially coherent transcriptional patches while retaining cell-level metadata and coordinates.  

Following pseudobinning, low-count cells were filtered, and the resulting objects were embedded using a standard **Scanpy** workflow involving principal component analysis (PCA), neighborhood graph construction, and **UMAP** embedding. Spatial domains were then identified by **unsupervised community detection** (Leiden or Louvain clustering) across a range of resolutions. Each resulting cluster—termed a **read-based domain (RBD)**—represents a transcriptionally defined, spatially contiguous compartment.  

These domains were subsequently annotated using the **Allen Brain Reference Atlas**, integrating both their spatial localization and transcriptional signatures.  

In total, we identified **27 distinct anatomical and transcriptional compartments** spanning both gray and white matter regions, including the **cortex (layers I–VI and subplate)**, **caudoputamen**, **thalamic nuclei**, **hypothalamus**, **hippocampal formation**, **dentate gyrus**, **olfactory areas**, **pallidum**, and **striatal ventral regions**. Additional compartments corresponded to **vascular and fiber tract–associated regions** (e.g., cortical vasculature, parenchymal vasculature, corticospinal tracts, and fiber tracts I–II), as well as border regions such as the **meningeal/glia limitans** and **ventricular system**. A small number of **unclassified (unknown)** domains likely represent transitional or low-density areas.  

By applying the same read-based analysis to tissues in which mtDNA double-strand breaks were induced at **postnatal day 21 (P21)** versus **postnatal day 60 (P60)**, we could directly compare the emergence and composition of these domains between the **developing (actively myelinating)** and **mature (myelin-maintaining)** oligodendrocyte environments. This revealed that the global spatial organization of transcriptional compartments remains largely conserved, although local cellular composition and stress gene enrichment vary markedly between developmental stages.

---

### **Figure 4. Spatial domain discovery and anatomical compartmentalization**

**(A)** Overview of the read-based domain (RBD) discovery workflow combining pseudobinning, dimensionality reduction, and community detection to resolve transcriptionally contiguous regions.  
**(B)** Spatial compartment map showing 27 distinct anatomical and transcriptional domains identified across Xenium mtDSB tissue sections. Major regions correspond to cortical layers, subcortical nuclei, fiber tracts, and vascular or meningeal compartments.  
![Spatial compartments](/results/figures/compartments.png)
**(C)** Comparison of P21 versus P60 mtDSB tissues highlights conserved large-scale organization but altered composition and stress gene enrichment within specific compartments.

