## 2. Unsupervised clustering and model-guided annotation reveal transcriptionally distinct oligodendrocyte states

To resolve the cellular composition of the Xenium **mtDNA double-strand break (mtDSB)** model, we performed unsupervised clustering of all captured cells following dimensionality reduction and neighborhood graph construction. **Leiden clustering** identified discrete transcriptional populations that were visualized using **UMAP** embeddings, revealing clear segregation of glial, neuronal, and stromal lineages (Fig. 2A).  

Automated cell-type annotation was achieved using a **marker-informed, large language model (LLM)-based approach**, in which cluster-specific marker genes were supplied to the `annotate_clusters()` function. The model (OpenAI `gpt-4o-mini`) integrated marker information, species, and brain context to assign biologically coherent labels, identifying the major expected CNS cell types including **oligodendrocytes**, **oligodendrocyte precursor cells (OPCs)**, **astrocytes**, **microglia**, **neurons**, **endothelial**, and **fibroblast-like** populations (Fig. 2B). The integration of automated, context-aware annotation accelerated curation while maintaining transparency via the underlying marker gene sets.  

To further resolve heterogeneity within the oligodendrocyte lineage, we subclustered the annotated oligodendrocyte population. This analysis revealed a continuum of transcriptional states spanning immature precursor-like to mature myelinating oligodendrocytes, as well as a distinct **damage-associated (DA) oligodendrocyte** cluster that emerged exclusively in mtDSB tissues (Fig. 2C).  

**Differential expression analysis** across these subclusters identified two dominant transcriptional trajectories:  
(1) a **differentiation axis**, progressing from OPC-like and immature oligodendrocyte states (Immature I–IV) toward fully myelinating cells, and  
(2) a **stress-adaptation axis**, diverging toward the DA-OL population characterized by integrated stress response (ISR) activation and mitochondrial stress signatures.  

Within the **DA-OL cluster**, genes associated with the ISR (*Atf4*, *Atf5*, *Trib3*, *Hspa9*), oxidative defense (*Gstp1*, *Mt2*), and immune activation (*H2-D1*, *B2m*, *Serpina3n*) were strongly upregulated, consistent with an adaptive but hypomyelinating phenotype. By contrast, **mature oligodendrocytes** showed robust expression of classical myelin and lipid synthesis genes (*Mbp*, *Plp1*, *Mag*, *Mog*, *Ugt8a*), indicative of intact myelin maintenance programs (Fig. 2D).  

Intermediate states between these extremes reflected distinct phases of oligodendrocyte maturation. **Immature I** cells expressed *Ptprz1*, *Pdgfra*, and *Sox4*, resembling differentiating OPCs with active migratory and axon–glia signaling. **Immature II** cells were enriched for *Sirt2*, *Cnp*, *Plp1*, and *Mbp*, marking the onset of myelin gene induction and biosynthetic upregulation. **Immature III–IV** states, while still expressing OPC-related transcripts, showed elevated expression of inhibitory and axon-guidance genes (*Lingo1*, *Rgma*, *Sema5a*), suggesting partial arrest or reactive remodeling under stress conditions.  

Overall, this transcriptional continuum delineates a **bifurcation between productive differentiation and stress adaptation**, in which mitochondrial DNA damage drives oligodendrocytes toward an ISR-dominant, damage-associated state rather than terminal myelination.  

---

**Figure 2.**  
![UMAP embedding](/results/figures/umap_cell_class.png)
*(A)* UMAP embedding of all Xenium mtDSB cells showing major annotated cell types.  
![Spatial maps](/results/figures/spatial_cell_class.png)
*(B)* Spatial maps of annotated populations across tissue sections.  
![Spatial maps](/results/figures/OL_sub_clust.png)
*(C)* Oligodendrocyte subclustering reveals immature, differentiating, and damage-associated states.  
![Dotplot](/results/figures/dotplot__OL_ANNO_cellclass_markers.png)
*(D)* Marker gene ranking highlighting lineage and stress-specific signatures.  