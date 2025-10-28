
## 6. Regional activation fingerprints reveal cell type– and age–specific stress responses
To assess how mitochondrial DNA double-strand breaks (mtDSBs) remodel gene-expression programs across the CNS, we derived **regional activation fingerprints** for each major glial population.  
Differential expression was aggregated per *cell type × anatomical compartment × age*, restricted to upregulated and sufficiently expressed genes (*log₂FC > 0.5*, *baseMean > 50*, *FDR < 0.05*).  
An **activation index**, defined as the product of the number of significant genes and their mean absolute log₂ fold change, quantified both the breadth and magnitude of activation per region.

Across glial cell types, the results revealed pronounced **age- and region-specific transcriptional remodeling**:

- **Oligodendrocytes** transitioned from modest, developmentally linked activation at P21 (notably in hypothalamus and fiber tracts; *Hcrt*, *Trh*, *Syt4*) to strong ISR- and immune-related activation at P60 (*Cdkn1a*, *Trib3*, *Serpina3n*, *Cd74*), particularly in thalamic and cortical regions.
- **Astrocytes** showed broad stress induction across both ages, with hypothalamic activation at P21 (*Hcrt*, *Pmch*, *Sst*) shifting to widespread inflammatory and complement signatures at P60 (*Thbs4*, *Gfap*, *C4b*, *Serpina3n*).
- **Microglia** exhibited minimal activation at P21 but a striking late response at P60, dominated by *Cd74*, *Cst7*, *Lpl*, and *Gfap*, consistent with a **disease-associated microglia (DAM-like)** phenotype in thalamic and hypothalamic compartments.
- **Oligodendrocyte precursor cells (OPCs)** remained largely quiescent, showing only mild activation at P60 (*Gfap*, *Hspa5*, *Mertk*), suggesting partial engagement of stress pathways without full arrest.

Collectively, these findings indicate that mtDSB-induced stress elicits a **regionally organized, cell type–specific activation hierarchy** that intensifies with cellular maturation.  
White-matter and metabolically active regions, particularly the **thalamus, hypothalamus, and fiber tracts**, emerge as hotspots of late-stage glial activation, reflecting increased susceptibility to mitochondrial genome perturbation in the adult CNS.

---

### **Figure 6. Regional activation fingerprints by age and cell type**
**(A)** Overview of the activation index computation and visualization approach.  
**(B)** Representative cell type–specific activation fingerprints comparing P21 and P60 mtDSB induction.  

**Oligodendrocytes**  
![Oligodendrocytes](/results/figures/region_fingerprints_by_age/Oligodendrocytes_activation_by_age.png)

**Astrocytes**  
![Astrocytes](/results/figures/region_fingerprints_by_age/Astrocytes_activation_by_age.png)

**Microglia**  
![Microglia](/results/figures/region_fingerprints_by_age/Microglia_activation_by_age.png)

**Oligodendrocyte Precursor Cells (OPCs)**  
![OPCs](/results/figures/region_fingerprints_by_age/Oligodendrocyte_Precursor_Cells_activation_by_age.png)

### Neural stem cells: indirect sensing of mitochondrial dysfunction in the ventricular niche

Although mitochondrial DNA damage was restricted to oligodendrocytes, **neural stem cells (NSCs)** within the **ventricular system** exhibited one of the most striking secondary transcriptional responses across all analyzed cell types.  
At **P21**, their activation was minimal and largely confined to stress-adaptive and barrier-associated transcripts (*Cldn10, Gpc4, Tsc22d3*), consistent with a **homeostatic response** to subtle metabolic perturbations in the surrounding tissue.  

By **P60**, however, the same ventricular niche showed a **robust and coordinated transcriptional activation** dominated by **cell cycle regulators and neurogenic transcription factors** (*Ccna2, Dlx2, Mki67, Myt1, Meis2, Dcx, Stmn4, Sox1–4*).  
This pattern is emblematic of a **latent progenitor reawakening or pseudo-neurogenic activation**, occurring despite the absence of direct mtDNA lesions in these cells.  

We interpret this as evidence that **adult NSCs are metabolically and paracrinely sensitive to glial mitochondrial stress**. Signals originating from mtDSB-affected oligodendrocytes—potentially involving lactate, cytokines, or ISR-related factors—appear to propagate to the ventricular niche, triggering a **reactive proliferation-like transcriptional state**.  
Thus, the NSC compartment behaves as a **sensor of the metabolic and inflammatory milieu**, translating glial mitochondrial dysfunction into a regenerative-like transcriptional response.

---

**Figure. Neural stem cells: regional activation fingerprints by age**  
![Neural stem cells](/results/figures/region_fingerprints_by_age/Neural_Stem_Cells_activation_by_age.png)