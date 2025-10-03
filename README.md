# oligo-mtDSB  
<p align="center">
  <img src="assets/logo.png" alt="oligo-mtDSB Logo" width="300"/>
</p>
This project investigates the consequences of double-stranded breaks (DSBs) in mitochondrial DNA (mtDNA) of oligodendrocytes using a genetically engineered mouse model. Xenium is applied to capture cellular and regional responses to mtDNA damage in the CNS.  

## Link to original publication and thesis
Link: https://pmc.ncbi.nlm.nih.gov/articles/PMC5647772/ and https://scholarship.miami.edu/view/pdfCoverPage?instCode=01UOML_INST&filePid=13386228880002976&download=true


## Project Overview  
- **Model**: Mouse with induced mtDNA double-strand breaks in oligodendrocytes  
- **Focus**: Linking mitochondrial genome instability to oligodendrocyte dysfunction, demyelination, and lesion formation  
- **Data**: 10x Genomics Xenium (spatial transcriptomics)  
- **Goal**: Define how mtDNA integrity impacts oligodendrocyte biology and contributes to neuroinflammation and MS-like pathology  

## Repository Structure  
- `notebooks/` – analysis workflows (QC, clustering, differential expression, spatial mapping)  

## Analysis Highlights  
- Cell type–resolved effects of mtDNA DSBs  
- Spatial mapping and proximity-dependent transcriptomic changes  
- Pathway-level insights into mitochondrial stress and myelin integrity  

## Preliminary findings:🧬 Mitochondrial DNA Damage and Oligodendrocyte Stress Programs

In this project, we investigate how mitochondrial DNA double-strand breaks (mtDSBs) affect oligodendrocytes (OLs) and the surrounding microenvironment. Our analysis reveals that OLs exposed to mtDNA damage do not simply die, but instead activate a **coordinated multi-pathway stress response**:

- 🔥 **Oxidative stress defenses** — *Mt2, Gstp1, Sqstm1* and antioxidant transcriptional control (*Nfe2l1*).  
- ⚡ **Mitochondrial and ER stress signaling** — integrated stress response (ISR) and unfolded protein response (UPRmt) (*Atf4, Jun, Hspa5, Hspd1, Hspa9*).  
- 🛡️ **Antigen presentation pathways** — upregulation of MHC-I molecules (*B2m, H2-D1, H2-K1*) and antigen processing enzymes (*Ctss*), suggesting stressed OLs become visible to the immune system.  
- 🚚 **Transport and myelin support disruption** — altered kinesin/dynein motors (*Kif5a, Kif5b, Dync1li1*) and adhesion molecules (*Itgb1, Mpzl1*).  
- 🌟 **Reactivity and glial crosstalk** — induction of *Ndrg2*, *Gfap*, and calcium-binding proteins (*S100a1, Calb2*).  
- 📣 **Inflammatory mediators** — *Nmu* (neuromedin U) and *Ccl3*, pointing to signaling that may recruit and activate immune cells.

**✨ Overall, these findings show that mitochondrial genome instability in oligodendrocytes drives a “stressed but not dead” phenotype**, characterized by oxidative stress, antigen presentation, and glial–immune interactions. This creates a permissive environment for **immune activation and tissue remodeling**, providing a mechanistic link to multiple sclerosis (MS)-like pathology.