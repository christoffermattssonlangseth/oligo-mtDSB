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

## 🔄 Updates
1. 🧩 **Cluster annotation** completed at *leiden* resolution 2.  
2. 🗺️ **Anatomical compartments** have been defined and annotated.  
3. 📊 **Differentially expressed genes (DEGs)** identified across both **condition** and **time**.  

## Preliminary findings:🧬 Mitochondrial DNA Damage and Oligodendrocyte Stress Programs
## 🧬 Preliminary findings: Mitochondrial DNA Damage and Oligodendrocyte Stress Programs

In this project, we investigate how mitochondrial DNA double-strand breaks (mtDSBs) affect oligodendrocytes (OLs) and their surrounding microenvironment. Our analysis reveals that OLs exposed to mtDNA damage do not simply undergo cell death but instead activate a **coordinated multi-pathway stress and communication program**, combining cell-intrinsic adaptation with extracellular signaling.

### 🔍 Core stress modules identified

- 🔥 **Oxidative stress defenses** — Induction of *Mt2*, *Gstp1*, *Sqstm1*, and transcriptional control via *Nfe2l1*, reflecting activation of antioxidant and proteostasis mechanisms.  
- ⚡ **Mitochondrial and ER stress signaling** — Hallmarks of the integrated stress response (ISR) and mitochondrial unfolded protein response (UPRmt), with activation of *Atf4*, *Jun*, *Atf5*, *Hspa5*, *Hspd1*, and *Hspa9*.  
- 🧩 **Mitokine-like secretome** — Upregulation of secreted stress mediators *Gdf15*, *Adm*, *Cst7*, *Igfbp3*, and *Serpina3n*, indicating that stressed OLs release paracrine signals capable of influencing neighboring glia and immune cells.  
- 🛡️ **Antigen presentation and immune visibility** — Upregulation of *B2m*, *H2-D1*, *H2-K1*, and *Ctss*, suggesting that metabolically stressed OLs may become visible to the adaptive immune system.  
- 🚚 **Transport and myelin support disruption** — Altered kinesin/dynein motors (*Kif5a*, *Kif5b*, *Dync1li1*) and adhesion molecules (*Itgb1*, *Mpzl1*) indicate structural and metabolic remodeling of the myelin unit.  
- 🌟 **Reactivity and glial crosstalk** — Increased *Ndrg2*, *Gfap*, and calcium-binding proteins (*S100a1*, *Calb2*) highlight astrocyte-like reactivity within the OL lineage.  
- 📣 **Inflammatory mediators** — Elevated *Nmu* (neuromedin U) and *Ccl3* point to recruitment and activation of immune cells.

**✨ Overall, these findings show that mitochondrial genome instability in oligodendrocytes drives a “stressed-but-not-dead” phenotype**, characterized by oxidative stress, integrated stress signaling, antigen presentation, and secretion of mitokine-like factors. This creates a permissive environment for **immune activation and tissue remodeling**, providing a mechanistic link to multiple sclerosis (MS)–like pathology.

### Differential gene expression reveals cell type–specific changes in integrated stress response, cytokine and inflammatory signaling, and metabolic/mitochondrial reprogramming

Mitochondrial DNA double-strand breaks (mtDSBs) in oligodendrocytes trigger a coordinated, cell type–specific transcriptional reprogramming across the CNS. Differential gene expression analysis revealed broad activation of the **integrated stress response (ISR)**, **cytokine and inflammatory signaling**, and **metabolic/mitochondrial remodeling**. Mature oligodendrocytes and microglia exhibited the strongest responses, characterized by upregulation of *Atf5*, *Trib3*, *Cdkn1a*, *Cst7*, *Cd40*, and *Gdf15*, reflecting a transition from mitochondrial stress adaptation to chronic immune activation. Astrocytes showed region-dependent engagement, with olfactory astrocytes mirroring the oligodendrocyte ISR signature, while telencephalic astrocytes and OPCs displayed more constrained transcriptional shifts.  

At later stages, oligodendrocytes developed a **mitokine-like secretome**, prominently featuring *Gdf15*, *Adm*, *Cst7*, and *Serpina3n*, suggesting active communication of mitochondrial distress to neighboring glia. This paracrine signaling likely underlies the widespread inflammatory and metabolic reprogramming observed in astrocytes and microglia, linking oligodendrocyte mitochondrial dysfunction to non–cell-autonomous glial activation and CNS remodeling. 📖 **Read more:** [Detailed on DGE analysis](results/preliminary-results/2-agexcondition-glia.md)

---

## 🧠 Region-independent effects (Oligodendrocytes)

### 21 weeks (early phase)
- **Overall pattern:** Mild transcriptional response dominated by early stress and differentiation shifts.  
- **Upregulated genes:** *Mt1*, *Trh*, *Nkx6-2*, *Cst7*, *Otp*  
- **Downregulated genes:** *Hoxa2*, *Hoxb3*, *Hoxb5*, *Col1a1*, *Bglap*, *Trim36*  

🧩 **Interpretation:**  
Early mtDNA damage elicits a subtle oxidative and mitochondrial stress response (*Mt1*), with suppression of developmental/myelin structural programs. This likely represents an initial compensatory or remodeling phase preceding full ISR activation.

---

### 60 weeks (late phase)
- **Overall pattern:** Broad activation of the integrated stress response (ISR), UPRmt, and inflammatory signaling.  
- **Upregulated genes:** *Atf5*, *Trib3*, *Cdkn1a*, *Gdf15*, *Adm*, *C4b*, *Cd40*, *Cst7*, *Gpnmb*, *Serpina3n*  
- **Downregulated genes:** *Pmch*, *Hcrt*, *Cercam*  

🧩 **Interpretation:**  
At later stages, oligodendrocytes show a robust ISR/UPRmt profile (*Atf5*, *Trib3*, *Gdf15*, *Cdkn1a*) together with glial activation (*Gfap*, *Gpnmb*) and complement/inflammatory pathways (*C4b*, *Cd40*).  
Notably, *Gdf15* was strongly induced (log₂FC ≈ 4.7, p ≈ 0.001), whereas *Fgf21*—another canonical mitokine—remained absent, highlighting a **CNS-adapted ISR^mt** centered on *Gdf15*-mediated paracrine signaling rather than systemic *Fgf21* induction. This pattern reflects a transition from **mitochondrial stress adaptation** to **chronic inflammatory signaling** and glial cross-communication.

---

## 🧬 Comparative summary of regional responses to mtDNA DSBs (Oligodendrocyte model)

| Category | **Hypothalamus** | **Thalamus** | **Pallidum + Corpus Callosum** |
|-----------|------------------|--------------|--------------------------------|
| **Neuronal (Hcrt / Pmch)** | 🔻 **Strong down** | 🔻 Mild | ⛔ None |
| **Stress (Atf5 / Trib3)** | 🔺 Moderate | 🔺 Strong | 🔺 Strong |
| **Complement / Inflammatory (C4b, Cst7, Cd74)** | 🔺 Strong | 🔺 Medium | 🔺 Medium |
| **Astro / Microglial activation (Gfap, Gpnmb)** | 🔺 Strong | 🔺 Medium | 🔺 Medium |
| **Metabolic adaptation (Gdf15, Adm)** | 🔺 Medium | 🔺 Medium | 🔺 Medium |
| **White-matter–specific response** | 🟠 Some | 🟠 Some | 🔴 Strong |

**Interpretation:**  
- The **hypothalamus** exhibits the most pronounced and multifaceted response—strong glial stress and inflammatory activation accompanied by loss of neuropeptidergic markers (*Hcrt*, *Pmch*).  
- The **thalamus** shows a robust glial ISR/stress pattern but minimal neuronal involvement.  
- The **pallidum + corpus callosum** displays a predominantly **cell-intrinsic oligodendrocyte stress response**, characteristic of white-matter regions.

---

### 🧠 Summary of mechanistic insight

Together, these data indicate that **mtDNA damage in oligodendrocytes activates an integrated ISR^mt–UPRmt axis coupled to a mitokine-like secretome**, dominated by *Gdf15* and *Adm*. This secretory program likely acts in a **paracrine manner**, modulating astrocyte and microglial states, and may represent a CNS-specific adaptation of mitochondrial stress signaling that links **oligodendrocyte mitochondrial dysfunction to broader neuroinflammatory remodeling**.
---