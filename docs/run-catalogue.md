# Run catalogue

One entry per analysis run of the oligo-mtDSB project, numbered once (cite as §N). Built 2026-10-23 from code, notebooks, result files and docs; feeds the analysis map in docs/analysis-map/.

### 1. Cell typing yields 15 cell classes, including a DA-OL state that is a quarter to two thirds of OLs in three mtDSB mice
*Notebook:* legacy/notebooks-01/00-batch-processing.ipynb, legacy/notebooks-02/01-clustering-take2.ipynb, legacy/notebooks-02/02-annotation-3-clusters.ipynb, legacy/notebooks-02/13-cytetype-annotation.ipynb, legacy/notebooks-04/06-subclustering-OLs.ipynb, legacy/notebooks-04/00-inspect-data.ipynb, legacy/notebooks-02/17-novae-oligo-mtDSB.ipynb (also legacy/notebooks-01/01-03, 01/17, 02/06) · *Results:* `results/cell_table/cells.csv.gz`, `results/cell_table/sections.csv`

**Rationale.** Every later analysis needs a cell type per cell. The labels were built in several rounds between 2025-09 and 2026-01 (git log), each one replacing the last.

**Technical approach.** The 12 Xenium output folders (slides 0060539 = P60, 0060541 = P21) were concatenated (01/00): 983,582 cells before QC. Cells with < 40 counts or < 15 genes were removed, leaving 980,474 cells x 5,101 genes. Counts were normalised to 100 per cell and log1p-transformed, then PCA, a 15-neighbour / 30-PC graph, UMAP and Leiden at resolutions 0.5-2 (01/00) and up to 3 (02/01). Three annotation rounds followed:
(a) manual names for leiden_2 clusters, 22 classes such as "Telencephalon astrocytes", "Olfactory astrocytes", "Mature oligodendrocytes (mic)" (01/17); the legacy notebooks 01/21-30 use these;
(b) top-20 t-test markers of the 61 leiden_3 clusters sent to an LLM (`mllmcelltype.annotate_clusters`, OpenAI gpt-4o-mini, tissue "brain") and collapsed to 12 coarse classes (02/02; manuscript sec. 2 describes this);
(c) CyteType on leiden_3 with a study-context prompt, giving a separate 61-category column `cytetype_annotation_leiden_3` (02/13).
OL-lineage cells (Oligodendrocytes + OPC classes) were then sub-clustered: PCA 50, 20 neighbours / 30 PCs, Louvain resolution 1.3 (18 clusters), with clusters named by marker genes (04/06, run on the novae3 file). An earlier version of this step (02/06) gave a different split. 04/00 adds `cell_class_updated` (Immature I-IV merged) to the 250-gene MANA file. Novae (02/17) gave zero-shot spatial domains (`novae-mouse-0`, radius 80, level 15). They do not feed `cell_class`.

**Results.** Cells per section after QC range from 69,913 (RB4630) to 135,378 (RB4403); RB4403 is later called over-segmented (README). The 04/06 sub-cluster counts are Oligodendrocytes 145,940, Immature OL I 22,087, DA-Oligodendrocytes 15,261, Immature OL II 3,119 and Immature OL III 1,584. They match the current `cell_class` in `results/cell_table/cells.csv.gz` exactly. That file has 15 classes: Excitatory Neurons 210,407, Neurons 160,462, Astrocytes 155,119, Oligodendrocytes 145,940, Endothelial 127,681, GABAergic 50,052, Medium Spiny Neurons 47,381, Microglia 28,916, Immature OL I-III, DA-OL, NSC 5,022, Endocrine 4,922 and Myocytes 2,521. DA-OL top markers (04/06) are Gstp1, Car2, Atf4, Atf5, Adamts1, Bcl2l1, H2-D1, Serpina3n and Hspa9. DA-OL share of OL + DA-OL per mouse (computed from cells.csv.gz for this catalogue): RB4405 64.7%, RB4627 42.6%, RB4401 21.4%, all others 0.26-2.3%.

**Reading.** The DA-OL cluster is the responder OL state found later (§10-11).

**Caveats.** Labels come from an LLM, CyteType and manual naming, with no reference mapping in the final labels. The Tangram / Allen notebooks 01/13-16 were not used. "Immature Oligodendrocytes I" has OPC markers (Tnr, Ptprz1, Pdgfra, Cspg5), and the current `cell_class` has no class called OPC. The `markers` printout in 04/06 cell 44 is stale (it repeats the 02/06 numbers), and the dotplot save in 04/06 fails with a NameError. README "Analysis history" item 1 cites 02/13 and 04/00, but the coarse labels come from 02/02 and the OL sub-labels from 04/06. Clustering pools all 12 mice, so a disease state can form its own cluster (it did: DA-OL).

### 2. Read-based compartments assign every cell to one of 12 anatomical regions
*Notebook:* legacy/notebooks-02/03-read-based-compartments.ipynb, legacy/notebooks-01/19-compartment-cell-driven-all-samples.ipynb, legacy/notebooks-01/20-compartment-read-based-all-sections.ipynb, legacy/notebooks-01/34-RBD-subset.ipynb · *Code:* `scripts/compartments/run_compartment_read_based.py`, `scripts/compartments/run_compartment_read_based.sh`, `scripts/compartments/submit_compartment_read_based.sub` · *Results:* `results/cell_table/cells.csv.gz`

**Rationale.** The regional DE (§3, §8) and all later region-level work need an anatomical label per cell.

**Technical approach.** A first attempt (01/19) clustered cell-type counts within 200 µm of each cell (Leiden) and named about 20 compartments by hand. It was superseded by read-based domains (RBD). The RBD script runs per section. It adds each cell's expression to that of its k nearest spatial neighbours (squidpy, `--neighbors` default 20, plus the cell itself). It then removes pooled profiles with < 4 counts and runs PCA (30), a 20-neighbour graph, UMAP and Leiden at 0.1, 0.2, 0.5 and 1.0. It ran on an HTCondor cluster (8 CPU, 128 GB) on `mtDNA_DSB_5k_clustered_manual_annotation.h5ad`. In 02/03 the `rbd_domain_1.0` clusters (about 50) were named by hand against the Allen atlas into `RBD_compartment`. The manuscript sec. 5 counts 27 named compartments. These were then collapsed to 12 regions in `RBD_compartment_simplified`. 01/34 repeats the step on a subset at resolution 0.5.

**Results.** Cells per region (all 12 sections): Cortex 259,386, Thalamus 130,938, Striatum 124,957, Hypothalamus 109,934, Fiber tracts 77,564, Olfactory areas 72,098, Meningeal border/glia limitans 62,048, Ventricular system 44,863, Hippocampus 35,364, Pallidum 33,663, Vasculature 8,085, Unknown 21,574. "Unknown" is excluded downstream.

**Reading.** The map is anatomically plausible and was used by every regional analysis up to analysis/14.

**Caveats.** All genes went into the pooled profiles, disease genes included. The input X of that file was log-normalised (01/00, 01/17); raw counts were attached only later (01/21), so the pooling summed log values. Naming was manual (the code spells it "Palladium" before simplification). **Corrected 2026-10-21 (§43):** in responders the disease moved pallidum into "Fiber tracts" (pallidum 0-0.7% of cells vs 3-7%). Use `region_new` from `results/region_map/region_transfer_labels.csv.gz` for new region-level work.

### 3. Earlier group-level DE reported a "P21 modest, P60 strong" OL stress response; that pattern came from responder composition
*Notebook:* legacy/notebooks-01/05-tissue-wide-expression-differences.ipynb, legacy/notebooks-01/06-predicting-sample-condition-OL.ipynb, legacy/notebooks-01/11-comparing-time-points.ipynb, legacy/notebooks-01/18-cell-type-log2fc.ipynb, legacy/notebooks-01/21-more-DGE.ipynb, legacy/notebooks-01/23-pertpy-OL.ipynb, legacy/notebooks-01/24-pertpy-agexcondition-oligo.ipynb.ipynb, legacy/notebooks-01/25-agexcondition-glia.ipynb, legacy/notebooks-01/26-agexconditionxregion-OLs.ipynb, legacy/notebooks-02/04-region-DGE.ipynb, legacy/notebooks-02/08-ISR-celltype-specific.ipynb · *Results:* `legacy/results/age/Mature_oligodendrocytes/simple_effect_age_21.csv`, `legacy/results/age/Mature_oligodendrocytes/simple_effect_age_60.csv` · *Figure:* `legacy/results/figures/pathway_ISR_UPR_age21.png`, `legacy/results/figures/pathway_ISR_UPR_age60.png`

**Rationale.** These runs were the first answer to "what does mtDSB do, and does it depend on induction age?". They fed manuscript sections 3, 4 and 6.

**Technical approach.** Methods changed from notebook to notebook:
(i) Blind per-section markers (01/05, 01/06): scanpy t-test of each section against all others, run before metadata was attached. The notebooks guessed conditions, e.g. control RB4282 / RB4350 as "DSB".
(ii) Cell-level t-tests and mean-count comparisons of mtDSB vs control OLs, overall, by condition x age (01/11) and at P60 only (01/21).
(iii) Pseudobulk pyDESeq2 via pertpy / decoupler. 01/23 pooled 10 OL sub-types x 12 sections into 120 "samples" with `~condition`. 01/24 used `~age + condition + age:condition` on mature OL / OPC pseudobulks.
(iv) pyDESeq2 per glia class per age (01/25 -> `legacy/results/age/`), and per class x region x age (01/26, 02/04). In 02/04, 4,432 of 5,101 genes passed gene QC.
(v) ISR / UPR gene panels plotted per cell type and age (02/08).

**Results.** Per-section top genes were stress-hormone / immediate-early genes (Sgk1, Arc, Ddit4, Fos, Tsc22d3 in RB4282), Kif5a / Map2 (RB4403) and Xist (female sections). RB4405 already showed Gfap, Serpina3n, Gstp1, Atf4 and Atf5. OL cell-level mtDSB vs control was led by control-high Il33 and Sez6l2. The pooled 01/23 model was led by Il33 (log2FC -2.13), Prl (-4.15), Fkbp5 (+0.69), Ccl3 (+2.60) and Atf5 (+0.81). In 01/24 the condition term had no gene with padj below 0.63 among those printed. The top interaction genes were the neuropeptides Hcrt (log2FC -8.44, adj p 1.1e-9) and Pmch (-6.01, 1.8e-5). The notebook nevertheless concluded that P60 OLs shift to "inflammatory and dedifferentiation genes (Irf8, Il1a, Hoxb3, Fcgr2b)".
Genes at padj < 0.05 in `legacy/results/age`:
- Mature OL: 4 at P21, 14 at P60
- Microglia: 1 and 3
- Olfactory astrocytes: 2 and 10
- OPC: 6 and 2
- Telencephalon astrocytes: 3 and 4

Mature-OL Trib3 was log2FC 3.89 ± 1.78 at P21 (padj NaN) and 5.10 ± 1.45 at P60 (padj 0.080). In 02/04 the OL regions ranked by summed -log10 q were Hypothalamus 110.4, Thalamus 38.0, Olfactory 8.7, Fiber tracts 7.2 and Striatum 4.6. 02/08 concluded "P21: ISR activation is modest and cell-restricted; P60: robust".

**Reading.** The stress genes were real, but group-level DE diluted them across non-responders. The age pattern reflects one responder at P21 against two at P60 (§7, §9, §11).

**Caveats.** Cell-level tests treat cells as replicates. 01/23 pools sub-types as replicates. Hypothalamic hits are neuropeptide spill-over (Hcrt, Pmch). The per-section top genes are section-level confounds: stress hormones, circadian genes (Dbp, Nr1d1), the odd section RB4403, and Il33 from one control. Age is confounded with slide.

**Corrected 2026-10-14 (§7):** the "P21 modest / P60 strong" claim is not supported. Note that ANALYSIS_REVIEW.md quotes "Trib3 log2FC 5.1 ± 1.4" as the P21 estimate; in `legacy/results/age` that value is the P60 estimate (P21 = 3.89 ± 1.78).

### 4. Hierarchical Bayesian case-control models found small but direction-consistent OL shifts; part of the PyMC output is not a condition effect
*Notebook:* legacy/notebooks-02/09-mtDSB-casecontrol-hierarchical-bayes.ipynb, legacy/notebooks-02/11-case-control-hierach-age-condition-21.ipynb, legacy/notebooks-02/11-case-control-hierach-age-condition-60.ipynb, legacy/notebooks-02/12-comparing-ages-posterior-mean.ipynb, legacy/notebooks-01/28-run-bayes-all-genes.ipynb, legacy/notebooks-01/29-run-bayes-for-potential-genes.ipynb, legacy/notebooks-01/30-run-bayes-cond-OL.ipynb, legacy/notebooks-01/30-run-bayes-cond-OL-OPC.ipynb · *Results:* `legacy/results/baysian_diff_condition.csv`

**Rationale.** These models were meant to account for the nested structure (cells within sections within condition), which the cell-level tests of §3 ignore.

**Technical approach.**
(a) 02/09, "EUCLID-inspired", JAX / NumPyro. The model fits one gene at a time on single cells: expression = section effect + per-cell-type baseline + per-cell-type condition shift + noise. It uses stochastic variational inference (AutoNormal, 800 steps). The run covered the highly variable genes (370 features started; 360 per cell type in the CSV). Output is a posterior mean shift, a PEP, a q-value and an FDR-0.05 flag. The section key equals the sample key, because there is one section per mouse.
(b) 02/11-12: the same model by age, on all 5,101 genes.
(c) 01/28-30: PyMC per gene on pseudobulk log1p-CPM of the 22-class annotation.

**Results.** (a) Comparing the 22 analysis/03 driver genes in the OL rows of `legacy/results/baysian_diff_condition.csv` (re-checked for this catalogue): all 22 agree in direction (17 up, 5 down), and 12 are FDR < 0.05. These are Mt2 +0.051, Gstp1 +0.044, B2m, H2-D1, Atf4, Atf5, Olig2, Cldn11 -0.167, Cnp -0.151, Mal, Mog and Ugt8a. The one-sided sign test is 0.5^22 = 2.4e-7 (two-sided 4.8e-7). Effects are tiny: median |posterior mean| across OL features is 0.0014. Serpina3n is +0.012 (q 0.10).
(b) The age-stratified run is incomplete. 02/12 read 470 genes for P21 and 155 for P60; the folders now hold 495 and 181 files of 5,101.
(c) 01/28 reported 40 divergences and r-hat > 1.01. The 01/30 model is `mu = alpha + beta[condition]`, with one beta per condition, and the notebook reports `beta[1]`. That is the mtDSB group level, not the mtDSB minus control difference. So "Serpina3n β = 2.145", and module "effects" that are all positive (myelin_identity +2.91, mitochondrial_ISR +2.47), are expression levels. Ranking genes by |beta[1]| put highly expressed myelin genes first (Cnp, Cldn11, Mal).

**Reading.** The pooled NumPyro model supports the direction of today's OL programme. Its effect sizes are diluted by the three non-responders (§7).

**Caveats.** The likelihood is at cell level and section is confounded with mouse. Variational inference ran for only 800 steps. The P21 / P60 comparison rests on partial gene sets. **Conflict:** ANALYSIS_REVIEW.md ("01/30 (β = 2.1)") and the README ("PyMC models (01/29-30): myelin identity down, Serpina3n up") cite 01/30 numbers that are not condition effects. The underlying Serpina3n and myelin findings are supported per animal elsewhere (§7, §15).

### 5. Earlier OL sub-state work found the DA-OL programme (ISR, redox, MHC-I, myelin loss) that later turned out to be the responder state
*Notebook:* legacy/notebooks-01/03-sub-cluster-OLs.ipynb, legacy/notebooks-01/06-predicting-sample-condition-OL.ipynb, legacy/notebooks-01/07-pseudotime.ipynb, legacy/notebooks-01/21-more-DGE.ipynb, legacy/notebooks-02/06-subclustering-OLs.ipynb, legacy/notebooks-04/06-subclustering-OLs.ipynb · *Results:* `results/cell_table/cells.csv.gz`

**Rationale.** These runs asked whether OLs split into normal and damaged states, and where the samples sit on a damage axis.

**Technical approach.**
- OL sub-clustering (Leiden in 01/03; Louvain 1.0-1.3 in 02/06 and 04/06), with marker t-tests and scanpy scores for prespecified OPC / COP / NFOL / myelin / DOL / ISR / IFN sets.
- A blinded sample-level pseudotime (01/07): diffusion pseudotime on per-section OL pseudobulks, rooted at RB4676, plus a stress score (Sgk1, Ddit4, Tsc22d3, Hspa1a/b, Atf4, Atf5, Serpina3n, ...) and a myelin score.
- "Functional modules" from a P60-only mtDSB vs control comparison of mean OL counts (01/21).

**Results.** 02/06 called DA-OL a "stress-adaptation axis" (11,778 cells, top Gstp1, Atf4, Car2, Atf5, Bcl2l1). The 04/06 rerun gives the current 15,261 DA-OLs. DA-OL share is 64.7% / 42.6% / 21.4% of OLs in RB4405 / RB4627 / RB4401 and ≤ 2.3% elsewhere (§1). In 01/07 the pseudotime correlated with the stress score (Spearman r 0.769, p 0.0034) but not with the myelin score (r -0.217, p 0.50). The myelin score was lowest in RB4405 (-1.83), RB4627 (-0.95) and RB4401 (-0.45). 01/21 at P60 found Trib3 log2FC 4.63, Cd40 4.30, Gdf15 4.20, Cdkn1a 3.18, Nmu 2.72, Gfap 1.18, Mt2 0.72 and Gstp1 0.60. It grouped them into oxidative-stress, ISR / UPRmt, MHC-I, transport and adhesion modules.

**Reading.** The programme is right. The claim that DA-OLs "emerged exclusively in mtDSB tissue" needs to read "in 3 of 6 mtDSB mice" (§7).

**Caveats.** The pseudotime order follows age / slide: the first five sections are P21 mice. The stress score includes section-level confound genes (Sgk1, Tsc22d3, Ddit4). The 01/21 modules use cell-level means without replicate statistics. ANALYSIS_REVIEW.md gives the responder DA-OL range as "22-65%"; the current labels give 21.4% for RB4401 (minor rounding).

### 6. Pathway and external explorations (MS serum, decoupler, sccellfie, MANA, EAE integration, NSC) produced claims that mostly did not hold
*Notebook:* legacy/notebooks-02/07-mtDSB-OL-mito-inflammatory-signature-vs-MS-serum.ipynb.ipynb, legacy/notebooks-01/22-decoupler.ipynb, legacy/notebooks-01/23-pertpy-OL.ipynb, legacy/notebooks-03/02-sccellfie-test.ipynb, legacy/notebooks-03/03-sccellfie-oligomtDSB.ipynb, legacy/notebooks-04/04-MANA.ipynb, legacy/notebooks-01/12-RR-chronic-mtDSB-integration.ipynb; results/paper/7-neural-stem-cells-mount-proliferative-reponse.md

**Rationale.** These runs looked for a mechanism and for disease relevance: a mito-inflammatory signature like pre-onset MS serum (Abdelhak et al. 2025), TF and pathway activity, metabolism, spatial niches, EAE similarity, and NSC reactivation.

**Technical approach.**
- 02/07: scanpy module scores averaged by condition over 187,991 OL-lineage cells, with no replicate-level test.
- 01/22-23: decoupler ULM with CollecTRI TFs, PROGENy and Hallmark on the pyDESeq2 Wald statistics of the pooled OL model (§3).
- 03/03: sccellfie metabolic tasks (211 genes, 296 reactions, 158 tasks).
- 04/04: MANA. 250 HVGs, scVI (30 epochs), a Delaunay graph, 3-layer neighbour aggregation and Leiden 0.25 (15 clusters). It writes the 250-gene MANA h5ad.
- 01/12: concatenation and clustering with RR-EAE and chronic-EAE 5k data (1,703,389 cells).
- Manuscript sec. 7: NSC "proliferative / neurogenic response at P60", from legacy region DE (Neural Stem Cells @ Ventricular system).

**Results.** 02/07 mean module scores, control vs mtDSB:

| Module | control | mtDSB |
|---|---|---|
| ISR / UPR | 0.593 | 0.735 |
| Mito / redox | 0.453 | 0.511 |
| Glycolysis | 0.504 | 0.543 |
| cGAS-STING / ISG | 0.070 | 0.080 |
| NF-κB | 0.108 | 0.112 |
| IL-3 axis | 0.066 | 0.069 |
| Myelin | 1.801 | 1.695 |

The notebook read this as an mtDNA leakage / innate-sensing signature mirroring pre-onset MS. PROGENy hits (01/23) were MAPK 5.55, EGFR 5.46, p53 4.38 and JAK-STAT 2.83; Hallmark hits were HYPOXIA 5.70, MYC V1 5.28, UPR 4.47, P53 4.44 and TNFA-NFKB 3.44. 28 TFs were at padj < 0.05, including Nupr1 4.70 and Atf5 3.26. sccellfie, MANA and the EAE integration printed no condition-level conclusion. Manuscript sec. 7 lists Mki67, Ccna2, Dlx2, Sox2 and Dcx as up at P60 with padj < 1e-10.

**Reading.** In the review (§7), cGAS-STING / IFN, NF-κB and glycolysis show no responder pattern; only Hk2 rises. The NSC signal does not follow responders and most likely reflects section plane. Later runs confirm both: §40 (no cGAS-STING) and §39 (NSC = plane). The p53 and UPR pathway hits fit the later ATF4 + p53 finding (§40).

**Caveats.** The decoupler input statistics were dominated by confound genes (Il33, Prl, Fkbp5). 02/07 has no unit-of-replication test. The MANA file holds only 250 genes and should not be used for DE. These runs are exploratory and record no fixed hypotheses.

### 7. A claim-by-claim review confirmed 10 earlier observations and attributed 11 others to responder composition or confounds
*Code:* `docs/ANALYSIS_REVIEW.md` · *Results:* `results/region_time/pseudobulk_counts.csv`, `results/ol_gene_drivers/tables/driver_genes.csv`, `legacy/results/baysian_diff_condition.csv`

**Rationale.** On 2026-10-14 the leave-one-section-out classifiers showed that only RB4405, RB4401 and RB4627 respond (§10-11). Every earlier analysis had compared all six mtDSB mice with all six controls, so each earlier claim needed a re-check.

**Technical approach.** Earlier notebooks, legacy result folders and manuscript sections were read for stated conclusions. Each claim was checked per animal on the novae4 OL (DA-OL folded in), microglia, astrocyte and NSC pseudobulks (`results/region_time/pseudobulk_*.csv`, CPM per section) and against the analysis/03 driver genes. Claims were sorted into holds, needs rewording, and no testable conclusion.

**Results.** 10 earlier claims confirmed and attributed to responders, 11 explained by responder composition or confounds, and about 35 notebooks with no testable conclusions.
Confirmed:
- per-section OL markers of RB4405 / RB4627 are today's driver genes
- the 01/21 modules overlap about 80% with the 34 up-regulated driver genes
- the DA-OL cluster is 22-65% of OLs in responders and 0.3-2.3% elsewhere
- myelin genes down
- ISR fold changes are responder-only (Trib3 77-214 vs 1-7 CPM)
- Gpnmb 31-52 vs 2-7 CPM
- microglial DAM: Cd74 876-2705 vs 40-228 CPM, Cst7 394-1100 vs 12-89
- thalamus / fiber-tract hotspots
- the Bayesian model: all 22 driver genes agree in direction, 12 at FDR < 0.05
- Serpina3n 242-620 vs 34-59 CPM

Not holding:
- the P21-vs-P60 progression (1 vs 2 responders; RB4627 has the highest Trib3 and Gdf15)
- Hspa5 / Sod2 as part of the ISR (Sod2 is an age effect)
- weak group statistics (dilution)
- Sgk1 / Tsc22d3 / Fkbp5 / Zbtb16 / Dnajb1 / Hsph1 / Ddit4 (section-level stress, circadian and immediate-early confounds)
- Kif5a / Map2 / Camk2a (RB4403)
- Il33 (control RB4620, about 20x the others) and Xist (sex)
- olfactory-astrocyte ISR
- the NSC proliferation (does not track responders; 141-701 cells; section plane)
- cGAS-STING / NF-κB
- glycolysis (only Hk2 25-32 vs 10-17 CPM)
- Irf8 / Il1a / Hoxb3 / Fcgr2b

The Bayesian numbers were re-checked for this catalogue and are correct (§4).

**Reading.** The earlier biology was largely real, but carried by half of the mtDSB mice. Manuscript actions: replace the age-progression narrative, reword DA-OL to "3 of 6", drop the olfactory-astrocyte and NSC claims.

**Caveats.** Checks are per-animal CPM comparisons with three responders, not formal tests. Conflicts found while cataloguing:
- "Trib3 log2FC 5.1 ± 1.4" is cited as P21 but is the P60 estimate (§3).
- "Serpina3n β = 2.1" (01/30) is a group level, not an effect (§4).
- RB4401's DA-OL share is 21.4%, not ≥ 22% (§5).

### 8. A single section-level pseudobulk pipeline tests mtDSB vs control in every cell type x region x age stratum; most significant genes are spill-over or NSC artefacts, and OL hits are few
*Notebook:* analysis/01-region-time-treatment-figures.ipynb (parts 0-3, 5-7) · *Code:* `utils/region_time_de.py`, `utils/region_time_plots.py`, `scripts/run_region_time_de.py`, `tests/test_region_time_de.py` · *Results:* `results/region_time/summary_by_region_age.csv`, `results/region_time/summary_table.md`, `results/region_time/region_specific_hits.csv`, `results/region_time/run_info.json`, `results/region_time_cell_class_updated/summary_by_region_age.csv`, `results/region_time_MANA250/summary_by_region_age.csv` · *Figure:* `results/region_time/figures/overview_n_sig.png`, `results/region_time/figures/heatmap_panel_Oligodendrocytes.png`, `results/region_time/figures/volcano_Oligodendrocytes_All_regions_P60.png`, `results/region_time/figures/heatmap_panel_Microglia.png`, `results/region_time/figures/region_specific_Oligodendrocytes_P60.png`, `results/region_time/figures/heatmap_top_Neural_Stem_Cells.png`

**Rationale.** Three legacy notebooks (01/25, 01/26, 02/04) held pasted copies of the same region DE code. This run replaced them with one tested module that uses the section as the replicate (commit 1ac2c50, 2026-09-30).

**Technical approach.** Raw counts are summed per (cell_class, RBD_compartment_simplified, section), plus an "All regions" whole-section pseudobulk per cell type. Pseudobulks from < 20 cells are dropped, and "Unknown" is excluded. In each stratum, pyDESeq2 0.5.4 fits `~ condition` (genes with total count ≥ 10) and runs a Wald test of mtDSB vs control. A stratum is tested only if each condition has ≥ 2 sections. "Significant" means padj < 0.05, |log2FC| ≥ 0.5 and baseMean ≥ 20. A region-vs-rest test compares the treatment log2FC in a region with that in the whole section minus the region (z-test on the difference, BH-adjusted). The activation index is the number of up-genes x their mean log2FC. Inputs: the full novae4 h5ad (980,474 cells x 5,101 genes; 1,624 s). A second run (`_cell_class_updated`) re-summed the pseudobulk with Immature OL I-III merged. A third (`_MANA250`) used the 250-gene MANA file as an example. A synthetic test with planted effects checks the code.

**Results.** 1,264 pseudobulk samples were kept and 585 dropped, mostly Immature OL III (96), NSC (91), DA-OL (79) and MSN (72). 197 of 336 strata could be tested.

Significant genes across all strata (2,024 in total):

| Cell type | P21 | P60 |
|---|---|---|
| NSC | 0 | 1,325 (Ventricular system 671, All regions 654) |
| Neurons | 107 | 222 |
| Astrocytes | 37 | 55 |
| Endothelial | 27 | 31 |
| Oligodendrocytes | 21 | 25 |
| Microglia | 15 | 17 |
| DA-OL | 1 | 2 |

OL whole-section hits at P60 were Mt2 (padj 0.0047), Gfap (0.0076), Gpnmb (0.011) and Atf5 (0.048). Large fold changes missed significance: Trib3 4.00 (padj 0.14), Cst7 3.62 (0.51), Gdf15 2.76 (0.44), Cdkn1a 2.52 (0.53). The same held for microglial Cd74 (4.16, padj 0.15) and Cst7 (4.44, 0.15). The OL hypothalamus hits at P21 (10 up) are Fezf1, Pdyn, Gfap and Nptx2.

Region-vs-rest found 257 specific genes, 126 of them NSC and 74 neurons. The 7 OL hits are Pmch ×3, Pdyn, Nptx2, Gfap and Ntsr1. The merged immature-OL run had 36 hits, including Cnp, Cldn11 and Mag up in hypothalamus at P60. The MANA250 example had 284 hits.

**Reading.** At group level, OL ISR genes show large but non-significant fold changes, as expected when 3 of 6 cases do not respond (§10). The many neuronal and hypothalamic hits are neuropeptide spill-over. The NSC hits reflect ventricle coverage and section plane (README caveats; §39).

**Caveats.** Groups are 3 vs 3 per age, with age confounded with slide. The regions come from the disease-distorted map (§2). **Corrected 2026-10-21 (§49):** on the disease-free map the totals barely move (2,024 to 2,007).
Other conflicts and limits:
- The notebook header command names the novae3 h5ad, while the commit message and README say novae4. `run_info.json` does not record the input path.
- The FOCUS list includes "Oligodendrocyte Precursor Cells", a class that no longer exists, so OPCs appear only as Immature OL I.
- The MANA250 run is only an illustration (250 HVGs, old OPC labels).
- The large long tables are gitignored and must be regenerated with the script.

### 9. The treatment effect does not change between P21 and P60 for any OL stress gene; the age-by-condition hits are neuropeptide spill-over and NSC plane effects
*Notebook:* analysis/01-region-time-treatment-figures.ipynb (part 4; age scatter of part 2) · *Code:* `utils/region_time_de.py`, `utils/region_time_plots.py` · *Results:* `results/region_time/interaction_hits.csv`, `results/region_time/de_interaction_by_region_log.csv` · *Figure:* `results/region_time/figures/age_scatter_Oligodendrocytes.png`, `results/region_time/figures/interaction_Oligodendrocytes_All_regions.png`, `results/region_time/figures/age_scatter_Microglia.png`

**Rationale.** This run formally tests the legacy claim that the response is modest at P21 and strong at P60 (§3). Comparing hit counts between ages is not such a test.

**Technical approach.** One joint pyDESeq2 model per (cell type, region): `~ age + condition + age:condition`, with sections as replicates. It needs ≥ 2 sections in every age x condition cell. A Wald test on the interaction coefficient asks whether the treatment effect differs between P60 and P21. The same thresholds as §8 apply. As a cross-check, a z-test on the difference of the two stratum log2FCs (`de_age_contrast_by_region.csv`) gives a second version.

**Results.** 95 of 168 (cell type, region) pairs could be fitted. There were 358 interaction hits:

| Cell type | Interaction hits |
|---|---|
| Neurons | 145 |
| NSC | 110 |
| Endothelial | 31 |
| Microglia | 28 |
| Astrocytes | 23 |
| Oligodendrocytes | 8 |
| Immature OL I | 5 |

The 8 OL hits are Pmch, Hcrt, Gabra2, C1ql1 and Pcdh20 (thalamus, hypothalamus and whole section). These are neuronal / neuropeptide transcripts.

No OL stress gene interacts with age (whole-section padj 1.0 for all of them):

| Gene | Interaction log2FC ± SE | P21 effect | P60 effect |
|---|---|---|---|
| Trib3 | +1.31 ± 1.84 | 2.69 | 4.00 |
| Gdf15 | +0.05 ± 1.57 | 2.72 | 2.76 |
| Atf5 | +0.37 ± 0.45 | 0.57 | 0.95 |
| Cdkn1a | +0.85 ± 1.12 | 1.67 | 2.52 |

Mt2, Gstp1, Mog and Cldn11 likewise show no age interaction. The delta-z version gave 496 hits; among its 8 OL hits only thalamus Cd74 (delta 6.12, padj 7e-5) is not a neuropeptide. Cd74 is a microglial / antigen-presentation gene.

**Reading.** There is no evidence that the OL response depends on induction age. The earlier impression came from responder composition: one responder of three at P21, two of three at P60 (§10-11; README finding 3).

**Caveats.** Each age x condition cell holds 3 sections, so power is low. Age and slide cannot be separated. Responder dilution affects both ages unequally. A non-significant interaction is absence of evidence, not proof of equality. **Corrected 2026-10-14:** the age question is properly answered by cross-cohort classifier transfer (§11 / §16).

### 10. A predefined ISR gene set splits the six mtDSB mice into three strong responders and three that look like controls
*Notebook:* analysis/02-ol-response-and-classifiers.ipynb (part 1) · *Code:* `utils/ol_response.py`, `scripts/run_region_time_de.py` · *Results:* `results/sample_classifier/tables/sections_with_responder_calls.csv`, `results/sample_classifier/tables/isr_per_gene_whole_section.csv`, `results/sample_classifier/tables/isr_group_tests.csv`, `results/sample_classifier/tables/isr_scores_by_section_region.csv` · *Figure:* `results/sample_classifier/figures/F1_isr_responders.png`

**Rationale.** Group-level DE and the earlier models (§3, §4, §8, §9) gave weak, age-dependent mtDSB effects, and the leave-one-section-out classifiers (§12) could not predict condition at all. The question was whether the mtDSB group is uniformly shifted or heterogeneous. A readout that uses neither the clustering nor any model was needed, so the call could not be circular with the DA-OL cluster.

**Technical approach.** Section-level oligodendrocyte pseudobulk from the region x age pipeline (§8; `results/region_time/pseudobulk_*.csv`). `DA-Oligodendrocytes` is folded into `Oligodendrocytes` ("general labels", `olr.general_pseudobulk`). ISR score = share of OL transcripts (CPM) from 10 canonical ATF4 / mitochondrial ISR targets fixed from the literature before inspecting per-section values (`ISR_GENES`: Atf4, Atf5, Ddit3, Trib3, Gdf15, Sesn2, Fgf21, Psat1, Cth, Eif4ebp1; Mt1/Mt2, Vegfa, Hspd1 deliberately left out). Rule (`responder_calls`, margin 0.2): an mtDSB section is a responder if its whole-section score is > 1.2x the highest control. Statistics on the raw prespecified score without any responder selection: one-sided Mann-Whitney mtDSB > control, one-sided F-test of variance ratio and Brown-Forsythe (Levene, median) on log ISR CPM. Unit: section = mouse (n = 6 vs 6).

**Results.** Responders RB4405 (4031 CPM, 1.57x highest control), RB4627 (3674, 1.43x), RB4401 (3477, 1.35x). Controls 2039-2568 CPM; non-responders RB4403 2347, RB4630 2269, RB4653 2059 (0.80-0.91x). Group tests: Mann-Whitney U = 25, p = 0.155; variance ratio F = 9.58, p = 0.013; Brown-Forsythe W = 38.3, p = 1.0e-4. Per gene (whole-section OL CPM), responder vs control ranges do not overlap for Atf4 (1255-1801 vs 637-950), Atf5 (503-741 vs 136-172), Ddit3 (256-377 vs 115-151), Trib3 (77-214 vs 0.9-7.0), Sesn2 (133-208 vs 36-65), Eif4ebp1 (13-27 vs 0.9-2.5) and also Gdf15 (12.7-65.4 vs 0.4-3.0). Psat1 goes the other way in two responders (RB4405 671, RB4627 782 vs controls 1061-1325).

**Reading.** The mtDSB group is heterogeneous, not uniformly shifted: about half of the animals carry a strong OL stress response and half none.

**Caveats.** The 1.2x margin is a description of 12 animals, not a test; the formal evidence is the variance test. FINDINGS.md says non-responders "sit inside the control range on every gene", but RB4403 is slightly above the control maximum for Atf5 (179.9 vs 171.6), Ddit3 (164.4 vs 150.6) and Cth (19.8 vs 19.1). FINDINGS counts six non-overlapping genes; the table shows seven (Gdf15 also separates). The score is dominated in absolute CPM by Atf4 and Psat1, and Psat1 falls in two responders, so the summed score understates the per-gene effects. The panel has no transgene probes, so induction cannot be checked. Age cohort = Xenium slide.

**Corrected 2026-10-15:** analysis/07 (§25) replaced the yes/no rule with P(responder): the three responders > 0.99; RB4630 0.17 (borderline), RB4653 and RB4403 < 0.01. **Corrected 2026-10-21:** analysis/15 (§45) showed responders also separate on microglia and on OL genes not used for the call.

**Next.** §11 checks the call with independent readouts; §25 attaches a probability to each call.

### 11. Three readouts built from different inputs single out the same three mice
*Notebook:* analysis/02-ol-response-and-classifiers.ipynb (part 2) · *Code:* `utils/ol_response.py`, `utils/sample_classifier.py`, `scripts/run_sample_classifier.py` · *Results:* `results/sample_classifier/tables/three_readouts_by_section.csv`, `results/sample_classifier/tables/three_readouts_spearman.csv`, `results/sample_classifier/tables/ol_classifier_scores_by_section_region.csv`, `results/sample_classifier/tables/da_cluster_share_by_section_region.csv` · *Figure:* `results/sample_classifier/figures/F2_three_readouts.png`

**Rationale.** A single prespecified gene set (§10) could be fooled by a handful of genes. If the responder split is real, readouts with different inputs should rank the same animals at the top.

**Technical approach.** Three per-section readouts. (1) ISR score from §10. (2) Classifier score: mean out-of-fold p(mtDSB) of a section's OLs from a pooled single-cell L2 logistic regression (C = 0.1, balanced classes, up to 2000 cells per training section, log-normalised and standardised genes, 5101 genes, general labels), trained leave-one-section-out on the other 11 sections' condition labels; also the fraction of OLs with p > 0.8. (3) Cluster share: fraction of OLs annotated `DA-Oligodendrocytes` (descriptive; depends on the clustering). Agreement: Spearman correlation across the 12 sections and overlap of each readout's top 3; chance that the top 3 of 12 are all mtDSB under random labels = 20/220.

**Results.** Top 3 by every readout = RB4405, RB4627, RB4401. Classifier mean p: responders 0.61-0.81, controls 0.38-0.60, non-responders 0.28-0.37. Fraction p > 0.8: responders 0.46-0.71, controls 0.18-0.43, non-responders 0.13-0.19. DA-OL cluster share: responders 0.215-0.647 vs 0.003-0.023 in all other sections. Spearman: ISR vs classifier 0.71, ISR vs cluster 0.56, classifier vs cluster 0.59. Chance for a single readout's top 3 all being mtDSB = 0.091.

**Reading.** The DA-OL cluster is the responder state, and a classifier that never saw the cluster finds the same animals.

**Caveats.** The readouts are not fully independent: all come from OL expression. The classifier is poorly calibrated (control OLs also score high; control means up to 0.60). The 0.091 chance figure applies to one readout and is not a combined test. Responders differ in severity (RB4401 mildest on the cluster share).

### 12. Plain mtDSB-vs-control classifiers fail because half the mtDSB labels are non-responders; sex, genotype, plane and slide do not explain the split
*Notebook:* analysis/02-ol-response-and-classifiers.ipynb (parts 3, 6) · *Code:* `utils/sample_classifier.py`, `scripts/run_sample_classifier.py` · *Results:* `results/sample_classifier/tables/pseudobulk_auc_whole_section.csv`, `results/sample_classifier/tables/pseudobulk_permp_whole_section.csv`, `results/sample_classifier/tables/ol_cell_model_section_predictions.csv`, `results/sample_classifier/tables/technical_checks.csv`, `results/sample_classifier/pseudobulk_general/metrics.csv` · *Figure:* `results/sample_classifier/figures/F3_classifiers.png`, `results/sample_classifier/figures/F6_technical_checks.png`

**Rationale.** The classifiers were first built to ask how separable age and condition are, and where the signal lives. Condition turned out not to be predictable. This run explains why, tests whether the responder signature is the same at both ages, and rules out technical explanations for non-response.

**Technical approach.** Pseudobulk level: one L2 logistic regression (kernel form, C = 10, genes centred not standardised, min 10 CPM) per (cell type, region) stratum on centred log-CPM, plus a composition stratum. Leave-one-section-out (LOSO) AUC with sample-level label permutations kept balanced (condition permuted within age and vice versa; 1000 permutations or exact). Tasks: age, condition, condition within each cohort, and cross-cohort transfer (train P21 -> test P60 and the reverse). Single-cell OL models (§11 settings) give section-level transfer predictions. Technical checks per section: Xist CPM vs recorded sex, % of cells in fibre tracts, cell count, genotype label.

**Results.** Age: 52 of 101 pseudobulk strata at BH-FDR < 0.1; single-cell models reach AUC 1.0 for 10 of 12 cell types. Condition (pooled): 0 of 101 strata at permutation p < 0.05 (minimum 0.05); OL pseudobulk condition AUC 0.44. Cross-age transfer of OL single-cell models: trained on P60, RB4627 is the top P21 section (0.65 vs 0.29-0.50); trained on P21, RB4405 and RB4401 are the top P60 sections (0.88, 0.85 vs 0.64-0.74). Sex labels match Xist in all 12 sections. Responders are 2 M + 1 F, as are non-responders. All six mtDSB mice are PlpTA:mtPst1; control RB4620 is mtPst1-only. Fibre-tract share: responders 8.8-10.8%, others 5.8-9.5%. RB4403 has 135,378 cells (others 69,913-88,718). Responders and non-responders share each slide (P60: RB4401, RB4405 vs RB4403; P21: RB4627 vs RB4630, RB4653).

**Reading.** Age is easy to predict and condition is not, because half of the "mtDSB" labels belong to animals that look like controls. The responder signature transfers across cohorts, so the earlier "P21 modest / P60 strong" reading (§3) reflects 1 vs 2 responders.

**Caveats.** Age cohort = Xenium slide (P60 on 0060539, P21 on 0060541), so the age signal cannot be separated from slide. The within-cohort and transfer tasks cannot reach p < 0.05 by construction (3 vs 3; the observed minimum p is 0.05-0.1), so "0 strata significant" is only informative for the pooled task. LOSO AUC well below 0.5 is a known artefact without signal. The 2/3 vs 1/3 responder split by age is uninformative (Fisher p = 1.0, stated in FINDINGS.md). Two facts in FINDINGS.md are not printed by the notebook: RB4403's 128,804 interior-stain cells (segmentation) and the mirror-image mounting of controls vs mtDSB sections. Sections being mirror-mounted is a possible confound for weak all-mtDSB signals (§14).

### 13. In responders every region is affected, most strongly thalamus and fibre tracts
*Notebook:* analysis/02-ol-response-and-classifiers.ipynb (parts 4, 5, 8, 9) · *Code:* `utils/ol_response.py`, `utils/spatial_orientation.py`, `scripts/export_ol_gene_scores.py`, `scripts/export_cell_table.py` · *Results:* `results/sample_classifier/tables/region_summary_by_group.csv`, `results/sample_classifier/tables/region_rank_in_responders.csv`, `results/sample_classifier/tables/ol_frac_above_control_p90.csv`, `results/sample_classifier/tables/ol_score_percentile_by_section.csv`, `results/sample_classifier/tables/ol_isr_cell_score_by_section.csv` · *Figure:* `results/sample_classifier/figures/F4_regions.png`, `results/sample_classifier/figures/F5_spatial_ol_scores.png`, `results/sample_classifier/figures/F8_ol_score_maps.png`, `results/sample_classifier/figures/F9_ol_isr_maps.png`

**Rationale.** Earlier work (§2, §8) and IHC had suggested hypothalamus / thalamus and white matter. With the responders known (§10), the question was whether the response is focal or brain-wide and which regions lead.

**Technical approach.** Region labels = original read-based compartments (`RBD_compartment_simplified`, §2). Per (section, region), each of the three readouts of §11 is divided by the mean of the six controls in that region. Regions with < 50 OLs are blanked. Regions are ranked within each responder (1 = highest) and the ranks averaged over the three readouts. Maps: (a) each OL's out-of-fold classifier p(mtDSB), flagged if above the control 90th percentile (0.979) or expressed as a percentile of pooled control OLs; (b) per-OL ISR score log1p(1e4 x ISR counts / total), averaged over the 15 nearest OLs and expressed as a percentile of smoothed control OLs on the same slide. All sections oriented dorsal up, midline right.

**Results.** In responders every region is up: ISR 1.23-1.88x controls, classifier high-score fraction 1.39-2.93x, DA cluster share 21-65% (controls 0.4-4.2%). Mean rank over readouts: fibre tracts 2.7, thalamus 2.7, vasculature 4.1, striatum 4.6, hypothalamus 4.8, cortex 5.0, ventricular 5.6, pallidum 6.7, meningeal 7.2, hippocampus 8.7, olfactory 8.9. OLs above the control 90th percentile of classifier p: responders 24-40%, controls 5-16%, non-responders 4-6%. Classifier percentile median: responders 65-85, controls 38-63, non-responders 23-35. Smoothed ISR above the control 90th percentile: responders 31-55% (RB4405 0.549, RB4627 0.537, RB4401 0.305), controls 7-13%, non-responders 8-14%.

**Reading.** The response is brain-wide in responders, with thalamus and white-matter tracts strongest; hypothalamus is affected but mid-ranked.

**Caveats.** Regional differences are smaller than the responder vs non-responder difference and rest on three animals. "Vasculature" ranks third on the classifier but holds few OLs. The original region map was itself distorted by disease (pallidum absorbed into fibre tracts in responders). Smoothing over 15 OLs is for display and inflates spatial coherence. The classifier probabilities are poorly calibrated, so they are shown as control percentiles.

**Corrected 2026-10-21:** analysis/14 (§43, §44) re-ran the regions on a disease-free map: thalamus (28-36% of OLs above the control 90th percentile) and fibre tracts (19-33%) stay at the top and pallidum joins them (20-35%, vs 6.7 mean rank here, where responders held only a few hundred pallidum OLs). Fibre-tract statements should read "white matter including the pallidum".

### 14. Outside OLs, single-cell classifiers give only weak, untested cross-cohort hints in astrocytes, endothelium and microglia
*Notebook:* analysis/02-ol-response-and-classifiers.ipynb (part 7) · *Code:* `utils/sample_classifier.py`, `scripts/run_sample_classifier.py` · *Results:* `results/sample_classifier/tables/cells_auc_whole_section.csv`, `results/sample_classifier/tables/cells_auc_glia_vascular_by_region.csv`, `results/sample_classifier/cells_general/predictions.csv` · *Figure:* `results/sample_classifier/figures/F7_cells_auc_all_cell_types.png`

**Rationale.** Once OLs carried the signal, the question was whether other cell types also distinguish mtDSB from control, which would point to a non-cell-autonomous response.

**Technical approach.** Same LOSO single-cell logistic regression as §11, one model per cell type (12 general cell types), tasks age, condition, condition within cohort and cross-cohort transfer. Section-level prediction = mean cell p per section and region. AUC only; no permutation test at single-cell level.

**Results.** Pooled condition AUC (whole section): astrocytes 0.50, endothelial 0.50, microglia 0.44, OLs 0.50; within-cohort AUCs 0.0-0.33 for most types. Cross-cohort transfer: astrocytes 0.78 / 0.78 (21->60 / 60->21), endothelial 0.78 / 1.0, microglia 0.67 / 0.89, OLs 0.78 / 0.89. Endothelial 60->21 ranks all three P21 mtDSB sections (0.44-0.51) above all three P21 controls (0.37-0.43), including the OL non-responders. In microglia and astrocytes the P21->P60 model ranks the two P60 responders highest (microglia RB4405 0.77, RB4401 0.72; astrocytes 0.72, 0.68). Per region, astrocyte and endothelial transfer AUCs reach 0.89-1.0 in some strata (e.g. endothelial vasculature 0.89; astrocyte meningeal and pallidum 60->21 1.0).

**Reading.** Listed as hypotheses only in FINDINGS.md: microglia echo the responder pattern; endothelial and astrocyte signals may follow condition rather than response.

**Caveats.** Margins are tiny (endothelial 0.44 vs 0.43), there is no permutation test, many cell types and regions were tested, and 3 vs 3 transfer has a floor of 1/20. Control and mtDSB sections were mounted mirror-image on both slides, so weak all-mtDSB signals could be technical.

**Corrected 2026-10-21:** analysis/12 (§37) found, after a spill-over model, real own changes in microglia (DAM, MHC-I), astrocytes (Serpina3n, C4b, Gfap, small IFN arm) and some endothelial changes in responders vs controls; most neuronal and immature-OL hits were neighbours' RNA.

### 15. Within each responder section, 34 genes rise and 9 myelin genes fall with an OL's stress score
*Notebook:* analysis/03-ol-stress-genes.ipynb (parts 1, 3, 4, 5) · *Code:* `scripts/run_ol_gene_drivers.py`, `utils/ol_response.py` · *Results:* `results/ol_gene_drivers/tables/driver_genes.csv`, `results/ol_gene_drivers/tables/all_genes_association_and_weights.csv`, `results/ol_gene_drivers/tables/driver_genes_log2fc_by_section.csv`, `results/ol_gene_drivers/tables/driver_programme_by_region.csv`, `results/ol_gene_drivers/tables/glucocorticoid_genes_by_section.csv` · *Figure:* `results/ol_gene_drivers/figures/G1_within_section_association.png`, `results/ol_gene_drivers/figures/G3_driver_genes_by_section.png`, `results/ol_gene_drivers/figures/G4_driver_programme_by_region.png`

**Rationale.** The OL classifier (§11) separates responders but its weights mix real biology with section-level confounds (batch, plane, time of sacrifice). Genes that separate stressed from unstressed OLs inside the same tissue section are immune to those.

**Technical approach.** For every OL (159,301 scored cells), Pearson correlation between each gene's log-normalised expression and the logit of the cell's out-of-fold p(mtDSB) (§11), computed within each section (`gene_score_association`). Driver rule fixed before inspecting genes: r > 0.1 in each of the three responders and mean r < 0.1 in both controls and non-responders (mirror rule for down genes). Programmes annotated manually. Per-section effect sizes: whole-section OL pseudobulk log2 fold change vs the six-control mean; per region, the mean log2 FC of up-drivers and of down-drivers vs that region's controls. Checks: glucocorticoid-response genes, neuron vs OL CPM ratio in controls (flag > 4x), Xist.

**Results.** 34 up, 9 down. Up: ATF4/ISR (Trib3 r 0.29, Atf5, Atf4, Ddit3, Cebpg, Sesn2, Yars, Vegfa), mitochondrial / protein quality control (Mt2, Hspa9, Gstp1, Sqstm1, Nfe2l1, Lamp2), cell cycle / survival (Cdkn1a, Ccnd1, Bcl2l1, Myc), immune / reactive (Cd40, Serpina3n, H2-D1, Adamts1, B2m), lineage (Olig2, Adgrg1), cytoskeleton (Flna, Mpzl1, Ptpra, Coro1c, Msn), lipid (Abca1, Slc6a8), other (Pkd2l1, Zfp281). Down: Cldn11 (r -0.28), Mog, Mal, Stmn4, Cnp, Aplp1, Ugt8a, Tppp, Mapk8ip1. Six of the ten prespecified ISR genes are recovered. Per animal (OL CPM): Cd40 1-4 -> 49-165, Cdkn1a 6-31 -> 147-235, Serpina3n 34-59 -> 242-620, Mog 3741-6376 -> 1570-3778. Region means in responders: up-programme thalamus 1.64, fibre tracts 1.60, hypothalamus 1.45 ... olfactory 0.83 (controls -0.58 to -0.05); myelin down -0.77 in thalamus and fibre tracts. Glucocorticoid genes (Sgk1 weight 0.43, Zbtb16 0.30, Plin4 0.26, Fkbp5 0.20) have r_responder -0.03 to 0.10 and Sgk1 is highest in a control (RB4282 4450 CPM); none is a driver. No driver is neuron-enriched (max 1.8x, Nfe2l1).

**Reading.** Stressed OLs switch on an ATF4 / p21 / immune-visible programme and lose mature-OL identity. Glucocorticoid genes track a section-level state and are not mtDSB drivers.

**Caveats.** The association is with the classifier score, which was trained on condition labels, so drivers are genes that co-vary with that score. Gene ranks rest on 3 animals and the programme labels are manual. Mog ranges touch (RB4401 3778 vs non-responder RB4403 3741). FINDINGS.md says the up-programme is "highest in the thalamus in all three responders"; G4 shows that is true only for the group mean (one responder is higher in fibre tracts, about 2.0, and in ventricular system, about 2.2, than in its thalamus). Regions use the original map (pallidum missing in responders; G4 has only two responder points for pallidum). G4 panels are lettered D/E although the figure stands alone.

**Corrected 2026-10-21:** analysis/13 (§40, §42) refined the programme (ATF4 + p53, no cGAS-STING; myelin RNA lost per OL without OL loss); analysis/12 (§37) confirmed OLs keep Atf5, Trib3, Cdkn1a after spill-over correction.

### 16. Classifiers trained on two disjoint sets of mice share their top genes far above chance
*Notebook:* analysis/03-ol-stress-genes.ipynb (part 2) · *Code:* `scripts/run_ol_gene_drivers.py`, `utils/ol_response.py`, `utils/sample_classifier.py` · *Results:* `results/ol_gene_drivers/tables/cohort_weight_replication.csv`, `results/ol_gene_drivers/coefficients.csv.gz` · *Figure:* `results/ol_gene_drivers/figures/G2_cohort_replication.png`

**Rationale.** With n = 3 responders, a replication with no animals in common is the strongest available check that the gene programme (§15) is not one mouse's quirk.

**Technical approach.** OL single-cell logistic regressions (C = 0.1, 2000 cells per section, standardised genes) trained only on the six P21 mice (condition 21->60) or only on the six P60 mice (condition 60->21). Compare their gene weights: genome-wide Spearman, and overlap of the top-k up and down genes (k = 25, 50, 100) against a hypergeometric null over 5101 genes.

**Results.** Genome-wide Spearman 0.035. Top-50 up overlap 10 (0.49 expected, p = 2.4e-11): A2m, Atf5, Bcl2l1, Mc3r, Meg3, Mt2, Plin4, Rnf144a, Trib3, Zbtb16. Top-50 down 4 (p = 0.0014): Cldn11, Il33, Opalin, Vstm2a. Top-25 up 5 (p = 9.2e-8); top-100 up 15 (p = 5.4e-10); top-25 down 1 (p = 0.12); top-100 down 6 (p = 0.013).

**Reading.** The strongest stress genes (Trib3, Atf5, Mt2, Bcl2l1) and myelin genes (Cldn11, Opalin) replicate across disjoint animals.

**Caveats.** Two of the ten shared top-50 up genes (Zbtb16, Plin4), and Fkbp5 in the top-100, are glucocorticoid-response genes that §15 excludes as section-level confounds; Il33 among the shared down genes was earlier flagged as a one-control confound (README finding 8). So part of the overlap reflects shared section-level state, not mtDSB alone. Each cohort model rests on 3 vs 3 animals including one (P21) or two (P60) responders.

### 17. Astrocytes and microglia next to stressed OLs are more activated, only in responders; endothelium is not
*Notebook:* analysis/04-neighbourhood-and-staging.ipynb (part 1) · *Code:* `utils/neighbourhood.py`, `utils/gene_sets.py`, `scripts/export_cell_gene_scores.py` · *Results:* `results/neighbourhood_staging/tables/coupling_summary.csv`, `results/neighbourhood_staging/tables/coupling_per_section.csv` · *Figure:* `results/neighbourhood_staging/figures/N1_neighbour_coupling.png`

**Rationale.** §14 hinted at glial signals and the legacy work (§5) reported a microglial DAM response. The question was whether glia respond locally to stressed OLs, measured inside sections so slide, plane and sacrifice time cannot create it.

**Technical approach.** Per-cell scores for prespecified literature gene sets (DAM, homeostatic, reactive astrocyte, endothelial activation, MHC-I, type-I IFN; `utils/gene_sets.py`). For each microglia, astrocyte or endothelial cell, local OL stress = mean ISR score of OLs within 50 um (at least 3 OLs). Both local stress and the readout are centred within (section, region); Pearson r per section (at least 50 cells). Responders vs controls by an exact one-sided permutation over all 84 relabellings of 3 + 6 mice (floor p = 0.012; homeostatic tested for lower). Repeated at 30 and 100 um. Controls are the null for baseline co-patterning.

**Results.** At 50 um, mean r responders vs controls: reactive astrocytes 0.056 vs -0.039 (p = 0.012); microglial MHC-I 0.031 vs 0.001 (p = 0.012); microglial DAM 0.082 vs 0.034 (p = 0.048); homeostatic microglia -0.048 vs -0.019 (p = 0.13); endothelial activation -0.015 vs -0.002 (p = 0.87); IFN in all three types p = 0.14-0.90. Astrocyte reactivity p = 0.012 at 30, 50 and 100 um; DAM p = 0.083 / 0.048 / 0.012 and MHC-I 0.071 / 0.012 / 0.012 at 30 / 50 / 100 um.

**Reading.** Glia couple to local OL stress only where OLs are stressed. Correlations are small because single-cell scores are noisy.

**Caveats.** The readout gene sets overlap the OL driver genes: MHC-I includes H2-D1 and B2m, and reactive astrocyte includes Serpina3n and H2-D1, all up in stressed OLs (§15). This cell-level coupling has no spill-over control, so part of it could be OL RNA picked up by neighbouring glia. Ten readouts x three radii are tested with no multiplicity correction. n = 3 responders. Region centring uses the original map. Co-location, not causation.

**Corrected 2026-10-21:** analysis/11 (§36) showed the extra microglia are spread through white matter rather than piled on the most stressed OLs (top-1% OLs have 5-25% more microglia within 50 um, as in controls). analysis/12 (§37) found microglial MHC-I and astrocyte Serpina3n / H2-D1 changes to be own (not spill-over) at the pseudobulk level, and analysis/15 (§46) reported microglial DAM and astrocyte reactivity holding with any one responder left out.

### 18. Stressed OLs come first, myelin loss follows: the share of stressed OLs that lost myelin rises with severity
*Notebook:* analysis/04-neighbourhood-and-staging.ipynb (part 2) · *Code:* `utils/ol_staging.py`, `utils/gene_sets.py`, `scripts/export_cell_gene_scores.py` · *Results:* `results/neighbourhood_staging/tables/ol_isr_myelin_states.csv`, `results/neighbourhood_staging/tables/myelin_by_isr_bin.csv` · *Figure:* `results/neighbourhood_staging/figures/N2_isr_myelin_staging.png`

**Rationale.** §15 showed stress genes up and myelin genes down in the same cells. With one time point (42 days), the only handle on order is a cross-sectional comparison of animals of different severity.

**Technical approach.** Per OL, prespecified ISR score (10 genes) and myelin score (13 genes: Mbp, Mog, Mag, Cldn11, Cnp, Mal, Opalin, Ugt8a, Aspa, Ermn, Myrf, Enpp6, Car2). Both smoothed over the 15 nearest OLs and converted to percentiles of control OLs of the same slide and region. ISR-high = > 90th percentile; myelin-low = < 10th. Per section: fractions in each state, P(myelin-low | ISR-high) and P(ISR-high | myelin-low); mean myelin percentile per ISR percentile bin. Responders vs controls by exact permutation (84 relabellings).

**Results.** ISR-high: RB4401 32%, RB4627 60%, RB4405 58% (controls 7-13%). Myelin-low: RB4401 21%, RB4627 65%, RB4405 86% (controls 4-15%). P(myelin-low | ISR-high): RB4401 0.40, RB4627 0.73, RB4405 0.90 vs controls 0.09-0.25 (difference 0.49, p = 0.012). ISR-high with normal myelin: responders - controls 0.054, p = 0.083 (RB4405 0.055, inside the control range). Myelin percentile falls with ISR bin in responders.

**Reading.** In the mildest responder most stressed OLs still have normal myelin genes; with severity a growing share also lose them, and in the most severe mouse myelin loss extends beyond stressed cells. This is consistent with stress preceding myelin loss.

**Caveats.** Not a time course; the severity ordering rests on three animals. Spatial smoothing makes neighbouring OLs share scores and couples the two percentiles regionally. Percentiles are per (slide, original region), and the original map put responder pallidum into fibre tracts (§43), which can shift myelin percentiles. Non-responder RB4403 has 27% myelin-low OLs (likely over-segmentation). The "ISR-high, myelin normal" excess is not significant.

**Corrected 2026-10-21:** analysis/13 (§42) counted myelin RNA without segmentation and confirmed Mog / Opalin fall 40-55% per OL while Sox10 and Olig2 stay, so myelin genes are lost, not OLs.

### 19. A microglia's own expression predicts how close it sits to a stressed OL, in responders only
*Notebook:* analysis/04-neighbourhood-and-staging.ipynb (part 3) · *Code:* `scripts/run_proximity_model.py`, `utils/neighbourhood.py` · *Results:* `results/proximity_model/metrics.csv`, `results/proximity_model/coefficients.csv.gz`, `results/neighbourhood_staging/tables/proximity_model_tests.csv`, `results/neighbourhood_staging/tables/proximity_model_summary.csv` · *Figure:* `results/neighbourhood_staging/figures/N3_proximity_model.png`

**Rationale.** The coupling in §17 uses small prespecified gene sets that overlap OL genes. An unbiased model with OL-related genes removed tests whether glia near stressed OLs have a distinct expression state of their own.

**Technical approach.** Target per microglia, astrocyte or endothelial cell: log1p distance to the nearest stressed OL (smoothed ISR > 90th percentile of same-slide, same-region control OLs), residualised within each section on region (one-hot), log distance to any OL and log OL count within 50 um. Features: the cell's log-normalised expression, standardised on training cells; "strict" set removes OL-enriched genes (> 2x CPM in OLs vs that type in controls: 366 / 282 / 401 genes), the 43 driver genes and the ISR genes. Ridge regression (alpha by efficient LOO, 1e1-1e5), up to 3000 cells per section. Held-out Pearson r per mouse: responders (train 2, test 1), controls (train 5, test 1; the null), responder model applied to controls and non-responders. Responders vs controls by exact permutation (84).

**Results.** Microglia, strict: responders 0.135-0.196 (mean 0.168) vs controls 0.076-0.126 (0.102), every responder above every control, p = 0.012; responder model on controls 0.074, on non-responders 0.047. All genes: 0.168 vs 0.111, p = 0.024. Astrocytes strict 0.101 vs 0.087 (p = 0.024; all genes p = 0.012). Endothelial strict 0.081 vs 0.094 (p = 0.77). Top microglial "closer" weights (strict): Scn4b, Cd68, Cd83, Hcn1, Speg, Pde10a, Adora1, Trem2, S1pr1.

**Reading.** Microglia near stressed OLs carry an activation / phagocytic signature (Cd68, Cd83, Trem2) that controls do not show.

**Caveats.** Controls are not at zero (r about 0.1, read through residual sub-regional and neuronal transcripts), so the evidence is the difference. Several of the largest microglial weights are neuronal genes (Scn4b, Hcn1, Pde10a, Adora1, Arc, Bsn), so the signature is not purely microglial. Responder models train on 2 mice and control models on 5. Effects are modest. Stressed OLs are defined on the original region map.

**Corrected 2026-10-21:** analysis/11 (§34-§36): responders have 1.5-3x the microglia density in thalamus and fibre tracts, spread through the region rather than clustered on stressed OLs.

### 20. With half the mtDSB mice responding, a classic group comparison is underpowered; an induction check makes 5-6 mice per group enough
*Notebook:* analysis/05-design-power.ipynb · *Code:* (inline in the notebook) · *Results:* `results/design_power/power_by_design.csv` · *Figure:* `results/design_power/figures/D1_power_by_design.png`

**Rationale.** §10 showed group means do not differ while variance does. The next cohort needs a design and a test that work when only some mice respond.

**Technical approach.** Inputs: per-mouse log whole-section OL ISR CPM (§10), minus each cohort's control mean (slide baseline). Parameters from 12 sections: controls mean 0, SD 0.026 (n = 6); responders mean 0.452, SD 0.112 (n = 3); non-responders -0.009, SD 0.068. Generative model: controls ~ Normal(control); each mtDSB mouse responds with probability p_resp (0.5 or 0.75) and is then drawn from the responder normal, otherwise like a control. Tests (one-sided, alpha 0.05): Mann-Whitney on means; Brown-Forsythe variance test (with larger mtDSB variance); responder count above the control maximum vs a 200-permutation null; Mann-Whitney after a perfect induction QC that drops non-responders. 600 simulations per cell, n = 3-12 per group.

**Results.** Responder shift = 17.1 control SDs. At p_resp 0.5, power at 6 / 12 mice per group: mean 0.45 / 0.69; variance 0.47 / 0.93; responder count 0.48 / 0.91; QC mean 0.88 / 1.00 (0.80 already at 5). Variance and count reach 0.82 at 10. At p_resp 0.75: mean 0.79 at 6; count 0.85 at 6. Recommendations: induction QC (mtPstI qPCR / mtDNA damage) per mouse; analyse and report per mouse; put both cohorts on every slide.

**Reading.** Without induction QC, per-mouse tests need about 10 mice per group; with QC that removes failed inductions, 5-6 suffice.

**Caveats.** Parameters come from 3 responders and 6 controls. The control SD is estimated after subtracting each cohort's mean from 3 mice, which shrinks it, so the 17-SD shift is optimistic (analysis/07 gives 11-19 control SDs). Non-responders are simulated like controls although their observed SD is 2.6x larger. QC power assumes non-response = failed induction and a perfect check. 600 simulations give Monte Carlo error of about 0.02. The variance-test power is non-monotone at p_resp 0.5 (0.41 at n = 4, 0.18 at n = 5), probably a small-sample artefact of the Levene test, not explained in the notebook. The responder-count test uses margin 0 although the markdown mentions a margin.

**Next.** analysis/07 (§25) adds a probabilistic responder rule for the next cohort (ISR >= 1.17x the slide's control mean gives P(responder) >= 0.9).

### 21. OL stress does not rise towards the ventricles or meninges: it is not border-driven
*Notebook:* analysis/06-borders-patches-repair-neurons.ipynb (part A) · *Code:* `utils/spatial_patterns.py`, `utils/ol_staging.py` · *Results:* `results/spatial_patterns/tables/border_coupling_summary.csv`, `results/spatial_patterns/tables/border_coupling_per_section.csv` · *Figure:* `results/spatial_patterns/figures/S1_border_distances.png`

**Rationale.** The collaborators' IHC found immune cells at the ventral meninges, around the ventricle and in the internal capsule / thalamus. If OL stress started at the CSF or meningeal surface (outside-in), responder OLs should be more stressed close to these borders. Follows the regional finding of §13 and the neighbour coupling of §17. Committed 2026-10-14 (cefdd05).

**Technical approach.** Input: per-cell gene-set scores and centroids (`results/cell_table/cell_gene_scores.csv.gz`, `cells.csv.gz`; cells with region "Unknown" dropped), old read-based region map (§2). OLs only (`cell_type == "Oligodendrocytes"`). Each OL's ISR score (prespecified gene set, `utils/gene_sets.py` / `utils/ol_response.py`) was smoothed over its 15 nearest OLs (`ol_staging.smooth_scores`, k = 15) and converted to a percentile among control OLs of the same slide and region. `border_distances`: distance of every cell to the nearest cell in region "Ventricular system", region "Meningeal border/glia limitans" and the nearest endothelial cell (vessel proxy), same section, um. `distance_coupling`: OLs inside the ventricle and meningeal regions excluded; smoothed ISR and log(1 + distance) centred within (section, region), then rank correlation per section. Unit of inference: the section (= mouse). Responders (3) vs controls (6) by exact permutation over all 84 relabellings, two-sided (floor p = 0.012). Panel D of S1 plots mean control-percentile ISR by distance-to-ventricle bin (0-100, 100-250, 250-500, 500-1000, 1000-5000 um).

**Results.** Mean per-section r (stress vs log distance; positive = more stress farther away):
- ventricle: controls 0.042, responders 0.109, p = 0.095
- meninges: controls 0.021, responders 0.087, p = 0.25
- nearest vessel: controls -0.045, responders -0.076, p = 0.083

The notebook concludes no concentration at ventricles or meninges; responder OLs are, if anything, more stressed farther from these borders. OLs near vessels are slightly more stressed in all groups (baseline anatomy). The README states responders sit at the 75th-80th control percentile at every distance from the ventricle (S1 D); that number is read off the figure and is not in a saved table.

**Reading.** The response is a regional shift spread through the tissue, not an outside-in front from CSF or meninges.

**Caveats.** Negative result with 3 responders. Distances depend on the old read-based region map; the ventricle and meningeal classes were later shown to be imperfect (§43). This analysis was not re-run on the disease-free map. The vessel proxy is the nearest endothelial-labelled cell. The smoothing (k = 15) mixes neighbouring OLs into each score.

### 22. OL stress is diffuse within regions, not focal patches
*Notebook:* analysis/06-borders-patches-repair-neurons.ipynb (part B) · *Code:* `utils/spatial_patterns.py` · *Results:* `results/spatial_patterns/tables/isr_correlogram_summary.csv`, `results/spatial_patterns/tables/isr_correlogram.csv` · *Figure:* `results/spatial_patterns/figures/S2_isr_correlogram.png`

**Rationale.** If mtDSB damage produced lesions, stressed OLs would cluster into patches. A first patch analysis on smoothed scores (`patches` in `utils/spatial_patterns.py`, stressed-label shuffle within section x region) was discarded because smoothing builds clustering in. The correlogram replaces it.

**Technical approach.** `correlogram` on the raw (unsmoothed) per-cell ISR score of OLs. For each OL and each distance band (0-25, 25-50, 50-100, 100-200, 200-400 um), the mean ISR of the other OLs in that ring (at least 3). Own score and ring mean are centred within region, then correlated per section (at least 50 cells). The unit is the section. Responders vs controls by exact two-sided permutation per band (84 relabellings).

**Results.** Mean per-section r (controls vs responders, p):
- 0-25 um: 0.047 vs 0.054, p = 0.71
- 25-50 um: 0.051 vs 0.056, p = 0.80
- 50-100 um: 0.055 vs 0.071, p = 0.50
- 100-200 um: 0.073 vs 0.072, p = 0.89
- 200-400 um: 0.035 vs 0.070, p = 0.083

Neighbouring OLs share stress no more in responders than in controls at any distance up to 400 um.

**Reading.** Within a region, responder stress is a uniform upward shift (thalamus, white matter), not focal lesions.

**Caveats.** The single-cell ISR score (about 10 genes) is noisy, so all correlations are small (r ≈ 0.03-0.07). Weak clustering could be missed. Region-centring removes any clustering that lines up with region boundaries, on the old region map (§2, corrected in §43). The 200-400 um band trends higher in responders (p = 0.083).

### 23. No repair response: immature OLs are not more common or more active near stressed OLs
*Notebook:* analysis/06-borders-patches-repair-neurons.ipynb (part C) · *Code:* `utils/spatial_patterns.py`, `utils/neighbourhood.py`, `utils/gene_sets.py` · *Results:* `results/spatial_patterns/tables/repair_coupling_summary.csv`, `results/spatial_patterns/tables/repair_density.csv`, `results/spatial_patterns/tables/repair_coupling_per_section.csv` · *Figure:* `results/spatial_patterns/figures/S3_repair.png`

**Rationale.** Damaged OLs might recruit OPCs or trigger new OL formation. The test is whether immature OLs gather near stressed OLs or switch on proliferation or differentiation genes there. Uses the neighbourhood framework of §17.

**Technical approach.** Immature OL types: "Immature Oligodendrocytes I" (OPC-like), "II" (pre-myelinating), "III" (`IMMATURE_OL`).
- Coupling: for each immature-OL cell, the mean raw ISR of OLs within 50 um (at least 3 OLs; `neighbourhood.local_ol_stress`). This was correlated with the cell's own PROLIF, OPC and NEWLY_FORMED_OL scores, centred within (section, region), per section (at least 50 cells; `neighbourhood.coupling`).
- Density: for each OL, the fraction of all cells within 50 um that are immature OLs (`local_fraction`). This was correlated with the OL's smoothed ISR after (section, region) centring. Whole-section immature-OL fractions were also compared.

Responders vs controls by exact one-sided permutation (84 relabellings). Type III did not have all 9 resp/control sections with enough cells and is absent from the summary.

**Results.**
- Local immature-OL fraction vs OL stress: responders r = -0.007 vs controls 0.020, p = 0.976.
- Section immature-OL fraction: responders 2.6-2.7 % vs controls 2.5-3.4 % (means about 2.7 % vs 2.9 %).
- Coupling in immature OL I (median 886 cells per section): all p ≥ 0.30.
- Coupling in immature OL II (median 169 cells per section): NEWLY_FORMED_OL r 0.124 vs -0.013 (p = 0.071), PROLIF 0.104 vs 0.013 (p = 0.071), OPC 0.074 vs -0.020 (p = 0.14).

**Reading.** No detectable repair attempt at 42 days. There is a weak, non-significant hint in pre-myelinating (type II) cells near stressed OLs.

**Caveats.** Few type II cells per section, so this is low power with 3 responders. The README cites only the PROLIF p = 0.07 and not the equal NEWLY_FORMED_OL p = 0.07. Readouts are short gene-set scores. Regions are from the old map (§2/§43). This is a snapshot, so later repair is not excluded.

### 24. No neuronal axon-injury or stress response near stressed OLs at 42 days
*Notebook:* analysis/06-borders-patches-repair-neurons.ipynb (part D) · *Code:* `utils/neighbourhood.py`, `utils/gene_sets.py` · *Results:* `results/spatial_patterns/tables/neuron_coupling_summary.csv`, `results/spatial_patterns/tables/neuron_coupling_per_section.csv` · *Figure:* `results/spatial_patterns/figures/S4_neurons.png`

**Rationale.** Do neurons whose axons run past stressed OLs show signs of injury (regeneration-associated genes) or their own ISR? This sets the baseline for later time points with demyelination (4-6 months).

**Technical approach.** Same framework as §23. For each neuron of the four neuron classes (Neurons, Excitatory, GABAergic, Medium Spiny), the mean ISR of OLs within 50 um was correlated with the neuron's AXON_INJURY and ISR scores. Both were centred within (section, region), per section (at least 50 cells). Responders vs controls by exact one-sided permutation (84 relabellings), 8 tests, no multiple-testing correction.

**Results.** AXON_INJURY: Excitatory p = 0.63, GABAergic p = 0.70, Medium Spiny p = 0.83. Generic "Neurons" r 0.041 in responders vs -0.014 in controls, **p = 0.036**. ISR: Excitatory p = 0.18, GABAergic p = 0.13, Medium Spiny p = 0.17, Neurons p = 0.54. The notebook and README conclude no neuronal axonal-injury or ISR response.

**Reading.** No consistent neuronal response. Together with §23, the 42-day response is confined to OLs and their glial neighbours.

**Caveats.** The one nominal hit (generic "Neurons" class, AXON_INJURY, p = 0.036 uncorrected among 8 tests, small r) is not mentioned in the notebook text or the README. At this level of testing it is most likely chance, but it is not reported. Spill-over of OL / glial RNA into neuron segments (shown later in §37) could mask or create small effects. 3 responders.

### 25. Each mouse's probability of responding: 3 mice P > 0.99, RB4630 open at P = 0.17
*Notebook:* analysis/07-responder-probability.ipynb · *Code:* `utils/responder_mixture.py` · *Results:* `results/responder_probability/tables/p_responder_primary.csv`, `results/responder_probability/tables/p_responder_by_readout.csv`, `results/responder_probability/tables/prior_sensitivity.csv`, `results/responder_probability/tables/operating_characteristics.csv`, `results/responder_probability/tables/calibration.csv`, `results/responder_probability/tables/next_cohort_rule.csv` · *Figure:* `results/responder_probability/figures/R1_responder_probability.png`

**Rationale.** The hard call (whole-section OL ISR above the highest control, `utils.ol_response.responder_calls`, §10) gives yes / no. With 3 controls per slide, the control range is itself uncertain. A posterior probability per mouse says how sure each call is. It also gives a rule for the next cohort (§20). Committed 2026-10-15 (efd979c).

**Technical approach.** The unit is the mouse (12 sections). Readouts were prespecified or label-free:
- primary: log whole-section OL ISR CPM (10 fixed ATF4 targets);
- the per-OL ISR mean;
- the logit of the out-of-fold control-vs-mtDSB classifier probability (never trained on responder labels).

Driver-gene scores were excluded as circular.

Model, a two-component Bayesian mixture:
- Readouts are centred on the slide's control mean.
- Non-responders follow a multivariate t with the pooled control covariance (df = 6 controls - 2 slides), inflated by 1 + 1/n for slide-mean uncertainty.
- Responders are shifted by delta (uniform grid 0-12 control SDs) with spread k x control (k in 1, 1.5, 2, 3). They must lie above the slide control mean in ISR (one-sided).
- pi ~ Beta(1, 1).

The posterior is exact: all 2^6 labelings are enumerated, pi is integrated out, and delta and k are summed over the grid. Checks:
- Null check: each control in turn is scored as unlabelled and left out of the reference.
- Prior sensitivity: 9 variants.
- Operating characteristics: 300 simulated cohorts of this design (2 slides x 3 controls + 3 mtDSB, 3 of 6 responders, shift = posterior mean 10.6 control SDs, spread = control).
- Posterior predictive rule for a new mouse with 3 slide controls.

**Results.** Primary readout (fold over slide controls, shift in control SDs, P):
- RB4627: 1.73x, 18.6 SD, P = 0.999
- RB4405: 1.61x, 16.2 SD, P = 0.999
- RB4401: 1.39x, 11.2 SD, P = 0.999
- RB4630: 1.07x, 2.3 SD, P = 0.172
- RB4653: 0.97x, P = 0.000
- RB4403: 0.94x, P = 0.000

Further results:
- Labelings: the most probable (3 responders) has posterior 0.825; the one adding RB4630 has 0.171.
- Responding fraction: E[pi] = 0.52, 90 % credible interval 0.23-0.80. Control SD of log ISR is 0.029.
- Null check: no control exceeds P = 0.082 (RB4498).
- Other readouts: the per-OL ISR mean is noisier (control RB4676 P = 0.925; RB4401 0.69). The classifier alone separates poorly (control RB4658 0.90; responders 0.31-0.63). ISR + classifier: responders 0.999, RB4630 0.064, max control 0.222.
- Prior sensitivity: RB4627, RB4405 and RB4401 stay ≥ 0.986 under every prior. RB4630 ranges 0.007-0.336, except 0.667 under "shift up to 6 SD", which the notebook dismisses because that prior cannot fit the observed 11-19 SD responders. The left-out control maximum reaches 0.47 under "shift up to 6 SD" and 0.31 with two-sided responders.
- Operating characteristics: true responders get P > 0.9 in 100 % of cases. Non-responders get P > 0.9 in 1.0 %, P < 0.1 in 86.7 % and 0.1-0.9 in 12.3 %. No simulated mouse with 0.1 < P ≤ 0.9 was a responder.
- Next-cohort rule: P ≥ 0.5 at ≥ 1.11x the slide-control mean, ≥ 0.9 at 1.17x, ≥ 0.95 at 1.19x.

**Reading.** Three mice almost certainly respond and two almost certainly do not. RB4630 is a possible weak responder and the first mouse to check at the bench. Six mice cannot distinguish "half respond" from "a quarter" or "most".

**Caveats.** The model uses the hard-call readout family (ISR) whose responders it confirms. That is circular for the existence of responders, but not for RB4630 or the control null. The simulation assumes the fitted shift (10.6 SD), so it is optimistic for weaker responders. Its intermediate band is conservative rather than calibrated (mean P 0.20 and 0.66 in bins where the observed responder rate was 0). The README states RB4630's range as "0.01-0.34 across priors" without noting the 0.67 variant. With 3 controls per slide, the rule rests on an SD estimated from 6 mice. **Corrected 2026-10-21 (analysis/15, §46):** dropping any single responder keeps the whole-section ISR result at worst p = 0.036.

**Next.** Apply `utils/responder_mixture.py` to the next cohort and report P per mouse.

### 26. T cells are real but rare: 242 co-expressing cells against about 4 in a shuffle null; Trbc1 and Ighm are off-target
*Notebook:* analysis/08-t-cells.ipynb (parts A, B) · *Code:* `utils/tcells.py`, `scripts/export_tcells.py` · *Results:* `results/tcells/panel_summary.csv`, `results/tcells/null_shuffle.csv`, `results/tcells/tables/threshold_grid.csv`, `results/tcells/tables/expression_profile.csv` · *Figure:* `results/tcells/figures/T1_probe_qc.png`, `results/tcells/figures/T2_calls_and_profile.png`

**Rationale.** No cluster or annotation (`cell_class`, `cell_type`, cytetype) contains lymphocytes, yet the IHC reported immune cells at the meninges, ventricle and internal capsule / thalamus. T cells (about 0.01-0.1 % of brain cells) are too rare to cluster. They must be called per cell, and the panel's lymphocyte probes must be checked first. Committed 2026-10-19 (6825e4c).

**Technical approach.** `scripts/export_tcells.py` (office Mac, `sopa` env, polars) read the novae4 h5ad and raw Xenium transcripts. It wrote:
- per-probe QC (`panel_summary.csv`);
- per-cell counts for every cell with ≥ 1 core transcript;
- a shuffle null: each core gene permuted independently across all cells within each section, 20 reps (`tcells.shuffle_null`). This keeps every gene's per-section total but breaks co-expression.

Call (`tcells.call_tcells`): ≥ 3 transcripts from the 9 core genes (Cd3d/e/g, Cd247, Cd2, Cd8a/b1, Itk, Cxcr6) spread over ≥ 2 genes. Trbc1 and Ighm were excluded as noisy; Trac and Nkg7 are not on the panel. Stability was checked on a threshold grid (2-5 transcripts x 1-3 genes). The expression profile of called cells was compared with other cells carrying ≥ 1 core transcript.

**Results.**
- Probe QC: Trbc1 is detected in 29,880 cells (70.4 % excitatory neurons) and Ighm in 67,574 cells (30.9 % excitatory neurons). Off-target. A Trbc1-including score calls about 6,000 "T cells", almost all cortical neurons (notebook text; no saved table for the 6,000).
- Core genes are detected in 331-1,954 cells each (Cd247 1,954, 65 % in "Neurons").
- Calls: 242 cells against a shuffled null of 3.9 (max 8, 20 reps), about 60x.
- Threshold grid: 3 transcripts / ≥ 2 genes gives 242. Other cells of the grid: 2 / ≥ 2 gives 320, 4 / ≥ 2 gives 198, 5 / ≥ 2 gives 166, and ≥ 3 genes gives 160-208.
- Profile of called cells: CD8 (Cd8a/b1) with Cxcr6, Ptprc, Lck, Lat, Il7r and some Gzma; Cd4 is low. Trbc1 is high in the called cells (2.3 per cell, notebook text).
- About 2.5 per 10k cells overall.

**Reading.** CD8 / Cxcr6 tissue-resident-type T cells are present at ordinary brain levels. CD4 T cells are rare or undetectable.

**Caveats.** The call depends on a hand-set threshold, though counts are stable around it. Real T cells with fewer than 3 core transcripts are missed. The "about 6,000 false T cells" and "Trbc1 2.3 per cell" figures are in notebook markdown only. The relaxed-rule null (about 20) is quoted in text but not saved.

### 27. T cells carry their neighbours' labels and sit mostly at meninges, ventricles and vessels
*Notebook:* analysis/08-t-cells.ipynb (parts C, D) · *Code:* `utils/tcells.py`, `scripts/export_tcells.py` · *Results:* `results/tcells/tables/labels.csv`, `results/tcells/zoom_cells.csv.gz`, `results/tcells/zoom_transcripts.csv.gz` · *Figure:* `results/tcells/figures/T3_labels.png`, `results/tcells/figures/T4_overview.png`, `results/tcells/figures/T5_polygon_zooms.png`

**Rationale.** If T cells exist (§26), what were they annotated as, and where are they? Are the calls real cells, or transcript smears from neighbours?

**Technical approach.** Label counts of the 242 called cells in `cell_class_updated`, `cell_type` and `cytetype_annotation_leiden_3`. Whole-section maps with T cells marked, with border defined as old read-based region (`RBD_compartment_simplified`) in {meninges, ventricular system, vasculature}. Close-ups (±45 um) of Xenium cell and nucleus polygons with T-cell transcripts (qv ≥ 20) for the two strongest T cells per section (T5), and for all 242 in `T5b_all_tcells_polygon_zooms.pdf`.

**Results.**
- `cell_class_updated`: Endothelial 133 (55 %), Astrocytes 44 (18 %), DA-Oligodendrocytes 37 (15 %), Oligodendrocytes 6, Myocytes 5, Microglia 4, Excitatory neurons 3, NSC 3.
- `cell_type`: Vascular leptomeningeal I 88, Mature oligodendrocytes V 40, Vascular endothelial I 34, Endothelial I 28.
- cytetype: "Stress-Responsive Oligodendrocyte (GFAP-elevated)" 39, Leptomeningeal 35, Perivascular astrocyte 34.
- No T cell carries a lymphocyte label at any level. Leiden at resolutions 0.5-3 splits them over the same 2-5 large clusters (notebook text).
- Polygons (visual, notebook text): small round cells about 8-12 um, with transcripts inside the polygon around the nucleus. The endothelial-labelled cells sit along vessels and meninges (possibly intravascular). The OL-labelled cells are isolated in tissue (striatum, hypothalamus, fibre tracts) among OLs.

**Reading.** T cells are too few to cluster, so they take the label of what surrounds them. The DA-OL cluster carries 15 % of them.

**Caveats.** Morphology is judged by eye. "Border" uses the old region map. RB4658's region map has no meningeal or ventricular cells under its edge T cells, so its "away from borders" fraction of 1.0 is an annotation artefact. Intravascular and parenchymal positions cannot be separated in 2D.

### 28. Responders do not have more T cells (3.4 vs 2.8 per 10k, p = 0.46)
*Notebook:* analysis/08-t-cells.ipynb (part E) · *Code:* `utils/tcells.py` · *Results:* `results/tcells/tables/per_section.csv`, `results/tcells/tables/responder_vs_control.csv` · *Figure:* `results/tcells/figures/T6_responders_vs_controls.png`

**Rationale.** Is the responder phenotype accompanied by T-cell infiltration, as in inflammatory demyelination? Also relevant to the IHC sites and to the EAE comparison (§31).

**Technical approach.** Per section: T cells per 10k cells, T cells away from borders per 10k, the fraction away from borders, and the fraction labelled oligodendrocyte (`cell_class_updated` contains "Oligodendrocytes"). Responders (3) vs controls (6) by exact two-sided permutation (84 relabellings).

**Results.** 6-33 T cells per section.

| Metric | Controls | Responders | Non-responders | p |
|---|---|---|---|---|
| Per 10k cells | 2.82 | 3.35 | 1.11 | 0.46 |
| In tissue, per 10k | 1.16 | 1.82 | 0.56 | 0.29 |
| Fraction away from borders | 0.44 | 0.53 | 0.44 | 0.70 |
| Fraction OL-labelled | 0.11 | 0.30 | 0.13 | 0.060 |

- Fraction OL-labelled in responders: RB4627 0.41, RB4405 0.39, RB4401 0.11. Controls range 0-0.19.
- RB4405 has 3.07 tissue T cells per 10k. Notebook text says "≤ 2.2 in any control", but RB4658 is 2.25 (its value is an artefact, §27).
- Non-responders have the fewest T cells (0.8-1.7 per 10k).

**Reading.** mtDSB does not raise T-cell numbers. In two responders more T cells sit in tissue among OLs. With n = 3 this is a lead to check against the IHC or flow data, not evidence of T-cell-driven damage.

**Caveats.** About 10-30 T cells per section, 3 responders. Border calls rely on the old region map. Detection is limited by the call threshold (§26).

### 29. Non-responders hide no responding OLs: flat upper tails overall and in every region
*Notebook:* analysis/09-mosaic-response.ipynb (parts A, B) · *Code:* `utils/mosaic.py` · *Results:* `results/mosaic/tables/tail_enrichment.csv`, `results/mosaic/tables/region_excess_q99_isr.csv`, `results/mosaic/tables/region_excess_q99_drivers_up.csv` · *Figure:* `results/mosaic/figures/M1_percentiles_and_tails.png`, `results/mosaic/figures/M2_region_tails.png`

**Rationale.** The leading explanation for the 3 non-responders (§10-11, §25) is failed or partial induction. If induction worked in only some OLs, a non-responder would carry a small set of strongly stressed OLs on top of a normal majority, invisible in section averages. Committed 2026-10-21 (e7ae7ee).

**Technical approach.** Input: `results/cell_table/ol_gene_scores.csv.gz` (159,301 OLs, 12 sections; old region map, "Unknown" removed). Scores:
- the prespecified ISR score (primary);
- the up driver genes (§15; data-derived, so circular for responders but informative for non-responders).

`loo_percentiles`: each OL's percentile among control OLs of the same slide and region, own section left out, ties broken at random, strata needing ≥ 50 reference cells. Uniform under the null. `tail_enrichment`: observed / expected share above the control 90th, 99th and 99.9th percentile, per section and per section x region (6 main regions).

**Results.** ISR, share in the control top 1 % (observed / expected):
- controls 0.88-1.21
- non-responders: RB4403 0.61, RB4653 0.95, RB4630 1.22
- responders: RB4401 5.0, RB4405 7.8, RB4627 8.3

Top 0.1 %: non-responders 0.47-1.70, controls 0.86-2.24, responders 7.1-14.0.

Driver genes, top 1 %:
- controls 0.35-2.72
- non-responders 0.52-1.67 (RB4630 0.75)
- responders 15.1-34.1

Regions:
- The highest non-responder region (outside hippocampus) is RB4630 fibre tracts at 1.6x for ISR and RB4653 cortex at 3.1x for drivers. Controls reach up to 2.3x and 3.3x outside hippocampus.
- Responders are elevated in every region: ISR 2.5-11.5x, strongest in thalamus and fibre tracts.
- The hippocampus has as few as 88 OLs per section, so its ratios swing (e.g. RB4630 2.6x, control RB4350 4.5x).

**Reading.** Non-responders are flat like controls at the single-cell level and in every region. If induction failed, it failed almost completely. RB4630 is at the top of the control range for ISR (1.2x) but not for drivers (0.7x).

**Caveats.** The ISR single-cell score is noisy. The region map is the old one: §44 (analysis/14 part 5) re-ran the regional tail on the disease-free map for responders vs controls (all regions hold), but the non-responder pocket check itself was not re-run. **Corrected 2026-10-21 (analysis/15, §46):** "OLs in control top 1 % (ISR)" responders vs controls p = 0.012, worst leave-one-responder-out p = 0.036.

### 30. A non-responder could hide at most about 3 % responding OLs; in responders at least 20-36 % respond
*Notebook:* analysis/09-mosaic-response.ipynb (parts C, D) · *Code:* `utils/mosaic.py` · *Results:* `results/mosaic/tables/spike_in_detection.csv`, `results/mosaic/tables/responding_fraction_lower_bound.csv` · *Figure:* `results/mosaic/figures/M3_detection_power.png`, `results/mosaic/figures/M4_responding_fraction.png`

**Rationale.** A flat tail (§29) means something only if a small responding share would have been visible. The same percentiles also bound what share of OLs respond in responders.

**Technical approach.**
- Spike-in (`spike_in_detection`): a crude responding component was taken from responder percentiles (histogram minus its floor, `responding_component`). Fractions f = 0-20 % of a control section's OLs were replaced with draws from it, 20 reps x 6 controls. Detection means the section's tail enrichment exceeds every other control's. The chance rate at f = 0 is 1/6.
- Lower bound (`responding_fraction_lb`): with F = (1 - pi) U + pi G, pi ≥ max_u (S(u) - (1 - u)) / u over u in 0.5-0.999. This is nonparametric.

**Results.** Detection by the ISR top-1 % rule:
- 1 %: 0.33
- 2 %: 0.65
- 3 %: 0.98
- 5 % and above: 1.0

ISR top-10 % reaches 1.0 at 5 %. Drivers top-1 % needs 10 % (0.89) to 20 % (1.0), because control RB4350 has a high driver tail.

Lower bound on the responding share, ISR:
- responders: RB4627 0.36, RB4405 0.33, RB4401 0.20
- controls: 0-0.07
- non-responders: 0-0.05

Lower bound, drivers:
- responders: 0.79, 0.68, 0.39
- controls: 0-0.16
- non-responders: 0.04-0.17

**Reading.** At least a fifth to a third of OLs respond in responders (plausibly 40-80 % by drivers), in every region. A non-responder could hide less than a tenth of that. The response is all-or-none per mouse, which fits a per-animal event (failed induction, or a whole-brain threshold) rather than a stochastic per-cell one.

**Caveats.** The 3 % sensitivity depends on the crude responding component and is approximate (notebook). Lower bounds from a noisy 10-gene score underestimate the true share. The driver score is circular for responders.

### 31. The OL stress (ISR) and p21 programme is mtDSB-specific; the immune-visible core is shared with EAE
*Notebook:* analysis/10-eae-comparison.ipynb (parts A, B) · *Code:* `scripts/export_eae_pseudobulk.py`, `utils/gene_sets.py`, `utils/ol_response.py` · *Results:* `results/eae_comparison/tables/ol_log2fc_mtdsb_vs_eae.csv`, `results/eae_comparison/tables/programme_log2fc.csv`, `results/eae_comparison/eae_pseudobulk_meta.csv` · *Figure:* `results/eae_comparison/figures/E1_gene_scatter.png`, `results/eae_comparison/figures/E2_programmes.png`

**Rationale.** Responder OLs switch on stress genes, become immune-visible and lose myelin genes (§15). Disease-associated OLs (DAOs) in neuroinflammation are also immune-visible. Is the mtDSB response specific, or the generic injured-OL response? The RR-map EAE atlas used the same 5k panel, so genes compare one to one. Committed 2026-10-21 (e7ae7ee).

**Technical approach.** EAE data: RR-map spinal-cord Xenium atlas (`RRMAP2 ... Anno2to4Updated.h5ad`, about 1.4M cells, read in chunks by `scripts/export_eae_pseudobulk.py` on the laptop). It covers 8 control (CFA) and 59 EAE mice, chronic (MOG) and relapsing-remitting (PLP), onset to remission. Pseudobulk raw counts were summed per mouse x L1 x L2 x curated niche.

Contrasts, as log2(1 + CPM) of OL pseudobulks per mouse, mean differences:
- mtDSB: responders minus same-slide controls, averaged over the two slides; Oligodendrocytes + DA-Oligodendrocytes pooled; region x time pseudobulks from §8.
- EAE OLs (MOL / DAO / NFOL) in lesion niches (active / progressive / mix) minus control OLs.
- EAE OLs in normal-appearing ("Physiological") niche minus control OLs.
- DAO minus MOL within EAE mice.

Genes were kept if > 20 CPM in mtDSB OLs (2,967 genes). Programmes are mean log2 FC over member genes passing that filter: driver programmes (§15) and prespecified sets.

**Results.** Genome-wide r of mtDSB responder log2 FC with EAE lesion 0.29, normal-appearing 0.10, DAO vs MOL 0.18.

Gene log2 FC (mtDSB / EAE lesion):
- Trib3 5.41 / 0.75
- Cdkn1a 3.97 / 0.03
- Atf5 2.07 / 0.62
- Sesn2 1.83 / 0.12
- Ddit3 1.22 / 0.09
- Ccnd1 2.64 / 0.86
- Cd40 5.09 / 0.87
- Nmu 4.74 / -0.28
- shared: Serpina3n 2.87 / 4.73, Cd74 1.76 / 4.04, H2-D1 1.39 / 3.64, B2m 0.72 / 2.60, Lgals1 1.51 / 1.88, Hmox1 1.12 / 2.01
- EAE-specific: Irf1 0.40 / 3.56, Psmb8 0.57 / 3.79

| Programme mean log2 FC | mtDSB | EAE lesion | EAE normal-appearing | DAO vs MOL |
|---|---|---|---|---|
| ATF4 / ISR drivers (8 genes) | 1.85 | 0.06 | -0.15 | 0.06 |
| Prespecified ISR (6 of 10 genes) | 1.83 | 0.21 | -0.02 | 0.02 |
| Cell cycle / survival | 2.52 | 0.67 | | |
| Immune / reactive | 2.25 | 2.41 | | |
| MHC I | 0.83 | 3.27 | | |
| IFN-gamma | 0.32 | 3.11 | | |
| Type I IFN | 0.19 | 1.27 | | |
| Myelin, prespecified | -0.49 | -0.15 | | |

**Reading.** mtDSB OLs carry an ISR + p21 / cell-cycle arm that EAE OLs do not. EAE OLs carry an interferon-gamma / MHC arm that mtDSB lacks. Both share an immune-visible core. A Trib3 / Gdf15 / p21 panel could mark mitochondrial-damage OLs in other tissues (pursued in §50-55).

**Caveats.** Brain vs spinal cord; one time point (42 d) vs a disease course. Lesion OL pools contain myeloid spill-over (Aif1, Ctss, Itgam) and surviving or newly formed OLs. Gdf15 and Eif4ebp1 fall below the 20 CPM filter: they are absent from the saved table and from the programme means and correlations. The notebook text still quotes Gdf15 +4.5 vs 0.0. Recomputed here from the saved pseudobulks, it is 4.47 vs 0.04 (pooled CPM 9.5), so the numbers are right but untraceable in the outputs. The prespecified ISR mean uses only 6 of 10 genes. Fold-change wording differs: README "up 3-40x" vs notebook summary "4-40x". Ddit3 is 2.3x (1.22 log2), so neither lower bound is exact. Mean differences only, with no per-gene test.

### 32. No EAE stage switches on the OL ISR, though the top EAE mouse overlaps the mildest responder on z
*Notebook:* analysis/10-eae-comparison.ipynb (part C) · *Code:* `scripts/export_eae_pseudobulk.py` · *Results:* `results/eae_comparison/tables/eae_per_mouse_ol_scores.csv`, `results/eae_comparison/tables/mtdsb_per_mouse_ol_isr.csv` · *Figure:* `results/eae_comparison/figures/E3_isr_by_stage.png`

**Rationale.** Pooled fold changes (§31) could hide a stage of EAE (e.g. peak) where OLs do switch on the ISR.

**Technical approach.** Per-cell scores were computed in the export script as log1p(1e4 x set counts / total) for the ISR set. The per-mouse mean was taken over OLs (`Anno_L1_curated == "Oligodendrocyte"`) in lesion or normal-appearing niche, keeping mouse x niche groups with ≥ 50 OLs. Each was z-scored against the mean and SD of control mouse x niche rows (12 rows: 8 normal-appearing, 4 "lesion"-niche rows in CFA controls). The values were plotted by stage. mtDSB: per-mouse mean OL ISR (`ol_gene_scores.csv.gz`) z-scored against the 3 controls of the same slide.

**Results.**
- EAE: max z 4.10 (C_P1_5, normal-appearing, PEAK1, chronic). Next highest: 3.15 (C_P1_4 lesion), 2.79, 2.53, 2.47.
- Mean z: EAE lesion -0.30 (n = 56), EAE normal-appearing +0.85 (n = 54).
- mtDSB responders z: RB4401 3.99, RB4627 6.37, RB4405 8.88. Non-responders 0.03-0.97.

The notebook reads this as EAE within about -4 to +4 SD (mostly -2 to +2) and responders at +4 to +9.

**Reading.** Combined with the gene-level fold changes (§31), inflammatory demyelination does not drive the OL ISR at any stage.

**Caveats.** The z scales are not comparable in precision (8 EAE controls vs 3 per mtDSB slide). The single highest EAE mouse (4.10) exceeds the mildest responder RB4401 (3.99), so on this per-mouse z the separation is not complete. The notebook text glosses over this. Normal-appearing EAE OLs sit slightly above controls on average (+0.85). CFA controls carry "lesion" niche labels, an annotation oddity, and these were pooled into the control reference.

### 33. Microglia and astrocytes react similarly in mtDSB and EAE, milder in mtDSB
*Notebook:* analysis/10-eae-comparison.ipynb (part D) · *Code:* `scripts/export_eae_pseudobulk.py`, `utils/gene_sets.py` · *Results:* `results/eae_comparison/tables/glia_programmes.csv`

**Rationale.** If the OL response is mtDSB-specific (§31), is the glial reaction around it also specific, or a common injury response? Follows the microglia / astrocyte coupling of §17.

**Technical approach.** Same pseudobulk contrasts as §31 for microglia (mtDSB "Microglia"; EAE L1 "Myeloid", L2 "Microglia") and astrocytes (mtDSB "Astrocytes"; EAE L1 "Astrocyte"). The comparison was mtDSB responders vs slide controls against EAE lesion niche vs controls. The genome-wide r used genes > 20 CPM in mtDSB. Prespecified sets: DAM, HOMEOSTATIC, MHC1, IFN, ASTRO_REACTIVE.

**Results.** Microglia (mtDSB vs EAE lesion):
- r = 0.61 across 3,228 genes
- DAM +2.65 vs +2.54
- homeostatic -0.20 vs -0.95
- MHC I 0.81 vs 2.05
- IFN 0.62 vs 1.41

Astrocytes (mtDSB vs EAE lesion):
- r = 0.43 across 3,304 genes
- reactive +0.83 vs +3.03
- MHC I 0.48 vs 2.43
- IFN 1.03 vs 1.94

**Reading.** The neighbouring glia react much more alike than the OLs do. Microglia show a DAM response of equal size, and EAE adds more interferon / MHC. mtDSB glia look like a milder version of a common injury response, while the OL carries the distinct signature.

**Caveats.** There is no figure for this part. Brain vs spinal cord. Pseudobulks of segmented cells carry spill-over from neighbours (addressed for mtDSB in §37, not here). The table has cross-type rows (ASTRO_REACTIVE in microglia) that are not meaningful. There is no statistical test.

### 34. Responders have 1.5-3x more microglia per mm2 in thalamus and fibre tracts, and their sections seemed to lack pallidum
*Notebook:* analysis/11-composition-microglia.ipynb (parts A, B) · *Code:* `utils/composition.py` · *Results:* `results/composition/tables/region_coverage_pct.csv`, `results/composition/tables/density_tests.csv`, `results/composition/tables/density_long.csv` · *Figure:* `results/composition/figures/C1_region_coverage.png`, `results/composition/figures/C2_density_heatmap.png`, `results/composition/figures/C3_density_dots.png`

**Rationale.** Earlier runs measured cell *states* (OL ISR §10-11, microglial DAM §17, §19). This run counts *cells*: do responders gain microglia or lose OLs, and in which regions? Part A first checks that every section samples the same regions, because density comparisons assume comparable anatomy.

**Technical approach.** 980,474 cells in 12 sections from `results/cell_table/cells.csv.gz`, with the original transcript-based region labels (`RBD_compartment_simplified`, §2). OLs and DA-OLs are counted together, since the DA-OL label partly captures responding OLs (§11). A: share of each section's cells per region. B: density = cells / area, where area = number of 50 um grid squares holding at least one cell of that region (`occupied_area`, `region_density`); 7 cell types x 6 regions, density and share. Unit = mouse. Responders (3) and non-responders (3) are each compared with controls (6) by an exact permutation test over all 84 relabellings (smallest two-sided p = 1/84 = 0.012); effect = log2 (group mean / control mean).

**Results.** A: pallidum is 0.04-0.68 % of cells in the three responders vs 3.1-7.3 % in every control and non-responder (p = 0.012); fibre tracts are larger in responders (9.8 vs 7.0 %, p = 0.012); other regions do not differ. B: microglia density in thalamus 92-153 /mm2 in responders vs 60-82 in controls (mean 119.7 vs 70.8, log2 +0.76, p = 0.012), in fibre tracts 71-219 vs 51-68 (134.4 vs 60.3, log2 +1.16, p = 0.012); hypothalamus 93.1 vs 64.7 (p = 0.060); cortex, striatum and hippocampus unchanged. Non-responders lie in the control range (thalamus p = 0.64, fibre tracts p = 0.15). OL (+ DA-OL) density is lower in responder cortex (246 vs 325 /mm2, p = 0.036) and fibre tracts (1,396 vs 1,674, p = 0.036) with overlapping ranges; immature OLs, astrocytes and endothelial cells unchanged. "Neurons" in responder fibre tracts are 3.7x controls (88 vs 24 /mm2, p = 0.024). The fibre-tract microglia range is wide: RB4405 219 /mm2, RB4401 71.

**Reading.** A white-matter microgliosis in responders only. The pallidum gap plus neurons inside "fibre tracts" suggested that grey-matter tissue was being labelled fibre tract, either by section plane or by disease distorting the region map.

**Caveats.** Region labels come from pooled neighbourhood transcripts and might shift with disease; this run could not tell plane from distortion (resolved in §43). Neuron densities are noisy. 3 vs 6 mice. The neuronal "excess" is likely a labelling artefact.

**Corrected 2026-10-21:** §43-44 showed the pallidum gap was the disease distorting the region map, not section plane; on the disease-free map thalamus microgliosis holds (p = 0.012), hypothalamus becomes significant (p = 0.036), and the fibre-tract excess weakens (p = 0.036, RB4401 inside the control range) because part of it was pallidum.
**Corrected 2026-10-21 (analysis/15, §47):** this 42-test density screen does not survive FDR on its own: only 6 tests have p < 0.05 and the two at the floor (microglia in thalamus and fibre tracts) reach q = 0.25. The microgliosis claim rests on the label-free replication (§35), the disease-free map (§44) and leave-one-responder-out (§46). The OL density drops are nominal.

**Next.** Iba1 counts in thalamus, pallidum and internal capsule (README finding 14).

### 35. The microglia excess sits in OL-rich tissue whatever the region label, with a modest rise in proliferating microglia
*Notebook:* analysis/11-composition-microglia.ipynb (parts C, D) · *Code:* `utils/composition.py` · *Results:* `results/composition/tables/microglia_density_by_ol_fraction.csv`, `results/composition/tables/microglia_proliferation.csv` · *Figure:* `results/composition/figures/C4_microglia_by_ol_fraction.png`, `results/composition/figures/C5_microglia_proliferation.png`

**Rationale.** §34 found the pallidum gap, which put every region-based result in doubt. Part C repeats the microglia result without region labels; part D asks whether the extra microglia are dividing.

**Technical approach.** C: each section is cut into 100 um squares (squares with >= 20 cells kept); each square is binned by its OL fraction (0-10 %, 10-30 %, 30-50 %, 50-100 %; a white-matter proxy) and microglia per mm2 is computed per section per bin (`ol_rich_bins`). Number of squares per section per bin is printed (e.g. 93-231 squares at 50-100 %). D: a microglia is "proliferating" if it carries >= 1 transcript of the PROLIF set (Mki67, Top2a, Cdk1, Pcna, ...). Because detection depends on sequencing depth and region, a binomial GLM (region dummies + log total counts) is fitted on control microglia only and gives each cell an expected probability; per mouse observed / expected is compared. Also proliferating microglia per mm2 of section. Unit = mouse; exact permutation tests, responders and non-responders vs controls.

**Results.** C: responders have more microglia in 10-30 % OL squares (p = 0.024) and 30-50 % OL squares (p = 0.012), same direction at 50-100 % (p = 0.083), none in OL-poor < 10 % squares (p = 0.75). Graded per mouse: RB4405 147-229 /mm2 in the 10-50 % bins vs controls 79-109; RB4627 103-123; RB4401 108-111, and 69 in the 50-100 % bin (control range 62-88). D: 15-17 % of responder microglia carry a proliferation transcript vs 11-15 % in controls (p = 0.012); observed / expected 1.17-1.23 vs 0.87-1.14 (p = 0.024 recomputed from the saved table; the notebook text and README say 0.012); proliferating microglia per mm2 of section 14.6-23.0 vs 9.8-13.2 (p = 0.012). Non-responders are around 1.0 observed / expected (p = 0.83).

**Reading.** The microgliosis is a white-matter phenomenon independent of how regions were drawn; local proliferation probably contributes.

**Caveats.** One transcript is a weak proliferation call (the panel lacks several cell-cycle genes); Ki67/Iba1 co-staining would settle it. The notebook says RB4401 is "barely above controls in the purest white matter", but in the 50-100 % bin it is 69 /mm2, inside the control range (62-88). The non-responder RB4403 has high label-free densities (128 and 137 /mm2 at 10-30 % and 30-50 %, above every control); analysis/15 attributes this to its over-segmentation (doubled cell count), so "none in non-responders" holds for region densities (§34) but needs that explanation here.

**Corrected 2026-10-21 (analysis/15, §46-47):** label-free OL-rich microglia density holds with any one responder left out (p = 0.036) and keeps q < 0.10; proliferating microglia are "consistent but weaker" (worst drop-one p = 0.071).

### 36. Microglia are spread through the affected white matter, not piled on the most stressed OLs
*Notebook:* analysis/11-composition-microglia.ipynb (parts E, F) · *Code:* `utils/composition.py`, `utils/mosaic.py` · *Results:* `results/composition/tables/microglia_near_ol_by_isr_band.csv` · *Figure:* `results/composition/figures/C6_microglia_near_stressed_ols.png`, `results/composition/figures/C7_maps.png`

**Rationale.** If microglia were recruited by individual stressed OLs (§17, §19 found weak local coupling), the extra microglia of §34-35 should cluster around the most stressed OLs.

**Technical approach.** Each OL (`results/cell_table/ol_gene_scores.csv.gz`, region "Unknown" dropped) gets an ISR percentile relative to same-slide, same-region controls with leave-one-out (`mosaic.loo_percentiles`, §29). OLs are binned (< 50, 50-90, 90-95, 95-99, > 99th percentile). Microglia within 50 um of each OL are counted (`neighbour_counts`). Within each section x region, the mean count per band is divided by the mean for OLs below the 50th percentile (cells >= 10 per band), then averaged over regions weighted by OL number. Controls get the same calculation for their own "top" OLs, which controls for expression co-varying with local cell crowding. F: whole-section maps of microglia over OLs for six sections.

**Results.** An OL above the control 99th percentile has 0.90-1.15x the microglia of a below-median OL in controls and 1.05-1.28x in responders (RB4627 1.05, RB4405 1.18, RB4401 1.28); responders vs controls on the > 99 band p = 0.14 (recomputed from the saved table; the figure title prints the same test). Non-responders 0.96-1.24. Maps: in RB4405 and RB4627 the extra microglia concentrate in thalamus and surrounding white matter (a dense thalamic patch in RB4405); RB4401 shows a smaller thalamic patch; whole-section density hides this (RB4401 98 vs a control 94 /mm2).

**Reading.** The microgliosis is regional and diffuse, matching the diffuse OL response (§22) and weak DAM-OL coupling (§17).

**Caveats.** 50 um radius only; a single-band ratio per mouse with 3 vs 6 mice has low power, so a small focal component cannot be excluded. Maps are qualitative.

### 37. After removing segmentation spill-over, only microglia, astrocytes and some endothelial cells keep their own response
*Notebook:* analysis/12-celltype-response-spillover.ipynb (parts A-E) · *Code:* `utils/spillover.py`, `utils/gene_sets.py` · *Results:* `results/celltype_de/tables/de_by_celltype.csv`, `results/celltype_de/tables/contamination_alpha.csv`, `results/celltype_de/tables/source_markers.csv`, `results/celltype_de/tables/own_changes_with_region_lfc.csv`, `results/celltype_de/tables/programmes_in_own_changes.csv` · *Figure:* `results/celltype_de/figures/D1_contamination.png`, `results/celltype_de/figures/D2_own_vs_spill.png`, `results/celltype_de/figures/D3_surviving_genes.png`

**Rationale.** §8 compared all mtDSB mice with controls; with the responder split (§10-11) each cell type can be tested for responders and non-responders separately. But when OLs lose myelin RNA and microglia switch on Cst7, a diluted copy of both changes appears in every neighbouring cell type through segmentation errors; this run separates own from spilled changes.

**Technical approach.** Input: whole-section pseudobulks per cell type and mouse (`results/region_time/pseudobulk_*`, from §8). A: pyDESeq2 per cell type, `~ age + group`, responders vs controls and non-responders vs controls; hit = padj < 0.05, |log2FC| > 0.5, baseMean > 20. B: linear mixing model on pseudobulk CPM (`spillover.py`): source-specific markers (CPM > 50 and >= 5x every other source type in controls; 43 OL, 133 microglia, 48 astrocyte, 126 endothelial, 5 neuron markers); contamination alpha[target, source] = median ratio of marker CPM in target vs source. Group means average the two slides (ages) equally. C: expected spill-over of each gene in each target = sum over other sources of alpha x the source's responder-control CPM change; a hit is "spill-over" if this predicts >= 50 % of the observed change (same sign), else "own". D: each own change must reappear (same sign, |log2FC| > 0.3, region pseudobulk with +5 CPM pseudocount) in >= 2 of 7 regions. E: surviving genes mapped to prespecified programmes (`GENE_SETS`).

**Results.** A: non-responders 0-4 hits per cell type; responders: microglia 152, neural stem cells 133, endothelial 59, astrocytes 45, OLs 43, neurons 38, others 1-9. B: contamination roughly 2-5 % of OL, microglial and endothelial signal and 6-14 % of astrocytic signal in other types; NSCs carry 32 % astrocytic signal; immature OL II/III 0.24-0.73 OL signal (partly mature OLs, so an upper bound). The neuron column (5 markers) is least reliable. C: own / spill-over: microglia 135 / 17, NSC 133 / 0, endothelial 46 / 13, astrocytes 32 / 13, OLs 33 / 10, neurons 12 / 26, immature OL I 1 / 8. Spilled hits come mainly from OLs (myelin genes) and microglia (Cst7 / Ccl3), plus astrocytes (Gfap / C4b). D: region-consistent own changes: microglia 135 / 135, astrocytes 32 / 32, OLs 31 / 33, endothelial 25 / 46, neurons 10 / 12; NSCs untestable (one region). E: microglia DAM / lipid handling (14 DAM genes incl. Cd74 +4.1, Axl +2.2; Csf1, Cxcl16, Plau / Plaur), MHC-I, P2ry12 -0.6, Mki67 +1.9; astrocytes reactive (Serpina3n, C4b, Gfap, A2m) and type-I interferon (Ifit1, Ifit3, Isg15, Oasl2, Usp18, Cxcl10 +2.4, Nlrc5); endothelial Ifit3, Gfap, Timp1; OLs ISR (Trib3 +4.4, Atf5 +1.2, Sesn2 +0.9).

**Reading.** Induction failure in non-responders is tissue-wide. In responders, most neuronal and immature-OL "responses" are neighbours' RNA; the genuine non-OL response is microglial DAM plus astrocyte reactivity and interferon. Surviving neuronal "down" genes (Foxd2, Gsx1, Fzd5, Fgl1) are nucleus markers, i.e. section plane.

**Caveats.** The mixing model is first order: it ignores region-dependent contamination and that more microglia means more microglial spill-over, so it is a filter, not a correction. Some own calls look like residual spill-over (e.g. astrocyte Mbp -0.57 is called own). The notebook text says all OL own changes pass the region check; the printout shows 31 of 33.

### 38. Astrocytes switch on a small, patchy type-I interferon programme that is not organised around stressed OLs
*Notebook:* analysis/12-celltype-response-spillover.ipynb (part F) · *Code:* `utils/composition.py`, `utils/neighbourhood.py` · *Results:* `results/celltype_de/tables/astro_ifn_positive_pct_by_region.csv`, `results/celltype_de/tables/astro_ifn_coupling.csv` · *Figure:* `results/celltype_de/figures/D4_astro_ifn_by_region.png`, `results/celltype_de/figures/D5_astro_ifn_coupling.png`

**Rationale.** §37 left an astrocyte type-I interferon programme (Isg15, Ifit1/3, Usp18, Cxcl10) in pseudobulk. A cGAS-STING-type response to released mtDNA would predict IFN in or near the damaged OLs; this checks it per cell.

**Technical approach.** Per astrocyte, IFN gene-set transcripts > 0 = "IFN-positive" (`cell_gene_scores.csv.gz`). Share per section and region (6 regions), responders and non-responders vs controls by exact permutation. Within each section: correlation between IFN-positivity and (i) microglia within 50 um and (ii) mean OL ISR within 50 um (`local_ol_stress`), both centred within region and with local total cell count regressed out of the predictor. Unit = mouse for the comparison.

**Results.** IFN-positive astrocytes: fibre tracts 12.8-19.4 % in responders vs 7.1-10.9 % in controls (p = 0.012), thalamus 8.6-9.2 % vs 6.6-7.7 % (p = 0.012); hypothalamus p = 0.083, cortex / striatum / hippocampus p >= 0.40. The non-responder RB4653 has the highest fibre-tract value (23.1 %). Within sections, r is within +/-0.03 everywhere: microglia r -0.030 to +0.032, OL ISR r -0.014 to +0.024. Responders vs controls on the OL-ISR coupling p = 0.26; on the microglia coupling p = 0.024 (recomputed; driven by RB4405 r = 0.032, the notebook text does not mention this p).

**Reading.** The interferon programme is real in pseudobulk and not spill-over, but small, patchy and not spatially tied to the OL lesion, unlike the IFN-gamma response in EAE (§33).

**Caveats.** "Any IFN transcript" is a coarse call; RB4653 shows region-level IFN-positivity is not responder-specific. Correlations are tiny in all mice.

### 39. The neural stem cell "proliferative response" is section plane: responders hold more lateral-ventricle wall, not more dividing NSCs
*Notebook:* analysis/12-celltype-response-spillover.ipynb (parts D, G) · *Code:* `scripts/export_cell_counts.py`, `utils/composition.py` · *Results:* `results/celltype_de/tables/nsc_niche_and_proliferation.csv`, `results/celltype_de/tables/own_changes_with_region_lfc.csv` · *Figure:* `results/celltype_de/figures/D6_nsc_niche.png`

**Rationale.** NSCs had 133 own hits in §37, all untestable for region consistency, and the manuscript (paper sec. 7, §6; questioned in §7) claims an NSC proliferative response. This run asks whether those hits are biology or which ventricular wall the section contains.

**Technical approach.** D: the NSC hits were read gene by gene. G: per-cell raw counts of selected genes in NSCs (`results/cell_table/nsc_ol_lineage_counts.csv.gz`, made by `scripts/export_cell_counts.py` from the h5ad). Per NSC: proliferating = any of Mki67, Top2a, Ccnb1, Cenpf, Pcna, Mcm2, Ccna2 > 0; lateral-ventricle wall = Foxg1 > 0; third ventricle / tanycyte = Tbx3 > 0. Per mouse: share Tbx3+, share proliferating, and share proliferating within Foxg1+ NSCs (same niche). Responders vs controls by exact permutation; P21 and P60 shown separately because they sit on different slides.

**Results.** D: NSC "down" genes are third-ventricle / tanycyte genes (Gpc3 -6.1, Tbx3 -4.2, Col23a1 -5.4, Frzb -6.2, Rara, Hes1); "up" genes are telencephalic / lateral-ventricle neurogenic genes (Foxg1 +1.2, Elavl2, Erbb4, Sema6c) with cell-cycle genes (Ccnb1 +1.6, Top2a +1.4, Ccna2, Cenpf, Knl1). G: at P60, 11.4-14.0 % of control NSCs are Tbx3+ vs 0.7 % (RB4401) and 0 % (RB4405); P21 sections have almost none except the non-responder RB4653 (25.7 %). Share proliferating does not differ (p = 0.45); within Foxg1+ NSCs responders sit in their age's control range (P60 57-58 % vs 52-61 %; P21 RB4627 73 % vs 44-66 %; p = 0.38). NSCs per section 214-1,041.

**Reading.** No evidence for an NSC proliferative response; the pseudobulk signal comes from sections holding different ventricular walls. The manuscript's NSC point should be dropped or rephrased (README "What we cannot claim").

**Caveats.** Single-transcript proliferation calls; 77-262 Foxg1+ NSCs per section; P21 RB4627 is slightly above its controls (73 vs 66 %), one mouse.

### 40. Inside responding OLs the programme is ATF4 / ISR plus p53, with no cGAS-STING, interferon, mitophagy or classic SASP
*Notebook:* analysis/13-ol-mechanism-myelin-rna.ipynb (part A) · *Code:* `utils/gene_sets.py`, `utils/spillover.py`, `utils/composition.py` · *Results:* `results/ol_mechanism/tables/pathway_scores.csv`, `results/ol_mechanism/tables/pathway_genes.csv` · *Figure:* `results/ol_mechanism/figures/O1_pathway_map.png`

**Rationale.** mtDNA double-strand breaks could act through several pathways: the ISR (found in §10, §15), p53 / DNA damage, mitochondrial UPR, mitophagy, or innate sensing of released mtDNA (cGAS-STING / type-I IFN). The earlier review (§7) found no responder pattern for cGAS-STING in legacy work. This run tests prespecified pathway sets in OLs, with a spill-over check for low-count IFN genes.

**Technical approach.** OL + DA-OL whole-section pseudobulk per mouse (from §8), CPM, log2(CPM + 1), centred on the mean of same-slide controls (slide = age). 11 prespecified sets (`MECHANISM_SETS`, 4-14 genes each: ATF4 / ISR, p53 targets, mito UPR, mitophagy, mtDNA / biogenesis, cGAS-STING / IFN-I, apoptosis, senescence / SASP, ER UPR, NRF2, cholesterol). Pathway score = mean gene log2FC per mouse; responders vs controls by exact permutation. Per gene, the spill-over fraction from §37's mixing model (share of the OL change predicted by other cell types' changes).

**Results.** Pathway scores (responder mean log2FC; control range; p): ATF4 / ISR +1.86 (-0.14 to 0.19; 0.012), p53 +1.03 (0.012), senescence / SASP +0.59, NRF2 +0.45, mito UPR +0.39, ER UPR +0.31, cGAS-STING / IFN-I +0.29 (-0.18 to 0.14; p = 0.024), mtDNA / biogenesis +0.27, apoptosis +0.20, mitophagy +0.07 (p = 0.036), cholesterol -0.42 (p = 0.036). Non-responders -0.12 to +0.17. Genes: Trib3 +5.4, Gdf15 +4.3, Eif4ebp1 +3.0, Atf5 +2.0, Aldh1l2 +1.9, Sesn2 +1.7, Psat1 -0.5; Cdkn1a +4.0, Bbc3 (PUMA) +2.1, Gadd45a +1.7, Pmaip1 and Fas flat; Bcl2l1 +1.0, caspases flat; Hmox1 +1.1, Hspa9 +0.8; Pink1, Prkn, Bnip3, Fundc1 flat, Sqstm1 +0.6; Serpine1 +2.4 and Igfbp3 +2.3 but Cdkn2a, Il1a, Cxcl1, Ccl2, Mmp3 flat; Cgas +0.48 from 1.1 CPM, Tmem173, Tbk1, Irf3, Ifnb1 flat; Ifit1 +0.94, Usp18 +0.54, Oasl2 +0.47, Isg15 +0.34 each with spill-over fraction 0.92-1.0; Ldlr -1.1.

**Reading.** A cell-intrinsic ISR + p53 programme with pro-survival BCL-XL and small adaptive arms; no evidence that cytosolic mtDNA triggers innate sensing inside OLs. For the bench: ISR (PERK / ISRIB) and p53 are the pathways to perturb, not STING.

**Caveats.** Pooled pseudobulk: with about half the OLs responding (§30) the per-cell change is larger; the ranking matters, not the size. Pathway scores near 0.3 are still "significant" by permutation (cGAS-STING p = 0.024, apoptosis p = 0.012) because controls are tight; "absent" rests on the core sensors being flat and the ISG rises being spill-over. Transcript level only: STING activation is post-translational.

**Corrected 2026-10-21 (analysis/15):** ATF4 and p53 scores hold with any one responder left out; all 11 pathway tests keep q < 0.10.

### 41. OPCs do not mount the response: their ISR signal is at the spill-over level of astrocytes
*Notebook:* analysis/13-ol-mechanism-myelin-rna.ipynb (part B) · *Code:* `scripts/export_cell_counts.py`, `utils/composition.py` · *Results:* `results/ol_mechanism/tables/isr_by_lineage_stage.csv` · *Figure:* `results/ol_mechanism/figures/O2_opc_isr.png`

**Rationale.** The construct is Plp1-driven, so myelinating OLs should respond but Pdgfra+ precursors perhaps not. If OPCs are spared, the precursor pool could replace damaged OLs (related to the no-repair finding §23).

**Technical approach.** Per-cell raw counts (`nsc_ol_lineage_counts.csv.gz`, `scripts/export_cell_counts.py`). Three groups: pure OPCs (Immature OL I with Pdgfra > 0 and Mog = 0), astrocytes (which should not respond, so they set the level expected from segmentation spill-over of stressed-OL RNA), mature OLs incl. DA-OLs. Per cell: ISR transcripts (Trib3, Atf5, Ddit3, Gdf15, Sesn2, Cebpg, Cdkn1a) per 1,000 counts; mean per mouse, divided by the control mean. Exact permutation tests.

**Results.** Responder mean (control max): mature OLs 4.21x (1.12; per mouse 3.18-5.02; p = 0.012); OPCs 1.19x (1.07; 1.11-1.36; p = 0.024); astrocytes 1.21x (1.08; 1.15-1.32; p = 0.012). Non-responders 1.01-1.06 (non-responder RB4403 OPCs 1.22, above two responders).

**Reading.** OPCs rise no more than astrocytes, i.e. to the spill-over level; the response is confined to myelinating OLs.

**Caveats.** The astrocyte level is a proxy for spill-over, not a model of it (OPCs may sit in different neighbourhoods). OPC ISR is mildly but consistently up (analysis/15: p = 0.024), so a small own response cannot be ruled out. Whether OPCs replace lost OLs needs a later time point.

### 42. Counting myelin RNA without segmentation: Mog and Opalin fall ~40-55 % per OL in every region while Sox10 stays, so OLs lose myelin genes rather than disappearing
*Notebook:* analysis/13-ol-mechanism-myelin-rna.ipynb (part C) · *Code:* `scripts/export_myelin_transcripts.py`, `utils/composition.py` · *Results:* `results/ol_mechanism/tables/myelin_tx_tests.csv`, `results/ol_mechanism/tables/myelin_tx_by_region.csv` · *Figure:* `results/ol_mechanism/figures/O3_myelin_tx_heatmap.png`, `results/ol_mechanism/figures/O4_myelin_tx_dots.png`

**Rationale.** Myelin-gene loss (finding 2, §18) was measured on segmented cells, but much myelin mRNA lies in sheaths and processes outside every cell outline. This run asks whether the loss is real when every transcript is counted, and whether OLs disappear or only lose the programme.

**Technical approach.** `scripts/export_myelin_transcripts.py` reads each section's raw Xenium transcripts.parquet (qv >= 20) for 13 genes (myelin: Mbp, Mog, Mag, Cnp, Opalin, Cldn11, Ermn; lineage: Sox10, Olig2; stress: Trib3, Cdkn1a; controls: Gfap, Snap25) and bins them on a 20 um grid, with counts inside and outside cells. Each bin gets the region of the nearest cell centroid (bins > 30 um from any cell dropped). Per section x region: transcripts per mm2 of region (50 um occupied area) and per OL (OL + DA-OL count). 5 regions x 13 genes x 3 metrics, responders vs controls by exact permutation. Snap25 (neuronal) is a negative control.

**Results.** Share of transcripts outside any cell: Mbp 0.77, Gfap 0.66, Snap25 0.39, Cnp 0.34, Mog 0.20, Trib3 0.18. Per OL (log2 responders / controls): Mog -0.72 to -0.91 (all 5 regions p <= 0.036), Opalin -0.88 to -1.15 (p 0.012-0.024; cortex 0.060), Cldn11 -0.70 to -0.78 (all p = 0.012), Ermn -0.39 to -0.66, Cnp -0.31 to -0.41, Mag -0.33 to -0.65 (mostly n.s.), Mbp -0.25 to -0.49 (p 0.036-0.12). Per mm2 larger: Mog -0.77 to -1.21, Opalin -0.90 to -1.41. Sox10 per OL -0.07 to +0.02 (p >= 0.42); Olig2 per OL +0.04 to +0.67 (fibre tracts p = 0.048). Gfap per mm2 +1.97 (fibre tracts), +2.48 (thalamus), +2.86 (hypothalamus). Snap25 per mm2 -0.03 to +0.54 (p >= 0.07). Trib3 per OL +1.65 to +4.32. Share of Mog outside cells rises slightly (log2 +0.24 to +0.35).

**Reading.** The cell-based myelin loss was not a segmentation artefact; OLs keep lineage identity and lose the late myelin programme. Mbp lags (long-lived in sheaths). Gfap shows strong astrocyte reactivity that cell-based scores understate.

**Caveats.** Region labels are the original map (§43). Sox10 per mm2 is down in cortex and fibre tracts (-0.39, -0.26, p = 0.024), consistent with the modest OL density drop of §34. Olig2 per OL is not strictly "unchanged" in fibre tracts (+0.67, p = 0.048). Snap25 per OL rises (cortex +0.46, p = 0.036) because the OL denominator falls, so the per-mm2 Snap25 is the clean control.

**Corrected 2026-10-21 (analysis/15):** Mog per OL holds with any responder left out; Opalin is consistent but weaker (worst drop-one p = 0.071); Mbp per OL is carried by RB4405 and RB4627 (RB4401 inside the control range). Transcripts-per-OL family: 45 of 65 tests keep q < 0.10.

### 43. The original region map was distorted by the disease: a map from 400 disease-free anatomy genes puts the responders' pallidum back
*Notebook:* analysis/14-region-map-disease-free.ipynb (parts 1-4, 6) · *Code:* `utils/region_transfer.py`, `scripts/export_cell_counts.py` · *Results:* `results/region_map/tables/anatomy_genes.csv`, `results/region_map/tables/held_out_accuracy.csv`, `results/region_map/tables/region_coverage_both_maps.csv`, `results/region_map/tables/original_fibre_tracts_relabelled.csv`, `results/region_map/tables/relabelled_share_by_region.csv`, `results/region_map/region_transfer_labels.csv.gz` · *Figure:* `results/region_map/figures/G1_coverage_both_maps.png`, `results/region_map/figures/G3_maps_both.png`

**Rationale.** §34 found pallidum nearly absent and fibre tracts enlarged in all three responders. The region labels (§2) come from clustering neighbourhood-pooled transcriptomes with all genes, so stressed OLs and microglia could pull a region into another cluster. Every region-level result (§8, §13, §34) rests on these labels.

**Technical approach.** 1: anatomy genes = top 400 by an F score (between-region variance / mean within-region variance of log2 CPM) in control region pseudobulks, mean CPM > 20, excluding 469 response genes (all programme sets, mechanism sets, microglial markers from §37, OL driver genes from §15, region-consistent own changes from §37). 2: raw counts of these genes for every cell (`export_cell_counts.py --genes-file`), summed with each cell's 20 nearest neighbours as in the original pipeline, log-normalised. 3: PCA (30 PCs, fitted on 150k control cells) + multinomial logistic regression trained on control cells with their original labels ("Unknown" excluded); checked by leaving one control section out; then applied to every section. 4: compare region shares, and where responder "fibre tract" cells go.

**Results.** Held-out accuracy 0.904-0.946 per control section. Pallidum share, original -> disease-free: RB4405 0.04 -> 3.87 %, RB4627 0.47 -> 3.34 %, RB4401 0.68 -> 1.78 %; controls 2.85-4.23 % on the new map; responders vs controls p 0.012 -> 0.31. Fibre tracts in responders 8.6-10.8 % -> 5.4-7.7 % (controls 5.9-8.3 %; p 0.012 -> 0.61). Of responder cells originally called fibre tracts, 65.8 % stay fibre tracts and 22.9 % become pallidum (controls 91.0 % and 4.2 %). Outside pallidum and fibre tracts, 1.2-6.0 % of cells are relabelled per main region in responders (olfactory areas 9.8 %), similar to controls (vasculature and meningeal border more; the notebook text says 2-6 %); non-responder RB4403 (over-segmented) differs more in hypothalamus and olfactory areas. Maps: in responders the original map paints pallidum (beside the internal capsule) as fibre tract; the new map restores it in place and size.

**Reading.** The disease, not section plane, moved the white-matter-rich pallidum into "fibre tracts". RB4401 keeps somewhat less pallidum, which may be plane. New region-level work should use `region_new`.

**Caveats.** The reference is the original control labels, so any error in them carries over; 5-10 % held-out disagreement sits mostly at region borders. Exclusion of response genes depends on earlier lists; a response gene missed by them could still enter. The anatomy count export needs the h5ad (gitignored, 83 MB).

### 44. On the disease-free map, OL stress stays highest in thalamus and white matter (pallidum joins the top) and thalamus microgliosis holds
*Notebook:* analysis/14-region-map-disease-free.ipynb (part 5) · *Code:* `utils/region_transfer.py`, `utils/composition.py`, `utils/mosaic.py` · *Results:* `results/region_map/tables/regional_results_both_maps.csv` · *Figure:* `results/region_map/figures/G2_regional_results_both_maps.png`

**Rationale.** §43 showed the region map was distorted in responders; the regional claims (OL stress by region §13, microgliosis §34) must be re-checked on the corrected labels.

**Technical approach.** For both maps, per mouse and region (8 regions): share of OLs above the same-slide control 90th ISR percentile (leave-one-out percentiles, §29) and microglia per mm2 (`region_density`). Responders vs controls by exact permutation; control min-max shown against each responder.

**Results.** OL share above the control 90th percentile, disease-free map: thalamus 28-36 % vs 8-11 % (p = 0.012), fibre tracts 19-33 % vs 8-12 % (p = 0.012), pallidum 20-35 % vs 8-13 % (p = 0.012), cortex, striatum, hypothalamus 17-28 % (all p = 0.012), hippocampus 17-21 % (p = 0.036), olfactory areas p = 0.071 (RB4401 8 %). Microglia per mm2: thalamus 90-154 vs 60-82 (p = 0.012 on both maps); hypothalamus 75-117 vs 50-80 (p 0.060 -> 0.036); fibre tracts 63-141 vs 51-70 (p 0.012 -> 0.036, RB4401 inside the control range); pallidum RB4405 212, RB4627 117, RB4401 62 vs 62-105 (p = 0.095); cortex, striatum, hippocampus, olfactory p >= 0.36. On the original map responder pallidum had fewer microglia than controls (53-67 vs 73-100, p = 0.024), an artefact of the few cells left there.

**Reading.** The regional OL ranking survives and pallidum joins the top; microgliosis is robust in thalamus, and "fibre tracts" should read "white matter including the pallidum". The label-free result (§35) is unaffected.

**Caveats.** The notebook summary says microgliosis is "strongest in thalamus, pallidum and hypothalamus"; pallidum is not significant (p = 0.095, RB4401 inside the control range), so that phrasing overstates it. Hypothalamus microgliosis is carried by two responders (analysis/15, worst drop-one p = 0.18). The region DE of §8 used the old labels; re-run in §49.

**Corrected 2026-10-21 (analysis/15):** regional results on the new map: 10 of 16 tests keep q < 0.10; thalamus microglia density (new map) holds with any one responder left out.

### 45. The three responders share a response in microglia, a cell type never used to call them, and in OL genes left out of the call
*Notebook:* analysis/15-robustness.ipynb (part A) · *Code:* `utils/robustness.py` · *Results:* `results/robustness/tables/heldout_signature_r.csv` · *Figure:* `results/robustness/figures/R1_heldout_signature.png`

**Rationale.** Responders were called from OL stress (ISR) genes (§10) and later described with OL driver genes (§15) and the ATF4 / p53 pathway sets (§40). A reviewer could argue that the responder group is defined by the very genes later reported as its signature. This run asks whether the three responders still look alike on data that played no part in the call.

**Technical approach.** Pseudobulk counts per mouse from the analysis/01 region x time run (§8; region "All regions"), per cell type: OLs (Oligodendrocytes + DA-Oligodendrocytes) with a "selection" gene set removed (prespecified ISR genes, all data-derived driver genes of §15, the ATF4 / ISR and p53 target sets), and microglia, astrocytes and endothelial cells with the same genes removed. Genes kept if mean CPM in controls > 20; log2(CPM + 1). For each mouse, the profile is its log2 expression minus the mean of the same-slide (same-age) controls, excluding the mouse itself from that baseline. The reference signature is the mean profile of the *other* responders. Statistic: Pearson r over genes between the two (`rb.heldout_signature_r`). Because controls share the control baseline, their r is not zero and they form the null. Responders vs controls by exact permutation test (3 vs 6 mice, minimum p = 0.012).

**Results.** Microglia: responders r = 0.57 (RB4627), 0.66 (RB4401), 0.82 (RB4405); controls -0.07 to 0.42; non-responders -0.21 to 0.27. Complete separation, p = 0.012. OLs with selection genes removed: responders 0.57, 0.37, 0.75; controls -0.03 to 0.37 (RB4401 ties the top control); p = 0.036. Astrocytes: responders 0.25-0.44, controls 0.07-0.35. Endothelial cells: responders 0.34-0.43, controls -0.06 to 0.52 (no separation).

**Reading.** The shared microglial profile is the cleanest evidence that the three responders share a real response that is not an artefact of how they were selected. In OLs, beyond the stress genes, RB4405 and RB4627 stand out and RB4401 is weak, consistent with RB4401 being the mildest responder.

**Caveats.** n = 3 responders; r is computed over thousands of correlated genes, so the per-mouse values are not independent tests. The pseudobulk uses the original (pre-§43) region labels, though "All regions" is unaffected by region calls. The microglia result could still reflect a section-level technical factor shared by the three responders; the slide (age) baseline handles slide effects only.

### 46. Most core claims hold with any one responder left out; Mbp loss, hypothalamus microgliosis and the OL non-ISR signature rest on two responders
*Notebook:* analysis/15-robustness.ipynb (part B) · *Code:* `utils/robustness.py` · *Results:* `results/robustness/tables/leave_one_responder_out.csv`, `results/robustness/tables/metric_z_by_mouse.csv` · *Figure:* `results/robustness/figures/R2_metric_matrix.png`

**Rationale.** All per-mouse results compare 3 responders with 6 controls, so one unusual responder could carry a claim. This run re-tests the key per-mouse metrics of the README with each responder dropped in turn.

**Technical approach.** 17 per-mouse metrics drawn from earlier runs: OL ISR score (§10, the selection metric), OL top-1 % excess (§29), OL and microglia held-out signature r (§45), ATF4 / ISR and p53 pathway scores (§40), microglial DAM and astrocyte reactive scores, microglia density in thalamus and hypothalamus on the disease-free map (§44), microglia density in 30-50 % OL squares (label-free, §35), proliferating microglia obs/exp (§35), Mog, Opalin and Mbp transcripts per OL (segmentation-free, §42), and two expected-null metrics (Sox10 per OL; OPC ISR, §41). For each: exact permutation p of 3 responders vs 6 controls; p with each responder dropped (2 vs 6, minimum 1/28 = 0.036), reporting the worst; whether each responder lies beyond the control range in the expected direction; non-responders vs controls p. For the heat map, each mouse's value as a z-score against same-slide controls using the pooled within-slide control SD (`rb.slide_z`).

**Results.** With any one responder dropped p stays at the 0.036 floor, and all three responders lie beyond the control range, for 10 of the 15 metrics with an expected change: OL ISR score, OL top-1 % excess, microglia held-out r, ATF4 / ISR score, p53 score, microglial DAM, astrocyte reactivity, thalamus microglia density, OL-rich-square microglia density, Mog per OL. Weaker: proliferating microglia (p 0.024, worst drop-one 0.071) and Opalin per OL (0.012, 0.071), with all three responders beyond controls. Carried by two responders: OL held-out r (0.036; 0.107; RB4401 inside the control range), Mbp per OL (0.048; 0.143; RB4401 inside) and hypothalamus microglia density (0.036; 0.179; RB4627 inside). Expected nulls: Sox10 per OL p = 0.893. OPC ISR p = 0.024 (worst drop-one 0.107), interpreted as spill-over level (§41). Non-responders vs controls: all p >= 0.274. The notebook notes RB4403's high label-free microglia density (z = 4.1) as an over-segmentation effect.

**Reading.** The core mouse claims do not hinge on a single animal. RB4401 is the mildest responder in nearly every readout.

**Caveats.** With n = 3 the test cannot go below p = 0.012, and with 2 vs 6 not below 0.036, so "holds" means "reaches the floor". Metrics are correlated (several are readouts of the same OL response), so 10 of 15 is not 10 independent confirmations. The "finding" labels in the table refer to README finding numbers, not to this catalogue.

### 47. Multiple-testing control: pathway, myelin-transcript, regional and key-metric families survive; the 42-test density screen does not
*Notebook:* analysis/15-robustness.ipynb (part C) · *Code:* `utils/robustness.py` · *Results:* `results/robustness/tables/multiple_testing.csv`

**Rationale.** Analyses 11-14 ran many exact tests, each with a floor of p = 0.012. The README had reported nominal p-values. This run applies false-discovery-rate control to see which claims survive.

**Technical approach.** Benjamini-Hochberg (`rb.bh`) within five test families taken from earlier result tables: cell densities, 7 types x 6 regions (§34; 42 tests); regional results on the disease-free map (§44; 16 tests); transcripts per OL, 13 genes x 5 regions (§42; 65 tests); OL pathway scores (§40; 11 tests); and the 17 key metrics of §46. Also BH across all 151 tests pooled (`q_all`). Unit: mouse-level exact p, responders vs controls.

**Results.** Tests with q < 0.10 within family: densities 0 of 42 (6 had p < 0.05; the two at the floor, microglia in thalamus and fibre tracts, reach q = 0.25); OL pathway scores 11 of 11; transcripts per OL 45 of 65 (44 with p < 0.05); regional results on the new map 10 of 16; key metrics 16 of 17. When all 151 tests are pooled, the six nominal density hits get q_all 0.032-0.066 (all < 0.10).

**Reading.** The microgliosis claim cannot rest on the §34 density screen alone. It rests on its replications in the label-free OL-rich squares (§35) and the disease-free map (§44), both q < 0.10, and on §46. The modest OL density drops in cortex and fibre tracts, and the neuron changes in that screen, should be treated as nominal.

**Caveats.** Results depend on how families are drawn: pooled BH rescues the density hits that within-family BH does not. The notebook reports the stricter within-family result as the conclusion. The families do not cover every test in the project (e.g. region x age DE has its own padj; §8).

### 48. Reproducibility audit: every analysis notebook ran without errors, and every gitignored input is listed with how to regenerate it
*Notebook:* analysis/15-robustness.ipynb (part D) · *Results:* `results/robustness/tables/reproducibility.csv`

**Rationale.** Many notebooks read large per-cell exports that are gitignored. A fresh clone must know which inputs are missing and how to rebuild them.

**Technical approach.** For each analysis notebook: count error outputs and code cells without an execution count; list gitignored files under results/ (`git ls-files --others --ignored`) whose names occur in the notebook's code; take the regeneration instruction from the comment line above the matching .gitignore pattern.

**Results.** 15 notebooks (01-15), 0 error outputs in each. 01, 03, 05 and 07 need no gitignored inputs. The others need some of: `cells.csv.gz`, `cell_gene_scores.csv.gz`, `ol_gene_scores.csv.gz` (scripts/export_cell_table.py), EAE scores (scripts/export_eae_pseudobulk.py), region x time pseudobulk counts (scripts/run_region_time_de.py), NSC / OL-lineage counts and anatomy counts (scripts/export_cell_counts.py), binned myelin transcripts (scripts/export_myelin_transcripts.py, office Mac, raw Xenium). analysis/15 itself showed 6 unexecuted code cells, because the scan ran while the notebook was executing.

**Reading.** A clone can open every notebook with its saved figures. Re-running analysis/02, 04, 06 and 08-15 needs the exports, most made from the h5ad or raw Xenium folders on the office Mac.

**Caveats.** "Unexecuted cells" for analysis/15 is a self-scan artefact, not a real gap. Analyses 16 and 17 were added later and are not in this table; their inputs (`results/human_ms/human_ol_*`, `results/human_xenium/cells.csv.gz`, `zones.csv.gz`) come from scripts/export_human_ms_ol.py and scripts/export_human_xenium_cells.py, and the human Xenium tables are gitignored as unpublished patient data. No check that re-running reproduces the same numbers was made.

### 49. Region x age DE re-run on the disease-free region map: totals and fold changes hold; striatum, fibre-tract and pallidum hit lists change
*Notebook:* analysis/15-robustness.ipynb (part E) · *Code:* `scripts/run_region_time_de.py`, `utils/region_time_de.py` · *Results:* `results/region_time_disease_free/summary_by_region_age.csv`, `results/region_time_disease_free/de_condition_by_region_age.csv`, `results/robustness/tables/region_de_old_vs_new.csv`, `results/robustness/tables/region_de_lfc_correlation.csv` · *Figure:* `results/robustness/figures/R3_region_de_old_vs_new.png`

**Rationale.** §43 showed that the original transcript-based region map was distorted by the disease (responder pallidum called fibre tracts). The region x age DE pipeline (§8) used that map, so its regional hit lists needed re-checking.

**Technical approach.** The §8 pipeline (pyDESeq2 0.5.4, pseudobulk per cell type x region x mouse, all six mtDSB mice vs controls per age) re-run with `--region-labels results/region_map/region_transfer_labels.csv.gz` (column region_new, from §43). Same settings: >= 20 cells per pseudobulk, "Unknown" excluded, padj < 0.05, |log2FC| >= 0.5, baseMean >= 20. 980,474 cells x 5,101 genes; 1,283 pseudobulk samples kept, 588 dropped (old run 1,264 / 585). Comparison: number of significant genes per cell type x region x age, and Pearson r of log2FC over shared genes per cell type x region.

**Results.** Total significant genes 2,024 (old) vs 2,007 (new). Median log2FC r across cell types: cortex 0.95, hypothalamus 0.95, thalamus 0.94, hippocampus 0.93, striatum 0.92, ventricular system 0.92, fibre tracts 0.87, pallidum 0.62, vasculature 0.58; all regions pooled 0.99. Largest shifts: astrocytes striatum P60 29 -> 1; endothelial striatum P60 10 -> 1; microglia hypothalamus P60 12 -> 24; OLs thalamus P60 4 -> 7; OLs hypothalamus P21 13 -> 10.

**Reading.** The region x age results of §8 stand overall. For striatum, fibre tracts and pallidum, use the new run; the old striatum hits partly came from mislabelled tissue.

**Caveats.** Like §8, this pools responders and non-responders as "mtDSB", so effects are diluted. `run_info.json` of the new run does not record that region labels were replaced (its keys still name RBD_compartment_simplified); only `run.log` does. The age contrast and interaction tables (§9) were also regenerated but not compared here.

### 50. Human MS OLs lack the mtDSB-specific markers but share a CHOP / p21 stress arm, alongside an EAE-like interferon-gamma programme
*Notebook:* analysis/16-human-ms-signature.ipynb (parts A, B) · *Code:* `utils/human_ms.py`, `scripts/export_human_ms_ol.py` · *Results:* `results/human_ms/tables/pathway_scores_lesion_vs_control.csv`, `results/human_ms/tables/pathway_scores_per_sample.csv`, `results/human_ms/tables/core_genes_log2fc.csv` · *Figure:* `results/human_ms/figures/H1_pathway_heatmap.png`, `results/human_ms/figures/H2_core_genes.png`

**Rationale.** In mice, the mtDSB OL response is an ATF4 / ISR programme with mtDSB-specific markers (Trib3, Gdf15, Atf5, Eif4ebp1) plus a p53 / p21 arm and little interferon-gamma, while EAE OLs show the reverse (§31, §40). The question is whether OLs in human MS lesions show the mtDSB pattern.

**Technical approach.** Four public MS single-nucleus RNA atlases, harmonised copies on the office Mac: Jäkel 2019, Schirmer 2019, Absinta 2021, Lerma-Martin 2024. `export_human_ms_ol.py` keeps OL and OPC nuclei (141,213: 130,897 OLs, 10,316 OPCs) and sums raw counts per (sample, lesion type) for human orthologs of the mouse panel genes. Orthologs by upper-casing plus a few exceptions (`to_human`). Gene sets: MTDSB_CORE (TRIB3, GDF15, ATF5, DDIT3, CDKN1A, EIF4EBP1, ALDH1L2, SESN2, BBC3), ISR, P53, driver genes up/down (§15), IFNG, IMMUNE_OL, MYELIN. OL pseudobulks with >= 30 nuclei, log2(CPM + 1), genes with mean CPM > 1; each sample minus the mean of that atlas's control samples; set score = mean over set genes. Lesion vs control samples by exact permutation test. Absinta lesion core (19 OL nuclei) dropped by the threshold. Mouse references from the §31 tables.

**Results.** Mouse reference (OL log2FC): mtDSB responders ATF4 / ISR 1.85, IFN-gamma 0.32; EAE lesion 0.06 and 3.11. MTDSB_CORE in lesion OLs: Schirmer +0.91 (A/CA, p 0.022) and +1.17 (CI, p 0.027); Absinta +0.69 to +1.04 (p 0.057-0.20); Lerma-Martin +1.16 (CA, p 0.004) and +0.82 (CI, p 0.024). IFNG rises as much or more: Lerma-Martin +1.76 / +1.87 (p 0.002 / 0.005), Absinta CA edge +1.91 (p 0.029), Schirmer A/CA +1.31 (p 0.011). Jäkel is inconsistent (core -0.92 in CA, -1.28 in CI). Gene by gene, the specific markers barely move in Absinta and Lerma-Martin: TRIB3 -0.10 to +0.31, GDF15 -0.11 to +0.27, EIF4EBP1 -0.28 to -0.02, ATF5 -0.80 to -0.18 (mouse responders +2.07 to +5.41). Schirmer: TRIB3 +0.85 / +1.14, ATF5 +0.90 / +1.32. Shared stress genes rise: DDIT3 +0.96 to +2.85 and CDKN1A +0.80 to +3.07 across Schirmer, Absinta and Lerma-Martin, with SESN2, GADD45A and CEBPG. IRF1 and PSMB8 rise as in EAE.

**Reading.** MS OLs look like inflammation-exposed OLs with a CHOP / p21 stress component, not like OLs with primary mtDNA damage. The core-set score rises mainly because the set includes general stress genes.

**Caveats.** snRNA (nuclear transcripts, low depth; TRIB3 and GDF15 are low-expressed), post-mortem tissue, mostly chronic lesions, 1-8 samples per lesion type. Samples, not donors, are the unit; whether one donor gives several lesion samples is not modelled. Per-lesion p-values are not corrected for multiple testing. The README finding-19 summary "within +/-0.3 log2" does not hold for ATF5 in Absinta (-0.54 to -0.80), though that is in the "not raised" direction.

**Next.** Test the markers in situ with the same technology as the mouse data (§53).

### 51. OLs carrying two or more mtDSB-specific transcripts are rare in all four MS atlases; DDIT3 + CDKN1A co-expressing OLs rise in lesions
*Notebook:* analysis/16-human-ms-signature.ipynb (part C) · *Code:* `scripts/export_human_ms_ol.py` · *Results:* `results/human_ms/tables/specific_marker_fraction_per_sample.csv`, `results/human_ms/tables/lerma_martin_subtypes.csv` · *Figure:* `results/human_ms/figures/H3_per_cell_fractions.png`

**Rationale.** Pseudobulks (§50) could dilute a small subpopulation of OLs with the mtDSB state, as in the mouse, where only a minority of responder OLs carry Trib3.

**Technical approach.** Per OL nucleus: "specific-positive" if the summed raw counts of TRIB3, GDF15, ATF5 and EIF4EBP1 are >= 2; "CHOP + p21" if DDIT3 > 0 and CDKN1A > 0. Fractions per (sample, lesion) with >= 30 nuclei. In Lerma-Martin, the published OL / OPC subtypes were scored with the gene-set scores (log1p of 10,000 x set counts / total) and their share per lesion type was tabulated.

**Results.** Specific-positive OLs (% per lesion type): Absinta 0.00-0.02 %, Lerma-Martin 0.06-0.09 %, Jäkel 0.00-0.20 % (control 0.20 %). Schirmer, the deepest libraries and including cortex: 1.82 % (A/CA) and 3.98 % (CI) vs 1.11 % in controls. CHOP + p21 OLs: Absinta 0.37-0.50 % vs 0.00 % in controls; Lerma-Martin 0.14-0.21 % vs 0.02 %; Schirmer CI 0.74 %, A/CA 0.08 % vs 0.00 %. Lerma-Martin subtypes: OL_Dis2 has the highest MTDSB_CORE score (0.431; next OPC_COP 0.401) and makes up 4.94 % of OL-lineage nuclei in CA, 2.09 % in CI and 1.27 % in controls. OL_Dis1 shifts more (22.7 % CA, 28.5 % CI vs 5.5 % controls) but has a lower core score (0.214).

**Reading.** No sizeable MS OL subpopulation carries the mtDSB-specific markers at snRNA depth. The core score of OL_Dis2 comes from shared stress genes.

**Caveats.** In mouse responders, Trib3+ OLs are about 1-3 % by Xenium (notebook text, citing analysis/13). snRNA detects such low-expressed genes less often, so a small population could be missed. Thresholds (>= 2 transcripts) were not varied. There are no statistical tests for the fractions or subtype shifts.

### 52. Genome-wide, MS lesion OL changes correlate only weakly with both mouse models, slightly more with mtDSB than with EAE
*Notebook:* analysis/16-human-ms-signature.ipynb (part D) · *Results:* `results/human_ms/tables/similarity_to_mouse_models.csv` · *Figure:* `results/human_ms/figures/H4_similarity.png`

**Rationale.** Gene sets (§50) test chosen genes. A genome-wide comparison asks which mouse model the MS OL response resembles overall.

**Technical approach.** Per atlas and lesion type: mean OL log2FC vs control samples (pseudobulk as in §50; genes with mean CPM > 5) correlated (Pearson r over genes) with three mouse OL profiles from §31: mtDSB responders, EAE lesion and EAE DAO vs MOL. Genes in common: 2,143-2,442 per atlas.

**Results.** r with mtDSB: Schirmer A/CA 0.170, CI 0.013; Absinta periplaque 0.117, CA edge 0.158, CI edge 0.149; Lerma-Martin CA 0.133, CI 0.095; Jäkel -0.044 to 0.079. r with EAE lesion: Schirmer 0.135 / -0.147; Absinta -0.000 to 0.101; Lerma-Martin 0.068 / 0.015; Jäkel 0.040-0.133. r with EAE DAO state: -0.336 to 0.133, mostly near 0.

**Reading.** In the three larger atlases MS OL changes resemble mtDSB slightly more than EAE: in Schirmer, Absinta and Lerma-Martin, r_mtDSB exceeds r_EAE lesion in all 7 lesion types, by 0.035-0.16. Per the notebook, the overlap with mtDSB is the shared CHOP / p21 arm and the EAE-like part is the interferon-gamma programme. All correlations are weak.

**Caveats.** There is no null distribution or confidence interval for the correlations. Species, technology (snRNA vs Xenium), tissue and time point all differ. Correlations over thousands of co-regulated genes overstate the effective sample size.

### 53. In human MS Xenium sections, OLs in demyelinated white matter raise TRIB3, EIF4EBP1, p21 and PUMA and lose MOG / OPALIN, while SOX10 stays flat
*Notebook:* analysis/17-human-xenium-lesion-ols.ipynb (parts A, B) · *Code:* `utils/lesion_zones.py`, `scripts/export_human_xenium_cells.py` · *Results:* `results/human_xenium/tables/cells_by_zone.csv`, `results/human_xenium/tables/log2_by_zone.csv`, `results/human_xenium/tables/cluster_markers.csv` · *Figure:* `results/human_xenium/figures/X1_zones.png`, `results/human_xenium/figures/X2_zone_heatmap.png`

**Rationale.** The snRNA atlases pool whole tissue blocks per sample (§50-§52), so a localised OL response could be diluted. An in-house human CNS Xenium 5k run (unpublished) allows a within-section comparison with the same technology as the mouse data.

**Technical approach.** Four sections, 274,472 cells: two MS sections (cross-sections that look like spinal cord or brainstem; called MS-1 and MS-2 here) and two Baló concentric sclerosis sections. `export_human_xenium_cells.py` exports per-cell raw counts of 38 marker genes, centroids and Leiden clusters. ATF5 is not on the human panel. OLs = Leiden (res. 1) clusters 0, 1, 2, 4 (SOX10+); OPCs = 18, 22 (PDGFRA+ OLIG2+), chosen from cluster marker means. There is no pathology annotation and no control tissue, so zones are drawn from the data (`lesion_zones.py`): per cell, MOG + MBP + MAG + CNP + MOBP counts per 1,000 counts among all cells within 100 µm; divided by the section's 90th percentile outside grey matter; < 0.25 demyelinated, < 0.6 edge, else myelinated. Cells with local SNAP25 > 1 per 1,000 are grey matter and excluded. Per section, OLs (and, as spill-over control, all non-OL / non-OPC "other" cells) are pooled per zone; expression per 10,000 counts; log2((zone + 0.05) / (myelinated + 0.05)); zones with >= 100 cells.

**Results.** OLs per zone: MS-1 3,235 demyelinated, 11,623 edge, 18,219 myelinated; MS-2 1,401 / 8,375 / 17,184. Baló-1: 6 / 26 / 198 (5,110 OLs in grey matter); Baló-2: 32 / 1,514 / 6,981. Demyelinated vs myelinated OLs, MS-1 / MS-2 (log2): TRIB3 +1.27 / +1.45, EIF4EBP1 +0.54 / +1.44, GDF15 +1.54 / +3.10, CDKN1A +1.80 / +2.27, BBC3 +1.66 / +2.40, SESN2 +0.64 / +1.04, GADD45A +0.53 / +0.93, DDIT3 -0.19 / +0.45, MOG -2.26 / -2.57, OPALIN -2.22 / -4.16, SOX10 -0.04 / +0.11, OLIG2 -1.37 / -2.45, IRF1 +1.02 / +1.18, CD74 +0.48 / +1.54. Edge zones are intermediate. Other cells in the same demyelinated zones: TRIB3 -0.19 / +0.40, EIF4EBP1 +0.20 / -0.07, CDKN1A +0.15 / +0.13, BBC3 +0.32 / +0.27, but GDF15 +2.83 / +1.87.

**Reading.** Within the same section, lesion OLs show the mtDSB pattern (specific markers and the p53 arm up, myelin genes down, SOX10 kept), and neighbours change TRIB3, EIF4EBP1, CDKN1A and BBC3 little. That argues the rises are OL-intrinsic. GDF15 also rises in neighbours (likely macrophages, per the notebook), so the OL GDF15 signal may partly be neighbours' RNA. Unlike mtDSB, lesion OLs also raise interferon / antigen-presentation genes.

**Caveats.** Two MS donors, no pathology annotation, no control tissue; Baló sections are mostly grey matter at this resolution. Zones are defined from myelin RNA, so the myelin-gene results for OLs are partly circular (stated in the module). Much of the demyelinated zone runs along the tissue border (see §54). OLIG2 falls in lesion OLs (-1.37 / -2.45), so "lineage kept" rests on SOX10 alone. OL clusters carry substantial GFAP signal in marker means (cluster 0: GFAP 49.1 vs MBP 28.2), and MOG / OPALIN also fall in "other" cells, so segmentation spill-over affects both groups. **Text vs table:** the notebook text (part B) gives MOG -2.2 / -2.3 and OPALIN -2.1 / -3.0; `log2_by_zone.csv` gives MOG -2.26 / -2.57 and OPALIN -2.22 / -4.16. The neighbour range "about -0.2 to +0.3" is slightly exceeded (TRIB3 +0.40 in MS-2, BBC3 +0.32 in MS-1).

### 54. Deep in the tissue, demyelinated OLs still carry p21 and PUMA about 3x and an mtDSB marker plus the p53 arm 4-10x more often; near the tissue border, myelinated OLs are just as high
*Notebook:* analysis/17-human-xenium-lesion-ols.ipynb (parts C, D) · *Code:* `utils/lesion_zones.py` · *Results:* `results/human_xenium/tables/rates_block_bootstrap.csv`, `results/human_xenium/tables/rates_by_zone_and_border_distance.csv`, `results/human_xenium/tables/rates_deep_tissue_block_bootstrap.csv` · *Figure:* `results/human_xenium/figures/X3_rates.png`, `results/human_xenium/figures/X3b_border.png`

**Rationale.** Zone pseudobulks (§53) cannot tell whether a few OLs express a lot or many express a little, and they have no uncertainty. Neighbouring cells are spatially correlated, so naive per-cell statistics would be overconfident. The demyelinated rim also lies along the tissue border, a possible artefact.

**Technical approach.** Per-cell events: >= 1 transcript of TRIB3, EIF4EBP1, GDF15, CDKN1A or BBC3, and a combined event (TRIB3 or EIF4EBP1 > 0 and CDKN1A or BBC3 > 0). Rate ratio demyelinated (or edge) vs myelinated within a section, with a 95 % CI from a spatial block bootstrap (`block_bootstrap_ratio`: 250 µm squares resampled with replacement, 300 resamples). This keeps neighbouring cells together, so spatial correlation is respected. Same for other cells as the spill-over control. Border distance per cell from a 25 µm occupancy grid (closed, holes filled). OLs binned < 150, 150-400 and > 400 µm from the border; ratios recomputed for OLs > 400 µm only.

**Results.** All OLs, demyelinated vs myelinated, MS-1 / MS-2: TRIB3+ 2.76x (CI 1.76-4.37) / 3.50x (1.83-6.45); CDKN1A+ 2.77x (2.44-3.25) / 2.85x (2.29-3.51); BBC3+ 3.37x (2.28-4.86) / 4.79x (3.28-7.35); EIF4EBP1+ 1.24x (0.92-1.65) / 2.17x (1.47-3.05); GDF15+ 7.9x / 24.5x; combined 3.02x (1.95-4.76) / 5.08x (2.57-9.29). Absolute rates are low, for example TRIB3+ 0.71 % vs 0.26 % of OLs in MS-1. Other cells: TRIB3 1.06x / 1.27x (CIs include 1), CDKN1A 1.29x / 1.08x, combined 1.92x (1.67-2.21) / 0.97x; GDF15 10.5x / 7.5x. Border: myelinated OLs within 400 µm of the border carry CDKN1A in 13.5-15.0 % (MS-1) vs 4.09 % deeper than 400 µm. Deep OLs only (> 400 µm; 1,107 and 578 demyelinated vs 16,889 and 15,178 myelinated): CDKN1A+ 3.14x (2.56-3.76) / 3.40x (2.44-4.24); BBC3+ 3.24x (1.81-5.56) / 3.33x (1.43-5.97); combined 3.71x (1.58-7.07) / 10.50x (3.05-19.65); TRIB3+ 2.18x (0.63-4.68) / 3.18x (0.69-6.60). Baló-2 edge vs myelinated OLs: CDKN1A 3.67x (2.87-4.58), EIF4EBP1 3.91x (2.17-6.37), BBC3 3.29x (1.27-9.93), combined 7.68x (1.91-29.93), TRIB3 0.92x (0-2.98).

**Reading.** The border is a real confound: either an edge artefact or subpial pathology, and these data cannot tell which. Removing it leaves a clear p21 / PUMA and combined-signature excess in lesion OLs. TRIB3 alone points the same way but is too sparse to be significant deep in the tissue.

**Caveats.** Two MS sections; each CI covers within-section spatial variation only, not donor-to-donor variation. 300 bootstrap resamples gives coarse CI tails. The block size (250 µm) and border cut-offs (150 / 400 µm) were not varied. The combined event in other cells is also up about 2x in MS-1. Baló numbers rest on few white-matter OLs (Baló-1 has 6 demyelinated OLs).

### 55. Away from the tissue border, the share of stressed OLs rises as local myelin RNA falls
*Notebook:* analysis/17-human-xenium-lesion-ols.ipynb (part E) · *Code:* `utils/lesion_zones.py` · *Results:* none saved (per-cell zones in the gitignored `results/human_xenium/zones.csv.gz`) · *Figure:* `results/human_xenium/figures/X4_gradient.png`, `results/human_xenium/figures/X5_maps.png`

**Rationale.** Three discrete zones (§53) hide whether the OL stress state tracks the degree of myelin loss. Restricting to deep tissue keeps the border confound (§54) out.

**Technical approach.** White-matter OLs > 400 µm from the border, binned by local myelin RNA relative to the section's myelinated level (bins < 0.15, 0.15-0.25, 0.25-0.4, 0.4-0.6, 0.6-0.8, 0.8-1, > 1); share of OLs with TRIB3+, CDKN1A+, BBC3+ and the combined event per bin, bins with >= 50 OLs, per section. Maps: zones with OLs carrying the combined event overlaid, for the two MS sections.

**Results.** No table or statistic was saved; the result exists only as figures. The notebook states that in both MS sections the share of OLs carrying TRIB3, CDKN1A or BBC3 rises as local myelin RNA falls, and that OLs carrying both arms cluster in the demyelinated rim and the adjacent edge zone, with few in the myelinated interior.

**Reading.** A graded relation between myelin loss and OL stress fits the mtDSB picture (§42, §40), though it does not show that mtDNA damage causes it.

**Caveats.** Purely descriptive: no trend test, no CI per bin. Myelin RNA defines the x-axis and is partly made by the same OLs (circular for OLs losing myelin RNA). Two donors. The maps include cells near the border, where §54 shows the stress transcripts are confounded.

**Next.** Pathology-annotated lesions (LFB / MOG IHC, lesion activity, pia vs cut edge) on more sections and donors, with control tissue.
