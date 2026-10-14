# Findings: oligodendrocyte response to mtDSB at 42 days after induction

Xenium 5k, 12 sections (one per animal): DOX removal at P21 or P60 x control / mtDSB x 3, all collected
42 days after induction. Status: 2026-10-14, full `novae4` data. Every number below is reproduced by
`analysis/02-ol-response-and-classifiers.ipynb`. **Fig** = `figures/`, **Tab** = `tables/` in this folder.

All models use *general* cell-type labels: the disease-defined cluster `DA-Oligodendrocytes` is folded into
`Oligodendrocytes`, so the treatment effect has to be found by the readouts, not by the clustering.

---

## Headline

**At 42 days after induction, mtDSB produces a strong oligodendrocyte stress response in 3 of the 6
mtDSB animals (RB4405, RB4401, RB4627) and no detectable response in the other 3 (RB4403, RB4630,
RB4653). In responders the response affects oligodendrocytes in every brain region, most strongly in
fiber tracts and the thalamus.**

---

## Claims and evidence

### 1. The mtDSB group splits into responders and non-responders on a predefined stress signature
*Fig F1, Tab `sections_with_responder_calls`, `isr_per_gene_whole_section`, `isr_group_tests`*

* Readout: share of OL transcripts from 10 canonical ATF4 / mitochondrial integrated-stress-response
  targets (Atf4, Atf5, Ddit3, Trib3, Gdf15, Sesn2, Fgf21, Psat1, Cth, Eif4ebp1). The gene set was fixed
  from the literature before looking at per-section values; it uses no clustering and no model.
* Three mtDSB sections score 1.35-1.57x the *highest* control (3477-4031 vs 2039-2568 CPM); the other three
  score within the control range (2059-2347).
* Per gene, responder and control ranges do not overlap for Atf4, Atf5, Ddit3, Trib3, Sesn2 and Eif4ebp1
  (e.g. Trib3 77-214 vs 1-7 CPM; Atf5 503-741 vs 136-172 CPM). Non-responders sit inside the control range
  on every gene.
* Test without any responder selection: the variance of the (log) ISR score is ~10x larger in mtDSB than
  in controls (Brown-Forsythe p = 1.0e-4; F-test p = 0.013), while the group means do not differ
  significantly (Mann-Whitney one-sided p = 0.15). *The mtDSB group is heterogeneous, not uniformly shifted.*

**Strength: strong.** Prespecified readout, non-overlapping ranges on six genes, formal variance test.

### 2. Two further readouts independently identify the same three animals
*Fig F2, Tab `three_readouts_by_section`, `three_readouts_spearman`*

* **Single-cell classifier** (oligodendrocytes, general labels; each section predicted by a model trained
  on the other 11): the three responders are the three highest-scoring sections of all 12 (mean held-out
  p(mtDSB) 0.61-0.81 vs controls 0.38-0.60 and non-responders 0.28-0.37).
* **Cluster share** (descriptive): 22-65% of responders' OLs fall in the `DA-Oligodendrocytes` cluster vs
  0.3-2.3% in every other section.
* Spearman correlation across the 12 sections: ISR vs classifier 0.71, ISR vs cluster 0.56, classifier vs
  cluster 0.59. Top 3 by each readout = the same three animals. (Chance that the top 3 of 12 are all mtDSB
  under random labels: 20/220 = 0.09 for a single readout.)

**Strength: strong for "same animals"; the readouts are not fully independent** (all are built from OL
expression), but they use different inputs (prior genes / labels of other sections / clustering).

### 3. The responder signature is the same at both ages
*Fig F3B, Tab `ol_cell_model_section_predictions`*

* An OL model trained only on P60 sections ranks RB4627 (the P21 responder) as the most mtDSB-like P21
  section (0.65 vs 0.29-0.50); a model trained only on P21 sections ranks RB4405 and RB4401 as the top two P60
  sections (0.88, 0.85 vs 0.64-0.74).
* Responders: 2 of 3 at P60, 1 of 3 at P21. With 3 animals per cohort this split is uninformative (Fisher
  exact p = 1.0); the earlier impression "stressed OLs are increased in the P60 cohort" is carried by RB4405
  and RB4401.

**Strength: moderate** (transfer AUC on 6 sections; 3 vs 3 permutation floor is p = 0.05).

### 4. Why plain "mtDSB vs control" prediction fails, while age is easy
*Fig F3A, F7, Tab `pseudobulk_auc_whole_section`, `pseudobulk_permp_whole_section`, `cells_auc_whole_section`*

* Age cohort (P21 vs P60) is predictable from almost every cell type: 52 of 101 cell type x region strata
  at BH-FDR < 0.1 (pseudobulk, exact permutation tests); single-cell models reach AUC 1.0 for 10 of 12
  cell types.
* Condition (mtDSB vs control) is not: 0 of 101 strata reach permutation p < 0.05 pooled, within cohort or
  across cohorts. This follows from claim 1: half of the "mtDSB" labels belong to animals whose tissue looks
  like control.
* Leave-one-out AUC far below 0.5 (e.g. 0.00 within cohort) is a known artefact when there is no signal
  in the training labels, not a reversed effect.

### 5. Where: brain-wide in responders, strongest in fiber tracts and thalamus
*Fig F4, F5, Tab `region_summary_by_group`, `region_rank_in_responders`, `ol_frac_above_control_p90`*

* In responders every region is affected: ISR score 1.23-1.88x controls, classifier high-score fraction
  1.39-2.93x, cluster share 21-65% (controls 0-4%).
* Ranking regions within each responder and averaging over the three readouts: **fiber tracts and thalamus
  tie for the top** (mean rank 2.7 of 11), then vasculature, striatum, **hypothalamus (4.8)**, cortex;
  hippocampus and olfactory areas are least affected (8.7, 8.9).
* Spatial maps (F5): OLs scoring above the 90th percentile of control OLs make up 24-40% of OLs in
  responders vs 5-16% in controls and 4-6% in non-responders, concentrated in thalamus and fiber-tract arcs.
* The earlier observation "hypothalamus and thalamus are enriched for stressed oligodendrocytes" holds for
  the **thalamus**; the hypothalamus is clearly affected but mid-ranked. It is not explained by neuronal
  neuropeptide spill-over (control hypothalamus has 0.8% cluster share).

**Strength: moderate.** Regional differences are smaller than the responder / non-responder difference and
rest on 3 animals.

### 6. Non-response is not explained by sex, labels, genotype or section plane
*Fig F6, Tab `technical_checks`*

* Sex: responders 2 M + 1 F, non-responders 2 M + 1 F. Sex labels agree with Xist expression in all 12
  sections (no sign of a sample swap).
* Genotype: all six mtDSB animals are PlpTA:mtPst1. Note that control RB4620 is mtPst1-only (no tTA); the
  other five controls are PlpTA-only.
* Section plane: responders contain somewhat more fiber tract (8.9-10.8% of cells vs 5.8-9.5% for the
  others), but the response is present *within* every region, including cortex and olfactory areas
  (claim 5), and the non-responder RB4403 has the largest fiber-tract OL population of any section.
* RB4403 has ~2x the cells of other sections (135 k); its proportions are normal, and it is a
  non-responder on all three readouts.
* The panel has no probes for the transgenes (tTA, mtPstI), so induction cannot be measured in this data.

### 7. What the responding oligodendrocytes do: 34 genes up, 9 down
*`results/ol_gene_drivers/`: Fig G1-G4, Tab `driver_genes`, `cohort_weight_replication`,
`driver_genes_log2fc_by_section`, `driver_programme_by_region`; notebook `analysis/03-ol-stress-genes.ipynb`*

* **Driver rule (fixed before inspecting genes):** within each responder section, across individual OLs,
  the gene correlates with the cell's out-of-fold mtDSB score (r > 0.1 in each of the three responders) and
  not in controls or non-responders (mean r < 0.1); mirror rule for down genes. Being within a section, this
  is immune to batch, section plane and time of sacrifice.
* **Up (34):** the ATF4 / integrated stress response (Atf4, Atf5, Ddit3, Trib3, Cebpg, Sesn2, Yars, Vegfa);
  mitochondrial and protein quality control (Hspa9 = mortalin, Sqstm1 = p62, Lamp2, Nfe2l1, Gstp1, Mt2);
  cell cycle / survival (Cdkn1a = p21, Ccnd1, Myc, Bcl2l1); immune / reactive (MHC-I H2-D1 and B2m, Cd40,
  Serpina3n, Adamts1); lineage / immature (Olig2, Adgrg1); cytoskeleton (Flna, Msn, Coro1c, Mpzl1, Ptpra);
  lipid (Abca1). Six of the ten prespecified ISR genes are recovered without being given to the analysis.
* **Down (9):** mature-OL / myelin genes (Cldn11, Mog, Mal, Cnp, Ugt8a, Tppp, Stmn4, Aplp1, Mapk8ip1). With
  Olig2 / Adgrg1 up, stressed OLs appear to lose part of their mature, myelinating identity.
* **Per animal:** in OL pseudobulk, Cd40 goes from 1-4 CPM (controls, non-responders) to 49-165 (responders);
  Cdkn1a 6-31 -> 147-235; Serpina3n 34-59 -> 242-620; Mog 3741-6376 -> 1570-3778. Non-responders are flat.
* **Replication:** the classifier trained only on the 6 P21 mice and the one trained only on the 6 P60 mice
  (no animals in common) share 10 of their top-50 up genes (0.5 expected; hypergeometric p = 2e-11),
  including Atf5, Trib3, Mt2, Bcl2l1, and 4 of their top-50 down genes (Cldn11, Opalin, Il33, Vstm2a; p = 0.001).
* **Region:** the up-programme is highest in the thalamus in all three responders (mean log2 FC 1.4-1.8),
  then fiber tracts and hypothalamus; the myelin genes drop most in fiber tracts and thalamus.
* **Checked and excluded:** glucocorticoid-response genes (Sgk1, Fkbp5, Zbtb16, Plin4) carry large
  classifier weights but no within-section signal and no responder pattern (Sgk1 is highest in a control),
  so they reflect a section-level state such as stress / time at sacrifice, not mtDSB. Xist has a negative
  weight (sex imbalance in training folds). No driver is neuron-enriched (max neuron / OL ratio 1.8x), so
  segmentation spill-over does not explain them.

**Strength: strong for the programme** (within-section design, replication across disjoint animals,
recovery of the prior ISR genes). Individual gene ranks rest on 3 animals; the programme annotation is manual.

### 8. Consistent with the earlier Bayesian and DE analyses
*`docs/ANALYSIS_REVIEW.md`; `legacy/results/baysian_diff_condition.csv`*

* The hierarchical Bayesian case-control model (`legacy/notebooks-02/09-12`, OLs, ages pooled) tested 22 of the
  43 driver genes; **all 22 have the same direction** (17 up, 5 down; one-sided sign test p = 2.4e-7) and 12
  reached its FDR < 0.05 (Mt2, Gstp1, B2m, H2-D1, Atf4, Atf5, Olig2 up; Cldn11, Cnp, Mal, Mog, Ugt8a down). Its
  posterior shifts were small (0.01-0.05) because each estimate averaged responders with non-responders.
* Per age: P21 13 / 13 tested drivers in the same direction; P60 13 / 18.
* The earliest per-sample marker lists (`legacy/notebooks-01/06`) already show the driver genes in RB4405 and
  RB4627; the DA-OL sub-clustering found the same programme. The earlier "P21 modest / P60 strong" reading is
  explained by responder composition (claim 3).

---

## Hypotheses only (not claims)

* **Astrocytes and endothelial cells** carry weak cross-cohort signals (single-cell transfer AUC 0.78-1.0
  on 6 sections; endothelial P60 -> P21 ranks all three P21 mtDSB above all three P21 controls, including
  the OL non-responders), strongest in fiber tracts, vasculature and the meningeal border. Small margins,
  no permutation test at this level, many cell types tested. *Tab `cells_auc_glia_vascular_by_region`.*
* **Microglia** echo the responder pattern (the two P60 responders score highest).
* **Batch:** if mtDSB and control sections were on different Xenium slides, the weak signals above could
  be slide effects. The OL responder split cannot be (it divides animals within the same condition).

## What we cannot claim

* That mtDSB "increases stressed OLs" as a group effect: the group means do not differ (claim 1).
* An age (P21 vs P60) difference in susceptibility: 2/3 vs 1/3 responders is uninformative.
* Anything about time since induction: it is 42 days for every Xenium section.

## Open questions and suggested next steps

1. **Why do three animals not respond?** Induction / recombination efficiency is the leading hypothesis.
   Bench test on remaining tissue (the unused 60-day fresh-frozen cohort includes RB4403, RB4401, RB4405):
   mtPstI transcript (qPCR) and mtDNA damage (copy number / deletion assay).
2. **Xenium slide layout** of the 12 sections, to rule batch in or out for the weak non-OL signals.
3. Follow up the immune-presentation genes on stressed OLs (MHC-I, Cd40) with IHC, and test whether the
   myelin-gene loss (claim 7) precedes the demyelination seen at 4-6 months.
4. Link to histology: are RB4627 / RB4401 / RB4405 the animals with immune infiltrates (IHC, flow), and does
   the thalamus / fiber-tract pattern precede the 4-6 month demyelination in corpus callosum and cerebellum?

## Reproduce

```bash
conda activate karospace311
python scripts/run_sample_classifier.py --pseudobulk results/region_time \
    --outdir results/sample_classifier/pseudobulk_general --merge DA-Oligodendrocytes=Oligodendrocytes
# office Mac (h5ad):
python scripts/run_sample_classifier.py --h5ad <novae4.h5ad> \
    --outdir results/sample_classifier/cells_general --merge DA-Oligodendrocytes=Oligodendrocytes
python scripts/export_cell_table.py --h5ad <novae4.h5ad> --outdir results/cell_table
python scripts/run_ol_gene_drivers.py --h5ad <novae4.h5ad> \
    --cell-scores results/sample_classifier/cells_general/cell_scores/Oligodendrocytes.csv.gz --outdir results/ol_gene_drivers
# figures + tables:
jupyter nbconvert --to notebook --execute --inplace analysis/02-ol-response-and-classifiers.ipynb
jupyter nbconvert --to notebook --execute --inplace analysis/03-ol-stress-genes.ipynb
```
