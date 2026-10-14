# Review of earlier analyses in light of the OL responder finding (2026-10-14)

**Why this document.** The leave-one-section-out classifiers and the label-free gene analysis
(`results/sample_classifier/FINDINGS.md`) showed that the oligodendrocyte (OL) response to mtDSB is carried by
**3 of the 6 mtDSB mice** (responders RB4405, RB4401, RB4627; non-responders RB4403, RB4630, RB4653 look like
controls). Every earlier analysis compared all six mtDSB mice with all six controls. This document goes
through the earlier notebooks, result folders and manuscript sections and records, claim by claim, what
holds, what does not, and why. Checks were made per animal on the OL (DA-OL folded in), microglia, astrocyte
and NSC pseudobulk of the current `novae4` annotation (`results/region_time/pseudobulk_*.csv`).

## Summary

The biology seen before is largely real: the earliest per-sample marker lists (`notebooks-01/06`) already
show today's driver genes in RB4405 and RB4627. What the earlier analyses could not see is that the effect
comes from half of the mtDSB animals. That single fact explains the weak group-level statistics, the
apparent "P21 modest / P60 strong" age dependence, and why section-level confounds (stress hormones,
circadian and immediate-early genes, one odd section) ended up at the top of several gene lists.

| | count |
|---|---|
| earlier claims confirmed (now attributed to responders) | 9 |
| claims explained by responder composition or confounds | 11 |
| notebooks without testable conclusions (processing, plotting, exploration) | ~35 |

---

## Holds up: seen before, confirmed per animal

| Earlier observation | Where | Per-animal check (current data) |
|---|---|---|
| Per-sample OL markers: RB4405 Gfap, Gstp1, Serpina3n, Atf4, Atf5, B2m, H2-D1, Mt2; RB4627 Gstp1, Atf5, Ccnd1, Olig2, Adgrg1, H2-D1 | `notebooks-01/06`, `01/05` | These are today's driver genes (`results/ol_gene_drivers/tables/driver_genes.csv`). |
| OL "functional modules": oxidative stress (Mt2, Gstp1, Sqstm1, Nfe2l1, Hspa9), ISR (Atf4), MHC-I / Cd40, adhesion (Ptpra, Mpzl1) | `notebooks-01/21` | ~80% overlap with the 34 within-section driver genes. |
| DA-OL cluster with ISR, redox and immune genes; differentiation vs stress-adaptation split | `notebooks-02/06`, `04/06`; paper sec. 2 | Confirmed without the cluster label; the cluster is the responder OL state (22-65% of OLs in responders, 0.3-2.3% elsewhere). |
| Myelin / mature-OL genes down | `notebooks-01/07` (myelin score), `01/30` (Bayes "myelin identity"), `02/07` | Confirmed (Mog, Cnp, Cldn11, Mal, Ugt8a down). The `01/07` myelin score was already lowest in RB4405, RB4627, RB4401. |
| Large ISR fold changes in OLs: Trib3, Gdf15, Cdkn1a, Cd40, Atf5, Eif4ebp1 | `notebooks-02/08`, `01/25` | Real; responder-only (e.g. Trib3 77-214 vs 1-7 CPM). |
| Late OL genes Gpnmb, Cst7, Igfbp3, Cd68, Hmox1 | paper sec. 3 | Real; responder-only (Gpnmb 31-52 vs 2-7 CPM). |
| Microglial DAM-like response: Cd74, Cst7, Cybb, H2-Ab1, Lpl | paper sec. 3, 6 | Real and strong, responder-only: Cd74 876-2705 vs 40-228 CPM, Cst7 394-1100 vs 12-89. The immune partner of the OL response. |
| Regional hotspots: thalamus, fiber tracts, hypothalamus | `notebooks-02/04`; paper sec. 6 | Same regions; thalamus and fiber tracts rank top, hypothalamus mid. |
| Serpina3n strongly up in OLs | `notebooks-01/30` (β = 2.1) | Real, responder-only (242-620 vs 34-59 CPM). |

## Does not hold, or needs rewording

| Earlier claim | Where | What explains it |
|---|---|---|
| **"P21 modest / transient, P60 broad / chronic" OL response** | `notebooks-01/24`, `02/08`, `02/12`; paper abstract, sec. 3, 4 | Responder composition: P21 has 1 responder of 3, P60 has 2 of 3. The P21 responder RB4627 is as strong as any (highest Trib3 and Gdf15). `02/08` shows large P21 fold changes with very wide SEs (Trib3 log2FC 5.1 ± 1.4): one responder among three animals, not a weak response. |
| Hspa5, Sod2 as part of the P60 ISR/UPR programme | `notebooks-02/08`; paper sec. 4 | Not responder-specific; Sod2 is higher in all P60 animals including controls (age effect). |
| Weak / null group-level DE (padj 0.4-1) and tiny Bayesian shifts (0.02-0.06) | `notebooks-01/23`-`25`, `02/09`, `02/12` | Dilution: half of the "cases" are unaffected. |
| Top "DSB vs control" genes Sgk1, Tsc22d3, Fkbp5, Zbtb16, Dnajb1, Hsph1, Ddit4 | `notebooks-01/06`, `01/05`, `01/30` | Section-level confounds: glucocorticoid / heat-shock genes vary between sections regardless of group (Sgk1 highest in control RB4282). `01/05` also shows immediate-early genes (Arc, Fos, Egr1) and circadian genes (Dbp, Nr1d1) in individual sections, i.e. differences in stress / time of day at sacrifice. No within-section signal. |
| Kif5a, Kif5b, Map2, Camk2a as mtDSB genes | `notebooks-01/06`, `01/21` | RB4403 neuronal / segmentation profile; RB4403 is a non-responder. |
| Control markers Il33, Sez6l2, Xist, Kdr | `notebooks-01/06`, `01/11` | Il33 comes from one control (RB4620, mtPst1-only genotype, ~20x others); Xist reflects sex imbalance. |
| Olfactory astrocytes induce Gdf15, Trib3, Cd40, Serpina3n | paper sec. 3 (old 5-glia annotation, `results/figures/*_21w_vs_60w.png`) | Near zero in olfactory astrocytes in every section with the current annotation. |
| NSC proliferative / neurogenic response at P60 | paper sec. 7 | Does not follow responders: Mki67, Dcx, Dlx2 high in two P21 controls and in non-responders RB4630, RB4403; low in all P60 controls and RB4653. Small cell numbers (141-701); most likely section plane (how much ventricular zone is captured). |
| cGAS-STING / interferon (Isg15, Ifit, Tmem173) and NF-kB (Nfkb2, Tnfaip3, Ccl2) activation | `notebooks-02/07` | No responder pattern; at most a hint in Ifit1. |
| Glycolytic shift (Ldha, Pkm, Aldoa, Pgk1, Slc16a1) | `notebooks-02/07`, `01/11` | Only Hk2 rises in responders (25-32 vs 10-17 CPM); the rest only in RB4405 or not at all. |
| P60 inflammatory / dedifferentiation genes Irf8, Il1a, Hoxb3, Fcgr2b | `notebooks-01/24` | Irf8 / Fcgr2b only in RB4405 (likely microglial contamination); Il1a, Hoxb3 flat. |

## Notebooks without testable conclusions

Processing, annotation, mapping, compartment calling, integration and plotting notebooks
(`notebooks-01/00-04, 08-10, 12-17, 19-20, 27, 31-34`; `notebooks-02/01-03, 05, 10-11, 13, 17, 32-34`;
`notebooks-03/*`; `notebooks-04/00, 04`) were checked for stated conclusions; they contain none that bear on
the condition effect. `notebooks-01/18` (cell-type proportion log2FC) shows only tiny changes, consistent
with the composition classifier finding no signal. `notebooks-02/05` (microglia sub-clustering) has no
condition-level conclusions.

## Why the complementary analysis is more robust

1. **Unit of response = the animal.** All earlier approaches averaged within groups (DE, Bayes, pseudotime)
   or within cells (clustering). Predicting held-out animals exposed the heterogeneity instead of averaging it.
2. **Confounds removed by design.** Driver genes must correlate with the OL score *within* a section, so
   section-level differences (stress hormones, circadian genes, plane, segmentation) cannot create them.
3. **Not dependent on clustering or group labels.** A literature ISR gene set, a classifier that never saw
   the DA label, and the cluster agree; models trained on disjoint animals (P21-only vs P60-only) share their
   top genes (10 / 50, p = 2e-11).
4. **Limits.** Three responder animals; the responder grouping comes from the data (the variance test,
   p = 1e-4, is the grouping-free evidence); gene-programme labels are manual. Whether non-response reflects
   induction failure is open (bench test proposed in FINDINGS.md).

## What this means for the manuscript

| Section | Action |
|---|---|
| Abstract, sec. 3, sec. 4 | Replace the P21-vs-P60 progression narrative with the responder finding; keep the gene programmes. |
| Sec. 2 | "DA-OL emerged exclusively in mtDSB tissue" -> "in 3 of 6 mtDSB animals". |
| Sec. 3 (astrocytes) | Drop or re-derive the olfactory-astrocyte panel with the current annotation. |
| Sec. 6 | Keep the regions; re-derive from `results/sample_classifier/` and `results/ol_gene_drivers/`. |
| Sec. 7 (NSC) | Do not use as an mtDSB effect without a section-plane control. |
| New | Responder / non-responder result, within-section driver genes, microglial DAM restricted to responders. |
