# Conflicts found while building the analysis map (2026-10-23)

Agents writing the run catalogue (`docs/run-catalogue.md`) and figure legends checked code, outputs, notebooks and
README against each other. Nothing below has been fixed yet. §N = catalogue entry.

## Worth acting on
- **Human case identifiers in public figures.** `results/human_xenium/figures/X1-X5` (PNG and PDF) print the MS section
  IDs S07-256 / S13-019; the repo is public and README embeds X5. Left out of the map; the catalogue says MS-1 / MS-2.
- **§16 replication p = 2e-11** counts Zbtb16 and Plin4 (glucocorticoid genes that analysis/03 excludes as confounds)
  among the 10 shared top-50 up genes; Il33 (README confound) is among the shared down genes.
- **§17 neighbour coupling** gene sets contain OL driver genes (MHC-I: H2-D1, B2m; reactive astrocyte: Serpina3n, H2-D1)
  with no spill-over control; analysis/12 addresses this only at pseudobulk level.
- **§4 legacy PyMC** (`legacy/notebooks-01/30`): the reported `beta[1]` is the mtDSB group level, not mtDSB - control.
  "Serpina3n beta = 2.1" (ANALYSIS_REVIEW) and "myelin identity down" are not condition effects (findings themselves are
  confirmed elsewhere).
- **docs/ANALYSIS_REVIEW.md** cites Trib3 log2FC 5.1 +/- 1.4 as P21; it is the P60 estimate (P21 3.89 +/- 1.78).
- **§32 EAE overlap:** the highest EAE mouse (z 4.10) is above the mildest responder RB4401 (3.99); notebook text says
  "EAE -4..+4 vs responders +4..+9". CFA controls with lesion-niche labels were pooled into the reference.
- **§31 EAE:** Gdf15 and Eif4ebp1 fail the 20 CPM filter, so the "prespecified 10-gene ISR" bar uses 6 genes; the
  quoted Gdf15 +4.5 vs 0.0 is right (recomputed) but not in any saved output.
- **§24 neurons:** generic Neurons axon-injury coupling p = 0.036 (1 of 8 uncorrected tests) is not mentioned.
- **§47 multiple testing:** pooling all 151 tests gives density hits q 0.032-0.066; within-family best q 0.25. The
  conclusion depends on family definition; say so.

## Number / wording mismatches
- FINDINGS.md: non-responders "inside the control range on every gene" - RB4403 exceeds control max for Atf5, Ddit3, Cth.
  Seven (not six) ISR genes have non-overlapping ranges (Gdf15). "Highest in thalamus in all three responders" holds
  only for the mean.
- §35 proliferating microglia: README/notebook p = 0.012, recomputed and printed in figure C5 p = 0.024.
- §8 OL P60 whole-section hits: 8 genes by the run's own rule (Pmch, Hcrt down; Gabra2, Trh, Gfap, Gpnmb, Rec8, Atf5 up);
  Mt2 fails the |log2FC| 0.5 cut-off.
- §25 README "RB4630 P 0.01-0.34 across priors" omits 0.67 under the "shift up to 6 SD" prior (dismissed in notebook).
- §53 analysis/17 text vs `log2_by_zone.csv`: MOG -2.26/-2.57, OPALIN -2.22/-4.16; OLIG2 falls -1.37/-2.45 in lesion
  OLs (text says lineage kept, on SOX10 alone); neighbour TRIB3 +0.40 in MS-2. Human Xenium OL clusters carry high GFAP.
- README finding 19: ATF5 in Absinta is -0.54 to -0.80, outside "+/-0.3".
- §42 Olig2 per OL +0.67 in fibre tracts (p 0.048); Sox10 per mm2 falls in cortex and fibre tracts (p 0.024).
- §44 pallidum microgliosis p = 0.095 on the new map (notebook summary lists pallidum among strongest); control range 8-13 %.
- §34 figure C2 stars Immature OL I in non-responder hypothalamus (+0.30), contradicting its panel title.
- §29 responder regional ISR range is 2.3-11.5x (RB4401 hippocampus), not 2.5-11.5x.
- §43 relabelled share per region 1.2-6.0 % (olfactory 9.8 %), notebook says 2-6 %.
- §28 control RB4658 tissue T-cell density 2.25 vs "<= 2.2 in any control".
- analysis/05 shrinks control SD (17.1-SD shift optimistic vs 11-19 in analysis/07).
- Figure names G1-G3 exist in both `results/region_map` and `results/ol_gene_drivers`; README says analysis/02 makes F1-F7
  (it makes F1-F9).
- README credits cell labels to CyteType / 04/00; coarse labels come from gpt-4o-mini (02/02), OL sub-labels from 04/06.
  No OPC class exists ("Immature OL I" carries OPC markers) though analysis/01 lists OPC.
- analysis/01 header names novae3; README/commit say novae4. `results/region_time_disease_free/run_info.json` does not
  record the `--region-labels` override.
- analysis/15 labels "13 OL pathway scores"; there are 11.
