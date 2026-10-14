# oligo-mtDSB
<p align="center">
  <img src="assets/logo.png" alt="oligo-mtDSB Logo" width="300"/>
</p>

Consequences of double-stranded breaks (DSBs) in mitochondrial DNA (mtDNA) of oligodendrocytes (OLs), in a
genetically engineered mouse model (PlpTA:mtPst1), measured with 10x Xenium 5k spatial transcriptomics.

**Design.** 12 sections, one per mouse: DOX removal at P21 or P60 x control / mtDSB x 3, all collected 42 days
after induction. Cells carry a cell type (`cell_class`) and an anatomical compartment
(`RBD_compartment_simplified`). The mouse is the biological replicate in every test.

---

## Headline (2026-10-14)

**At 42 days after induction, mtDSB produces a strong oligodendrocyte stress response in 3 of the 6 mtDSB mice
(RB4405, RB4401, RB4627) and no detectable response in the other 3 (RB4403, RB4630, RB4653).**

* **How we found it.** Leave-one-mouse-out classifiers: each mouse's group is predicted by a model that never
  saw it, per cell type and region, using general cell-type labels (the disease-defined DA-OL cluster is
  folded into Oligodendrocytes). Age is predicted almost perfectly; mtDSB vs control is not, because half of
  the mtDSB mice look like controls.
* **Three independent readouts agree on the same three mice:** a predefined ISR gene set, the label-free OL
  classifier, and the DA-OL cluster. On the prespecified ISR score the mtDSB group is ~10x more variable than
  the controls (p = 1e-4) while the group means do not differ: heterogeneous, not uniformly shifted.
* **What responding OLs do** (34 genes up, 9 down; defined *within* tissue sections, so batch, section plane
  and time of sacrifice cannot create them): ATF4 / integrated stress response (Atf4, Atf5, Ddit3, Trib3,
  Sesn2), mitochondrial and protein quality control (Hspa9, Sqstm1/p62, Gstp1, Mt2), p21 / cell cycle,
  immune-visibility genes (MHC-I H2-D1 / B2m, Cd40, Serpina3n), and loss of myelin genes (Mog, Cnp, Cldn11,
  Mal). The top genes replicate between models trained on disjoint mice (P21-only vs P60-only; 10 / 50
  shared, p = 2e-11). Microglia in the same three mice show a DAM-like response (Cd74, Cst7 10-30x).
* **Where.** Throughout the brain in responders, strongest in the thalamus and fiber tracts.
* **Not explained by** sex, genotype label, sample swap (sex labels match Xist) or white-matter content.
  Open: induction efficiency (no transgene probe on the panel).

<p align="center">
  <img src="results/sample_classifier/figures/F2_three_readouts.png" alt="three readouts" width="900"/>
  <br><em>Three readouts, same three mice (filled orange = responders, open = non-responders, blue = controls).</em>
</p>

Full claims with evidence and strength: **[`results/sample_classifier/FINDINGS.md`](results/sample_classifier/FINDINGS.md)**.

### How this confirms and explains the earlier analyses

The new analysis does not overturn the earlier work; it shows that the earlier signal was real but carried by
half of the mtDSB mice, which also explains why it looked weak or age-dependent.

* **Hierarchical Bayesian case-control model** (`legacy/notebooks-02/09-12`, `legacy/results/baysian_diff_condition.csv`,
  OLs, ages pooled): all 22 of today's driver genes that it tested have the **same direction** (17 up, 5 down;
  sign test p < 1e-6), and 12 reached its FDR < 0.05 (Mt2, Gstp1, B2m, H2-D1, Atf4, Atf5, Olig2 up; Cldn11, Cnp,
  Mal, Mog, Ugt8a down). Its posterior shifts were small (0.01-0.05) because every estimate averaged
  responders with non-responders.
* **PyMC models** (`legacy/notebooks-01/29-30`): "myelin identity" down and Serpina3n up (beta = 2.1), as now.
* **Per-sample markers** (`legacy/notebooks-01/06`): RB4405 and RB4627 already showed Atf4, Atf5, Gstp1,
  Serpina3n, H2-D1, B2m, Mt2, Olig2, i.e. today's driver genes.
* **DA-OL sub-clustering** (`legacy/notebooks-02/06`, `legacy/notebooks-04/06`): same ISR / redox / immune programme;
  the cluster is the responder OL state.
* **Region x time DE** (`results/region_time/`, paper sec. 6): same hotspots (thalamus, fiber tracts, hypothalamus).
* **What changes:** the "P21 modest / P60 strong" progression is responder composition (1 of 3 vs 2 of 3 mice
  respond; the P21 responder is as strong as any); top genes such as Sgk1, Fkbp5, Zbtb16 (stress hormones),
  Kif5a (one odd section) and Il33 (one control) were section-level confounds. Claim-by-claim review of every
  earlier notebook: **[`docs/ANALYSIS_REVIEW.md`](docs/ANALYSIS_REVIEW.md)**.

---

## Start here

| I want to ... | Go to |
|---|---|
| read the findings and their evidence | `results/sample_classifier/FINDINGS.md` |
| see the figures (F1-F8 responders / regions / score maps / checks; G1-G4 genes) | `results/sample_classifier/figures/`, `results/ol_gene_drivers/figures/` |
| check whether an earlier notebook's conclusion still holds | `docs/ANALYSIS_REVIEW.md` |
| regenerate figures and tables (laptop, no h5ad needed) | `analysis/02-ol-response-and-classifiers.ipynb`, `analysis/03-ol-stress-genes.ipynb` |
| region x time-point DE figures | `analysis/01-region-time-treatment-figures.ipynb` |
| understand the methods | docstrings of `utils/sample_classifier.py`, `utils/ol_response.py`, `utils/region_time_de.py` |
| check the code works without the data | `python tests/test_sample_classifier.py`, `tests/test_ol_response.py`, `tests/test_region_time_de.py` |
| edit the A0 poster (KI style, editable SVG; regenerate with `make_poster.py`) | `docs/poster/poster.svg` |
| read the manuscript draft (with review notes) | `results/paper/manuscript.md` (built from `results/paper/[0-9]*-*.md`) |
| see how the annotation / compartments were made | `legacy/notebooks-02/01-03`, `legacy/notebooks-02/13`, `legacy/notebooks-04/00`, `scripts/compartments/` |

## Repository layout

```
analysis/        current notebooks (figures + tables from the result folders; no h5ad needed)
  01-region-time-treatment-figures.ipynb   region x time-point DE
  02-ol-response-and-classifiers.ipynb     responders, classifiers, regions, checks   -> F1-F7
  03-ol-stress-genes.ipynb                 driver genes of the OL response            -> G1-G4
utils/           maintained modules: sample_classifier, ol_response, region_time_de, region_time_plots, spatial_utils
scripts/         batch runs (need the h5ad): run_sample_classifier, run_ol_gene_drivers, run_region_time_de,
                 export_cell_table; compartments/ = read-based compartment cluster job
tests/           synthetic tests for the maintained modules
results/
  sample_classifier/   FINDINGS.md, figures F1-F7, tables, classifier outputs (per-cell scores gitignored)
  ol_gene_drivers/     gene weights, within-section associations, driver genes, figures G1-G4
  region_time/         region x time DE tables (+ _cell_class_updated, _MANA250 example)
  cell_table/          sections.csv (age, condition, sex, genotype); cells.csv.gz gitignored
  paper/               manuscript sections + build_manuscript.py
docs/            ANALYSIS_REVIEW.md; poster/ (A0 poster SVG + generator)
legacy/          earlier work, kept for the record (not maintained; relative paths inside no longer resolve)
  notebooks-01 ... notebooks-04, utils_builder.ipynb
  utils/       auto-harvested helper modules
  results/     age, region, bayes, casecontrol_*, figures*, cell_types_*, compartments_all, rbd_highlight_*
assets/  data/ (gitignored)  logs/
```

---

## Reproducing the analysis

### Data location

The annotated h5ad files are **not in the repository**. They live on an external SSD (`/Volumes/processing2`)
normally attached to the office Mac (reachable over Tailscale; here `<office-mac>`). Copying over Tailscale
runs at ~2 MB/s (10 GB takes ~1.5 h), so run the h5ad steps on the office Mac or plug the disk in.

```
/Volumes/processing2/oligo-mtDSB/data/
    mtDNA_DSB_5k_clustered_annotation_with_rbd_2_cytetype_brain_novae4.h5ad   <- current: 980 474 cells x 5 101 genes,
                                                                                 layers['counts'] raw (10.4 GB)
    mtDNA_DSB_5k_clustered_annotation_MANA.h5ad                                <- 250 genes only; NOT for DE
    ... older intermediates
```

### Environments

| conda env | use for |
|---|---|
| `karospace311` (laptop) | all maintained code, tests, `analysis/` notebooks (scanpy, scikit-learn, pydeseq2 0.5.4) |
| `cellcharter-env` (office Mac) | the h5ad steps below (scikit-learn, anndata, pydeseq2) |

### Steps

```bash
H5AD=/Volumes/processing2/oligo-mtDSB/data/mtDNA_DSB_5k_clustered_annotation_with_rbd_2_cytetype_brain_novae4.h5ad
# office Mac (h5ad):
python scripts/run_region_time_de.py --h5ad $H5AD --outdir results/region_time --layer counts --n-cpus 12   # 27 min
python scripts/run_sample_classifier.py --h5ad $H5AD --outdir results/sample_classifier/cells_general \
    --merge DA-Oligodendrocytes=Oligodendrocytes                                                            # ~30 min, resumable
python scripts/run_ol_gene_drivers.py --h5ad $H5AD --outdir results/ol_gene_drivers \
    --cell-scores results/sample_classifier/cells_general/cell_scores/Oligodendrocytes.csv.gz               # 1-2 min
python scripts/export_cell_table.py --h5ad $H5AD --outdir results/cell_table                               # 1 min
# laptop (no h5ad):
python scripts/run_sample_classifier.py --pseudobulk results/region_time \
    --outdir results/sample_classifier/pseudobulk_general --merge DA-Oligodendrocytes=Oligodendrocytes     # 3 min
jupyter nbconvert --to notebook --execute --inplace analysis/0[1-3]-*.ipynb
```

### Methods in brief

* **Classifiers** (`utils/sample_classifier.py`): L2 logistic regression, balanced classes. Pseudobulk level:
  one model per (cell type, region) on centered log-CPM, exact kernel solution, sample-level permutation tests
  (floors: 1/400 pooled, 1/20 within-cohort / transfer). Single-cell level: one model per cell type, 2000
  cells per training section, every cell of the held-out section scored. Tasks: `age`, `condition`, `joint`,
  `condition@21/@60`, `condition 21->60 / 60->21` (train on one cohort, predict the other).
* **Responders and driver genes** (`utils/ol_response.py`): ISR score from 10 prespecified ATF4 targets;
  driver genes = within-section correlation with the out-of-fold OL score, r > 0.1 in each responder and
  < 0.1 in controls and non-responders (fixed before looking at genes).
* **Region x time DE** (`utils/region_time_de.py`): section-level pseudobulk, pyDESeq2 per (cell type,
  region, age), joint `~ age + condition + age:condition` models, region-vs-rest tests.

### Caveats when reading `results/region_time/`

* **Neuropeptide spill-over:** Hcrt, Pmch, Sst, Pdyn appear as top "glial" hits in the hypothalamus; they are
  neuronal transcripts picked up by neighbouring segmented cells.
* **Neural Stem Cells @ Ventricular system @ P60** (671 hits) reflects section plane / ventricle coverage,
  not regulated genes.
* Group-level mtDSB vs control tests average responders with non-responders (see Headline).

---

## Analysis history

1. Cluster annotation, refined with cytetype and manual merging (`legacy/notebooks-02/13`, `legacy/notebooks-04/00`).
2. Read-based anatomical compartments (`scripts/compartments/`, `legacy/notebooks-02/03`).
3. DEGs across condition and time, whole tissue and per region (`legacy/notebooks-01/25-26`, `legacy/notebooks-02/04`).
4. Hierarchical Bayesian case-control models (`legacy/notebooks-01/28-30`, `legacy/notebooks-02/09-12`).
5. 2026-09: consolidated region x time-point DE pipeline (`utils/region_time_de.py`, `analysis/01`).
6. 2026-10-14: leave-one-mouse-out classifiers, OL responder finding, within-section driver genes
   (`analysis/02-03`); all earlier analyses re-checked in `docs/ANALYSIS_REVIEW.md`; repository split into
   current (`analysis/`, `results/`) and `legacy/`.

## Background references
- https://pmc.ncbi.nlm.nih.gov/articles/PMC5647772/
- https://scholarship.miami.edu/view/pdfCoverPage?instCode=01UOML_INST&filePid=13386228880002976&download=true
