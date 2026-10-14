#!/usr/bin/env python
"""Predict age cohort / condition / both per section (leave-one-section-out).

Two levels, see ``utils/sample_classifier.py``:

* ``--pseudobulk DIR``: reuse the pseudobulk written by
  ``run_region_time_de.py`` (``DIR/pseudobulk_*.csv``). Fast (minutes), runs
  anywhere, includes permutation tests.
* ``--h5ad FILE``: single-cell models per cell type on the annotated h5ad
  (needs the full data, i.e. the office Mac or the SSD). Writes one
  out-of-fold score per cell to OUTDIR/cell_scores/<cell type>.csv.gz and
  refreshes predictions.csv / metrics.csv / progress.json after every cell
  type; rerunning the same command resumes where it stopped.

Examples
--------
    conda activate karospace311
    python scripts/run_sample_classifier.py --pseudobulk results/region_time \
        --outdir results/sample_classifier/pseudobulk

    python scripts/run_sample_classifier.py \
        --h5ad /Volumes/processing2/oligo-mtDSB/data/mtDNA_DSB_5k_clustered_annotation_with_rbd_2_cytetype_brain_novae4.h5ad \
        --outdir results/sample_classifier/cells_general --merge "DA-Oligodendrocytes=Oligodendrocytes"

Use ``--merge FROM=TO`` to fold disease-defined clusters into their general
class (here DA-Oligodendrocytes -> Oligodendrocytes), so the classifier, not
the clustering, has to find the treatment effect.
"""

from __future__ import annotations

import argparse
import logging
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from utils import sample_classifier as sc  # noqa: E402
from utils.region_time_de import Keys, Pseudobulk, merge_cell_types, read_minimal  # noqa: E402


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--pseudobulk", help="directory with pseudobulk_counts.csv / pseudobulk_meta.csv")
    src.add_argument("--h5ad", help="annotated AnnData with raw counts (single-cell models)")
    p.add_argument("--outdir", required=True)
    p.add_argument("--tasks", nargs="+", default=None, help="subset of task names (default: all, see default_tasks)")
    p.add_argument("--C", type=float, default=None, help="inverse L2 strength (default 10 pseudobulk, 0.1 cells)")
    # pseudobulk
    p.add_argument("--n-perm", type=int, default=1000)
    p.add_argument("--min-cpm", type=float, default=10.0)
    p.add_argument("--merge", nargs="+", default=[], metavar="FROM=TO",
                   help='relabel cell types before classifying, e.g. "DA-Oligodendrocytes=Oligodendrocytes" (both levels)')
    # cells
    p.add_argument("--layer", default="counts")
    p.add_argument("--cell-type-key", default="cell_class")
    p.add_argument("--region-key", default="RBD_compartment_simplified")
    p.add_argument("--cell-types", nargs="+", default=None)
    p.add_argument("--cells-per-sample", type=int, default=2000)
    p.add_argument("--min-cells-per-sample", type=int, default=50)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--no-resume", action="store_true", help="refit cell types whose scores are already in OUTDIR/cell_scores/")
    return p.parse_args(argv)


def main(argv=None):
    a = parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    keys = Keys(cell_type=a.cell_type_key, region=a.region_key)
    tasks = sc.default_tasks(keys)
    if a.tasks:
        tasks = [t for t in tasks if t.name in set(a.tasks)]
    mapping = dict(m.split("=", 1) for m in a.merge)
    t0 = time.time()
    info = {"tasks": sc.tasks_info(tasks), "keys": vars(keys) if hasattr(keys, "__dict__") else str(keys)}

    if a.pseudobulk:
        pb = Pseudobulk.load(a.pseudobulk)
        if mapping:
            pb = merge_cell_types(pb, mapping)
        C = 10.0 if a.C is None else a.C
        tables = sc.classify_pseudobulk(pb, keys, tasks, C=C, min_cpm=a.min_cpm, n_perm=a.n_perm, seed=a.seed)
        info.update(level="pseudobulk", source=a.pseudobulk, C=C, min_cpm=a.min_cpm, n_perm=a.n_perm,
                    cell_type_map=mapping)
    else:
        adata = read_minimal(a.h5ad, layer=a.layer, keys=keys)
        logging.info("loaded %d cells x %d genes (%.0f s)", adata.n_obs, adata.n_vars, time.time() - t0)
        C = 0.1 if a.C is None else a.C
        tables = sc.classify_cells(adata, keys, layer=a.layer, tasks=tasks, cell_types=a.cell_types,
                                   cells_per_sample=a.cells_per_sample, min_cells_per_sample=a.min_cells_per_sample,
                                   C=C, seed=a.seed, outdir=a.outdir, resume=not a.no_resume, cell_type_map=mapping)
        info.update(level="cells", source=a.h5ad, cell_type_map=mapping, C=C, cells_per_sample=a.cells_per_sample,
                    min_cells_per_sample=a.min_cells_per_sample, n_cells=int(adata.n_obs), n_genes=int(adata.n_vars))

    info["seconds"] = round(time.time() - t0, 1)
    sc.save_results(a.outdir, tables, info)
    logging.info("wrote %s (%.0f s)", a.outdir, time.time() - t0)


if __name__ == "__main__":
    main()
