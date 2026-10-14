#!/usr/bin/env python
"""Which genes drive the oligodendrocyte mtDSB classifier? (needs the h5ad)

Two outputs in ``--outdir``:

1. ``coefficients/Oligodendrocytes.csv.gz`` - gene weights of every fold model
   of the OL single-cell classifier (general labels) for three tasks:
   ``condition`` (12 leave-one-section-out folds), ``condition 21->60``
   (trained on the 6 P21 animals only) and ``condition 60->21`` (trained on the
   6 P60 animals only). The two cohort models share no animals, so agreement
   between their weights is a replication.
2. ``association.csv.gz`` - within-section Pearson correlation between every
   gene's log-normalized expression and the out-of-fold OL score
   ``p_pos[condition]`` from the main run (``--cell-scores``), per section.

Example (office Mac)
--------------------
    python scripts/run_ol_gene_drivers.py \
        --h5ad /Volumes/processing2/oligo-mtDSB/data/mtDNA_DSB_5k_clustered_annotation_with_rbd_2_cytetype_brain_novae4.h5ad \
        --cell-scores results/sample_classifier/cells_general/cell_scores/Oligodendrocytes.csv.gz \
        --outdir results/ol_gene_drivers
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from utils import ol_response as olr  # noqa: E402
from utils import sample_classifier as sc  # noqa: E402
from utils.region_time_de import Keys, read_minimal  # noqa: E402

TASKS = ("condition", "condition 21->60", "condition 60->21")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--h5ad", required=True)
    p.add_argument("--cell-scores", required=True, help="Oligodendrocytes.csv.gz from the main cells_general run")
    p.add_argument("--outdir", required=True)
    p.add_argument("--layer", default="counts")
    p.add_argument("--cells-per-sample", type=int, default=2000)
    p.add_argument("--C", type=float, default=0.1)
    a = p.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    t0 = time.time()
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    keys = Keys()
    adata = read_minimal(a.h5ad, layer=a.layer, keys=keys)
    logging.info("loaded %d cells (%.0f s)", adata.n_obs, time.time() - t0)

    # 1. fold-model weights
    tasks = [t for t in sc.default_tasks(keys) if t.name in TASKS]
    res = sc.classify_cells(adata, keys, layer=a.layer, tasks=tasks, cell_types=[olr.OL],
                            cells_per_sample=a.cells_per_sample, C=a.C, outdir=out / "model",
                            resume=False, cell_type_map=olr.GENERAL_LABELS, save_coefficients=True)
    res["coefficients"].to_csv(out / "coefficients.csv.gz", index=False, compression="gzip")
    logging.info("coefficients written (%.0f s)", time.time() - t0)

    # 2. within-section association with the published out-of-fold scores
    scores = pd.read_csv(a.cell_scores, index_col=0)
    col = "p_pos[condition]"
    scores = scores.loc[scores[col].notna()]
    pos = adata.obs_names.get_indexer(scores.index)
    if (pos < 0).any():
        raise KeyError(f"{(pos < 0).sum()} scored cells not in the h5ad")
    X = sparse.csr_matrix(adata.layers[a.layer] if a.layer else adata.X)
    Xn, target = sc._lognorm_cells(X)  # same normalization as the classifier (median total over all cells)
    assoc = olr.gene_score_association(Xn[pos], scores["sample_id"].to_numpy(), scores[col].to_numpy(),
                                       adata.var_names)
    assoc.to_csv(out / "association.csv.gz", index=False, compression="gzip")
    info = {"h5ad": a.h5ad, "cell_scores": a.cell_scores, "tasks": list(TASKS), "C": a.C,
            "cells_per_sample": a.cells_per_sample, "norm_target": target, "n_scored_cells": int(len(scores)),
            "seconds": round(time.time() - t0, 1)}
    (out / "run_info.json").write_text(json.dumps(info, indent=2))
    logging.info("wrote %s (%.0f s)", out, time.time() - t0)


if __name__ == "__main__":
    main()
