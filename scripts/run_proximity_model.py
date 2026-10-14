#!/usr/bin/env python
"""Can a neighbour cell's own expression tell how close it is to a stressed oligodendrocyte? (needs the h5ad)

Target per cell (``results/proximity_model/targets.csv.gz``, from ``utils.neighbourhood.proximity_targets``):
log distance to the nearest stressed OL, residualised within each section on region, distance to any OL and
local OL density. Features: log-normalised expression of the cell, standardised on the training cells.
Model: ridge regression (alpha chosen by efficient leave-one-out on the training cells).

Feature sets (to rule out transcript spill-over from the neighbouring OL itself):
    all       every gene
    strict    minus OL-enriched genes (> 2x CPM in OLs vs this cell type in controls), minus the OL driver genes
              and the ISR genes

Evaluation, per cell type and feature set (Pearson r between predicted and true target in the held-out mouse):
    resp_loso   train on 2 responders, test on the 3rd
    ctrl_loso   train on 5 controls, test on the 6th (the null)
    resp_to_ctrl / resp_to_nonresp   train on all responders, test on each control / non-responder

Writes ``metrics.csv`` and ``coefficients.csv.gz`` (weights of the model fitted on all responders) to --outdir.
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

from utils import sample_classifier as sc  # noqa: E402
from utils.ol_response import ISR_GENES  # noqa: E402
from utils.region_time_de import Keys, read_minimal  # noqa: E402


def fit_predict(Xtr, ytr, Xte_list):
    from sklearn.linear_model import RidgeCV
    mu, sd = Xtr.mean(0), Xtr.std(0)
    sd[sd == 0] = 1
    m = RidgeCV(alphas=np.logspace(1, 5, 9)).fit((Xtr - mu) / sd, ytr)
    return m, [m.predict((X - mu) / sd) for X in Xte_list], (mu, sd)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--h5ad", required=True)
    p.add_argument("--targets", default=str(REPO / "results/proximity_model/targets.csv.gz"))
    p.add_argument("--ol-enriched", default=str(REPO / "results/proximity_model/ol_enriched_genes.json"))
    p.add_argument("--drivers", default=str(REPO / "results/ol_gene_drivers/tables/driver_genes.csv"))
    p.add_argument("--outdir", default=str(REPO / "results/proximity_model"))
    p.add_argument("--cells-per-section", type=int, default=3000)
    p.add_argument("--seed", type=int, default=0)
    a = p.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    t0 = time.time()
    rng = np.random.default_rng(a.seed)
    tg = pd.read_csv(a.targets, index_col=0)
    adata = read_minimal(a.h5ad, layer="counts", keys=Keys())
    Xn, _ = sc._lognorm_cells(sparse.csr_matrix(adata.layers["counts"]))
    genes = np.asarray(adata.var_names)
    pos = pd.Series(np.arange(adata.n_obs), index=adata.obs_names)
    ol_enriched = json.loads(Path(a.ol_enriched).read_text())
    drivers = set(pd.read_csv(a.drivers, index_col=0).index) | set(ISR_GENES)
    metrics, coefs = [], []
    for ct, d in tg.groupby("cell_type"):
        # subsample per section once, reuse across schemes
        idx = []
        for sid, g in d.groupby("sample_id"):
            take = g.index if len(g) <= a.cells_per_section else rng.choice(g.index, a.cells_per_section, replace=False)
            idx.extend(take)
        d = d.loc[idx]
        Xall = Xn[pos.loc[d.index].to_numpy()].toarray().astype(np.float32)
        excl = set(ol_enriched.get(ct, [])) | drivers
        fsets = {"all": np.ones(len(genes), bool), "strict": ~np.isin(genes, list(excl))}
        sections = d.groupby("sample_id")["group"].first()
        R = list(sections.index[sections == "responder"])
        C = list(sections.index[sections == "control"])
        N = list(sections.index[sections == "non-responder"])
        for fs, mask in fsets.items():
            X = Xall[:, mask]
            y = d["y_resid"].to_numpy()
            sid = d["sample_id"].to_numpy()

            def run(train, tests, scheme):
                tr = np.isin(sid, train)
                m, preds, _ = fit_predict(X[tr], y[tr], [X[sid == t] for t in tests])
                for t, pr in zip(tests, preds):
                    yt = y[sid == t]
                    metrics.append({"cell_type": ct, "features": fs, "scheme": scheme, "test_section": t,
                                    "test_group": sections[t], "train": "+".join(train), "n_test": len(yt),
                                    "r": float(np.corrcoef(pr, yt)[0, 1]), "alpha": float(m.alpha_)})
                return m

            for t in R:
                run([s for s in R if s != t], [t], "resp_loso")
            for t in C:
                run([s for s in C if s != t], [t], "ctrl_loso")
            m = run(R, C, "resp_to_ctrl")
            run(R, N, "resp_to_nonresp")
            coefs.append(pd.DataFrame({"cell_type": ct, "features": fs, "gene": genes[mask], "coef": m.coef_}))
            logging.info("%s / %s done (%.0f s)", ct, fs, time.time() - t0)
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(metrics).to_csv(out / "metrics.csv", index=False)
    pd.concat(coefs).to_csv(out / "coefficients.csv.gz", index=False, compression="gzip")
    logging.info("wrote %s (%.0f s)", out, time.time() - t0)


if __name__ == "__main__":
    main()
