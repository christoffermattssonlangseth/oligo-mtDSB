#!/usr/bin/env python
"""Slim per-cell table of the human CNS Xenium 5k run (MS + Baló sections) for analysis/17.

Reads the clustered h5ad (raw counts in ``layers/counts``) and writes one row per cell: sample, Leiden clusters,
centroid, area, total counts and raw counts of ``GENES`` (mtDSB markers, shared stress, interferon, myelin,
lineage and cell-type markers).

Example (office Mac, env cellcharter-env)
-----------------------------------------
    python scripts/export_human_xenium_cells.py \
        --h5ad /Volumes/moldiassd/human_MS_Balo/processed_data/baloMS_clustered.h5ad --out results/human_xenium/cells.csv.gz
"""

from __future__ import annotations

import argparse

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

GENES = ["TRIB3", "GDF15", "EIF4EBP1", "DDIT3", "CDKN1A", "SESN2", "BBC3", "GADD45A", "ATF4", "VEGFA", "HMOX1", "BNIP3",
         "SLC2A1", "HIF1A", "IRF1", "STAT1", "SERPINA3", "CD74", "HLA-B", "B2M",
         "MOG", "MBP", "MAG", "CNP", "MOBP", "OPALIN", "SOX10", "OLIG2", "PDGFRA",
         "GFAP", "AQP4", "P2RY12", "CX3CR1", "CD68", "SNAP25", "CLDN5", "PTPRC", "CD3E"]


def obs_col(f, k):
    g = f["obs"][k]
    if isinstance(g, h5py.Group):
        cats = np.array([c.decode() for c in g["categories"][:]], dtype=object)
        return cats[g["codes"][:]]
    return g[:]


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--h5ad", required=True)
    p.add_argument("--out", required=True)
    a = p.parse_args(argv)
    with h5py.File(a.h5ad, "r") as f:
        v = f["var"]; genes = pd.Index([x.decode() for x in v[v.attrs["_index"]][:]])
        out = pd.DataFrame({k: obs_col(f, k) for k in ["sample_id", "cell_id", "leiden_1", "leiden_2", "x_centroid", "y_centroid", "cell_area"]})
        X = f["layers/counts"]
        X = sparse.csr_matrix((X["data"][:], X["indices"][:], X["indptr"][:]), shape=(len(out), len(genes)))
    out["total_counts"] = np.asarray(X.sum(1)).ravel()
    use = [g for g in GENES if g in genes]
    out[use] = X[:, genes.get_indexer(use)].toarray().astype(np.int32)
    out.to_csv(a.out, index=False, compression="gzip")
    print(f"wrote {a.out}: {len(out):,} cells; not on panel: {[g for g in GENES if g not in genes]}")


if __name__ == "__main__":
    main()
