#!/usr/bin/env python
"""Per-cell scores for every prespecified gene set (utils/gene_sets.py), all cells (needs the h5ad).

Score = ``log1p(1e4 * counts in the set / total counts of the cell)``. Cell types use general labels
(DA-Oligodendrocytes folded into Oligodendrocytes). Output columns: cell_id, sample_id, cell_type, region,
total_counts, one column per gene set.

Example (office Mac)
--------------------
    python scripts/export_cell_gene_scores.py --h5ad <novae4.h5ad> --out results/cell_table/cell_gene_scores.csv.gz
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from utils import ol_response as olr  # noqa: E402
from utils.gene_sets import GENE_SETS  # noqa: E402
from utils.region_time_de import Keys, read_minimal  # noqa: E402


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--h5ad", required=True)
    p.add_argument("--out", default=str(REPO / "results" / "cell_table" / "cell_gene_scores.csv.gz"))
    p.add_argument("--layer", default="counts")
    a = p.parse_args(argv)
    adata = read_minimal(a.h5ad, layer=a.layer, keys=Keys())
    X = sparse.csr_matrix(adata.layers[a.layer] if a.layer else adata.X)
    genes = pd.Index(adata.var_names)
    total = np.asarray(X.sum(1)).ravel()
    out = pd.DataFrame({
        "sample_id": adata.obs["sample_id"].astype(str).to_numpy(),
        "cell_type": adata.obs["cell_class"].astype(str).map(lambda c: olr.GENERAL_LABELS.get(c, c)).to_numpy(),
        "region": adata.obs["RBD_compartment_simplified"].astype(str).to_numpy(),
        "total_counts": total}, index=adata.obs_names)
    out.index.name = "cell_id"
    denom = np.where(total > 0, total, 1)
    for name, gs in GENE_SETS.items():
        idx = genes.get_indexer([g for g in gs if g in genes])
        out[name] = np.round(np.log1p(1e4 * np.asarray(X[:, idx].sum(1)).ravel() / denom), 4)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(a.out, compression="gzip")
    print(f"wrote {a.out}: {len(out)} cells x {len(GENE_SETS)} gene sets")


if __name__ == "__main__":
    main()
