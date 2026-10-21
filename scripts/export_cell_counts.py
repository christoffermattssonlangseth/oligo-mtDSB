#!/usr/bin/env python
"""Raw counts of selected genes for selected cell types, one row per cell (needs the h5ad).

For per-cell questions that the gene-set scores cannot answer, e.g. analysis/12 G (proliferation within
lateral-ventricle NSCs), analysis/13 (ISR genes in OPCs) and analysis/14 (anatomy genes in every cell). Output
columns: cell_id, sample_id, cell_type, region, total_counts, then one column per gene. Cell types use the
``cell_class`` labels as in results/cell_table. Without ``--cell-types`` all cells are kept. ``--genes-file`` reads one
gene per line. An ``--out`` ending in ``.npz`` writes a sparse matrix (cells x genes, ``scipy.sparse.save_npz``) plus
``<out>.obs.csv.gz`` and ``<out>.genes.txt`` instead of a dense table.

Example (office Mac)
--------------------
    python scripts/export_cell_counts.py --h5ad <novae4.h5ad> --out results/cell_table/nsc_opc_counts.csv.gz \
        --cell-types "Neural Stem Cells" "Immature Oligodendrocytes I" --genes Foxg1 Tbx3 Mki67 Trib3
    python scripts/export_cell_counts.py --h5ad <novae4.h5ad> --out anatomy_counts.npz --genes-file anatomy_genes.txt
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

from utils.region_time_de import Keys, read_minimal  # noqa: E402


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--h5ad", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--cell-types", nargs="+")
    p.add_argument("--genes", nargs="+", default=[])
    p.add_argument("--genes-file")
    a = p.parse_args(argv)
    if a.genes_file:
        a.genes += [g.strip() for g in open(a.genes_file) if g.strip()]
    adata = read_minimal(a.h5ad, layer="counts", keys=Keys())
    keep = adata.obs["cell_class"].astype(str).isin(a.cell_types).to_numpy() if a.cell_types else np.ones(adata.n_obs, bool)
    X = sparse.csr_matrix(adata.layers["counts"])[keep]
    genes = pd.Index(adata.var_names)
    missing = [g for g in a.genes if g not in genes]
    if missing:
        print("not on the panel:", missing)
    use = [g for g in a.genes if g in genes]
    obs = adata.obs[keep]
    out = pd.DataFrame({"sample_id": obs["sample_id"].astype(str).to_numpy(), "cell_type": obs["cell_class"].astype(str).to_numpy(),
                        "region": obs["RBD_compartment_simplified"].astype(str).to_numpy(),
                        "total_counts": np.asarray(X.sum(1)).ravel()}, index=obs.index)
    out.index.name = "cell_id"
    M = X[:, genes.get_indexer(use)]
    if a.out.endswith(".npz"):
        sparse.save_npz(a.out, sparse.csr_matrix(M, dtype=np.int32))
        out.to_csv(a.out + ".obs.csv.gz", compression="gzip")
        Path(a.out + ".genes.txt").write_text("\n".join(use) + "\n")
    else:
        out[use] = M.toarray().astype(np.int32)
        out.to_csv(a.out, compression="gzip")
    print(f"wrote {a.out}: {len(out):,} cells x {len(use)} genes")


if __name__ == "__main__":
    main()
