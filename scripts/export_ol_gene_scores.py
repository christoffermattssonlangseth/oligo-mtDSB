#!/usr/bin/env python
"""Per-cell gene-set scores for every oligodendrocyte (needs the h5ad).

For each OL (general label: DA-Oligodendrocytes folded into Oligodendrocytes) and each gene set, the score
is ``log1p(1e4 * counts in the set / total counts of the cell)``:

* ``isr``          - the 10 prespecified ISR genes (``utils.ol_response.ISR_GENES``)
* ``drivers_up``   - the up driver genes (``results/ol_gene_drivers/tables/driver_genes.csv``)
* ``myelin_down``  - the down driver genes (mature-OL / myelin)

Output: ``--out`` (csv.gz) with cell_id, sample_id, region, total_counts and the three scores.

Example (office Mac)
--------------------
    python scripts/export_ol_gene_scores.py --h5ad <novae4.h5ad> --out results/cell_table/ol_gene_scores.csv.gz
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
from utils.region_time_de import Keys, read_minimal  # noqa: E402


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--h5ad", required=True)
    p.add_argument("--out", default=str(REPO / "results" / "cell_table" / "ol_gene_scores.csv.gz"))
    p.add_argument("--drivers", default=str(REPO / "results" / "ol_gene_drivers" / "tables" / "driver_genes.csv"))
    p.add_argument("--layer", default="counts")
    a = p.parse_args(argv)

    adata = read_minimal(a.h5ad, layer=a.layer, keys=Keys())
    ct = adata.obs["cell_class"].astype(str).map(lambda c: olr.GENERAL_LABELS.get(c, c))
    sel = np.flatnonzero((ct == olr.OL).to_numpy())
    X = sparse.csr_matrix(adata.layers[a.layer] if a.layer else adata.X)[sel]
    genes = pd.Index(adata.var_names)
    drivers = pd.read_csv(a.drivers, index_col=0)
    sets = {"isr": list(olr.ISR_GENES),
            "drivers_up": drivers.index[drivers.direction == "up"].tolist(),
            "myelin_down": drivers.index[drivers.direction == "down"].tolist()}
    total = np.asarray(X.sum(1)).ravel()
    out = pd.DataFrame({"sample_id": adata.obs["sample_id"].astype(str).to_numpy()[sel],
                        "region": adata.obs["RBD_compartment_simplified"].astype(str).to_numpy()[sel],
                        "total_counts": total}, index=adata.obs_names[sel])
    out.index.name = "cell_id"
    for name, gs in sets.items():
        idx = genes.get_indexer([g for g in gs if g in genes])
        s = np.asarray(X[:, idx].sum(1)).ravel()
        out[name] = np.log1p(1e4 * s / np.where(total > 0, total, 1))
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(a.out, compression="gzip")
    print(f"wrote {a.out}: {len(out)} OLs; sets: " + ", ".join(f"{k}={len(v)}" for k, v in sets.items()))


if __name__ == "__main__":
    main()
