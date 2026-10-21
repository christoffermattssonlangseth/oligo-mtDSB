#!/usr/bin/env python
"""Pseudobulk + per-cell gene-set scores of glia in the RR-map EAE spinal-cord Xenium atlas (same 5k panel).

For analysis/10 (mtDSB vs EAE). Reads the annotated RR-map h5ad in row chunks (1.4M cells; never loads the whole
matrix) and keeps oligodendrocytes, myeloid cells and astrocytes (``Anno_L1_curated``). Writes to ``--outdir``:

* ``eae_pseudobulk_counts.csv.gz``  summed raw counts per (mouse, L1, L2, curated niche state)
* ``eae_pseudobulk_meta.csv``       those groups: mouse (``sample_name``), condition, model, stage, n_cells
* ``eae_cell_scores.csv.gz``        per cell: metadata + ``log1p(1e4 * set counts / total)`` for every gene set in
                                    ``utils/gene_sets.py`` and the mtDSB driver genes (up / down)

Example (laptop, ~10 min)
-------------------------
    python scripts/export_eae_pseudobulk.py --h5ad ~/Downloads/RRMAP2_xenium_all_samples...Anno2to4Updated.h5ad
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from utils.gene_sets import GENE_SETS  # noqa: E402

KEEP = ["Oligodendrocyte", "Myeloid", "Astrocyte"]
OBS = ["sample_name", "sample_id", "condition", "model", "stage", "day_of_sacrifice", "score_sacrifice", "region",
       "Anno_L1_curated", "Anno_L2", "Curated_niche_state", "Global_anatomical_region"]
GROUP = ["sample_name", "Anno_L1_curated", "Anno_L2", "Curated_niche_state"]


def read_obs(f) -> pd.DataFrame:
    def col(k):
        g = f["obs"][k]
        if isinstance(g, h5py.Group):
            cats = np.array([c.decode() if isinstance(c, bytes) else c for c in g["categories"][:]], dtype=object)
            codes = g["codes"][:]
            return np.where(codes >= 0, cats[np.clip(codes, 0, None)], None)
        return g[:]
    return pd.DataFrame({k: col(k) for k in OBS})


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--h5ad", required=True)
    p.add_argument("--outdir", default=str(REPO / "results" / "eae_comparison"))
    p.add_argument("--chunk", type=int, default=100_000)
    a = p.parse_args(argv)
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)

    drivers = pd.read_csv(REPO / "results" / "ol_gene_drivers" / "tables" / "driver_genes.csv", index_col=0)
    sets = dict(GENE_SETS, drivers_up=drivers.index[drivers.direction == "up"].tolist(),
                myelin_down=drivers.index[drivers.direction == "down"].tolist())

    with h5py.File(a.h5ad, "r") as f:
        genes = pd.Index([g.decode() if isinstance(g, bytes) else g for g in f["var"][f["var"].attrs["_index"]][:]])
        obs = read_obs(f)
        keep = obs.Anno_L1_curated.isin(KEEP).to_numpy()
        obs = obs.assign(_row=np.arange(len(obs)))
        grp_key = obs[GROUP].astype(str).agg("|".join, axis=1)
        groups = pd.Index(sorted(grp_key[keep].unique()))
        gcode = groups.get_indexer(grp_key)
        S = np.zeros((len(sets), len(genes)), np.float32)   # gene-set membership
        for i, gs in enumerate(sets.values()):
            S[i, genes.get_indexer([g for g in gs if g in genes])] = 1
        pb = np.zeros((len(groups), len(genes)), np.int64)
        scores, totals = [], []
        g = f["layers/counts"]
        indptr = g["indptr"][:]
        n = len(obs)
        for start in range(0, n, a.chunk):
            stop = min(start + a.chunk, n)
            k = np.flatnonzero(keep[start:stop])
            if len(k) == 0:
                continue
            lo, hi = indptr[start], indptr[stop]
            X = sparse.csr_matrix((g["data"][lo:hi], g["indices"][lo:hi], indptr[start:stop + 1] - lo),
                                  shape=(stop - start, len(genes)))[k]
            codes = gcode[start + k]
            M = sparse.csr_matrix((np.ones(len(k)), (codes, np.arange(len(k)))), shape=(len(groups), len(k)))
            pb += (M @ X).toarray().astype(np.int64)
            tot = np.asarray(X.sum(1)).ravel()
            scores.append(np.asarray(X @ S.T))
            totals.append(tot)
            print(f"rows {stop:,}/{n:,}", flush=True)

    meta = pd.DataFrame([k.split("|") for k in groups], columns=GROUP, index=pd.Index(range(len(groups)), name="pb_id"))
    first = obs[keep].groupby("sample_name")[["condition", "model", "stage", "day_of_sacrifice", "score_sacrifice"]].first()
    meta = meta.join(first, on="sample_name")
    meta["n_cells"] = np.bincount(gcode[keep], minlength=len(groups))
    meta.to_csv(out / "eae_pseudobulk_meta.csv")
    pd.DataFrame(pb, columns=genes, index=meta.index).to_csv(out / "eae_pseudobulk_counts.csv.gz", compression="gzip")

    sc = np.vstack(scores)
    tot = np.concatenate(totals)
    cells = obs[keep].drop(columns="_row").reset_index(drop=True)
    cells["total_counts"] = tot
    denom = np.where(tot > 0, tot, 1)
    for i, name in enumerate(sets):
        cells[name] = np.round(np.log1p(1e4 * sc[:, i] / denom), 4)
    cells.to_csv(out / "eae_cell_scores.csv.gz", index=False, compression="gzip")
    print(f"wrote {out}: {len(groups)} pseudobulks, {len(cells):,} cells")


if __name__ == "__main__":
    main()
