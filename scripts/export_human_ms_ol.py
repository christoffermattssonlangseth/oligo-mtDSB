#!/usr/bin/env python
"""OL-lineage nuclei of human MS single-nucleus atlases: per-cell gene-set scores and per-sample pseudobulks (analysis/16).

Reads harmonised h5ad files (``obs``: ``cell_type_coarse``, ``condition_original``, ``sample``, optional subtype
column; raw counts in ``layers/counts``) in row chunks and keeps Oligodendrocyte and OPC nuclei. Writes to ``--out``:

* ``human_ol_cells.csv.gz``      per nucleus: dataset, sample, lesion, cell type, subtype, total counts, the scores
                                 ``log1p(1e4 * set counts / total)`` of utils/human_ms.gene_sets, and raw counts of the
                                 mtDSB core genes
* ``human_ol_pseudobulk.csv.gz`` summed counts per (dataset, sample, lesion, cell type) for human orthologs of the
                                 mouse panel genes expressed in OLs (analysis/10 table)
* ``human_ol_pseudobulk_meta.csv``

Example (office Mac, env cellcharter-env; data on /Volumes/jamboree)
--------------------------------------------------------------------
    python scripts/export_human_ms_ol.py --data-dir /Volumes/jamboree/OligoC4b_public/processed \
        --out /Volumes/moldiassd/oligo-mtDSB/human_ms
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

from utils.human_ms import gene_sets, to_human  # noqa: E402

DATASETS = {"Jakel2019": ("Jakel2019_MS_human_snRNA.h5ad", "cell_type_original"),
            "Schirmer2019": ("Schirmer2019_MS_human_snRNA.h5ad", "cell_type_original"),
            "Absinta2021": ("Absinta2021_MS_human_snRNA.h5ad", "cell_type_original"),
            "LermaMartin2024": ("LermaMartin2024_MS_human_snRNA.h5ad", "cell_type_original")}
KEEP = ["Oligodendrocyte", "OPC"]


def obs_col(f, k):
    g = f["obs"][k]
    if isinstance(g, h5py.Group):
        cats = np.array([c.decode() if isinstance(c, bytes) else str(c) for c in g["categories"][:]], dtype=object)
        codes = g["codes"][:]
        return np.where(codes >= 0, cats[np.clip(codes, 0, None)], None)
    v = g[:]
    return np.array([x.decode() if isinstance(x, bytes) else x for x in v], dtype=object)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--data-dir", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--chunk", type=int, default=50_000)
    a = p.parse_args(argv)
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    sets = gene_sets(REPO)
    core = sets["MTDSB_CORE"]
    # genes expressed in mtDSB / EAE OLs on the mouse 5k panel (analysis/10), as human orthologs
    panel = to_human(pd.read_csv(REPO / "results" / "eae_comparison" / "tables" / "ol_log2fc_mtdsb_vs_eae.csv", index_col=0).index)
    cells, pbs, metas = [], [], []
    for name, (fname, sub_col) in DATASETS.items():
        with h5py.File(Path(a.data_dir) / fname, "r") as f:
            v = f["var"]; genes = pd.Index([x.decode() if isinstance(x, bytes) else x for x in v[v.attrs["_index"]][:]])
            obs = pd.DataFrame({"sample": obs_col(f, "sample"), "lesion": obs_col(f, "condition_original"),
                                "cell_type": obs_col(f, "cell_type_coarse"), "subtype": obs_col(f, sub_col)})
            keep = obs.cell_type.isin(KEEP).to_numpy()
            S = np.zeros((len(sets), len(genes)), np.float32)
            for i, gs in enumerate(sets.values()):
                S[i, genes.get_indexer([g for g in gs if g in genes])] = 1
            pidx = genes.get_indexer([g for g in panel if g in genes]); pnames = genes[pidx]
            cidx = genes.get_indexer([g for g in core if g in genes]); cnames = genes[cidx]
            key = obs.loc[keep, ["sample", "lesion", "cell_type"]].astype(str).agg("|".join, axis=1)
            groups = pd.Index(sorted(key.unique())); gcode = np.full(len(obs), -1); gcode[keep] = groups.get_indexer(key)
            pb = np.zeros((len(groups), len(pidx)), np.int64)
            g = f["layers/counts"]; indptr = g["indptr"][:]
            sc, tot, cc = [], [], []
            for start in range(0, len(obs), a.chunk):
                stop = min(start + a.chunk, len(obs))
                k = np.flatnonzero(keep[start:stop])
                if len(k) == 0:
                    continue
                lo, hi = indptr[start], indptr[stop]
                X = sparse.csr_matrix((g["data"][lo:hi], g["indices"][lo:hi], indptr[start:stop + 1] - lo), shape=(stop - start, len(genes)))[k]
                M = sparse.csr_matrix((np.ones(len(k)), (gcode[start + k], np.arange(len(k)))), shape=(len(groups), len(k)))
                pb += (M @ X[:, pidx]).toarray().round().astype(np.int64)
                sc.append(np.asarray(X @ S.T)); tot.append(np.asarray(X.sum(1)).ravel()); cc.append(X[:, cidx].toarray())
        c = obs[keep].reset_index(drop=True).assign(dataset=name)
        t = np.concatenate(tot); c["total_counts"] = t
        scores = np.vstack(sc); denom = np.where(t > 0, t, 1)
        for i, s in enumerate(sets):
            c[s] = np.round(np.log1p(1e4 * scores[:, i] / denom), 4)
        c[list(cnames)] = np.vstack(cc).round().astype(np.int32)
        cells.append(c)
        meta = pd.DataFrame([k.split("|") for k in groups], columns=["sample", "lesion", "cell_type"]).assign(dataset=name)
        meta["n_cells"] = np.bincount(gcode[keep], minlength=len(groups))
        metas.append(meta); pbs.append(pd.DataFrame(pb, columns=pnames))
        print(name, f"{len(c):,} OL-lineage nuclei, {len(groups)} pseudobulks, {len(pnames)} panel genes", flush=True)
    meta = pd.concat(metas, ignore_index=True); meta.index.name = "pb_id"
    meta.to_csv(out / "human_ol_pseudobulk_meta.csv")
    pd.concat(pbs, ignore_index=True).fillna(0).astype(np.int64).set_index(meta.index).to_csv(out / "human_ol_pseudobulk.csv.gz", compression="gzip")
    pd.concat(cells, ignore_index=True).to_csv(out / "human_ol_cells.csv.gz", index=False, compression="gzip")
    print("wrote", out)


if __name__ == "__main__":
    main()
