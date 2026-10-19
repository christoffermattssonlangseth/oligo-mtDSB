#!/usr/bin/env python
"""T-cell calls plus the polygons and transcripts around them, for analysis/08 (needs the h5ad + raw Xenium).

Writes to ``--outdir`` (default results/tcells):

* ``panel_summary.csv``      lymphocyte probes: counts, cells detected, most common cell class / region
* ``marker_counts.csv.gz``   per-cell panel counts + labels for every cell with >= 1 core T-cell transcript
* ``null_shuffle.csv``       calls after shuffling each core gene within section (``utils.tcells.shuffle_null``)
* ``zoom_cells.csv.gz``      cell polygons (+ labels) within ``--half-width`` um of every T cell
* ``zoom_nuclei.csv.gz``     nucleus polygons of the T cells
* ``zoom_transcripts.csv.gz`` T-cell gene transcripts (qv >= 20) in those windows

Transcripts are read with polars: pyarrow (<= 22) fails on most row groups of these transcripts.parquet
files ("Unexpected end of stream") although polars reads them cleanly. Run in an env with polars + anndata
(office Mac: ``sopa``).

Example (office Mac)
--------------------
    python scripts/export_tcells.py --h5ad <novae4.h5ad> --raw-dir /Volumes/moldiassd/oligo-mtDSB/raw_xenium/run4
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl
from scipy import sparse

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from utils import tcells as tc  # noqa: E402
from utils.region_time_de import Keys, read_minimal  # noqa: E402

LABELS = ["cell_class_updated", "cell_type", "cytetype_annotation_leiden_3"]
EXTRA_OBS = ["run", "cell_id", "x_centroid", "y_centroid", *LABELS]
TX_GENES = ["Cd3d", "Cd3e", "Cd3g", "Cd247", "Cd8a", "Cd8b1", "Cxcr6", "Lck", "Itk", "Cd2", "Ptprc"]


def in_windows(x, y, centres: np.ndarray, w: float) -> np.ndarray:
    keep = np.zeros(len(x), bool)
    for cx, cy in centres:
        keep |= (np.abs(x - cx) < w) & (np.abs(y - cy) < w)
    return keep


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--h5ad", required=True)
    p.add_argument("--raw-dir", required=True, help="folder holding the output-XETG... run folders")
    p.add_argument("--outdir", default=str(REPO / "results" / "tcells"))
    p.add_argument("--half-width", type=float, default=60.0)
    p.add_argument("--n-null", type=int, default=20)
    a = p.parse_args(argv)
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)

    adata = read_minimal(a.h5ad, layer="counts", keys=Keys(), extra_obs=EXTRA_OBS)
    obs = adata.obs.copy()
    for c in obs.columns:
        if isinstance(obs[c].dtype, pd.CategoricalDtype):
            obs[c] = obs[c].astype(str)
    genes = pd.Index(adata.var_names)
    panel = [g for g in tc.PANEL if g in genes]
    print("absent from panel:", [g for g in tc.PANEL if g not in genes])
    X = sparse.csc_matrix(adata.layers["counts"][:, genes.get_indexer(panel)])
    counts = pd.DataFrame(X.toarray().astype(np.int32), columns=panel, index=obs.index)
    total = np.asarray(adata.layers["counts"].sum(1)).ravel()
    del adata

    rows = []
    for g in panel:
        det = counts[g] > 0
        cls = obs.loc[det, "cell_class_updated"].value_counts(normalize=True)
        reg = obs.loc[det, "RBD_compartment_simplified"].value_counts(normalize=True)
        rows.append({"gene": g, "role": "core" if g in tc.CORE else "noisy" if g in tc.NOISY else "extra",
                     "total_counts": int(counts[g].sum()), "cells_ge1": int(det.sum()), "cells_ge2": int((counts[g] >= 2).sum()),
                     "top_class": cls.index[0] if len(cls) else "", "top_class_frac": round(cls.iloc[0], 3) if len(cls) else np.nan,
                     "top_region": reg.index[0] if len(reg) else "", "top_region_frac": round(reg.iloc[0], 3) if len(reg) else np.nan})
    pd.DataFrame(rows).to_csv(out / "panel_summary.csv", index=False)

    called = tc.call_tcells(counts)
    null = tc.shuffle_null(counts, obs["sample_id"], n_rep=a.n_null)
    pd.DataFrame({"rep": range(len(null)), "n_called": null}).to_csv(out / "null_shuffle.csv", index=False)
    print(f"called {called.sum()} T cells; shuffled null {null.mean():.1f} (max {null.max()})")

    keep = counts[tc.CORE].sum(1) > 0
    mc = obs.loc[keep, ["sample_id", "age", "condition", "RBD_compartment_simplified", *EXTRA_OBS]].join(counts[keep])
    mc.insert(0, "total_counts", total[keep.to_numpy()])
    mc["tcell"] = called[keep]
    mc.to_csv(out / "marker_counts.csv.gz", compression="gzip")

    labels = obs.set_index(["sample_id", "cell_id"])[LABELS]
    cells, nuclei, txs = [], [], []
    for smp, t in obs[called].groupby("sample_id"):
        run = Path(a.raw_dir) / t["run"].iloc[0]
        centres = t[["x_centroid", "y_centroid"]].to_numpy()
        cb = pl.read_parquet(run / "cell_boundaries.parquet").to_pandas()
        hit = cb.loc[in_windows(cb.vertex_x.to_numpy(), cb.vertex_y.to_numpy(), centres, a.half_width + 10), "cell_id"].unique()
        cb = cb[cb.cell_id.isin(hit)].assign(sample_id=smp)
        cells.append(cb.join(labels.loc[smp], on="cell_id"))
        nb = pl.read_parquet(run / "nucleus_boundaries.parquet").to_pandas()
        nuclei.append(nb[nb.cell_id.isin(set(t.cell_id))].assign(sample_id=smp))
        tx = (pl.scan_parquet(run / "transcripts.parquet")
              .select(["cell_id", "feature_name", "x_location", "y_location", "qv"])
              .filter(pl.col("feature_name").is_in(TX_GENES) & (pl.col("qv") >= 20)).collect().to_pandas())
        txs.append(tx[in_windows(tx.x_location.to_numpy(), tx.y_location.to_numpy(), centres, a.half_width)].assign(sample_id=smp))
        print(smp, len(t), "T cells,", len(hit), "cells in windows")
    pd.concat(cells).round(2).to_csv(out / "zoom_cells.csv.gz", index=False)
    pd.concat(nuclei).round(2).to_csv(out / "zoom_nuclei.csv.gz", index=False)
    pd.concat(txs).drop(columns="qv").round(2).to_csv(out / "zoom_transcripts.csv.gz", index=False)
    print("wrote", out)


if __name__ == "__main__":
    main()
