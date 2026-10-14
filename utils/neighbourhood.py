"""Non-cell-autonomous effects: are cells next to stressed oligodendrocytes (OLs) more activated?

For every cell of a target type (microglia, astrocytes, endothelial cells), ``local_ol_stress`` averages the
ISR score of all OLs within ``radius`` micrometres in the same section. ``coupling`` then asks, within each
section, whether the target cell's own readout (e.g. a DAM score) rises with that local OL stress.

Both variables are centred within (section, region) before correlating, so the association cannot come from
regions that are simply high in both (e.g. thalamus) or from section-level differences (slide, plane). The
controls are the null: their OLs are not stressed, so any residual coupling there is baseline co-patterning;
the test of interest is responders vs controls.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree


def local_ol_stress(cells: pd.DataFrame, target_types, radius: float = 50.0, score: str = "ISR",
                    x: str = "x_centroid", y: str = "y_centroid", ol_type: str = "Oligodendrocytes",
                    min_ols: int = 3) -> pd.DataFrame:
    """Mean OL ``score`` within ``radius`` of each target cell (same section).

    ``cells`` needs ``sample_id, cell_type, region``, coordinates and the score column. Returns the target
    cells with ``ol_stress`` (NaN if fewer than ``min_ols`` OLs in range) and ``n_ols``.
    """
    out = []
    for sid, sec in cells.groupby("sample_id", observed=True):
        ol = sec[sec.cell_type == ol_type]
        tg = sec[sec.cell_type.isin(target_types)]
        if ol.empty or tg.empty:
            continue
        tree = cKDTree(ol[[x, y]].to_numpy())
        nbrs = tree.query_ball_point(tg[[x, y]].to_numpy(), r=radius)
        s = ol[score].to_numpy()
        n = np.array([len(i) for i in nbrs])
        m = np.array([s[i].mean() if len(i) >= min_ols else np.nan for i in nbrs])
        out.append(tg.assign(ol_stress=m, n_ols=n))
    return pd.concat(out)


def _centre(df: pd.DataFrame, cols, by=("sample_id", "region")) -> pd.DataFrame:
    return df[list(cols)] - df.groupby(list(by), observed=True)[list(cols)].transform("mean")


def coupling(local: pd.DataFrame, readouts, group: pd.Series, min_cells: int = 50) -> pd.DataFrame:
    """Within-section, region-centred Pearson r between ``ol_stress`` and each readout, per (cell type, section).

    ``group`` maps sample_id -> control / non-responder / responder. Also returns the slope (readout units per
    unit of OL ISR score).
    """
    rows = []
    for (ct, sid), d in local.dropna(subset=["ol_stress"]).groupby(["cell_type", "sample_id"], observed=True):
        if len(d) < min_cells:
            continue
        for r in readouts.get(ct, []):
            c = _centre(d, ["ol_stress", r])
            xv, yv = c["ol_stress"].to_numpy(), c[r].to_numpy()
            if xv.std() == 0 or yv.std() == 0:
                continue
            rr = float(np.corrcoef(xv, yv)[0, 1])
            slope = float(np.polyfit(xv, yv, 1)[0])
            rows.append({"cell_type": ct, "readout": r, "sample_id": sid, "group": group.get(sid), "n_cells": len(d),
                         "r": rr, "slope": slope})
    return pd.DataFrame(rows)


def binned_response(local: pd.DataFrame, readout: str, cell_type: str, group: pd.Series, n_bins: int = 10) -> pd.DataFrame:
    """Readout (centred within section x region) by decile of region-centred local OL stress, per section."""
    d = local[(local.cell_type == cell_type)].dropna(subset=["ol_stress"]).copy()
    c = _centre(d, ["ol_stress", readout])
    d["stress_c"], d["readout_c"] = c["ol_stress"], c[readout]
    d["bin"] = d.groupby("sample_id", observed=True)["stress_c"].transform(
        lambda v: pd.qcut(v.rank(method="first"), n_bins, labels=False))
    b = d.groupby(["sample_id", "bin"], observed=True)["readout_c"].mean().reset_index()
    b["group"] = b.sample_id.map(group)
    return b
