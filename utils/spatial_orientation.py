"""Put every section in the same anatomical orientation before plotting.

The Xenium sections were imaged in different rotations, and controls and mtDSB sections are mirror images
of each other (mounted the other way round or opposite hemispheres). Image coordinates also have y pointing
down. Orientation is fitted from anatomy rather than by eye:

1. Region centroids (cell-count weighted mean position of each ``RBD_compartment_simplified`` label) are
   computed per section.
2. A reference section is put in the canonical frame: **dorsal up with a vertical midline** (the
   hypothalamus sits directly below the thalamus along the midline, so the hypothalamus -> thalamus vector
   points up) and **midline to the right** (hypothalamus to the right of striatum).
3. Every section is aligned to the reference by an orthogonal Procrustes fit of its region centroids
   (rotation, with a mirror flip allowed; no scaling), after centring on its cell centroid.

``fit_orientation`` returns one row per section with the 2x2 matrix, the centre, the rotation angle, whether
it was mirrored, and the centroid RMSD of the fit (a check: large values mean a poor match).
``apply_orientation`` adds ``x_oriented`` / ``y_oriented`` to a cell table.
"""

from __future__ import annotations

from typing import Sequence

import numpy as np
import pandas as pd

DEFAULT_REGIONS = ("Cortex", "Thalamus", "Hypothalamus", "Striatum", "Hippocampus", "Fiber tracts",
                   "Olfactory areas", "Ventricular system", "Pallidum")


def region_centroids(cells: pd.DataFrame, regions: Sequence[str] = DEFAULT_REGIONS, sample: str = "sample_id",
                     region: str = "RBD_compartment_simplified", x: str = "x_centroid", y: str = "y_centroid") -> pd.DataFrame:
    """Mean (x, y) per (section, region); y is flipped to point up."""
    c = cells.loc[cells[region].isin(regions), [sample, region, x, y]].copy()
    c[y] = -c[y]
    return c.groupby([sample, region], observed=True)[[x, y]].mean()


def _canonical(ref: pd.DataFrame) -> np.ndarray:
    """2x2 matrix putting the reference dorsal-up with a vertical midline, midline on the right (after centring)."""
    up = ref.loc["Thalamus"].to_numpy() - ref.loc["Hypothalamus"].to_numpy()
    ang = np.arctan2(up[1], up[0])
    rot = np.pi / 2 - ang
    R = np.array([[np.cos(rot), -np.sin(rot)], [np.sin(rot), np.cos(rot)]])
    p = (ref.to_numpy() - ref.to_numpy().mean(0)) @ R.T
    pr = pd.DataFrame(p, index=ref.index)
    if pr.loc["Hypothalamus", 0] < pr.loc["Striatum", 0]:
        R = np.diag([-1.0, 1.0]) @ R
    return R


def fit_orientation(cells: pd.DataFrame, reference: str, regions: Sequence[str] = DEFAULT_REGIONS,
                    sample: str = "sample_id", region: str = "RBD_compartment_simplified",
                    x: str = "x_centroid", y: str = "y_centroid") -> pd.DataFrame:
    cen = region_centroids(cells, regions, sample, region, x, y)
    centres = cells.assign(_y=-cells[y]).groupby(sample, observed=True)[[x, "_y"]].mean()
    ref = cen.loc[reference]
    R0 = _canonical(ref)
    target = (ref.to_numpy() - centres.loc[reference].to_numpy()) @ R0.T
    target = pd.DataFrame(target, index=ref.index)
    rows = []
    for sid in centres.index:
        src = cen.loc[sid]
        common = src.index.intersection(target.index)
        A = src.loc[common].to_numpy() - centres.loc[sid].to_numpy()
        B = target.loc[common].to_numpy()
        U, _, Vt = np.linalg.svd(A.T @ B)
        R = (U @ Vt).T  # maps a (column) vector: B ~ A @ R.T
        rmsd = float(np.sqrt(((A @ R.T - B) ** 2).sum(1).mean()))
        det = float(np.linalg.det(R))
        rows.append({sample: sid, "r11": R[0, 0], "r12": R[0, 1], "r21": R[1, 0], "r22": R[1, 1],
                     "cx": centres.loc[sid, x], "cy": centres.loc[sid, "_y"], "mirrored": det < 0,
                     "angle_deg": float(np.degrees(np.arctan2(R[1, 0] * np.sign(det), R[0, 0]))),
                     "rmsd": rmsd, "n_regions": len(common), "reference": reference})
    return pd.DataFrame(rows).set_index(sample)


def apply_orientation(cells: pd.DataFrame, orient: pd.DataFrame, sample: str = "sample_id",
                      x: str = "x_centroid", y: str = "y_centroid") -> pd.DataFrame:
    """Return a copy with ``x_oriented`` / ``y_oriented`` (same units, y up)."""
    out = cells.copy()
    out["x_oriented"] = np.nan
    out["y_oriented"] = np.nan
    for sid, idx in out.groupby(sample, observed=True).groups.items():
        if sid not in orient.index:
            continue
        o = orient.loc[sid]
        R = np.array([[o.r11, o.r12], [o.r21, o.r22]])
        P = np.column_stack([out.loc[idx, x].to_numpy() - o.cx, -out.loc[idx, y].to_numpy() - o.cy])
        Q = P @ R.T
        out.loc[idx, "x_oriented"] = Q[:, 0]
        out.loc[idx, "y_oriented"] = Q[:, 1]
    return out
