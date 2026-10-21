"""Lesion zones in human MS Xenium sections from local myelin RNA, and cell-level rates with spatial uncertainty.

There is no pathology annotation for these sections, so zones are drawn from the transcripts themselves:

* ``local_density``  per cell, the summed counts of a gene set in all cells within ``radius`` um, per 1,000 counts
* ``assign_zones``   myelin density relative to the section's well-myelinated tissue (its 90th percentile):
                     < ``lo`` demyelinated, < ``hi`` edge, else myelinated. Neuron-rich neighbourhoods (local SNAP25
                     above ``gm_max``) are grey matter and get no zone, so low myelin there is not called a lesion
* ``border_distance``  distance of each cell to the tissue border (occupancy grid, closed and hole-filled), because
                     cells near a section edge can show stress transcripts for reasons unrelated to lesions
* ``block_bootstrap_ratio``  rate of a per-cell event (e.g. >= 1 TRIB3 transcript) in one zone over another, with a CI
                     from resampling square spatial blocks, which respects the spatial correlation of neighbouring cells

The zones use myelin genes, so myelin genes of the OLs themselves are partly circular; the stress genes are not.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import ndimage
from scipy.spatial import cKDTree


def local_density(cells: pd.DataFrame, genes, radius: float = 100.0, x: str = "x_centroid", y: str = "y_centroid",
                  total: str = "total_counts") -> np.ndarray:
    """Within one section: counts of ``genes`` summed over neighbours within ``radius``, per 1,000 neighbour counts."""
    P = cells[[x, y]].to_numpy()
    nb = cKDTree(P).query_ball_point(P, radius)
    g = cells[list(genes)].sum(1).to_numpy(float)
    t = cells[total].to_numpy(float)
    return np.array([g[i].sum() / max(t[i].sum(), 1.0) * 1e3 for i in nb])


def assign_zones(myelin_local: np.ndarray, gm_local: np.ndarray, lo: float = 0.25, hi: float = 0.6, gm_max: float = 1.0,
                 ref_q: float = 90) -> np.ndarray:
    wm = gm_local <= gm_max
    ref = np.percentile(myelin_local[wm], ref_q) if wm.any() else np.nan
    r = myelin_local / ref
    z = np.select([r < lo, r < hi], ["demyelinated", "edge"], "myelinated").astype(object)
    z[~wm] = "grey matter"
    return z


def border_distance(cells: pd.DataFrame, bin: float = 25.0, close_iter: int = 3, x: str = "x_centroid",
                    y: str = "y_centroid") -> np.ndarray:
    """Within one section: distance (um) from each cell to the nearest empty grid square outside the tissue."""
    gx = (cells[x] // bin).astype(int).to_numpy(); gy = (cells[y] // bin).astype(int).to_numpy()
    occ = np.zeros((gx.max() + 3, gy.max() + 3), bool)
    occ[gx + 1, gy + 1] = True
    occ = ndimage.binary_fill_holes(ndimage.binary_closing(occ, iterations=close_iter))
    return ndimage.distance_transform_edt(occ)[gx + 1, gy + 1] * bin


def block_bootstrap_ratio(cells: pd.DataFrame, event: str, zone: str, num: str, den: str, block: float = 250.0,
                          n_boot: int = 500, seed: int = 0, x: str = "x_centroid", y: str = "y_centroid") -> dict:
    """Rate of boolean ``event`` among cells in zone ``num`` divided by the rate in zone ``den`` (one section),
    with a 95 % CI from resampling ``block`` um squares with replacement."""
    rng = np.random.default_rng(seed)
    b = (cells[x] // block).astype(int).astype(str) + "_" + (cells[y] // block).astype(int).astype(str)
    d = pd.DataFrame({"b": b.to_numpy(), "z": cells[zone].to_numpy(), "e": cells[event].to_numpy(bool)})
    agg = d[d.z.isin([num, den])].groupby(["b", "z"]).e.agg(["sum", "size"]).unstack(fill_value=0)
    blocks = agg.index.to_numpy()

    def ratio(a):
        rn = a[("sum", num)].sum() / max(a[("size", num)].sum(), 1)
        rd = a[("sum", den)].sum() / max(a[("size", den)].sum(), 1)
        return rn / rd if rd > 0 else np.nan, rn, rd

    est, rn, rd = ratio(agg)
    boots = [ratio(agg.loc[rng.choice(blocks, len(blocks), replace=True)])[0] for _ in range(n_boot)]
    lo, hi = np.nanpercentile(boots, [2.5, 97.5])
    return {"rate_num": rn, "rate_den": rd, "ratio": est, "ci_lo": lo, "ci_hi": hi,
            "n_num": int(agg[("size", num)].sum()), "n_den": int(agg[("size", den)].sum())}
