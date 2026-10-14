"""Spatial structure of the oligodendrocyte (OL) stress response.

* ``border_distances``: distance of every cell to the nearest ventricle cell, meningeal-border cell and endothelial
  cell (a proxy for the nearest vessel), same section, micrometres.
* ``distance_coupling``: within each section, does OL stress fall with distance to a border? Both stress and
  log distance are centred within (section, region), so the test is not driven by which regions lie near borders.
* ``patches``: are stressed OLs clustered into patches? Stressed OLs closer than ``link`` um are joined; patch sizes
  are compared with a null in which the stressed label is shuffled among OLs of the same section and region
  (which keeps regional enrichment, e.g. thalamus, and tests clustering *within* regions).
* ``local_fraction``: fraction of cells of given types within ``radius`` of each centre cell.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

VENTRICLE = "Ventricular system"
MENINGES = "Meningeal border/glia limitans"


def border_distances(cells: pd.DataFrame, x: str = "x_centroid", y: str = "y_centroid") -> pd.DataFrame:
    out = []
    for sid, sec in cells.groupby("sample_id", observed=True):
        P = sec[[x, y]].to_numpy()
        d = pd.DataFrame(index=sec.index)
        for name, mask in [("d_ventricle", sec.region == VENTRICLE), ("d_meninges", sec.region == MENINGES),
                           ("d_vessel", sec.cell_type == "Endothelial Cells")]:
            ref = sec.loc[mask, [x, y]].to_numpy()
            d[name] = cKDTree(ref).query(P)[0] if len(ref) else np.nan
        out.append(d)
    return pd.concat(out)


def distance_coupling(ols: pd.DataFrame, score: str, dist_cols, group: pd.Series, exclude_regions=(VENTRICLE, MENINGES)) -> pd.DataFrame:
    """Per section: Spearman-like (rank) correlation of region-centred stress vs region-centred log distance."""
    d = ols[~ols.region.isin(exclude_regions)].copy()
    rows = []
    for col in dist_cols:
        d["_ld"] = np.log1p(d[col])
        c = d[[score, "_ld"]] - d.groupby(["sample_id", "region"], observed=True)[[score, "_ld"]].transform("mean")
        for sid, idx in d.groupby("sample_id", observed=True).groups.items():
            a, b = c.loc[idx, score].rank(), c.loc[idx, "_ld"].rank()
            rows.append({"sample_id": sid, "distance": col, "group": group.get(sid), "r": float(np.corrcoef(a, b)[0, 1]),
                         "n_OLs": len(idx)})
    return pd.DataFrame(rows)


def _components(P: np.ndarray, link: float) -> np.ndarray:
    if len(P) == 0:
        return np.array([], dtype=int)
    pairs = cKDTree(P).query_pairs(r=link, output_type="ndarray")
    n = len(P)
    A = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(n, n)) if len(pairs) else coo_matrix((n, n))
    return connected_components(A, directed=False)[1]


def patch_stats(P: np.ndarray, flag: np.ndarray, link: float, big: int = 10) -> dict:
    lab = _components(P[flag], link)
    sizes = np.bincount(lab) if len(lab) else np.array([0])
    n = int(flag.sum())
    return {"n_stressed": n, "n_patches": int((sizes > 0).sum()), "largest_patch": int(sizes.max()) if n else 0,
            "frac_in_big_patches": float(sizes[sizes >= big].sum() / n) if n else np.nan,
            "mean_patch_size_per_cell": float((sizes ** 2).sum() / n) if n else np.nan}


def patches(ols: pd.DataFrame, flag_col: str, link: float = 30.0, n_perm: int = 50, seed: int = 0,
            x: str = "x_centroid", y: str = "y_centroid", big: int = 10) -> pd.DataFrame:
    """Observed patch statistics per section vs label-shuffled null (within section x region)."""
    rng = np.random.default_rng(seed)
    rows = []
    for sid, sec in ols.groupby("sample_id", observed=True):
        P = sec[[x, y]].to_numpy()
        flag = sec[flag_col].to_numpy().astype(bool)
        obs = patch_stats(P, flag, link, big)
        reg = sec.region.to_numpy()
        null = []
        for _ in range(n_perm):
            f = flag.copy()
            for r in np.unique(reg):
                i = np.flatnonzero(reg == r)
                f[i] = rng.permutation(f[i])
            null.append(patch_stats(P, f, link, big))
        null = pd.DataFrame(null)
        row = {"sample_id": sid, **obs}
        for k in ("mean_patch_size_per_cell", "frac_in_big_patches", "largest_patch"):
            row[f"{k}_null"] = float(null[k].mean())
            row[f"{k}_ratio"] = obs[k] / row[f"{k}_null"] if row[f"{k}_null"] else np.nan
            row[f"{k}_perm_p"] = float((1 + (null[k] >= obs[k]).sum()) / (1 + n_perm))
        rows.append(row)
    return pd.DataFrame(rows).set_index("sample_id")


def local_fraction(cells: pd.DataFrame, centres: pd.DataFrame, target_types, radius: float = 50.0,
                   x: str = "x_centroid", y: str = "y_centroid") -> pd.Series:
    """For each centre cell: fraction of all cells within ``radius`` (same section) that are of ``target_types``."""
    out = []
    for sid, cen in centres.groupby("sample_id", observed=True):
        sec = cells[cells.sample_id == sid]
        tree = cKDTree(sec[[x, y]].to_numpy())
        is_t = sec.cell_type.isin(target_types).to_numpy()
        nb = tree.query_ball_point(cen[[x, y]].to_numpy(), r=radius)
        out.append(pd.Series([is_t[i].mean() if len(i) else np.nan for i in nb], index=cen.index))
    return pd.concat(out)


def correlogram(ols: pd.DataFrame, score: str, bands=((0, 25), (25, 50), (50, 100), (100, 200), (200, 400)),
                x: str = "x_centroid", y: str = "y_centroid", min_nb: int = 3) -> pd.DataFrame:
    """Spatial correlogram of a raw (unsmoothed) per-cell score among OLs.

    For each OL and distance band, the mean score of the *other* OLs in that annulus; the correlation (per section)
    between a cell's own score and its annulus mean, after centring both within (section, region). High values at
    larger distances = larger spatially coherent domains (patches); a fast decay = cell-by-cell (salt-and-pepper).
    Because no smoothing is applied, the clustering is not built in.
    """
    rows = []
    for sid, sec in ols.groupby("sample_id", observed=True):
        P = sec[[x, y]].to_numpy()
        s = sec[score].to_numpy()
        tree = cKDTree(P)
        reg = sec.region.to_numpy()
        prev = None
        outer_sets = {}
        for lo, hi in bands:
            outer_sets[hi] = tree.query_ball_point(P, r=hi)
        for lo, hi in bands:
            inner = outer_sets.get(lo) if lo > 0 else [[i] for i in range(len(P))]
            means = np.full(len(P), np.nan)
            for i in range(len(P)):
                ring = set(outer_sets[hi][i]) - set(inner[i]) - {i}
                if len(ring) >= min_nb:
                    means[i] = s[list(ring)].mean()
            df = pd.DataFrame({"own": s, "nb": means, "region": reg}).dropna()
            c = df[["own", "nb"]] - df.groupby("region")[["own", "nb"]].transform("mean")
            r = float(np.corrcoef(c["own"], c["nb"])[0, 1]) if len(c) > 50 else np.nan
            rows.append({"sample_id": sid, "band": f"{lo}-{hi}", "band_mid": (lo + hi) / 2, "r": r, "n": len(df)})
    return pd.DataFrame(rows)
