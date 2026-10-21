"""Tissue composition per mouse: cell densities per region, region-free white-matter bins, local neighbour counts.

Densities use the area a region covers in the section: the number of ``bin`` x ``bin`` um grid squares holding at
least one cell of that region. Shares (% of cells) are also reported, but a share moves when *other* cell types
change, so density is the primary readout. Because region labels come from pooled neighbourhood transcripts and
could themselves shift with disease, ``ol_rich_bins`` gives a label-free alternative: grid squares stratified by
their oligodendrocyte fraction (a white-matter proxy).

Group comparisons are exact permutation tests over mice (``exact_perm``), as elsewhere in the repo: with 3 vs 6
mice there are 84 relabellings and the smallest two-sided p is 1/84 = 0.012.
"""

from __future__ import annotations

import itertools

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree


def exact_perm(values: pd.Series, a, b) -> tuple[float, float]:
    """Mean(a) - mean(b) and its two-sided exact permutation p over all relabellings of a + b."""
    v = values.loc[list(a) + list(b)].to_numpy(float)
    na = len(a)
    obs = v[:na].mean() - v[na:].mean()
    null = np.array([v[list(k)].mean() - np.delete(v, list(k)).mean() for k in itertools.combinations(range(len(v)), na)])
    return float(obs), float((np.abs(null) >= abs(obs) - 1e-12).mean())


def occupied_area(cells: pd.DataFrame, by=("sample_id", "region"), bin: float = 50.0,
                  x: str = "x_centroid", y: str = "y_centroid") -> pd.Series:
    """Area (mm2) per ``by`` group: number of occupied ``bin`` um grid squares x bin^2."""
    sq = pd.DataFrame({**{k: cells[k] for k in by}, "bx": (cells[x] // bin).astype(int), "by_": (cells[y] // bin).astype(int)})
    return sq.drop_duplicates().groupby(list(by)).size() * bin * bin / 1e6


def region_density(cells: pd.DataFrame, cell_type: str = "cell_type", region: str = "region", bin: float = 50.0) -> pd.DataFrame:
    """Long table per (sample, region, cell type): n, share of the region's cells (%) and density (cells / mm2)."""
    area = occupied_area(cells, ("sample_id", region), bin).rename("area_mm2")
    n = cells.groupby(["sample_id", region, cell_type], observed=True).size().rename("n").reset_index()
    n = n.join(area, on=["sample_id", region])
    n["share"] = 100 * n.n / n.groupby(["sample_id", region]).n.transform("sum")
    n["density"] = n.n / n.area_mm2
    return n


def ol_rich_bins(cells: pd.DataFrame, is_ol: pd.Series, is_target: pd.Series, bin: float = 100.0, min_cells: int = 20,
                 edges=(0, 0.1, 0.3, 0.5, 1.0), x: str = "x_centroid", y: str = "y_centroid") -> pd.DataFrame:
    """Target-cell density (per mm2) in grid squares stratified by OL fraction, per sample. Region labels unused."""
    d = pd.DataFrame({"sample_id": cells.sample_id, "bx": (cells[x] // bin).astype(int), "by": (cells[y] // bin).astype(int),
                      "ol": is_ol.to_numpy(), "tg": is_target.to_numpy()})
    b = d.groupby(["sample_id", "bx", "by"]).agg(n=("ol", "size"), ol=("ol", "sum"), tg=("tg", "sum")).reset_index()
    b = b[b.n >= min_cells]
    labels = [f"{lo:.0%}-{hi:.0%}" for lo, hi in zip(edges[:-1], edges[1:])]
    b["ol_frac"] = pd.cut(b.ol / b.n, list(edges), labels=labels, include_lowest=True)
    out = b.groupby(["sample_id", "ol_frac"], observed=True).agg(n_bins=("tg", "size"), tg=("tg", "sum"))
    out["density"] = out.tg / out.n_bins / (bin * bin / 1e6)
    return out.reset_index()


def neighbour_counts(query: pd.DataFrame, target: pd.DataFrame, radius: float = 50.0,
                     x: str = "x_centroid", y: str = "y_centroid") -> np.ndarray:
    """Number of ``target`` cells within ``radius`` um of each ``query`` cell (same section assumed)."""
    if target.empty:
        return np.zeros(len(query), int)
    return cKDTree(target[[x, y]].to_numpy()).query_ball_point(query[[x, y]].to_numpy(), radius, return_length=True)
