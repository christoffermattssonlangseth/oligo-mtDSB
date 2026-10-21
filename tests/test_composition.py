"""Tests for utils.composition on small synthetic sections.

Run with pytest, or directly:  python tests/test_composition.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from utils import composition as co  # noqa: E402


def test_exact_perm_extremes():
    v = pd.Series([5, 6, 7, 1, 1, 2, 1, 2, 1], index=list("abcdefghi"))
    d, p = co.exact_perm(v, "abc", "defghi")
    assert d > 4 and abs(p - 1 / 84) < 1e-9   # complete separation: only the observed split reaches |d|


def test_density_and_bins():
    rng = np.random.default_rng(0)
    # one 1 x 1 mm section: left half region A (OL-rich, 10 % microglia), right half region B (no OLs)
    n = 4000
    x, y = rng.uniform(0, 1000, n), rng.uniform(0, 1000, n)
    left = x < 500
    ct = np.where(left, rng.choice(["OL", "MG"], n, p=[0.9, 0.1]), "N")
    cells = pd.DataFrame({"sample_id": "s", "region": np.where(left, "A", "B"), "cell_type": ct, "x_centroid": x, "y_centroid": y})
    d = co.region_density(cells).set_index(["region", "cell_type"])
    assert abs(d.loc[("A", "MG"), "share"] - 10) < 2
    assert abs(d.loc[("A", "OL"), "area_mm2"] - 0.5) < 0.03          # half the section
    assert abs(d.loc[("A", "MG"), "density"] - 0.1 * n / 2 / 0.5) < 60
    b = co.ol_rich_bins(cells, cells.cell_type == "OL", cells.cell_type == "MG").set_index("ol_frac")
    assert b.loc["50%-100%", "density"] > 300 and b.loc["0%-10%", "density"] < 50
    q = cells.iloc[:5]
    assert (co.neighbour_counts(q, cells, radius=1e4) == n).all()


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("ok", name)
