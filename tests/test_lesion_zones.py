"""Tests for utils.lesion_zones on a synthetic section with a demyelinated core.

Run with pytest, or directly:  python tests/test_lesion_zones.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from utils import lesion_zones as lz  # noqa: E402


def make(seed=0, n=6000):
    """1 x 1 mm section: myelin RNA everywhere except a 300 um core; core cells have a 5x higher event rate."""
    rng = np.random.default_rng(seed)
    x, y = rng.uniform(0, 1000, n), rng.uniform(0, 1000, n)
    core = np.hypot(x - 500, y - 500) < 150
    mbp = rng.poisson(np.where(core, 1, 30))
    return pd.DataFrame({"x_centroid": x, "y_centroid": y, "MBP": mbp, "SNAP25": 0, "total_counts": 300,
                         "core": core, "event": rng.random(n) < np.where(core, 0.25, 0.05)})


def test_zones_and_bootstrap():
    c = make()
    c["my"] = lz.local_density(c, ["MBP"], radius=50)
    c["gm"] = lz.local_density(c, ["SNAP25"], radius=50)
    c["zone"] = lz.assign_zones(c.my.to_numpy(), c.gm.to_numpy())
    assert (c[c.core & (np.hypot(c.x_centroid - 500, c.y_centroid - 500) < 80)].zone == "demyelinated").mean() > 0.95
    assert (c[~c.core & (np.hypot(c.x_centroid - 500, c.y_centroid - 500) > 300)].zone == "myelinated").mean() > 0.95
    r = lz.block_bootstrap_ratio(c, "event", "zone", "demyelinated", "myelinated", n_boot=200)
    assert 3 < r["ratio"] < 8 and r["ci_lo"] > 2 and r["ci_lo"] < r["ratio"] < r["ci_hi"]
    d = lz.border_distance(c)
    centre = np.hypot(c.x_centroid - 500, c.y_centroid - 500) < 50
    assert d[centre].min() > 400 and np.percentile(d, 1) < 60          # centre deep, outermost cells near the border
    c["SNAP25"] = 10                                   # everything neuron-rich -> grey matter, no lesion call
    c["gm"] = lz.local_density(c, ["SNAP25"], radius=50)
    assert (lz.assign_zones(c.my.to_numpy(), c.gm.to_numpy()) == "grey matter").all()


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("ok", name)
