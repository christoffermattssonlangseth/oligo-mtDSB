"""Tests for utils.mosaic on synthetic sections with a planted responding subset.

Run with pytest, or directly:  python tests/test_mosaic.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from utils import mosaic as mo  # noqa: E402


def make(pi=0.2, n=4000, seed=0):
    """One slide, 3 controls + 1 mosaic section (fraction ``pi`` of cells shifted up by 3 SD), two regions."""
    rng = np.random.default_rng(seed)
    rows = []
    for s in ["c1", "c2", "c3", "m1"]:
        for r in ["A", "B"]:
            x = rng.normal(0 if r == "A" else 2, 1, n)
            if s == "m1":
                x[: int(pi * n)] += 3
            rows += [(s, r, v) for v in x]
    df = pd.DataFrame(rows, columns=["sample_id", "region", "score"])
    group = pd.Series({"c1": "control", "c2": "control", "c3": "control", "m1": "mtDSB"})
    slide = pd.Series("s", index=group.index)
    return df, group, slide


def test_controls_uniform_and_mosaic_detected():
    df, group, slide = make()
    u = mo.loo_percentiles(df, "score", group, slide)
    for s in ["c1", "c2", "c3"]:
        v = u[df.sample_id == s]
        assert abs(v.mean() - 0.5) < 0.03 and 0.6 < mo.tail_enrichment(v, 0.99) < 1.6
    m = u[df.sample_id == "m1"]
    assert mo.tail_enrichment(m, 0.99) > 5
    assert 0.12 < mo.responding_fraction_lb(m) <= 0.22   # a lower bound, close to the planted 0.2


def test_spike_in_power_rises_with_fraction():
    df, group, slide = make(pi=0.3)
    u = mo.loo_percentiles(df, "score", group, slide)
    comp = mo.responding_component(u[df.sample_id == "m1"].to_numpy())
    ctrl = {s: u[df.sample_id == s].to_numpy() for s in ["c1", "c2", "c3"]}
    det = mo.spike_in_detection(ctrl, comp, [0.0, 0.05, 0.2], n_rep=10)
    assert det[0.0] < 0.5 and det[0.2] == 1.0 and det[0.05] >= det[0.0]


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("ok", name)
