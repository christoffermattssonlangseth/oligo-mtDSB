"""Tests for utils.robustness.

Run with pytest, or directly:  python tests/test_robustness.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from utils import robustness as rb  # noqa: E402

CTRL = ["c1", "c2", "c3", "c4", "c5", "c6"]
RESP = ["r1", "r2", "r3"]
SLIDE = pd.Series({**{c: "a" if i < 3 else "b" for i, c in enumerate(CTRL)}, "r1": "a", "r2": "b", "r3": "b"})


def test_slide_z_and_loro():
    v = pd.Series({"c1": 1, "c2": 2, "c3": 3, "c4": 11, "c5": 12, "c6": 13, "r1": 6, "r2": 20, "r3": 12.5})
    z = rb.slide_z(v, SLIDE, CTRL)
    assert abs(z["c2"]) < 1e-9 and abs(z["r1"] - 4) < 1e-9 and abs(z["r2"] - 8) < 1e-9   # slide offset removed
    out = rb.loro(v, RESP, CTRL)
    assert out["beyond_r2"] and not out["beyond_r1"] and out["n_beyond_control_range"] == 1


def test_heldout_signature_and_bh():
    rng = np.random.default_rng(0)
    sig = rng.normal(0, 1, 500)
    X = pd.DataFrame(rng.normal(0, 0.3, (9, 500)), index=CTRL + RESP)
    X.loc[RESP] += sig
    r = rb.heldout_signature_r(X, SLIDE, RESP, CTRL)
    assert r[RESP].min() > 0.8 > r[CTRL].max()
    adj = rb.bh([0.01, 0.04, 0.03, 0.5])
    assert np.allclose(adj, [0.04, 0.04 * 4 / 3, 0.04 * 4 / 3, 0.5])


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("ok", name)
