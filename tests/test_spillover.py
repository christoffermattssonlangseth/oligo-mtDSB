"""Tests for utils.spillover: recover a planted contamination fraction and separate own from spilled changes.

Run with pytest, or directly:  python tests/test_spillover.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from utils import spillover as sp  # noqa: E402


def make(alpha=0.05):
    """Two source types (OL, MG) and a target (AS) that picks up ``alpha`` of each neighbour's RNA."""
    genes = [f"ol{i}" for i in range(10)] + [f"mg{i}" for i in range(10)] + ["own_as", "induced"]
    own = pd.DataFrame(0.0, index=["OL", "MG", "AS"], columns=genes)
    own.loc["OL", [g for g in genes if g.startswith("ol")]] = 1000
    own.loc["MG", [g for g in genes if g.startswith("mg")]] = 1000
    own.loc["AS", "own_as"] = 100
    mix = pd.DataFrame(np.eye(3), index=own.index, columns=own.index)
    mix.loc["AS", ["OL", "MG"]] = alpha
    return own, mix


def test_alpha_recovered_and_changes_classified():
    own, mix = make()
    ctrl = mix @ own
    m = sp.source_markers(ctrl, ["OL", "MG"])
    assert set(m["OL"]) == {f"ol{i}" for i in range(10)} and "induced" not in m["MG"]   # silent in controls
    a = sp.contamination(ctrl, m)
    assert abs(a.loc["AS", "OL"] - 0.05) < 1e-9 and abs(a.loc["AS", "MG"] - 0.05) < 1e-9
    case = own.copy()
    case.loc["MG", "induced"] = 2000      # induced in microglia only -> appears in AS by spill-over
    case.loc["AS", "own_as"] = 300           # genuine astrocyte change
    delta = mix @ case - ctrl
    f = sp.spill_fraction(delta, sp.expected_spill(a, delta))
    assert f.loc["AS", "induced"] > 0.99 and f.loc["AS", "own_as"] < 0.01


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("ok", name)
