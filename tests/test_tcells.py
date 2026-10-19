"""Tests for utils.tcells on synthetic counts with planted T cells.

Run with pytest, or directly:  python tests/test_tcells.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from utils import tcells as tc  # noqa: E402


def make(n=20000, n_t=40, seed=0):
    """Sparse ambient noise in every core gene plus ``n_t`` planted T cells with 2-4 transcripts in each of 3 genes."""
    rng = np.random.default_rng(seed)
    m = rng.poisson(0.003, size=(n, len(tc.CORE)))
    planted = rng.choice(n, n_t, replace=False)
    for i in planted:
        m[i, rng.choice(len(tc.CORE), 3, replace=False)] += rng.integers(2, 5, 3)
    groups = pd.Series(np.repeat(["s1", "s2"], n // 2))
    return pd.DataFrame(m, columns=tc.CORE), groups, planted


def test_call_rules():
    df = pd.DataFrame([[3, 0, 0, 0, 0, 0, 0, 0, 0], [2, 1, 0, 0, 0, 0, 0, 0, 0], [1, 1, 0, 0, 0, 0, 0, 0, 0]], columns=tc.CORE)
    assert tc.call_tcells(df).tolist() == [False, True, False]   # one gene only / 3 tx in 2 genes / too few
    assert tc.call_tcells(df, min_tx=2).tolist() == [False, True, True]


def test_recovers_planted_and_null_is_small():
    df, groups, planted = make()
    called = tc.call_tcells(df)
    assert set(planted) <= set(np.flatnonzero(called))
    assert called.sum() <= len(planted) + 2
    null = tc.shuffle_null(df, groups, n_rep=5)
    assert null.max() < len(planted) / 4   # co-expression broken -> few calls


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("ok", name)
