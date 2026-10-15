"""Tests for utils.responder_mixture on synthetic cohorts with planted responders.

Run with pytest, or directly:  python tests/test_responder_mixture.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from utils import responder_mixture as rm  # noqa: E402


def make(shift=6.0, responders=("m21_0", "m60_1"), cohort_offset=3.0, seed=0, d=1):
    """2 cohorts x (3 controls + 3 mtDSB); responders shifted by ``shift`` SDs in every readout."""
    rng = np.random.default_rng(seed)
    rows = []
    for coh in ("21", "60"):
        for cond in ("c", "m"):
            for r in range(3):
                mid = f"{cond}{coh}_{r}"
                v = rng.normal(0, 1, d) + (cohort_offset if coh == "60" else 0)
                if mid in responders:
                    v = v + shift
                rows.append((mid, coh, cond == "c", *v))
    df = pd.DataFrame(rows, columns=["id", "cohort", "ctrl", *[f"r{j}" for j in range(d)]]).set_index("id")
    return df[[f"r{j}" for j in range(d)]], df.cohort, df.ctrl


def test_t_density_matches_scipy():
    S = np.array([[1.0, 0.4], [0.4, 2.0]])
    x = np.random.default_rng(1).normal(size=(7, 2))
    ref = stats.multivariate_t(np.zeros(2), S * 1.3, df=5).logpdf(x)
    assert np.allclose(rm._logpdf_t(x, S, 5, 1.3), ref)


def test_planted_responders_found():
    x, coh, ctrl = make()
    res = rm.fit(x, coh, ctrl)
    p = res["p_responder"]
    assert p[["m21_0", "m60_1"]].min() > 0.9
    assert p.drop(["m21_0", "m60_1"]).max() < 0.5
    dx = res["pi_grid"][1] - res["pi_grid"][0]
    assert abs(res["pi_density"].sum() * dx - 1) < 0.02


def test_cohort_offset_removed():
    # a large cohort offset must not turn the whole second cohort into responders
    x, coh, ctrl = make(cohort_offset=20.0)
    p = rm.fit(x, coh, ctrl)["p_responder"]
    assert p.drop(["m21_0", "m60_1"]).max() < 0.5


def test_null_gives_no_confident_responders():
    hits = []
    for seed in range(10):
        x, coh, ctrl = make(responders=(), seed=seed)
        hits.append((rm.fit(x, coh, ctrl)["p_responder"] > 0.9).sum())
        hits.append((rm.control_check(x, coh, ctrl) > 0.9).sum())
    assert np.mean(hits) < 0.5


def test_two_readouts_and_predictive_monotone():
    x, coh, ctrl = make(d=2)
    res = rm.fit(x, coh, ctrl, rm.Prior(n_delta=25))
    assert res["p_responder"][["m21_0", "m60_1"]].min() > 0.9
    z = np.c_[np.linspace(0.1, 8, 20), np.linspace(0.1, 8, 20)]
    pr = rm.predictive_p(res, z)
    # rises from "control-like" to "responder"; t tails allow a tiny dip far beyond the responders
    assert np.all(np.diff(pr) >= -0.01) and pr[0] < 0.5 < pr[-1]


if __name__ == "__main__":
    for name, f in list(globals().items()):
        if name.startswith("test_"):
            f()
            print("ok", name)
