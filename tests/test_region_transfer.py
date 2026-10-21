"""Tests for utils.region_transfer: a region absorbed by disease genes is recovered from disease-free genes.

Run with pytest, or directly:  python tests/test_region_transfer.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from utils import region_transfer as rt  # noqa: E402


def section(name, rng, n=1500):
    """Left third region A, middle B, right C; genes 0-9 mark A, 10-19 B, 20-29 C (Poisson)."""
    x, y = rng.uniform(0, 300, n), rng.uniform(0, 100, n)
    reg = np.select([x < 100, x < 200], ["A", "B"], "C")
    lam = np.full((n, 30), 0.2)
    for j, r in enumerate("ABC"):
        lam[reg == r, 10 * j:10 * j + 10] = 3
    return pd.DataFrame({"sample_id": name, "region": reg, "x": x, "y": y}), rng.poisson(lam)


def test_neighbour_pool_and_transfer():
    rng = np.random.default_rng(0)
    parts = [section(s, rng) for s in ["c1", "c2", "c3", "d1"]]
    obs = pd.concat([p[0] for p in parts], ignore_index=True)
    X = sparse.csr_matrix(np.vstack([p[1] for p in parts]))
    P = rt.neighbour_pool(X, obs.sample_id.to_numpy(), obs[["x", "y"]].to_numpy(), k=5)
    assert np.allclose(P.sum(), X.sum() * 6, rtol=0.01) and P.shape == X.shape
    # in the "disease" section the original map called region B as C
    obs["orig"] = np.where((obs.sample_id == "d1") & (obs.region == "B"), "C", obs.region)
    N = rt.log_normalise(P)
    m = rt.RegionTransfer(n_pcs=5)
    ref = obs.sample_id.isin(["c1", "c2", "c3"]).to_numpy()
    Z = m.embed_fit(N, ref)
    acc = m.loco_accuracy(Z, obs.orig.to_numpy(), obs.sample_id.to_numpy(), ["c1", "c2", "c3"], ref)
    assert (acc > 0.9).all()
    lab, _ = m.fit(Z, obs.orig.to_numpy(), ref).predict(Z)
    d = (obs.sample_id == "d1") & (obs.region == "B")
    assert (lab[d] == "B").mean() > 0.9


def test_select_anatomy_genes_excludes():
    rng = np.random.default_rng(1)
    idx = pd.MultiIndex.from_product([["A", "B"], ["s1", "s2", "s3"]])
    L = pd.DataFrame(rng.normal(8, 0.1, (6, 4)), index=idx, columns=["g1", "g2", "g3", "g4"])
    L.loc["A", "g1"] = 12; L.loc["A", "g2"] = 11
    F = rt.select_anatomy_genes(L, pd.Series(idx.get_level_values(0)), exclude={"g1"}, n=2)
    assert F.index[0] == "g2" and "g1" not in F.index


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("ok", name)
