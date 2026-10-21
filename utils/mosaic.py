"""Is the OL response mosaic? Per-cell control percentiles, tail enrichment and a bound on the responding fraction.

If induction worked in only some oligodendrocytes, a "non-responder" would hide a small set of strongly stressed
OLs on top of a control-like majority. Each OL is placed at its percentile among the same-slide, same-region control
OLs (leaving its own section out, so controls are scored fairly). Under the null the percentiles are uniform; a
responding subset appears as an excess in the upper tail.

* ``loo_percentiles``    per-cell control percentile, ties broken at random
* ``tail_enrichment``    observed / expected share of cells above the q-th control percentile
* ``responding_fraction_lb``  nonparametric lower bound on the responding share: with F = (1 - pi) U + pi G,
                              S(u) <= (1 - pi)(1 - u) + pi, so pi >= max_u (S(u) - (1 - u)) / u
* ``spike_in_detection`` how often a responding share f, planted into control sections, exceeds every other control
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def loo_percentiles(df: pd.DataFrame, score: str, group: pd.Series, slide: pd.Series, by: str = "region",
                    min_ref: int = 50, seed: int = 0) -> pd.Series:
    """Percentile of each cell among control cells of the same slide and ``by`` stratum, own section left out.

    ``df`` needs ``sample_id``, ``by`` and ``score``; ``group`` and ``slide`` map sample_id -> group / slide.
    Strata with fewer than ``min_ref`` reference cells give NaN.
    """
    rng = np.random.default_rng(seed)
    out = pd.Series(np.nan, index=df.index)
    for _, g in df.groupby([df.sample_id.map(slide), by]):
        ctr = g[g.sample_id.map(group) == "control"]
        for s, gs in g.groupby("sample_id"):
            ref = np.sort(ctr.loc[ctr.sample_id != s, score].dropna().to_numpy())
            if len(ref) < min_ref:
                continue
            v = gs[score].to_numpy()
            lo, hi = np.searchsorted(ref, v, "left"), np.searchsorted(ref, v, "right")
            out.loc[gs.index] = (lo + rng.uniform(size=len(v)) * (hi - lo)) / len(ref)
    return out


def tail_enrichment(u, q: float) -> float:
    u = np.asarray(u, float)
    u = u[~np.isnan(u)]
    return float((u > q).mean() / (1 - q))


def responding_fraction_lb(u, grid=np.linspace(0.5, 0.999, 200)) -> float:
    u = np.asarray(u, float)
    u = u[~np.isnan(u)]
    S = np.array([(u > g).mean() for g in grid])
    return float(np.max((S - (1 - grid)) / grid))


def responding_component(u_resp, bins: int = 50):
    """Histogram of the responding component in responder percentiles: density minus its floor (the uniform part).

    Crude (assumes the least-populated bin holds no responding cells), so it is used only to plant realistic
    responding cells in ``spike_in_detection``. Returns (bin edges, probabilities).
    """
    h, e = np.histogram(np.asarray(u_resp)[~np.isnan(u_resp)], bins=bins, range=(0, 1))
    d = h / h.sum() * bins
    comp = np.clip(d - d.min(), 0, None)
    return e, comp / comp.sum()


def spike_in_detection(u_by_control: dict, component, fractions, q: float = 0.99, n_rep: int = 20, seed: int = 0) -> pd.Series:
    """Share of (control, replicate) pairs in which planting fraction f of responding cells pushes the section's
    tail enrichment above the maximum of the *other* controls."""
    rng = np.random.default_rng(seed)
    e, p = component
    width = e[1] - e[0]
    te = {s: tail_enrichment(v, q) for s, v in u_by_control.items()}
    out = {}
    for f in fractions:
        hits = []
        for s, v in u_by_control.items():
            v = np.asarray(v)[~np.isnan(v)]
            thr = max(t for o, t in te.items() if o != s)
            for _ in range(n_rep):
                w = v.copy()
                k = int(round(f * len(w)))
                b = rng.choice(len(p), k, p=p)
                w[rng.choice(len(w), k, replace=False)] = e[b] + rng.uniform(size=k) * width
                hits.append(tail_enrichment(w, q) > thr)
        out[f] = float(np.mean(hits))
    return pd.Series(out, name=f"detected_q{q}")
