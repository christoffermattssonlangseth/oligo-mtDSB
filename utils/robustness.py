"""Robustness checks for per-mouse results (analysis/15): one responder at a time, circularity, multiple testing.

* ``slide_z``       per-mouse value as a z-score against the controls of the same slide (age), using the pooled
                    within-slide control SD, so the two slides' baselines do not inflate the spread
* ``loro``          leave-one-responder-out: the exact permutation p with each responder dropped (2 vs 6 mice,
                    smallest p = 1/28 = 0.036) and whether every responder lies beyond the control range
* ``heldout_signature_r``  circularity check: correlation of a mouse's log2FC profile (genes used to call responders
                    removed) with the mean profile of the *other* responders. Control mice give the null, since
                    both profiles share the same control baseline
* ``bh``            Benjamini-Hochberg adjusted p-values
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from utils.composition import exact_perm


def slide_z(values: pd.Series, slide: pd.Series, controls) -> pd.Series:
    v = values.astype(float)
    ctrl = [c for c in controls if c in v.index and pd.notna(v[c])]
    mu = v[ctrl].groupby(slide[ctrl]).mean()
    resid = v[ctrl] - slide[ctrl].map(mu)
    sd = np.sqrt((resid ** 2).sum() / max(len(ctrl) - mu.size, 1))
    return (v - slide[v.index].map(mu)) / (sd if sd > 0 else np.nan)


def loro(values: pd.Series, responders, controls, direction: int = 1) -> dict:
    """``direction`` +1 if responders are expected above controls, -1 below."""
    v = values.astype(float)
    _, p_all = exact_perm(v, responders, controls)
    p_drop = {r: exact_perm(v, [x for x in responders if x != r], controls)[1] for r in responders}
    lim = v[list(controls)].max() if direction > 0 else v[list(controls)].min()
    beyond = {r: bool((v[r] - lim) * direction > 0) for r in responders}
    return {"p_all": p_all, "p_worst_drop_one": max(p_drop.values()), "dropped_for_worst": max(p_drop, key=p_drop.get),
            "n_beyond_control_range": sum(beyond.values()), **{f"beyond_{r}": b for r, b in beyond.items()}}


def heldout_signature_r(log_expr: pd.DataFrame, slide: pd.Series, responders, controls) -> pd.Series:
    """``log_expr``: mice x genes (log2 CPM). Returns r per mouse (see module docstring)."""
    def profile(m, ref):
        cc = [c for c in ref if slide[c] == slide[m]]
        return log_expr.loc[m] - log_expr.loc[cc].mean()
    out = {}
    for m in log_expr.index:
        ref = [c for c in controls if c != m]
        sig = pd.concat([profile(r, ref) for r in responders if r != m], axis=1).mean(1)
        out[m] = float(np.corrcoef(profile(m, ref), sig)[0, 1])
    return pd.Series(out, name="heldout_signature_r")


def bh(p) -> np.ndarray:
    p = np.asarray(p, float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(adj, 0, 1)
    return out
