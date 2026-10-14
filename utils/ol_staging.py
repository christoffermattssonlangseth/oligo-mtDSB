"""Order of events in responding oligodendrocytes: does the stress response come before myelin loss?

Each OL has a prespecified ISR score and a prespecified myelin score (``utils/gene_sets.py``). Both are
smoothed over the ``k`` nearest OLs in the section (single cells detect few transcripts) and converted to
percentiles of the controls on the same slide. Each OL is then classed as

    ISR-high   ISR percentile    > ``hi`` (default 90)
    myelin-low myelin percentile < ``lo`` (default 10)

If stress precedes myelin loss, responder tissue should contain many ISR-high OLs with normal myelin, while
myelin-low OLs should mostly be ISR-high (myelin loss rarely occurs without stress); the share of ISR-high OLs
that have also lost myelin should rise with the severity of the animal. This is a cross-sectional inference,
not a time course.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree


def smooth_scores(ols: pd.DataFrame, cols, k: int = 15, x: str = "x_centroid", y: str = "y_centroid") -> pd.DataFrame:
    out = ols.copy()
    for c in cols:
        out[f"{c}_smooth"] = np.nan
    for sid, idx in ols.groupby("sample_id", observed=True).groups.items():
        sub = ols.loc[idx]
        _, nbr = cKDTree(sub[[x, y]].to_numpy()).query(sub[[x, y]].to_numpy(), k=min(k, len(sub)))
        for c in cols:
            out.loc[idx, f"{c}_smooth"] = sub[c].to_numpy()[nbr].mean(1)
    return out


def control_percentile(df: pd.DataFrame, col: str, controls, by: str = "slide") -> pd.Series:
    """Percentile of ``col`` against control cells of the same ``by`` group."""
    pct = pd.Series(np.nan, index=df.index)
    for _, idx in df.groupby(by, observed=True).groups.items():
        sub = df.loc[idx]
        ref = np.sort(sub.loc[sub.sample_id.isin(controls), col].to_numpy())
        pct.loc[idx] = 100 * np.searchsorted(ref, sub[col].to_numpy(), side="right") / len(ref)
    return pct


def state_table(df: pd.DataFrame, isr_pct: str = "ISR_pct", myelin_pct: str = "MYELIN_pct",
                hi: float = 90, lo: float = 10) -> pd.DataFrame:
    """Per section: fractions in the four ISR x myelin states and the two conditional probabilities."""
    d = df.assign(isr_high=df[isr_pct] > hi, myelin_low=df[myelin_pct] < lo)
    rows = []
    for sid, g in d.groupby("sample_id", observed=True):
        ih, ml = g.isr_high, g.myelin_low
        rows.append({"sample_id": sid, "n_OLs": len(g),
                     "frac_isr_high": ih.mean(), "frac_myelin_low": ml.mean(),
                     "isr_high_myelin_normal": (ih & ~ml).mean(), "isr_high_myelin_low": (ih & ml).mean(),
                     "isr_normal_myelin_low": (~ih & ml).mean(),
                     "p_myelin_low_given_isr_high": (ih & ml).sum() / max(ih.sum(), 1),
                     "p_isr_high_given_myelin_low": (ih & ml).sum() / max(ml.sum(), 1)})
    return pd.DataFrame(rows).set_index("sample_id")


def myelin_by_isr_bin(df: pd.DataFrame, bins=(0, 50, 75, 90, 95, 99, 100.01), isr_pct: str = "ISR_pct",
                      myelin: str = "MYELIN_pct") -> pd.DataFrame:
    """Mean myelin percentile per ISR-percentile bin, per section."""
    d = df.assign(isr_bin=pd.cut(df[isr_pct], bins=list(bins), right=False))
    return d.groupby(["sample_id", "isr_bin"], observed=True)[myelin].agg(["mean", "size"]).reset_index()
