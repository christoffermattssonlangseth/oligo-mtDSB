"""How much of a cell type's expression change could be RNA spilled over from neighbouring cell types?

Xenium segmentation assigns some transcripts of neighbouring cells to the wrong cell, so a strong change in one cell
type (myelin genes down in OLs, Cst7 up in microglia) also appears, diluted, in every cell type next to it. A linear
mixing model makes this quantitative on pseudobulk CPM:

    observed_t(g) = own_t(g) + sum_s alpha[t, s] * own_s(g)

* ``source_markers``  genes expressed almost only in one source type (>= ``fold`` x every other source type, in
                      controls). Their level in a target type can only be contamination.
* ``contamination``   alpha[t, s] = median over the markers of s of CPM_t / CPM_s (controls).
* ``expected_spill``  change in target t predicted from the changes in the *other* source types:
                      sum_{s != t} alpha[t, s] * delta_s(g).

A target-type change is called own when spill-over predicts less than half of it (same sign). The model is first
order: it ignores that contamination differs by region and that disease changes cell density (more microglia means
more microglial spill-over), so it is a filter, not a correction.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def source_markers(cpm_ctrl: pd.DataFrame, sources, fold: float = 5.0, min_cpm: float = 50.0) -> dict:
    """``cpm_ctrl``: cell types x genes (control mean CPM). Marker genes per source type."""
    out = {}
    for s in sources:
        others = cpm_ctrl.loc[[o for o in sources if o != s]]
        ok = (cpm_ctrl.loc[s] > min_cpm) & (cpm_ctrl.loc[s] >= fold * others.max())
        out[s] = cpm_ctrl.columns[ok].tolist()
    return out


def contamination(cpm_ctrl: pd.DataFrame, markers: dict) -> pd.DataFrame:
    """alpha: target types (rows) x source types (columns); 1 on the diagonal for source types."""
    return pd.DataFrame({s: cpm_ctrl[g].div(cpm_ctrl.loc[s, g], axis=1).median(1) for s, g in markers.items() if g})


def expected_spill(alpha: pd.DataFrame, delta: pd.DataFrame) -> pd.DataFrame:
    """``delta``: cell types x genes (CPM change). Expected change per target type from the other source types."""
    a = alpha.copy()
    for t in a.index.intersection(a.columns):
        a.loc[t, t] = 0.0
    return pd.DataFrame(a.to_numpy() @ delta.loc[a.columns].to_numpy(), index=a.index, columns=delta.columns)


def spill_fraction(observed: pd.DataFrame, expected: pd.DataFrame) -> pd.DataFrame:
    """Share of the observed change that spill-over predicts (0 when the signs differ, capped at 1)."""
    with np.errstate(divide="ignore", invalid="ignore"):
        f = expected / observed
    return f.clip(lower=0, upper=1).fillna(0)
