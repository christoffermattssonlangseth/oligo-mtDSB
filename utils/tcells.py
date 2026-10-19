"""T cells in the Xenium sections: marker co-expression calls and a co-expression null.

T cells are far too rare (~2-4 per 10k cells) to form a cluster, so they are called per cell from the
co-expression of core T-cell genes. ``Trbc1`` and ``Ighm`` are excluded: both are detected in tens of
thousands of cortical neurons (off-target signal), and ``Trac`` / ``Nkg7`` are not in the panel.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

CORE = ["Cd3d", "Cd3e", "Cd3g", "Cd247", "Cd2", "Cd8a", "Cd8b1", "Itk", "Cxcr6"]
"""Core T-cell genes counted for the call (TCR complex, CD8, T-cell kinases, tissue-resident Cxcr6)."""
EXTRA = ["Ptprc", "Lck", "Lat", "Zap70", "Cd6", "Cd5", "Cd28", "Il7r", "Cd4", "Gzma", "Gzmb", "Gzmk", "Prf1",
         "Ifng", "Klrk1", "Klrd1", "Pdcd1", "Ctla4", "Icos", "Foxp3", "Cd69", "Ccr7", "Sell"]
"""Supporting T / NK / activation genes, reported but not used for the call."""
NOISY = ["Trbc1", "Ighm"]
"""Lymphocyte probes with off-target signal in neurons; shown in the panel QC only."""
PANEL = CORE + EXTRA + NOISY + ["Cd19", "Ms4a1", "Cd79a", "Cd79b", "Jchain", "Mzb1"]

BORDER_REGIONS = {"Meningeal border/glia limitans", "Ventricular system", "Vasculature"}


def call_tcells(counts: pd.DataFrame, min_tx: int = 3, min_genes: int = 2, genes=CORE) -> pd.Series:
    """True for cells with >= ``min_tx`` transcripts from ``genes`` spread over >= ``min_genes`` genes."""
    m = counts[list(genes)].to_numpy()
    return pd.Series(((m.sum(1) >= min_tx) & ((m > 0).sum(1) >= min_genes)), index=counts.index)


def shuffle_null(counts: pd.DataFrame, groups: pd.Series, n_rep: int = 20, seed: int = 0, **call_kw) -> np.ndarray:
    """Number of calls after permuting each gene independently across cells within each group.

    Keeps every gene's per-section total and every cell's identity but breaks co-expression, so it gives
    the number of calls expected from scattered (ambient / off-target) transcripts alone. ``counts`` must
    include the zero rows (all cells of each group), otherwise the null is inflated.
    """
    rng = np.random.default_rng(seed)
    genes = call_kw.pop("genes", CORE)
    m = counts[list(genes)].to_numpy()
    blocks = [np.flatnonzero(groups.to_numpy() == g) for g in pd.unique(groups)]
    out = np.empty(n_rep, int)
    for r in range(n_rep):
        p = m.copy()
        for ii in blocks:
            for j in range(p.shape[1]):
                p[ii, j] = m[rng.permutation(ii), j]
        out[r] = call_tcells(pd.DataFrame(p, columns=genes), genes=genes, **call_kw).sum()
    return out
