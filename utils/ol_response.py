"""Per-animal oligodendrocyte (OL) response to mtDSB: who responds, how strongly, and where.

Three readouts of the same question, deliberately built from different inputs
so that agreement between them is informative:

1. **ISR score** (``isr_scores``): share of OL transcripts from a fixed set of
   canonical ATF4 / mitochondrial integrated-stress-response targets
   (``ISR_GENES``), per (section, region), from the pseudobulk. The gene set
   was fixed from the literature *before* looking at the per-section results
   and uses neither the clustering nor any model.
2. **Classifier score** (``high_score_fraction``): fraction of OLs whose
   out-of-fold p(mtDSB) from ``sample_classifier.classify_cells`` exceeds a
   threshold, per (section, region). Uses the condition labels of the *other*
   sections, never the cluster labels (DA-Oligodendrocytes are merged into
   Oligodendrocytes before fitting).
3. **Cluster share** (``cluster_share``): fraction of OLs annotated as
   ``DA-Oligodendrocytes`` per (section, region). Descriptive only; depends on
   the clustering.

``responder_calls`` combines the section-level ISR score with the control
range: a mtDSB section is a responder when its whole-section ISR score exceeds
the highest control by more than ``margin`` (relative). This is a description
of 12 animals, not a test; it is reported together with every per-gene value.
"""

from __future__ import annotations

from typing import Sequence

import numpy as np
import pandas as pd

from utils.region_time_de import TISSUE_WIDE, Pseudobulk, merge_cell_types

ISR_GENES = ["Atf4", "Atf5", "Ddit3", "Trib3", "Gdf15", "Sesn2", "Fgf21", "Psat1", "Cth", "Eif4ebp1"]
"""Canonical ATF4 / mitochondrial ISR targets present on the 5k panel (fixed a priori).

Left out although on the panel: Mt1/Mt2 (metal response), Vegfa (hypoxia),
Hspd1 (chaperone; mitochondrial UPR rather than ISR).
"""

OL = "Oligodendrocytes"
DA_OL = "DA-Oligodendrocytes"
GENERAL_LABELS = {DA_OL: OL}
"""Fold the disease-defined cluster into its general class."""


def general_pseudobulk(pb: Pseudobulk) -> Pseudobulk:
    return merge_cell_types(pb, GENERAL_LABELS)


def isr_scores(pb: Pseudobulk, genes: Sequence[str] = ISR_GENES, cell_type: str = OL) -> pd.DataFrame:
    """ISR set CPM and per-gene CPM per (section, region) for one cell type.

    ``pb`` should already have general labels (``general_pseudobulk``).
    Returns long-format rows with ``sample_id, region, n_cells, isr_cpm`` and
    one column per gene (CPM).
    """
    missing = [g for g in genes if g not in pb.counts.columns]
    if missing:
        raise KeyError(f"ISR genes not in panel: {missing}")
    sel = (pb.meta["cell_type"] == cell_type).to_numpy()
    c = pb.counts.loc[sel]
    lib = c.sum(1)
    out = pb.meta.loc[sel, ["sample_id", "region", "n_cells"]].copy()
    out["isr_cpm"] = c[list(genes)].sum(1) / lib * 1e6
    for g in genes:
        out[g] = c[g] / lib * 1e6
    return out.reset_index(drop=True)


def responder_calls(isr: pd.DataFrame, sections: pd.DataFrame, margin: float = 0.2, control: str = "control") -> pd.DataFrame:
    """Classify each section as control / responder / non-responder from the whole-section ISR score.

    ``sections`` is indexed by ``sample_id`` with at least ``condition``.
    Adds ``isr_cpm``, ``isr_vs_control_max`` (ratio to the highest control) and ``group``.
    """
    tw = isr.loc[isr["region"] == TISSUE_WIDE].set_index("sample_id")["isr_cpm"]
    s = sections.copy()
    s["isr_cpm"] = tw.reindex(s.index)
    ctrl_max = s.loc[s["condition"] == control, "isr_cpm"].max()
    s["isr_vs_control_max"] = s["isr_cpm"] / ctrl_max
    s["group"] = np.where(
        s["condition"] == control, "control",
        np.where(s["isr_vs_control_max"] > 1 + margin, "responder", "non-responder"),
    )
    return s


def high_score_fraction(cell_scores: pd.DataFrame, task: str = "condition", threshold: float = 0.8,
                        cell_type: str = OL, min_cells: int = 50) -> pd.DataFrame:
    """Fraction of a cell type's cells with out-of-fold p(mtDSB) > ``threshold``, per (section, region).

    Includes region ``All regions``. Groups with fewer than ``min_cells`` cells are dropped.
    """
    col = f"p_pos[{task}]"
    cs = cell_scores.loc[(cell_scores["cell_type"] == cell_type) & cell_scores[col].notna(), ["sample_id", "region", col]]
    cs = pd.concat([cs, cs.assign(region=TISSUE_WIDE)])
    g = cs.assign(high=cs[col] > threshold).groupby(["sample_id", "region"])
    out = g.agg(n_cells=("high", "size"), frac_high=("high", "mean"), mean_p=(col, "mean")).reset_index()
    return out.loc[out["n_cells"] >= min_cells].reset_index(drop=True)


def cluster_share(pb_raw: Pseudobulk, cluster: str = DA_OL, parent: str = OL, min_cells: int = 50) -> pd.DataFrame:
    """Share of ``parent`` + ``cluster`` cells annotated as ``cluster``, per (section, region).

    Uses kept and dropped pseudobulk rows so small counts are not lost. Takes the
    pseudobulk with the *original* labels.
    """
    meta = pd.concat([pb_raw.meta, pb_raw.dropped]) if len(pb_raw.dropped) else pb_raw.meta
    meta = meta.assign(n_cells=meta["n_cells"].astype(int))
    meta = meta.loc[meta["cell_type"].isin([cluster, parent])]
    tab = meta.pivot_table(index=["sample_id", "region"], columns="cell_type", values="n_cells", aggfunc="sum", fill_value=0)
    for c in (cluster, parent):
        if c not in tab:
            tab[c] = 0
    tot = tab[cluster] + tab[parent]
    out = pd.DataFrame({"n_cells": tot, "n_cluster": tab[cluster], "share": tab[cluster] / tot.where(tot > 0)}).reset_index()
    return out.loc[out["n_cells"] >= min_cells].reset_index(drop=True)


def region_by_group(values: pd.DataFrame, sections: pd.DataFrame, value: str) -> pd.DataFrame:
    """Mean of per-section ``value`` per region and group, plus the responder / control ratio."""
    v = values.merge(sections[["group"]], left_on="sample_id", right_index=True)
    t = v.pivot_table(index="region", columns="group", values=value, aggfunc="mean")
    if {"responder", "control"} <= set(t.columns):
        t["responder / control"] = t["responder"] / t["control"]
    return t


def region_rank_per_animal(values: pd.DataFrame, sections: pd.DataFrame, value: str, group: str = "responder") -> pd.DataFrame:
    """Rank of each region within each animal of ``group`` (1 = highest ``value``), excluding ``All regions``."""
    ids = sections.index[sections["group"] == group]
    v = values.loc[values["sample_id"].isin(ids) & (values["region"] != TISSUE_WIDE)]
    w = v.pivot(index="region", columns="sample_id", values=value)
    r = w.rank(ascending=False)
    r["mean rank"] = r.mean(1)
    return r.sort_values("mean rank")
