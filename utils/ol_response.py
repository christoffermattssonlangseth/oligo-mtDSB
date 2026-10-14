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


# ----------------------------------------------------------------------------
# Which genes drive the classifier?
# ----------------------------------------------------------------------------


def summarize_coefficients(coefs: pd.DataFrame) -> pd.DataFrame:
    """Per (task, gene): fold-mean weight, SD, sign consistency and mean rank.

    ``coefs`` is the long table from ``classify_cells(save_coefficients=True)``.
    ``sign_consistency`` = fraction of folds whose weight has the sign of the
    mean; ``mean_rank`` ranks genes by weight within each fold (1 = most
    positive, i.e. most mtDSB-like) and averages over folds.
    """
    c = coefs.copy()
    c["rank"] = c.groupby(["task", "fold"])["coef"].rank(ascending=False)
    g = c.groupby(["task", "gene"])
    out = g.agg(mean_coef=("coef", "mean"), sd_coef=("coef", "std"), n_folds=("coef", "size"), mean_rank=("rank", "mean"))
    c["same_sign"] = np.sign(c["coef"]) == np.sign(c.groupby(["task", "gene"])["coef"].transform("mean"))
    out["sign_consistency"] = c.groupby(["task", "gene"])["same_sign"].mean()
    return out.reset_index()


def gene_score_association(X, sample_of: np.ndarray, score: np.ndarray, genes: Sequence[str],
                           min_cells: int = 200, eps: float = 1e-4) -> pd.DataFrame:
    """Within-section Pearson correlation of every gene with a per-cell score.

    ``X`` (cells x genes, sparse or dense) holds log-normalized expression for
    the cells of one type, ``score`` their out-of-fold p(mtDSB); the score is
    logit-transformed. Computing the correlation *within* each section removes
    section-level (batch, plane) differences: a driver gene must separate
    high- from low-scoring cells inside the same section.

    Returns long rows ``sample_id, gene, r, mean_expr, frac_expr, n_cells``.
    """
    from scipy import sparse as sp

    genes = np.asarray(genes)
    y_all = np.log(np.clip(score, eps, 1 - eps) / (1 - np.clip(score, eps, 1 - eps)))
    rows = []
    for sid in pd.unique(sample_of):
        idx = np.flatnonzero(sample_of == sid)
        if len(idx) < min_cells:
            continue
        Xs = X[idx]
        Xs = Xs.toarray() if sp.issparse(Xs) else np.asarray(Xs)
        Xs = Xs.astype(np.float64)
        y = y_all[idx] - y_all[idx].mean()
        mu = Xs.mean(0)
        Xc = Xs - mu
        sx = np.sqrt((Xc**2).sum(0))
        sy = np.sqrt((y**2).sum())
        with np.errstate(invalid="ignore", divide="ignore"):
            r = (Xc.T @ y) / (sx * sy)
        rows.append(pd.DataFrame({"sample_id": sid, "gene": genes, "r": r, "mean_expr": mu,
                                  "frac_expr": (Xs > 0).mean(0), "n_cells": len(idx)}))
    return pd.concat(rows, ignore_index=True)


def summarize_association(assoc: pd.DataFrame, sections: pd.DataFrame) -> pd.DataFrame:
    """Mean within-section r per gene and group (responder / non-responder / control).

    Also ``min_r_responders`` (the weakest responder), so a gene driven by one
    animal does not rank high.
    """
    a = assoc.merge(sections[["group"]], left_on="sample_id", right_index=True)
    wide = a.pivot_table(index="gene", columns="group", values="r", aggfunc="mean")
    wide.columns = [f"r_{c}" for c in wide.columns]
    resp = a.loc[a["group"] == "responder"].pivot(index="gene", columns="sample_id", values="r")
    wide["min_r_responders"] = resp.min(1)
    wide["r_resp_minus_ctrl"] = wide.get("r_responder") - wide.get("r_control")
    return wide.reset_index()
