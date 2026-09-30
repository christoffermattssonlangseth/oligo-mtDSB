"""Figures for the region x time-point x treatment comparison.

Every function takes the tidy tables produced by :mod:`region_time_de`
(``load_results`` / ``run_all``) and returns a matplotlib ``Figure``. Nothing
here touches the AnnData object, so the notebook can build all figures from the
saved CSVs on a laptop without the h5ad.

Colour rules (one job per palette):
* magnitude (hit counts, activation index) -> single blue ramp
* polarity (log2 fold change, delta) -> blue / grey / red diverging
* identity (significance category on scatter plots) -> fixed categorical order
"""

from __future__ import annotations

import textwrap
from typing import Sequence

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

try:
    from .region_time_de import TISSUE_WIDE, Keys, gene_stratum_matrix, significant
except ImportError:  # utils/ itself is on sys.path
    from region_time_de import TISSUE_WIDE, Keys, gene_stratum_matrix, significant

# --- palette -----------------------------------------------------------------
SEQ_BLUE = ["#fcfcfb", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
CMAP_SEQ = LinearSegmentedColormap.from_list("seq_blue", SEQ_BLUE)
CMAP_DIV = LinearSegmentedColormap.from_list("div_blue_red", ["#0d366b", "#3987e5", "#f0efec", "#e34948", "#7a1414"])

CAT = {  # fixed order, never cycled
    "both": "#4a3aa7",
    "early only": "#eb6834",
    "late only": "#2a78d6",
    "neither": "#c9c8c2",
}
TEXT = "#0b0b0b"
TEXT_2 = "#52514e"
GRID = "#e6e5e1"


def _style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=TEXT_2, labelsize=8)
    ax.title.set_color(TEXT)


def _region_order(regions: Sequence[str]) -> list[str]:
    regs = sorted(set(regions) - {TISSUE_WIDE})
    if TISSUE_WIDE in set(regions):
        regs.append(TISSUE_WIDE)
    return regs


# --- 1. overview: hits per region x age, one panel per cell type --------------


def plot_hits_overview(
    summary: pd.DataFrame,
    value: str = "n_sig",
    cell_types: Sequence[str] | None = None,
    keys: Keys = Keys(),
    ncols: int = 4,
    annotate: bool = True,
    vmax: float | None = None,
) -> plt.Figure:
    """Small multiples: regions (rows) x age (columns) heatmap per cell type.

    ``value`` is one of ``n_sig, n_up, n_down, activation_index``. The colour
    scale is capped at ``vmax`` (default: the 95th percentile of the non-zero
    values) so that one extreme stratum does not wash out the rest; the
    annotated numbers are always the true values.
    """
    cell_types = list(cell_types or sorted(summary["cell_type"].unique()))
    regions = _region_order(summary["region"])
    ages = [str(a) for a in keys.ages]
    if vmax is None:
        nonzero = summary.loc[summary[value] > 0, value]
        vmax = float(np.percentile(nonzero, 95)) if len(nonzero) else 1.0
    vmax = max(vmax, 1.0)
    n = len(cell_types)
    ncols = min(ncols, n)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(2.1 * ncols + 1.5, 0.28 * len(regions) * nrows + 1.6), squeeze=False, layout="constrained"
    )
    for ax in axes.ravel()[n:]:
        ax.axis("off")
    for ax, ct in zip(axes.ravel(), cell_types):
        sub = summary[summary["cell_type"] == ct].pivot_table(index="region", columns="age", values=value, aggfunc="first")
        sub = sub.reindex(index=regions, columns=ages)
        im = ax.imshow(sub.to_numpy(dtype=float), cmap=CMAP_SEQ, vmin=0, vmax=vmax, aspect="auto")
        ax.set_xticks(range(len(ages)))
        ax.set_xticklabels([f"P{a}" for a in ages])
        ax.set_yticks(range(len(regions)))
        ax.set_yticklabels(regions if ax is axes[np.unravel_index(list(axes.ravel()).index(ax), axes.shape)[0], 0] else [])
        ax.set_title("\n".join(textwrap.wrap(ct, 20)), fontsize=9, loc="left")
        ax.tick_params(length=0)
        for s in ax.spines.values():
            s.set_visible(False)
        if annotate:
            arr = sub.to_numpy(dtype=float)
            for i in range(arr.shape[0]):
                for j in range(arr.shape[1]):
                    v = arr[i, j]
                    if np.isnan(v):
                        ax.text(j, i, "–", ha="center", va="center", fontsize=7, color=TEXT_2)
                    else:
                        txt = f"{v:.0f}" if value.startswith("n_") else f"{v:.1f}"
                        ax.text(j, i, txt, ha="center", va="center", fontsize=7, color="white" if v > 0.6 * vmax else TEXT)
    cbar = fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.5, pad=0.02)
    label = {"n_sig": "significant genes", "n_up": "up-regulated genes", "n_down": "down-regulated genes", "activation_index": "activation index"}.get(value, value)
    if float(summary[value].max()) > vmax:
        label += f" (colour capped at {vmax:.0f})"
    cbar.set_label(label, color=TEXT_2)
    cbar.outline.set_visible(False)
    fig.suptitle("Treatment effect (mtDSB vs control) by region and time point", ha="left", x=0.02, fontsize=11, color=TEXT)
    return fig


# --- 2. LFC at P21 vs LFC at P60 for one region ------------------------------


def _sig_category(g: pd.DataFrame, keys: Keys, padj_thr: float) -> pd.Series:
    e = g[f"padj_age{keys.early}"] < padj_thr
    l = g[f"padj_age{keys.late}"] < padj_thr
    return pd.Series(np.select([e & l, e & ~l, ~e & l], ["both", "early only", "late only"], "neither"), index=g.index)


def plot_age_scatter(
    age_contrast: pd.DataFrame,
    cell_type: str,
    region: str,
    keys: Keys = Keys(),
    padj_thr: float = 0.05,
    label_top: int = 12,
    forced_genes: Sequence[str] = (),
    ax: plt.Axes | None = None,
) -> plt.Figure:
    """Scatter of per-age treatment log2FC; colour = where the gene is significant.

    Genes with a significant late-minus-early difference (interaction z-test)
    get a dark ring. Labels: ``label_top`` genes by smallest interaction
    p-value plus ``forced_genes``.
    """
    g = age_contrast[(age_contrast["cell_type"] == cell_type) & (age_contrast["region"] == region)].copy()
    if ax is None:
        fig, ax = plt.subplots(figsize=(4.6, 4.4))
    else:
        fig = ax.figure
    if g.empty:
        ax.text(0.5, 0.5, "no data", ha="center", va="center", transform=ax.transAxes, color=TEXT_2)
        ax.set_title(f"{cell_type} - {region}", fontsize=9)
        return fig
    x, y = f"log2FC_age{keys.early}", f"log2FC_age{keys.late}"
    g["cat"] = _sig_category(g, keys, padj_thr)
    lim = float(np.nanmax(np.abs(g[[x, y]].to_numpy()))) * 1.05 or 1.0
    ax.axhline(0, color=GRID, lw=1, zorder=0)
    ax.axvline(0, color=GRID, lw=1, zorder=0)
    ax.plot([-lim, lim], [-lim, lim], color=GRID, lw=1, ls="--", zorder=0)
    for cat in ["neither", "early only", "late only", "both"]:
        s = g[g["cat"] == cat]
        ax.scatter(s[x], s[y], s=10 if cat == "neither" else 22, color=CAT[cat], alpha=0.5 if cat == "neither" else 0.9,
                   lw=0, label=f"{cat} (n={len(s)})", zorder=2, rasterized=cat == "neither")
    inter = g[g["padj"] < padj_thr]
    if len(inter):
        ax.scatter(inter[x], inter[y], s=60, facecolor="none", edgecolor=TEXT, lw=0.9, zorder=3,
                   label=f"P{keys.late} ≠ P{keys.early} (n={len(inter)})")
    to_label = pd.concat([g.sort_values("pvalue").head(label_top), g[g["gene"].isin(forced_genes)]]).drop_duplicates("gene")
    for _, r in to_label.iterrows():
        ax.annotate(r["gene"], (r[x], r[y]), xytext=(3, 3), textcoords="offset points", fontsize=7, color=TEXT)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel(f"log2FC mtDSB vs control, P{keys.early}", color=TEXT_2, fontsize=8)
    ax.set_ylabel(f"log2FC mtDSB vs control, P{keys.late}", color=TEXT_2, fontsize=8)
    ax.set_title(f"{cell_type} - {region}", fontsize=9, loc="left")
    ax.legend(frameon=False, fontsize=7, loc="upper left")
    _style(ax)
    return fig


def plot_age_scatter_grid(age_contrast: pd.DataFrame, cell_type: str, regions: Sequence[str] | None = None, ncols: int = 4, **kw) -> plt.Figure:
    """One :func:`plot_age_scatter` per region for a cell type."""
    regions = list(regions or _region_order(age_contrast.loc[age_contrast["cell_type"] == cell_type, "region"]))
    nrows = int(np.ceil(len(regions) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.9 * ncols, 3.8 * nrows), squeeze=False)
    for ax in axes.ravel()[len(regions):]:
        ax.axis("off")
    for ax, reg in zip(axes.ravel(), regions):
        plot_age_scatter(age_contrast, cell_type, reg, ax=ax, **kw)
        ax.get_legend().set_visible(ax is axes[0, 0])
    fig.tight_layout()
    return fig


# --- 3. gene x (region @ age) heatmap ------------------------------------------


def plot_gene_heatmap(
    results: pd.DataFrame,
    cell_type: str,
    genes: Sequence[str] | None = None,
    mask_padj: float | None = 0.05,
    top_n: int = 15,
    regions: Sequence[str] | None = None,
    keys: Keys = Keys(),
    vmax: float | None = None,
    padj_thr: float = 0.05,
    lfc_thr: float = 0.5,
    min_base_mean: float = 20.0,
) -> plt.Figure:
    """log2FC heatmap; cells that are not significant are blank.

    ``padj_thr / lfc_thr / min_base_mean`` define "significant" when the gene
    list is chosen automatically (``genes=None``).
    """
    mat = gene_stratum_matrix(
        results, cell_type, genes, "log2FC", mask_padj, top_n, padj_thr, lfc_thr, min_base_mean, regions=regions, keys=keys
    )
    if mat.empty:
        fig, ax = plt.subplots(figsize=(4, 2))
        ax.text(0.5, 0.5, f"{cell_type}: no significant genes", ha="center", va="center", color=TEXT_2)
        ax.axis("off")
        return fig
    vmax = vmax or float(np.nanmax(np.abs(mat.to_numpy()))) or 1.0
    fig, ax = plt.subplots(figsize=(0.42 * mat.shape[1] + 2.5, 0.22 * mat.shape[0] + 1.5))
    cmap = CMAP_DIV.copy()
    cmap.set_bad("#ffffff")
    im = ax.imshow(np.ma.masked_invalid(mat.to_numpy(dtype=float)), cmap=cmap, norm=TwoSlopeNorm(0, -vmax, vmax), aspect="auto")
    ax.set_xticks(range(mat.shape[1]))
    ax.set_xticklabels([c.replace(" @ ", "\nP") for c in mat.columns], fontsize=7, rotation=90)
    ax.set_yticks(range(mat.shape[0]))
    ax.set_yticklabels(mat.index, fontsize=7, style="italic")
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_xticks(np.arange(-0.5, mat.shape[1]), minor=True)
    ax.set_yticks(np.arange(-0.5, mat.shape[0]), minor=True)
    ax.grid(which="minor", color="#ffffff", lw=2)
    cbar = fig.colorbar(im, ax=ax, shrink=0.4, pad=0.02)
    cbar.set_label("log2FC mtDSB vs control", color=TEXT_2, fontsize=8)
    cbar.outline.set_visible(False)
    ax.set_title(f"{cell_type}" + (f"  (blank: padj ≥ {mask_padj})" if mask_padj else ""), fontsize=9, loc="left")
    return fig


# --- 4. interaction / region-specific hits --------------------------------------


def plot_interaction_lollipop(
    interaction: pd.DataFrame,
    cell_type: str,
    region: str,
    keys: Keys = Keys(),
    top_n: int = 20,
    padj_thr: float = 0.05,
    ax: plt.Axes | None = None,
) -> plt.Figure:
    """Top genes by age x condition interaction for one (cell type, region).

    Each gene shows the treatment effect at both time points; the dot colour
    follows the fixed early/late palette and a line joins them. Genes are
    ordered by the interaction estimate.
    """
    sub = interaction[(interaction["cell_type"] == cell_type) & (interaction["region"] == region)]
    inter = sub[sub["term"] == "age_x_condition"].sort_values("pvalue").head(top_n)
    if ax is None:
        fig, ax = plt.subplots(figsize=(4.5, 0.28 * max(len(inter), 4) + 1.2))
    else:
        fig = ax.figure
    if inter.empty:
        ax.text(0.5, 0.5, "no model (insufficient sections)", ha="center", va="center", transform=ax.transAxes, color=TEXT_2)
        ax.axis("off")
        return fig
    e = sub[sub["term"] == f"condition_at_{keys.early}"].set_index("gene")["log2FC"]
    l = sub[sub["term"] == f"condition_at_{keys.late}"].set_index("gene")["log2FC"]
    inter = inter.assign(lfc_e=inter["gene"].map(e), lfc_l=inter["gene"].map(l)).sort_values("log2FC")
    y = np.arange(len(inter))
    ax.axvline(0, color=GRID, lw=1, zorder=0)
    ax.hlines(y, inter["lfc_e"], inter["lfc_l"], color=GRID, lw=2, zorder=1)
    ax.scatter(inter["lfc_e"], y, color=CAT["early only"], s=28, zorder=2, label=f"P{keys.early}")
    ax.scatter(inter["lfc_l"], y, color=CAT["late only"], s=28, zorder=2, label=f"P{keys.late}")
    sig = inter["padj"] < padj_thr
    ax.set_yticks(y)
    ax.set_yticklabels([f"{g}{' *' if s else ''}" for g, s in zip(inter["gene"], sig)], fontsize=7, style="italic")
    ax.set_xlabel("log2FC mtDSB vs control", fontsize=8, color=TEXT_2)
    ax.set_title(f"{cell_type} - {region}\nage x condition interaction, * padj < {padj_thr}", fontsize=9, loc="left")
    ax.legend(frameon=False, fontsize=7, loc="lower right")
    _style(ax)
    return fig


def plot_region_specificity(
    region_vs_rest: pd.DataFrame,
    cell_type: str,
    age: str,
    padj_thr: float = 0.05,
    label_top: int = 10,
    ax: plt.Axes | None = None,
) -> plt.Figure:
    """Region log2FC vs rest-of-tissue log2FC, all regions in one panel.

    Points far from the diagonal are region-specific responses; those with a
    significant region-minus-rest difference are drawn opaque with a ring.
    """
    g = region_vs_rest[(region_vs_rest["cell_type"] == cell_type) & (region_vs_rest["age"].astype(str) == str(age))].copy()
    if ax is None:
        fig, ax = plt.subplots(figsize=(4.8, 4.6))
    else:
        fig = ax.figure
    if g.empty:
        ax.text(0.5, 0.5, "no data", ha="center", va="center", transform=ax.transAxes, color=TEXT_2)
        return fig
    lim = float(np.nanmax(np.abs(g[["log2FC_rest", "log2FC_region"]].to_numpy()))) * 1.05 or 1.0
    ax.plot([-lim, lim], [-lim, lim], color=GRID, lw=1, ls="--", zorder=0)
    ax.axhline(0, color=GRID, lw=1, zorder=0)
    ax.axvline(0, color=GRID, lw=1, zorder=0)
    ns = g[g["padj"] >= padj_thr]
    ax.scatter(ns["log2FC_rest"], ns["log2FC_region"], s=8, color=CAT["neither"], alpha=0.4, lw=0, rasterized=True, zorder=1)
    sig = g[g["padj"] < padj_thr]
    ax.scatter(sig["log2FC_rest"], sig["log2FC_region"], s=26, color=CAT["both"], edgecolor=TEXT, lw=0.6, zorder=2)
    for _, r in sig.sort_values("pvalue").head(label_top).iterrows():
        ax.annotate(f"{r['gene']} ({r['region']})", (r["log2FC_rest"], r["log2FC_region"]), xytext=(3, 3),
                    textcoords="offset points", fontsize=7, color=TEXT)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel("log2FC, rest of tissue", fontsize=8, color=TEXT_2)
    ax.set_ylabel("log2FC, region", fontsize=8, color=TEXT_2)
    ax.set_title(f"{cell_type}, P{age}: region-specific treatment effects (n sig = {len(sig)})", fontsize=9, loc="left")
    _style(ax)
    return fig


# --- 5. volcano --------------------------------------------------------------


def plot_volcano(
    res: pd.DataFrame,
    title: str = "",
    padj_thr: float = 0.05,
    lfc_thr: float = 0.5,
    min_base_mean: float = 20.0,
    label_top: int = 15,
    forced_genes: Sequence[str] = (),
    ax: plt.Axes | None = None,
) -> plt.Figure:
    """Standard volcano on one result table (one stratum or one term)."""
    g = res.dropna(subset=["padj"]).copy()
    g["sig"] = significant(g, padj_thr, lfc_thr, min_base_mean)
    g["nlp"] = -np.log10(g["padj"].clip(lower=1e-300))
    if ax is None:
        fig, ax = plt.subplots(figsize=(4.6, 4.2))
    else:
        fig = ax.figure
    ns = g[~g["sig"]]
    ax.scatter(ns["log2FC"], ns["nlp"], s=8, color=CAT["neither"], alpha=0.5, lw=0, rasterized=True)
    up = g[g["sig"] & (g["log2FC"] > 0)]
    dn = g[g["sig"] & (g["log2FC"] < 0)]
    ax.scatter(up["log2FC"], up["nlp"], s=16, color="#e34948", lw=0, label=f"up (n={len(up)})")
    ax.scatter(dn["log2FC"], dn["nlp"], s=16, color="#2a78d6", lw=0, label=f"down (n={len(dn)})")
    ax.axhline(-np.log10(padj_thr), color=GRID, ls="--", lw=1)
    for x in (-lfc_thr, lfc_thr):
        ax.axvline(x, color=GRID, ls="--", lw=1)
    lab = pd.concat([g[g["sig"]].sort_values("padj").head(label_top), g[g["gene"].isin(forced_genes)]]).drop_duplicates("gene")
    for _, r in lab.iterrows():
        ax.annotate(r["gene"], (r["log2FC"], r["nlp"]), xytext=(3, 2), textcoords="offset points", fontsize=7, color=TEXT)
    ax.set_xlabel("log2FC mtDSB vs control", fontsize=8, color=TEXT_2)
    ax.set_ylabel("-log10 padj", fontsize=8, color=TEXT_2)
    ax.set_title(title, fontsize=9, loc="left")
    ax.legend(frameon=False, fontsize=7)
    _style(ax)
    return fig


def save(fig: plt.Figure, path, formats: Sequence[str] = ("png", "pdf"), dpi: int = 200) -> None:
    """Save one figure in several formats next to each other."""
    from pathlib import Path

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    for ext in formats:
        fig.savefig(path.with_suffix(f".{ext}"), dpi=dpi, bbox_inches="tight", facecolor="white")
