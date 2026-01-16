"""
Auto-generated utilities for spatial_utils.
Do not edit by hand without moving changes back into notebooks.

Each function below was extracted from exported analysis notebooks.
"""

from typing import *
import os

"""
Spatial plotting + ordering utilities for RRMap / EAE projects.
"""

import warnings
warnings.filterwarnings("ignore")

from typing import Iterable, Optional, Union

import numpy as np
import pandas as pd
import scanpy as sc

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib import colors as mcolors

from pandas.api.types import is_numeric_dtype, is_categorical_dtype
from scipy.sparse import issparse
from anndata import AnnData


# ------------------------------------------------------------------
# 1) Global ordering: course + region
# ------------------------------------------------------------------

COURSE_ORDER = [
    "MOG CFA",
    "early onset",
    "chronic peak",
    "chronic long",

    "PLP CFA",
    "non symptomatic",
    "monophasic",

    "onset I",
    "onset II",
    "peak I",
    "remitt I",
    "peak II",
    "remitt II",
    "peak III",
]

REGION_ORDER = ["L", "T", "C"]


def apply_course_region_and_sample_order(
    ad: AnnData,
    sample_key: str = "sample_id",
    course_order: Optional[Iterable[str]] = None,
    region_order: Optional[Iterable[str]] = None,
) -> AnnData:
    """
    Make ad.obs['course'], ad.obs['region'], and ad.obs[sample_key]
    ordered categoricals based on (course, region, sample_id).

    Stores the ordered sample IDs in:
        ad.uns[f"{sample_key}_order_by_course_region"].

    Parameters
    ----------
    ad
        AnnData object with obs columns 'course', 'region', and `sample_key`.
    sample_key
        Name of the sample column in ad.obs.
    course_order
        Optional custom ordering for 'course'. If None, uses COURSE_ORDER.
    region_order
        Optional custom ordering for 'region'. If None, uses REGION_ORDER.

    Returns
    -------
    AnnData
        The same object, modified in-place and also returned for convenience.
    """
    if course_order is None:
        course_order = COURSE_ORDER
    if region_order is None:
        region_order = REGION_ORDER

    # --- course ---
    if "course" in ad.obs.columns:
        ad.obs["course"] = ad.obs["course"].astype(
            pd.CategoricalDtype(categories=list(course_order), ordered=True)
        )

    # --- region ---
    if "region" in ad.obs.columns:
        ad.obs["region"] = ad.obs["region"].astype(
            pd.CategoricalDtype(categories=list(region_order), ordered=True)
        )

    # --- sample_id ordered by (course, region, sample_id) ---
    if (
        sample_key in ad.obs.columns
        and "course" in ad.obs.columns
        and "region" in ad.obs.columns
    ):
        tmp = (
            ad.obs[[sample_key, "course", "region"]]
            .drop_duplicates()
            .dropna(subset=["course", "region"])
            .copy()
        )

        tmp["course"] = tmp["course"].astype(
            pd.CategoricalDtype(categories=list(course_order), ordered=True)
        )
        tmp["region"] = tmp["region"].astype(
            pd.CategoricalDtype(categories=list(region_order), ordered=True)
        )

        tmp = tmp.sort_values(["course", "region", sample_key])
        sample_order = tmp[sample_key].tolist()

        ad.obs[sample_key] = ad.obs[sample_key].astype(
            pd.CategoricalDtype(categories=sample_order, ordered=True)
        )

        ad.uns[f"{sample_key}_order_by_course_region"] = sample_order

    return ad


# ------------------------------------------------------------------
# 2) Compact spatial plot (obs or gene, grouped by sample)
# ------------------------------------------------------------------

def plot_spatial_compact_fast(
    ad: AnnData,
    color: str = "leiden_2",  # obs column *or* gene name
    groupby: str = "sample_id",
    spot_size: float = 8,
    cols: int = 3,
    height: float = 8,
    legend_col_width: float = 1.2,
    palette: Optional[Union[dict, list, str]] = None,
    rasterized: bool = True,
    invert_y: bool = True,
    dpi: int = 120,
    highlight: Optional[Union[str, Iterable[str]]] = None,
    group_order: Optional[Iterable[str]] = None,
    background: str = "white",
    grey_alpha: float = 0.2,     # alpha for non-highlighted categories
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    cmap_name: str = "viridis",
    shared_scale: bool = False,  # if True: vmin/vmax from whole `ad` (not per subset)
):
    """
    Compact faceted spatial plot, one panel per group (e.g. sample_id),
    for either an obs column or a gene.

    - Handles continuous (gene / numeric obs) and categorical obs.
    - Respects ordered categoricals in `groupby` if present.
    - Can share color scale across all panels (shared_scale=True).
    - Writes color settings into `ad.uns[f"{color}_continuous"]` (continuous)
      or `ad.uns[f"{color}_colors"]` (categorical).

    Parameters
    ----------
    ad
        AnnData with `ad.obsm['spatial']`.
    color
        Name of an obs column (categorical or continuous) or a gene.
    groupby
        obs column to facet by (e.g. 'sample_id').
    spot_size
        Size of scatter markers.
    cols
        Number of columns in the panel grid.
    height
        Overall figure height (in inches).
    legend_col_width
        Width allocated for the legend / colorbar column.
    palette
        For categorical obs: dict {category -> color}, list of colors,
        or name of a matplotlib colormap. If None, uses default Scanpy palette.
    rasterized
        Rasterize the scatter artists (good for large ST plots).
    invert_y
        If True, invert y-axis to match image coordinates.
    dpi
        Figure DPI.
    highlight
        For categorical obs: value or list of values to highlight.
        Non-highlighted categories are greyed out (alpha=grey_alpha).
    group_order
        Optional manual ordering of groups in `groupby`.
    background
        Figure and axes background color.
    grey_alpha
        Alpha for non-highlighted categories in categorical mode.
    vmin, vmax
        If provided, overrides automatic min/max for continuous color.
    cmap_name
        Name of colormap for continuous mode (or if `palette` is str).
    shared_scale
        If True, compute vmin/vmax globally across the full AnnData
        for `color`. If False, based on current values only.

    Returns
    -------
    None
        Shows the plot.
    """
    # pick background + text color
    fig_face = background
    ax_face = background
    text_color = "white" if background in ("black", "#000000", "k") else "black"

    # ----- 0) Preconditions -----
    if "spatial" not in ad.obsm:
        raise ValueError("ad.obsm['spatial'] not found")

    if groupby not in ad.obs.columns:
        raise KeyError(f"groupby {groupby!r} not in ad.obs")

    # decide where 'color' comes from: obs vs var (gene)
    if color in ad.obs.columns:
        color_source = "obs"
        col_series = ad.obs[color]
    elif color in ad.var_names:
        color_source = "var"   # gene expression
        col_series = None
    else:
        raise KeyError(
            f"{color!r} not found in ad.obs.columns or ad.var_names "
            "(expected an obs column or a gene name)."
        )

    coords = np.asarray(ad.obsm["spatial"])[:, :2]

    # Detect continuous vs categorical
    if color_source == "var":
        # genes: always continuous
        is_continuous = True
    else:
        if is_categorical_dtype(col_series):
            is_continuous = False
        else:
            is_continuous = is_numeric_dtype(col_series)

    # ----------------------------------------------------
    # 1) Build colors_arr differently for cont vs cat
    # ----------------------------------------------------
    if is_continuous:
        # ===== CONTINUOUS MODE =====
        if color_source == "obs":
            vals = col_series.to_numpy(dtype=float)
        else:
            # from var / gene expression
            gene_idx = ad.var_names.get_loc(color)
            x = ad.X[:, gene_idx]
            if issparse(x):
                vals = x.toarray().ravel()
            else:
                vals = np.asarray(x).ravel()

        # choose colormap
        if palette is None:
            cmap = plt.get_cmap(cmap_name)
        elif isinstance(palette, str):
            cmap = plt.get_cmap(palette)
        else:
            cmap = palette

        # ----- determine vmin/vmax -----
        if shared_scale:
            # global min/max across the *full AnnData* passed to this function
            full = ad
            if color_source == "var":
                gx = full.X[:, full.var_names.get_loc(color)]
                if issparse(gx):
                    full_vals = gx.toarray().ravel()
                else:
                    full_vals = np.asarray(gx).ravel()
            else:
                full_vals = pd.to_numeric(full.obs[color], errors="coerce").to_numpy()

            finite_full = np.isfinite(full_vals)
            if finite_full.sum() == 0:
                raise ValueError(f"All values for '{color}' are NaN or non-finite.")
            vmin_use = float(np.min(full_vals[finite_full]))
            vmax_use = float(np.max(full_vals[finite_full]))
        else:
            finite_mask = np.isfinite(vals)
            if finite_mask.sum() == 0:
                raise ValueError(f"All values for '{color}' are NaN or non-finite.")
            vmin_use = float(np.min(vals[finite_mask]))
            vmax_use = float(np.max(vals[finite_mask]))

        # user overrides everything
        if vmin is not None:
            vmin_use = float(vmin)
        if vmax is not None:
            vmax_use = float(vmax)

        # avoid zero-range
        if vmin_use == vmax_use:
            vmin_use -= 1.0
            vmax_use += 1.0

        norm = mcolors.Normalize(vmin=vmin_use, vmax=vmax_use)

        colors_arr = np.zeros((vals.size, 4), dtype=float)
        finite_mask = np.isfinite(vals)
        colors_arr[finite_mask] = cmap(norm(vals[finite_mask]))
        colors_arr[~finite_mask] = (0, 0, 0, 0)

        # store continuous settings in uns (also fine for genes)
        ad.uns[f"{color}_continuous"] = {
            "vmin": float(vmin_use),
            "vmax": float(vmax_use),
            "cmap": cmap.name if hasattr(cmap, "name") else str(cmap),
        }

        cat_names = None
        cat_codes = None

    else:
        # ===== CATEGORICAL (only for obs) =====
        # preserve existing categorical order if present
        if is_categorical_dtype(col_series):
            cats = col_series.cat.remove_unused_categories()
        else:
            cats = col_series.astype("category")

        cat_names = cats.cat.categories
        cat_codes = cats.cat.codes.to_numpy()

        # palette handling
        if isinstance(palette, dict):
            col_list = [palette[c] for c in cat_names]
        elif isinstance(palette, (list, tuple)):
            if len(palette) < len(cat_names):
                raise ValueError("Palette shorter than number of categories.")
            col_list = list(palette)[:len(cat_names)]
        elif f"{color}_colors" in ad.uns:
            col_list = list(ad.uns[f"{color}_colors"])
            if len(col_list) != len(cat_names):
                raise ValueError(f"{color}_colors length != categories.")
        else:
            base = list(
                sc.pl.palettes.default_64
                if hasattr(sc.pl.palettes, "default_64")
                else sc.pl.palettes.default_102
            )
            reps = int(np.ceil(len(cat_names) / len(base)))
            col_list = (base * reps)[:len(cat_names)]

        ad.uns[f"{color}_colors"] = col_list

        rgba = np.array([mcolors.to_rgba(c) for c in col_list], dtype=float)

        colors_arr = np.empty((cat_codes.size, 4), dtype=float)
        colors_arr[cat_codes >= 0] = rgba[cat_codes[cat_codes >= 0]]
        colors_arr[cat_codes < 0] = (0, 0, 0, 0)

        # ----- highlighting logic -----
        if highlight is not None:
            # allow single value or list/tuple/array
            if not isinstance(highlight, (list, tuple, set, np.ndarray)):
                highlight = [highlight]
            # convert to string for robust matching
            highlight_str = {str(h) for h in highlight}

            cat_name_str = np.array([str(c) for c in cat_names])
            keep_cat_mask = np.isin(cat_name_str, list(highlight_str))  # per-category

            # use user-defined alpha for greyed-out categories
            grey_rgba = (0.8, 0.8, 0.8, float(grey_alpha))

            valid = cat_codes >= 0
            keep_flag = np.zeros_like(cat_codes, dtype=bool)
            keep_flag[valid] = keep_cat_mask[cat_codes[valid]]

            # grey out all non-highlighted cells
            colors_arr[valid & ~keep_flag] = grey_rgba

            # also grey in legend
            col_list = [
                col_list[k] if keep_cat_mask[k] else mcolors.to_hex(grey_rgba)
                for k in range(len(cat_names))
            ]

    # ----------------------------------------------------
    # 2) Precompute group indices (RESPECT ORDERED CATEGORICAL)
    # ----------------------------------------------------
    gser = ad.obs[groupby]

    if group_order is not None:
        group_order = [str(g) for g in group_order]
        present = set(gser.astype(str))
        uniq_groups = [g for g in group_order if g in present]
    else:
        # if groupby is an ordered categorical, respect its category order
        if is_categorical_dtype(gser) and gser.cat.ordered:
            cats = list(gser.cat.categories)
            present = set(gser.astype(str))
            uniq_groups = [str(c) for c in cats if str(c) in present]
        else:
            # fallback: sorted unique strings
            uniq_groups = sorted(gser.astype(str).unique())

    gvals = gser.astype(str).to_numpy()
    gid_to_idx = {g: i for i, g in enumerate(uniq_groups)}
    gcodes = np.array([gid_to_idx.get(g, -1) for g in gvals], dtype=int)

    group_indices = [np.flatnonzero(gcodes == gi) for gi in range(len(uniq_groups))]

    # ----------------------------------------------------
    # 3) Figure layout
    # ----------------------------------------------------
    n = len(uniq_groups)
    rows = int(np.ceil(n / cols))
    panel_w = height * cols * 0.6 / rows
    fig_w = panel_w + legend_col_width

    plt.ioff()
    fig = plt.figure(figsize=(fig_w, height), dpi=dpi, constrained_layout=False)

    # background for figure
    fig.patch.set_facecolor(fig_face)

    gs = GridSpec(
        rows, cols + 1, figure=fig,
        width_ratios=[1]*cols + [legend_col_width / (fig_w - legend_col_width)],
        wspace=0.02, hspace=0.02
    )

    # ----------------------------------------------------
    # 4) Panels
    # ----------------------------------------------------
    for i, sid in enumerate(uniq_groups):
        r, c = divmod(i, cols)
        ax = fig.add_subplot(gs[r, c])

        # panel background
        ax.set_facecolor(ax_face)

        idx = group_indices[i]
        if idx.size:
            xy = coords[idx]
            ax.scatter(
                xy[:, 0], -xy[:, 1],
                c=colors_arr[idx],
                s=spot_size,
                marker='o',
                linewidths=0,
                rasterized=rasterized
            )

        meta_strings = []

        if "region" in ad.obs.columns:
            region_vals = (
                ad.obs.loc[ad.obs[groupby] == sid, "region"]
                .dropna().astype(str).unique()
            )
            if len(region_vals) == 0:
                meta_strings.append("Region: unknown")
            elif len(region_vals) == 1:
                meta_strings.append(f"Region: {region_vals[0]}")
            else:
                meta_strings.append("Region: mixed")

        if "course" in ad.obs.columns:
            course_vals = (
                ad.obs.loc[ad.obs[groupby] == sid, "course"]
                .dropna().astype(str).unique()
            )
            if len(course_vals) == 0:
                meta_strings.append("Course: unknown")
            elif len(course_vals) == 1:
                meta_strings.append(f"Course: {course_vals[0]}")
            else:
                meta_strings.append("Course: mixed")

        if meta_strings:
            title = f"{sid}\n[{ ' | '.join(meta_strings) }]"
        else:
            title = str(sid)

        ax.set_title(title, fontsize=5, pad=2, color=text_color)
        ax.set_aspect("equal")
        if invert_y:
            ax.invert_yaxis()
        ax.set_axis_off()

    # blank unused panels
    for j in range(n, rows * cols):
        r, c = divmod(j, cols)
        ax = fig.add_subplot(gs[r, c])
        ax.set_facecolor(ax_face)
        ax.axis("off")

    # ----------------------------------------------------
    # 5) Legend / Colorbar
    # ----------------------------------------------------
    ax_leg = fig.add_subplot(gs[:, -1])
    ax_leg.set_facecolor(ax_face)
    ax_leg.axis("off")

    if is_continuous:
        sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
        sm.set_array([])
        cbar = fig.colorbar(sm, cax=ax_leg)
        cbar.set_label(color, rotation=90, color=text_color)
        cbar.ax.yaxis.set_tick_params(color=text_color)
        plt.setp(plt.getp(cbar.ax.axes, "yticklabels"), color=text_color)
    else:
        handles = [
            Line2D(
                [0], [0], marker="o", color="w",
                markerfacecolor=col_list[k], markersize=7, label=str(cat)
            )
            for k, cat in enumerate(cat_names)
        ]
        leg = ax_leg.legend(
            handles=handles,
            title=color,
            frameon=False,
            loc="center left",
            labelcolor=text_color,
            title_fontsize=10
        )
        leg.get_title().set_color(text_color)
        for text in leg.get_texts():
            text.set_color(text_color)

    fig.subplots_adjust(
        left=0.01, right=0.98, top=0.98, bottom=0.02,
        wspace=0.02, hspace=0.02
    )

    plt.ion()
    plt.show()


def plot_spatial_compact(
    ad,
    color="leiden_2",
    groupby="sample_id",
    spot_size=18,
    cols=3,
    height=8,
    legend_col_width=1.2,
    palette=None,          # dict {category: "#hex"} or list in desired order
):
    # --- categories from the FULL object (freeze order) ---
    cats = ad.obs[color].astype("category").cat.categories

    # --- build a shared color list in that order ---
    if isinstance(palette, dict):
        cols_list = [palette[c] for c in cats]
    elif isinstance(palette, (list, tuple)):
        if len(palette) < len(cats):
            raise ValueError("Palette list shorter than number of categories.")
        cols_list = list(palette)[:len(cats)]
    elif f"{color}_colors" in ad.uns:
        cols_list = list(ad.uns[f"{color}_colors"])
        if len(cols_list) != len(cats):
            raise ValueError(f"{color}_colors length does not match categories.")
    else:
        base = sc.pl.palettes.default_64 if hasattr(sc.pl.palettes, "default_64") else sc.pl.palettes.default_102
        reps = int(np.ceil(len(cats) / len(base)))
        cols_list = (base * reps)[:len(cats)]

    # store for reuse
    ad.uns[f"{color}_colors"] = cols_list

    # --- layout ---
    sids = list(ad.obs[groupby].unique())
    n = len(sids)
    rows = int(np.ceil(n / cols))

    panel_w = height * cols * 0.6 / rows
    fig_w = panel_w + legend_col_width
    fig = plt.figure(figsize=(fig_w, height), constrained_layout=False)

    gs = GridSpec(
        rows, cols + 1, figure=fig,
        width_ratios=[1]*cols + [legend_col_width / (fig_w - legend_col_width)],
        wspace=0.02, hspace=0.02
    )

    # --- panels ---
    for i, sid in enumerate(sids):
        r, c = divmod(i, cols)
        ax = fig.add_subplot(gs[r, c])

        ad_sub = ad[ad.obs[groupby] == sid].copy()

        # force identical categories & palette on each subset
        ad_sub.obs[color] = ad_sub.obs[color].astype("category")
        ad_sub.obs[color] = ad_sub.obs[color].cat.set_categories(cats)
        ad_sub.uns[f"{color}_colors"] = cols_list

        sc.pl.spatial(
            ad_sub, color=color, spot_size=spot_size,
            show=False, ax=ax, legend_loc=None, frameon=False, title=str(sid)
        )
        ax.set_xlabel(""); ax.set_ylabel("")
        ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)

    # blank unused
    for j in range(n, rows*cols):
        r, c = divmod(j, cols)
        fig.add_subplot(gs[r, c]).axis("off")

    # --- legend ---
    ax_leg = fig.add_subplot(gs[:, -1])
    ax_leg.axis("off")
    handles = [
        Line2D([0], [0], marker='o', color='w',
               markerfacecolor=cols_list[k], markersize=7, label=str(cat))
        for k, cat in enumerate(cats)
    ]
    ax_leg.legend(handles=handles, title=color, frameon=False, loc="center left")

    fig.subplots_adjust(left=0.01, right=0.98, top=0.98, bottom=0.02,
                        wspace=0.02, hspace=0.02)

    return fig


def plot_spatial_highlight(
    ad,
    color,                      # obs column with categories
    highlight,                  # one category to highlight
    groupby="sample_id",
    spot_size=12,
    cols=3,
    height=7,
    flip_y=True,
    fix_limits=True,            # same axes for all panels
    base_color="#B0B0B0",       # de-emphasized others
    base_alpha=0.15,
    palette=None,               # dict {cat:"#hex"} or list
    rasterized=True,
    dpi=120,
):
    if "spatial" not in ad.obsm:
        raise ValueError("ad.obsm['spatial'] is required.")
    coords = np.asarray(ad.obsm["spatial"])[:, :2]

    # categories & palette
    cats = ad.obs[color].astype("category")
    cat_names = list(cats.cat.categories)
    if isinstance(palette, dict):
        col_list = [palette.get(c, "#377eb8") for c in cat_names]
    elif isinstance(palette, (list, tuple)) and len(palette) >= len(cat_names):
        col_list = list(palette)[:len(cat_names)]
    elif f"{color}_colors" in ad.uns and len(ad.uns[f"{color}_colors"]) == len(cat_names):
        col_list = list(ad.uns[f"{color}_colors"])
    else:
        from scanpy.pl import palettes
        base = palettes.default_64 if hasattr(palettes, "default_64") else palettes.default_102
        reps = int(np.ceil(len(cat_names) / len(base)))
        col_list = (base * reps)[:len(cat_names)]

    if highlight not in cat_names:
        raise ValueError(f"'{highlight}' not found in ad.obs['{color}'] categories.")
    h_idx = cat_names.index(highlight)
    hl_hex = col_list[h_idx]
    hl_rgba = np.array(mcolors.to_rgba(hl_hex), float)

    base_rgba = np.array(mcolors.to_rgba(base_color), float)
    base_rgba[3] = base_alpha

    cat_codes = cats.cat.codes.to_numpy()  # -1 for NaN
    is_highlight = (cat_codes == h_idx)

    colors_arr = np.empty((cat_codes.size, 4), dtype=float)
    colors_arr[:] = base_rgba
    colors_arr[is_highlight] = hl_rgba
    colors_arr[cat_codes < 0] = (0, 0, 0, 0)  # transparent for NA

    # group indices
    gseries = ad.obs[groupby]
    if pd.api.types.is_categorical_dtype(gseries):
        uniq_groups = list(gseries.cat.categories)
        uniq_groups = [g for g in uniq_groups if (gseries == g).any()]
        gcodes = pd.Categorical(gseries, categories=uniq_groups, ordered=True).codes
    else:
        uniq_groups = list(pd.unique(gseries))
        mapping = {g: i for i, g in enumerate(uniq_groups)}
        gcodes = gseries.map(mapping).to_numpy()

    group_indices = [np.flatnonzero(gcodes == i) for i in range(len(uniq_groups))]

    # common limits
    if fix_limits:
        x_min, y_min = coords.min(axis=0)
        x_max, y_max = coords.max(axis=0)

    # layout
    n = len(uniq_groups)
    rows = int(np.ceil(n / cols))
    panel_w = height * cols * 0.6 / rows
    fig_w = panel_w
    plt.ioff()
    fig = plt.figure(figsize=(fig_w, height), dpi=dpi, constrained_layout=False)
    gs = GridSpec(rows, cols, figure=fig, wspace=0.02, hspace=0.02)

    axes = []
    for i, sid in enumerate(uniq_groups):
        r, c = divmod(i, cols)
        ax = fig.add_subplot(gs[r, c])
        idx = group_indices[i]
        if idx.size:
            xy = coords[idx]
            ax.scatter(
                xy[:, 0], xy[:, 1],
                c=colors_arr[idx],
                s=spot_size,
                marker='o',
                linewidths=0,
                rasterized=rasterized
            )
        ax.set_title(f"{sid} — {highlight}", fontsize=9, pad=2)
        ax.set_aspect("equal")
        if flip_y:
            ax.invert_yaxis()
        if fix_limits:
            ax.set_xlim(x_min, x_max)
            if flip_y:
                ax.set_ylim(y_max, y_min)
            else:
                ax.set_ylim(y_min, y_max)
        ax.set_axis_off()
        axes.append(ax)

    # blank unused
    for j in range(n, rows * cols):
        r, c = divmod(j, cols)
        fig.add_subplot(gs[r, c]).axis("off")

    plt.ion()
    return fig, axes


def plot_spatial_genes_compact(
    ad,
    genes,
    groupby="sample_id",
    spot_size=18,
    cols=3,
    height=8,
    cmap="magma",          # any matplotlib colormap
    log1p=False,           # log1p-transform before plotting
    robust=True,           # use robust 2–98% limits instead of min/max
    vmin=None, vmax=None,  # override global color limits if you want
    nan_color="lightgray", # background for zero/NaN
    title_suffix="",
):
    """
    Plot per-sample spatial maps for one or multiple genes (continuous).
    Uses a shared color scale across all samples for each gene.
    NOTE: Requires Normalize and mcolors, which are imported at top.
    """

    genes = [genes] if isinstance(genes, str) else list(genes)
    sids  = list(ad.obs[groupby].unique())
    nS, nG = len(sids), len(genes)

    cols_total = nG + 1
    rows = int(np.ceil(nS))
    fig = plt.figure(
        figsize=(height * cols_total * 0.6, height),
        constrained_layout=False
    )
    gs = GridSpec(
        rows, cols_total, figure=fig,
        width_ratios=[1]*nG + [0.25],
        wspace=0.03, hspace=0.03
    )

    for gi, g in enumerate(genes):
        if g not in ad.var_names:
            print(f"⚠️ Gene '{g}' not in var_names; skipping.")
            continue

        vals = ad[:, g].X
        vals = vals.toarray().ravel() if hasattr(vals, "toarray") else np.asarray(vals).ravel()
        if log1p:
            vals = np.log1p(vals)

        if vmin is None or vmax is None:
            if robust:
                lo, hi = np.nanpercentile(vals, [2, 98])
            else:
                lo, hi = np.nanmin(vals), np.nanmax(vals)
        else:
            lo, hi = vmin, vmax
        if hi <= lo:
            hi = lo + 1e-9

        norm = Normalize(vmin=lo, vmax=hi)

        for si, sid in enumerate(sids):
            ax = fig.add_subplot(gs[si, gi]) if rows > 1 else fig.add_subplot(gs[0, gi])
            ad_sub = ad[ad.obs[groupby] == sid].copy()

            v = ad_sub[:, g].X
            v = v.toarray().ravel() if hasattr(v, "toarray") else np.asarray(v).ravel()
            if log1p:
                v = np.log1p(v)

            cmap_obj = plt.get_cmap(cmap)
            rgba = cmap_obj(norm(v))

            mask = ~np.isfinite(v) | (v <= 0)
            rgba[mask, :3] = mcolors.to_rgb(nan_color)
            rgba[mask, 3]  = 1.0

            hex_colors = [mcolors.to_hex(c) for c in rgba]
            ad_sub.obs[f"__{g}__colors"] = hex_colors

            sc.pl.spatial(
                ad_sub,
                color=f"__{g}__colors",
                spot_size=spot_size,
                show=False,
                ax=ax,
                title=f"{sid}" + (f" • {title_suffix}" if title_suffix else ""),
                frameon=False,
                legend_loc=None,
            )
            ax.set_xlabel(""); ax.set_ylabel("")
            ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)

        # colorbar column
        cb_ax = fig.add_subplot(gs[:, -1]) if rows > 1 else fig.add_subplot(gs[0, -1])
        if gi == nG - 1:
            cb_ax.axis("off")
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cbar = fig.colorbar(sm, ax=cb_ax, fraction=0.8)
            label_txt = f"{'log1p(' if log1p else ''}{genes[-1]}{')' if log1p else ''}"
            cbar.set_label(label_txt, rotation=90)

    fig.subplots_adjust(left=0.02, right=0.98, top=0.96, bottom=0.02,
                        wspace=0.02, hspace=0.02)
    plt.show()


def plot_spatial_any(ad, feature, **kwargs):
    """
    Convenience dispatcher:
      - If `feature` is a column in .obs (categorical), call plot_spatial_compact
      - If `feature` is a gene (in .var_names), call plot_spatial_genes_compact
    """
    if feature in ad.obs.columns and pd.api.types.is_categorical_dtype(ad.obs[feature]):
        return plot_spatial_compact(ad, color=feature, **kwargs)
    elif feature in ad.var_names:
        return plot_spatial_genes_compact(ad, genes=feature, **kwargs)
    else:
        raise ValueError(f"'{feature}' not found as categorical obs or gene.")


def _flip_y_inplace(ad):
    """Flip Y in .obsm['spatial'] so 'up' is up."""
    XY = ad.obsm["spatial"].copy()
    XY[:, 1] = XY[:, 1].max() - XY[:, 1]
    ad.obsm["spatial"] = XY
    return ad


def plot_compartment_bar(df: pd.DataFrame, top_n=None, title=None, ax=None):
    """
    Horizontal barplot of composition (percent).
    Accepts the DataFrame returned by `compartment_composition`.
    """
    if top_n is not None:
        d = df.head(top_n)
    else:
        d = df

    if ax is None:
        h = max(3, 0.35 * len(d))
        fig, ax = plt.subplots(figsize=(7, h))

    ax.barh(d.index[::-1], d["percent"][::-1])
    for i, (y, v) in enumerate(zip(d.index[::-1], d["percent"][::-1])):
        ax.text(v + 0.5, i, f"{v:.1f}%", va="center", fontsize=9)

    comp = df.attrs.get("compartment", "")
    groupby = df.attrs.get("groupby", "")
    total = df.attrs.get("total", None)
    ttl = title or f"{comp} • {groupby} composition"
    if total is not None:
        ttl += f"  (n={total})"

    ax.set_title(ttl)
    ax.set_xlabel("Percent of cells")
    ax.set_ylabel(groupby)
    ax.set_xlim(0, max(100, float(d["percent"].max()) + 5))
    ax.spines[['top','right']].set_visible(False)
    plt.tight_layout()
    return ax


def extract_dotplot_data_from_adata(
    adata,
    genes,
    groupby="rbd_annotation",
    layer=None,               # e.g. "log1p", else uses .X
    thresh=0.0,               # expression > thresh counts as "expressing"
    standard_scale=None,      # None or "var" (min–max per gene)
    return_tidy=True,
    to_csv_prefix=None
):
    """
    Compute per-group mean expression and fraction expressing for selected genes.
    Returns:
      {"mean": DataFrame[genes x groups],
       "frac": DataFrame[genes x groups],
       "mean_scaled": (optional),
       "tidy": (optional)}
    """
    present = [g for g in genes if g in adata.var_names]
    missing = sorted(set(genes) - set(present))
    if missing:
        print(f"[warn] {len(missing)} genes not found and will be skipped: {missing}")
    if not present:
        raise ValueError("None of the requested genes are in adata.var_names.")

    groups = adata.obs[groupby].astype("category")
    group_names = list(groups.cat.categories)

    X_all = adata.layers[layer] if layer is not None else adata.X
    is_sparse = sp.issparse(X_all)

    mean_mat = np.zeros((len(present), len(group_names)), dtype=float)
    frac_mat = np.zeros((len(present), len(group_names)), dtype=float)

    idx_cols = adata.var_names.get_indexer(present)

    for j, g in enumerate(group_names):
        mask = (groups.values == g)
        if mask.sum() == 0:
            continue
        Xg = X_all[mask][:, idx_cols]

        if is_sparse:
            means = np.asarray(Xg.mean(axis=0)).ravel()
            frac  = np.asarray((Xg > thresh).mean(axis=0)).ravel()
        else:
            means = Xg.mean(axis=0)
            frac  = (Xg > thresh).mean(axis=0)

        mean_mat[:, j] = means
        frac_mat[:, j] = frac

    mean_df = pd.DataFrame(mean_mat, index=present, columns=group_names)
    frac_df = pd.DataFrame(frac_mat, index=present, columns=group_names)

    out = {"mean": mean_df, "frac": frac_df}

    if standard_scale == "var":
        mm = (mean_df.T - mean_df.min(axis=1)) / (
            mean_df.max(axis=1) - mean_df.min(axis=1) + 1e-12
        )
        mean_scaled = mm.T
        out["mean_scaled"] = mean_scaled

    if return_tidy:
        base = out.get("mean_scaled", mean_df)
        left = (
            base.reset_index()
                .rename(columns={"index": "gene"})
                .melt(id_vars="gene", var_name=groupby, value_name="mean_scaled_or_mean")
        )
        right = (
            frac_df.reset_index()
                   .rename(columns={"index": "gene"})
                   .melt(id_vars="gene", var_name=groupby, value_name="frac")
        )
        tidy = left.merge(right, on=["gene", groupby], how="left")
        out["tidy"] = tidy

    if to_csv_prefix:
        mean_df.to_csv(f"{to_csv_prefix}.mean.csv")
        frac_df.to_csv(f"{to_csv_prefix}.frac.csv")
        if "mean_scaled" in out:
            out["mean_scaled"].to_csv(f"{to_csv_prefix}.mean_scaled.csv")
        if "tidy" in out:
            out["tidy"].to_csv(f"{to_csv_prefix}.tidy.csv", index=False)

    return out


def plot_spatial_grid_by_condition_age(
    adata,
    gene,
    *,
    layer=None,
    use_raw=False,
    spot_size=30,
    cmap="coolwarm",
    percentile=99,
    conditions=None,
    ages=None,
    sample_index=0,
    condition_col="condition",
    age_col="age",
    sample_col="sample_id",
    title_fmt="{cond} | {age} | {sid}",
    vmin=0.0,
    figsize_scale=5.0,
    tight=True,
):
    """
    Plot a grid of spatial maps: rows = condition, cols = age.
    Uses a shared color scale computed from the chosen matrix/layer (or .raw).
    Returns fig, axes.
    """
    def _to_dense(x):
        return x.toarray() if sp.issparse(x) else np.asarray(x)

    def _get_gene_values(adata_local, gene_local, layer=None, use_raw=False):
        if use_raw:
            if adata_local.raw is None:
                raise ValueError("use_raw=True but adata.raw is None.")
            if gene_local not in adata_local.raw.var_names:
                raise ValueError(f"{gene_local} not in adata.raw.var_names.")
            mat = adata_local.raw[:, gene_local].X
            return _to_dense(mat).ravel()
        if gene_local not in adata_local.var_names:
            raise ValueError(f"{gene_local} not in adata.var_names.")
        if layer is None:
            mat = adata_local[:, gene_local].X
        else:
            if layer not in adata_local.layers:
                raise ValueError(f"Layer '{layer}' not in adata.layers.")
            mat = adata_local[:, gene_local].layers[layer]
        return _to_dense(mat).ravel()

    for col in (sample_col, condition_col, age_col):
        if col not in adata.obs.columns:
            raise ValueError(f"'{col}' not found in adata.obs.")
    if "spatial" not in adata.obsm_keys():
        raise ValueError("No 'spatial' coords in adata.obsm.")

    cond_series = adata.obs[condition_col].astype(str)
    age_series  = adata.obs[age_col].astype(str)

    if ages is None:
        try:
            uniq_ages = sorted(pd.to_numeric(age_series.unique()))
            ages = list(map(str, uniq_ages))
        except Exception:
            ages = sorted(age_series.unique(), key=lambda x: (len(x), x))
    else:
        ages = list(map(str, ages))

    if conditions is None:
        conditions = sorted(cond_series.unique(), key=lambda x: (len(x), x))
    else:
        conditions = list(map(str, conditions))

    df_keys = adata.obs[[sample_col, condition_col, age_col]].copy()
    df_keys[condition_col] = df_keys[condition_col].astype(str)
    df_keys[age_col]       = df_keys[age_col].astype(str)

    comb2samples = {}
    for sid, grp in df_keys.groupby(sample_col):
        c = grp[condition_col].iloc[0]
        a = grp[age_col].iloc[0]
        comb2samples.setdefault((c, a), []).append(sid)

    vals = _get_gene_values(adata, gene, layer=layer, use_raw=use_raw)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        raise ValueError(f"No finite values for {gene} in selected matrix/layer.")
    vmax = np.percentile(vals, percentile)

    n_rows, n_cols = len(conditions), len(ages)
    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(figsize_scale * n_cols, figsize_scale * n_rows),
        squeeze=False
    )

    mappable = None
    for i, cond in enumerate(conditions):
        for j, age in enumerate(ages):
            ax = axes[i, j]
            sids = comb2samples.get((cond, age), [])
            if not sids:
                ax.axis("off")
                ax.set_title(f"{cond} | {age}\n(no sample)", fontsize=10)
                continue

            idx = min(sample_index, len(sids) - 1)
            sid = sids[idx]
            ad_int = adata[adata.obs[sample_col] == sid]

            sc.pl.spatial(
                ad_int,
                color=gene,
                layer=layer,
                spot_size=spot_size,
                vmin=vmin,
                vmax=vmax,
                cmap=cmap,
                show=False,
                ax=ax,
            )
            ax.set_title(title_fmt.format(cond=cond, age=age, sid=sid), fontsize=10)

            if mappable is None:
                ims = ax.get_images()
                if ims:
                    mappable = ims[0]

    if mappable is not None:
        cbar = fig.colorbar(mappable, ax=axes, fraction=0.02, pad=0.02)
        cbar.set_label(gene)

    if tight:
        plt.tight_layout()
    plt.show()
    return fig, axes


def plot_spatial_grid_basic(
    adata,
    *,
    condition_col="condition",
    age_col="age",
    sample_col="sample_id",
    conditions=None,
    ages=None,
    sample_index=0,
    figsize_scale=5.0,
    spot_size=30,
    # highlight mode
    highlight_col=None,
    highlight_value=None,
    highlight_color="tab:red",
    background_grey="#B0B0B0",
    background_alpha=1.0,
    title_fmt="{cond} | {age} | {sid}",
    invert_y=True,
):
    for col in (condition_col, age_col, sample_col):
        if col not in adata.obs.columns:
            raise ValueError(f"'{col}' not found in adata.obs.")
    if "spatial" not in adata.obsm_keys():
        raise ValueError("No 'spatial' coordinates in adata.obsm.")

    cond_series = adata.obs[condition_col].astype(str)
    age_series  = adata.obs[age_col].astype(str)

    if ages is None:
        try:
            uniq_ages = sorted(pd.to_numeric(age_series.unique()))
            ages = list(map(str, uniq_ages))
        except Exception:
            ages = sorted(age_series.unique(), key=lambda x: (len(x), x))
    else:
        ages = list(map(str, ages))

    if conditions is None:
        conditions = sorted(cond_series.unique(), key=lambda x: (len(x), x))

    key_df = adata.obs[[sample_col, condition_col, age_col]].astype(str)
    comb2samples = {}
    for sid, grp in key_df.groupby(sample_col):
        c = grp[condition_col].iloc[0]
        a = grp[age_col].iloc[0]
        comb2samples.setdefault((c, a), []).append(sid)

    n_rows, n_cols = len(conditions), len(ages)
    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(figsize_scale * n_cols, figsize_scale * n_rows),
        squeeze=False
    )

    for i, cond in enumerate(conditions):
        for j, age in enumerate(ages):
            ax = axes[i, j]
            sids = comb2samples.get((cond, age), [])
            if not sids:
                ax.axis("off")
                ax.set_title(f"{cond} | {age}\n(no sample)", fontsize=10)
                continue

            sid = sids[min(sample_index, len(sids) - 1)]
            ad_sub = adata[adata.obs[sample_col] == sid]
            xy = ad_sub.obsm["spatial"]
            if xy is None or len(xy) == 0:
                ax.axis("off")
                ax.set_title(f"{cond} | {age} | {sid}\n(no coords)", fontsize=10)
                continue

            # background: all cells in grey
            ax.scatter(
                xy[:, 0], xy[:, 1],
                s=spot_size,
                c=background_grey,
                alpha=background_alpha,
                edgecolors="none",
                zorder=1,
            )

            # overlay highlight
            if highlight_col is not None and highlight_value is not None:
                if highlight_col not in ad_sub.obs.columns:
                    ax.set_title(
                        f"{cond} | {age} | {sid}\n(no '{highlight_col}')",
                        fontsize=10
                    )
                else:
                    mask = ad_sub.obs[highlight_col].astype(str).values == str(highlight_value)
                    if mask.any():
                        ax.scatter(
                            xy[mask, 0], xy[mask, 1],
                            s=spot_size,
                            c=highlight_color,
                            edgecolors="black",
                            linewidths=0.2,
                            zorder=2,
                        )

            ax.set_title(title_fmt.format(cond=cond, age=age, sid=sid), fontsize=10)
            ax.set_aspect("equal")
            ax.axis("off")
            if invert_y:
                ax.invert_yaxis()

    plt.tight_layout()
    plt.show()
    return fig, axes