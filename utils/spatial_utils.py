"""
Auto-generated utilities for spatial_utils.
Do not edit by hand without moving changes back into notebooks.

Each function below was extracted from exported analysis notebooks.
"""

from typing import *
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
from scipy import sparse

def plot_spatial_compact_fast(
    ad,
    color="leiden_2",          # obs column (categorical) OR a gene name
    groupby="sample_id",
    spot_size=8,
    cols=3,
    height=8,
    legend_col_width=1.2,
    palette=None,              # for categorical only: dict {cat:"#hex"} or list
    rasterized=True,
    invert_y=True,
    dpi=120,
    # --- extra knobs (used when color is a gene) ---
    cmap="viridis",
    vmin=None,
    vmax=None,
    na_alpha=0.0               # alpha for NaN gene values (0 = fully transparent)
):
    """
    If `color` is an obs column -> categorical plot with legend.
    If `color` is a gene present in ad.var_names (or ad.raw.var_names) -> continuous plot with colorbar.
    """
    # 0) Preconditions ---------------------------------------------------------
    if "spatial" not in ad.obsm:
        raise ValueError("ad.obsm['spatial'] is required.")

    coords = np.asarray(ad.obsm["spatial"])[:, :2]

    has_raw = hasattr(ad, "raw") and (ad.raw is not None)
    is_obs  = color in ad.obs.columns
    is_gene = (color in ad.var_names) or (has_raw and color in ad.raw.var_names)

    if not (is_obs or is_gene):
        raise KeyError(f"'{color}' not found as obs column or gene name.")

    # 1) Grouping layout -------------------------------------------------------
    gvals = ad.obs[groupby].astype(str).to_numpy()
    uniq_groups, gcodes = np.unique(gvals, return_inverse=True)
    group_indices = [np.flatnonzero(gcodes == gi) for gi in range(len(uniq_groups))]

    n = len(uniq_groups)
    rows = int(np.ceil(n / cols))
    panel_w = height * cols * 0.6 / max(rows, 1)
    fig_w = panel_w + legend_col_width

    plt.ioff()
    fig = plt.figure(figsize=(fig_w, height), dpi=dpi, constrained_layout=False)
    gs = GridSpec(
        rows, cols + 1, figure=fig,
        width_ratios=[1]*cols + [legend_col_width / max(1e-9, (fig_w - legend_col_width))],
        wspace=0.02, hspace=0.02
    )

    # 2A) Categorical path -----------------------------------------------------
    if is_obs:
        cats = ad.obs[color].astype("category")
        cat_names = cats.cat.categories
        cat_codes = cats.cat.codes.to_numpy()  # -1 for NaN

        # palette
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
            base = (mpl.cm.get_cmap("tab20").colors
                    if hasattr(mpl.cm.get_cmap("tab20"), "colors")
                    else list(__import__("scanpy").pl.palettes.default_64))
            base = list(base)
            reps = int(np.ceil(len(cat_names) / len(base)))
            col_list = (base * reps)[:len(cat_names)]

        # store for consistency
        ad.uns[f"{color}_colors"] = col_list

        # fast RGBA lookup
        rgba = np.array([mcolors.to_rgba(c) for c in col_list], dtype=float)
        colors_arr = np.empty((cat_codes.size, 4), dtype=float)
        mask_valid = cat_codes >= 0
        colors_arr[mask_valid] = rgba[cat_codes[mask_valid]]
        colors_arr[~mask_valid] = (0, 0, 0, 0)  # transparent for NA

        # panels
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
            ax.set_title(str(sid), fontsize=9, pad=2)
            ax.set_aspect("equal")
            if invert_y: ax.invert_yaxis()
            ax.set_axis_off()

        # blank unused
        for j in range(n, rows * cols):
            r, c = divmod(j, cols)
            fig.add_subplot(gs[r, c]).axis("off")

        # legend column
        ax_leg = fig.add_subplot(gs[:, -1])
        ax_leg.axis("off")
        handles = [
            Line2D([0], [0], marker='o', color='w',
                   markerfacecolor=col_list[k], markersize=7, label=str(cat))
            for k, cat in enumerate(cat_names)
        ]
        ax_leg.legend(handles=handles, title=color, frameon=False, loc="center left")

    # 2B) Gene (continuous) path ----------------------------------------------
    else:
        # pull expression vector
        expr = (ad.raw[:, color].X if has_raw and color in ad.raw.var_names else ad[:, color].X)
        expr = np.asarray(expr).squeeze()  # dense 1D
        # colormap normalization
        finite = np.isfinite(expr)
        vmin_eff = np.nanmin(expr[finite]) if vmin is None else vmin
        vmax_eff = np.nanmax(expr[finite]) if vmax is None else vmax
        norm = mpl.colors.Normalize(vmin=vmin_eff, vmax=vmax_eff)
        cmap_obj = mpl.cm.get_cmap(cmap)

        # panels
        for i, sid in enumerate(uniq_groups):
            r, c = divmod(i, cols)
            ax = fig.add_subplot(gs[r, c])
            idx = group_indices[i]
            if idx.size:
                xy = coords[idx]
                vals = expr[idx]
                # mask NaNs -> transparent by alpha
                alphas = np.where(np.isfinite(vals), 1.0, na_alpha)
                # build RGBA from cmap + alpha
                col_rgba = cmap_obj(norm(np.nan_to_num(vals, nan=vmin_eff)))
                col_rgba[:, 3] = alphas
                ax.scatter(
                    xy[:, 0], xy[:, 1],
                    c=col_rgba,
                    s=spot_size,
                    marker='o',
                    linewidths=0,
                    rasterized=rasterized
                )
            ax.set_title(str(sid), fontsize=9, pad=2)
            ax.set_aspect("equal")
            if invert_y: ax.invert_yaxis()
            ax.set_axis_off()

        # blank unused
        for j in range(n, rows * cols):
            r, c = divmod(j, cols)
            fig.add_subplot(gs[r, c]).axis("off")

        # colorbar in legend column
        cax = fig.add_subplot(gs[:, -1])
        cax.set_visible(True)
        sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap_obj)
        sm.set_array([])
        cb = fig.colorbar(sm, cax=cax)
        cb.set_label(color)

    # finish
    fig.subplots_adjust(left=0.01, right=0.98, top=0.98, bottom=0.02, wspace=0.02, hspace=0.02)
    plt.ion()
    plt.show()
    return fig

def plot_spatial_compact(
    ad,
    color="leiden_2",
    groupby="sample_id",
    spot_size=18,
    cols=3,
    height=8,
    legend_col_width=1.2,
    palette=None,          # <- NEW: dict {category: "#hex"} or list in desired order
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
        # reuse existing palette if present on the full object
        cols_list = list(ad.uns[f"{color}_colors"])
        if len(cols_list) != len(cats):
            raise ValueError(f"{color}_colors length does not match categories.")
    else:
        base = sc.pl.palettes.default_64 if hasattr(sc.pl.palettes, "default_64") else sc.pl.palettes.default_102
        reps = int(np.ceil(len(cats) / len(base)))
        cols_list = (base * reps)[:len(cats)]

    # store on the full object too (handy if you plot elsewhere later)
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

        # CRUCIAL: force identical categories & palette on each subset
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

    fig.subplots_adjust(left=0.01, right=0.98, top=0.98, bottom=0.02, wspace=0.02, hspace=0.02)

    # RETURN instead of showing
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
    palette=None,               # dict {cat:"#hex"} or list; if None try ad.uns[f"{color}_colors"]
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
        # fallback: simple cycle
        from scanpy.pl import palettes
        base = palettes.default_64 if hasattr(palettes, "default_64") else palettes.default_102
        reps = int(np.ceil(len(cat_names) / len(base)))
        col_list = (base * reps)[:len(cat_names)]

    # pick highlight color
    if highlight not in cat_names:
        raise ValueError(f"'{highlight}' not found in ad.obs['{color}'] categories.")
    h_idx = cat_names.index(highlight)
    hl_hex = col_list[h_idx]
    hl_rgba = np.array(mcolors.to_rgba(hl_hex), float)

    # precompute color array: highlight vs base
    base_rgba = np.array(mcolors.to_rgba(base_color), float)
    base_rgba[3] = base_alpha

    cat_codes = cats.cat.codes.to_numpy()  # -1 for NaN
    is_highlight = (cat_codes == h_idx)
    # initialize all as base, then overwrite highlights
    colors_arr = np.empty((cat_codes.size, 4), dtype=float)
    colors_arr[:] = base_rgba
    colors_arr[is_highlight] = hl_rgba
    # transparent for NaN
    colors_arr[cat_codes < 0] = (0, 0, 0, 0)

    # group indices (preserve observed order if possible)
    gseries = ad.obs[groupby]
    if pd.api.types.is_categorical_dtype(gseries):
        uniq_groups = list(gseries.cat.categories)
        # keep only those present
        uniq_groups = [g for g in uniq_groups if (gseries == g).any()]
        gcodes = pd.Categorical(gseries, categories=uniq_groups, ordered=True).codes
    else:
        uniq_groups = list(pd.unique(gseries))  # preserves first-appearance order
        mapping = {g:i for i,g in enumerate(uniq_groups)}
        gcodes = gseries.map(mapping).to_numpy()
    group_indices = [np.flatnonzero(gcodes == i) for i in range(len(uniq_groups))]

    # common limits (global) if requested
    if fix_limits:
        x_min, y_min = coords.min(axis=0)
        x_max, y_max = coords.max(axis=0)

    # layout
    n = len(uniq_groups)
    rows = int(np.ceil(n / cols))
    panel_w = height * cols * 0.6 / rows
    fig_w = panel_w  # legend-free; keep compact
    plt.ioff()
    fig = plt.figure(figsize=(fig_w, height), dpi=dpi, constrained_layout=False)
    gs = GridSpec(rows, cols, figure=fig, wspace=0.02, hspace=0.02)

    # panels
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
            ax.set_ylim(y_max, y_min) if flip_y else ax.set_ylim(y_min, y_max)
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
    """

    # ensure genes is a list
    genes = [genes] if isinstance(genes, str) else list(genes)
    sids  = list(ad.obs[groupby].unique())
    nS, nG = len(sids), len(genes)

    # grid: rows = samples, cols = genes + 1 legend column
    cols_total = nG + 1
    rows = int(np.ceil(nS))
    fig = plt.figure(figsize=(height * cols_total * 0.6, height), constrained_layout=False)
    gs  = GridSpec(rows, cols_total, figure=fig, width_ratios=[1]*nG + [0.25],
                   wspace=0.03, hspace=0.03)

    for gi, g in enumerate(genes):
        if g not in ad.var_names:
            print(f"⚠️ Gene '{g}' not in var_names; skipping.")
            continue

        # collect values across ALL samples to set shared limits
        vals = ad[:, g].X
        vals = vals.toarray().ravel() if hasattr(vals, "toarray") else np.asarray(vals).ravel()
        if log1p: vals = np.log1p(vals)

        if vmin is None or vmax is None:
            if robust:
                lo, hi = np.nanpercentile(vals, [2, 98])
            else:
                lo, hi = np.nanmin(vals), np.nanmax(vals)
        else:
            lo, hi = vmin, vmax
        if hi <= lo:  # degenerate case
            hi = lo + 1e-9
        norm = Normalize(vmin=lo, vmax=hi)

        # draw per-sample panels for this gene
        for si, sid in enumerate(sids):
            ax = fig.add_subplot(gs[si, gi]) if rows > 1 else fig.add_subplot(gs[0, gi])
            ad_sub = ad[ad.obs[groupby] == sid].copy()

            # precompute color values so Scanpy respects global norm
            v = ad_sub[:, g].X
            v = v.toarray().ravel() if hasattr(v, "toarray") else np.asarray(v).ravel()
            if log1p: v = np.log1p(v)
            # map to colors
            cmap_obj = plt.get_cmap(cmap)
            rgba = cmap_obj(norm(v))
            # make zeros / NaNs pale
            mask = ~np.isfinite(v) | (v <= 0)
            rgba[mask, :3] = plt.matplotlib.colors.to_rgb(nan_color)  # set color
            rgba[mask, 3]  = 1.0                                     # opaque

            # Scanpy trick: provide precomputed colors via obs (as hex)
            import matplotlib.colors as mcolors
            hex_colors = [mcolors.to_hex(c) for c in rgba]
            ad_sub.obs[f"__{g}__colors"] = hex_colors

            sc.pl.spatial(
                ad_sub,
                color=f"__{g}__colors",           # plot with our precomputed colors
                spot_size=spot_size,
                show=False,
                ax=ax,
                title=f"{sid}" + (f" • {title_suffix}" if title_suffix else ""),
                frameon=False,
                legend_loc=None,
            )
            ax.set_xlabel(""); ax.set_ylabel("")
            ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)

        # add a colorbar column only once per gene
        cb_ax = fig.add_subplot(gs[:, -1]) if rows > 1 else fig.add_subplot(gs[0, -1])
        if gi == nG - 1:  # only draw once (last gene) to avoid redraw
            cb_ax.axis("off")
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cbar = fig.colorbar(sm, ax=cb_ax, fraction=0.8)
            cbar.set_label(f"{'log1p(' if log1p else ''}{genes[-1]}{')' if log1p else ''}", rotation=90)

    fig.subplots_adjust(left=0.02, right=0.98, top=0.96, bottom=0.02, wspace=0.02, hspace=0.02)
    plt.show()

def plot_spatial_any(
    ad,
    feature,
    **kwargs
):
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
    Horizontal barplot of composition (percent). Accepts the DataFrame returned by `compartment_composition`.
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
    to_csv_prefix=None        # e.g. "dotplot"; will write 2 CSVs if set
):
    """
    Compute per-group mean expression and fraction expressing for selected genes.
    Returns:
      {"mean": DataFrame[genes x groups],
       "frac": DataFrame[genes x groups],
       "mean_scaled": (optional) DataFrame if standard_scale="var",
       "tidy": (optional) long-form DataFrame}
    """
    # 1) restrict to genes present
    present = [g for g in genes if g in adata.var_names]
    missing = sorted(set(genes) - set(present))
    if missing:
        print(f"[warn] {len(missing)} genes not found and will be skipped: {missing}")

    if not present:
        raise ValueError("None of the requested genes are in adata.var_names.")

    groups = adata.obs[groupby].astype("category")
    group_names = list(groups.cat.categories)

    # 2) choose matrix
    X_all = adata.layers[layer] if layer is not None else adata.X
    is_sparse = sp.issparse(X_all)

    # 3) build result holders
    mean_mat = np.zeros((len(present), len(group_names)), dtype=float)
    frac_mat = np.zeros((len(present), len(group_names)), dtype=float)

    # 4) compute per group
    for j, g in enumerate(group_names):
        mask = (groups.values == g)
        if mask.sum() == 0:
            continue
        # slice rows by group, columns by genes
        idx_cols = adata.var_names.get_indexer(present)  # integer indices
        Xg = X_all[mask][:, idx_cols]

        # means
        if is_sparse:
            # sparse mean keeps shape (1, n) — convert to 1D
            means = np.asarray(Xg.mean(axis=0)).ravel()
            # fraction expressing
            frac  = np.asarray((Xg > thresh).mean(axis=0)).ravel()
        else:
            means = Xg.mean(axis=0)
            frac  = (Xg > thresh).mean(axis=0)

        mean_mat[:, j] = means
        frac_mat[:, j] = frac

    mean_df = pd.DataFrame(mean_mat, index=present, columns=group_names)
    frac_df = pd.DataFrame(frac_mat, index=present, columns=group_names)

    out = {"mean": mean_df, "frac": frac_df}

    # 5) optional standard_scale="var" → min–max per gene across groups
    if standard_scale == "var":
        mm = (mean_df.T - mean_df.min(axis=1)) / (mean_df.max(axis=1) - mean_df.min(axis=1) + 1e-12)
        mean_scaled = mm.T
        out["mean_scaled"] = mean_scaled

    # 6) optional tidy long format
    # 6) optional tidy long format  ✅ fixed merge keys
    if return_tidy:
        base = out.get("mean_scaled", mean_df)  # what you'd map to color
        left = (
            base.reset_index()                       # has columns: ["gene", groups...]
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

    # 7) optional CSVs
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
    layer=None,           # e.g. "log1p", "counts"; None -> .X
    use_raw=False,        # use adata.raw
    spot_size=30,
    cmap="coolwarm",
    percentile=99,        # global vmax percentile across ALL cells
    conditions=None,      # optional subset/list to order conditions (rows)
    ages=None,            # optional subset/list to order ages (cols)
    sample_index=0,       # which sample to show if multiple per (cond, age)
    condition_col="condition",
    age_col="age",
    sample_col="sample_id",
    title_fmt="{cond} | {age} | {sid}",
    vmin=0.0,
    figsize_scale=5.0,    # figure size scaling per row/col
    tight=True,
):
    """
    Plot a grid of spatial maps: rows = condition, cols = age.
    Uses a shared color scale computed from the chosen matrix/layer (or .raw).

    Returns
    -------
    fig, axes : matplotlib Figure and 2D axes array
    """
    # ---- helpers ----
    def _to_dense(x):
        return x.toarray() if sp.issparse(x) else np.asarray(x)

    def _get_gene_values(adata, gene, layer=None, use_raw=False):
        if use_raw:
            if adata.raw is None:
                raise ValueError("use_raw=True but adata.raw is None.")
            if gene not in adata.raw.var_names:
                raise ValueError(f"{gene} not in adata.raw.var_names.")
            mat = adata.raw[:, gene].X
            return _to_dense(mat).ravel()
        if gene not in adata.var_names:
            raise ValueError(f"{gene} not in adata.var_names.")
        if layer is None:
            mat = adata[:, gene].X
        else:
            if layer not in adata.layers:
                raise ValueError(f"Layer '{layer}' not in adata.layers: {list(adata.layers.keys())}")
            mat = adata[:, gene].layers[layer]
        return _to_dense(mat).ravel()

    # ---- basic sanity ----
    if sample_col not in adata.obs.columns:
        raise ValueError(f"'{sample_col}' not found in adata.obs.")
    if condition_col not in adata.obs.columns:
        raise ValueError(f"'{condition_col}' not found in adata.obs.")
    if age_col not in adata.obs.columns:
        raise ValueError(f"'{age_col}' not found in adata.obs.")

    # ---- labels & ordering ----
    cond_series = adata.obs[condition_col].astype(str)
    age_series  = adata.obs[age_col].astype(str)

    # ages: try numeric sort if possible
    if ages is None:
        try:
            uniq_ages = sorted(pd.to_numeric(age_series.unique()))
            ages = list(map(str, uniq_ages))
        except Exception:
            ages = sorted(age_series.unique(), key=lambda x: (len(x), x))
    else:
        ages = list(map(str, ages))

    # conditions: simple sorted unless user supplied
    if conditions is None:
        conditions = sorted(cond_series.unique(), key=lambda x: (len(x), x))
    else:
        conditions = list(map(str, conditions))

    # ---- map (cond, age) -> list of sample_ids
    df_keys = adata.obs[[sample_col, condition_col, age_col]].copy()
    df_keys[condition_col] = df_keys[condition_col].astype(str)
    df_keys[age_col]       = df_keys[age_col].astype(str)

    comb2samples = {}
    for sid, grp in df_keys.groupby(sample_col):
        c = grp[condition_col].iloc[0]
        a = grp[age_col].iloc[0]
        comb2samples.setdefault((c, a), []).append(sid)

    # ---- shared color scale from matrix values ----
    vals = _get_gene_values(adata, gene, layer=layer, use_raw=use_raw)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        raise ValueError(f"No finite values for {gene} in selected matrix/layer.")
    vmax = np.percentile(vals, percentile)

    # ---- figure grid ----
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

            # choose sample
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

    # ---- single global colorbar ----
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
    highlight_col=None,        # e.g. "celltype_merged"
    highlight_value=None,      # e.g. "Microglia"
    highlight_color="tab:red",
    background_grey="#B0B0B0", # darker grey so it shows
    background_alpha=1.0,      # 1.0 = fully opaque background
    title_fmt="{cond} | {age} | {sid}",
    invert_y=True,             # Visium-style orientation
):
    # --- sanity
    for col in (condition_col, age_col, sample_col):
        if col not in adata.obs.columns:
            raise ValueError(f"'{col}' not found in adata.obs.")
    if "spatial" not in adata.obsm_keys():
        raise ValueError("No 'spatial' coordinates in adata.obsm.")

    # order axes
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

    # map (cond, age) -> sample_ids
    key_df = adata.obs[[sample_col, condition_col, age_col]].astype(str)
    comb2samples = {}
    for sid, grp in key_df.groupby(sample_col):
        c = grp[condition_col].iloc[0]; a = grp[age_col].iloc[0]
        comb2samples.setdefault((c, a), []).append(sid)

    # figure
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
                ax.axis("off"); ax.set_title(f"{cond} | {age}\n(no sample)", fontsize=10)
                continue

            sid = sids[min(sample_index, len(sids) - 1)]
            ad_sub = adata[adata.obs[sample_col] == sid]
            xy = ad_sub.obsm["spatial"]
            if xy is None or len(xy) == 0:
                ax.axis("off"); ax.set_title(f"{cond} | {age} | {sid}\n(no coords)", fontsize=10)
                continue

            # --- background: ALL cells in grey (opaque unless you change background_alpha)
            ax.scatter(
                xy[:, 0], xy[:, 1],
                s=spot_size,
                c=background_grey,
                alpha=background_alpha,
                edgecolors="none",
                zorder=1,
            )

            # --- overlay highlight (optional)
            if highlight_col is not None and highlight_value is not None:
                if highlight_col not in ad_sub.obs.columns:
                    ax.set_title(f"{cond} | {age} | {sid}\n(no '{highlight_col}')", fontsize=10)
                else:
                    mask = ad_sub.obs[highlight_col].astype(str).values == str(highlight_value)
                    if mask.any():
                        ax.scatter(
                            xy[mask, 0], xy[mask, 1],
                            s=spot_size,
                            c=highlight_color,
                            edgecolors="black",  # thin outline helps pop
                            linewidths=0.2,
                            zorder=2,
                        )

            ax.set_title(title_fmt.format(cond=cond, age=age, sid=sid), fontsize=10)
            ax.set_aspect("equal")
            ax.axis("off")
            if invert_y:
                ax.invert_yaxis()  # match Scanpy/Visium orientation

    plt.tight_layout()
    plt.show()
    return fig, axes
