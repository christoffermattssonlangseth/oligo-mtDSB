#!/usr/bin/env python
# coding: utf-8

# In[ ]:





# In[2]:


# define functions
import scanpy as sc
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib import colors as mcolors
import warnings
warnings.filterwarnings("ignore")


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib import colors as mcolors
import matplotlib as mpl

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

import os
import numpy as np
import pandas as pd
from anndata import AnnData
from pydeseq2.dds import DeseqDataSet, DefaultInference
from pydeseq2.ds import DeseqStats
from statsmodels.stats.multitest import multipletests

import pandas as pd
import numpy as np
import scipy.sparse as sp
from collections import defaultdict

def _to_dense(X):
    return X.toarray() if sp.issparse(X) else np.asarray(X)

def pseudobulk_by(
    adata,
    groupby=["cell_class", "region", "age", "condition"],
    layer=None,
    min_cells=30,
):
    """
    Create pseudobulk counts per unique combination of the groupby factors.

    Returns:
        pb_dict[cell_class][(region, age)] = (counts_df, meta_df)
        - counts_df: genes x pseudobulk-samples
        - meta_df:   rows = pseudobulk-samples (index = pb_key), cols = groupby fields
    """
    X = _to_dense(adata.layers[layer] if layer else adata.X)
    obs = adata.obs.copy()
    genes = adata.var_names

    # Build a stable pseudobulk key
    obs["pb_key"] = obs[groupby].astype(str).agg("__".join, axis=1)

    # Sum counts per pseudobulk key
    summed = (
        pd.DataFrame(X, index=obs["pb_key"], columns=genes)
        .groupby(level=0)
        .sum()
    )

    # Group metadata (indexed by pb_key already)
    meta = obs[groupby + ["pb_key"]].drop_duplicates("pb_key").set_index("pb_key")

    pb = defaultdict(dict)
    for key, row in meta.iterrows():
        cell_class = row["cell_class"]
        reg = row["RBD_compartment_simplified"]
        age = str(row["age"])

        # one pseudobulk column for this (cell_class, region, age, condition)
        counts_col = summed.loc[[key]].T
        counts_col.columns = [key]

        # pull the corresponding one-row metadata frame (indexed by pb_key)
        meta_row = meta.loc[[key]]

        if (reg, age) in pb[cell_class]:
            counts_df_all, meta_df = pb[cell_class][(reg, age)]
            counts_df_all = pd.concat([counts_df_all, counts_col], axis=1)
            meta_df = pd.concat([meta_df, meta_row], axis=0)
        else:
            counts_df_all = counts_col
            meta_df = meta_row

        pb[cell_class][(reg, age)] = (counts_df_all, meta_df)

    # Optional: drop group slots with too few cells total
    # (counts of cells in obs, not number of pseudobulk samples)
    to_prune = []
    for cl in list(pb.keys()):
        for (reg, age), (counts_df, meta_df) in list(pb[cl].items()):
            n_cells = (
                (obs["cell_class"] == cl)
                & (obs["RBD_compartment_simplified"] == reg)
                & (obs["age"].astype(str) == age)
            ).sum()
            if n_cells < min_cells:
                to_prune.append((cl, (reg, age)))
    for cl, key in to_prune:
        del pb[cl][key]

    return pb

# ---------- helpers from earlier ----------
def _align(counts_df, meta_df):
    if not counts_df.columns.equals(meta_df.index):
        common = counts_df.columns.intersection(meta_df.index)
        if len(common)==0:
            raise ValueError("No overlapping sample IDs")
        counts_df = counts_df.loc[:, common]
        meta_df   = meta_df.loc[common]
    return counts_df, meta_df

def _to_anndata(counts_df, meta_df):
    X = counts_df.T.astype(np.int64)
    ad = AnnData(X=X.values)
    ad.obs_names = X.index.astype(str)
    ad.var_names = counts_df.index.astype(str)
    ad.obs = meta_df.copy()
    return ad

def _prep_condition(mdf, condition_col="condition", ref_condition="control"):
    mdf = mdf.copy()
    mdf[condition_col] = mdf[condition_col].astype(str).str.strip()
    uniq = list(pd.unique(mdf[condition_col]))
    if ref_condition not in uniq:
        m = [u for u in uniq if u.lower()==ref_condition.lower()]
        if not m:
            raise ValueError(f"Ref '{ref_condition}' not present. Found: {uniq}")
        ref_condition = m[0]
    mdf[condition_col] = pd.Categorical(mdf[condition_col],
                                        categories=[ref_condition]+[u for u in uniq if u!=ref_condition],
                                        ordered=True)
    return mdf, ref_condition

def simple_effect_at_subset(counts_df, meta_df, mask, condition_col="condition",
                            ref_condition="control", n_cpus=8, label=""):
    """Generic: run mtDSB vs control inside boolean `mask` on meta_df."""
    counts_df, meta_df = _align(counts_df, meta_df)
    if mask.sum() < 2:
        raise ValueError(f"Not enough samples for subset: {label}")
    cdf = counts_df.loc[:, mask]
    mdf = meta_df.loc[mask].copy()
    mdf, ref_condition = _prep_condition(mdf, condition_col=condition_col, ref_condition=ref_condition)

    # detect case level
    cats = mdf[condition_col].cat.categories.tolist()
    case = [c for c in cats if c != ref_condition][0] if len(cats)>1 else None
    if case is None:
        raise ValueError(f"No case level in subset: {label}")

    ad = _to_anndata(cdf, mdf)
    dds = DeseqDataSet(adata=ad, design_factors=[condition_col],
                       refit_cooks=True, inference=DefaultInference(n_cpus=n_cpus))
    dds.deseq2()
    st = DeseqStats(dds, contrast=[condition_col, case, ref_condition], inference=DefaultInference(n_cpus=n_cpus))
    st.summary()
    res = st.results_df.copy()
    # standardize column names
    res = res.rename(columns={
        "log2FoldChange": f"log2FC_{case}_vs_{ref_condition}",
        "lfcSE":          f"lfcSE_{case}_vs_{ref_condition}",
        "pvalue":         f"pvalue_{case}_vs_{ref_condition}",
        "padj":           f"padj_{case}_vs_{ref_condition}",
    })
    res.attrs = {"subset_label": label, "case": case, "ref": ref_condition}
    return res

def combine_delta(res_A, res_B, lfc_prefix="log2FC_", se_prefix="lfcSE_",
                  p_colname="p_delta", q_colname="q_delta"):
    """Generic Δ (B − A) with z-test using SEs; works for ages or regions."""
    lfc_A = [c for c in res_A.columns if c.startswith(lfc_prefix)][0]
    se_A  = [c for c in res_A.columns if c.startswith(se_prefix)][0]
    lfc_B = [c for c in res_B.columns if c.startswith(lfc_prefix)][0]
    se_B  = [c for c in res_B.columns if c.startswith(se_prefix)][0]
    df = pd.DataFrame(index=res_A.index.union(res_B.index))
    df[lfc_A] = res_A[lfc_A]; df[se_A] = res_A[se_A]
    df[lfc_B] = res_B[lfc_B]; df[se_B] = res_B[se_B]
    df = df.dropna(subset=[lfc_A, se_A, lfc_B, se_B]).copy()
    df["delta_log2FC"] = df[lfc_B] - df[lfc_A]
    df["se_delta"] = np.sqrt(np.square(df[se_A]) + np.square(df[se_B]))
    df["se_delta"] = df["se_delta"].replace(0, np.nan)
    df = df.dropna(subset=["se_delta"])
    df["z"] = df["delta_log2FC"] / df["se_delta"]
    from scipy.stats import norm
    df[p_colname] = 2*(1 - norm.cdf(np.abs(df["z"])))
    df[q_colname] = multipletests(df[p_colname], method="fdr_bh")[1]
    return df

# ---------- NEW: region-aware analysis ----------
def run_by_region(pb,                      # dict: ct -> (counts_df, meta_df)
                  condition_col="condition",
                  age_col="age",
                  region_col="region",     # <- your anatomical label column
                  min_per_group=2,         # min samples per condition in a subset
                  n_cpus=8):
    """
    For each cell type & region:
      - per-age simple effect: mtDSB vs control within (region, age)
      - region impact score per age (intensity of effect)
      - region-vs-region ΔLFC at same age (interaction proxy)
    Returns nested dict: results[ct][region] -> {...}
    """
    results = {}
    for ct, (counts_df, meta_df) in pb.items():
        counts_df, meta_df = _align(counts_df, meta_df)
        # normalize fields
        md = meta_df.copy()
        md[age_col] = md[age_col].astype(str).str.strip()
        md[region_col] = md[region_col].astype(str).str.strip()
        md[condition_col] = md[condition_col].astype(str).str.strip()

        regions = list(pd.unique(md[region_col]))
        ages    = list(pd.unique(md[age_col]))
        ct_out = {}
        for reg in regions:
            reg_out = {"simple_effects": {}, "impact": {}, "age_list": ages, "region": reg}
            for a in ages:
                mask = (md[region_col]==reg) & (md[age_col]==a)
                # require at least min_per_group per condition
                if mask.sum() >= 2 and all((md.loc[mask, condition_col].value_counts() >= min_per_group).reindex(md[condition_col].unique(), fill_value=0)[:2]):
                    try:
                        res = simple_effect_at_subset(
                            counts_df, md, mask,
                            condition_col=condition_col, ref_condition="control",
                            n_cpus=n_cpus, label=f"{reg}@{a}"
                        )
                        reg_out["simple_effects"][a] = res
                        # impact score: count of sig genes & sum -log10(FDR) among sig
                        lfc_col = [c for c in res.columns if c.startswith("log2FC_")][0]
                        padj_col = [c for c in res.columns if c.startswith("padj_")][0]
                        sig = (res[padj_col] < 0.05) & (res[lfc_col].abs() >= 1.0)
                        reg_out["impact"][a] = {
                            "n_sig": int(sig.sum()),
                            "sum_neglog10q": float((-np.log10(res.loc[sig, padj_col].clip(lower=1e-300))).sum())
                        }
                    except Exception as e:
                        reg_out["simple_effects"][a] = None
                        reg_out["impact"][a] = {"n_sig": 0, "sum_neglog10q": 0.0}
                else:
                    reg_out["simple_effects"][a] = None
                    reg_out["impact"][a] = {"n_sig": 0, "sum_neglog10q": 0.0}
            ct_out[reg] = reg_out

        # region-vs-region ΔLFC at same age (interaction proxy)
        # For each age, compare each pair of regions’ mtDSB effects
        for regA, regB in [(r1, r2) for i, r1 in enumerate(regions) for r2 in regions[i+1:]]:
            for a in ages:
                resA = ct_out[regA]["simple_effects"].get(a)
                resB = ct_out[regB]["simple_effects"].get(a)
                if resA is None or resB is None:
                    continue
                dtab = combine_delta(resA, resB, p_colname="p_delta_region", q_colname="q_delta_region")
                if "deltas_region" not in ct_out[regA]:
                    ct_out[regA]["deltas_region"] = {}
                if "deltas_region" not in ct_out[regB]:
                    ct_out[regB]["deltas_region"] = {}
                ct_out[regA]["deltas_region"][(a, regB)] = dtab  # Δ = regB − regA at age a
                # also store the opposite direction (flip sign)
                dtab_flip = dtab.copy()
                dtab_flip["delta_log2FC"] = -dtab_flip["delta_log2FC"]
                dtab_flip["z"] = -dtab_flip["z"]
                ct_out[regB]["deltas_region"][(a, regA)] = dtab_flip

        results[ct] = ct_out
    return results

# ---------- ranking regions by "affectedness" ----------
def summarize_region_impact(results_by_region, ct, metric="sum_neglog10q", age=None):
    """
    Build a ranking table of regions by impact metric.
    metric: "n_sig" or "sum_neglog10q"
    age: specific age or None to aggregate across ages (sum).
    """
    rows = []
    for reg, d in results_by_region[ct].items():
        imp = d.get("impact", {})
        if age is None:
            vals = [imp[a][metric] for a in imp.keys()]
            score = np.nansum(vals)
        else:
            score = imp.get(str(age), {}).get(metric, 0.0) if str(age) in imp else 0.0
        rows.append({"region": reg, "score": score})
    out = pd.DataFrame(rows).sort_values("score", ascending=False).reset_index(drop=True)
    return out

from pydeseq2.dds import DeseqDataSet, DefaultInference
from pydeseq2.ds import DeseqStats
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from scipy.stats import norm

def _align(counts_df, meta_df):
    common = counts_df.columns.intersection(meta_df.index)
    return counts_df.loc[:, common], meta_df.loc[common]

def simple_effect_at_subset(counts_df, meta_df, mask, condition_col="condition",
                            ref_condition="control", n_cpus=8, label=""):
    counts_df, meta_df = _align(counts_df, meta_df)
    if mask.sum() < 2:
        raise ValueError(f"Not enough samples for subset {label}")
    cdf = counts_df.loc[:, mask]
    mdf = meta_df.loc[mask].copy()
    mdf[condition_col] = pd.Categorical(
        mdf[condition_col].astype(str),
        categories=[ref_condition] + [x for x in mdf[condition_col].unique() if x != ref_condition],
        ordered=True
    )
    adata = cdf.T
    inference = DefaultInference(n_cpus=n_cpus)
    dds = DeseqDataSet(adata=adata, clinical=mdf, design_factors=[condition_col],
                       refit_cooks=True, inference=inference)
    dds.deseq2()
    st = DeseqStats(dds, contrast=[condition_col, "mtDSB", ref_condition], inference=inference)
    st.summary()
    return st.results_df.copy()

def run_by_region(pb_class, condition_col="condition", age_col="age",
                  region_col="region", min_per_group=2, n_cpus=8):
    results = {}
    for cl, (counts_df, meta_df) in pb_class.items():
        md = meta_df.copy()
        md[age_col] = md[age_col].astype(str)
        md[region_col] = md[region_col].astype(str)
        md[condition_col] = md[condition_col].astype(str)
        ages = md[age_col].unique()
        regions = md[region_col].unique()

        cl_results = {}
        for reg in regions:
            reg_out = {"simple_effects": {}, "impact": {}}
            for age in ages:
                mask = (md[region_col]==reg) & (md[age_col]==age)
                if mask.sum() < 2: continue
                try:
                    res = simple_effect_at_subset(
                        counts_df, md, mask,
                        condition_col=condition_col,
                        ref_condition="control",
                        n_cpus=n_cpus,
                        label=f"{cl}_{reg}_{age}"
                    )
                    lfc_col = [c for c in res.columns if "log2FoldChange" in c][0]
                    padj_col = [c for c in res.columns if "padj" in c][0]
                    sig = (res[padj_col] < 0.05) & (res[lfc_col].abs() > 1)
                    reg_out["simple_effects"][age] = res
                    reg_out["impact"][age] = {
                        "n_sig": int(sig.sum()),
                        "sum_neglog10q": float((-np.log10(res.loc[sig, padj_col].clip(lower=1e-300))).sum())
                    }
                except Exception:
                    pass
            cl_results[reg] = reg_out
        results[cl] = cl_results
    return results

def summarize_region_impact(results_by_region, cell_class, metric="sum_neglog10q", age=None):
    rows = []
    for reg, d in results_by_region[cell_class].items():
        imp = d["impact"]
        if not imp: continue
        if age is None:
            vals = [v[metric] for v in imp.values()]
            score = np.nansum(vals)
        else:
            score = imp.get(str(age), {}).get(metric, 0.0)
        rows.append({"region": reg, "score": score})
    out = pd.DataFrame(rows).sort_values("score", ascending=False).reset_index(drop=True)
    return out
import numpy as np
import pandas as pd
from anndata import AnnData
from pydeseq2.dds import DeseqDataSet, DefaultInference
from pydeseq2.ds import DeseqStats

def _as_anndata(counts_df, meta_df):
    ad = AnnData(counts_df.T)
    ad.obs = meta_df.copy()
    ad.var_names = counts_df.index
    ad.obs_names = counts_df.columns
    return ad

def run_by_region(
    pb_class,
    condition_col="condition",
    age_col="age",
    region_col="region",
    ref_condition="control",
    min_per_condition=2,
    n_cpus=8,
    verbose=True
):
    results = {}
    for cl, (counts_df, meta_df) in pb_class.items():
        counts_df, meta_df = counts_df.copy(), meta_df.copy()
        md = meta_df.copy()
        md[age_col] = md[age_col].astype(str)
        md[region_col] = md[region_col].astype(str)
        md[condition_col] = md[condition_col].astype(str)
        ages = sorted(md[age_col].unique(), key=lambda x: (len(x), x))
        regions = sorted(md[region_col].unique())
        inference = DefaultInference(n_cpus=n_cpus)

        cl_out = {}
        for reg in regions:
            reg_out = {"simple_effects": {}, "impact": {}}
            for age in ages:
                mask = (md[region_col]==reg) & (md[age_col]==age)
                sub = md.loc[mask]
                if sub.empty: 
                    continue
                vc = sub[condition_col].value_counts()
                if len(vc) < 2 or any(vc < min_per_condition):
                    if verbose: print(f"[{cl}] skip {reg}@{age} (counts {vc.to_dict()})")
                    continue
                case = [x for x in vc.index if x != ref_condition][0]
                if verbose: print(f"[{cl}] running {reg}@{age} ({vc.to_dict()})")

                ad = _as_anndata(counts_df.loc[:, sub.index], sub)
                dds = DeseqDataSet(
                    adata=ad,
                    design_factors=[condition_col],
                    refit_cooks=True,
                    inference=inference
                )
                dds.deseq2()
                st = DeseqStats(dds, contrast=[condition_col, case, ref_condition], inference=inference)
                st.summary()
                res = st.results_df.copy()

                lfc_col = "log2FoldChange"
                padj_col = "padj"
                sig = (res[padj_col] < 0.05) & (res[lfc_col].abs() >= 1.0)
                reg_out["simple_effects"][age] = res
                reg_out["impact"][age] = {
                    "n_sig": int(sig.sum()),
                    "sum_neglog10q": float((-np.log10(res.loc[sig, padj_col].clip(lower=1e-300))).sum())
                }

            cl_out[reg] = reg_out
        results[cl] = cl_out
    return results
import numpy as np
import matplotlib.pyplot as plt

def volcano_clean(
    res_df,
    lfc_col="log2FoldChange",
    padj_col="padj",
    title="mtDSB vs control",
    lfc_thr=1.0,
    q=0.05,
    label_top_fdr=10,           # top genes by FDR to label
    label_top_abs_lfc=10,       # top by |LFC|
    max_labels=20,              # absolute label cap
    highlight_genes=None,       # list of genes to force label
    figsize=(9,7),
    fontsize=13
):
    """Pretty volcano plot for DESeq2 / PyDESeq2 results."""

    df = res_df.copy().replace([np.inf, -np.inf], np.nan).dropna(subset=[lfc_col, padj_col])
    df["neglog10q"] = -np.log10(df[padj_col].clip(lower=1e-300))
    sig = (df[padj_col] < q) & (df[lfc_col].abs() >= lfc_thr)

    # pick label set
    idx_fdr = df.sort_values(padj_col).head(label_top_fdr).index
    idx_lfc = df.reindex(df[lfc_col].abs().sort_values(ascending=False).head(label_top_abs_lfc).index).index
    pick = set(idx_fdr).union(set(idx_lfc))
    if highlight_genes:
        pick.update(set(df.index.intersection(highlight_genes)))
    if len(pick) > max_labels:
        pick = set(df.loc[list(pick)].sort_values("neglog10q", ascending=False).head(max_labels).index)

    # figure
    plt.figure(figsize=figsize, dpi=150)
    plt.scatter(df[lfc_col], df["neglog10q"], s=16, alpha=0.35, color="lightgrey", label="not sig")
    up = df[sig & (df[lfc_col] > 0)]
    dn = df[sig & (df[lfc_col] < 0)]
    plt.scatter(up[lfc_col], up["neglog10q"], s=28, color="#d73027", label="Up (sig)")
    plt.scatter(dn[lfc_col], dn["neglog10q"], s=28, color="#4575b4", label="Down (sig)")

    # thresholds
    plt.axvline(+lfc_thr, ls="--", color="k", lw=1)
    plt.axvline(-lfc_thr, ls="--", color="k", lw=1)
    plt.axhline(-np.log10(q), ls="--", color="k", lw=1)

    plt.xlabel("log₂ fold change (mtDSB vs control)", fontsize=fontsize)
    plt.ylabel("-log₁₀(FDR)", fontsize=fontsize)
    plt.title(title, fontsize=fontsize+2, pad=10)

    # label genes
    texts = []
    try:
        from adjustText import adjust_text
        for g in pick:
            x, y = df.at[g, lfc_col], df.at[g, "neglog10q"]
            t = plt.text(x, y, g, fontsize=fontsize-3, va="bottom", ha="center")
            texts.append(t)
            plt.scatter([x],[y], s=40, edgecolor="black", linewidth=0.6, zorder=4)
        adjust_text(
            texts,
            only_move={'points':'y', 'text':'xy'},
            expand_points=(1.2, 1.4),
            expand_text=(1.1, 1.2),
            arrowprops=dict(arrowstyle="-", lw=0.8, color="0.25"),
        )
    except ImportError:
        for g in pick:
            x, y = df.at[g, lfc_col], df.at[g, "neglog10q"]
            plt.annotate(
                g, xy=(x, y), xytext=(x+0.2*np.sign(x), y+0.3),
                fontsize=fontsize-3,
                arrowprops=dict(arrowstyle="-", lw=0.8, color="0.25"),
                ha="center", va="bottom"
            )
            plt.scatter([x],[y], s=40, edgecolor="black", linewidth=0.6, zorder=4)

    plt.legend(frameon=False, fontsize=fontsize-3)
    plt.tight_layout()
    plt.show()


# In[3]:


adata = sc.read_h5ad('../data/rbd_annotated_monod_annotated.h5ad')
adata.obs_names_make_unique()


# In[4]:


adata_raw = sc.read_h5ad('../data/mtDNA_DSB_5k_clustered_manual_annotation_with_raw.h5ad')
adata_raw.obs_names_make_unique()


# In[5]:


adata_anno_new = sc.read_h5ad('../data/mtDNA_DSB_5k_clustered_higher_res_LLM_anno.h5ad')
adata_anno_new.obs_names_make_unique()


# In[6]:


cols = ['leiden_0.5', 'leiden_1', 'leiden_2', 'leiden_1.5',
        'leiden_2.1', 'leiden_2.3', 'leiden_2.5', 'leiden_3', 'cell_class']

for col in cols:
    map_dict = dict(zip(adata_anno_new.obs.index, adata_anno_new.obs[col]))
    adata.obs[col] = adata.obs.index.map(map_dict)


# In[393]:


adata[adata.obs.cell_class == 'Neural Stem Cells'].obs.condition.value_counts()


# In[392]:


adata[adata.obs.cell_class == 'Neural Stem Cells'].obs.age.value_counts()


# In[281]:


import numpy as np
import pandas as pd
import scanpy as sc

def compute_gene_qc(adata):
    """
    Returns a DataFrame with per-gene stats:
      n_cells_detected
      detect_rate
      total_counts
    Works for dense or sparse X.
    """
    X = adata.X

    # number of cells where gene is expressed (>0)
    if hasattr(X, "toarray"):  # sparse
        n_cells_detected = np.asarray((X > 0).sum(axis=0)).ravel()
        total_counts = np.asarray(X.sum(axis=0)).ravel()
    else:  # dense
        n_cells_detected = (X > 0).sum(axis=0)
        total_counts = X.sum(axis=0)

    detect_rate = n_cells_detected / adata.n_obs

    qc_df = pd.DataFrame({
        "gene": adata.var_names,
        "n_cells_detected": n_cells_detected,
        "detect_rate": detect_rate,
        "total_counts": total_counts,
    }).set_index("gene")

    return qc_df

qc_df = compute_gene_qc(adata)
qc_df.head()


# In[286]:


qc_df.sort_values(by = 'total_counts', ascending = False).loc['Nr5a1']


# In[290]:


# tune thresholds here
min_detect_rate = 0.01        # keep genes seen in ≥0.5% of cells
min_cells =1000                 # or seen in at least 30 cells total
min_total_counts = 2000          # or at least 50 UMIs total across all cells

informative_mask = (
    (qc_df["detect_rate"] >= min_detect_rate) |
    (qc_df["n_cells_detected"] >= min_cells) |
    (qc_df["total_counts"] >= min_total_counts) 
)

genes_to_keep = qc_df.index[informative_mask].tolist()
genes_to_drop = qc_df.index[~informative_mask].tolist()

print(f"Keeping {len(genes_to_keep)} / {qc_df.shape[0]} genes")
print("Example drops:", genes_to_drop[:20])


# In[296]:


adata = adata[:, genes_to_keep].copy()
adata


# In[7]:


adata.layers['counts'] = adata_raw.layers['counts'].copy()


# In[162]:


adata.write('../data/mtDNA_DSB_5k_clustered_annotation_with_rbd.h5ad')


# In[8]:


sc.pl.umap(adata, color = 'RBD_compartment_simplified')


# In[12]:


sc.pl.umap(adata, color = 'cell_class')


# In[11]:


del adata.uns['cell_class_colors']


# In[13]:


for domain in ['cell_class']:
    if domain in adata.obs.columns:
        print(f"Plotting {domain} …")
        plot_spatial_compact_fast(
            adata,
            color=domain,
            groupby="sample_id",
            spot_size=0.3,
            cols=6,
            height=8,
            legend_col_width=1.0,
            rasterized=True,   # big speedup for large data
            dpi=100             # lower = faster preview
        )
    else:
        print(f"Skipping {domain} — not in adata.obs")


# In[14]:


for domain in ['RBD_compartment_simplified']:
    if domain in adata.obs.columns:
        print(f"Plotting {domain} …")
        plot_spatial_compact_fast(
            adata,
            color=domain,
            groupby="sample_id",
            spot_size=0.3,
            cols=6,
            height=8,
            legend_col_width=1.0,
            rasterized=True,   # big speedup for large data
            dpi=100             # lower = faster preview
        )
    else:
        print(f"Skipping {domain} — not in adata.obs")


# In[15]:


adata.X = adata.layers['counts']
adata.X = adata.X.astype(int)


# # pseudobulking

# In[297]:


pb = pseudobulk_by(
    adata,
    groupby=["cell_class", "RBD_compartment_simplified", "age", "condition", "sample_id"],  # 👈 add sample_id
    layer="counts" if "counts" in adata.layers else None,
    min_cells=30,
)


# In[412]:


import pandas as pd

# Get expression matrix as DataFrame (cells × genes)
df_expr = adata.to_df()
df_expr['cell_class'] = adata.obs['cell_class'].values

# 1️⃣ Fraction of cells expressing each gene (>0 counts)
expr_frac = (
    (df_expr.drop(columns='cell_class') > 0)
    .assign(cell_class=df_expr['cell_class'])
    .groupby('cell_class')
    .mean()
)

# 2️⃣ Mean expression across all cells (including zeros)
expr_mean = (
    df_expr
    .groupby('cell_class')
    .mean(numeric_only=True)
)

# 3️⃣ Combine into a single tidy DataFrame
expr_summary = (
    expr_mean
    .stack()
    .to_frame('mean_expr')
    .join(expr_frac.stack().to_frame('frac_expressing'))
    .reset_index()
    .rename(columns={'level_1': 'gene'})
)

expr_summary.head()


# In[415]:


expr_summary[expr_summary.cell_class == 'Oligodendrocytes']


# In[419]:


import matplotlib.pyplot as plt
import seaborn as sns

def plot_expression_fraction_vs_mean(expr_summary, cell_class=None, top_n=15):
    """
    Visualize expression prevalence (fraction of cells) vs mean expression per gene.
    """
    if cell_class:
        df = expr_summary.query("cell_class == @cell_class").copy()
        title = f"{cell_class}: Fraction vs Mean Expression"
    else:
        df = expr_summary.copy()
        title = "Fraction vs Mean Expression (all cell classes)"

    plt.figure(figsize=(7, 6))
    sns.scatterplot(
        data=df,
        x='frac_expressing',
        y='mean_expr',
        s=20,
        alpha=0.6,
        edgecolor='none'
    )

    # highlight top genes
    top = df.nlargest(top_n, 'mean_expr')
    for _, r in top.iterrows():
        plt.text(
            r['frac_expressing'],
            r['mean_expr'],
            r['gene'],
            fontsize=8,
            color='black',
            ha='left',
            va='bottom'
        )

    plt.xlabel("Fraction of cells expressing gene")
    plt.ylabel("Mean expression (counts or normalized)")
    plt.title(title)
    plt.xscale('log')
    plt.yscale('log')
    plt.grid(True, linestyle='--', alpha=0.3)
    plt.tight_layout()
    plt.axvline(0.05, color='red', linestyle='--', label='min fraction (5%)')
    plt.axhline(0.5, color='blue', linestyle='--', label='min mean expr (0.5)')
    plt.legend()
    plt.show()

# Example usage:
plot_expression_fraction_vs_mean(expr_summary, cell_class="Oligodendrocytes")


# In[420]:


plot_expression_fraction_vs_mean(expr_summary, cell_class="Astrocytes")


# In[423]:


# heuristic combined score: geometric mean of fraction and mean expression
expr_summary['expr_score'] = (expr_summary['frac_expressing'] * expr_summary['mean_expr']) ** 0.5


# In[426]:


expr_summary_sub = expr_summary[expr_summary.cell_class == 'Oligodendrocytes']


# In[436]:


expr_summary_sub[expr_summary_sub.gene.isin(['Atf4','Atf5','Pmch','Nmu',])]


# In[439]:


import numpy as np

def pick_celltype_gene_universe(expr_summary, percentile=50):
    """
    Returns: dict[cell_class] = list of genes to keep for that cell_class
    Uses expr_score percentile per class as cutoff.
    """
    keep = {}
    for cc, sub in expr_summary.groupby("cell_class"):
        cutoff = np.percentile(sub['expr_score'], percentile)
        gene_list = sub.loc[sub['expr_score'] >= cutoff, 'gene'].unique().tolist()
        keep[cc] = gene_list
    return keep

gene_universe = pick_celltype_gene_universe(expr_summary, percentile=50)


# In[440]:


gene_universe


# In[438]:


import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def plot_expr_score_distribution(expr_summary, cell_class, bins=50):
    # subset for that class
    df_cc = expr_summary.query("cell_class == @cell_class").copy()

    # guard in case it's empty
    if df_cc.empty:
        print(f"No data for {cell_class}")
        return

    # basic histogram / KDE
    plt.figure(figsize=(6,4))
    sns.histplot(
        df_cc['expr_score'],
        bins=bins,
        kde=True,
        edgecolor='none',
        alpha=0.7
    )

    plt.xlabel("expr_score = sqrt(frac_expressing * mean_expr)")
    plt.ylabel("Number of genes")
    plt.title(f"{cell_class}: expr_score distribution")

    # optional: suggest a heuristic cutoff
    cutoff = np.percentile(df_cc['expr_score'], 75)  # top quartile, for example
    plt.axvline(cutoff, color='red', linestyle='--', linewidth=1)
    plt.text(
        cutoff,
        plt.ylim()[1]*0.9,
        f"75th pct = {cutoff:.3f}",
        color='red',
        rotation=90,
        va='top',
        ha='right',
        fontsize=8
    )

    plt.tight_layout()
    plt.show()

    return df_cc, cutoff

# example:
df_ol, cutoff_ol = plot_expr_score_distribution(expr_summary, "Astrocytes")


# In[298]:


pb_class = {}
for cl, regdict in pb.items():
    all_counts = []
    all_meta = []
    for (reg, age), (counts_df, meta_df) in regdict.items():
        all_counts.append(counts_df)
        all_meta.append(meta_df)
    counts_df_all = pd.concat(all_counts, axis=1)
    meta_df_all = pd.concat(all_meta)
    pb_class[cl] = (counts_df_all, meta_df_all)


# ## run DGE

# In[299]:


results_by_region = run_by_region(
    pb_class,
    condition_col="condition",
    age_col="age",
    region_col="RBD_compartment_simplified",
    ref_condition="control",
    min_per_condition=2,
    n_cpus=8,
    verbose=True
)


# In[300]:


import numpy as np
import pandas as pd

def summarize_region_impact(
    results_by_region, cell_class, metric="sum_neglog10q",
    age=None, agg="sum", fill_value=0.0
):
    """
    Aggregate region impact for a cell class.
    - metric: key inside d['impact'][age], e.g. 'sum_neglog10q'
    - age: None -> aggregate across ages; else str/int for a single age
    - agg: 'sum' | 'mean' | 'max' for combining ages when age=None
    """
    rows = []
    items = results_by_region.get(cell_class, {})
    for reg, d in items.items():
        imp = d.get("impact", {})
        if age is None:
            vals = [v.get(metric, np.nan) for v in imp.values() if isinstance(v, dict)]
            vals = np.array(vals, dtype=float)
            if vals.size == 0:
                score = fill_value
            else:
                if agg == "mean":
                    score = np.nanmean(vals)
                elif agg == "max":
                    score = np.nanmax(vals)
                else:  # 'sum'
                    score = np.nansum(vals)
        else:
            score = imp.get(str(age), {}).get(metric, fill_value)
        rows.append({"region": reg, "score": float(score)})
    return (pd.DataFrame(rows)
              .sort_values("score", ascending=False)
              .reset_index(drop=True))


# In[301]:


region_dict = {}
for cell in adata.obs.cell_class.unique():
    print(cell)
    rank_regions = summarize_region_impact(results_by_region,cell)
    print(rank_regions.head(10))
    region_dict[cell] = rank_regions


# In[302]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Pick a score column automatically if user doesn’t pass one
def _pick_score_col(df, preferred=("score","n_sig_genes","n_genes","n_regions",
                                   "sum_neglog10FDR","mean_neglog10FDR","mean_LFC")):
    num_cols = [c for c in df.columns if np.issubdtype(df[c].dtype, np.number)]
    for c in preferred:
        if c in num_cols: 
            return c
    # fallback: first numeric column
    return num_cols[0] if num_cols else None

def plot_region_impact_for_glia(
    results_by_region,
    glia_classes,
    summarize_region_impact,     # your function(region_summary = summarize_region_impact(...))
    top_k=10,
    score_col=None,
    figsize_bar=(8,4),
    figsize_heatmap=(10,0.6),
    cmap="viridis",
    save_prefix=None
):
    """
    For each glia class:
      1) calls summarize_region_impact(results_by_region, glia)
      2) plots top-K regions by selected score (barh)
    Also builds a glia x region heatmap of the selected score.
    """
    all_rows = []
    per_glia_tables = {}

    # --- per-glia bars
    for glia in glia_classes:
        df = summarize_region_impact(results_by_region, glia)
        if df is None or len(df)==0:
            print(f"[skip] {glia}: empty result")
            continue

        # ensure region is index
        if "region" in df.columns:
            df = df.set_index("region")

        sc = score_col or _pick_score_col(df)
        if sc is None:
            print(f"[skip] {glia}: no numeric score column found")
            continue

        # store for heatmap
        per_glia_tables[glia] = df.copy()
        tmp = df[[sc]].copy()
        tmp["glia"] = glia
        tmp["region"] = tmp.index
        all_rows.append(tmp.reset_index(drop=True))

        # plot top-K bar
        top = df.sort_values(sc, ascending=False).head(top_k).iloc[::-1]
        plt.figure(figsize=figsize_bar)
        sns.barplot(x=top[sc], y=top.index, color="#5B8DEF")
        plt.title(f"{glia} — top {top_k} regions ({sc})")
        plt.xlabel(sc); plt.ylabel("region")
        plt.tight_layout()
        if save_prefix:
            plt.savefig(f"{save_prefix}_{glia}_top_regions.png", dpi=300, bbox_inches="tight")
        plt.show()

    if not all_rows:
        print("No data to build heatmap.")
        return

    # --- glia x region heatmap
    long = pd.concat(all_rows, ignore_index=True)
    sc = score_col or long.columns[0]  # already picked above
    mat = long.pivot_table(index="glia", columns="region", values=sc, aggfunc="mean")

    # z-score per glia (optional but helps comparability)
    mat_z = (mat - mat.mean(axis=1, keepdims=True)) / mat.std(axis=1, keepdims=True)

    h = max(3, figsize_heatmap[1]*mat.shape[0])  # scale height by #glia
    plt.figure(figsize=(figsize_heatmap[0], h))
    sns.heatmap(mat_z, cmap=cmap, center=0, linewidths=.4, linecolor="white",
                cbar_kws={"label": f"{sc} (row z-score)"})
    plt.title(f"Region impact heatmap across glia (score={sc})")
    plt.xlabel("region"); plt.ylabel("glia class")
    plt.tight_layout()
    if save_prefix:
        plt.savefig(f"{save_prefix}_heatmap.png", dpi=300, bbox_inches="tight")
    plt.show()

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def plot_integrated_region_scores(region_scores_dict, cmap="mako", normalize=True, save=None):
    """
    region_scores_dict: dict {glia_class_name: df with columns ['region','score']}
    Produces integrated plots across glial classes.
    """
    # --- combine all into one long DataFrame ---
    df_all = []
    for glia, df in region_scores_dict.items():
        tmp = df.copy()
        tmp["glia"] = glia
        df_all.append(tmp)
    df_all = pd.concat(df_all, ignore_index=True)

    # --- pivot into glia x region matrix ---
    mat = df_all.pivot_table(index="glia", columns="region", values="score", aggfunc="mean", fill_value=0)

    # optional normalization (row-wise)
    if normalize:
        mat = mat.div(mat.max(axis=1), axis=0)

    # --- plot heatmap ---
    plt.figure(figsize=(1*mat.shape[1], 0.8*mat.shape[0]+2))
    sns.heatmap(mat, cmap=cmap, annot=False, linewidths=0.4, linecolor="white", cbar_kws={"label": "Relative score"})
    plt.title("Regional impact across glial classes")
    plt.xlabel("Region")
    plt.ylabel("Glial cell class")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    if save:
        plt.savefig(f"{save}_heatmap.png", dpi=300, bbox_inches="tight")
    plt.show()

    # --- summarize top regions across all glia ---
    region_mean = df_all.groupby("region")["score"].mean().sort_values(ascending=False)
    plt.figure(figsize=(6,4))
    sns.barplot(x=region_mean.values[:10], y=region_mean.index[:10], palette="crest")
    plt.title("Top 10 regions (mean score across glia)")
    plt.xlabel("Mean score"); plt.ylabel("")
    plt.tight_layout()
    if save:
        plt.savefig(f"{save}_region_bar.png", dpi=300, bbox_inches="tight")
    plt.show()

    # --- per-glia line plot (optional) ---
    plt.figure(figsize=(9,5))
    for glia in mat.index:
        plt.plot(mat.columns, mat.loc[glia], marker="o", label=glia)
    plt.legend(bbox_to_anchor=(1.05,1), loc="upper left")
    plt.title("Regional profiles per glia (normalized)")
    plt.ylabel("Relative score")
    plt.xlabel("Region")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    if save:
        plt.savefig(f"{save}_profiles.png", dpi=300, bbox_inches="tight")
    plt.show()

    return mat, df_all


# In[303]:


mat, df_all = plot_integrated_region_scores(
    region_dict,
    cmap="rocket_r",
    normalize=False,

    #save="figs/glia_region_impact"
)


# In[304]:


import numpy as np
import matplotlib.pyplot as plt
from adjustText import adjust_text

def volcano_pretty(
    df,
    title="",
    alpha=0.05,
    fc=1.0,
    max_labels=20,
    s_range=(10, 120),
    cmap_up="Reds",
    cmap_down="Blues",
    dpi=140,
    ylim=None,          # ✅ new argument
    xlim=None,          # also handy if you want consistency
):
    """
    Publication-quality volcano:
      • point size ~ baseMean (log10 scaled)
      • color by up/down direction
      • annotated top genes
      • optional axis limits (ylim, xlim)
    Expects columns: log2FoldChange, padj, pvalue, baseMean, and gene/index.
    """
    # --- extract arrays ---
    df = df.copy()
    x = df["log2FoldChange"].to_numpy()
    y = -np.log10(np.clip(df["pvalue"].to_numpy(), 1e-300, 1.0))
    bm = np.log10(df["baseMean"].clip(lower=1) + 1.0).to_numpy()

    # --- size scaling ---
    bm_range = np.ptp(bm) if np.ptp(bm) > 0 else 1
    bm_norm = (bm - bm.min()) / bm_range
    sizes = s_range[0] + bm_norm * (s_range[1] - s_range[0])

    # --- significance ---
    sig = (df["padj"] < alpha) & (np.abs(x) >= fc)
    up = sig & (x > 0)
    down = sig & (x < 0)

    # --- colormaps ---
    cu = plt.cm.get_cmap(cmap_up)
    cd = plt.cm.get_cmap(cmap_down)

    # --- plot ---
    plt.ioff()
    fig, ax = plt.subplots(figsize=(6, 5), dpi=dpi)

    # background (non-sig)
    ax.scatter(
        x[~sig], y[~sig], s=sizes[~sig], color="lightgrey",
        alpha=0.5, linewidths=0, rasterized=True
    )
    # significant up/down
    ax.scatter(
        x[up], y[up], s=sizes[up], color=cu(0.8),
        alpha=0.9, edgecolors="none", rasterized=True
    )
    ax.scatter(
        x[down], y[down], s=sizes[down], color=cd(0.8),
        alpha=0.9, edgecolors="none", rasterized=True
    )

    # --- guides ---
    ax.axvline(fc, ls="--", lw=0.8, c="grey")
    ax.axvline(-fc, ls="--", lw=0.8, c="grey")
    ax.axhline(-np.log10(alpha), ls="--", lw=0.8, c="grey")

    # --- limits (✅ new) ---
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)

    # --- labeling ---
    gene_names = (
        df["gene"].astype(str).to_numpy()
        if "gene" in df.columns
        else df.index.astype(str).to_numpy()
    )
    score = -np.log10(df["padj"].clip(1e-300)).to_numpy() * np.abs(x)
    up_idx = np.where(up)[0]
    down_idx = np.where(down)[0]
    top_up = up_idx[np.argsort(score[up_idx])[::-1][:max_labels // 2]]
    top_down = down_idx[np.argsort(score[down_idx])[::-1][:max_labels // 2]]
    label_idx = np.concatenate([top_up, top_down])

    texts = []
    for i in label_idx:
        texts.append(
            ax.text(
                x[i], y[i], gene_names[i],
                fontsize=8, ha="center", va="bottom", color="black"
            )
        )

    adjust_text(
        texts, ax=ax,
        expand=(1.05, 1.2),
        arrowprops=dict(arrowstyle="-", lw=0.5, color="grey")
    )

    # --- aesthetics ---
    ax.set_xlabel("log2 fold change", fontsize=11)
    ax.set_ylabel("-log10(p-value)", fontsize=11)
    ax.set_title(title, fontsize=12, pad=10)
    ax.spines[['top','right']].set_visible(False)
    plt.tight_layout()
    plt.ion()
    plt.show()

    return fig, ax


# In[305]:


import numpy as np
import pandas as pd

def summarize_results(
    results,
    table="simple_effects",
    alpha=0.05,
    lfc=1.0,               # abs log2FC threshold
    require_cols=("log2FoldChange","pvalue","padj","baseMean"),
):
    """
    Walk results[cell_type][region][table][age] and compute per-slice summary stats.
    Returns a tidy DataFrame with one row per (cell_type, region, age).
    """
    rows = []
    for ct, reg_dict in results.items():
        if not isinstance(reg_dict, dict): 
            continue
        for reg, tab_dict in reg_dict.items():
            if not isinstance(tab_dict, dict) or table not in tab_dict: 
                continue
            age_dict = tab_dict[table]
            if not isinstance(age_dict, dict): 
                continue
            for age, df in age_dict.items():
                if not isinstance(df, pd.DataFrame) or df.empty:
                    continue
                # sanity columns
                if any(c not in df.columns for c in require_cols):
                    # skip or relax here if needed
                    missing = [c for c in require_cols if c not in df.columns]
                    # you can `print(f"Skip {ct}/{reg}/{age}: missing {missing}")`
                    continue

                x = df["log2FoldChange"].to_numpy()
                p = np.clip(df["pvalue"].to_numpy(), 1e-300, 1.0)
                padj = np.clip(df["padj"].to_numpy(), 1e-300, 1.0)
                bm = df["baseMean"].to_numpy()

                sig = (padj < alpha) & (np.abs(x) >= lfc)
                n = df.shape[0]
                n_sig = int(sig.sum())
                n_up  = int((sig & (x > 0)).sum())
                n_dn  = int((sig & (x < 0)).sum())
                frac_sig = n_sig / n if n else 0.0

                # robust stats (handle slices with no sig)
                if n_sig > 0:
                    abs_lfc_sig = np.abs(x[sig])
                    med_abs_lfc = float(np.median(abs_lfc_sig))
                    mean_abs_lfc = float(abs_lfc_sig.mean())
                    max_abs_lfc = float(abs_lfc_sig.max())
                    # expression-weighted magnitude (optional)
                    w = bm[sig]
                    w = w / w.sum() if w.sum() > 0 else np.full_like(w, 1/len(w))
                    w_mean_abs_lfc = float((np.abs(x[sig]) * w).sum())
                    min_padj = float(padj[sig].min())
                else:
                    med_abs_lfc = mean_abs_lfc = max_abs_lfc = w_mean_abs_lfc = 0.0
                    min_padj = 1.0

                rows.append(dict(
                    cell_type=ct, region=reg, age=str(age),
                    n_genes=n, n_sig=n_sig, n_up=n_up, n_down=n_dn, frac_sig=frac_sig,
                    median_abs_lfc=med_abs_lfc, mean_abs_lfc=mean_abs_lfc,
                    max_abs_lfc=max_abs_lfc, w_mean_abs_lfc=w_mean_abs_lfc,
                    best_padj=min_padj
                ))
    return pd.DataFrame(rows)


# In[306]:


# --- deps & vector-friendly output ---
from IPython.display import display, HTML
import pandas as pd, numpy as np, os
import matplotlib as mpl, matplotlib.pyplot as plt
from adjustText import adjust_text

mpl.rcParams['svg.fonttype'] = 'none'   # text stays as text in SVG
mpl.rcParams['pdf.fonttype'] = 42       # embed editable fonts in PDF
mpl.rcParams['text.usetex'] = False

# --------------------------------------------
# Output folder
# --------------------------------------------
outdir = "../results/figures"
os.makedirs(outdir, exist_ok=True)

# --------------------------------------------
# Helpers to access your results structure
# --------------------------------------------
def _get_age_df(results_by_region, *, cell_type, region, age, table="simple_effects"):
    df = results_by_region[cell_type][region][table][age]
    return df.set_index("gene") if "gene" in df.columns else df

def _get_expr_series(results_by_region, *, cell_type, region,
                     stats_table="gene_stats", expr_col="mean_cpm_all"):
    pack = results_by_region[cell_type][region]
    if stats_table in pack and isinstance(pack[stats_table], pd.DataFrame):
        gs = pack[stats_table]
        if "gene" in gs.columns and expr_col in gs.columns:
            return gs.set_index("gene")[expr_col].squeeze()
    # fallback: average baseMean over ages if available
    se = pack.get("simple_effects", {})
    pieces = []
    for age, df in se.items():
        d = df.set_index("gene") if "gene" in df.columns else df
        if "baseMean" in d.columns:
            pieces.append(d[["baseMean"]].rename(columns={"baseMean": age}))
    if pieces:
        return pd.concat(pieces, axis=1).mean(axis=1).rename("mean_expr_fallback")
    return None

def _apply_expr_filter(res_a, res_b, expr_series, min_cpm=80, forced=None):
    forced = set(forced or [])
    genes = res_a.index.intersection(res_b.index)
    A, B = res_a.loc[genes].copy(), res_b.loc[genes].copy()
    if expr_series is not None:
        e = expr_series.reindex(genes).fillna(0)
        keep = (e >= float(min_cpm)) | pd.Index(genes).isin(forced)
        return A.loc[keep], B.loc[keep], e.loc[keep]
    return A, B, pd.Series(np.nan, index=genes)

# --------------------------------------------
# Display a small, clean top table
# --------------------------------------------
def _display_top_table(res, top, cell, age, head=30):
    if res is None or not isinstance(res, pd.DataFrame) or len(top) == 0:
        display(HTML(f"<h4 style='color:gray;margin:0.6em 0;'>⚠️ {cell} — {age}: no data</h4>"))
        return
    idx = res.index.intersection(pd.Index(top))
    if len(idx) == 0:
        display(HTML(f"<h4 style='color:gray;margin:0.6em 0;'>⚠️ {cell} — {age}: none of the requested genes found</h4>"))
        return
    sub = res.loc[idx].copy()
    num = sub.select_dtypes(include="number").columns
    sub[num] = sub[num].round(3)
    display(HTML(f"<h3 style='color:#2a4d69;margin:0.6em 0 0.2em 0;'>{cell} — Top {min(head,len(idx))} genes at {age}</h3>"))
    display(sub.head(head).sort_values(by = 'log2FoldChange', ascending = False))

# --------------------------------------------
# Pretty ΔLFC scatter (21w vs 60w) with adjustable axes
# --------------------------------------------
def compare_lfc_scatter_colored_labels(
    res_early, res_late,
    lfc_col="log2FoldChange",
    padj_col="padj",
    title=None,
    highlight_genes=None,
    expr=None,
    size_min=20, size_max=180,
    use_log_expr=True,
    q=0.05,
    figsize=(8, 8),
    fontsize=12,
    top_n_delta=15,
    save_base=None,
    save_formats=("png", "svg", "pdf"),
    dpi=300,
    xlim=None, ylim=None,
    symmetric_axes=True
):
    genes = res_early.index.intersection(res_late.index)
    df = (
        res_early.loc[genes, [lfc_col, padj_col]]
        .rename(columns={lfc_col: "lfc_21w", padj_col: "padj_21w"})
        .join(
            res_late.loc[genes, [lfc_col, padj_col]]
            .rename(columns={lfc_col: "lfc_60w", padj_col: "padj_60w"})
        )
        .replace([np.inf, -np.inf], np.nan)
        .dropna(subset=["lfc_21w", "lfc_60w"])
    )

    df["deltaLFC"] = df["lfc_60w"] - df["lfc_21w"]
    df["sig"] = (df["padj_21w"] < q) | (df["padj_60w"] < q)

    # point size ~ expression
    if expr is not None:
        df["expr"] = expr.reindex(df.index).fillna(0)
        base = np.log1p(df["expr"]) if use_log_expr else df["expr"]
        lo, hi = base.min(), base.max()
        rng = (hi - lo) if hi > lo else 1.0
        df["size"] = size_min + (base - lo) / rng * (size_max - size_min)
    else:
        df["expr"] = 0.0
        df["size"] = size_min

    # plot
    fig, ax = plt.subplots(figsize=figsize, dpi=150)
    sc = ax.scatter(df["lfc_21w"], df["lfc_60w"],
                    c=df["deltaLFC"], s=df["size"],
                    cmap="coolwarm", alpha=0.7, edgecolor="none")

    ax.scatter(df.loc[df["sig"], "lfc_21w"], df.loc[df["sig"], "lfc_60w"],
               c=df.loc[df["sig"], "deltaLFC"], cmap="coolwarm",
               s=df.loc[df["sig"], "size"] * 1.1,
               edgecolor="black", linewidth=0.3, alpha=0.9, zorder=3)

    # limits
    if xlim is None or ylim is None:
        lim_auto = float(np.nanmax(np.abs(df[["lfc_21w", "lfc_60w"]].to_numpy()))) * 1.05
        if symmetric_axes:
            xlim = xlim or (-lim_auto, lim_auto)
            ylim = ylim or (-lim_auto, lim_auto)
        else:
            xlim = xlim or (df["lfc_21w"].min()-0.5, df["lfc_21w"].max()+0.5)
            ylim = ylim or (df["lfc_60w"].min()-0.5, df["lfc_60w"].max()+0.5)
    ax.set_xlim(xlim); ax.set_ylim(ylim)

    # guides
    ax.axhline(0, color="k", ls="--", lw=1)
    ax.axvline(0, color="k", ls="--", lw=1)
    ax.plot([xlim[0], xlim[1]], [ylim[0], ylim[1]], color="k", ls=":", lw=1)

    ax.set_xlabel("log₂FC (21 weeks)", fontsize=fontsize)
    ax.set_ylabel("log₂FC (60 weeks)", fontsize=fontsize)
    ax.set_title(title or "ΔLFC comparison", fontsize=fontsize+2)

    # colorbar & size legend
    cb = fig.colorbar(sc, ax=ax); cb.set_label("Δ log₂FC (60 − 21)", fontsize=fontsize-1)
    if expr is not None:
        qs = np.percentile(df["expr"], [25, 50, 90])
        handles, labels = [], []
        base = np.log1p(df["expr"]) if use_log_expr else df["expr"]
        lo, hi = base.min(), base.max(); rng = (hi - lo) if hi > lo else 1.0
        for v in qs:
            vv = np.log1p(v) if use_log_expr else v
            s = size_min + (vv - lo) / rng * (size_max - size_min)
            handles.append(plt.scatter([], [], s=s, color="grey", alpha=0.5))
            labels.append(f"{v:.1f}")
        leg = ax.legend(handles, labels, title="Mean CPM", loc="lower right", frameon=False, fontsize=fontsize-2)
        ax.add_artist(leg)

    # labels (top ΔLFC + any provided)
    picks = set(highlight_genes or [])
    picks.update(df["deltaLFC"].abs().sort_values(ascending=False).head(top_n_delta).index)
    texts = []
    for g in picks:
        if g in df.index:
            texts.append(ax.text(df.at[g, "lfc_21w"], df.at[g, "lfc_60w"],
                                 g, fontsize=fontsize-3, ha="center", va="center"))
    if texts:
        adjust_text(texts, ax=ax, expand_points=(1.2, 1.3),
                    arrowprops=dict(arrowstyle="-", lw=0.7, color="0.4", alpha=0.8))

    plt.tight_layout()
    if save_base:
        for fmt in save_formats:
            fig.savefig(f"{save_base}.{fmt}", dpi=dpi, bbox_inches="tight")
        print(f"💾 Saved: {save_base}.{{{', '.join(save_formats)}}}")
    plt.show()
    return df

# --------------------------------------------
# Driver that pulls ages, filters, and plots
# --------------------------------------------
def run_age_contrast(
    results_by_region,
    *,
    cell_type, region,
    age_early="21", age_late="60",
    min_cpm=80,
    forced_genes=('Trib3','Gdf15','Fgf21','Atf4'),
    lfc_col="log2FoldChange", padj_col="padj",
    xlim=(-5,5), ylim=(-5,5),
    save=True
):
    res_early = _get_age_df(results_by_region, cell_type=cell_type, region=region, age=age_early)
    res_late  = _get_age_df(results_by_region, cell_type=cell_type, region=region, age=age_late)
    expr      = _get_expr_series(results_by_region, cell_type=cell_type, region=region)

    res_e_f, res_l_f, expr_f = _apply_expr_filter(res_early, res_late, expr, min_cpm=min_cpm, forced=forced_genes)

    # simple “top by |ΔLFC|” for the tables
    merged = (
        res_e_f[[lfc_col, padj_col]].rename(columns={lfc_col:"lfc_21w", padj_col:"padj_21w"})
        .join(res_l_f[[lfc_col, padj_col]].rename(columns={lfc_col:"lfc_60w", padj_col:"padj_60w"}))
    ).dropna()
    merged["deltaLFC"] = merged["lfc_60w"] - merged["lfc_21w"]
    top_delta = merged["deltaLFC"].abs().sort_values(ascending=False).head(25).index.tolist()
    labels = list(dict.fromkeys(top_delta + list(forced_genes)))

    safe = f"{cell_type}_{region}_{age_early}vs{age_late}".replace(" ", "_").replace("/", "_")
    save_base = os.path.join(outdir, safe) if save else None

    df_plot = compare_lfc_scatter_colored_labels(
        res_e_f, res_l_f,
        lfc_col=lfc_col, padj_col=padj_col,
        title=f"{cell_type} — {region}: {age_early} vs {age_late} (min CPM={min_cpm})",
        highlight_genes=labels,
        expr=expr_f,
        xlim=xlim, ylim=ylim,
        save_base=save_base
    )

    _display_top_table(res_e_f, top_delta, f"{cell_type} — {region}", f"{age_early} w")
    _display_top_table(res_l_f, top_delta, f"{cell_type} — {region}", f"{age_late} w")
    return df_plot


# In[307]:


for compartment in list(adata.obs.RBD_compartment_simplified.unique()): 
    try: 
        _ = run_age_contrast(
            results_by_region,
            cell_type="Mature oligodendrocytes",
            region=compartment,
            age_early="21",
            age_late="60",
            forced_genes=[],#('Trib3', 'Gdf15', 'Fgf21', 'Atf4'),
            min_cpm=50,
            xlim=(-7,7), ylim=(-7,7),
            save=False,

        )
    except KeyError:
        continue


# In[308]:


for compartment in list(adata.obs.RBD_compartment_simplified.unique()): 
    try: 
        _ = run_age_contrast(
            results_by_region,
            cell_type="Microglia",
            region=compartment,
            age_early="21",
            age_late="60",
            forced_genes=[],#('Trib3', 'Gdf15', 'Fgf21', 'Atf4'),
            min_cpm=50,
            xlim=(-8,8), ylim=(-8,8),
            save=False,

        )
    except KeyError:
        continue


# In[310]:


for compartment in list(adata.obs.RBD_compartment_simplified.unique()):
    try:
        _ = run_age_contrast(
            results_by_region,
            cell_type="Astrocytes",
            region=compartment,
            age_early="21",
            age_late="60",
            forced_genes=[],
            min_cpm=50,
            xlim=(-8,8), ylim=(-8,8),
            save=False,
        )
    except (KeyError, IndexError) as e:
        print(f"Skipping {compartment}: {type(e).__name__} — {e}")
        continue


# In[311]:


for compartment in list(adata.obs.RBD_compartment_simplified.unique()):
    try:
        _ = run_age_contrast(
            results_by_region,
            cell_type="Oligodendrocyte Precursor Cells",
            region=compartment,
            age_early="21",
            age_late="60",
            forced_genes=[],
            min_cpm=50,
            xlim=(-8,8), ylim=(-8,8),
            save=False,
        )
    except (KeyError, IndexError) as e:
        print(f"Skipping {compartment}: {type(e).__name__} — {e}")
        continue


# In[312]:


for compartment in list(adata.obs.RBD_compartment_simplified.unique()):
    try:
        _ = run_age_contrast(
            results_by_region,
            cell_type="Endothelial Cells",
            region=compartment,
            age_early="21",
            age_late="60",
            forced_genes=[],
            min_cpm=50,
            xlim=(-8,8), ylim=(-8,8),
            save=False,
        )
    except (KeyError, IndexError) as e:
        print(f"Skipping {compartment}: {type(e).__name__} — {e}")
        continue


# In[313]:


for compartment in list(adata.obs.RBD_compartment_simplified.unique()):
    try:
        _ = run_age_contrast(
            results_by_region,
            cell_type="Medium Spiny Neurons",
            region=compartment,
            age_early="21",
            age_late="60",
            forced_genes=[],
            min_cpm=50,
            xlim=(-8,8), ylim=(-8,8),
            save=False,
        )
    except (KeyError, IndexError) as e:
        print(f"Skipping {compartment}: {type(e).__name__} — {e}")
        continue


# In[314]:


for compartment in list(adata.obs.RBD_compartment_simplified.unique()):
    try:
        _ = run_age_contrast(
            results_by_region,
            cell_type="GABAergic Neurons",
            region=compartment,
            age_early="21",
            age_late="60",
            forced_genes=[],
            min_cpm=50,
            xlim=(-8,8), ylim=(-8,8),
            save=False,
        )
    except (KeyError, IndexError) as e:
        print(f"Skipping {compartment}: {type(e).__name__} — {e}")
        continue


# In[385]:


for compartment in list(adata.obs.RBD_compartment_simplified.unique()):
    try:
        _ = run_age_contrast(
            results_by_region,
            cell_type="Neural Stem Cells",
            region=compartment,
            age_early="21",
            age_late="60",
            forced_genes=[],
            min_cpm=30,
            xlim=(-8,8), ylim=(-8,8),
            save=True,
        )
    except (KeyError, IndexError) as e:
        print(f"Skipping {compartment}: {type(e).__name__} — {e}")
        continue


# In[316]:


for compartment in list(adata.obs.RBD_compartment_simplified.unique()):
    try:
        _ = run_age_contrast(
            results_by_region,
            cell_type="Excitatory Neurons",
            region=compartment,
            age_early="21",
            age_late="60",
            forced_genes=[],
            min_cpm=50,
            xlim=(-8,8), ylim=(-8,8),
            save=False,
        )
    except (KeyError, IndexError) as e:
        print(f"Skipping {compartment}: {type(e).__name__} — {e}")
        continue


# In[317]:


import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib import colors as mcolors

def plot_spatial_compact_fast(
    ad,
    color="leiden_2",          # obs column (categorical) OR a gene name (continuous)
    groupby="sample_id",
    spot_size=8,
    cols=3,
    height=8,
    legend_col_width=1.2,
    palette=None,              # for categorical only: dict {cat:"#hex"} or list
    rasterized=True,
    invert_y=True,
    dpi=120,
    # --- gene plotting params ---
    cmap="inferno",
    vmin=None,
    vmax=None,
    robust=True,               # if True and vmin/vmax are None, use robust percentiles
    robust_pct=(1, 99),        # lower/upper percentiles for robust scaling
    na_alpha=0.0,              # alpha for NaN gene values (0 = transparent)
    # --- disambiguation when name exists in both obs & var ---
    prefer="gene"              # one of {"gene","obs"}; default: prefer gene
):
    """
    If `color` is an obs column -> categorical plot with legend.
    If `color` is a gene in ad.var_names -> continuous plot with colorbar.
    If name exists in BOTH, `prefer` determines which path is used (default: gene).
    """
    # 0) Preconditions ---------------------------------------------------------
    if "spatial" not in ad.obsm:
        raise ValueError("ad.obsm['spatial'] is required.")

    coords = np.asarray(ad.obsm["spatial"])[:, :2]

    # membership checks (avoid array truthiness by using sets)
    obs_cols  = set(map(str, ad.obs.columns))
    var_names = set(map(str, ad.var_names))

    in_obs  = color in obs_cols
    in_var  = color in var_names

    if not (in_obs or in_var):
        raise KeyError(f"'{color}' not found as obs column or gene name in AnnData.")

    # choose path deterministically to avoid accidental doubles
    if in_obs and in_var:
        path = "gene" if prefer == "gene" else "obs"
    elif in_var:
        path = "gene"
    else:
        path = "obs"

    # 1) Grouping layout -------------------------------------------------------
    gvals = ad.obs[groupby].astype(str).to_numpy()
    uniq_groups, gcodes = np.unique(gvals, return_inverse=True)
    group_indices = [np.flatnonzero(gcodes == gi) for gi in range(len(uniq_groups))]

    n = len(uniq_groups)
    rows = int(np.ceil(n / max(cols, 1))) or 1

    panel_w = height * cols * 0.6 / rows
    fig_w = panel_w + legend_col_width

    plt.ioff()
    fig = plt.figure(figsize=(fig_w, height), dpi=dpi, constrained_layout=False)
    gs = GridSpec(
        rows, cols + 1, figure=fig,
        width_ratios=[1]*cols + [legend_col_width / max(1e-9, (fig_w - legend_col_width))],
        wspace=0.02, hspace=0.02
    )

    # 2A) Categorical path -----------------------------------------------------
    if path == "obs":
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
            base = list(mpl.cm.get_cmap("tab20").colors)
            reps = int(np.ceil(len(cat_names) / len(base))) if base else 1
            col_list = (base * max(reps,1))[:len(cat_names)]

        # store for consistency elsewhere
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
    else:  # path == "gene"
        # pull expression vector from ad[:, gene].X and make it dense 1-D float
        expr_raw = ad[:, color].X
        if hasattr(expr_raw, "A1"):            # scipy sparse .A1 (csr/csc)
            expr = expr_raw.A1.astype(float, copy=False)
        elif hasattr(expr_raw, "toarray"):     # generic sparse
            expr = expr_raw.toarray().ravel().astype(float, copy=False)
        else:
            expr = np.asarray(expr_raw).ravel().astype(float, copy=False)

        finite = np.isfinite(expr)
        if vmin is None or vmax is None:
            if robust and finite.any():
                lo, hi = robust_pct
                vmin_eff = float(np.nanpercentile(expr[finite], lo)) if vmin is None else vmin
                vmax_eff = float(np.nanpercentile(expr[finite], hi)) if vmax is None else vmax
            else:
                vmin_eff = (float(np.nanmin(expr[finite])) if finite.any() else None) if vmin is None else vmin
                vmax_eff = (float(np.nanmax(expr[finite])) if finite.any() else None) if vmax is None else vmax
        else:
            vmin_eff, vmax_eff = vmin, vmax

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
                # NaNs -> alpha
                alphas = np.where(np.isfinite(vals), 1.0, na_alpha)
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
        sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap_obj)
        cb = fig.colorbar(sm, cax=cax)
        cb.set_label(color)

    # finish
    fig.subplots_adjust(left=0.01, right=0.98, top=0.98, bottom=0.02, wspace=0.02, hspace=0.02)
    plt.ion()
    plt.show()
    return fig


# In[318]:


import re
import numpy as np
import pandas as pd
from scipy.stats import combine_pvalues

def _pick_cols(df, lfc_col=None, pval_col=None):
    """Pick LFC and pval columns from df; returns (lfc,pval) or (None,None)."""
    if df is None or not isinstance(df, pd.DataFrame) or df.empty:
        return None, None

    cols = {c.lower(): c for c in df.columns}
    # if user provided, respect when present
    if lfc_col in df.columns and pval_col in df.columns:
        return lfc_col, pval_col

    # common names (order = preference)
    lfc_candidates = [
        "log2fc_mtdsb_vs_control", "log2fc", "log2foldchange", "lfc", "log2_fc",
    ]
    pval_candidates = [
        "pvalue_mtdsb_vs_control", "padj", "p_adj", "pvalue", "pval",
    ]

    def find_one(cands, pattern=None):
        for nm in cands:
            if nm in cols: 
                return cols[nm]
        if pattern:
            hits = [orig for lower, orig in cols.items() if re.search(pattern, lower)]
            if hits:
                # prefer padj over pvalue if both appear
                hits_sorted = sorted(hits, key=lambda x: (not re.search("padj|p_adj", x.lower()), x))
                return hits_sorted[0]
        return None

    lfc = find_one(lfc_candidates, pattern=r"log2.*(fc|fold)")
    pvl = find_one(pval_candidates, pattern=r"p(adj|value|val)")
    return lfc, pvl

def _ensure_gene_index(df):
    """Make gene names the index if there's a 'gene' column."""
    if df is None or not isinstance(df, pd.DataFrame) or df.empty:
        return df
    if df.index.name is None and "gene" in df.columns:
        df = df.set_index("gene")
    return df

import numpy as np
import pandas as pd
from scipy.stats import combine_pvalues
from scipy.stats import norm

def summarize_condition_across_regions_flex(
    results_by_region,
    min_regions=2,
    lfc_col=None,
    pval_col=None,
    expr_col="baseMean",      # NEW: expression proxy to carry through
    combine="fisher",         # "fisher" or "stouffer"
    per_celltype=True,
    verbose=True,
):
    """
    Aggregate simple_effects across ALL regions (and ages) per cell type.

    We expect:
        results_by_region[cell_type][region]['simple_effects'][age] -> DataFrame

    Each DF should have:
        - log2FC column (lfc_col or autodetected)
        - p-value  column (pval_col or autodetected)
        - optionally an expression column (expr_col, e.g. 'baseMean')

    We return:
        dict[cell_type] -> DataFrame indexed by gene with:
            n_regions       (# distinct regions contributing)
            n_rows          (# total rows/region-age combos contributing)
            mean_lfc
            median_lfc
            mean_pval
            p_combined      (Fisher or Stouffer across regions)
            mean_expr       (avg expr_col across contributing regions/ages)
            abs_mean_lfc
    """

    def _ensure_gene_index(df):
        # if gene names are in a column instead of index, try to fix
        if not isinstance(df.index, pd.Index) or df.index.dtype == "int64":
            # heuristic: if there's a column literally called 'gene' or 'Gene'
            for cand in ["gene", "Gene", "symbol", "Symbol"]:
                if cand in df.columns:
                    df = df.set_index(cand)
                    break
        return df

    def _pick_cols(df, lfc_col=None, pval_col=None):
        # try user-specified first
        if lfc_col is not None and pval_col is not None:
            if lfc_col in df.columns and pval_col in df.columns:
                return lfc_col, pval_col

        # otherwise guess
        guess_lfc = [c for c in df.columns if "log2" in c.lower() or "lfc" in c.lower()]
        guess_pv  = [c for c in df.columns if "pval" in c.lower() or c.lower() == "pvalue"]

        lfc_name = guess_lfc[0] if guess_lfc else None
        p_name   = guess_pv[0]  if guess_pv  else None
        return lfc_name, p_name

    out = {}
    any_rows = 0

    for cell, by_region in results_by_region.items():
        if not isinstance(by_region, dict):
            if verbose:
                print(f"⚠️ Skip cell {cell}: not a dict")
            continue

        rows = []

        for region, payload in by_region.items():
            if not isinstance(payload, dict):
                if verbose:
                    print(f"  ⚠️ Skip region {region} in {cell}: payload not dict")
                continue

            se = payload.get("simple_effects", None)
            if not isinstance(se, dict) or not se:
                if verbose:
                    print(f"  ⚠️ {cell}/{region}: no 'simple_effects'")
                continue

            # iterate over ages (e.g. "21", "60")
            for age, df in se.items():
                df = _ensure_gene_index(df)

                lfc_name, p_name = _pick_cols(df, lfc_col=lfc_col, pval_col=pval_col)
                if df is None or lfc_name is None or p_name is None:
                    if verbose:
                        have = list(df.columns) if isinstance(df, pd.DataFrame) else None
                        print(
                            f"    ⚠️ {cell}/{region}/{age}: cannot find LFC/P columns (have={have})"
                        )
                    continue

                cols_to_keep = [lfc_name, p_name]
                has_expr = False
                if expr_col in df.columns:
                    cols_to_keep.append(expr_col)
                    has_expr = True

                sub = df[cols_to_keep].copy()
                sub.columns = ["lfc", "pval"] + (["expr"] if has_expr else [])

                sub["region"] = region
                sub["age"] = str(age)

                # clean
                sub = sub.replace([np.inf, -np.inf], np.nan).dropna(subset=["lfc", "pval"])
                sub = sub[(sub["pval"] >= 0) & (sub["pval"] <= 1)]
                if sub.empty:
                    if verbose:
                        print(f"    ⚠️ {cell}/{region}/{age}: empty after clean")
                    continue

                rows.append(sub)

        if not rows:
            if verbose:
                print(f"❕ No rows aggregated for {cell}")
            continue

        all_df = pd.concat(rows, axis=0)
        any_rows += len(all_df)

        # group by gene across regions+ages
        def _combine_func(df_sub):
            # how many distinct regions contributed this gene?
            n_reg = df_sub["region"].nunique()
            n_tot = df_sub.shape[0]

            mean_lfc = df_sub["lfc"].mean()
            med_lfc  = df_sub["lfc"].median()
            mean_pv  = df_sub["pval"].mean()

            # expression summary if present
            mean_expr = df_sub["expr"].mean() if "expr" in df_sub.columns else np.nan

            # combine p-values across regions/ages
            ps = df_sub["pval"].tolist()
            ps = [p for p in ps if np.isfinite(p) and 0 <= p <= 1]

            if len(ps) < min_regions:
                comb = np.nan
            else:
                if combine == "stouffer":
                    zs = [norm.isf(p) for p in ps if p > 0]
                    if len(zs) == 0:
                        comb = np.nan
                    else:
                        z_sum = np.sum(zs) / np.sqrt(len(zs))
                        comb = float(norm.sf(z_sum))
                else:  # Fisher
                    comb = float(combine_pvalues(ps, method="fisher")[1])

            return pd.Series({
                "n_regions": n_reg,
                "n_rows": n_tot,
                "mean_lfc": mean_lfc,
                "median_lfc": med_lfc,
                "mean_pval": mean_pv,
                "p_combined": comb,
                "mean_expr": mean_expr,          # <--- NEW
                "abs_mean_lfc": abs(mean_lfc),
            })

        agg = all_df.groupby(all_df.index).apply(_combine_func)

        # basic filter: require gene seen in >= min_regions
        agg = agg[agg["n_regions"] >= min_regions]

        # rank
        agg = agg.sort_values(
            ["abs_mean_lfc", "p_combined"],
            ascending=[False, True]
        )

        out[cell] = agg

    if verbose:
        if any_rows == 0:
            print(
                "🚫 Aggregation produced no rows. "
                "Check structure/column names.\n"
                "   Expecting results_by_region[cell][region]['simple_effects'][age] -> DataFrame"
            )
        else:
            print(f"✅ Aggregated {any_rows} rows across cell types.")

    return out


# In[319]:


summary_regions_by_age = region_activity_summary_by_age(
    results_by_region,
    lfc_col="log2FoldChange",
    pval_col="pvalue",
    p_thresh=0.05,
    min_expr=30,  # expression floor in mtDSB cells
    expr_col="mean_cpm_cond_mtDSB",
    lfc_direction="up"
)


# In[260]:


for ct, df in summary_by_cell.items():
    hits = df.query("p_combined < 0.05").sort_values("abs_mean_lfc", ascending=False)
    print(f"\n🔥 {ct}: {len(hits)} significant genes (FDR<0.05)")
    display(hits.head(20))
    print(list(hits.head(20).index))


# In[261]:


import numpy as np
import pandas as pd

def region_activity_summary_by_age(
    results_by_region,
    lfc_col="log2FoldChange",
    pval_col="pvalue",
    p_thresh=0.05,
    min_expr=None,          # e.g. 50 to require baseMean ≥ 50
    expr_col="baseMean",
    lfc_direction="up",     # "up", "down", or "both"
    verbose=True
):
    """
    Summarize region-level activation for each age group (e.g., 21 vs 60).
    Returns a tidy DataFrame with one row per (cell_type, region, age).
    """

    rows_out = []

    for cell_type, region_dict in results_by_region.items():
        if not isinstance(region_dict, dict):
            continue

        for region_name, payload in region_dict.items():
            se = payload.get("simple_effects", None)
            if not isinstance(se, dict) or not se:
                continue

            # --- iterate over age labels ("21", "60", etc.)
            for age_label, df in se.items():
                if not isinstance(df, pd.DataFrame) or df.empty:
                    continue

                df_here = df.copy()

                # ensure gene names are index
                if (df_here.index.name is None or df_here.index.dtype == "int64"):
                    for cand in ["gene", "Gene", "symbol", "Symbol"]:
                        if cand in df_here.columns:
                            df_here = df_here.set_index(cand)
                            break

                # check required columns
                if lfc_col not in df_here.columns or pval_col not in df_here.columns:
                    continue

                keep_cols = [lfc_col, pval_col]
                if expr_col in df_here.columns:
                    keep_cols.append(expr_col)

                sub = df_here[keep_cols].copy()
                sub = sub.replace([np.inf, -np.inf], np.nan).dropna(subset=[lfc_col, pval_col])
                sub = sub[(sub[pval_col] >= 0) & (sub[pval_col] <= 1)]

                # optional expression filter
                if (min_expr is not None) and (expr_col in sub.columns):
                    sub = sub[sub[expr_col] >= min_expr]

                # LFC direction
                if lfc_direction == "up":
                    sub = sub[sub[lfc_col] > 0]
                elif lfc_direction == "down":
                    sub = sub[sub[lfc_col] < 0]

                # significance mask
                sig_mask = sub[pval_col] < p_thresh
                sub_sig = sub[sig_mask]

                if sub_sig.empty:
                    n_sig = 0
                    mean_abs_lfc_sig = 0.0
                    top5_genes = []
                else:
                    collapsed = (
                        sub_sig
                        .groupby(sub_sig.index)
                        .agg(
                            mean_abs_lfc = (lfc_col, lambda x: np.mean(np.abs(x))),
                            mean_lfc     = (lfc_col, "mean"),
                            min_p        = (pval_col, "min"),
                            n_obs        = (pval_col, "size"),
                        )
                        .sort_values("mean_abs_lfc", ascending=False)
                    )
                    n_sig = collapsed.shape[0]
                    mean_abs_lfc_sig = collapsed["mean_abs_lfc"].mean()
                    top5_genes = list(collapsed.head(5).index)

                rows_out.append({
                    "cell_type": cell_type,
                    "region": region_name,
                    "age": str(age_label),
                    "n_sig_genes": n_sig,
                    "mean_abs_lfc_sig": mean_abs_lfc_sig,
                    "activation_index": n_sig * mean_abs_lfc_sig,
                    "top5_genes": top5_genes
                })

    summary_df = pd.DataFrame(rows_out)
    if not summary_df.empty:
        summary_df = summary_df.sort_values(
            ["cell_type", "age", "activation_index"],
            ascending=[True, True, False]
        )

    return summary_df


# In[263]:


summary_regions_by_age = region_activity_summary_by_age(
    results_by_region,
    lfc_col="log2FoldChange",
    pval_col="pvalue",
    p_thresh=0.05,
    min_expr=30,  # expression floor in mtDSB cells
    expr_col="mean_cpm_cond_mtDSB",
    lfc_direction="up"
)


# In[264]:


summary_regions_by_age


# In[265]:


ct  = 'Oligodendrocytes' 
_summary = summary_regions_by_age[
    summary_regions_by_age["cell_type"] == ct
]

import seaborn as sns
import matplotlib.pyplot as plt

plt.figure(figsize=(8, 5))
sns.barplot(
    data=_summary,
    x="region",
    y="activation_index",
    hue="age",
    palette="coolwarm"
)
plt.xticks(rotation=90)
plt.title("Oligodendrocyte activation index across regions (P21 vs P60)")
plt.tight_layout()
plt.show()


# In[266]:


def expand_summary_genes(summary_df):
    """Expand top5 gene lists into long format for plotting."""
    rows = []
    for _, r in summary_df.iterrows():
        for g in r["top5_genes"]:
            rows.append({
                "cell_type": r["cell_type"],
                "region": r["region"],
                "age": r["age"],
                "gene": g,
                "n_sig_genes": r["n_sig_genes"],
                "mean_abs_lfc_sig": r["mean_abs_lfc_sig"],
                "activation_index": r["activation_index"],
            })
    return pd.DataFrame(rows)


# In[267]:


summary_genes_long = expand_summary_genes(summary_regions_by_age)


# In[268]:


def plot_region_activation_fingerprint(df_cell, cell_type_name="cell type"):
    # sort by activation index descending
    df_plot = df_cell.sort_values("activation_index", ascending=True)

    y = df_plot["region"]
    x = df_plot["activation_index"]
    sizes = df_plot["n_sig_genes"]

    fig, ax = plt.subplots(figsize=(10,4))

    bars = ax.barh(
        y=y,
        width=x,
        color=plt.cm.Reds(np.interp(sizes, (sizes.min(), sizes.max()), (0.3, 1.0))),
        edgecolor="k",
        linewidth=0.5
    )

    # annotate with top genes + direction
    for bar, genes_info in zip(bars, df_plot["top5_genes_with_lfc"]):
        ax.text(
            bar.get_width() + 0.02 * x.max(),
            bar.get_y() + bar.get_height()/2,
            ", ".join([f"{g} ({lfc:+.1f})" for g, lfc in genes_info]),
            va="center",
            ha="left",
            fontsize=7,
        )

    ax.set_xlabel("Activation index (breadth × magnitude)")
    ax.set_title(f"{cell_type_name}: regional activation fingerprint")
    plt.tight_layout()
    plt.show()


# In[330]:


def filter_results_by_expression(
    results_by_region,
    expr_col="baseMean",
    min_expr=100,
    copy=True,
    verbose=True,
):
    """
    Recursively filters results_by_region so that only rows with
    expr_col >= min_expr are kept in each simple_effects DataFrame.

    Returns a new dict (unless copy=False, then modifies in place).
    """
    import copy as cp

    results_out = cp.deepcopy(results_by_region) if copy else results_by_region

    for cell_type, region_dict in results_out.items():
        if not isinstance(region_dict, dict):
            continue
        for region_name, payload in region_dict.items():
            se_dict = payload.get("simple_effects", None)
            if not isinstance(se_dict, dict):
                continue
            for age_label, df in se_dict.items():
                if not isinstance(df, pd.DataFrame):
                    continue
                if expr_col not in df.columns:
                    if verbose:
                        print(f"⚠️ {cell_type}/{region_name}/{age_label}: missing {expr_col}")
                    continue

                n_before = df.shape[0]
                df_filtered = df[df[expr_col] >= min_expr].copy()
                n_after = df_filtered.shape[0]

                results_out[cell_type][region_name]['simple_effects'][age_label] = df_filtered

                if verbose:
                    print(f"{cell_type}/{region_name}/{age_label}: {n_after}/{n_before} genes kept (>{min_expr} {expr_col})")

    return results_out


# In[332]:


results_by_region_filtered = filter_results_by_expression(
    results_by_region,
    expr_col="baseMean",   # or "mean_cpm_cond_mtDSB" if you want only active expression
    min_expr=50,          # threshold
    copy=True,             # create a new dict, keep original safe
    verbose=True,
)


# In[372]:


import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd

def _prep_age_df(df_age):
    """
    Sort regions so that highest activation_index ends up at the top visually.
    We sort ascending because barh plots from bottom to top.
    """
    return df_age.sort_values("activation_index", ascending=True).copy()

def _plot_single_age_panel(ax, df_plot, age_label, cmap="Reds", annotate_genes=True):
    """
    Draw one panel for a single age.

    df_plot must have:
      region
      activation_index
      n_sig_genes
      top_genes_with_lfc  (list[(gene, lfc), ...])
    """
    if df_plot.empty:
        ax.set_title(f"Age {age_label} (no signal)")
        ax.set_xlabel("Activation index\n(breadth × magnitude)")
        # Remove redundant y-label completely
        ax.set_ylabel(None)
        ax.tick_params(axis="y", labelsize=9)
        ax.axis("off")
        return None, None

    regions = df_plot["region"].tolist()
    act_idx = df_plot["activation_index"].to_numpy()
    n_sig   = df_plot["n_sig_genes"].to_numpy()

    # map n_sig_genes -> color
    vmin = np.min(n_sig) if len(n_sig) else 0
    vmax = np.max(n_sig) if len(n_sig) else 1
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    colors = sm.to_rgba(n_sig)

    bars = ax.barh(
        regions,
        act_idx,
        color=colors,
        edgecolor="black",
        linewidth=0.5
    )

    # annotate genes to the right of each bar
    if annotate_genes and "top_genes_with_lfc" in df_plot.columns:
        xmax = act_idx.max() if len(act_idx) else 1
        for bar, gene_list in zip(bars, df_plot["top_genes_with_lfc"]):
            if gene_list and isinstance(gene_list, list):
                label_txt = ", ".join(
                    [f"{g} ({lfc:+.1f})" for (g, lfc) in gene_list[:4]]
                )
                ax.text(
                    bar.get_width() + 0.02 * xmax,
                    bar.get_y() + bar.get_height()/2,
                    label_txt,
                    va="center",
                    ha="left",
                    fontsize=8,
                )

    ax.set_title(f"Age {age_label}", fontsize=12)
    ax.set_xlabel("Activation index\n(breadth × magnitude)")
    #ax.set_ylabel("Region")
    ax.grid(axis="x", linestyle=":", alpha=0.4)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    return sm, norm


def plot_region_activation_fingerprint_by_age(
    summary_regions_by_age,
    cell_type_name,
    ages=("21", "60"),
    figsize=(14,7),
    cmap="Reds",
    annotate_genes=True,
    share_x=True,
):
    """
    Side-by-side panels comparing activation fingerprints across regions
    for two ages (e.g. 21 vs 60) in a given cell type.
    Colorbar is placed horizontally under both panels to avoid overlap.
    """

    # slice this cell type
    df_ct = summary_regions_by_age[
        summary_regions_by_age["cell_type"] == cell_type_name
    ].copy()

    # slice data for each requested age
    dfs = []
    for age in ages:
        dfa = df_ct[df_ct["age"] == str(age)].copy()
        dfa = _prep_age_df(dfa)
        dfs.append(dfa)

    df_left, df_right = dfs

    # make figure: 2 columns for ages, shared y? we keep share_x and adjust manually
    fig, axes = plt.subplots(
        1, 2,
        figsize=figsize,
        sharex=share_x,
        gridspec_kw={"wspace": 1}
    )

    # left panel
    sm_left, norm_left = _plot_single_age_panel(
        axes[0],
        df_left,
        age_label=str(ages[0]),
        cmap=cmap,
        annotate_genes=annotate_genes,
    )

    # right panel
    sm_right, norm_right = _plot_single_age_panel(
        axes[1],
        df_right,
        age_label=str(ages[1]),
        cmap=cmap,
        annotate_genes=annotate_genes,
    )

    # sync x-axis limits across both if share_x
    if share_x:
        xmax = 0.0
        if not df_left.empty:
            xmax = max(xmax, df_left["activation_index"].max())
        if not df_right.empty:
            xmax = max(xmax, df_right["activation_index"].max())
        if xmax <= 0:
            xmax = 1.0
        for ax in axes:
            ax.set_xlim(0, xmax * 1.2)

    # ---- construct shared colorbar data from both ages ----
    # pull n_sig_genes from both dataframes
    n_all = []
    for d in [df_left, df_right]:
        if "n_sig_genes" in d.columns and not d.empty:
            n_all.extend(d["n_sig_genes"].tolist())
    if len(n_all) == 0:
        vmin, vmax = 0, 1
    else:
        vmin, vmax = min(n_all), max(n_all)

    norm_all = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    sm_all = mpl.cm.ScalarMappable(norm=norm_all, cmap=cmap)
    sm_all.set_array([])

    # ---- ADD HORIZONTAL COLORBAR UNDER BOTH AXES ----
    # We'll create a new axes that spans the bottom
        # ---- ADD HORIZONTAL COLORBAR UNDER BOTH AXES ----
    # [left, bottom, width, height] in figure fraction
    # narrower (width=0.4), lower (bottom=0.02), and thin (height=0.02)
    cbar_ax = fig.add_axes([0.3, 0.02, 0.4, 0.02])

    cbar = plt.colorbar(
        sm_all,
        cax=cbar_ax,
        orientation="horizontal",
    )
    cbar.set_label("# sig. upregulated genes (FDR<0.05)", fontsize=10)
    cbar.ax.tick_params(labelsize=8, length=3)

    # main title
    fig.suptitle(
        f"{cell_type_name}: Regional activation fingerprints by age",
        fontsize=14,
        y=0.98,
    )

    # increase bottom spacing slightly to avoid any tick overlap
    fig.subplots_adjust(top=0.88, bottom=0.17, left=0.08, right=0.97)

    plt.show()

    return fig, axes

############################################
# 3. example usage
############################################
# 0. Filter each simple_effects table in results_by_region
results_by_region_filtered = filter_results_by_expression(
    results_by_region,
    expr_col="baseMean",   # we trust baseMean as a general "this gene isn't garbage"
    min_expr=50,           # keep genes with baseMean >= 50
    copy=True,             # don't overwrite the original
    verbose=True,
)

# 1. Summarize per cell type / region / age,
#    using STRICT biologically meaningful criteria
summary_regions_by_age = region_activity_summary_by_age_clean(
    results_by_region_filtered,
    lfc_col="log2FoldChange",
    pval_col="pvalue",
    expr_col_mt="baseMean",  # expression in mtDSB condition
    p_thresh=0.05,        # must be significant
    lfc_min=0.5,          # must be meaningfully up (~1.4x+)
    expr_min_cpm=30,      # must actually be expressed in mtDSB cells in that region
    top_n_genes=4,
    verbose=True,
)


# In[383]:


results_by_region_filtered['Neural Stem Cells']['Ventricular system']['simple_effects']['60'].sort_values(by = 'log2FoldChange', ascending = False).head(50)


# In[375]:


import os

# --- Define output directory ---
out_dir = "../results/figures/region_fingerprints_by_age"
os.makedirs(out_dir, exist_ok=True)

# --- Loop across all cell types ---
for ct in summary_regions_by_age.cell_type.unique():
    print(f"Plotting {ct} …")

    # Create the figure
    fig, axes = plot_region_activation_fingerprint_by_age(
        summary_regions_by_age,
        cell_type_name=ct,
        ages=("21", "60"),
        figsize=(14, 6),
        cmap="Reds",
        annotate_genes=True,
        share_x=True,
    )

    # --- Save the figure ---
    base = os.path.join(out_dir, f"{ct.replace(' ', '_')}_activation_by_age")
    for ext in ["png"]:
        fig.savefig(f"{base}.{ext}", bbox_inches="tight", dpi=300)
    plt.close(fig)

print(f"✅ Saved all figures to: {out_dir}")


# In[367]:


fig = plot_top_genes_by_celltype(summary_by_cell, p_cutoff=0.05, top_n=10, n_cols=3)
#fig.savefig("results/figures/top_genes_by_celltype.png", dpi=300, bbox_inches="tight")


# In[368]:


def summarize_global_hits(summary_by_cell, p_cutoff=0.05):
    """
    Find genes that are recurrently significant across cell types.
    """
    records = []
    for ct, df in summary_by_cell.items():
        df_ct = df.copy()
        filt = df_ct.query("p_combined < @p_cutoff").copy()
        if filt.empty:
            continue
        for gene, row in filt.iterrows():
            records.append({
                "gene": gene,
                "cell_type": ct,
                "mean_lfc": row["mean_lfc"],
                "abs_mean_lfc": row["abs_mean_lfc"],
                "p_combined": row["p_combined"],
            })

    if not records:
        return pd.DataFrame()

    long_df = pd.DataFrame(records)
    summary = (
        long_df
        .groupby("gene")
        .agg(
            in_how_many_celltypes=("cell_type", "nunique"),
            celltypes=("cell_type", lambda xs: ",".join(sorted(set(xs)))),
            mean_abs_lfc=("abs_mean_lfc", "mean"),
            max_abs_lfc=("abs_mean_lfc", "max"),
        )
        .sort_values(
            ["in_how_many_celltypes", "mean_abs_lfc"],
            ascending=[False, False]
        )
    )
    return summary

