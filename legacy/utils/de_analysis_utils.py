"""
Auto-generated utilities for de_analysis_utils.
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

def screen_condition_ols(df_long: pd.DataFrame, cell_type: str,
                         condition_col="condition",
                         value_col="value",
                         sample_col="sample"):
    """
    Fast OLS screen: value ~ 1 + condition (aggregated per sample to avoid pseudorep).
    Returns tidy DF with columns:
      ['gene','cell_class','ols_effect','ols_se','ols_p','n_samples','n_groups']
    """
    d = df_long[df_long["cell_class"] == cell_type].copy()
    if d.empty:
        return pd.DataFrame(columns=["gene","cell_class","ols_effect","ols_se","ols_p","n_samples","n_groups"])

    # aggregate to per-sample means (prevents many per-sample replicates)
    agg = (d
           .groupby(["gene", sample_col, condition_col], as_index=False)[value_col]
           .mean())

    genes = sorted(agg["gene"].unique())
    rows = []
    for g in tqdm(genes, desc=f"Screening (OLS) in {cell_type}"):
        sub = agg[agg["gene"] == g].copy()
        # drop genes with <2 unique samples or single condition
        if sub[sample_col].nunique() < 2 or sub[condition_col].nunique() < 2:
            continue
        X = pd.get_dummies(sub[condition_col], drop_first=True)  # 0/1 for case vs control
        # If your positive level should be mtDSB, ensure categories ordered
        # e.g., sub[condition_col] = pd.Categorical(sub[condition_col], ["control","mtDSB"])
        X = sm.add_constant(X)
        y = sub[value_col].values
        try:
            fit = sm.OLS(y, X).fit()
            # take the single condition coefficient (last column)
            coef_name = [c for c in X.columns if c != "const"][-1]
            beta = fit.params[coef_name]
            se   = fit.bse[coef_name]
            pval = fit.pvalues[coef_name]
            rows.append({
                "gene": g,
                "cell_class": cell_type,
                "ols_effect": float(beta),
                "ols_se": float(se),
                "ols_p": float(pval),
                "n_samples": int(sub[sample_col].nunique()),
                "n_groups": len(sub)
            })
        except Exception:
            # e.g., singular matrix; skip
            continue

    return pd.DataFrame(rows)

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

def _pick_score_col(df, preferred=("score","n_sig_genes","n_genes","n_regions",
                                   "sum_neglog10FDR","mean_neglog10FDR","mean_LFC")):
    num_cols = [c for c in df.columns if np.issubdtype(df[c].dtype, np.number)]
    for c in preferred:
        if c in num_cols: 
            return c
    # fallback: first numeric column
    return num_cols[0] if num_cols else None

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

def compare_lfc_scatter_colored_labels(
    res_early, res_late,
    lfc_col="log2FC_mtDSB_vs_control",
    padj_col="pvalue_mtDSB_vs_control",
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
    xlim=None, ylim=None,   # 👈 NEW optional args
    symmetric_axes=True     # 👈 NEW toggle to enforce symmetric limits
):
    """ΔLFC scatter with colorbar and size legend (now with adjustable axes)."""
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

    # --- Scale size by expression if provided ---
    if expr is not None:
        df["expr"] = expr.reindex(df.index).fillna(0)
        if use_log_expr:
            df["size"] = np.interp(np.log1p(df["expr"]),
                                   (np.log1p(df["expr"]).min(), np.log1p(df["expr"]).max()),
                                   (size_min, size_max))
        else:
            df["size"] = np.interp(df["expr"], 
                                   (df["expr"].min(), df["expr"].max()),
                                   (size_min, size_max))
    else:
        df["size"] = size_min

    # --- Plot ---
    fig, ax = plt.subplots(figsize=figsize, dpi=150)
    sc = ax.scatter(df["lfc_21w"], df["lfc_60w"],
                    c=df["deltaLFC"], s=df["size"],
                    cmap="coolwarm", alpha=0.7,
                    edgecolor="none")

    # outline for significant
    ax.scatter(df.loc[df["sig"], "lfc_21w"],
               df.loc[df["sig"], "lfc_60w"],
               c=df.loc[df["sig"], "deltaLFC"],
               cmap="coolwarm", s=df.loc[df["sig"], "size"] * 1.1,
               edgecolor="black", linewidth=0.3, alpha=0.9, zorder=3)

    # --- Axis limits ---
    if xlim is None or ylim is None:
        lim_auto = np.nanmax(np.abs(df[["lfc_21w", "lfc_60w"]])) * 1.05
        if symmetric_axes:
            xlim = xlim or (-lim_auto, lim_auto)
            ylim = ylim or (-lim_auto, lim_auto)
        else:
            xlim = xlim or (df["lfc_21w"].min() - 0.5, df["lfc_21w"].max() + 0.5)
            ylim = ylim or (df["lfc_60w"].min() - 0.5, df["lfc_60w"].max() + 0.5)

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    # diagonal + reference lines
    ax.axhline(0, color="k", ls="--", lw=1)
    ax.axvline(0, color="k", ls="--", lw=1)
    ax.plot([xlim[0], xlim[1]], [ylim[0], ylim[1]], color="k", ls=":", lw=1)

    ax.set_xlabel("log₂FC (21 weeks)", fontsize=fontsize)
    ax.set_ylabel("log₂FC (60 weeks)", fontsize=fontsize)
    ax.set_title(title or "ΔLFC comparison", fontsize=fontsize + 2)

    # --- Colorbar (ΔLFC legend) ---
    cb = fig.colorbar(sc, ax=ax)
    cb.set_label("Δ log₂FC (60 w – 21 w)", fontsize=fontsize - 1)

    # --- Size legend (expression proxy) ---
    if expr is not None:
        handles = []
        for val in np.percentile(df["expr"], [25, 50, 90]):
            handles.append(plt.scatter([], [], s=np.interp(np.log1p(val),
                                                           (np.log1p(df["expr"]).min(), np.log1p(df["expr"]).max()),
                                                           (size_min, size_max)),
                                       color="grey", alpha=0.5))
        labels = [f"{v:.1f}" for v in np.percentile(df["expr"], [25, 50, 90])]
        legend = ax.legend(handles, labels, title="Mean CPM", 
                           loc="lower right", frameon=False, fontsize=fontsize - 2)
        ax.add_artist(legend)

    # --- Label top genes ---
    picks = set(highlight_genes or [])
    picks.update(df["deltaLFC"].abs().sort_values(ascending=False).head(top_n_delta).index)
    texts = []
    for g in picks:
        if g not in df.index:
            continue
        x, y = df.at[g, "lfc_21w"], df.at[g, "lfc_60w"]
        texts.append(ax.text(x, y, g, fontsize=fontsize - 3, ha="center", va="center"))
    adjust_text(texts, ax=ax, expand_points=(1.2, 1.3),
                arrowprops=dict(arrowstyle="-", lw=0.7, color="0.4", alpha=0.8))

    plt.tight_layout()

    # --- Save ---
    if save_base:
        for fmt in save_formats:
            fig.savefig(f"{save_base}.{fmt}", dpi=dpi, bbox_inches="tight")
        print(f"💾 Saved: {save_base}.{{{', '.join(save_formats)}}}")

    plt.show()
    return df

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

def summarize_de_by_region(
    results_by_region,
    cell_class,
    regions=None,          # list or None
    ages=None,             # list or None
    padj_thr=0.05,
    lfc_thr=1.0,
    top_n=500,             # how many rows to pull from extract_top_genes
    top_k_display=30,      # how many rows to show in the styled table
    caption="Top genes by number of regions detected",
    show=True,             # display styled table in notebook
    save_csv_path=None,    # e.g. "summary_by_region.csv"
):
    """
    Run extract_top_genes() and summarize per gene across regions.
    Returns (summary_df, styled_table_or_None).
    """
    # 1) collect significant genes across region/age
    top_genes = extract_top_genes(
        results_by_region=results_by_region,
        cell_class=cell_class,
        regions=regions,
        ages=ages,
        padj_thr=padj_thr,
        lfc_thr=lfc_thr,
        top_n=top_n,
    )
    if top_genes is None or top_genes.empty:
        print("No genes passed thresholds.")
        return pd.DataFrame(), None

    # 2) summarize per gene (index is gene name)
    # robust to any index name by grouping on index level 0
    summary = (
        top_genes.groupby(level=0)
        .agg(
            n_regions=("region", "nunique"),
            mean_LFC=("log2FC", "mean"),
            mean_neglog10FDR=("neglog10FDR", "mean"),
        )
        .sort_values("n_regions", ascending=False)
    )

    if save_csv_path:
        summary.to_csv(save_csv_path)

    # 3) optional pretty styling for slides
    styled = (
        summary.head(top_k_display)
        .style
        .background_gradient(subset="n_regions", cmap="Purples")
        .bar(subset="mean_LFC", color="lightblue")
        .format({"mean_LFC": "{:.2f}", "mean_neglog10FDR": "{:.1f}"})
        .set_caption(caption)
        .set_table_styles([{'selector': 'caption',
                            'props': [('font-size', '16px'), ('font-weight', 'bold')]}])
    )

    if show:
        display(styled)

    return summary, styled

def plot_de_summary_bars(summary: pd.DataFrame, top_k=30, save=None, title=None):
    """
    Two-panel figure:
      (A) n_regions per gene (horizontal bars)
      (B) mean_LFC per gene (colored by sign, value labels optional)
    """
    df = summary.head(top_k).copy()
    df = df.iloc[::-1]  # plot top at top
    fig, axes = plt.subplots(1, 2, figsize=(12, 0.35*len(df)+2), gridspec_kw=dict(wspace=0.25))

    # Panel A: n_regions
    sns.barplot(ax=axes[0], y=df.index, x="n_regions", data=df, color="#6C5CE7")
    axes[0].set_title("(A) Regions per gene")
    axes[0].set_xlabel("# regions"); axes[0].set_ylabel("")

    # Panel B: mean_LFC with sign color
    palette = df["mean_LFC"].apply(lambda v: "#2ECC71" if v>0 else "#E74C3C")
    axes[1].barh(df.index, df["mean_LFC"], color=palette)
    axes[1].axvline(0, color="k", lw=0.7, alpha=0.6)
    axes[1].set_title("(B) Mean log2FC")
    axes[1].set_xlabel("mean log2FC"); axes[1].set_ylabel("")

    if title:
        fig.suptitle(title, y=1.02, fontsize=14, fontweight="bold")
    plt.tight_layout()
    if save:
        fig.savefig(save, dpi=300, bbox_inches="tight")
    plt.show()

def plot_de_summary_bubbles(summary: pd.DataFrame, top_k=30, save=None, title=None):
    """
    Bubble plot: x = mean_LFC, y = gene, size = n_regions, color = mean_neglog10FDR
    """
    df = summary.head(top_k).copy().sort_values("mean_LFC")
    plt.figure(figsize=(8, 0.35*len(df)+2))
    sc = plt.scatter(
        df["mean_LFC"], np.arange(len(df)),
        s=60 + 60*df["n_regions"],    # bigger = found in more regions
        c=df["mean_neglog10FDR"], cmap="coolwarm", edgecolor="k", linewidth=0.3
    )
    plt.yticks(np.arange(len(df)), df.index)
    plt.xlabel("mean log2FC"); plt.ylabel("")
    plt.axvline(0, color="k", lw=0.7, alpha=0.6)
    cbar = plt.colorbar(sc); cbar.set_label("mean -log10(FDR)")
    if title:
        plt.title(title)
    plt.tight_layout()
    if save:
        plt.savefig(save, dpi=300, bbox_inches="tight")
    plt.show()

def plot_region_presence_heatmap(top_genes: pd.DataFrame, top_k=30, save=None, title=None):
    """
    Heatmap of region presence/sign for the top_k recurrent genes.
    Expects columns: ['region','log2FC'] in top_genes (long format).
    Cell values: sign(log2FC) ∈ {-1, 0, +1}
    """
    genes = (top_genes.index.value_counts().sort_values(ascending=False).head(top_k).index)
    df = top_genes.loc[genes]
    mat = (df.assign(sign=np.sign(df["log2FC"]))
             .pivot_table(index=df.index, columns="region", values="sign", aggfunc="mean")
             .fillna(0.0)
             .sort_index())

    plt.figure(figsize=(0.6*mat.shape[1]+4, 0.35*mat.shape[0]+2))
    sns.heatmap(mat, cmap=sns.color_palette(["#E74C3C", "#EEEEEE", "#2ECC71"], as_cmap=True),
                vmin=-1, vmax=1, linewidths=0.5, linecolor="white", cbar=False)
    plt.xlabel("Region"); plt.ylabel("")
    if title:
        plt.title(title)
    plt.tight_layout()
    if save:
        plt.savefig(save, dpi=300, bbox_inches="tight")
    plt.show()

def region_gene_matrix(results_by_region, cell_class, genes, ages=None, padj_thr=0.05):
    """
    Build region×gene matrix of log2FC for a given cell_class.
    Rows = "region @ age", Cols = genes.
    Non-significant entries (padj >= padj_thr) are set to NaN.
    """
    rows = []
    res_cl = results_by_region[cell_class]

    for reg, block in res_cl.items():
        se = block.get("simple_effects", {})
        if not se:
            continue

        # Map available age keys to strings so we can match robustly
        age_key_map = {str(k): k for k in se.keys()}
        ages_iter = [str(a) for a in (ages if ages is not None else se.keys())]

        for age_str in ages_iter:
            key = age_key_map.get(age_str)
            if key is None:
                continue
            res = se.get(key)
            if res is None or res.empty:
                continue
            # keep only requested genes that exist
            present = res.index.intersection(genes)
            if len(present) == 0:
                continue

            # standard DESeq2 columns: log2FoldChange, padj
            sub = res.loc[present, ["log2FoldChange", "padj"]].copy()
            # mask non-significant to NaN
            sub.loc[sub["padj"] >= padj_thr, "log2FoldChange"] = np.nan

            row = sub["log2FoldChange"]
            row.name = f"{reg} @ {age_str}"
            rows.append(row)

    if not rows:
        return pd.DataFrame(columns=genes)

    mat = pd.DataFrame(rows).reindex(columns=genes)
    return mat

def heatmap(mat, title):
    if mat.empty:
        print("No data to plot for:", title)
        return
    arr = mat.values.astype(float)
    vmax = np.nanmax(np.abs(arr))
    vmax = 1.0 if not np.isfinite(vmax) or vmax == 0 else vmax
    plt.figure(figsize=(0.5*mat.shape[1]+4, 0.4*mat.shape[0]+3), dpi=150)
    im = plt.imshow(arr, aspect="auto", vmin=-vmax, vmax=+vmax)
    plt.colorbar(im, label="log2FC (mtDSB vs control)")
    plt.yticks(range(mat.shape[0]), mat.index)
    plt.xticks(range(mat.shape[1]), mat.columns, rotation=60, ha="right")
    plt.title(title)
    plt.tight_layout()
    plt.show()

def make_shortlist(screen_df: pd.DataFrame, top_k=300, by="abs_effect",
                   effect_col=None, p_col=None):
    """
    Create a shortlist of genes from a screen DataFrame.
    - by='abs_effect' ranks by |effect|
    - by='fdr' ranks by BH-FDR on p-values
    Auto-detects column names if not provided.
    """
    if screen_df is None or screen_df.empty:
        return []

    df = screen_df.copy()

    # Infer column names if not provided
    if effect_col is None:
        for cand in ["ols_effect","effect","beta","coef","estimate"]:
            if cand in df.columns:
                effect_col = cand; break
    if p_col is None:
        for cand in ["ols_p","p","pval","p_value"]:
            if cand in df.columns:
                p_col = cand; break

    if effect_col is None:
        raise KeyError("No effect column found (looked for ols_effect/effect/beta/coef/estimate).")
    if p_col is None and by == "fdr":
        raise KeyError("No p-value column found (looked for ols_p/p/pval/p_value).")

    # Compute FDR if requested and available
    if p_col and p_col in df.columns:
        df["ols_fdr"] = multipletests(df[p_col].values, method="fdr_bh")[1]
    else:
        df["ols_fdr"] = np.nan

    # Ranking
    if by == "abs_effect":
        df = df.reindex(df[effect_col].abs().sort_values(ascending=False).index)
    elif by == "fdr":
        if "ols_fdr" not in df or df["ols_fdr"].isna().all():
            raise ValueError("FDR selected but p-values were not available to compute it.")
        df = df.sort_values("ols_fdr", ascending=True)
    else:
        raise ValueError("by must be 'abs_effect' or 'fdr'.")

    # Keep top_k unique genes
    out = df.drop_duplicates("gene").head(top_k)
    return out["gene"].tolist()

def screen_condition_ols_robust(df_long: pd.DataFrame, cell_type: str,
                                positive="mtDSB", reference="control",
                                min_total_samples=2, min_conditions=2):
    """
    Robust fast screen of condition effect within a cell type.
    Returns columns: ['gene','cell_class','ols_effect','ols_se','ols_p','n_samples','n_groups','method']
    """
    std = _detect_cols(df_long)
    cols_out = ["gene","cell_class","ols_effect","ols_se","ols_p","n_samples","n_groups","method"]

    # subset cell type (be exact on the string!)
    d = df_long[df_long[std["cell_class"]] == cell_type].copy()
    if d.empty:
        print(f"[screen] No rows for cell type: {cell_type}")
        return pd.DataFrame(columns=cols_out)

    # enforce condition coding (reference vs positive)
    try:
        d = _coerce_condition(d, std["condition"], positive=positive, reference=reference)
    except Exception as e:
        print("[screen] Condition coercion failed:", e)
        return pd.DataFrame(columns=cols_out)

    # per-sample mean to avoid pseudoreplication
    agg = (d.groupby([std["gene"], std["sample"], std["condition"]], as_index=False)[std["value"]]
             .mean()
             .rename(columns={std["value"]:"y"}))

    # quick diagnostics
    print(f"[screen] {cell_type}: rows={len(d)}, agg_rows={len(agg)}")
    print("[screen] condition counts:", d[std["condition"]].value_counts().to_dict())
    print("[screen] unique genes:", agg[std["gene"]].nunique())

    out = []
    genes = agg[std["gene"]].unique().tolist()
    for g in tqdm(genes, desc=f"Screening in {cell_type}"):
        sub = agg[agg[std["gene"]] == g]
        n_total = sub[std["sample"]].nunique()
        n_cond  = sub[std["condition"]].nunique()
        if (n_total < min_total_samples) or (n_cond < min_conditions):
            continue

        # design: y ~ 1 + I(condition == positive)
        X = pd.get_dummies(sub[std["condition"]], drop_first=True)  # column for positive vs reference
        if X.shape[1] == 0:
            # only one level present (shouldn't happen after checks)
            continue
        X = sm.add_constant(X)
        y = sub["y"].values

        # OLS fit; if it fails, fall back to Welch t-test
        try:
            fit = sm.OLS(y, X).fit()
            coef_name = [c for c in X.columns if c != "const"][-1]
            beta = float(fit.params.get(coef_name, np.nan))
            se   = float(fit.bse.get(coef_name, np.nan))
            pval = float(fit.pvalues.get(coef_name, np.nan))
            method = "OLS"
        except Exception:
            # fall back: Welch t-test on per-sample means
            g0 = sub[sub[std["condition"]].astype(str) == reference]["y"].values
            g1 = sub[sub[std["condition"]].astype(str) == positive]["y"].values
            if len(g0) < 1 or len(g1) < 1:
                continue
            beta = float(np.nanmean(g1) - np.nanmean(g0))
            # SE from two-sample stats (rough)
            se = float(np.sqrt(np.nanvar(g0, ddof=1)/max(len(g0),1) + np.nanvar(g1, ddof=1)/max(len(g1),1)))
            pval = float(stats.ttest_ind(g1, g0, equal_var=False, nan_policy="omit").pvalue)
            method = "welch_t"

        out.append({
            "gene": g,
            "cell_class": cell_type,
            "ols_effect": beta,
            "ols_se": se,
            "ols_p": pval,
            "n_samples": int(n_total),
            "n_groups": int(len(sub)),
            "method": method
        })

    res = pd.DataFrame(out, columns=cols_out)
    if res.empty:
        print("[screen] Result is empty after fitting. Check column names and condition balance.")
    return res

def plot_pathway_panels(
    cts,
    pathway_genes,
    module_name="Pathway",
    results_root="../results/age",
    age_key="mean_cpm_age_60",
    expr_min=30,
    p_cutoff=0.1,
    n_cols=3,
    shared_xlim=None,       # e.g. (-1, 6) or None for auto
    save_dir=None,          # <-- NEW: optional save directory
    save_formats=("png",),  # <-- NEW: save formats
    dpi=300,                # <-- NEW: resolution
    show=True,              # <-- NEW: whether to display figure
):
    """
    Multi-panel barplots of fold-change for a gene module across cell types.

    pathway_genes : list of str
        Genes to include for this pathway/module.
    """

    import numpy as np
    import matplotlib.pyplot as plt
    import math

    n_ct = len(cts)
    n_rows = math.ceil(n_ct / n_cols)
    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(n_cols * 5.5, n_rows * 5.5),
        sharex=True
    )
    axes = axes.flatten() if isinstance(axes, np.ndarray) else [axes]

    # infer which file to load
    age_suffix = "60" if "60" in age_key else "21"

    for i, ct in enumerate(cts):
        ax = axes[i]
        fn = os.path.join(results_root, ct, f"simple_effect_age_{age_suffix}.csv")

        if not os.path.exists(fn):
            print(f"Skipping {ct}: file not found at {fn}")
            ax.axis("off")
            continue

        df = pd.read_csv(fn, index_col=0)
        df_t = df.T  # rows = metrics, cols = genes

        # intersect this module with available genes
        genes_here = [g for g in pathway_genes if g in df_t.columns]
        if not genes_here:
            print(f"Skipping {ct}: none of {module_name} genes found.")
            ax.axis("off")
            continue

        # slice to pathway genes
        df_mod = df_t[genes_here]

        # confirm required rows
        required_rows = ["log2FC_mtDSB_vs_control", "pvalue_mtDSB_vs_control", age_key]
        if any(r not in df_mod.index for r in required_rows):
            print(f"{ct}: missing required rows in table, skipping.")
            ax.axis("off")
            continue

        # filter low expression
        df_filt = df_mod.T[df_mod.T[age_key] > expr_min].T
        if df_filt.shape[1] == 0:
            print(f"{ct}: no {module_name} genes pass expr_min={expr_min} for {age_key}")
            ax.axis("off")
            continue

        # tidy table
        plot_df = pd.DataFrame({
            "gene": df_filt.columns,
            "log2FC": df_filt.loc["log2FC_mtDSB_vs_control"].values,
            "pval": df_filt.loc["pvalue_mtDSB_vs_control"].values,
            "expr": df_filt.loc[age_key].values,
        })

        plot_df["sig_flag"] = plot_df["pval"] < p_cutoff
        plot_df = plot_df.sort_values("log2FC", ascending=True)

        colors = ["#b30000" if sig else "#cccccc" for sig in plot_df["sig_flag"]]

        ax.barh(
            plot_df["gene"],
            plot_df["log2FC"],
            color=colors,
            edgecolor="black",
            linewidth=0.6,
        )

        ax.axvline(0, color="k", lw=0.8, ls="--")

        ax.set_title(ct.replace("_", " "), fontsize=11)
        ax.set_xlabel("log₂FC (mtDSB vs control)")
        ax.set_ylabel("Gene" if i % n_cols == 0 else "")

        # consistent axis
        if isinstance(shared_xlim, tuple):
            ax.set_xlim(shared_xlim)
            x_text = shared_xlim[1] - 0.5
        else:
            x_max = max(plot_df["log2FC"].max(), 0)
            ax.set_xlim(-1, x_max + 1.5)
            x_text = x_max + 0.4

        # annotate
        for _, row in plot_df.iterrows():
            fc_mult = 2 ** row["log2FC"]
            label = f"{fc_mult:>4.1f}× | {row['expr']:>5.0f}"
            ax.text(
                x_text,
                row["gene"],
                label,
                va="center",
                ha="left",
                fontsize=7,
                fontfamily="monospace",
            )

    for j in range(i + 1, len(axes)):
        axes[j].axis("off")

    fig.suptitle(
        f"{module_name} across cell types\n"
        f"(label = fold-change × | {age_key.replace('_',' ')} ; p<{p_cutoff})",
        fontsize=13,
        y=1.02,
    )
    plt.tight_layout()

    # --- NEW: optional saving ---
    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
        base = os.path.join(save_dir, f"{module_name}_{age_suffix}")
        for ext in save_formats:
            fig.savefig(f"{base}.{ext}", bbox_inches="tight", dpi=dpi)
        print(f"✅ Saved figure to {base}.[{', '.join(save_formats)}]")

    if show:
        plt.show()
    else:
        plt.close(fig)

    return fig

def plot_region_bubbles(df_all, top_n_label=5):
    """
    df_all columns:
      region, gene, baseMean, log2FoldChange, pvalue
    """

    # clean / transform
    df_all = df_all.copy()
    df_all["neglog10p"] = -np.log10(df_all["pvalue"].clip(lower=1e-300))
    df_all["baseMean"] = df_all["baseMean"].astype(float)
    df_all["log2FoldChange"] = df_all["log2FoldChange"].astype(float)

    regions = df_all["region"].unique()
    n_regions = len(regions)

    ncols = 3
    nrows = int(np.ceil(n_regions / ncols))

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(5*ncols, 4*nrows),
        sharex=False, sharey=False
    )
    axes = np.atleast_1d(axes).flatten()

    for ax, region in zip(axes, regions):
        sub = df_all[df_all["region"] == region].copy()

        # pick which genes to label
        # rank = high log2FC and also decent baseMean
        sub["rank_metric"] = sub["log2FoldChange"].abs() * np.log10(sub["baseMean"]+1)
        sub = sub.sort_values("rank_metric", ascending=False)

        # scatter
        sc = ax.scatter(
            sub["log2FoldChange"],
            sub["baseMean"],
            c=sub["neglog10p"],
            s=np.clip(sub["baseMean"], 30, 400),  # size ~ abundance
            cmap="Reds",
            edgecolor="black",
            linewidth=0.4,
            alpha=0.8,
        )

        # annotate top_n_label
        for _, row in sub.head(top_n_label).iterrows():
            ax.text(
                row["log2FoldChange"],
                row["baseMean"],
                row["gene"],
                fontsize=8,
                ha="left",
                va="bottom"
            )

        ax.axvline(0, color="gray", lw=0.6, ls="--")
        ax.set_title(region, fontsize=11)
        ax.set_xlabel("log2FC (mtDSB vs control)")
        ax.set_ylabel("baseMean (avg expression)")

    # hide any empty axes
    for j in range(len(regions), len(axes)):
        axes[j].axis("off")

    # shared colorbar
    cbar = fig.colorbar(sc, ax=axes.tolist(), shrink=0.6)
    cbar.set_label("-log10(p value)")

    fig.suptitle("Regional driver genes: effect size vs abundance", fontsize=14)
    fig.tight_layout(rect=[0,0,1,0.95])
    plt.show()

def summarize_region_signature(df_region, top_n=5):
    """
    df_region: subset of df_all for one region
    returns:
      activation_score (float)
      signature_text (str)
    """
    # choose "important" genes: high positive log2FC, decent baseMean
    tmp = df_region.copy()
    tmp = tmp[tmp["log2FoldChange"] > 0]  # we care about upregulation here
    tmp["rank_metric"] = tmp["log2FoldChange"] * np.log10(tmp["baseMean"]+1)
    tmp = tmp.sort_values("rank_metric", ascending=False)

    sig_genes = []
    for _, row in tmp.head(top_n).iterrows():
        sig_genes.append(f"{row['gene']} (+{row['log2FoldChange']:.1f})")

    # simple activation score = sum of positive log2FC for the labeled ones
    activation_score = tmp.head(top_n)["log2FoldChange"].sum()

    return activation_score, ", ".join(sig_genes)

def plot_region_signature_bars_by_age(
    df_all_celltype,
    top_n_genes_per_region=5,
    figsize=(14,6),
    xlim=None,
    sort_regions_by="max_lfc"  # or "sum_lfc" etc.
):
    """
    df_all_celltype: output from build_df_all_by_age(...), but filtered to ONE cell_type.
                     must contain columns: ['gene','region','age','log2FoldChange','pvalue','baseMean']
    We'll generate 2 panels (age 21 and 60).
    """

    # helper: build plotting frame for one age
    def prep_one_age(df_age):
        if df_age.empty:
            return pd.DataFrame(columns=["region_gene","log2FoldChange","region","gene"])

        # rank genes within each region by log2FC
        df_ranked = (
            df_age
            .sort_values(["region","log2FoldChange"], ascending=[True, False])
            .groupby("region")
            .head(top_n_genes_per_region)
            .copy()
        )

        # optionally reorder regions by some score (like strongest single gene)
        region_strength = (
            df_ranked.groupby("region")["log2FoldChange"].max()
            if sort_regions_by == "max_lfc"
            else df_ranked.groupby("region")["log2FoldChange"].sum()
        )

        # sort regions strongest → weakest
        ordered_regions = region_strength.sort_values(ascending=False).index.tolist()
        df_ranked["region"] = pd.Categorical(df_ranked["region"], categories=ordered_regions, ordered=True)

        # build y labels as "Region: gene"
        df_ranked["region_gene"] = df_ranked.apply(
            lambda r: f"{r['region']}: {r['gene']}", axis=1
        )

        # final sort by region category order then log2FC
        df_ranked = df_ranked.sort_values(["region","log2FoldChange"], ascending=[True, False])
        return df_ranked

    # slice each age
    ages = ["21","60"]
    dfs_age = {age: prep_one_age(df_all_celltype[df_all_celltype["age"] == age]) for age in ages}

    # figure
    fig, axes = plt.subplots(1, 2, figsize=figsize, sharex=True)

    for ax, age in zip(axes, ages):
        dfp = dfs_age[age]

        if dfp.empty:
            ax.set_title(f"Age {age} (no data)")
            ax.axis("off")
            continue

        # plot bars
        y_positions = np.arange(len(dfp))[::-1]  # top-to-bottom
        ax.barh(
            y_positions,
            dfp["log2FoldChange"].values,
            color="crimson",
            edgecolor="black",
            linewidth=0.6,
        )

        # annotate baseMean (expression) on each bar
        for y0, (lfc, base, gene) in enumerate(
            zip(dfp["log2FoldChange"].values,
                dfp["baseMean"].values,
                dfp["gene"].values)
        ):
            ax.text(
                lfc + 0.1 if lfc >= 0 else lfc - 0.1,
                y_positions[y0],
                f"{lfc:+.2f}  | baseMean={base:.0f}",
                va="center",
                ha="left" if lfc >= 0 else "right",
                fontsize=7,
                fontfamily="monospace",
            )

        ax.set_yticks(y_positions)
        ax.set_yticklabels(dfp["region_gene"].values, fontsize=8)
        ax.set_xlabel("log₂FC (mtDSB vs control)")
        ax.set_title(f"Age {age}")

        ax.grid(axis="x", linestyle=":", alpha=0.4)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    # shared x limit
    if xlim is not None:
        for ax in axes:
            ax.set_xlim(xlim)
    else:
        # auto: get max abs value across both ages and pad
        max_abs = 0
        for age in ages:
            if not dfs_age[age].empty:
                max_abs = max(
                    max_abs,
                    np.nanmax(np.abs(dfs_age[age]["log2FoldChange"].values))
                )
        for ax in axes:
            ax.set_xlim(0, max_abs * 1.3 if max_abs > 0 else 1)

    fig.suptitle(
        f"{df_all_celltype['cell_type'].iloc[0]}: top induced genes by region and age",
        fontsize=14
    )
    fig.tight_layout(rect=[0,0,1,0.95])
    plt.show()

    return fig, axes

def combine_delta_lfc(res_a, res_b):
    lfc_col_a = [c for c in res_a.columns if c.startswith("log2FC_")][0]
    se_col_a  = [c for c in res_a.columns if c.startswith("lfcSE_")][0]
    lfc_col_b = [c for c in res_b.columns if c.startswith("log2FC_")][0]
    se_col_b  = [c for c in res_b.columns if c.startswith("lfcSE_")][0]

    df = pd.DataFrame(index=res_a.index.union(res_b.index))
    df[lfc_col_a] = res_a[lfc_col_a]
    df[se_col_a]  = res_a[se_col_a]
    df[lfc_col_b] = res_b[lfc_col_b]
    df[se_col_b]  = res_b[se_col_b]

    df = df.dropna(subset=[lfc_col_a, se_col_a, lfc_col_b, se_col_b]).copy()

    df["delta_log2FC"] = df[lfc_col_b] - df[lfc_col_a]
    df["se_delta"] = np.sqrt(np.square(df[se_col_a]) + np.square(df[se_col_b]))
    df["se_delta"] = df["se_delta"].replace(0, np.nan)
    df = df.dropna(subset=["se_delta"])
    from scipy.stats import norm
    df["z"] = df["delta_log2FC"] / df["se_delta"]
    df["p_delta"] = 2 * (1 - norm.cdf(np.abs(df["z"])))
    df["q_delta"] = multipletests(df["p_delta"], method="fdr_bh")[1]
    return df

def run_all(pb, condition_col="condition", outdir=None, n_cpus=8, save_csv=False):
    """
    For each cell type:
      - per-age simple effects (model vs control)
      - pairwise delta-LFC interaction proxy
      - gene_stats computed once and JOINED onto every table
    """
    if save_csv and outdir:
        os.makedirs(outdir, exist_ok=True)

    results = {}
    for ct, (counts_df, meta_df) in pb.items():
        # ensure alignment & compute gene stats once
        counts_df, meta_df = _align(counts_df, meta_df)
        gene_stats = _compute_gene_stats(counts_df, meta_df)

        ages = list(pd.unique(meta_df["age"].astype(str).str.strip()))
        per_age = {}

        # simple effects (join stats)
        for a in ages:
            try:
                res = simple_effect_at_age(counts_df, meta_df, a,
                                           condition_col=condition_col,
                                           ref_condition="control",
                                           n_cpus=n_cpus)
                res = res.join(gene_stats, how="left")
                per_age[a] = res
                if save_csv and outdir:
                    res.to_csv(os.path.join(outdir, f"{ct.replace(' ','_')}_simple_effect_age_{a}.csv"))
            except Exception as e:
                per_age[a] = None
                print(f"[{ct}] simple-effect failed at age={a}: {e}")

        # deltas (join stats)
        deltas = {}
        age_pairs = []
        for i in range(len(ages)):
            for j in range(i+1, len(ages)):
                a, b = ages[i], ages[j]
                if per_age.get(a) is None or per_age.get(b) is None:
                    continue
                dtab = combine_delta_lfc(per_age[a], per_age[b])
                dtab = dtab.join(gene_stats, how="left")
                dtab.attrs = {"age_A": a, "age_B": b,
                              "note": "delta_log2FC = LFC(age_B) - LFC(age_A)"}
                deltas[(a,b)] = dtab
                age_pairs.append((a,b))
                if save_csv and outdir:
                    dtab.to_csv(os.path.join(outdir, f"{ct.replace(' ','_')}_deltaLFC_{a}_vs_{b}.csv"))

        results[ct] = {
            "simple_effects": per_age,
            "deltas": deltas,
            "ages": ages,
            "age_pairs": age_pairs,
            "gene_stats": gene_stats,     # <--- stored for convenience
        }

    return results

def plot_lfc_across_ages_grid(
    results,
    ct,
    genes,
    q=0.05,
    ncols=4,
    figsize=(16, 9),
    sharey=True,
    fontsize_base=13
):
    """
    Plot log2FC across ages for multiple genes in subplots for a given cell type.
    Red dots = padj < q, grey = non-significant.
    Shared y-axis and increased text sizes for clarity.
    """
    ages = results[ct]["ages"]
    try:
        ages = sorted(ages, key=lambda x: float(x))
    except Exception:
        ages = sorted(ages)

    n_genes = len(genes)
    nrows = int(np.ceil(n_genes / ncols))

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=figsize,
        sharex=True,
        sharey=sharey,
        constrained_layout=True
    )
    axes = np.ravel(axes)

    for i, gene in enumerate(genes):
        ax = axes[i]
        lfc, padj = [], []

        for a in ages:
            res = results[ct]["simple_effects"].get(a)
            if res is not None and gene in res.index:
                lfc_col = next(c for c in res.columns if c.startswith("log2FC_"))
                padj_col = next(c for c in res.columns if c.startswith("padj_"))
                lfc.append(res.loc[gene, lfc_col])
                padj.append(res.loc[gene, padj_col])
            else:
                lfc.append(np.nan)
                padj.append(np.nan)

        sig_mask = np.array(padj) < q
        ax.plot(ages, lfc, lw=1.5, c="black", zorder=2)
        ax.scatter(
            ages, lfc,
            c=["red" if s else "grey" for s in sig_mask],
            s=60, zorder=3
        )
        ax.axhline(0, ls="--", c="0.4", lw=1.0, zorder=1)

        # larger fonts for readability
        ax.set_title(gene, fontsize=fontsize_base, pad=6)
        ax.tick_params(axis='x', labelrotation=45, labelsize=fontsize_base - 2)
        ax.tick_params(axis='y', labelsize=fontsize_base - 2)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    # remove unused subplots
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    # global labels — offset + larger text
    fig.suptitle(f"{ct}: log₂FC (mtDSB vs control) across ages",
                 fontsize=fontsize_base + 3, y=1.05, fontweight="bold")
    fig.text(0.5, -0.03, "Age (weeks)", ha="center", fontsize=fontsize_base + 1)
    fig.text(-0.04, 0.5, "log₂FC", va="center", rotation="vertical", fontsize=fontsize_base + 1)

    plt.subplots_adjust(bottom=0.12, left=0.10, right=0.98, top=0.93, wspace=0.25, hspace=0.35)
    plt.show()

def pick_top_per_age(
    res21, res60,
    lfc_col="log2FC_mtDSB_vs_control",
    padj_col="padj_mtDSB_vs_control",  # can also be "pvalue_mtDSB_vs_control"
    q=0.05,
    top_n_each=20,
    alpha_lfc=1.0,     # weight for |LFC|
    beta_p=0.0,        # weight for -log10(p); set 0 for pure LFC
    force_lfc_min=0.0, # require |LFC| >= this (0 disables)
    require_sig=False, # if True, require p<q for that age
    age_specific=False,# if True, drop genes significant in both ages
    max_labels=None    # optional cap on size of union 'picks'
):
    """
    Rank genes per age by a linear score:
        score_age = alpha_lfc * |LFC_age| + beta_p * (-log10(p_age))
    Then take top_n_each per age and return their union (+ ordering preserved).
    """

    # build combined table
    genes = res21.index.intersection(res60.index)
    A = res21.loc[genes, [lfc_col, padj_col]].rename(columns={lfc_col:"lfc21", padj_col:"p21"})
    B = res60.loc[genes, [lfc_col, padj_col]].rename(columns={lfc_col:"lfc60", padj_col:"p60"})
    df = A.join(B)
    # clean
    df = df.replace([np.inf, -np.inf], np.nan)
    df = df.dropna(subset=["lfc21","lfc60"])  # keep rows with both LFCs

    # choose p-values/FDRs; clip to avoid -log10(0)
    p21 = df["p21"].astype(float).clip(lower=1e-300)
    p60 = df["p60"].astype(float).clip(lower=1e-300)

    # basic masks
    sig21 = p21 < q
    sig60 = p60 < q
    lfc_ok_21 = df["lfc21"].abs() >= float(force_lfc_min)
    lfc_ok_60 = df["lfc60"].abs() >= float(force_lfc_min)

    # pools for ranking
    if require_sig:
        pool21 = df.loc[sig21 & lfc_ok_21].copy()
        pool60 = df.loc[sig60 & lfc_ok_60].copy()
    else:
        pool21 = df.loc[lfc_ok_21].copy()
        pool60 = df.loc[lfc_ok_60].copy()

    # fallbacks if empty
    if pool21.empty:
        pool21 = df.copy()
    if pool60.empty:
        pool60 = df.copy()

    # scores
    s21 = alpha_lfc * pool21["lfc21"].abs() + beta_p * (-np.log10(p21.loc[pool21.index]))
    s60 = alpha_lfc * pool60["lfc60"].abs() + beta_p * (-np.log10(p60.loc[pool60.index]))

    top21_idx = s21.sort_values(ascending=False).head(top_n_each).index
    top60_idx = s60.sort_values(ascending=False).head(top_n_each).index

    if age_specific:
        # keep genes uniquely significant in that age (by q)
        top21_idx = [g for g in top21_idx if (g in df.index) and (sig21.loc[g] and not sig60.loc[g])]
        top60_idx = [g for g in top60_idx if (g in df.index) and (sig60.loc[g] and not sig21.loc[g])]

    # union with order preserved: first 21, then add any new from 60
    picks_ordered = list(dict.fromkeys(list(top21_idx) + list(top60_idx)))

    if max_labels is not None and len(picks_ordered) > max_labels:
        # prioritize by max score across ages
        sc_all = pd.Series(0.0, index=df.index)
        sc_all.loc[pool21.index] = np.maximum(sc_all.loc[pool21.index], s21)
        sc_all.loc[pool60.index] = np.maximum(sc_all.loc[pool60.index], s60)
        picks_ordered = list(
            pd.Index(picks_ordered)
              .to_series()
              .loc[pd.Index(picks_ordered)]
              .sort_values(key=lambda idx: sc_all.loc[idx].values, ascending=False)
              .head(max_labels)
              .index
        )

    return picks_ordered, list(top21_idx), list(top60_idx)

def save_all_results(results, outdir="../../results"):
    """
    Save all DE results from the nested `results` dict into organized subfolders.
    Structure:
        outdir/
          └── <cell_type>/
                ├── simple_effect_age_<age>.csv
                ├── deltaLFC_<ageA>_vs_<ageB>.csv
    """
    os.makedirs(outdir, exist_ok=True)

    for ct, rdict in results.items():
        ct_dir = os.path.join(outdir, ct.replace(" ", "_"))
        os.makedirs(ct_dir, exist_ok=True)

        # --- Simple effects (mtDSB vs control per age) ---
        for age, df in rdict.get("simple_effects", {}).items():
            if df is not None:
                fname = f"simple_effect_age_{age}.csv"
                df.to_csv(os.path.join(ct_dir, fname))

        # --- Delta-LFC (interaction proxy between ages) ---
        for (a1, a2), df in rdict.get("deltas", {}).items():
            if df is not None:
                fname = f"deltaLFC_{a1}_vs_{a2}.csv"
                df.to_csv(os.path.join(ct_dir, fname))

        print(f"✅ Saved {ct} results to {ct_dir}")
