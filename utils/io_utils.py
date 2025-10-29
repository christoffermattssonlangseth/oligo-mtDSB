"""
Auto-generated utilities for io_utils.
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

def filter_low_genes(df_long, min_mean=MIN_MEAN, min_detect=MIN_DET):
    kept = []
    for ct, dfc in df_long.groupby("cell_class"):
        stats = (dfc
                 .groupby("gene")
                 .agg(mean_expr=("value","mean"),
                      n_detect=("value", lambda x: (x>0).sum()))
                 .reset_index())
        good = stats[(stats.mean_expr >= min_mean) & (stats.n_detect >= min_detect)].gene
        kept.append(dfc[dfc.gene.isin(good)])
    return pd.concat(kept, ignore_index=True)

def fit_many_genes_condition(df_long, genes, celltypes,
                             *, advi=False, draws=800, tune=800, chains=2, cores=2,
                             target_accept=0.9, seed=42):
    rows = []
    for ct in celltypes:
        df_ct = df_long[df_long["cell_class"] == ct]
        for g in genes:
            sub = df_ct[df_ct["gene"] == g]
            if sub.empty or sub["condition"].nunique() < 2:
                continue
            idata, meta = fit_gene_pymc_condition(
                sub, advi=advi, draws=draws, tune=tune, chains=chains, cores=cores,
                target_accept=target_accept, seed=seed
            )
            eff = summarize_condition_effect(idata, meta)
            # keep only the mtDSB row if present
            eff = eff[eff["condition"].str.contains("mtDSB", case=False, na=False)]
            eff.insert(0, "gene", g)
            eff.insert(0, "cell_class", ct)
            rows.append(eff)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()

def fit_many_genes_condition_parallel(df_long, genes, celltypes, n_jobs=8, **kwargs):
    out = []
    for ct in celltypes:
        res = Parallel(n_jobs=n_jobs, backend="loky")(
            delayed(_fit_one_gene_cond)(
                df_long, g, ct,
                advi=True,
                draws=kwargs.get("draws",600),
                tune=0, chains=1, cores=1,             # ADVI ignores cores
                target_accept=kwargs.get("target_accept",0.9),
                aggregate_by_sample=kwargs.get("aggregate_by_sample",True),
                sample_col=kwargs.get("sample_col","sample"),   # set to "sample_id" if needed
                progressbar=False
            )
            for g in genes
        )
        ct_df = [r for r in res if r is not None]
        if ct_df:
            out.append(pd.concat(ct_df, ignore_index=True))
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()

def spatial_neighbourhoods_fast(
    ad,
    cluster_label="leiden_0.5",
    max_distance=200.0,
    x_col="x_centroid",
    y_col="y_centroid",
):
    """
    Returns an AnnData where X is (cells x clusters) with counts of neighbors
    within `max_distance` per cluster. Uses sparse ops (fast + memory-light).
    """
    # --- coordinates ---
    coords = ad.obs[[x_col, y_col]].to_numpy()

    # Sparse adjacency (N x N) with 1 if within radius, excluding self
    A = radius_neighbors_graph(
        coords, radius=max_distance, mode="connectivity",
        include_self=False, n_jobs=-1
    ).tocsr()

    # --- cluster one-hot (N x K) ---
    cats = ad.obs[cluster_label].astype("category")
    cat_codes = cats.cat.codes.to_numpy()
    N = ad.n_obs
    K = len(cats.cat.categories)
    rows = np.arange(N, dtype=int)
    C = sp.csr_matrix((np.ones(N, dtype=np.int32), (rows, cat_codes)), shape=(N, K))

    # --- neighbor counts per cluster: F = A @ C  (N x K) ---
    F = A @ C  # integer counts

    # Build output AnnData
    out = sc.AnnData(
        X=F,
        obs=ad.obs.copy(),
        var=pd.DataFrame(index=[str(c) for c in cats.cat.categories])
    )
    # Keep the adjacency if you want it later (in obsp, not layers)
    out.obsp["neighbors_within_radius"] = A
    out.uns["spatial_neighbourhoods"] = {
        "cluster_label": cluster_label,
        "max_distance": float(max_distance),
        "x_col": x_col,
        "y_col": y_col,
    }
    return out

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

def _to_anndata(counts_df, meta_df):
    X = counts_df.T.astype(np.int64)  # samples x genes
    ad = AnnData(X=X.values)
    ad.obs_names = X.index.astype(str)
    ad.var_names = counts_df.index.astype(str)
    ad.obs = meta_df.copy()
    return ad

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

def extract_top_genes(
    results_by_region,
    cell_class,
    regions=None,
    ages=None,
    padj_thr=0.05,
    lfc_thr=1.0,
    top_n=50
):
    rows = []
    res_cl = results_by_region[cell_class]
    if regions is None:
        regions = list(res_cl.keys())

    for reg in regions:
        reg_data = res_cl.get(reg, {})
        effects = reg_data.get("simple_effects", {})
        if not effects:
            continue

        # normalize available age keys to strings
        avail_age_keys = {str(k): k for k in effects.keys()}
        ages_to_use = [str(a) for a in (ages if ages is not None else effects.keys())]

        for age_str in ages_to_use:
            key = avail_age_keys.get(age_str, None)
            if key is None:
                continue
            res = effects.get(key, None)
            if res is None or res.empty:
                continue

            # standard DESeq2/PyDESeq2 columns
            if ("padj" not in res.columns) or ("log2FoldChange" not in res.columns):
                # skip if table doesn't look like a DE result
                continue

            sig = res[(res["padj"] < padj_thr) & (res["log2FoldChange"].abs() >= lfc_thr)].copy()
            if sig.empty:
                continue

            sig["region"] = reg
            sig["age"] = age_str
            sig["neglog10FDR"] = -np.log10(sig["padj"].clip(lower=1e-300))
            sig.rename(columns={"log2FoldChange": "log2FC"}, inplace=True)

            rows.append(sig[["log2FC", "padj", "neglog10FDR", "region", "age"]])

    if not rows:
        print("No genes passed the thresholds; try relaxing padj_thr or lfc_thr.")
        return pd.DataFrame()

    df = pd.concat(rows).sort_values("neglog10FDR", ascending=False)
    return df.head(top_n)

def load_and_pad(paths, bg="white", scale=0.5):
    imgs, sizes = [], []
    for p in paths:
        im = Image.open(p).convert("RGB")
        # ↓ downscale early for smaller memory footprint
        if scale != 1.0:
            w, h = im.size
            im = im.resize((int(w * scale), int(h * scale)), Image.Resampling.LANCZOS)
        imgs.append(im)
        sizes.append(im.size)

    max_w = max(w for w, h in sizes)
    max_h = max(h for w, h in sizes)

    padded = []
    for im in imgs:
        w, h = im.size
        canvas = Image.new("RGB", (max_w, max_h), color=bg)
        canvas.paste(im, ((max_w - w)//2, (max_h - h)//2))
        padded.append(np.asarray(canvas))
    return padded

def fit_many_genes(df_long, genes, celltypes,
                   *, celltype_col="cell_class",
                   advi=False, draws=1000, tune=1000, chains=4,
                   target_accept=0.9, seed=42):
    rows = []
    for ct in celltypes:
        df_ct = df_long[df_long[celltype_col] == ct]
        for g in genes:
            sub = df_ct[df_ct["gene"] == g]
            if sub.empty:
                continue
            idata, meta = fit_gene_pymc_ac(sub, advi=advi, draws=draws, tune=tune,
                                           chains=chains, target_accept=target_accept, seed=seed)
            eff = summarize_effects_ac(idata, meta)
            eff.insert(0, "gene", g)
            eff.insert(0, celltype_col, ct)   # keep the same column name
            rows.append(eff)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()

def run_shortlisted_fits(df_filtered, shortlist_csv, out_dir="./",
                         *, celltype_col="cell_class",
                         advi=False, draws=1000, tune=1000, chains=4,
                         target_accept=0.9, seed=42):
    from pathlib import Path
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)

    shortlists = load_shortlists(shortlist_csv)
    shortlists = shortlist_present(df_filtered, shortlists, celltype_col=celltype_col)

    all_summaries = []
    for ct, genes in shortlists.items():
        if not genes:
            continue
        print(f"\nFitting {ct} — {len(genes)} shortlisted genes")

        sct = fit_many_genes(df_filtered, genes, [ct],
                             celltype_col=celltype_col,
                             advi=advi, draws=draws, tune=tune, chains=chains,
                             target_accept=target_accept, seed=seed)
        if not sct.empty:
            sct.to_csv(out_dir / f"bayes_mtDSB_{ct.replace(' ', '_')}_nosex.csv",
                       index=False)
            all_summaries.append(sct)

    if all_summaries:
        all_df = pd.concat(all_summaries, ignore_index=True)
        all_df.to_csv(out_dir / "bayes_mtDSB_ALL_nosex.csv", index=False)
        return all_df
    return pd.DataFrame()

def plot_interaction_caterpillar(df, cell_type, top=25, save=None):
    d = df[(df["cell_class"]==cell_type) & (df["effect_clean"]=="age interaction (60–21)")].copy()
    if d.empty:
        print(f"No interaction rows for {cell_type}")
        return

    d = d.sort_values("mean", ascending=False)

    # build 'keep' safely without duplicates
    half = max(1, top // 2)
    keep = pd.concat([d.head(half), d.tail(half)]).drop_duplicates(subset=["gene"]).sort_values("mean")

    # order genes explicitly (avoid pandas Categorical)
    gene_order = keep["gene"].tolist()

    plt.figure(figsize=(9, max(6, 0.35*len(keep))))
    ax = sns.pointplot(data=keep, y="gene", x="mean", order=gene_order, join=False, color="C0")
    for _, r in keep.iterrows():
        plt.plot([r["hdi_2.5%"], r["hdi_97.5%"]], [r["gene"], r["gene"]], color="C0", lw=2, alpha=.7)
    plt.axvline(0, ls="--", c="k", lw=1)
    plt.title(f"{cell_type} — Age interaction (60–21): effect size with 95% HDIs")
    plt.xlabel("Posterior mean (log scale)"); plt.ylabel("gene")
    plt.tight_layout()
    if save: plt.savefig(save, dpi=200, bbox_inches="tight")
    plt.show()

def plot_condition_caterpillar(df, cell_type, top=30):
    d = df[df["cell_class"]==cell_type].copy()
    if d.empty:
        print(f"No rows for {cell_type}")
        return

    d = d[d["condition"].str.contains("mtDSB", case=False)]
    d = d.sort_values("mean", ascending=False)
    keep = pd.concat([d.head(top//2), d.tail(top//2)]).drop_duplicates(subset="gene").sort_values("mean")

    plt.figure(figsize=(9, max(6, 0.35*len(keep))))
    sns.pointplot(data=keep, y="gene", x="mean", join=False, color="C0")
    for _, r in keep.iterrows():
        plt.plot([r["hdi_2.5%"], r["hdi_97.5%"]], [r["gene"], r["gene"]], color="C0", lw=2, alpha=.7)
    plt.axvline(0, ls="--", c="k", lw=1)
    plt.title(f"{cell_type} — mtDSB effect (condition-only model)")
    plt.xlabel("Posterior mean (log scale)"); plt.ylabel("gene")
    plt.tight_layout()
    plt.show()

def format_data_neighs_colapse(
    adata: AnnData,
    spatial_key: str = "spatial",
    neighs: int = 10
) -> AnnData:
    """Pseudobin each cell’s expression by summing over its spatial neighbors."""
    if spatial_key not in adata.obsm:
        coords = adata.obs[['X','Y']].values.astype(float)
        adata.obsm[spatial_key] = coords

    ad = adata.copy()
    print(f"    → Computing spatial neighbors (k={neighs})…", end="", flush=True)
    sq.gr.spatial_neighbors(ad, spatial_key=spatial_key, n_neighs=neighs)
    print("done")

    C = ad.obsp['spatial_connectivities'].astype(bool).astype(int)
    C = C + sp.eye(C.shape[0], format='csr')
    X = ad.X.copy() if not isinstance(ad.X, np.ndarray) else ad.X
    Xc = C.dot(X)

    print("    → Pseudobinning complete")
    return AnnData(
        X    = Xc,
        obs  = ad.obs.copy(),
        var  = ad.var.copy(),
        obsm = ad.obsm.copy(),
        uns  = ad.uns.copy()
    )

def domains_by_rbd(
    adata: AnnData,
    hyperparams: dict,
    sample_key: str = 'sample',
    cell_id_key: str = None,
    checkpoint_dir: str = './checkpoints',
    force_rerun: bool = False
) -> (AnnData, AnnData):
    """
    Read-based domains with checkpointed saves, force-rerun option, and printouts.
    """
    os.makedirs(checkpoint_dir, exist_ok=True)

    # A) ensure unique IDs
    if cell_id_key is None:
        print("1) Ensuring unique cell IDs…", end="", flush=True)
        adata.obs_names_make_unique()
        cell_id_key = 'cell_id'
        adata.obs[cell_id_key] = adata.obs_names
        print("done\n")

    def load_or_compute(path, compute_fn, step_name):
        if (not force_rerun) and os.path.exists(path):
            print(f"{step_name}: Loading from checkpoint…", end="", flush=True)
            ad = sc.read_h5ad(path)
            print("loaded\n")
            return ad, False
        else:
            print(f"{step_name}: Computing…", flush=True)
            ad = compute_fn()
            ad.write_h5ad(path)
            print(f"{step_name}: Saved to {path}\n")
            return ad, True

    # 1) PSEUDOBINNING
    pseudobin_path = os.path.join(checkpoint_dir, 'pseudobinned.h5ad')
    def compute_pseudobin():
        parts = []
        for samp in adata.obs[sample_key].unique():
            print(f"  Sample '{samp}':")
            sub = adata[adata.obs[sample_key] == samp].copy()
            parts.append(format_data_neighs_colapse(
                sub, spatial_key='spatial', neighs=hyperparams['neighbors']
            ))
        return sc.concat(parts, merge='same')

    pseudobinned, _ = load_or_compute(
        pseudobin_path, compute_pseudobin, "Step 1) Pseudobinning"
    )

    # 2) FILTERING
    filtered_path = os.path.join(checkpoint_dir, 'filtered.h5ad')
    def compute_filtered():
        ad_neigh = pseudobinned.copy()
        ad_neigh.X = np.nan_to_num(ad_neigh.X)
        sc.pp.filter_cells(ad_neigh, min_counts=4)
        ad_neigh.raw = ad_neigh
        return ad_neigh

    adata_neigh, _ = load_or_compute(
        filtered_path, compute_filtered, "Step 2) Filtering"
    )

    # 3) UMAP
    umap_path = os.path.join(checkpoint_dir, 'umap.h5ad')
    def compute_umap():
        nn  = hyperparams.get('n_neighbors', 20)
        pcs = hyperparams.get('n_pcs', None)
        if pcs is None:
            sc.pp.neighbors(adata_neigh, n_neighbors=nn, use_rep='X')
        else:
            sc.tl.pca(adata_neigh, n_comps=pcs, svd_solver='arpack')
            sc.pp.neighbors(adata_neigh, n_neighbors=nn, n_pcs=pcs)
        sc.tl.umap(adata_neigh, min_dist=hyperparams.get('min_dist',0.1))
        return adata_neigh

    adata_neigh, _ = load_or_compute(
        umap_path, compute_umap, "Step 3) UMAP"
    )

    # 4) CLUSTERING at one or many resolutions
    algo = hyperparams['clustering_algorithm'].lower()
    res_list = hyperparams.get('resolutions',
               [hyperparams.get('resolution')])

    for res in res_list:
        key = f"rbd_domain_{res}"
        chk = os.path.join(checkpoint_dir, f"clustered_{res}.h5ad")

        def compute_cluster():
            if algo == 'leiden':
                sc.tl.leiden(adata_neigh, resolution=res, key_added=key)
            else:
                sc.tl.louvain(adata_neigh, resolution=res, key_added=key)
            return adata_neigh

        adata_neigh, _ = load_or_compute(
            chk, compute_cluster,
            f"Step 4) Clustering (res={res})"
        )

        # map back to original
        print(f"  → Mapping rbd_domain_{res} back to original…", end="", flush=True)
        mapping = dict(zip(
            adata_neigh.obs[cell_id_key],
            adata_neigh.obs[key]
        ))
        adata.obs[key] = adata.obs[cell_id_key].map(mapping).astype(str)
        print("done\n")

    print("All steps complete! Returning results.")
    return adata, adata_neigh

def compartment_composition(
    adata,
    compartment,
    anno_key="rbd_annotation",
    groupby="cell_class",
    sort=True
) -> pd.DataFrame:
    """
    Return counts and percentages of `groupby` within one compartment.

    Parameters
    ----------
    adata : AnnData
    compartment : str                      # e.g. 'Meningeal–parenchymal border domain'
    anno_key : str                         # column in .obs that holds the compartment labels
    groupby : str                          # column in .obs to summarize (e.g. cell_class)
    sort : bool                            # sort by descending percentage

    Returns
    -------
    DataFrame with columns: ['count', 'percent'] indexed by group
    """
    mask = (adata.obs[anno_key] == compartment)
    sub = adata.obs.loc[mask, groupby].value_counts()
    total = int(sub.sum())
    df = pd.DataFrame({"count": sub.astype(int)})
    df["percent"] = df["count"] / total * 100.0
    if sort:
        df = df.sort_values("percent", ascending=False)
    df.attrs["total"] = total
    df.attrs["compartment"] = compartment
    df.attrs["groupby"] = groupby
    return df

def compare_compartments_heatmap(
    adata,
    compartments,
    anno_key="rbd_annotation",
    groupby="cell_class",
    top_k=None,
    figsize=(10, 6),
    cmap="viridis"
):
    """
    Heatmap of percentages across multiple compartments.
    If `top_k` is set, uses the globally most abundant K classes across all chosen compartments.
    """
    mats = []
    for comp in compartments:
        df = compartment_composition(adata, comp, anno_key=anno_key, groupby=groupby, sort=False)
        mats.append(df["percent"])

    M = pd.concat(mats, axis=1)
    M.columns = list(compartments)
    M = M.fillna(0.0)

    if top_k is not None and top_k < M.shape[0]:
        order = M.sum(1).sort_values(ascending=False).index[:top_k]
        M = M.loc[order]

    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(M.values, aspect="auto", cmap=cmap)
    ax.set_xticks(range(M.shape[1]))
    ax.set_xticklabels(M.columns, rotation=0)
    ax.set_yticks(range(M.shape[0]))
    ax.set_yticklabels(M.index)
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("% within compartment")
    ax.set_title(f"{groupby} composition across compartments")
    plt.tight_layout()
    return ax, M

def pseudobulk_to_long_with_celltype(pb_norm, groups_df, var_names,
                                     sample_col="sample_id", age_col="age",
                                     cond_col="condition", celltype_col="cell_class",
                                     value_name="value"):
    df_wide = pd.DataFrame(pb_norm, columns=var_names)
    df = pd.concat([groups_df.reset_index(drop=True), df_wide], axis=1)

    # Ensure these exist
    for col in [sample_col, age_col, cond_col, celltype_col]:
        if col not in df.columns:
            raise KeyError(f"Missing column in groups_df: {col}")

    # Melt ONLY the gene columns (var_names)
    long = df.melt(
        id_vars=[sample_col, age_col, cond_col, celltype_col],
        value_vars=list(var_names),               # <- critical
        var_name="gene",
        value_name=value_name
    )

    # Standardize column names
    long = long.rename(columns={
        sample_col: "sample",
        age_col: "age",
        cond_col: "condition",
        celltype_col: "cell_type"
    })

    # Coerce numeric values and drop rows where value is NA
    long[value_name] = pd.to_numeric(long[value_name], errors="coerce")
    long = long.dropna(subset=[value_name])

    return long

def fit_all_genes(df_long, genes=None, draws=1200, tune=1200, target_accept=0.9):
    if genes is None:
        genes = sorted(df_long["gene"].unique().tolist())
    out = []
    for g in genes:
        sub = df_long[df_long["gene"] == g]
        # require both conditions present overall
        labels = sub["condition"].astype(str).str.lower().unique().tolist()
        if not any(l in ("control", "ctrl") for l in labels) or "mtdsb" not in labels:
            continue
        idata, D = fit_gene_pymc(sub, draws=draws, tune=tune, target_accept=target_accept, seed=42)
        out.append(summarize_gene(idata, D, g))
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()

def filter_genes_in_data(data, genes):
    """
    Return only genes from the input list that are present in the given data.

    Parameters
    ----------
    data : AnnData | pandas.DataFrame
        The data object containing gene names either in:
        - adata.var_names (for AnnData)
        - df.index (for DataFrame)
    genes : list of str
        List of gene names to filter.

    Returns
    -------
    list of str
        Filtered list containing only genes found in the data.
    """
    import pandas as pd
    from anndata import AnnData

    # Determine available gene names
    if isinstance(data, AnnData):
        available = set(data.var_names)
        src = "adata.var_names"
    elif isinstance(data, pd.DataFrame):
        available = set(data.index)
        src = "df.index"
    else:
        raise TypeError("Input must be AnnData or pandas DataFrame")

    found = [g for g in genes if g in available]
    missing = [g for g in genes if g not in available]

    print(f"✅ {len(found)} genes found in {src}.")
    if missing:
        print(f"⚠️ {len(missing)} genes not found and removed: {', '.join(missing)}")

    return found

def _load_transcripts_polars(files, nucleus_only, allowed_categories):
    """
    Fast path: read transcripts with polars scan_parquet over one or many files.
    Returns pandas DataFrame with ['cell_id','feature_name'].
    """
    import polars as pl

    lazy_list = []
    for f in files:
        lf = pl.scan_parquet(f).select(
            "cell_id",
            "feature_name",
            "codeword_category",
            "overlaps_nucleus",
        )
        lazy_list.append(lf)

    lf_all = pl.concat(lazy_list)

    filt = (pl.col("codeword_category").is_in(list(allowed_categories)))
    if nucleus_only:
        filt = filt & (pl.col("overlaps_nucleus") == 1)

    lf_filt = (
        lf_all
        .filter(filt)
        .select(["cell_id", "feature_name"])
    )

    df = lf_filt.collect().to_pandas()
    df = df.dropna(subset=["cell_id", "feature_name"])
    return df

def _load_transcripts_pyarrow(files, nucleus_only, allowed_categories):
    """
    Fallback: use pyarrow to read + concat shards.
    Returns pandas DataFrame with ['cell_id','feature_name'].
    """
    import pyarrow.parquet as pq
    import pyarrow as pa

    cols = ["cell_id", "feature_name", "codeword_category", "overlaps_nucleus"]
    tables = []
    for f in files:
        tbl = pq.read_table(f, columns=cols)
        tables.append(tbl)

    big_tbl = pa.concat_tables(tables, promote=True)
    df = big_tbl.to_pandas()

    mask = df["codeword_category"].isin(list(allowed_categories))
    if nucleus_only:
        mask = mask & (df["overlaps_nucleus"] == 1)

    df = df.loc[mask, ["cell_id", "feature_name"]].copy()
    df = df.dropna(subset=["cell_id", "feature_name"])
    return df

def _load_transcripts_pandas(files, nucleus_only, allowed_categories):
    """
    Last-resort fallback: pandas.read_parquet on each file.
    Returns pandas DataFrame with ['cell_id','feature_name'].
    """
    dfs = []
    for f in files:
        tdf = pd.read_parquet(f)
        # filter
        good = tdf["codeword_category"].isin(list(allowed_categories))
        if nucleus_only and "overlaps_nucleus" in tdf.columns:
            good = good & (tdf["overlaps_nucleus"] == 1)
        tdf = tdf.loc[good, ["cell_id", "feature_name"]].copy()
        tdf = tdf.dropna(subset=["cell_id", "feature_name"])
        dfs.append(tdf)

    if len(dfs) == 0:
        return pd.DataFrame(columns=["cell_id","feature_name"])
    return pd.concat(dfs, axis=0, ignore_index=True)

def _read_transcripts_fast(transcripts_path, nucleus_only, allowed_categories):
    """
    Robust loader for transcripts.
    Tries Polars -> PyArrow -> Pandas.
    Returns pandas DF with ['cell_id','feature_name'].
    """
    files = _list_parquet_files(transcripts_path)

    last_err = None
    # try polars
    try:
        df = _load_transcripts_polars(files, nucleus_only, allowed_categories)
        if df.shape[0] > 0:
            return df
        print(f"[warn] polars returned 0 rows from {transcripts_path}, trying pyarrow...")
    except Exception as e:
        last_err = e
        print(f"[warn] polars failed on {transcripts_path}: {e}")

    # try pyarrow
    try:
        df = _load_transcripts_pyarrow(files, nucleus_only, allowed_categories)
        if df.shape[0] > 0:
            return df
        print(f"[warn] pyarrow returned 0 rows from {transcripts_path}, trying pandas...")
    except Exception as e:
        last_err = e
        print(f"[warn] pyarrow failed on {transcripts_path}: {e}")

    # fallback pandas
    df = _load_transcripts_pandas(files, nucleus_only, allowed_categories)
    if df.shape[0] == 0:
        raise RuntimeError(
            f"no transcripts passed filters in {transcripts_path}; last_err={last_err}"
        )
    return df

def _load_cells_any(cells_path, x_key, y_key, extra_obs_cols):
    """
    Load cells.parquet quickly.
    We'll try pyarrow w/ column subset first; fallback to pandas.
    Returns a pandas DF indexed by cell_id.
    """
    cols = ["cell_id", x_key, y_key, *extra_obs_cols]
    # dedupe but keep order
    cols = list(dict.fromkeys(cols))

    # try pyarrow
    try:
        import pyarrow.parquet as pq
        tbl = pq.read_table(cells_path, columns=cols)
        cells_df = tbl.to_pandas()
    except Exception as e:
        print(f"[warn] pyarrow cells loader failed: {e}")
        # pandas fallback
        full_df = pd.read_parquet(cells_path)
        # keep only cols that actually exist
        use_cols = [c for c in cols if c in full_df.columns]
        cells_df = full_df[use_cols].copy()

    if "cell_id" not in cells_df.columns:
        raise KeyError(f"'cell_id' not found in {cells_path}")

    cells_df = cells_df.set_index("cell_id")
    return cells_df

def load_xenium_to_anndata(
    transcripts_path,
    cells_path,
    min_counts_per_cell=100,
    min_cells_per_gene=100,
    nucleus_only=True,
    allowed_categories=("predesigned_gene", "custom_gene"),
    x_key="x_centroid",
    y_key="y_centroid",
    extra_obs_cols=("nucleus_area",),
    layer_name="counts",
):
    """
    Build an AnnData from raw Xenium parquet exports for ONE sample.
    Uses fast loaders (polars/pyarrow) with graceful fallback.
    """

    # ---------------------------
    # 1. read / filter transcripts
    # ---------------------------
    df_tx = _read_transcripts_fast(
        transcripts_path,
        nucleus_only=nucleus_only,
        allowed_categories=allowed_categories,
    )

    if df_tx.empty:
        raise ValueError(f"No transcripts passed filters in {transcripts_path}")

    # ---------------------------
    # 2. count transcripts per (cell_id × gene)
    # ---------------------------
    counts = (
        df_tx.groupby(["cell_id", "feature_name"])
        .size()
        .reset_index(name="count")
    )

    count_matrix = (
        counts.pivot(index="cell_id", columns="feature_name", values="count")
        .fillna(0)
        .astype(int)
    )

    # sparse CSR
    X_sparse = sparse.csr_matrix(count_matrix.values)

    # ---------------------------
    # 3. create AnnData
    # ---------------------------
    adata = ad.AnnData(X=X_sparse)

    adata.obs_names = count_matrix.index.astype(str)
    adata.var_names = count_matrix.columns.astype(str)

    adata.obs["cell_id"] = adata.obs_names
    adata.var["gene"] = adata.var_names

    # ---------------------------
    # 4. attach spatial / morphology metadata
    # ---------------------------
    cells_df = _load_cells_any(cells_path, x_key, y_key, extra_obs_cols)

    for col in [x_key, y_key, *extra_obs_cols]:
        if col in cells_df.columns:
            adata.obs[col] = adata.obs_names.map(cells_df[col])
        else:
            adata.obs[col] = np.nan

    # spatial coords for squidpy/plotting
    if x_key in adata.obs.columns and y_key in adata.obs.columns:
        adata.obsm["spatial"] = adata.obs[[x_key, y_key]].to_numpy()

    # ---------------------------
    # 5. stash raw counts layer
    # ---------------------------
    adata.layers[layer_name] = adata.X.copy()

    # ---------------------------
    # 6. QC + basic filters
    # ---------------------------
    sc.pp.calculate_qc_metrics(
        adata,
        inplace=True,
        percent_top=None,
        log1p=False,
    )

    sc.pp.filter_cells(adata, min_counts=min_counts_per_cell)
    sc.pp.filter_genes(adata, min_cells=min_cells_per_gene)

    return adata

def load_all_xenium_runs(
    base_dir,
    min_counts_per_cell=100,
    min_cells_per_gene=100,
    nucleus_only=True,
    allowed_categories=("predesigned_gene", "custom_gene"),
    x_key="x_centroid",
    y_key="y_centroid",
    extra_obs_cols=("nucleus_area",),
    layer_name="counts",
):
    """
    Walk `base_dir`, find each sample subfolder, build AnnData,
    and concatenate them.

    Returns
    -------
    adata_all : AnnData (concatenated union-of-genes, CSR sparse)
    adata_list : list[AnnData] (per sample)
    """

    adata_list = []

    for sample in os.listdir(base_dir):
        sample_path = os.path.join(base_dir, sample)
        if not os.path.isdir(sample_path):
            continue  # skip stray files

        transcripts_path = os.path.join(sample_path, "transcripts.parquet")
        cells_path = os.path.join(sample_path, "cells.parquet")

        # tolerate either file or dir for transcripts
        # cells is usually 1 parquet file
        if not (os.path.exists(transcripts_path) or os.path.isdir(transcripts_path)):
            continue
        if not os.path.exists(cells_path):
            continue

        print(f"→ loading {sample} ...")
        try:
            _adata_ = load_xenium_to_anndata(
                transcripts_path=transcripts_path,
                cells_path=cells_path,
                min_counts_per_cell=min_counts_per_cell,
                min_cells_per_gene=min_cells_per_gene,
                nucleus_only=nucleus_only,
                allowed_categories=allowed_categories,
                x_key=x_key,
                y_key=y_key,
                extra_obs_cols=extra_obs_cols,
                layer_name=layer_name,
            )

            # guess sample_id from folder name (your old split('__')[2])
            parts = sample.split("__")
            if len(parts) >= 3:
                sample_id = parts[2]
            else:
                sample_id = sample
            _adata_.obs["sample_id"] = sample_id

            adata_list.append(_adata_)

        except Exception as e:
            print(f"[warn] skipping {sample}: {e}")
            continue

    if len(adata_list) == 0:
        raise RuntimeError(f"No valid samples loaded from {base_dir}")

    # concat (outer join on vars; fill missing with 0)
    adata_all = ad.concat(
        adata_list,
        join="outer",
        label="sample_id",
        keys=[a.obs["sample_id"].unique()[0] for a in adata_list],
        fill_value=0,
    )

    # enforce sparse CSR after concat (scanpy may densify)
    if not sparse.issparse(adata_all.X):
        adata_all.X = sparse.csr_matrix(adata_all.X)

    return adata_all, adata_list

def load_ensembl_gene_map(gtf_path):
    """
    Parse an Ensembl GTF and return a dict: {ensembl_gene_id: gene_symbol}.
    Keeps only 'gene' features; strips version suffixes.
    """
    gene_map = {}
    pat_id   = re.compile(r'gene_id "([^"]+)"')
    pat_name = re.compile(r'gene_name "([^"]+)"')
    with open(gtf_path, "r") as fh:
        for line in fh:
            if line.startswith("#"): 
                continue
            # Only keep 'gene' feature lines to avoid huge memory use
            # GTF columns: seqname, source, feature, start, end, score, strand, frame, attributes
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            attrs = parts[8]
            m_id = pat_id.search(attrs)
            m_nm = pat_name.search(attrs)
            if not (m_id and m_nm):
                continue
            gid = m_id.group(1).split('.')[0]   # strip version, e.g. ENSMUSG... .1 → base
            gnm = m_nm.group(1)
            # keep first occurrence (usually fine). If you want, prefer protein_coding using parts[1]/attrs
            if gid not in gene_map:
                gene_map[gid] = gnm
    return gene_map

def build_df_all_by_age(results_by_region, cell_type):
    """
    For a given cell type (e.g. 'Oligodendrocytes'),
    grab top genes per region for each age ('21','60'), and
    return one tidy DataFrame with columns:
    ['gene','baseMean','log2FoldChange','pvalue','region','cell_type','age']
    """
    rows = []
    for region, payload in results_by_region[cell_type].items():
        # walk over both ages if present
        for age in ["21", "60"]:
            try:
                df_age = (
                    payload["simple_effects"][age]
                    .sort_values(by="log2FoldChange", ascending=False)
                    .head(10)[["baseMean", "log2FoldChange", "pvalue"]]
                    .copy()
                )
            except Exception:
                continue  # skip if that age / region combo doesn't exist

            df_age = df_age.reset_index()  # gene out of index
            df_age = df_age.rename(columns={"index": "gene"})
            df_age["region"] = region
            df_age["cell_type"] = cell_type
            df_age["age"] = age
            rows.append(df_age)

    if len(rows) == 0:
        return pd.DataFrame(columns=[
            "gene","baseMean","log2FoldChange","pvalue","region","cell_type","age"
        ])

    out = pd.concat(rows, ignore_index=True)

    # enforce numeric types
    out["baseMean"] = pd.to_numeric(out["baseMean"], errors="coerce")
    out["log2FoldChange"] = pd.to_numeric(out["log2FoldChange"], errors="coerce")
    out["pvalue"] = pd.to_numeric(out["pvalue"], errors="coerce")
    out = out[out.baseMean > 300]
    return out

def get_DEGs(adata, n_genes=100):
    res = adata.uns['rank_genes_groups']
    groups = res['names'].dtype.names
    dfs = []
    for g in groups:
        df = pd.DataFrame({
            "gene": res['names'][g],
            "logfoldchange": res['logfoldchanges'][g],
            "pval_adj": res['pvals_adj'][g],
            "group": g
        })
        dfs.append(df.head(n_genes))
    return pd.concat(dfs)

def simple_effect_at_age(counts_df, meta_df, age_value,
                         condition_col="condition", ref_condition="control",
                         n_cpus=8):
    counts_df, meta_df = _align(counts_df, meta_df)
    mask = meta_df["age"].astype(str).str.strip() == str(age_value)
    if mask.sum() < 2:
        raise ValueError(f"Not enough samples at age={age_value}")
    cdf = counts_df.loc[:, mask]
    mdf = meta_df.loc[mask].copy()

    mdf, ref_condition = _prep_condition(mdf, condition_col=condition_col, ref_condition=ref_condition)
    case_label = _detect_case(mdf, condition_col=condition_col, ref_condition=ref_condition)

    ad = _to_anndata(cdf, mdf)
    dds = DeseqDataSet(
        adata=ad,
        design_factors=[condition_col],
        refit_cooks=True,
        inference=DefaultInference(n_cpus=n_cpus),
    )
    dds.deseq2()

    st = DeseqStats(dds, contrast=[condition_col, case_label, ref_condition], inference=DefaultInference(n_cpus=n_cpus))
    st.summary()
    res = st.results_df.copy()
    res = res.rename(columns={
        "log2FoldChange": f"log2FC_{case_label}_vs_{ref_condition}",
        "lfcSE":          f"lfcSE_{case_label}_vs_{ref_condition}",
        "pvalue":         f"pvalue_{case_label}_vs_{ref_condition}",
        "padj":           f"padj_{case_label}_vs_{ref_condition}",
    })
    res.attrs = {"contrast": f"{case_label} vs {ref_condition} @ age {age_value}",
                 "age": str(age_value),
                 "case": case_label,
                 "ref": ref_condition}
    return res

def spatial_neigbourshoods(anndata,
                           cluster_label = 'leiden_0.5',
                           max_distance_allowed = 200):

    import numpy as np
    import scanpy as sc
    import pandas as pd
    from sklearn.metrics.pairwise import euclidean_distances

    distances_input=np.array([anndata.obs['x_centroid'],anndata.obs['y_centroid']])
    din=distances_input.transpose()
    distances=euclidean_distances(din, din)
    dist_df=pd.DataFrame(distances)
    max_distance_allowed=max_distance_allowed
    dist_binary=((dist_df<max_distance_allowed)*1)*((dist_df!=
                                                     0)*1)
    np.sum(np.sum(dist_binary))
    dist_binary['name']=list(anndata.obs[cluster_label])
    distbinsum=dist_binary.groupby('name').sum()
    adata=sc.AnnData(distbinsum.transpose())
    adata.obs=anndata.obs

    return adata

def download_file(
        src_url,
        dst_path,
        force_download=False):
    """
    Download a file.

    Parameters
    ----------
    src_url:
        URL of the file to be downloaded
    dst_path:
        path where the file will be saved
    force_download:
        if True and dst_path exists, overwrite
        (otherwise, assume dst_path is the file you want)
    """
    dst_path = pathlib.Path(dst_path)
    if force_download and dst_path.exists():
        dst_path.unlink()

    if not dst_path.is_file():
        args = [
            "wget",
            src_url,
            "-O",
            str(dst_path),
            "-q"
        ]
        p = subprocess.Popen(args)
        exit_status = p.wait()
        if exit_status != 0:
            if metadata_path.exists():
                metadata_path.unlink()
            raise RuntimeError("Failure downloading metadata file")
        print(f"SUCCESSFULLY DOWNLOADED {dst_path}")
    else:
        print(f"{dst_path} ALREADY EXISTS; NO ACTION TAKEN")

def h5ad_from_subset_of_hdf5(
        hdf5_path,
        row_idx,
        obs,
        gene_list,
        dst_path,
        chunk_size=20000):
    """
    Split off a subset of the cell-by-gene HDF5 file as an h5ad file

    Parameters
    ----------
    hdf5_path:
        path to the cell-by-gene hdf5
    row_idx:
        np.array of integers indicating the rows of the CSV to be
        saved to his h5ad file
    obs:
        obs dataframe for this h5ad file
    gene_list:
        list of gene identifiers for this h5ad file
    dst_path:
        path to the h5ad file being written
    chunk_size:
        number of cells to read/write at a time
    """

    n_cells = len(obs)
    n_genes = len(gene_list)

    last_printed = 0
    r0 = 0
    with h5py.File(hdf5_path, 'r') as src:
        with h5py.File(dst_path, 'w') as dst:
            x_dataset = dst.create_dataset(
                'X',
                shape=(n_cells, n_genes),
                dtype=src['data/counts'].dtype,
                chunks=True,
                compression='gzip',
                compression_opts=4
            )
            x_dataset.attrs.create(
                name='encoding-type',
                data='array'
            )
            x_dataset.attrs.create(
                name='encoding-version',
                data='0.2.0'
            )
            for i0 in range(0, len(row_idx), chunk_size):
                idx_chunk = row_idx[i0:i0+chunk_size]
                chunk = src['data/counts'][:, idx_chunk].transpose()
                x_dataset[r0: r0+chunk.shape[0], :] = chunk
                r0 += chunk.shape[0]
                if i0 >= (last_printed + n_cells//10):
                    last_printed = i0
                    print(f'    wrote {i0+chunk_size:.3e} of {n_cells:.3e} cells')

    # write metadata dataframes to h5ad file
    var = pd.DataFrame(
        [{'gene': gene} for gene in gene_list]
    ).set_index('gene')

    anndata_utils.write_df_to_h5ad(
        h5ad_path=dst_path,
        df_name='var',
        df_value=var)
    anndata_utils.write_df_to_h5ad(
        h5ad_path=dst_path,
        df_name='obs',
        df_value=obs
    )
    print(f'WROTE {dst_path}')
    print(f'contains {n_cells} cells and {n_genes} genes')

def h5ad_from_subset_of_csv(
        csv_path,
        row_idx,
        obs,
        gene_list,
        dst_path):
    """
    Split off a subset of the human cell-by-gene CSV as an h5ad file

    Parameters
    ----------
    csv_path:
        path to the human cell-by-gene CSV
    row_idx:
        np.array of integers indicating the rows of the CSV to be
        saved to his h5ad file
    obs:
        obs dataframe for this h5ad file
    gene_list:
        list of gene identifiers for this h5ad file
    dst_path:
        path to the h5ad file being written
    """

    row_idx_set = set(row_idx)
    n_cells = len(obs)
    n_genes = len(gene_list)

    (descriptor,
     tmp_path) = tempfile.mkstemp(
         suffix='.csv'
     )
    os.close(descriptor)

    tmp_path = pathlib.Path(
        tmp_path
    )

    try:
        # copy chosen rows of CSV to a new CSV at tmp_path
        dtype = dict()
        with open(csv_path, 'r') as src:
            header = src.readline()

            header_params = header.strip().split(',')
            dtype[header_params[0]] = str
            for hh in header_params[1:]:
                dtype[hh] = int

            with open(tmp_path, 'w') as dst:
                dst.write(header)
                for idx, line in enumerate(src):
                    if idx in row_idx_set:
                        dst.write(line)

        # use pd.read_csv with non-zero chunksize to iterate over
        # the tmp CSV one chunk at a time
        chunk_iterator = pd.read_csv(
            tmp_path,
            chunksize=5000,
            dtype=dtype
        )

        with h5py.File(dst_path, 'w') as dst:
            x_dataset = dst.create_dataset(
                'X',
                shape=(n_cells, n_genes),
                dtype=int,
                chunks=True,
                compression='gzip',
                compression_opts=4
            )
            x_dataset.attrs.create(
                name='encoding-type',
                data='array'
            )
            x_dataset.attrs.create(
                name='encoding-version',
                data='0.2.0'
            )
            r0 = 0
            for chunk in chunk_iterator:
                chunk = chunk.set_index('sample_name')
                chunk = chunk.to_numpy()
                x_dataset[r0:r0+chunk.shape[0], :] = chunk
                r0 += chunk.shape[0]

        # write metadata dataframes to h5ad file
        var = pd.DataFrame(
            [{'gene': gene} for gene in gene_list]
        ).set_index('gene')

        anndata_utils.write_df_to_h5ad(
            h5ad_path=dst_path,
            df_name='var',
            df_value=var)
        anndata_utils.write_df_to_h5ad(
            h5ad_path=dst_path,
            df_name='obs',
            df_value=obs
        )

    finally:
        if tmp_path.exists():
            print(f'deleting {tmp_path}')
            tmp_path.unlink()
