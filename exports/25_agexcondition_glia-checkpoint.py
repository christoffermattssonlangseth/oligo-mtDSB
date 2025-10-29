#!/usr/bin/env python
# coding: utf-8

# # mtDSB x Age DE analysis (pydeseq2, LRT-free)

# In[1]:


import scanpy as sc
import decoupler as dc
import pertpy as pt
sc.set_figure_params(figsize=(3, 3), frameon=False)
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
import warnings
warnings.filterwarnings("ignore")

def make_pseudobulk(adata, celltype, ct_col="celltype_merged",
                    sample_cols=("mouse_id","genotype","age"),
                    layer="counts"):
    ad = adata[adata.obs[ct_col] == celltype].copy()
    if layer in ad.layers:
        X = ad.layers[layer]
    else:
        X = ad.X
    if sparse.issparse(X):
        X = X.tocsr()

    # define samples as unique combinations of sample_cols
    meta = ad.obs.loc[:, sample_cols].astype(str)
    meta["sample_id"] = meta.apply(lambda r: "_".join(r.values), axis=1)

    # aggregate counts by sample_id
    sample_ids = meta["sample_id"].values
    uniq = pd.Index(np.unique(sample_ids), name="sample_id")
    # map each cell to column index
    col_index = pd.Series(pd.Categorical(sample_ids, categories=uniq)).cat.codes.values

    # build counts: genes x samples
    n_genes = ad.n_vars
    n_samp  = len(uniq)
    M = np.zeros((n_genes, n_samp), dtype=np.int64)
    # fast aggregate with group-by via numpy
    for s in range(n_samp):
        idx = np.where(col_index == s)[0]
        if len(idx)==0: 
            continue
        Xi = X[idx]
        M[:, s] = Xi.sum(axis=0).A.ravel() if sparse.issparse(Xi) else Xi.sum(axis=0)

    # sample metadata frame
    smd = uniq.to_frame(index=False)
    # split back columns
    parts = smd["sample_id"].str.split("_", expand=True)
    for i, c in enumerate(sample_cols):
        smd[c] = parts[i].values
    smd["condition"] = smd["condition"].astype("category")
    # age: keep both categorical and numeric versions handy
    smd["age"] = smd["age"].astype(str)
    smd["age_num"] = pd.to_numeric(smd["age"], errors="coerce")

    # gene names
    genes = ad.var_names.to_list()

    return pd.DataFrame(M, index=genes, columns=smd["sample_id"]), smd.set_index("sample_id")

# =======================
# mtDSB x Age DE analysis (pydeseq2, LRT-free)
# - Per-age simple effects: model vs control
# - Interaction proxy: delta-LFC between ages with z-test
# Requirements: pydeseq2, anndata, numpy, pandas
# Assumes: pb[ct] = (counts_df [genes x samples], meta_df indexed by sample_id with columns "condition","age")
# =======================

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from anndata import AnnData
from pydeseq2.dds import DeseqDataSet, DefaultInference
from pydeseq2.ds import DeseqStats

# ---------- helpers ----------
def _align(counts_df, meta_df):
    if not counts_df.columns.equals(meta_df.index):
        common = counts_df.columns.intersection(meta_df.index)
        if len(common) == 0:
            raise ValueError("No overlapping sample IDs between counts_df.columns and meta_df.index")
        counts_df = counts_df.loc[:, common]
        meta_df   = meta_df.loc[common]
    return counts_df, meta_df

def _to_anndata(counts_df, meta_df):
    X = counts_df.T.astype(np.int64)  # samples x genes
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
        m = [u for u in uniq if u.lower() == ref_condition.lower()]
        if not m:
            raise ValueError(f"Reference '{ref_condition}' not present. Found: {uniq}")
        ref_condition = m[0]
    mdf[condition_col] = pd.Categorical(
        mdf[condition_col],
        categories=[ref_condition] + [u for u in uniq if u != ref_condition],
        ordered=True
    )
    return mdf, ref_condition

def _detect_case(mdf, condition_col="condition", ref_condition="control"):
    cats = list(mdf[condition_col].cat.categories)
    case_levels = [c for c in cats if c != ref_condition]
    if len(case_levels) == 0:
        raise ValueError("No case level found (only control present).")
    return case_levels[0]

def _compute_gene_stats(counts_df, meta_df):
    """
    Compute gene-level stats to judge expression 'meaningfulness'.
    Returns DataFrame indexed by gene with:
      - mean_cpm_all, detect_rate_all
      - mean_cpm_cond_<cond> (per-condition)
      - mean_cpm_age_<age>   (per-age)
    """
    # CPM
    lib = counts_df.sum(axis=0).replace(0, np.nan)
    cpm = counts_df.divide(lib, axis=1) * 1e6

    gs = pd.DataFrame(index=counts_df.index)
    gs["mean_cpm_all"] = cpm.mean(axis=1)
    gs["detect_rate_all"] = (cpm > 0).mean(axis=1)

    if "condition" in meta_df.columns:
        for cond, cols in meta_df.groupby("condition").groups.items():
            cols = list(cols)
            if cols:
                gs[f"mean_cpm_cond_{cond}"] = cpm[cols].mean(axis=1)

    if "age" in meta_df.columns:
        for age, cols in meta_df.groupby("age").groups.items():
            cols = list(cols)
            if cols:
                gs[f"mean_cpm_age_{age}"] = cpm[cols].mean(axis=1)

    return gs

# ---------- DE pieces ----------
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

# ---------- main runner with gene stats joined ----------
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



# ## Read data

# In[20]:


adata = sc.read_h5ad('../data/mtDNA_DSB_5k_clustered_manual_annotation_with_raw.h5ad')


# In[21]:


pd.DataFrame(adata.var.index).to_csv('markers_mouse_xenium.csv')


# In[22]:


adata.X = adata.layers['counts']
adata.X = adata.X.astype(int)


# In[5]:


glia_classes = [
    # --- Astrocytes ---
    "Telencephalon astrocytes",
    "Olfactory astrocytes",

    # --- Oligodendrocyte lineage ---
    "Mature oligodendrocytes",
    "Oligodendrocytes precursor cells",

    # --- Microglia ---
    "Microglia",
]


# In[6]:


adata = adata[adata.obs.cell_class.isin(glia_classes)]


# In[23]:


adata


# ## Perform pseudobulk

# In[24]:


celltypes = adata.obs.cell_class.unique()
pb = {}   # dict: celltype -> (count_df, sample_metadata)
for ct in celltypes:
    pb[ct] = make_pseudobulk(adata, celltype=ct, ct_col = 'cell_class', sample_cols = ('sample_id', 'condition', 'age'))


# ## Perform DGE

# In[ ]:


results = run_all(pb, condition_col="condition", outdir=None, n_cpus=8, save_csv=False)


# In[26]:


import numpy as np
import matplotlib.pyplot as plt

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


# In[27]:


results.keys()


# In[28]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings
import math

def plot_expression_by_condition_grid(
    pb,
    ct,
    genes,
    norm=True,                 # CPM-normalize per sample
    q=None,                    # optional: highlight sig ages (not used here; keep for parity)
    ncols=4,                   # number of subplot columns
    sharey=True,               # share y-axis across genes for magnitude comparison
    figsize_per_panel=(4.0, 3.2),
    fontsize_base=13,
    suppress_warnings=True
):
    """
    Multi-panel plot of mean expression (pseudobulk) across ages for each condition.

    pb: dict-like, pb[ct] = (counts_df [genes x samples], meta_df [samples x meta])
    ct: cell class key present in pb
    genes: list of gene symbols to plot
    norm: if True, CPM-normalize per-sample library size before group averaging
    sharey: if True, use shared y-axis across all gene subplots
    """

    if suppress_warnings:
        warnings.filterwarnings("ignore")

    counts_df, meta_df = pb[ct]
    # keep only genes that exist
    genes_present = [g for g in genes if g in counts_df.index]
    missing = sorted(set(genes) - set(genes_present))
    if missing:
        print(f"Skipping missing genes ({len(missing)}): {', '.join(missing[:8])}" + ("..." if len(missing)>8 else ""))

    if len(genes_present) == 0:
        print("No requested genes found in counts_df.")
        return

    # Build long table for requested genes
    # sample-wise counts and metadata
    dat = counts_df.loc[genes_present].T  # samples x genes
    df = dat.stack().reset_index()
    df.columns = ["sample", "gene", "counts"]

    # attach meta
    meta = meta_df.copy()
    meta = meta.loc[df["sample"].unique()]  # ensure order/coverage
    # basic hygiene
    meta = meta.assign(
        age=meta["age"].astype(str).str.strip(),
        condition=meta["condition"].astype(str).str.strip()
    )
    df = df.merge(meta[["age", "condition"]], left_on="sample", right_index=True, how="left")

    # CPM normalize if requested
    if norm:
        lib_sizes = counts_df.sum(axis=0)  # per-sample library size
        df["cpm"] = df.apply(lambda r: (r["counts"] / lib_sizes.loc[r["sample"]]) * 1e6, axis=1)
        val_col = "cpm"
    else:
        val_col = "counts"

    # aggregate mean ± sem per condition × age × gene
    agg = (
        df.groupby(["gene", "condition", "age"], observed=True)[val_col]
          .agg(mean="mean", sem=lambda x: x.std(ddof=1)/np.sqrt(len(x)) if len(x)>1 else 0.0)
          .reset_index()
    )

    # sort ages numerically if possible
    try:
        age_order = sorted(agg["age"].unique(), key=lambda x: float(x))
    except Exception:
        age_order = sorted(agg["age"].unique())

    # figure layout
    n = len(genes_present)
    nrows = math.ceil(n / ncols)
    fig_w = max(ncols * figsize_per_panel[0], 6.0)
    fig_h = max(nrows * figsize_per_panel[1], 3.5)

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(fig_w, fig_h),
        sharex=True,
        sharey=sharey,
        constrained_layout=True
    )
    axes = np.ravel(axes)

    # If sharey=True, compute global y-limits from all panels for consistency
    ymins, ymaxs = [], []
    if sharey:
        for g in genes_present:
            sub = agg[agg["gene"] == g]
            if sub.empty:
                continue
            ymins.append((sub["mean"] - sub["sem"]).min())
            ymaxs.append((sub["mean"] + sub["sem"]).max())
        if ymins and ymaxs:
            y_min = 0 if norm else min(0, np.nanmin(ymins))
            y_max = np.nanmax(ymaxs) * 1.08
        else:
            y_min, y_max = None, None

    # per-plot drawing
    for i, g in enumerate(genes_present):
        ax = axes[i]
        sub = agg[agg["gene"] == g].copy()
        if sub.empty:
            ax.set_visible(False)
            continue

        # draw each condition line with error bars
        for cond, grp in sub.groupby("condition", observed=True):
            grp = grp.set_index("age").reindex(age_order).reset_index()
            ax.errorbar(
                grp["age"], grp["mean"], yerr=grp["sem"],
                marker="o", capsize=4, lw=2.0, markersize=6, label=str(cond)
            )

        ax.axhline(0, ls="--", color="0.7", lw=1)  # helpful baseline (works for CPM too)
        ax.set_title(g, fontsize=fontsize_base, pad=6)
        ax.tick_params(axis="x", labelrotation=45, labelsize=fontsize_base-2)
        ax.tick_params(axis="y", labelsize=fontsize_base-2)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        # set shared y limits if requested
        if sharey and (y_min is not None) and (y_max is not None):
            ax.set_ylim(y_min, y_max)

        # lightweight legend per first subplot only
        if i == 0:
            ax.legend(frameon=False, fontsize=fontsize_base-2, title="Condition", title_fontsize=fontsize_base-2)
        else:
            ax.legend().remove()

    # hide unused panels
    for j in range(i+1, len(axes)):
        axes[j].set_visible(False)

    # global labels offset away from plots
    units = "CPM" if norm else "counts"
    fig.suptitle(f"{ct}: mean {units} across ages by condition",
                 fontsize=fontsize_base+3, y=1.04, fontweight="bold")
    fig.text(0.5, -0.035, "Age (weeks)", ha="center", fontsize=fontsize_base+1)
    fig.text(-0.035, 0.5, f"Mean {units}", va="center", rotation="vertical", fontsize=fontsize_base+1)

    # extra margins so labels never overlap
    plt.subplots_adjust(bottom=0.14, left=0.11, right=0.99, top=0.93, wspace=0.25, hspace=0.35)
    plt.show()


# In[29]:


import numpy as np
import pandas as pd

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

# --- helper: filter genes by mean expression (CPM) but always keep 'forced'
def apply_expr_filter(res21, res60, expr_series, min_cpm=0.5, forced=None):
    """
    Keep genes with mean expression >= min_cpm (based on expr_series),
    intersected with genes present in both res tables.
    Always keep 'forced' genes if provided.
    """
    if expr_series is None:
        # nothing to filter on
        keep = res21.index.intersection(res60.index)
        return res21.loc[keep], res60.loc[keep], None

    # ensure alignment
    common = res21.index.intersection(res60.index).intersection(expr_series.index)
    if len(common) == 0:
        raise ValueError("No overlapping genes among res21, res60, and expr_series.")

    # main gate
    keep = set(expr_series.loc[common][expr_series.loc[common] >= min_cpm].index)

    # always-keep
    if forced:
        keep |= set([g for g in forced if g in common])

    keep = sorted(keep)
    # pass back a filtered expression series too (nice for point sizes)
    return res21.loc[keep], res60.loc[keep], expr_series.loc[keep]


# In[30]:


from IPython.display import display, HTML
import pandas as pd
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from adjustText import adjust_text

# --- Keep vector text editable for Illustrator ---
mpl.rcParams['svg.fonttype'] = 'none'   # text stays as text
mpl.rcParams['pdf.fonttype'] = 42       # embed fonts in PDF
mpl.rcParams['text.usetex'] = False

# --- Output directory ---
outdir = "../results/figures"
os.makedirs(outdir, exist_ok=True)


def _display_top_table(res, top, cell, age, head=20):
    """Display a nicely formatted subset of DE results."""
    if res is None or not isinstance(res, pd.DataFrame) or len(top) == 0:
        display(HTML(f"<h4 style='color:gray;margin:0.6em 0;'>⚠️ {cell} — {age}: no data</h4>"))
        return
    idx = res.index.intersection(pd.Index(top))
    if len(idx) == 0:
        display(HTML(f"<h4 style='color:gray;margin:0.6em 0;'>⚠️ {cell} — {age}: none of the requested genes found</h4>"))
        return

    sub = res.loc[idx].copy()
    num_cols = sub.select_dtypes(include="number").columns
    sub[num_cols] = sub[num_cols].round(3)

    display(HTML(f"<h3 style='color:#2a4d69;margin:0.6em 0 0.2em 0;'>{cell} — Top {len(idx)} genes at {age}</h3>"))
    display(sub.head(head))


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


# ## Compare LFC across age and cell type

# In[31]:


import matplotlib.pyplot as plt
import numpy as np

cutoff = 60
cells = list(results.keys())
n = len(cells)

# Choose grid layout
ncols = 2
nrows = int(np.ceil(n / ncols))

fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows), constrained_layout=True)

for ax, cell in zip(axes.flat, cells):
    df_hist = results[cell]['simple_effects'].get('60')
    if df_hist is None or 'baseMean' not in df_hist.columns:
        ax.set_visible(False)
        continue

    # Plot histogram
    ax.hist(np.log10(df_hist['baseMean'] + 1), bins=80, color='slateblue', edgecolor='black', alpha=0.7)
    ax.axvline(np.log10(cutoff), color='red', linestyle='--', lw=1.2, label=f'cutoff={cutoff}')
    ax.set_title(cell, fontsize=10)
    ax.set_xlabel('log10(baseMean + 1)')
    ax.set_ylabel('Number of genes')
    ax.legend(fontsize=8)

# Hide any unused axes
for i in range(len(cells), nrows * ncols):
    axes.flat[i].set_visible(False)

plt.suptitle("Distribution of average expression (baseMean) across cell types", fontsize=14)
plt.show()


# In[33]:


# === Main loop ===
forced = ['Trib3', 'Gdf15','Fgf21','Atf4']
MIN_CPM = 80  # <-- tweak threshold here

for cell, pack in results.items():
    se = pack.get("simple_effects", {})
    if "21" not in se or "60" not in se or se["21"] is None or se["60"] is None:
        print(f"Skipping {cell}: missing 21/60 results.")
        continue

    res21 = se["21"]
    res60 = se["60"]

    gstats = pack.get("gene_stats", None)
    expr_series = gstats["mean_cpm_all"] if (gstats is not None and "mean_cpm_all" in gstats.columns) else None
    if expr_series is not None and not isinstance(expr_series, pd.Series):
        # if it's a column from a DF, make sure it's a Series with gene index
        expr_series = expr_series.squeeze()

    # --- NEW: expression filter (keeps 'forced' regardless of CPM)
    res21_f, res60_f, expr_f = apply_expr_filter(res21, res60, expr_series, min_cpm=MIN_CPM, forced=forced)
    if cell=='Mature oligodendrocytes':
        top_n_each = 30
    else: 
        top_n_each = 20
    # pick tops on the FILTERED sets
    picks, top21, top60 = pick_top_per_age(
        res21_f, res60_f,
        lfc_col='log2FC_mtDSB_vs_control',
        padj_col='pvalue_mtDSB_vs_control',
        q=0.25, top_n_each=top_n_each,
        alpha_lfc=1.0, beta_p=0.0,
        force_lfc_min=0.0, max_labels=100
    )

    # labels = selected + forced (dedup) but only for genes that survived filter
    labels = [g for g in dict.fromkeys(list(picks) + forced) if g in res21_f.index]

    safe_name = cell.replace(" ", "_").replace("/", "_")
    save_base = os.path.join(outdir, f"{safe_name}_21w_vs_60w")

    _ = compare_lfc_scatter_colored_labels(
        res21_f, res60_f,
        highlight_genes=labels,
        expr=expr_f,  # sized points reflect filtered mean CPM
        title=f"{cell}: 21 w vs 60 w (LFC prioritized; size ~ expression, min_cpm={MIN_CPM})",
        xlim=(-5, 5), ylim=(-5, 5),
        save_base=save_base
    )

    _display_top_table(res21_f, top21, cell, "21 w")
    _display_top_table(res60_f, top60, cell, "60 w")


# In[ ]:





# ## mitokine expression

# In[54]:


mitokine_like = [
    "Gdf15", "Adm", "Cst7", "Igfbp3", "Serpina3n",
    "Cdkn1a", "Maff", "Eif4ebp1", "Aldh1l2"
]


# In[55]:


import pandas as pd

records = []
for cell, pack in results.items():
    se = pack.get("simple_effects", {})
    for age in ("21", "60"):
        if age in se and isinstance(se[age], pd.DataFrame):
            res = se[age]
            for gene in mitokine_like:
                if gene in res.index:
                    r = res.loc[gene]
                    records.append({
                        "cell_class": cell,
                        "age": age,
                        "gene": gene,
                        "log2FC": r.get("log2FC_mtDSB_vs_control", None),
                        "padj": r.get("pvalue_mtDSB_vs_control", None),
                        "mean_cpm": r.get("mean_cpm_all", None)
                    })
df = pd.DataFrame(records)


# In[57]:


import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# --- Abbreviations for cell class labels ---
cell_abbrev = {
    "Mature oligodendrocytes": "OL",
    "Oligodendrocytes precursor cells": "OPCs",
    "Olfactory astrocytes": "Ast (olf)",
    "Telencephalon astrocytes": "Ast (tel)",
    "Microglia": "MG",
}

# --- Prepare data ---
bars = (
    df[df["gene"].isin(mitokine_like)]
      .pivot_table(index=["gene", "cell_class", "age"], values="log2FC")
      .reset_index()
)

# Map abbreviations
bars["cell_class"] = bars["cell_class"].map(cell_abbrev).fillna(bars["cell_class"])

# Define consistent order
x_order_full = ["OL", "OPCs", "Ast (olf)", "Ast (tel)", "MG"]
x_order = [x for x in x_order_full if x in bars["cell_class"].unique()]

# --- Plot ---
sns.set_style("whitegrid")
g = sns.catplot(
    data=bars, kind="bar",
    x="cell_class", y="log2FC", hue="age",
    col="gene", col_wrap=4,
    order=x_order,
    palette=["#5DA5DA", "#F15854"],
    height=3.2,
    sharex=False, sharey=True
)

# --- Force draw so ticks exist ---
plt.draw()

# --- Customize each facet ---
for ax in g.axes.flat:
    ax.axhline(0, ls=":", c="0.6", lw=0.8)
    ax.set_xlabel("")
    ax.set_ylabel("log₂FC", fontsize=10)
    ax.set_xticks(np.arange(len(x_order)))
    ax.set_xticklabels(x_order, rotation=90, fontsize=15, fontweight="medium")
    ax.tick_params(axis='x', which='both', bottom=True, top=False, length=3, width=0.8)

    # Italicize gene titles
    title = ax.get_title()
    gene = title.split(" = ", 1)[-1] if " = " in title else title
    ax.set_title(rf"$\it{{{gene}}}$", fontsize=11)

# --- Move legend fully outside ---
g._legend.set_bbox_to_anchor((1.04, 0.5))
g._legend.set_title("Age")
g._legend._legend_box.align = "left"

# --- Figure title & layout ---
g.fig.suptitle(
    "Mitokine-like genes: age-dependent response by glial population",
    y=1.03,
    fontsize=13
)
plt.tight_layout(rect=[0, 0, 0.93, 1])  # extra room for legend

# --- Optional: SVG export (editable text in Illustrator) ---
g.fig.savefig("../results/figures/mitokine_bars_facets.png", bbox_inches="tight", format="png")

plt.show()


# ## Save results

# In[58]:


import os

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


# In[59]:


save_all_results(results, outdir="../results/age")


# In[60]:


import scanpy as sc


# In[61]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import scipy.sparse as sp

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


# In[466]:


sc.pp.normalize_total(adata, target_sum=1000)
sc.pp.log1p(adata)


# In[471]:





# In[736]:


mitokine_like = [
    "Gdf15",
]


# In[737]:


for gene in mitokine_like:
    print(gene)
    fig, axes = plot_spatial_grid_by_condition_age(
        adata,
        gene=gene,
        layer=None,        # or "log1p"/"counts"
        use_raw=False,
        cmap="coolwarm",
        percentile=99,
        spot_size=30,
        condition_col="condition",
        age_col="age",
        sample_col="sample_id",
    )


# In[501]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

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


# In[505]:


to_keep = ['Telencephalon astrocytes',
 'Olfactory astrocytes',
 'Mature oligodendrocytes',
 'Oligodendrocytes precursor cells']


# In[ ]:


for cell in to_keep:
    print(cell)
    fig, axes = plot_spatial_grid_basic(
        adata,
        condition_col="condition",
        age_col="age",
        sample_col="sample_id",
        highlight_col="cell_class",
        highlight_value=cell,
        highlight_color="tab:red",
        spot_size=5,
        background_alpha=0.1,

    )


# In[630]:


adata.X.max()


# In[ ]:





# In[ ]:




