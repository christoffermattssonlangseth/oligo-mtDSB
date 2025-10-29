"""
Auto-generated utilities for misc_utils.
Do not edit by hand without moving changes back into notebooks.

Each function below was extracted from exported analysis notebooks.
"""

"""
Auto-generated utilities for misc_utils.
Do not edit by hand without moving changes back into notebooks.

Each function below was extracted from exported analysis notebooks.
"""

from typing import *
import os
import re
import glob
from collections import defaultdict

import numpy as np
import pandas as pd

import scanpy as sc
import anndata as ad
from anndata import AnnData

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec
from matplotlib.colors import Normalize
import matplotlib.colors as mcolors

import seaborn as sns
from adjustText import adjust_text

from scipy import sparse
import scipy.sparse as sp
from scipy.spatial import Delaunay
from scipy.stats import norm, combine_pvalues
from statsmodels.stats.multitest import multipletests

from sklearn.neighbors import radius_neighbors_graph

from skimage.draw import polygon2mask
from skimage.morphology import binary_closing
from skimage.measure import find_contours
from shapely.geometry import Polygon  # for generate_smoothed_alpha_polygon

import pymc as pm
import arviz as az

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.inference import DefaultInference

def _as_csr(X):
    if sp.issparse(X):
        return X.tocsr()
    # dense -> csr
    return sp.csr_matrix(np.asarray(X), copy=False)

def pseudobulk_by_groups(adata, group_cols, layer=None, keep_celltypes=None):
    """
    Sum raw counts per (celltype, sample, age, condition).
    Returns:
      pb_counts : dense float32 array (n_groups x n_genes)
      groups_df : dataframe with group labels (n_groups x len(group_cols))
      var_names : pandas Index of gene names
    """
    # 0) pull matrix & obs
    if layer is None:
        X = _as_csr(adata.X)
    else:
        if layer not in adata.layers:
            raise ValueError(f"Layer '{layer}' not found in adata.layers")
        X = _as_csr(adata.layers[layer])

    obs = adata.obs.copy()

    # 1) optional filter of cell types
    if keep_celltypes is not None:
        obs = obs.loc[obs[CELLTYPE_COL].isin(keep_celltypes)]
        X   = X[obs.index, :]   # keep matrix aligned

    # 2) make sure group columns exist and are strings (avoid categorical gotchas)
    for c in group_cols:
        if c not in obs.columns:
            raise ValueError(f"Missing obs column: {c}")
    gdf = obs.loc[:, group_cols].astype(str).copy()

    # 3) build an integer code per group
    key = pd.MultiIndex.from_frame(gdf, names=group_cols)
    codes, uniques = pd.factorize(key, sort=True)

    # 4) aggregate counts by group using a sparse trick
    n_groups = len(uniques)
    n_genes  = X.shape[1]
    ones     = np.ones(X.shape[0], dtype=np.int64)
    row_idx  = codes  # which group each cell belongs to

    # make a sparse "group-by" matrix G (cells -> groups), then G @ X = group sums
    G = sp.csr_matrix((ones, (row_idx, np.arange(X.shape[0]))), shape=(n_groups, X.shape[0]))
    pb_sparse = G @ X
    pb_counts = np.asarray(pb_sparse.todense(), dtype=np.float32)

    # 5) unpack group labels
    groups_df = uniques.to_frame(index=False)  # columns=group_cols
    groups_df.columns = group_cols

    return pb_counts, groups_df.reset_index(drop=True), adata.var_names.copy()

def normalize_pseudobulk(pb_counts, method="cpm_log1p", axis=1):
    """
    Normalize GxP counts (G=groups, P=genes).
    method="cpm_log1p": counts per million (per group) then log1p
    method="log1p_sizefactor": sizefactor = library_size / median(library_size), log1p(counts/sizefactor)
    """
    pb = pb_counts.astype(np.float64, copy=True)
    lib = pb.sum(axis=axis, keepdims=True)  # group library sizes
    if method == "cpm_log1p":
        cpm = (pb / np.clip(lib, 1.0, None)) * 1e6
        return np.log1p(cpm)
    elif method == "log1p_sizefactor":
        size_factor = lib / np.median(lib)
        return np.log1p(pb / np.clip(size_factor, 1e-12, None))
    else:
        raise ValueError(f"Unknown method: {method}")

def pseudobulk_to_long_with_celltype(pb_norm, groups_df, var_names,
                                     sample_col=SAMPLE_COL, age_col=AGE_COL,
                                     cond_col=COND_COL, celltype_col=CELLTYPE_COL,
                                     value_name="value"):
    """
    Returns tidy df: [sample, age, condition, cell_type, gene, value]
    """
    n_groups, n_genes = pb_norm.shape
    assert n_groups == len(groups_df), "shape mismatch"

    # build a dataframe carefully
    df = pd.DataFrame(pb_norm, columns=var_names)
    df.insert(0, sample_col, groups_df[sample_col].values)
    df.insert(1, age_col,    groups_df[age_col].values)
    df.insert(2, cond_col,   groups_df[cond_col].values)
    df.insert(3, celltype_col, groups_df[celltype_col].values)

    # wide -> long
    long = df.melt(
        id_vars=[sample_col, age_col, cond_col, celltype_col],
        var_name="gene", value_name=value_name
    )
    # ensure expected dtypes
    long[age_col] = long[age_col].astype(str)
    long[cond_col] = long[cond_col].astype(str)
    long[celltype_col] = long[celltype_col].astype(str)
    long[sample_col] = long[sample_col].astype(str)

    return long

def fit_gene_pymc_ac(df_gene, *, advi=False, draws=1000, tune=1000, chains=4,
                     target_accept=0.9, seed=42):
    """
    value ~ age + condition + age:condition + (1|sample)
    df_gene is one (cell_type, gene).
    """
    df = df_gene.copy()
    ages    = sorted(df["age"].astype(str).unique())
    conds   = sorted(df["condition"].astype(str).unique())
    samples = sorted(df["sample_id"].astype(str).unique())

    A = len(ages); C = len(conds); S = len(samples)

    age_idx    = df["age"].astype(str).apply(ages.index).to_numpy()
    cond_idx   = df["condition"].astype(str).apply(conds.index).to_numpy()
    sample_idx = df["sample_id"].astype(str).apply(samples.index).to_numpy()
    y          = df["value"].to_numpy()

    with pm.Model() as m:
        sigma         = pm.HalfNormal("sigma", 0.7)
        sigma_sample  = pm.HalfNormal("sigma_sample", 0.5)
        sample_offset = pm.Normal("sample_offset", 0.0, sigma_sample, shape=S)

        beta_age         = pm.Normal("beta_age", 0.0, 0.5, shape=A)
        beta_cond        = pm.Normal("beta_cond", 0.0, 0.5, shape=C)
        beta_interaction = pm.Normal("beta_interaction", 0.0, 0.5, shape=(A, C))

        mu = (
            beta_age[age_idx]
          + beta_cond[cond_idx]
          + beta_interaction[age_idx, cond_idx]
          + sample_offset[sample_idx]
        )

        pm.Normal("obs", mu, sigma, observed=y)

        if advi:
            approx = pm.fit(20_000, random_seed=seed, progressbar=False)
            idata  = approx.sample(draws)
        else:
            idata = pm.sample(draws=draws, tune=tune, chains=chains,
                              target_accept=target_accept, random_seed=seed,
                              return_inferencedata=True, progressbar=True)

    meta = {"ages": ages, "conds": conds}
    return idata, meta

def _summarize_1d(samples):
    q = np.quantile(samples, [0.025, 0.975])
    return float(samples.mean()), float(q[0]), float(q[1])

def summarize_effects_ac(idata, meta):
    """
    Report:
      - mtDSB_at_age  = (mtDSB - control) at each age
      - age_interaction(old-young) = difference of those contrasts between oldest and youngest age
    """
    ages  = meta["ages"]
    conds = [c.lower() for c in meta["conds"]]

    # indices for control & mtDSB
    i_ctrl = conds.index("control") if "control" in conds else 0
    i_mt   = conds.index("mtdsb")   if "mtdsb"   in conds else (1 if len(conds)>1 else 0)

    beta_cond = idata.posterior["beta_cond"].values.reshape(-1, len(conds))
    beta_int  = idata.posterior["beta_interaction"].values.reshape(-1, len(ages), len(conds))

    rows = []

    # mtDSB effect at each age
    for a, age in enumerate(ages):
        diff = (beta_cond[:, i_mt] - beta_cond[:, i_ctrl]) + (beta_int[:, a, i_mt] - beta_int[:, a, i_ctrl])
        mean, lo, hi = _summarize_1d(diff)
        rows.append({
            "effect": "mtDSB_at_age", "age": age,
            "mean": mean, "hdi_2.5%": lo, "hdi_97.5%": hi,
            "term": f"beta[{age}] (mtDSB effect @ age)"
        })

    # age interaction: oldest – youngest
    if len(ages) >= 2:
        a_lo, a_hi = 0, len(ages)-1
        d_hi = (beta_cond[:, i_mt] - beta_cond[:, i_ctrl]) + (beta_int[:, a_hi, i_mt] - beta_int[:, a_hi, i_ctrl])
        d_lo = (beta_cond[:, i_mt] - beta_cond[:, i_ctrl]) + (beta_int[:, a_lo, i_mt] - beta_int[:, a_lo, i_ctrl])
        ai  = d_hi - d_lo
        mean, lo, hi = _summarize_1d(ai)
        rows.append({
            "effect": "age_interaction", "age": f"{ages[a_hi]}–{ages[a_lo]}",
            "mean": mean, "hdi_2.5%": lo, "hdi_97.5%": hi,
            "term": f"age_interaction({ages[a_hi]} - {ages[a_lo]})"
        })
    return pd.DataFrame(rows)

def fit_gene_pymc_condition(
    df_gene,
    aggregate_by_sample=True,
    sample_col="sample",          # or "sample_id"
    advi=False,
    draws=800, tune=800, chains=2, cores=2,
    target_accept=0.9, seed=42,
    progressbar=False
):
    dg = df_gene.copy()

    # be forgiving about sample column name
    if aggregate_by_sample and sample_col not in dg.columns:
        if "sample_id" in dg.columns:
            sample_col = "sample_id"
        else:
            aggregate_by_sample = False  # no sample column → skip aggregation

    if aggregate_by_sample:
        dg = (dg.groupby([sample_col, "condition"], as_index=False)
                .agg(value=("value","mean")))

    cond = _ensure_condition_cats(dg["condition"])
    cond_idx = pd.Categorical(cond).codes
    C = len(pd.Categorical(cond).categories)
    y = dg["value"].values

    with pm.Model() as m:
        alpha = pm.Normal("alpha", 0, 1)
        sigma = pm.HalfNormal("sigma", 0.5)
        beta  = pm.Normal("beta", 0, 1, shape=C)   # condition main effects
        mu = alpha + beta[cond_idx]
        pm.Normal("obs", mu, sigma, observed=y)

        if advi:
            approx = pm.fit(
                20_000,
                random_seed=seed,
                callbacks=[pm.callbacks.CheckParametersConvergence()],
                progressbar=progressbar,   # OK here
            )
            idata = approx.sample(draws, random_seed=seed)  # ← NO progressbar kw here
        else:
            idata = pm.sample(
                draws=draws, tune=tune, chains=chains, cores=cores,
                target_accept=target_accept, random_seed=seed,
                return_inferencedata=True, progressbar=progressbar
            )

    return idata, {"categories": list(pd.Categorical(cond).categories)}

def make_shortlist(screen_df, top_k=300, by="abs_effect"):
    if screen_df is None or screen_df.empty:
        return []
    df = screen_df.copy()
    df["ols_fdr"] = multipletests(df["ols_p"].values, method="fdr_bh")[1]
    if by == "abs_effect":
        df = df.reindex(df["ols_effect"].abs().sort_values(ascending=False).index)
    elif by == "fdr":
        df = df.sort_values("ols_fdr")
    else:
        raise ValueError("by must be 'abs_effect' or 'fdr'")
    return df.drop_duplicates("gene").head(top_k)["gene"].tolist()

def summarize_condition_effect(idata, meta):
    post = az.summary(idata, var_names=["beta"], hdi_prob=0.95).reset_index()
    post = post.rename(columns={"index": "param"})
    post["condition"] = meta["categories"]
    # Keep the mtDSB effect row only (if your categories are ['control','mtDSB'])
    return post

def _ensure_condition_cats(series):
    if "control" in series.unique() and "mtDSB" in series.unique():
        return pd.Categorical(series, categories=["control","mtDSB"])
    return pd.Categorical(series)

def _fit_one_gene_cond(df_long, gene, cell_type, **kwargs):
    sub = df_long[(df_long["cell_class"]==cell_type) & (df_long["gene"]==gene)]
    if sub.empty or sub["condition"].nunique() < 2:
        return None
    # jitter seed per gene so workers don’t collide
    kwargs = {**kwargs, "seed": (kwargs.get("seed", 42) + (hash(gene) % 10_000))}
    try:
        idata, meta = fit_gene_pymc_condition(sub, **kwargs)
        eff = summarize_condition_effect(idata, meta)
        if eff.empty:
            return None
        eff.insert(0, "gene", gene)
        eff.insert(0, "cell_class", cell_type)
        return eff
    except Exception as e:
        # optional: print(gene, "failed:", e)
        return None

def _p(x,p): return np.nanpercentile(x, p)

def _to_dense(X):
    return X.toarray() if sp.issparse(X) else np.asarray(X)

def _align(counts_df, meta_df):
    if not counts_df.columns.equals(meta_df.index):
        common = counts_df.columns.intersection(meta_df.index)
        if len(common) == 0:
            raise ValueError("No overlapping sample IDs between counts_df.columns and meta_df.index")
        counts_df = counts_df.loc[:, common]
        meta_df   = meta_df.loc[common]
    return counts_df, meta_df

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

def _get_age_df(results_by_region, *, cell_type, region, age, table="simple_effects"):
    df = results_by_region[cell_type][region][table][age]
    return df.set_index("gene") if "gene" in df.columns else df

def _apply_expr_filter(res_a, res_b, expr_series, min_cpm=80, forced=None):
    forced = set(forced or [])
    genes = res_a.index.intersection(res_b.index)
    A, B = res_a.loc[genes].copy(), res_b.loc[genes].copy()
    if expr_series is not None:
        e = expr_series.reindex(genes).fillna(0)
        keep = (e >= float(min_cpm)) | pd.Index(genes).isin(forced)
        return A.loc[keep], B.loc[keep], e.loc[keep]
    return A, B, pd.Series(np.nan, index=genes)

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

def _ensure_gene_index(df):
    """Make gene names the index if there's a 'gene' column."""
    if df is None or not isinstance(df, pd.DataFrame) or df.empty:
        return df
    if df.index.name is None and "gene" in df.columns:
        df = df.set_index("gene")
    return df

def save_all_region_results(results_by_region, outdir="../../results_regionwise"):
    """
    Save all DE results from `results_by_region` into organized subfolders.
    Structure:
        outdir/
          └── <cell_class>/
                ├── <region>/
                │     ├── simple_effect_age_<age>.csv
                │     └── ...
                └── summary_region_impact.csv
    """
    os.makedirs(outdir, exist_ok=True)

    for cell_class, region_dict in results_by_region.items():
        cell_dir = os.path.join(outdir, cell_class.replace(" ", "_"))
        os.makedirs(cell_dir, exist_ok=True)
        print(f"\n📂 Saving {cell_class} ...")

        region_summary_rows = []

        for region, reg_data in region_dict.items():
            reg_dir = os.path.join(cell_dir, region.replace(" ", "_"))
            os.makedirs(reg_dir, exist_ok=True)

            # --- Simple effects per age (mtDSB vs control) ---
            simple_effects = reg_data.get("simple_effects", {})
            for age, df in simple_effects.items():
                if df is not None and not df.empty:
                    fname = f"simple_effect_age_{age}.csv"
                    df.to_csv(os.path.join(reg_dir, fname))

            # --- Region-level impact summary ---
            imp = reg_data.get("impact", {})
            for age, metrics in imp.items():
                if metrics:
                    region_summary_rows.append({
                        "region": region,
                        "age": age,
                        "n_sig": metrics.get("n_sig", np.nan),
                        "sum_neglog10q": metrics.get("sum_neglog10q", np.nan),
                    })

        # --- Save region summary table per cell_class ---
        if region_summary_rows:
            summary_df = pd.DataFrame(region_summary_rows)
            summary_df.to_csv(os.path.join(cell_dir, "summary_region_impact.csv"), index=False)

        print(f"✅ Saved regionwise results for {cell_class} to {cell_dir}")

def label_by_max_score(ad, group_key="micro_leiden_1.0", panels=panels):
    scores = {}
    for name, genes in panels.items():
        gs = [g for g in genes if g in ad.var_names]
        if not gs: 
            continue
        sc.tl.score_genes(ad, gs, score_name=f"score_{name}", use_raw=False)
        scores[name] = f"score_{name}"
    labels = {}
    for cl in ad.obs[group_key].cat.categories:
        sub = ad[ad.obs[group_key] == cl]
        means = {name: float(sub.obs[scol].mean()) for name, scol in scores.items()}
        labels[cl] = max(means, key=means.get) if means else "Unknown"
    return labels

def natural_key(s):
    m = re.search(r'(\d+)', os.path.basename(s))
    return int(m.group(1)) if m else s

def shortlist_present(df_filtered, shortlists, celltype_col="cell_class"):
    present = {}
    for ct, genes in shortlists.items():
        have = set(df_filtered.loc[df_filtered[celltype_col] == ct, "gene"].unique())
        keep = [g for g in genes if g in have]
        if not keep:
            print(f"[warn] No shortlisted genes found in df for '{ct}' (skipping).")
        else:
            print(f"{ct}: using {len(keep)}/{len(genes)} shortlisted genes present in data")
        present[ct] = keep
    return present

def label_effect(row):
    e = row["effect"]
    if e.startswith("mtDSB_at_age"):
        return f"mtDSB @ age {row['age']}"
    if "age_interaction" in e or "interaction" in e:
        return "age interaction (60–21)"
    return e

def bootstrap_ci(data1, data2, n=1000):
    rng = np.random.default_rng(42)
    diffs = []
    for _ in range(n):
        s1 = rng.choice(data1, size=len(data1), replace=True)
        s2 = rng.choice(data2, size=len(data2), replace=True)
        diffs.append(s2.mean() - s1.mean())
    return np.percentile(diffs, [2.5, 97.5])

def _detect_cols(df):
    """Map whatever you have to standard keys we need."""
    cols = {c.lower(): c for c in df.columns}
    def pick(candidates, required=True):
        for c in candidates:
            if c in cols: 
                return cols[c]
        if required:
            raise KeyError(f"Could not find any of columns: {candidates}")
        return None

    return {
        "gene":       pick(["gene"]),
        "cell_class": pick(["cell_class","celltype","cell_type"]),
        "condition":  pick(["condition","cond","group"]),
        "value":      pick(["value","expr","expression","y"]),
        "sample":     pick(["sample","sample_id","donor","subject"])
    }

def _coerce_condition(d, condition_col, positive="mtDSB", reference="control"):
    # ensure only two levels, ordered [reference, positive]
    d = d.copy()
    d[condition_col] = d[condition_col].astype(str).str.strip()
    levels = sorted(d[condition_col].unique().tolist())
    if reference not in levels or positive not in levels:
        raise ValueError(
            f"Condition levels must include '{reference}' and '{positive}'. "
            f"Found: {levels}"
        )
    d[condition_col] = pd.Categorical(d[condition_col], categories=[reference, positive], ordered=True)
    return d

def score_signatures(
    adata,
    signatures: dict,
    prefix: str = "score_",
    ctrl_size: int = 50,
    n_bins: int = 25,
    use_raw: bool = False,
    verbose: bool = True,
):
    """
    Adds one column per signature to adata.obs: f'{prefix}{name}'.
    Auto-filters genes to those present in adata.var_names.
    """
    var_set = set(adata.var_names.astype(str))
    added = []
    for name, genes in signatures.items():
        genes = [g for g in genes if g in var_set]
        if len(genes) == 0:
            if verbose:
                print(f"⚠️  Skipping '{name}': no genes found in var_names.")
            continue
        if verbose:
            print(f"→ Scoring {name:12s} with {len(genes):2d} genes.")
        sc.tl.score_genes(
            adata,
            gene_list=genes,
            score_name=f"{prefix}{name}",
            ctrl_size=ctrl_size,
            n_bins=n_bins,
            use_raw=use_raw
        )
        added.append(f"{prefix}{name}")
    return added

def z(x):
    x = np.asarray(x, float)
    mu, sd = np.nanmean(x), np.nanstd(x)
    return (x - mu) / (sd if sd > 0 else 1)

def _get_counts_matrix(adata, layer=None):
    """Return (matrix, var_names) where matrix is CSR counts (cells x genes)."""
    if layer is not None and layer in adata.layers:
        M = adata.layers[layer]
    else:
        M = adata.X
    if sp.issparse(M):
        M = M.tocsr()
    else:
        M = sp.csr_matrix(M)  # keep memory reasonable downstream
    return M, adata.var_names.to_numpy()

def _as_category(series):
    cat = pd.Categorical(series)
    return cat.codes.astype(int), list(cat.categories)

def _encode_inputs(df_gene):
    # condition (control -> 0, mtDSB -> 1)
    cond = df_gene["condition"].astype(str).str.lower().str.strip()
    cond_map = {"control": 0, "ctrl": 0, "mtdsb": 1, "dsb": 1}
    cond_idx = cond.map(cond_map)
    if cond_idx.isna().any():
        bad = sorted(cond[cond_idx.isna()].unique().tolist())
        raise ValueError(f"Unknown condition labels: {bad}")

    age_idx, age_levels = _as_category(df_gene["age"].astype(str))
    sample_idx, sample_levels = _as_category(df_gene["sample"].astype(str))

    return {
        "y": df_gene["value"].to_numpy(dtype=float),
        "age_idx": age_idx,
        "cond_idx": cond_idx.to_numpy(dtype=int),
        "sample_idx": sample_idx,
        "A": len(age_levels),
        "S": len(sample_levels),
        "age_levels": age_levels,
        "sample_levels": sample_levels,
    }

def fit_gene_pymc(df_gene, draws=1200, tune=1200, target_accept=0.9, seed=1):
    """
    y ~ Normal(mu, sigma)
    mu = alpha_age[age] + beta_age[age] * cond + gamma_sample[sample]
    Priors:
        alpha_age[a] ~ Normal(0, 2)
        beta_age[a]  ~ Normal(0, 1)
        gamma[s]     ~ Normal(0, sigma_gamma)
        sigma_gamma  ~ HalfNormal(1)
        sigma        ~ HalfNormal(1)
    """
    D = _encode_inputs(df_gene)

    with pm.Model() as m:
        alpha_age   = pm.Normal("alpha_age", mu=0.0, sigma=2.0, shape=D["A"])
        beta_age    = pm.Normal("beta_age",  mu=0.0, sigma=1.0, shape=D["A"])
        sigma_gamma = pm.HalfNormal("sigma_gamma", sigma=1.0)
        gamma       = pm.Normal("gamma", mu=0.0, sigma=sigma_gamma, shape=D["S"])
        sigma       = pm.HalfNormal("sigma", sigma=1.0)

        mu = alpha_age[D["age_idx"]] + beta_age[D["age_idx"]] * D["cond_idx"] + gamma[D["sample_idx"]]
        pm.Normal("obs", mu=mu, sigma=sigma, observed=D["y"])

        idata = pm.sample(
            draws=draws, tune=tune, chains=4, target_accept=target_accept,
            random_seed=seed, return_inferencedata=True, progressbar=False
        )
    return idata, D

def summarize_gene(idata, D, gene_name):
    post = idata.posterior
    betas = post["beta_age"].stack(sample=("chain", "draw")).values  # A x N
    ages = D["age_levels"]

    def hdi(v, prob=0.95):
        lo, hi = az.hdi(v, hdi_prob=prob)
        return float(lo), float(hi)

    rows = []
    for a in range(len(ages)):
        v = betas[a, :]
        rows.append({
            "gene": gene_name,
            "term": f"beta[{ages[a]}] (mtDSB effect @ age)",
            "age": ages[a],
            "mean": float(v.mean()),
            "hdi_2.5%": hdi(v)[0],
            "hdi_97.5%": hdi(v)[1],
        })
    if len(ages) >= 2:
        diff = betas[-1, :] - betas[0, :]
        lo, hi = hdi(diff)
        rows.append({
            "gene": gene_name,
            "term": f"age_interaction({ages[-1]} - {ages[0]})",
            "age": f"{ages[-1]}-{ages[0]}",
            "mean": float(diff.mean()),
            "hdi_2.5%": lo,
            "hdi_97.5%": hi,
        })
    # diagnostics (quick overview)
    rhat = az.rhat(idata).to_array().mean().item()
    ess  = az.ess(idata).to_array().mean().item()
    for r in rows:
        r["rhat_mean"] = rhat
        r["ess_mean"] = ess
    return pd.DataFrame(rows)

def _list_parquet_files(path):
    """
    Handle both:
    - a single .parquet file
    - a directory containing many parquet shards (*.parquet)

    Returns list of parquet file paths.
    """
    if os.path.isdir(path):
        files = glob.glob(os.path.join(path, "*.parquet"))
        if len(files) == 0:
            raise FileNotFoundError(f"Directory {path} contains no .parquet files")
        return sorted(files)
    elif os.path.isfile(path):
        return [path]
    else:
        raise FileNotFoundError(f"{path} is not a file or directory")

def map_ids_to_symbols(ids, gene_map):
    # vectorized mapping via pandas for speed and NaN handling
    s = pd.Series(ids, dtype="string")
    base = s.str.replace(r"\.\d+$", "", regex=True)  # strip version suffixes
    symbols = base.map(gene_map).astype("string")
    return symbols.tolist()

def safe_mean(df, genes):
    """Mean across genes present in df (samples x genes)."""
    g = [g for g in genes if g in df.columns]
    if not g:
        return pd.Series(0.0, index=df.index)
    return df[g].mean(axis=1)

def score(mat, genes):
    g = [x for x in genes if x in mat.columns]
    return mat[g].mean(axis=1) if g else pd.Series(0, index=mat.index)

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

def wipe_neighbors(adata, key="neighbors"):
    if key in adata.uns:
        u = adata.uns[key]
        for kk in ("distances_key","connectivities_key"):
            if kk in u and u[kk] in adata.obsp:
                del adata.obsp[u[kk]]
        del adata.uns[key]
    for k in list(adata.obsp.keys()):
        if k.endswith(("distances","connectivities")) or k in ("distances","connectivities"):
            del adata.obsp[k]

def create_confusion_matrix(
        row_labels,
        col_labels):
    """
    Create a confusion matrix in which the numerical values are the Jaccard index between the corresponding
    cell type labels.

    Parameters
    ----------
    row_labels:
        an (n_cells, ) array listing the labels of each cell in the cell types that will be along the rows of
        the confusion matrix
    col_abels:
        an (n_cells, ) array listing the labels of each cell in the cell types that will be along the columns
        of the confusion matrix

    Returns
    -------
    A dict
        {'array': the (n_rows, n_cols) confusion matrix (jaccard index)
         'ct': the (n_rows, n_cols) count matrix (keeps track of how many cells are in each grid point)
         'rows': the label of each row (in the order they occur in the array)
         'cols': the label of each column (in the order they occur in the array)
         }
    """
    row_values = sorted(set(row_labels))
    col_values = sorted(set(col_labels))
    nrows = len(row_values)
    ncols = len(col_values)
    arr = np.zeros((nrows, ncols), dtype=float)
    ct_arr = np.zeros((nrows, ncols), dtype=float)
    row_to_idx = {r:ii for ii, r in enumerate(row_values)}
    col_to_idx = {c:ii for ii, c in enumerate(col_values)}
    row_idx = np.array([row_to_idx[r] for r in row_labels])
    col_idx = np.array([col_to_idx[c] for c in col_labels])
    for ir in range(nrows):
        r_set = set(np.where(row_idx==ir)[0])
        for ic in range(ncols):
            c_set = set(np.where(col_idx==ic)[0])
            jj = len(r_set.intersection(c_set))/len(r_set.union(c_set))
            arr[ir, ic] = jj
            ct_arr[ir, ic] = len(r_set.intersection(c_set))
    assert ct_arr.sum() == len(row_labels)
    assert len(row_labels) == len(col_labels)
    return {
        "array": arr,
        "ct": ct_arr,
        "rows": row_values,
        "cols": col_values
    }

def get_confusion_order(input_arr, ct_arr, rng, n_iter=500, inverse_temp=10.0, order_wgt=0.1):
    """
    Re-order the rows and columns of a confusion matrix so that it appears
    as diagnoal as possible. Use the Metropolis-Hastings algorithm to try
    to optimize a cost function.

    Parameters
    ----------
    input_arr:
        the confusion matrix. Probably the Jaccard index between
        two cel types
    ct_arr:
        an array the same shape as the confusion matrix showing how many
        cells are in each (row, col) pair of input_arr
    rng:
        a numpy random number generator
    n_iter:
        number of random steps to take with the Metropolis-Hastings
        algorithm
    inverse_temp:
        1/Temperature to use when deciding whether or not
        to accept a random step
    order_wgt:
        how important is ordering cell types by population.
        Higher values of order_wgt will privilege ordering cell types
        by n_cells over keeping the confusion matrix diagonal.

    Returns
    -------
    a dict
        {'rows': np.array of indexes of rows as they should be ordered
         'cols': np.array of indexes of cols as they should be ordered
         }
    """
    print("OPTIMIZING ORDER OF ROWS/COLUMNS IN CONFUSION MATRIX")
    def cost_fn(arr, dst, ct, order_wgt):

        val = (arr*dst).sum()
        row_sum = ct.sum(axis=1)
        row_sorted_idx = np.argsort(row_sum)[-1::-1]
        col_sum = ct.sum(axis=0)
        col_sorted_idx = np.argsort(col_sum)[-1::-1]

        desired_row = np.arange(arr.shape[0])
        row_term = ((desired_row-row_sorted_idx)**2).sum()
        desired_col = np.arange(arr.shape[1])
        col_term =((desired_col-col_sorted_idx)**2).sum()

        norm = ct.shape[0]+ct.shape[1]
        col_wgt = ct.shape[0]/norm
        row_wgt = ct.shape[1]/norm

        val += order_wgt*(col_wgt*col_term+row_wgt*row_term)

        return val

    baseline = np.copy(input_arr)
    candidate = np.copy(input_arr)

    candidate_ct = np.copy(ct_arr)
    baseline_ct = np.copy(ct_arr)

    grid = np.meshgrid(
        np.arange(input_arr.shape[0]),
        np.arange(input_arr.shape[1]),
        indexing='ij'
    )
    row_grid = grid[0]
    col_grid = grid[1]
    dst_grid = (row_grid-col_grid)**2
    print(f'raw cost {cost_fn(arr=input_arr, dst=dst_grid, ct=ct_arr, order_wgt=order_wgt):.2e}')

    rows = np.arange(input_arr.shape[0], dtype=int)
    cols = np.arange(input_arr.shape[1], dtype=int)

    # just try to put brightest pixels at top
    #row_max = baseline.max(axis=1)
    #assert row_max.shape == (baseline.shape[0], )
    #sorted_dex = np.argsort(row_max)[-1::-1]

    row_sum = ct_arr.sum(axis=1)
    sorted_dex = np.argsort(row_sum)[-1::-1]

    rows = rows[sorted_dex]
    baseline = baseline[rows, :]
    candidate = candidate[rows, :]
    baseline_ct = baseline_ct[rows, :]
    candidate_ct = candidate_ct[rows, :]
    print(f'cost after row shuffler {cost_fn(arr=baseline, dst=dst_grid, ct=baseline_ct, order_wgt=order_wgt):.2e}')

    #col_max = baseline.max(axis=0)
    #assert col_max.shape == (baseline.shape[1], )
    #sorted_dex = np.argsort(col_max)[-1::-1]

    col_sum = ct_arr.sum(axis=0)
    sorted_dex = np.argsort(col_sum)[-1::-1]

    cols = cols[sorted_dex]

    baseline = baseline[:, cols]
    candidate = candidate[:, cols]
    baseline_ct = baseline_ct[:, cols]
    candidate_ct = candidate_ct[:, cols]

    n_rows = input_arr.shape[0]
    n_cols = input_arr.shape[1]

    best_cost = cost_fn(arr=baseline, dst=dst_grid, ct=baseline_ct, order_wgt=order_wgt)
    actual_best_cost = best_cost
    actual_rows = np.copy(rows)
    actual_cols = np.copy(cols)
    print(f'base cost {best_cost}')
    for i_iter in range(n_iter):
        c0 = None
        c1 = None
        r0 = None
        r1 = None
        row_or_col = rng.integers(0, 2)
        if row_or_col == 0:
            to_swap = rng.choice(np.arange(n_rows), 2, replace=False)
            r0 = baseline[to_swap[0], :]
            r1 = baseline[to_swap[1], :]
            candidate[to_swap[1], :] = r0
            candidate[to_swap[0], :] = r1

            r0 = baseline_ct[to_swap[0], :]
            r1 = baseline_ct[to_swap[1], :]
            candidate_ct[to_swap[1], :] = r0
            candidate_ct[to_swap[0], :] = r1

        else:
            to_swap = rng.choice(np.arange(n_cols), 2, replace=False)
            c0 = baseline[:, to_swap[0]]
            c1 = baseline[:, to_swap[1]]
            candidate[:, to_swap[1]] = c0
            candidate[:, to_swap[0]] = c1

            c0 = baseline_ct[:, to_swap[0]]
            c1 = baseline_ct[:, to_swap[1]]
            candidate_ct[:, to_swap[1]] = c0
            candidate_ct[:, to_swap[0]] = c1

        candidate_cost = cost_fn(arr=candidate, dst=dst_grid, ct=candidate_ct, order_wgt=order_wgt)
        accept = False
        if candidate_cost < best_cost:
            accept = True
        else:
            roll = rng.random()
            delta = inverse_temp*(candidate_cost-best_cost)
            if np.exp(-0.5*(delta)) > roll:
                accept = True

        if accept:
            best_cost = candidate_cost
            if row_or_col == 0:
                baseline[to_swap[0], :] = candidate[to_swap[0], :]
                baseline[to_swap[1], :] = candidate[to_swap[1], :]
                baseline_ct[to_swap[0], :] = candidate_ct[to_swap[0], :]
                baseline_ct[to_swap[1], :] = candidate_ct[to_swap[1], :]

                r0 = rows[to_swap[0]]
                r1 = rows[to_swap[1]]
                rows[to_swap[1]] = r0
                rows[to_swap[0]] = r1
            else:
                baseline[:, to_swap[0]] = candidate[:, to_swap[0]]
                baseline[:, to_swap[1]] = candidate[:, to_swap[1]]

                baseline_ct[:, to_swap[0]] = candidate_ct[:, to_swap[0]]
                baseline_ct[:, to_swap[1]] = candidate_ct[:, to_swap[1]]

                c0 = cols[to_swap[0]]
                c1 = cols[to_swap[1]]
                cols[to_swap[1]] = c0
                cols[to_swap[0]] = c1
            if best_cost < actual_best_cost:
                actual_best_cost = best_cost
                actual_rows = np.copy(rows)
                actual_cols = np.copy(cols)
        else:
            if row_or_col == 0:
                candidate[to_swap[0], :] = baseline[to_swap[0], :]
                candidate[to_swap[1], :] = baseline[to_swap[1], :]

                candidate_ct[to_swap[0], :] = baseline_ct[to_swap[0], :]
                candidate_ct[to_swap[1], :] = baseline_ct[to_swap[1], :]
            else:
                candidate[:, to_swap[0]] = baseline[:, to_swap[0]]
                candidate[:, to_swap[1]] = baseline[:, to_swap[1]]

                candidate_ct[:, to_swap[0]] = baseline_ct[:, to_swap[0]]
                candidate_ct[:, to_swap[1]] = baseline_ct[:, to_swap[1]]

        if i_iter % (n_iter//10) == 0:
            print(f'iteration {i_iter} -- best_cost {actual_best_cost:.2e}')

    print(f'best_cost {actual_best_cost:.2e}')
    test = np.copy(input_arr)
    test = test[actual_rows, :]
    test = test[:, actual_cols]
    test_ct = np.copy(ct_arr)
    test_ct = test_ct[actual_rows, :]
    test_ct = test_ct[:, actual_cols]
    validate = cost_fn(arr=test, dst=dst_grid, ct=test_ct, order_wgt=order_wgt)
    print(f'validating {validate} -- {actual_best_cost}')

    return {'rows': actual_rows, 'cols': actual_cols}

def de_condition_within_age(ad, age_val):
    """mtDSB vs control at a given age."""
    ad_sub = ad[ad.obs['age'] == int(age_val)].copy()
    ad_sub = ad_sub[ad_sub.obs['condition'].isin(['control','mtdsb'])].copy()
    ad_sub.obs['condition'] = pd.Categorical(
        ad_sub.obs['condition'], categories=['control','mtdsb']
    )
    sc.tl.rank_genes_groups(
        ad_sub,
        groupby='condition',
        reference='control',
        method='wilcoxon',
        use_raw=True
    )
    df = sc.get.rank_genes_groups_df(ad_sub, group=None)
    df['contrast'] = f"mtDSB_vs_control@{age_val}"
    return df

def alpha_shape(points, alpha, only_outer=True):
    assert points.shape[0] > 3, "Need at least four points"

    def add_edge(edges, i, j):
        if (i, j) in edges or (j, i) in edges:
            if only_outer:
                edges.remove((j, i) if (j, i) in edges else (i, j))
            return
        edges.add((i, j))

    tri = Delaunay(points)
    edges = set()
    for ia, ib, ic in tri.simplices:
        pa, pb, pc = points[ia], points[ib], points[ic]
        a = np.linalg.norm(pa - pb)
        b = np.linalg.norm(pb - pc)
        c = np.linalg.norm(pc - pa)
        s = (a + b + c) / 2.0
        area = max(np.sqrt(s * (s - a) * (s - b) * (s - c)), 1e-12)
        circum_r = a * b * c / (4.0 * area)
        if circum_r < alpha:
            add_edge(edges, ia, ib)
            add_edge(edges, ib, ic)
            add_edge(edges, ic, ia)
    return edges

def stitch_boundaries(edges):
    edge_set = edges.copy()
    boundary_lst = []
    while len(edge_set) > 0:
        boundary = []
        edge0 = edge_set.pop()
        boundary.append(edge0)
        last_edge = edge0
        while True:
            i, j = last_edge
            j_first = [n for x, n in edge_set if x == j]
            j_second = [n for n, x in edge_set if x == j]
            next_j = j_first or j_second
            if not next_j:
                break
            k = next_j[0]
            edge_set.remove((j, k) if (j, k) in edge_set else (k, j))
            boundary.append((j, k))
            last_edge = (j, k)
            if boundary[0][0] == last_edge[1]:
                break
        boundary_lst.append(boundary)
    return boundary_lst

def generate_smoothed_alpha_polygon(coords_filtered, alpha=100, closing_radius=20, img_size=(4096, 4096)):
    edges = alpha_shape(coords_filtered, alpha=alpha, only_outer=True)
    stitched = stitch_boundaries(edges)
    if not stitched:
        return None

    largest_boundary = max(stitched, key=len)
    points = np.array([coords_filtered[i] for i, _ in largest_boundary])
    points = np.vstack([points, points[0]])  # close loop

    x_min, y_min = coords_filtered.min(axis=0)
    x_max, y_max = coords_filtered.max(axis=0)
    scale_x = img_size[1] / (x_max - x_min)
    scale_y = img_size[0] / (y_max - y_min)

    scaled_points = np.column_stack([
        (points[:, 1] - y_min) * scale_y,
        (points[:, 0] - x_min) * scale_x
    ])
    mask = polygon2mask(img_size, scaled_points)
    closed_mask = binary_closing(mask, iterations=closing_radius)
    contours = find_contours(closed_mask.astype(float), 0.5)
    if not contours:
        return None

    contour = contours[0]
    poly_points = np.column_stack([
        contour[:, 1] / scale_x + x_min,
        contour[:, 0] / scale_y + y_min
    ])
    return Polygon(poly_points)

def print_column_info(df):

    for c in df.columns:
        grouped = df[[c]].groupby(c).count()
        members = ''
        if len(grouped) < 30:
            members = str(list(grouped.index))
        print("Number of unique %s = %d %s" % (c, len(grouped), members))

def despine(ax):
    for sp in ["top","right","left","bottom"]:
        ax.spines[sp].set_visible(False)
    ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)

def _report_dups(name, idx):
    dup = pd.Index(idx)
    d = dup[dup.duplicated()].unique().tolist()
    if d:
        print(f"[WARN] {name}: {len(d)} duplicated names e.g. {d[:5]}")
    else:
        print(f"[OK] {name}: all unique")

def attach_counts_intersection_safe(adata, adata_raw, layer_name="counts", collapse_var=False):
    # 1) report duplicates
    _report_dups("adata.obs_names", adata.obs_names)
    _report_dups("adata.var_names", adata.var_names)
    _report_dups("adata_raw.obs_names", adata_raw.obs_names)
    _report_dups("adata_raw.var_names", adata_raw.var_names)

    # 2) handle duplicates
    if collapse_var:
        # collapse duplicated var names by summing columns (good for counts)
        def _collapse(A, names):
            df = pd.DataFrame.sparse.from_spmatrix(A) if sparse.issparse(A) else pd.DataFrame(A)
            df.columns = names
            return df.groupby(level=0, axis=1).sum()
        # collapse both objects on var (genes)
        A1 = _collapse(adata.X, adata.var_names)
        A2 = _collapse(adata_raw.X, adata_raw.var_names)
        # keep obs in original order
        A1.index = adata.obs_names
        A2.index = adata_raw.obs_names
        # replace X and var_names with collapsed
        adata = adata[:, []].copy()
        adata.X = A1.values
        adata.var_names = A1.columns.astype(str)
        adata.obs_names = A1.index.astype(str)

        adata_raw = adata_raw[:, []].copy()
        adata_raw.X = A2.values
        adata_raw.var_names = A2.columns.astype(str)
        adata_raw.obs_names = A2.index.astype(str)
    else:
        # just force uniqueness by appending suffixes
        adata.var_names_make_unique()
        adata.obs_names_make_unique()
        adata_raw.var_names_make_unique()
        adata_raw.obs_names_make_unique()

    # 3) intersect axes
    common_obs = adata.obs_names.intersection(adata_raw.obs_names)
    common_var = adata.var_names.intersection(adata_raw.var_names)

    ad_sub  = adata[common_obs, common_var].copy()
    raw_sub = adata_raw[common_obs, common_var]

    # 4) attach raw into a layer
    Xraw = raw_sub.X
    if sparse.issparse(Xraw):
        Xraw = Xraw.copy()
    else:
        Xraw = np.asarray(Xraw).copy()
    ad_sub.layers[layer_name] = Xraw

    print(f"[attach_counts] kept {ad_sub.n_obs} cells and {ad_sub.n_vars} genes "
          f"(intersected). Raw stored in .layers['{layer_name}'].")
    return ad_sub
