#!/usr/bin/env python
# coding: utf-8

# # mtDSB x Age DE analysis (pydeseq2, LRT-free)

# In[100]:


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
        reg = row["region"]
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
                & (obs["region"] == reg)
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
import scanpy as sc


# ## Read data

# In[101]:


comp_annos = sc.read_h5ad('../data/mtDNA_DSB_5k_clustered_manual_annotation_compartment.h5ad')


# In[102]:


adata = sc.read_h5ad('../data/mtDNA_DSB_5k_clustered_manual_annotation_with_raw.h5ad')


# In[103]:


adata.obs['region'] = adata.obs.index.map(dict(zip(comp_annos.obs.index, comp_annos.obs.compartment)))


# In[104]:


adata.obs["compartment"] = (
    adata.obs["region"]
    .replace({"Olfactory areas": "Olfactory areas + Hippocampal formation"})
    .astype("category")
)


# In[105]:


adata.obs


# In[106]:


adata.X = adata.layers['counts']
adata.X = adata.X.astype(int)


# In[154]:


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


# In[108]:


adata = adata[adata.obs.cell_class.isin(glia_classes)]


# In[109]:


adata


# In[110]:


for run in adata.obs['sample_id'].unique():
    print(run)
    ad_int = adata[adata.obs['sample_id'] == run]

    with plt.rc_context({'figure.figsize': (20, 10)}):
        sc.pl.spatial(ad_int, spot_size=20, color = 'compartment')
    plt.show()



# ## Perform pseudobulk

# In[111]:


pb = pseudobulk_by(
    adata,
    groupby=["cell_class", "region", "age", "condition", "sample_id"],  # 👈 add sample_id
    layer="counts" if "counts" in adata.layers else None,
    min_cells=30,
)

# Inspect one example
cell_type = "Mature oligodendrocytes"
region_age = list(pb[cell_type].keys())[0]
counts_df, meta_df = pb[cell_type][region_age]

print(region_age)
print("Counts shape:", counts_df.shape)
print(meta_df.head())


# In[112]:


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


# In[113]:


counts_df, meta_df = pb_class["Mature oligodendrocytes"]
print("columns match index:", counts_df.columns.equals(meta_df.index))


# In[114]:


results_by_region = run_by_region(
    pb_class,
    condition_col="condition",
    age_col="age",
    region_col="region",
    ref_condition="control",
    min_per_condition=2,
    n_cpus=8,
    verbose=True
)


# In[115]:


results_by_region


# In[167]:


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


# In[168]:


region_dict = {}
for glia in glia_classes:
    print(glia)
    rank_regions = summarize_region_impact(results_by_region,glia)
    print(rank_regions.head(10))
    region_dict[glia] = rank_regions


# In[169]:


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


# In[178]:


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


# In[179]:


region_dict


# In[195]:


top_genes_dotplot = ['Serpina3n', 'Gfap', 'Hcrt', 'Mt2', 'Sst', 'Mt2', 'Dlk1', 'Cd74', 'Pmch',
       'Ecel1', 'Syt4', 'Fezf1', 'Pdyn', 'Mt2', 'Glra1', 'C4b', 'Foxd2',
       'Igfbp3', 'Mt1', 'Ecel1']


# In[191]:


adata_mTDSB = adata[adata.obs.condition == 'mtDSB']
sc.tl.rank_genes_groups(adata_mTDSB, groupby='compartment', method='t-test')
# See top 5 marker genes per cluster
sc.pl.rank_genes_groups(adata_mTDSB, n_genes=25, sharey=False)

marker_genes = pd.DataFrame({
    group: adata_mTDSB.uns['rank_genes_groups']['names'][group][:10]
    for group in adata_mTDSB.uns['rank_genes_groups']['names'].dtype.names
})
marker_genes.head()
#marker_genes.to_csv('../data/broad_markers_leiden0-5.csv')
for col in marker_genes.columns: 
    print(col)
    genes = marker_genes[col].tolist()
    print(" ".join(genes))
    print(' ')


# In[187]:


mat, df_all = plot_integrated_region_scores(
    region_dict,
    cmap="rocket_r",
    normalize=False,
    #save="figs/glia_region_impact"
)


# In[173]:


plot_region_impact_for_glia(
    results_by_region=results_by_region,
    glia_classes=glia_classes,
    summarize_region_impact=summarize_region_impact,  # your existing function
    top_k=10,
    # score_col="sum_neglog10FDR",   # optionally force a column name
    #save_prefix="figs/region_impact"
)


# In[118]:


for region in rank_regions.head(10).region:
    res_hypo_60 = results_by_region["Olfactory astrocytes"][region]["simple_effects"]["60"]
    volcano_clean(res_hypo_60, title=region + " (60) — mtDSB vs Control", figsize=(5, 5))


# In[119]:


for region in rank_regions.head(10).region:
    res_hypo_60 = results_by_region["Mature oligodendrocytes"][region]["simple_effects"]["60"]
    volcano_clean(res_hypo_60, title=region + " (60) — mtDSB vs Control", figsize=(5, 5))


# In[121]:


top_regions = rank_regions.head(15)["region"].tolist()
print("Top affected regions:", top_regions)


# In[122]:


import numpy as np
import pandas as pd

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


# In[123]:


rank_regions


# In[189]:





# In[125]:


top_regions = rank_regions["region"].tolist()
top_genes = extract_top_genes(results_by_region, "Mature oligodendrocytes", regions=top_regions)
top_genes = top_genes.reset_index()

summary = (
    top_genes.groupby('index')
    .agg(
        n_regions=("region", "nunique"),          # number of regions where gene significant
        mean_LFC=("log2FC", "mean"),              # mean log2 fold change across them
        mean_neglog10FDR=("neglog10FDR", "mean")  # average -log10(FDR)
    )
    .sort_values("n_regions", ascending=False)
)

# your summary df is called `summary` (index=gene; cols: n_regions, mean_LFC, mean_neglog10FDR)
summary = summary.rename(columns={
    summary.columns[0]: "n_regions",
    summary.columns[1]: "mean_LFC",
    summary.columns[2]: "mean_neglog10FDR"
})


from IPython.display import display

summary_styled = (
    summary.head(30)
    .style
    .background_gradient(subset="n_regions", cmap="Purples")
    .bar(subset="mean_LFC", color="lightblue")
    .format({
        "mean_LFC": "{:.2f}",
        "mean_neglog10FDR": "{:.1f}"
    })
    .set_caption("Top 30 genes by number of regions detected")
    .set_table_styles([
        {'selector': 'caption', 'props': [('font-size', '16px'), ('font-weight', 'bold')]}
    ])
)

display(summary_styled)


# In[51]:


top_regions


# In[52]:


top_genes


# In[131]:


from IPython.display import display
import pandas as pd

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


# In[136]:


ct = "Microglia"
summary, styled = summarize_de_by_region(
    results_by_region,
    cell_class=ct,
    regions=top_regions,          # or None for all
    ages=["21","60"],
    padj_thr=1.0,
    lfc_thr=0.0,
    top_n=500,
    top_k_display=30,
    caption=ct+": top 30 genes by number of regions detected",
    #save_csv_path="mature_oligodendrocytes_summary.csv",
)


# In[138]:


ct = "Olfactory astrocytes"
summary, styled = summarize_de_by_region(
    results_by_region,
    cell_class=ct,
    regions=top_regions,          # or None for all
    ages=["21","60"],
    padj_thr=1.0,
    lfc_thr=0.0,
    top_n=500,
    top_k_display=30,
    caption=ct+": top 30 genes by number of regions detected",
    #save_csv_path="mature_oligodendrocytes_summary.csv",
)


# In[139]:


ct = "Mature oligodendrocytes"
summary, styled = summarize_de_by_region(
    results_by_region,
    cell_class=ct,
    regions=top_regions,          # or None for all
    ages=["21","60"],
    padj_thr=1.0,
    lfc_thr=0.0,
    top_n=500,
    top_k_display=30,
    caption=ct+": top 30 genes by number of regions detected",
    #save_csv_path="mature_oligodendrocytes_summary.csv",
)


# In[140]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

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


# In[148]:


ct = "Mature oligodendrocytes"
# After you run summarize_de_by_region(...)
summary, _ = summarize_de_by_region(
    results_by_region,
    cell_class=ct,
    regions=top_regions,
    caption=ct + ': top genes by number of regions detected', 
    ages=["21","60"],
    padj_thr=1.0, lfc_thr=0.0, top_n=500, top_k_display=30
)

# 1) Bars (regions + mean LFC)
plot_de_summary_bars(summary, top_k=30,
                     title=ct + ": recurrent DE genes across regions",
                     #save="summary_bars.png"
                    )

# 2) Bubble plot (nice single-panel summary)
plot_de_summary_bubbles(summary, top_k=30,
                        title=ct + ": recurrent genes: effect size vs recurrence",
                        #save="summary_bubbles.png
                        )

# 3) Heatmap of presence/sign by region (needs the long table)
top_genes = extract_top_genes(results_by_region, ct,
                              regions=top_regions, ages=["21","60"],
                              padj_thr=1.0, lfc_thr=0.0, top_n=500)
plot_region_presence_heatmap(top_genes, top_k=30,
                             title=ct + ": sign of change per region",
                             #save="region_presence_heatmap.png"
                            )


# In[147]:


ct = "Microglia"
# After you run summarize_de_by_region(...)
summary, _ = summarize_de_by_region(
    results_by_region,
    cell_class=ct,
    regions=top_regions,
    caption=ct + ': top genes by number of regions detected', 
    ages=["21","60"],
    padj_thr=1.0, lfc_thr=0.0, top_n=500, top_k_display=30
)

# 1) Bars (regions + mean LFC)
plot_de_summary_bars(summary, top_k=30,
                     title=ct + ": recurrent DE genes across regions",
                     #save="summary_bars.png"
                    )

# 2) Bubble plot (nice single-panel summary)
plot_de_summary_bubbles(summary, top_k=30,
                        title=ct + ": recurrent genes: effect size vs recurrence",
                        #save="summary_bubbles.png
                        )

# 3) Heatmap of presence/sign by region (needs the long table)
top_genes = extract_top_genes(results_by_region, ct,
                              regions=top_regions, ages=["21","60"],
                              padj_thr=1.0, lfc_thr=0.0, top_n=500)
plot_region_presence_heatmap(top_genes, top_k=30,
                             title=ct + ": sign of change per region",
                             #save="region_presence_heatmap.png"
                            )


# In[146]:


ct = "Olfactory astrocytes"
# After you run summarize_de_by_region(...)
summary, _ = summarize_de_by_region(
    results_by_region,
    cell_class=ct,
    regions=top_regions,
    caption=ct + ': top genes by number of regions detected', 
    ages=["21","60"],
    padj_thr=1.0, lfc_thr=0.0, top_n=500, top_k_display=30
)

# 1) Bars (regions + mean LFC)
plot_de_summary_bars(summary, top_k=30,
                     title=ct + ": recurrent DE genes across regions",
                     #save="summary_bars.png"
                    )

# 2) Bubble plot (nice single-panel summary)
plot_de_summary_bubbles(summary, top_k=30,
                        title=ct + ": recurrent genes: effect size vs recurrence",
                        #save="summary_bubbles.png
                        )

# 3) Heatmap of presence/sign by region (needs the long table)
top_genes = extract_top_genes(results_by_region, ct,
                              regions=top_regions, ages=["21","60"],
                              padj_thr=1.0, lfc_thr=0.0, top_n=500)
plot_region_presence_heatmap(top_genes, top_k=30,
                             title=ct + ": sign of change per region",
                             #save="region_presence_heatmap.png"
                            )


# In[141]:


# After you run summarize_de_by_region(...)
summary, _ = summarize_de_by_region(
    results_by_region,
    cell_class="Mature oligodendrocytes",
    regions=top_regions,
    ages=["21","60"],
    padj_thr=1.0, lfc_thr=0.0, top_n=500, top_k_display=30
)

# 1) Bars (regions + mean LFC)
plot_de_summary_bars(summary, top_k=30,
                     title="Mature OLs — recurrent DE genes across regions",
                     save="summary_bars.png")

# 2) Bubble plot (nice single-panel summary)
plot_de_summary_bubbles(summary, top_k=30,
                        title="Recurrent genes: effect size vs recurrence",
                        save="summary_bubbles.png")

# 3) Heatmap of presence/sign by region (needs the long table)
top_genes = extract_top_genes(results_by_region, "Mature oligodendrocytes",
                              regions=top_regions, ages=["21","60"],
                              padj_thr=1.0, lfc_thr=0.0, top_n=500)
plot_region_presence_heatmap(top_genes, top_k=30,
                             title="Sign of change per region",
                             save="region_presence_heatmap.png")


# In[91]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

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


# In[92]:


mat60 = region_gene_matrix(results_by_region, "Mature oligodendrocytes", gene_panel, ages=["60"], padj_thr=0.05)
mat21 = region_gene_matrix(results_by_region, "Mature oligodendrocytes", gene_panel, ages=["21"], padj_thr=0.05)

heatmap(mat60, "Region × gene (60 wks) — significant log2FC")
heatmap(mat21, "Region × gene (21 wks) — significant log2FC")


# In[93]:


import os
import pandas as pd

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


# In[94]:


save_all_region_results(results_by_region, outdir="../results/region")


# In[122]:


adata.obs


# In[ ]:




