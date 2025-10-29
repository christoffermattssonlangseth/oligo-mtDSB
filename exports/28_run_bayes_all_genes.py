#!/usr/bin/env python
# coding: utf-8

# In[3]:


import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp


# In[4]:


adata = sc.read_h5ad('../data/mtDNA_DSB_5k_clustered_manual_annotation_with_raw.h5ad')


# In[22]:


# ==== CONFIG: set your obs column names here ====
AGE_COL       = "age"         # e.g., 21 / 60 or "21" / "60"
COND_COL      = "condition"   # values like "control" / "mtDSB"
SAMPLE_COL    = "sample_id"   # section/mouse/slide identifier
LAYER_COUNTS  = "counts"      # or None to use .X
NORM_METHOD   = "cpm_log1p"   # "cpm_log1p" (recommended) or "log1p_sizefactor"

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

def pseudobulk_by_groups(adata, group_cols=(SAMPLE_COL, AGE_COL, COND_COL), layer=LAYER_COUNTS):
    """
    Sum counts per unique (sample, age, condition). Sparse-safe and reasonably fast.
    Returns:
      pb_counts: (G x genes) dense float32 array
      groups_df: DataFrame with columns [SAMPLE_COL, AGE_COL, COND_COL]
      var_names: gene names (np.array)
    """
    M, var_names = _get_counts_matrix(adata, layer=layer)
    # Build a single categorical key for grouping rows
    key_df = adata.obs[list(group_cols)].astype(str).copy()
    key_df["_group_key"] = key_df.apply(lambda r: "||".join(r.values.tolist()), axis=1)
    codes, groups = pd.factorize(key_df["_group_key"], sort=True)
    n_groups = len(groups)

    # Sum rows per group into a new (n_groups x n_genes) matrix
    # Efficient loop over groups (using CSR slicing and summation)
    pb = np.zeros((n_groups, M.shape[1]), dtype=np.float64)
    for g in range(n_groups):
        idx = np.where(codes == g)[0]
        if idx.size == 0:
            continue
        block_sum = M[idx].sum(axis=0)  # 1 x genes
        pb[g, :] = np.asarray(block_sum).ravel()

    # Decode group columns back out
    groups_df = pd.DataFrame([x.split("||") for x in groups], columns=group_cols)
    return pb, groups_df.reset_index(drop=True), var_names

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

# FIXED: melt with cell type kept as id_var, and only gene columns melted.
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


# In[23]:


import pymc as pm
import arviz as az

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


# In[ ]:


adata


# In[11]:


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


# In[12]:


adata = adata[adata.obs.cell_class.isin(glia_classes)]


# In[24]:


CELLTYPE_COL = "cell_class"

# If you grouped by (cell_class, sample_id, age, condition):
pb_counts, groups_df, var_names = pseudobulk_by_groups(
    adata, group_cols=("cell_class", "sample_id", "age", "condition"),
    layer=LAYER_COUNTS
)
pb_norm = normalize_pseudobulk(pb_counts, method="cpm_log1p", axis=1)

df_long = pseudobulk_to_long_with_celltype(
    pb_norm, groups_df, var_names,
    sample_col="sample_id", age_col="age",
    cond_col="condition", celltype_col="cell_class",
    value_name="value"
)


# In[32]:


celltypes


# In[ ]:


celltypes = sorted(df_long["cell_type"].unique())

results_all = []
for ct in ['Mature oligodendrocytes']:# celltypes:
    sub_ct = df_long[df_long["cell_type"] == ct]


    print(f"Running model for {ct}...")
    summary_ct = fit_all_genes(
        sub_ct,
        genes=None,
        draws=1000,
        tune=1000,
        target_accept=0.9
    )
    summary_ct["cell_type"] = ct
    results_all.append(summary_ct)

summary_all = pd.concat(results_all, ignore_index=True)


# In[ ]:


summary_all[summary_all.gene == 'Slc16a1']


# In[ ]:


credible = summary_all[
    (summary_all["term"].str.contains("age_interaction")) &
    (summary_all["hdi_2.5%"] * summary_all["hdi_97.5%"] > 0)
]
credible.sort_values("mean", ascending=False).head(20)


# In[ ]:


print('DONE')


# In[ ]:




