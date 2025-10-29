#!/usr/bin/env python
# coding: utf-8

# # processing eae 5k

# ## load packages

# In[1]:


import warnings
warnings.filterwarnings('ignore')
import os
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.cluster import KMeans


# ## data retrival

# In[5]:


import os
import scanpy as sc
import pandas as pd

base_dir = '/Volumes/processing2/20250917__105547__GoncaloLeslie_5kMouse_run4/'

# Grab only folders with "RB" in the name
runs = [d for d in os.listdir(base_dir) if d.startswith("output-") and os.path.isdir(os.path.join(base_dir, d))]

ad_list = []

for run in runs:
    run_path = os.path.join(base_dir, run)

    # Define file paths directly in each run folder
    h5_path = os.path.join(run_path, 'cell_feature_matrix.h5')
    cell_info_path = os.path.join(run_path, 'cells.csv.gz')

    if not (os.path.exists(h5_path) and os.path.exists(cell_info_path)):
        print(f"Skipping {run_path} (missing required files)")
        continue

    print(f"Loading: {run}")
    ad_int = sc.read_10x_h5(h5_path)
    cell_info = pd.read_csv(cell_info_path, index_col=0)

    ad_int.obs = cell_info
    ad_int.obs['run'] = run

    ad_list.append(ad_int)


# In[6]:


ad = sc.concat(ad_list)


# In[7]:


ad.obs


# ## write raw data

# In[13]:


ad.write('../data/mtDNA_DSB_5k_raw.h5ad')


# In[16]:


ad.obs = ad.obs.reset_index()


# In[14]:


sc.pp.calculate_qc_metrics(ad, percent_top=None, log1p=False, inplace=True)


# In[34]:


ad.obs['sample_id'] = ad.obs.run.str.split('__', expand = True)[2]


# In[24]:


import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# --- Load your Xenium dataset ---
# --- Example: basic summary table ---
summary = (
    ad.obs.groupby("sample_id")
    .agg(
        n_cells=("cell_id", "count"),
        mean_counts=("total_counts", "mean"),
        mean_genes=("n_genes_by_counts", "mean"),
    )
    .sort_values("n_cells", ascending=False)
)
display(summary)

# --- Plot 1: Cell count per sample ---
plt.figure(figsize=(8,4))
sns.barplot(y=summary.index, x=summary["n_cells"], palette="crest")
plt.title("Number of segmented cells per Xenium experiment")
plt.xlabel("Cell count")
plt.ylabel("Sample")
plt.tight_layout()
plt.show()

# --- Plot 2: Mean detected genes per sample ---
plt.figure(figsize=(8,4))
sns.barplot(y=summary.index, x=summary["mean_genes"], palette="mako")
plt.title("Mean number of genes per cell")
plt.xlabel("Mean n_genes_by_counts")
plt.tight_layout()
plt.show()

# --- Plot 3: Distributions of total counts ---
plt.figure(figsize=(6,4))
sns.histplot(ad.obs["total_counts"], bins=50, kde=True)
plt.title("Distribution of total counts per cell")
plt.xlabel("total_counts")
plt.tight_layout()
plt.show()

# --- Plot 4: Composition of cell types ---
if "cell_types" in ad.obs:
    ct_counts = ad.obs["cell_types"].value_counts()
    plt.figure(figsize=(8,4))
    sns.barplot(y=ct_counts.index[:15], x=ct_counts.values[:15], palette="viridis")
    plt.title("Top 15 cell types (by abundance)")
    plt.xlabel("Number of cells")
    plt.tight_layout()
    plt.show()

# --- Plot 5: Cell type composition per sample ---
if {"cell_types","sample_id"}.issubset(ad.obs.columns):
    comp = (
        ad.obs.groupby(["sample_id","cell_types"])
        .size()
        .groupby(level=0)
        .apply(lambda x: x / x.sum())
        .unstack(fill_value=0)
    )
    plt.figure(figsize=(10,5))
    comp.plot(kind="bar", stacked=True, colormap="tab20", figsize=(10,5))
    plt.ylabel("Fraction of cells")
    plt.title("Cell type composition per Xenium experiment")
    plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    plt.show()


# In[28]:


# Pretty print with formatting
styled = (
    summary.style
    .background_gradient(subset=["n_cells"], cmap="Purples")
    .background_gradient(subset=["counts_mean"], cmap="Blues")
    .background_gradient(subset=["genes_mean"], cmap="Greens")
    .format({
        "n_cells": "{:,.0f}",
        "counts_mean": "{:.1f}",
        "counts_median": "{:.0f}",
        "counts_p10": "{:.0f}",
        "counts_p90": "{:.0f}",
        "genes_mean": "{:.1f}",
        "genes_median": "{:.0f}",
        "genes_p10": "{:.0f}",
        "genes_p90": "{:.0f}",
    })
    .set_caption("📊 Xenium QC Summary per Sample")
)
styled


# In[27]:


# --- Xenium QC dashboard (drop-in) ---
import os, numpy as np, pandas as pd, scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

# ====== config (rename if your column names differ) ======
RUN_COL     = "run"            # batch/slide/run id
SAMPLE_COL  = "sample_id"      # sample or library id
CT_COL      = "cell_types"     # cell-type column (optional)
COUNTS_COL  = "total_counts"   # transcripts per cell
NGENES_COL  = "n_genes_by_counts"
SAVE_DIR    = "xenium_qc"
os.makedirs(SAVE_DIR, exist_ok=True)

# ====== 1) summary table (per run) ======
cols = [c for c in [RUN_COL, SAMPLE_COL, CT_COL, COUNTS_COL, NGENES_COL] if c in ad.obs.columns]
assert RUN_COL in ad.obs, f"Missing column '{RUN_COL}' in ad.obs"

def _p(x,p): return np.nanpercentile(x, p)

agg_dict = {"cell_id":("cell_id","count")} if "cell_id" in ad.obs else {"index":("index","count")}
if COUNTS_COL in ad.obs:
    agg_dict |= {
        "counts_mean":(COUNTS_COL,"mean"),
        "counts_median":(COUNTS_COL,"median"),
        "counts_p10":(COUNTS_COL,lambda x:_p(x,10)),
        "counts_p90":(COUNTS_COL,lambda x:_p(x,90)),
    }
if NGENES_COL in ad.obs:
    agg_dict |= {
        "genes_mean":(NGENES_COL,"mean"),
        "genes_median":(NGENES_COL,"median"),
        "genes_p10":(NGENES_COL,lambda x:_p(x,10)),
        "genes_p90":(NGENES_COL,lambda x:_p(x,90)),
    }

summary = (
    ad.obs.reset_index()
      .groupby(SAMPLE_COL)
      .agg(**agg_dict)
      .rename(columns={"cell_id":"n_cells","index":"n_cells"})
      .sort_values("n_cells", ascending=False)
)
display(summary)
summary.to_csv(f"{SAVE_DIR}/summary_by_run.csv")

# ====== 2) plots: cell counts & QC per run ======
plt.figure(figsize=(9,4.5))
sns.barplot(y=summary.index, x=summary["n_cells"], palette="Set3")
plt.title("Cells per Xenium run")
plt.xlabel("# cells"); plt.ylabel("Run")
plt.tight_layout(); plt.savefig(f"{SAVE_DIR}/cells_per_run_bar.png", dpi=200); plt.show()

if NGENES_COL in ad.obs:
    plt.figure(figsize=(10,4.5))
    sns.violinplot(data=ad.obs, x=SAMPLE_COL, y=NGENES_COL, inner="quartile", palette="rocket")
    plt.title("n_genes_by_counts per run"); plt.xlabel("Run"); plt.ylabel("n_genes_by_counts")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout(); plt.savefig(f"{SAVE_DIR}/ngenes_violin.png", dpi=200); plt.show()

if COUNTS_COL in ad.obs:
    plt.figure(figsize=(10,4.5))
    sns.violinplot(data=ad.obs, x=SAMPLE_COL, y=COUNTS_COL, inner="quartile", palette="mako")
    plt.title("total_counts per run"); plt.xlabel("Run"); plt.ylabel("total_counts")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout(); plt.savefig(f"{SAVE_DIR}/counts_violin.png", dpi=200); plt.show()

# counts vs genes hex
if set([COUNTS_COL, NGENES_COL]).issubset(ad.obs.columns):
    plt.figure(figsize=(6,5))
    plt.hexbin(ad.obs[COUNTS_COL], ad.obs[NGENES_COL], gridsize=50, mincnt=1)
    plt.xlabel("total_counts"); plt.ylabel("n_genes_by_counts")
    plt.title("Counts vs genes (all runs)")
    plt.tight_layout(); plt.savefig(f"{SAVE_DIR}/counts_vs_genes_hex.png", dpi=200); plt.show()

# ====== 3) cell-type composition ======
if CT_COL in ad.obs:
    ct_counts = ad.obs[CT_COL].value_counts()
    plt.figure(figsize=(9,5))
    sns.barplot(y=ct_counts.index[:20], x=ct_counts.values[:20], palette="Spectral")
    plt.title("Top cell types (all runs)"); plt.xlabel("# cells"); plt.ylabel("cell type")
    plt.tight_layout(); plt.savefig(f"{SAVE_DIR}/celltypes_top20.png", dpi=200); plt.show()

    # per run stacked fractions (collapse tail to 'Other' for readability)
    topK = 14
    ct_top = ct_counts.index[:topK].tolist()
    comp = (
        ad.obs.assign(ct_plot=lambda d: d[CT_COL].where(d[CT_COL].isin(ct_top), other="Other"))
              .groupby([SAMPLE_COL, "ct_plot"]).size()
              .groupby(level=0).apply(lambda x: x/x.sum())
              .unstack(fill_value=0)
    )
    ax = comp.plot(kind="bar", stacked=True, figsize=(10,5), colormap="tab20")
    plt.ylabel("fraction of cells"); plt.title("Cell-type composition per run")
    plt.legend(bbox_to_anchor=(1.02,1), loc="upper left", ncol=1)
    plt.tight_layout(); plt.savefig(f"{SAVE_DIR}/celltypes_per_run_stacked.png", dpi=200); plt.show()

# ====== 4) panel coverage (gene detection rate) ======
# fraction of cells with X>0 per gene (overall & per run)
try:
    import scipy.sparse as sp
    X = ad.X
    is_sparse = sp.issparse(X)
    # overall detection
    if is_sparse:
        detected = (X > 0).astype(np.int8)
        det_overall = np.array(detected.sum(axis=0)).ravel() / ad.n_obs
    else:
        det_overall = (X > 0).sum(axis=0) / ad.n_obs
    det_overall = pd.Series(det_overall, index=ad.var_names, name="fraction_cells")
    det_overall.sort_values(ascending=False).to_csv(f"{SAVE_DIR}/gene_detection_overall.csv")
    # top 30 genes by detection
    top30 = det_overall.sort_values(ascending=False).head(30)
    plt.figure(figsize=(8,5))
    sns.barplot(y=top30.index, x=top30.values, palette="coolwarm")
    plt.xlabel("fraction of cells detected"); plt.ylabel("gene"); plt.title("Panel coverage: top 30 genes")
    plt.tight_layout(); plt.savefig(f"{SAVE_DIR}/gene_detection_top30.png", dpi=200); plt.show()

    # per-run detection heatmap (top N informative genes)
    if SAMPLE_COL in ad.obs:
        # compute per run
        runs = ad.obs[SAMPLE_COL].astype(str).values
        run_idx = {r: np.where(runs == r)[0] for r in np.unique(runs)}
        det_run = {}
        for r, idx in run_idx.items():
            if len(idx)==0: continue
            if is_sparse:
                sub = X[idx,:]
                frac = np.array((sub > 0).sum(axis=0)).ravel() / len(idx)
            else:
                sub = X[idx,:]
                frac = (sub > 0).sum(axis=0) / len(idx)
            det_run[r] = frac
        det_df = pd.DataFrame(det_run, index=ad.var_names).T  # runs × genes

        # pick top genes by overall variance
        var_rank = det_df.var(axis=0).sort_values(ascending=False)
        sel = var_rank.head(40).index
        plt.figure(figsize=(min(12, len(sel)*0.3+4), max(4, len(det_df)*0.35+2)))
        sns.heatmap(det_df[sel], cmap="rocket", vmin=0, vmax=1, cbar_kws={"label":"fraction detected"})
        plt.title("Gene detection per run (top variable genes)")
        plt.xlabel("gene"); plt.ylabel("run")
        plt.tight_layout(); plt.savefig(f"{SAVE_DIR}/gene_detection_heatmap_runs.png", dpi=200); plt.show()
except Exception as e:
    print(f"[panel coverage] skipped: {e}")

# ====== 5) nice-to-haves (only if columns exist) ======
# Size/area QC (Xenium often provides cell area)
for cand in ["cell_area_um2","cell_area","area","nucleus_area_um2"]:
    if cand in ad.obs:
        plt.figure(figsize=(7,4))
        sns.violinplot(data=ad.obs, x=SAMPLE_COL, y=cand, inner="quartile", palette="PuBuGn")
        plt.title(f"{cand} per run"); plt.xlabel("Run"); plt.xticks(rotation=45, ha="right")
        plt.tight_layout(); plt.savefig(f"{SAVE_DIR}/{cand}_violin.png", dpi=200); plt.show()
        break

print(f"✓ saved outputs to {SAVE_DIR}/")


# ## preprocessing

# ### calculate qc metrics and filter cells for counts and number of genes

# In[14]:


ad.obs['segmentation_method'].value_counts()


# In[15]:


sc.pp.calculate_qc_metrics(ad, percent_top=None, log1p=False, inplace=True)
sc.pp.filter_cells(ad,min_counts=40)
sc.pp.filter_cells(ad,min_genes=15)


# ### normalizing and transforming

# In[16]:


sc.pp.normalize_total(ad, inplace=True,target_sum=100)
sc.pp.log1p(ad)
#sc.pp.scale(ad, )#max_value=10)


# ### pca and neighbors

# In[17]:


plt.rcdefaults()
sc.tl.pca(ad)
sc.pl.pca_variance_ratio(ad, n_pcs=50, log=True)
sc.pp.neighbors(ad, n_neighbors=15, n_pcs=30)


# ### umap

# In[18]:


sc.tl.umap(ad, min_dist=0.1)


# ### clustering

# In[42]:


resolutions = [0.5, 1,1.5, 2]

for resolution in resolutions:
    key = f'leiden_{resolution}'

    if key in ad.obs.columns:
        print(f"Skipping {resolution}: {key} already exists.")
    else:
        print(f"Clustering at resolution {resolution}...")
        sc.tl.leiden(ad, resolution=resolution, key_added=key)
        print("Done.")

    # plot UMAP
    sc.pl.umap(ad, color=key, legend_loc='on data', frameon=False)


# ### write clustered data

# In[39]:


ad.write('../data/mtDNA_DSB_5k_clustered.h5ad')


# In[29]:


ad = sc.read('../data/mtDNA_DSB_5k_clustered.h5ad')


# ### add spatial information

# In[30]:


spatial = np.array(ad.obs[['x_centroid','y_centroid']])
ad.obsm['spatial'] = spatial


# In[52]:


import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec

def plot_spatial_compact(
    ad,
    color="leiden_2",
    groupby="sample_id",
    spot_size=18,
    cols=3,
    height=8,                    # figure height (inches)
    legend_col_width=1.2,        # width (inches) reserved for legend
):
    sids = list(ad.obs[groupby].unique())
    n = len(sids)
    rows = int(np.ceil(n / cols))

    # compute figure width: panels + skinny legend column
    panel_w = height * cols * 0.6 / rows   # keeps a pleasing aspect
    fig_w = panel_w + legend_col_width
    fig = plt.figure(figsize=(fig_w, height), constrained_layout=False)

    # Grid: left (panels), right (legend)
    gs = GridSpec(
        rows, cols + 1, figure=fig,
        width_ratios=[1]*cols + [legend_col_width / (fig_w - legend_col_width)],
        wspace=0.02, hspace=0.02
    )

    # colors/categories (respect ad.uns if present)
    cats = ad.obs[color].astype("category").cat.categories
    if f"{color}_colors" in ad.uns:
        cols_list = ad.uns[f"{color}_colors"]
    else:
        base = sc.plotting.palettes.default_20
        cols_list = (base * int(np.ceil(len(cats)/len(base))))[:len(cats)]

    # panels
    for i, sid in enumerate(sids):
        r, c = divmod(i, cols)
        ax = fig.add_subplot(gs[r, c])
        ad_sub = ad[ad.obs[groupby] == sid].copy()

        sc.pl.spatial(
            ad_sub,
            color=color,
            spot_size=spot_size,
            show=False,
            ax=ax,
            legend_loc=None,
            frameon=False,     # no frame
            title=sid,         # short title; comment out to save more space
        )
        # remove axis labels Scanpy adds ("spatial1/2")
        ax.set_xlabel(""); ax.set_ylabel("")
        ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)

    # blank unused axes (if any)
    for j in range(n, rows*cols):
        r, c = divmod(j, cols)
        fig.add_subplot(gs[r, c]).axis("off")

    # legend in the skinny column
    ax_leg = fig.add_subplot(gs[:, -1])
    ax_leg.axis("off")
    handles = [
        Line2D([0],[0], marker='o', color='w',
               markerfacecolor=cols_list[k], markersize=7, label=str(cat))
        for k, cat in enumerate(cats)
    ]
    ax_leg.legend(handles=handles, title=color, frameon=False, loc="center left")

    # squeeze outer margins
    fig.subplots_adjust(left=0.01, right=0.98, top=0.98, bottom=0.02, wspace=0.02, hspace=0.02)
    plt.show()



# In[50]:


ad


# In[54]:


# usage
plot_spatial_compact(
    ad,
    color="leiden_2",
    groupby="sample_id",
    spot_size=20,
    cols=6,
    height=8,
    legend_col_width=1.0,
)


# ### plot clusters on basis of coordinates

# In[41]:


for run in ad.obs['run'].unique():
    print(run)
    ad_int = ad[ad.obs['run'] == run]
    with plt.rc_context({'figure.figsize': (20, 10)}):
        sc.pl.spatial(ad_int, spot_size=15, color = 'leiden_2')
    plt.show()



# In[ ]:




