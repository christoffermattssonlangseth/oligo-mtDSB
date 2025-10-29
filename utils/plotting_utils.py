"""
Auto-generated utilities for plotting_utils.
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

def plot_condition_effect_scatter(summary_df, cell_type, top_n=50, save=None):
    # filter one cell type
    d = summary_df[summary_df["cell_class"] == cell_type].copy()
    if d.empty:
        print(f"No data for {cell_type}")
        return

    # keep only mtDSB condition if others are present
    if "condition" in d.columns and d["condition"].nunique() > 1:
        d = d[d["condition"] == "mtDSB"]

    # rank by absolute mean effect
    d["abs_mean"] = d["mean"].abs()
    d = d.sort_values("abs_mean", ascending=False)

    plt.figure(figsize=(6, max(4, 0.25*len(d))))
    sns.scatterplot(data=d, x="mean", y="gene", s=40, color="C0", edgecolor="k", linewidth=0.3)
    # add HDI intervals as horizontal bars
    plt.hlines(y=d["gene"], xmin=d["hdi_2.5%"], xmax=d["hdi_97.5%"],
               color="gray", alpha=0.5, linewidth=2)

    plt.axvline(0, color="k", ls="--", lw=1)
    plt.title(f"{cell_type}\nPosterior mean ± 95% HDI (mtDSB vs control)")
    plt.xlabel("Posterior mean effect")
    plt.ylabel("")
    plt.tight_layout()

    if save:
        plt.savefig(save, dpi=200, bbox_inches="tight")
    plt.show()

def plot_top_genes_caterpillar(summary_df, cell_class, genes, save=None):
    d = summary_df[(summary_df["cell_class"] == cell_class) & (summary_df["gene"].isin(genes))].copy()
    if d.empty:
        print(f"No data for {cell_class}")
        return

    # sort by mean effect
    d = d.sort_values("mean", ascending=True)
    d["gene"] = pd.Categorical(d["gene"], categories=d["gene"], ordered=True)

    plt.figure(figsize=(3, 0.1 * len(d) + 2))
    ax = sns.pointplot(
        data=d, y="gene", x="mean", join=False, color="C0",
        errorbar=None, scale=0.6
    )
    plt.hlines(y=d["gene"], xmin=d["hdi_2.5%"], xmax=d["hdi_97.5%"],
               color="gray", lw=2, alpha=0.6)
    plt.axvline(0, color="k", ls="--", lw=1)
    plt.title(f"{cell_class} — mtDSB effect (posterior mean ± 95% HDI)")
    plt.xlabel("Posterior mean effect")
    plt.tight_layout()
    if save:
        plt.savefig(save, dpi=200, bbox_inches="tight")
    plt.show()

def post_merge(df, labels, post_merge_cutoff, linkage_method='single', 
               linkage_metric='correlation', fcluster_criterion='distance', name='', save=False):
    Z = scipy.cluster.hierarchy.linkage(df.T, method=linkage_method, metric=linkage_metric)
    merged_labels_short = scipy.cluster.hierarchy.fcluster(Z, post_merge_cutoff, criterion=fcluster_criterion)

    # map cluster name to merged label
    label_conversion = dict(zip(df.columns, merged_labels_short))
    new_labels = [label_conversion[i] for i in labels]

    # Plot
    fig, ax = plt.subplots(figsize=(20, 10))
    scipy.cluster.hierarchy.dendrogram(Z, labels=df.columns, color_threshold=post_merge_cutoff, ax=ax)
    ax.hlines(post_merge_cutoff, 0, ax.get_xlim()[1])
    ax.set_title('Merged clusters')
    ax.set_ylabel(linkage_metric, fontsize=20)
    ax.set_xlabel('Pre-merge cluster labels', fontsize=20)
    ax.tick_params(axis="x", labelrotation=90, labelsize=10)  # <-- rotate tick labels
    plt.show()
    if save:
        fig.savefig(f"{name}.svg", dpi=500)

    return new_labels

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

def _prep_age_df(df_age):
    """
    Sort regions so that highest activation_index ends up at the top visually.
    We sort ascending because barh plots from bottom to top.
    """
    return df_age.sort_values("activation_index", ascending=True).copy()

def plot_mtDSB_by_age(summary_df, cell_type, top_n=20, ascending=False):
    sub = summary_df[(summary_df["cell_type"]==cell_type) & (summary_df["effect"]=="mtDSB_at_age")]
    # pick the older age row per gene for ranking
    rank_age = sorted(sub["age"].astype(str).unique())[-1]
    ranker = sub[sub["age"]==rank_age].sort_values("mean", ascending=ascending).head(top_n).gene
    keep = sub[sub.gene.isin(ranker)]
    plt.figure(figsize=(10,4))
    ax = sns.pointplot(data=keep, x="gene", y="mean", hue="age", errorbar=None)
    for _, r in keep.iterrows():
        ax.plot([r.name%1+r["gene"], r.name%1+r["gene"]],
                [r["hdi_2.5%"], r["hdi_97.5%"]], c="k", lw=1, alpha=0.6)
    ax.axhline(0, ls="--", c="k", lw=1)
    plt.title(f"{cell_type}: mtDSB effect by age (top {top_n} by age={rank_age})")
    plt.ylabel("Effect (log scale)")
    plt.xticks(rotation=65, ha="right")
    plt.tight_layout()

def heatmap_across_celltypes(summary_df, effect="age_interaction", genes=None):
    dfp = summary_df[summary_df["effect"]==effect]
    if genes is not None:
        dfp = dfp[dfp["gene"].isin(genes)]
    mat = dfp.pivot(index="gene", columns="cell_type", values="mean").fillna(0)
    plt.figure(figsize=(7,10))
    sns.heatmap(mat, cmap="coolwarm", center=0)
    plt.title(f"{effect} (posterior means)")
    plt.tight_layout()

def plot_effect_heatmap(df, cell_type, genes=None, save=None):
    d = df[df["cell_class"]==cell_type].copy()
    d = d[d["effect_clean"].isin(["mtDSB @ age 21", "mtDSB @ age 60", "age interaction (60–21)"])]
    if genes:
        d = d[d["gene"].isin(genes)]
    # wide matrix genes × effects
    mat = d.pivot_table(index="gene", columns="effect_clean", values="mean", aggfunc="mean").fillna(0)
    if mat.empty:
        print(f"No rows to plot for {cell_type}")
        return
    mat = mat.loc[sorted(mat.index), ["mtDSB @ age 21","mtDSB @ age 60","age interaction (60–21)"]]

    plt.figure(figsize=(8, max(6, 0.35*len(mat))))
    sns.heatmap(mat, cmap="vlag", center=0, cbar_kws={"label":"mean effect (log)"})
    plt.title(f"{cell_type} — effects heatmap")
    plt.xlabel("effect"); plt.ylabel("gene")
    plt.tight_layout()
    if save: plt.savefig(save, dpi=200, bbox_inches="tight")
    plt.show()

def plot_raw_strip(df_long, cell_type, gene, save=None):
    d = df_long[(df_long["cell_class"]==cell_type) & (df_long["gene"]==gene)].copy()
    if d.empty:
        print(f"No raw rows for {cell_type} / {gene}")
        return

    plt.figure(figsize=(6,4))
    ax = sns.stripplot(
        data=d, x="age", y="value", hue="condition",
        dodge=True, alpha=0.7
    )

    # overlay group means
    means = d.groupby(["age","condition"])["value"].mean().reset_index()
    sns.pointplot(
        data=means, x="age", y="value", hue="condition",
        dodge=0.4, ci=None, markers="D", linestyles="",
        palette=sns.color_palette()[:2], ax=ax
    )

    plt.title(f"{gene} — raw pseudobulk (log1p-CPM)\n{cell_type}")
    plt.ylabel("expression")
    plt.xlabel("age")

    # move legend outside (right side)
    handles, labels = ax.get_legend_handles_labels()
    plt.legend(handles[:2], labels[:2], title="condition",
               bbox_to_anchor=(1.05, 0.5), loc="center left", borderaxespad=0.)

    plt.tight_layout(rect=[0, 0, 0.85, 1])  # leave room for legend
    if save:
        plt.savefig(save, dpi=200, bbox_inches="tight")
    plt.show()

def plot_21_vs_60(df, cell_type, label_top=15, save=None):
    # pick only the per-age mtDSB effects
    d = df[(df["cell_class"] == cell_type) & (df["effect"].str.startswith("mtDSB_at_age"))].copy()
    if d.empty:
        print(f"No mtDSB-per-age rows for {cell_type}")
        return

    # pivot: genes × age
    wide = d.pivot(index="gene", columns="age", values="mean")
    if wide.empty:
        print(f"Pivot produced empty table for {cell_type}")
        return

    # normalize column names to strings without spaces, then rename to age21/age60
    colmap = {}
    for c in wide.columns:
        cs = str(c).strip().lower().replace("age", "").strip()  # -> '21' or '60'
        if cs in {"21", "60"}:
            colmap[c] = f"age{cs}"
    wide = wide.rename(columns=colmap)

    # if still missing, try direct integer keys
    if "age21" not in wide.columns and 21 in wide.columns:  wide = wide.rename(columns={21: "age21"})
    if "age60" not in wide.columns and 60 in wide.columns:  wide = wide.rename(columns={60: "age60"})
    if "age21" not in wide.columns and "21" in wide.columns: wide = wide.rename(columns={"21": "age21"})
    if "age60" not in wide.columns and "60" in wide.columns: wide = wide.rename(columns={"60": "age60"})

    # guardrails
    if not {"age21", "age60"}.issubset(wide.columns):
        print("Available columns after rename:", list(wide.columns))
        raise ValueError("Could not find both age21 and age60 columns after pivot/rename.")

    wide = wide.dropna(subset=["age21", "age60"])

    # plot
    import seaborn as sns
    import matplotlib.pyplot as plt
    plt.figure(figsize=(7, 7))
    ax = sns.scatterplot(data=wide, x="age21", y="age60", s=60, edgecolor="k", linewidth=.3)

    lim = [min(wide.min()) - 0.5, max(wide.max()) + 0.5]
    plt.plot(lim, lim, ls="--", c="k", lw=1)  # diagonal
    plt.axhline(0, ls=":", c="gray"); plt.axvline(0, ls=":", c="gray")
    plt.xlim(lim); plt.ylim(lim)
    plt.title(f"{cell_type} — mtDSB effect: 21w vs 60w")
    plt.xlabel("mtDSB effect @ 21w (log)"); plt.ylabel("mtDSB effect @ 60w (log)")

    # label top genes by |Δ|
    wide["delta"] = wide["age60"] - wide["age21"]
    top_idx = wide["delta"].abs().sort_values(ascending=False).head(label_top).index
    for g in top_idx:
        x, y = wide.loc[g, "age21"], wide.loc[g, "age60"]
        plt.text(x, y, g, fontsize=9, ha="left", va="bottom")

    plt.tight_layout()
    if save:
        plt.savefig(save, dpi=200, bbox_inches="tight")
    plt.show()

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

def plot_lesion_dynamics_by_model(comp_frac, MODELS, palette):
    n_models = len(MODELS)
    fig, axes = plt.subplots(1, n_models, figsize=(6*n_models, 5), sharey=True)
    if n_models == 1:
        axes = [axes]

    for ax, (model, cfg) in zip(axes, MODELS.items()):
        # order timepoints for this model
        time_order = [cfg["baseline"]] + cfg["courses"]
        idx = [t for t in time_order if t in comp_frac.index]
        if not idx:
            ax.set_title(f"{model} (no matching courses)")
            continue

        df = comp_frac.loc[idx, lesion_cols]  # <-- restrict to lesions only

        # plot each lesion as line
        for comp in df.columns:
            ax.plot(df.index, df[comp], marker="o",
                    label=comp, color=palette.get(comp, "#BBBBBB"),
                    lw=2)

        ax.set_title(f"{model}")
        ax.set_xlabel("Course")
        ax.set_ylabel("Fraction of cells")
        ax.set_xticks(range(len(idx)))
        ax.set_xticklabels(idx, rotation=45, ha="right")

    # single legend outside
    handles, labels = axes[-1].get_legend_handles_labels()
    fig.legend(handles, labels, title="Lesion", bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.tight_layout()
    plt.show()

def assess_mapping(
        truth,
        mapping,
        taxonomy_level,
        fontsize=15):
    """
    Generate plots assessing mapping quality

    Parameters
    ----------
    truth:
        a dict mapping cell_label to the correct assignment at each level of the taxonomy
    mapping:
        the 'results' from the mapping run
    taxonomy_level:
        the level at which we are mapping
    """
    truth_arr = np.array(
        [truth[cell['cell_id']][taxonomy_level] for cell in mapping]
    )
    mapping_arr = np.array(
        [cell[taxonomy_level]['assignment'] for cell in mapping]
    )
    prob = np.array(
        [cell[taxonomy_level]['aggregate_probability'] for cell in mapping]
    )
    corr = np.array(
        [cell[taxonomy_level]['avg_correlation'] for cell in mapping]
    )

    is_true = (mapping_arr == truth_arr)

    corr_grid = np.linspace(0.0, 1.0, 20)
    prob_grid = np.linspace(0.0, 1.0, 20)

    n_rows = 1
    n_cols = 2
    fig = plt.figure(figsize=(n_cols*5, n_rows*5))
    axis_list = [fig.add_subplot(n_rows, n_cols, ii+1) for ii in range(n_rows*n_cols)]

    axis_list[0].set_xlabel('metric cut', fontsize=fontsize)
    axis_list[0].set_ylabel('F1', fontsize=fontsize)
    axis_list[1].set_xlabel('False positive rate', fontsize=fontsize)
    axis_list[1].set_ylabel('True positive rate', fontsize=fontsize)

    for metric_grid, metric, metric_name in [(corr_grid, corr, 'metric = Correlation'),
                                             (prob_grid, prob, 'metric = Probability')]:

        tp = np.zeros(metric_grid.shape, dtype=int)
        fp = np.zeros(metric_grid.shape, dtype=int)
        fn = np.zeros(metric_grid.shape, dtype=int)
        for ii, metric_value in enumerate(metric_grid):
            considered_true = (metric >= metric_value)
            tp[ii] = np.logical_and(
                is_true,
                considered_true).sum()
            fp[ii] = np.logical_and(
                np.logical_not(is_true),
                considered_true).sum()
            fn[ii] = np.logical_and(
                is_true,
                np.logical_not(considered_true)).sum()
            fn[ii] += np.logical_not(is_true).sum()

        f1 = tp/(tp+0.5*(fp+fn))
        axis_list[0].plot(metric_grid, f1, label=metric_name)
        axis_list[1].plot(fp/len(mapping_arr), tp/len(mapping_arr), label=metric_name)

    axis_list[0].set_title(taxonomy_level, fontsize=fontsize)
    for axis in axis_list:
        axis.legend(loc=0, fontsize=fontsize)
    fig.tight_layout()

def plot_confusion_matrix(
        row_assignments,
        col_assignments,
        row_axis_label,
        col_axis_label,
        label_elements=True,
        fontsize=15,
        colorbar_title='Jaccard index',
        col_cut=None,
        row_cut=None,
        order_wgt=0.1,
        ordering_iterations=300000):
    """
    Create a confusion matrix plot colored by the Jaccard index
    to show the correspondence between cell types in two taxonomies.

    Parameters
    ----------
    row_assignments:
        a list of str. One for each cell. The assignment of the cells to the
        cell types listed along the rows of the confusion matrix
    col_assignments:
        as above. The assignment of the cells to the cell types
        listed along the columns of the confusion matrix
    row_axis_label:
        label for the row axis (the vertical axis) of the
        confusion matrix
    col_axis_label:
        label for the column axis (the horizontal axis)
        of the confusion matrix
    label_elements:
        a boolean. If True, label the rows and columns of the confusion
        matrix. If False, do not (probably because there are too many of them)
    color_bar_title:
        title to list along the color bar of the confusion matrix
    col_cut:
        an optional int. If non-None, only show the first {col_cut}
        columns of the confusion matrix
    row_cut:
        an optional int. If non-None, only show the first {row_cut}
        rows of the founction matrix
    order_wgt:
        a float. Prioritize how important it is that cell types be in
        sorted order according to how many cells are assigned to those
        types when ordering rows and columns of the cofnusion matrix.
    ordering_iterations:
        how many random steps to take when trying to find the most
        aesthetically pleasing order for rows and columns of the
        confusion matrix
    """
    confusion = create_confusion_matrix(
        row_labels=row_assignments,
        col_labels=col_assignments
    )
    confusion_matrix = confusion['array']
    row_labels = confusion['rows']
    col_labels = confusion['cols']
    swapping = get_confusion_order(
        input_arr=confusion_matrix,
        ct_arr=confusion['ct'],
        rng=np.random.default_rng(8711121),
        n_iter=ordering_iterations,
        inverse_temp=10.0,
        order_wgt=order_wgt)

    row_labels = np.array(row_labels)[swapping['rows']]
    col_labels = np.array(col_labels)[swapping['cols']]

    fig = plt.figure(figsize=(20, 20))
    to_plot = np.copy(confusion_matrix)
    to_plot = to_plot[swapping['rows'], :]
    to_plot = to_plot[:, swapping['cols']]

    if col_cut is not None:
        col_labels = col_labels[:col_cut]
        to_plot = to_plot[:, :col_cut]
    if row_cut is not None:
        row_labels = row_labels[:row_cut]
        to_plot = to_plot[:row_cut, :]

    print(f'to_plot shape {to_plot.shape}')
    heatmap_axis = fig.add_subplot(2,1,1)
    heatmap_axis.set_title(f'{row_axis_label} vs {col_axis_label}', fontsize=fontsize)
    img = heatmap_axis.imshow(np.ma.masked_array(to_plot, to_plot==0.0), cmap=matplotlib.cm.Reds)

    divider = make_axes_locatable(heatmap_axis)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(
        img,
        ax=heatmap_axis,
        cax=cax,
        label=colorbar_title)

    cbar.ax.tick_params(axis='both', which='both', labelsize=fontsize)
    cbar.ax.set_ylabel(ylabel=colorbar_title, fontdict={'fontsize': fontsize})

    heatmap_axis.set_xlabel(col_axis_label, fontsize=fontsize)
    heatmap_axis.set_ylabel(row_axis_label, fontsize=fontsize)

    if label_elements:
        row_ticks = [ii+0.5 for ii in range(len(row_labels))]
        heatmap_axis.set_yticks(row_ticks, labels=row_labels, va='bottom', ha='right')

        col_ticks = [ii-0.5 for ii in range(len(col_labels))]
        heatmap_axis.set_xticks(col_ticks, labels=col_labels, va='top', ha='left')

        heatmap_axis.tick_params(
            which='both',
            axis='both',
            labelsize=fontsize,
            pad=0,
            size=10
        )
        heatmap_axis.tick_params(
            axis='x',
            labelrotation=-90.0)

    row_hist_axis = fig.add_subplot(2, 2, 3)
    col_hist_axis = fig.add_subplot(2, 2, 4)
    for axis, assignments, labels in [(row_hist_axis, row_assignments, row_labels),
                                      (col_hist_axis, col_assignments, col_labels)]:
        l_to_idx = {l: ii+0.5 for ii, l in enumerate(labels)}
        assn_idx = [l_to_idx[a] for a in assignments if a in l_to_idx]
        _ = axis.hist(assn_idx, bins=np.arange(len(labels)+1))
        if label_elements:
            ticks = [ii for ii in range(len(labels))]
            axis.set_xticks(ticks, labels, va='top', ha='left')
            axis.tick_params(
                which='both',
                axis='both',
                labelsize=fontsize
            )
            axis.tick_params(
                axis='x',
                labelrotation=-90.0
            )
    row_hist_axis.set_xlabel(row_axis_label, fontsize=fontsize)
    col_hist_axis.set_xlabel(col_axis_label, fontsize=fontsize)
    row_hist_axis.set_ylabel('N cells', fontsize=fontsize)
    col_ylim = col_hist_axis.get_ylim()
    row_ylim = row_hist_axis.get_ylim()
    if col_ylim[1] > row_ylim[1]:
        row_hist_axis.set_ylim(col_ylim)
    else:
        col_hist_axis.set_ylim(row_ylim)
    fig.tight_layout()

def plot_umap(xx, yy, cc=None, val=None, fig_width=8, fig_height=8, cmap=None):

    fig, ax = plt.subplots()
    fig.set_size_inches(fig_width, fig_height)

    if cmap is not None :
        plt.scatter(xx, yy, s=0.5, c=val, marker='.', cmap=cmap)
    elif cc is not None :
        plt.scatter(xx, yy, s=0.5, color=cc, marker='.')

    ax.axis('equal')
    ax.set_xlim(-18, 27)
    ax.set_ylim(-18, 27)
    ax.set_xticks([])
    ax.set_yticks([])

    return fig, ax
