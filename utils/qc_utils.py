"""
Auto-generated utilities for qc_utils.
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

def filter_metabolic(de_df, gene_lists):
    metabolic_genes = set().union(*gene_lists.values())
    return de_df[de_df["names"].isin(metabolic_genes)].copy()
