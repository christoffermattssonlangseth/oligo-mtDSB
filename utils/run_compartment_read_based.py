#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Run read-based domain discovery with checkpointing and robust cell_id handling.

Example:
python run_compartment_read_based.py \
  --data-base /date/gcb/gcb_CML/oligo-mtDSB \
  --adata data/mtDNA_DSB_5k_clustered_manual_annotation.h5ad \
  --results-subdir results/rbd_runs/TEST \
  --checkpoints-subdir results/rbd_checkpoints_mtDSB \
  --sample-key sample_id \
  --cell-id-key cell_id \
  --threads 8 \
  --force-rerun 0 \
  --tag monod
"""

import os
import time
import json
import argparse
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
import scipy.sparse as sp
from anndata import AnnData


# ------------------------------ helpers ------------------------------ #

def _ensure_spatial_in_obsm(adata: AnnData, spatial_key: str = "spatial"):
    """Ensure adata.obsm[spatial_key] exists, falling back to obs columns."""
    if spatial_key in adata.obsm:
        return
    for xy in (("X", "Y"), ("x", "y"), ("px_x", "px_y"), ("imagecol", "imagerow")):
        if set(xy).issubset(adata.obs.columns):
            coords = adata.obs[list(xy)].to_numpy(dtype=float)
            adata.obsm[spatial_key] = coords
            return
    raise KeyError(
        f"Could not find {spatial_key} in .obsm and no fallback XY columns present "
        f"in .obs. Add coordinates or rename to one of "
        f"[('X','Y'), ('x','y'), ('px_x','px_y'), ('imagecol','imagerow')]."
    )


def _normalize_obs_and_cell_id(ad: AnnData, cell_id_key: str) -> None:
    """
    Ensure:
      - obs_names are unique
      - obs[cell_id_key] == obs_names (as str)
      - obs.index.name != cell_id_key (clear if so)
    """
    ad.obs_names_make_unique()
    ad.obs[cell_id_key] = ad.obs_names.astype(str)
    if ad.obs.index.name == cell_id_key:
        ad.obs.index.name = None


def _load_h5ad_with_heal(path: str, cell_id_key: str) -> AnnData:
    """Load an .h5ad, normalize cell IDs, and persist the fix to disk."""
    ad = sc.read_h5ad(path)
    _normalize_obs_and_cell_id(ad, cell_id_key)
    ad.write_h5ad(path)
    return ad


def format_data_neighs_colapse(
    adata: AnnData,
    spatial_key: str = "spatial",
    neighs: int = 10
) -> AnnData:
    """Pseudobin each cell’s expression by summing over its spatial neighbors."""
    _ensure_spatial_in_obsm(adata, spatial_key=spatial_key)

    ad = adata.copy()
    print(f"    → Computing spatial neighbors (k={neighs})…", end="", flush=True)
    sq.gr.spatial_neighbors(ad, spatial_key=spatial_key, n_neighs=neighs)
    print("done")

    # adjacency (include self)
    C = ad.obsp["spatial_connectivities"].astype(bool).astype(int)
    C = C + sp.eye(C.shape[0], format="csr")

    X = ad.X
    if isinstance(X, np.ndarray):
        X = sp.csr_matrix(X)

    Xc = C.dot(X)

    print("    → Pseudobinning complete")
    out = AnnData(
        X=Xc,
        obs=ad.obs.copy(),
        var=ad.var.copy(),
        obsm=ad.obsm.copy(),
        uns=ad.uns.copy()
    )
    return out


def domains_by_rbd(
    adata: AnnData,
    hyperparams: dict,
    sample_key: str = "sample",
    cell_id_key: str = "cell_id",
    checkpoint_dir: str = "./checkpoints",
    force_rerun: bool = False
):
    """
    Read-based domains with checkpointed saves, force-rerun option, and printouts.
    Returns: (adata_with_domains, adata_pseudobinned_processed)
    """
    os.makedirs(checkpoint_dir, exist_ok=True)

    # Always normalize input (paranoid/safe)
    _normalize_obs_and_cell_id(adata, cell_id_key)

    def load_or_compute(path, compute_fn, step_name):
        if (not force_rerun) and os.path.exists(path):
            print(f"{step_name}: Loading from checkpoint…", end="", flush=True)
            ad = _load_h5ad_with_heal(path, cell_id_key)
            print("loaded\n")
            return ad, False
        else:
            print(f"{step_name}: Computing…", flush=True)
            ad = compute_fn()
            _normalize_obs_and_cell_id(ad, cell_id_key)
            ad.write_h5ad(path)
            print(f"{step_name}: Saved to {path}\n")
            return ad, True

    # 1) PSEUDOBINNING
    pseudobin_path = os.path.join(checkpoint_dir, "pseudobinned.h5ad")

    def compute_pseudobin():
        parts = []
        for samp in adata.obs[sample_key].unique():
            print(f"  Sample '{samp}':")
            sub = adata[adata.obs[sample_key] == samp].copy()
            _normalize_obs_and_cell_id(sub, cell_id_key)
            pb = format_data_neighs_colapse(
                sub, spatial_key="spatial", neighs=hyperparams["neighbors"]
            )
            _normalize_obs_and_cell_id(pb, cell_id_key)
            parts.append(pb)
        out = sc.concat(parts, merge="same")
        _normalize_obs_and_cell_id(out, cell_id_key)
        return out

    pseudobinned, _ = load_or_compute(
        pseudobin_path, compute_pseudobin, "Step 1) Pseudobinning"
    )

    # 2) FILTERING
    filtered_path = os.path.join(checkpoint_dir, "filtered.h5ad")

    def compute_filtered():
        ad_neigh = pseudobinned.copy()
        ad_neigh.X = np.nan_to_num(ad_neigh.X)
        sc.pp.filter_cells(ad_neigh, min_counts=4)
        ad_neigh.raw = ad_neigh
        _normalize_obs_and_cell_id(ad_neigh, cell_id_key)
        return ad_neigh

    adata_neigh, _ = load_or_compute(
        filtered_path, compute_filtered, "Step 2) Filtering"
    )

    # 3) NEIGHBORS/UMAP
    umap_path = os.path.join(checkpoint_dir, "umap.h5ad")

    def compute_umap():
        nn = hyperparams.get("n_neighbors", 20)
        pcs = hyperparams.get("n_pcs", None)
        if pcs is None:
            sc.pp.neighbors(adata_neigh, n_neighbors=nn, use_rep="X")
        else:
            sc.tl.pca(adata_neigh, n_comps=pcs, svd_solver="arpack")
            sc.pp.neighbors(adata_neigh, n_neighbors=nn, n_pcs=pcs)
        sc.tl.umap(adata_neigh, min_dist=hyperparams.get("min_dist", 0.1))
        _normalize_obs_and_cell_id(adata_neigh, cell_id_key)
        return adata_neigh

    adata_neigh, _ = load_or_compute(umap_path, compute_umap, "Step 3) UMAP")

    # 4) CLUSTERING at one or many resolutions
    algo = hyperparams["clustering_algorithm"].lower()
    res_list = hyperparams.get("resolutions", [hyperparams.get("resolution")])

    for res in res_list:
        key = f"rbd_domain_{res}"
        chk = os.path.join(checkpoint_dir, f"clustered_{res}.h5ad")

        def compute_cluster():
            if algo == "leiden":
                sc.tl.leiden(adata_neigh, resolution=res, key_added=key)
            else:
                sc.tl.louvain(adata_neigh, resolution=res, key_added=key)
            _normalize_obs_and_cell_id(adata_neigh, cell_id_key)
            return adata_neigh

        adata_neigh, _ = load_or_compute(
            chk, compute_cluster, f"Step 4) Clustering (res={res})"
        )

        # map back to original cells
        print(f"  → Mapping '{key}' back to original…", end="", flush=True)
        if cell_id_key not in adata_neigh.obs.columns:
            raise KeyError(
                f"Expected '{cell_id_key}' in adata_neigh.obs but it is missing."
            )
        mapping = dict(zip(adata_neigh.obs[cell_id_key], adata_neigh.obs[key]))
        _normalize_obs_and_cell_id(adata, cell_id_key)
        adata.obs[key] = adata.obs[cell_id_key].map(mapping).astype(str)
        print("done\n")

    print("All steps complete! Returning results.")
    return adata, adata_neigh


# ------------------------------ main ------------------------------ #

def parse_args():
    p = argparse.ArgumentParser(description="Run read-based domain discovery.")
    p.add_argument("--data-base", required=True, help="Base directory")
    p.add_argument("--adata", required=True, help="Path relative to base, e.g. data/foo.h5ad")
    p.add_argument("--results-subdir", required=True, help="Relative results dir to write outputs")
    p.add_argument("--checkpoints-subdir", required=True, help="Relative checkpoints dir")
    p.add_argument("--sample-key", default="sample_id")
    p.add_argument("--cell-id-key", default="cell_id")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--force-rerun", type=int, default=0, choices=[0, 1])
    p.add_argument("--tag", default="", help="Suffix tag for output filenames")
    # hyperparams (override if you like)
    p.add_argument("--neighbors", type=int, default=20)
    p.add_argument("--n-pcs", type=int, default=30)
    p.add_argument("--n-neighbors", type=int, default=20)
    p.add_argument("--min-dist", type=float, default=0.5)
    p.add_argument("--cluster", choices=["leiden", "louvain"], default="leiden")
    p.add_argument("--resolutions", default="0.1,0.2,0.5,1.0",
                   help="Comma-separated list, e.g. '0.1,0.2,0.5,1.0'")
    return p.parse_args()


def main():
    args = parse_args()

    # threads
    os.environ["OMP_NUM_THREADS"] = str(args.threads)
    os.environ["OPENBLAS_NUM_THREADS"] = str(args.threads)
    os.environ["MKL_NUM_THREADS"] = str(args.threads)
    sc.settings.n_jobs = args.threads

    # paths
    base = os.path.abspath(args.data_base)
    adata_path = os.path.join(base, args.adata)
    results_dir = os.path.join(base, args.results_subdir)
    ckpt_dir = os.path.join(base, args.checkpoints_subdir)
    os.makedirs(results_dir, exist_ok=True)
    os.makedirs(ckpt_dir, exist_ok=True)

    # read data (heal immediately)
    print(f"Reading AnnData: {adata_path}")
    ad = sc.read_h5ad(adata_path)
    _normalize_obs_and_cell_id(ad, args.cell_id_key)

    # hyperparameters
    res_list = [float(r.strip()) for r in args.resolutions.split(",") if r.strip()]
    hyper = {
        "neighbors": int(args.neighbors),
        "clustering_algorithm": args.cluster,
        "resolutions": res_list,
        "n_neighbors": int(args.n_neighbors),
        "n_pcs": int(args.n_pcs) if args.n_pcs else None,
        "min_dist": float(args.min_dist),
    }

    # run
    t0 = time.time()
    ad_annot, ad_pseudo = domains_by_rbd(
        adata=ad,
        hyperparams=hyper,
        sample_key=args.sample_key,
        cell_id_key=args.cell_id_key,
        checkpoint_dir=str(ckpt_dir),
        force_rerun=bool(args.force_rerun),
    )

    # save (normalize again to avoid any residual mismatch)
    tag = f"_{args.tag}" if args.tag else ""
    out_annot = os.path.join(results_dir, f"rbd_annotated{tag}.h5ad")
    out_pseudo = os.path.join(results_dir, f"rbd_pseudobinned{tag}.h5ad")
    print(f"💾 Writing: {out_annot}")
    _normalize_obs_and_cell_id(ad_annot, args.cell_id_key)
    ad_annot.write_h5ad(out_annot, compression="gzip")
    print(f"💾 Writing: {out_pseudo}")
    _normalize_obs_and_cell_id(ad_pseudo, args.cell_id_key)
    ad_pseudo.write_h5ad(out_pseudo, compression="gzip")

    # metadata
    meta = {
        "host": os.uname().nodename,
        "start_time": time.strftime("%F %T"),
        "elapsed_min": round((time.time() - t0) / 60.0, 2),
        "adata_path": adata_path,
        "results_dir": results_dir,
        "checkpoint_dir": ckpt_dir,
        "sample_key": args.sample_key,
        "cell_id_key": args.cell_id_key,
        "threads": args.threads,
        "force_rerun": bool(args.force_rerun),
        "hyperparams": hyper,
        "n_cells_in": int(ad.n_obs),
        "n_vars_in": int(ad.n_vars),
        "n_cells_annot": int(ad_annot.n_obs),
        "n_vars_annot": int(ad_annot.n_vars),
        "n_cells_pseudobinned": int(ad_pseudo.n_obs),
        "n_vars_pseudobinned": int(ad_pseudo.n_vars),
    }
    meta_path = os.path.join(results_dir, f"rbd_run_meta{tag}.json")
    with open(meta_path, "w") as f:
        json.dump(meta, f, indent=2)
    print(f"✅ Done in {meta['elapsed_min']} min | meta: {meta_path}")


if __name__ == "__main__":
    main()