#!/usr/bin/env python
"""Run the treatment x region x time-point DE pipeline on the full data set.

The annotated h5ad lives on the processing server, so this script is meant to
be run where the data is mounted (laptop with /Volumes/processing2, or the
cluster). It writes tidy CSV tables to ``--outdir`` that the notebook
``notebooks-05/01-region-time-treatment-figures.ipynb`` turns into figures.

Example
-------
    conda activate karospace311     # needs scanpy + pydeseq2 >= 0.5
    python scripts/run_region_time_de.py \
        --h5ad /Volumes/processing2/oligo-mtDSB/data/mtDNA_DSB_5k_clustered_annotation_with_rbd_2_cytetype_brain_novae3.h5ad \
        --outdir results/region_time \
        --layer counts --n-cpus 8

Use ``--cell-types "Oligodendrocytes" "Astrocytes"`` to restrict, and
``--no-region-vs-rest`` / ``--no-interaction`` to skip the slower tests.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from utils.region_time_de import Keys, read_minimal, run_all  # noqa: E402


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--h5ad", required=True, help="annotated AnnData with raw counts")
    p.add_argument("--outdir", default=str(REPO / "results" / "region_time"))
    p.add_argument("--layer", default="counts", help="layer with raw integer counts; use 'X' for adata.X")
    p.add_argument("--cell-type-key", default="cell_class")
    p.add_argument("--region-key", default="RBD_compartment_simplified")
    p.add_argument("--sample-key", default="sample_id")
    p.add_argument("--age-key", default="age")
    p.add_argument("--condition-key", default="condition")
    p.add_argument("--case", default="mtDSB")
    p.add_argument("--control", default="control")
    p.add_argument("--ages", nargs=2, default=["21", "60"], metavar=("EARLY", "LATE"))
    p.add_argument("--cell-types", nargs="*", default=None)
    p.add_argument("--regions", nargs="*", default=None)
    p.add_argument("--exclude-regions", nargs="*", default=["Unknown"])
    p.add_argument("--min-cells-per-sample", type=int, default=20)
    p.add_argument("--min-per-condition", type=int, default=2)
    p.add_argument("--no-region-vs-rest", action="store_true")
    p.add_argument("--no-interaction", action="store_true")
    p.add_argument("--n-cpus", type=int, default=None)
    p.add_argument("--padj", type=float, default=0.05)
    p.add_argument("--lfc", type=float, default=0.5)
    p.add_argument("--min-base-mean", type=float, default=20.0)
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    layer = None if args.layer == "X" else args.layer
    keys = Keys(
        cell_type=args.cell_type_key,
        region=args.region_key,
        sample=args.sample_key,
        age=args.age_key,
        condition=args.condition_key,
        case=args.case,
        control=args.control,
        ages=tuple(args.ages),
    )
    # only obs columns + the counts matrix are read; X and other layers stay on disk
    adata = read_minimal(args.h5ad, layer, keys)
    logging.info("loaded %d cells x %d genes from %s", adata.n_obs, adata.n_vars, args.h5ad)
    tables = run_all(
        adata,
        args.outdir,
        keys=keys,
        layer=layer,
        min_cells_per_sample=args.min_cells_per_sample,
        exclude_regions=args.exclude_regions,
        cell_types=args.cell_types,
        regions=args.regions,
        min_per_condition=args.min_per_condition,
        region_vs_rest=not args.no_region_vs_rest,
        interaction=not args.no_interaction,
        n_cpus=args.n_cpus,
        padj_thr=args.padj,
        lfc_thr=args.lfc,
        min_base_mean=args.min_base_mean,
    )
    summ = tables["summary"]
    print(f"\nwrote {args.outdir}")
    print(summ.groupby(["cell_type", "age"])["n_sig"].sum().unstack("age").to_string())


if __name__ == "__main__":
    main()
