#!/usr/bin/env python
"""Export a light per-cell table (metadata + spatial centroids) from the annotated h5ad.

The figure notebooks need cell positions and section metadata but not the
10 GB expression matrix, so this writes only obs columns:

    cell_id, sample_id, age, condition, sex, genotype, run, segmentation_method, cell_class,
    RBD_compartment_simplified, x_centroid, y_centroid

plus ``sections.csv`` with one row per section (metadata, Xenium slide, cell count, cells per segmentation method).

Example
-------
    python scripts/export_cell_table.py \
        --h5ad /Volumes/processing2/oligo-mtDSB/data/mtDNA_DSB_5k_clustered_annotation_with_rbd_2_cytetype_brain_novae4.h5ad \
        --outdir results/cell_table
"""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import pandas as pd
from anndata.io import read_elem

COLS = ["sample_id", "age", "condition", "sex", "genotype", "run", "segmentation_method", "cell_class", "RBD_compartment_simplified",
        "x_centroid", "y_centroid"]


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--h5ad", required=True)
    p.add_argument("--outdir", required=True)
    a = p.parse_args(argv)
    with h5py.File(a.h5ad, "r") as f:
        obs = f["obs"]
        df = pd.DataFrame({c: read_elem(obs[c]) for c in COLS}, index=pd.Index(read_elem(obs[obs.attrs["_index"]]), name="cell_id"))
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / "cells.csv.gz", compression="gzip")
    sections = df.groupby("sample_id", observed=True).agg(
        age=("age", "first"), condition=("condition", "first"), sex=("sex", "first"),
        genotype=("genotype", "first"), run=("run", "first"), n_cells=("age", "size"))
    sections["xenium_slide"] = sections["run"].str.extract(r"__(\d{7})__", expand=False)
    seg = df.groupby(["sample_id", "segmentation_method"], observed=True).size().unstack(fill_value=0)
    seg.columns = [f"n_seg[{c.split(' (')[0]}]" for c in seg.columns]
    sections = sections.join(seg)
    sections.to_csv(out / "sections.csv")
    print(sections.to_string())


if __name__ == "__main__":
    main()
