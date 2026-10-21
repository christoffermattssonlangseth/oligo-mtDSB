#!/usr/bin/env python
"""Myelin and OL-lineage transcripts binned on a grid, inside and outside segmented cells, for analysis/13.

About three quarters of Mbp transcripts fall outside every cell polygon (myelin sheaths and processes), so cell-based
myelin scores miss most of the signal. This reads each section's transcripts.parquet (qv >= 20) for ``GENES`` and
writes one table to ``--out``:

* ``myelin_tx_bins.csv.gz``  sample_id, bx, by (``--bin`` um grid), gene, n (all), n_unassigned (no cell),
                             n_nucleus (overlaps a nucleus)

Regions are attached later from the nearest cell centroid (results/cell_table/cells.csv.gz), so this needs only
the raw Xenium folders. Read with polars (pyarrow <= 22 fails on these files); office Mac env ``sopa``, ~5 min.

Example (office Mac)
--------------------
    python scripts/export_myelin_transcripts.py --raw-dir /Volumes/moldiassd/oligo-mtDSB/raw_xenium/run4 \
        --out /Volumes/moldiassd/oligo-mtDSB/myelin_tx
"""

from __future__ import annotations

import argparse
from pathlib import Path

import polars as pl

GENES = ["Mbp", "Mog", "Mag", "Cnp", "Opalin", "Cldn11", "Ermn", "Sox10", "Olig2", "Trib3", "Cdkn1a", "Gfap", "Snap25"]


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--raw-dir", required=True, help="folder holding the output-XETG... run folders")
    p.add_argument("--out", required=True)
    p.add_argument("--bin", type=float, default=20.0)
    a = p.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    parts = []
    for run in sorted(Path(a.raw_dir).glob("output-*")):
        sample = run.name.split("__")[2]
        tx = (pl.scan_parquet(run / "transcripts.parquet")
              .filter(pl.col("feature_name").is_in(GENES) & (pl.col("qv") >= 20))
              .select(pl.lit(sample).alias("sample_id"),
                      (pl.col("x_location") // a.bin).cast(pl.Int32).alias("bx"),
                      (pl.col("y_location") // a.bin).cast(pl.Int32).alias("by"),
                      pl.col("feature_name").cast(pl.Utf8).alias("gene"),
                      (pl.col("cell_id") == "UNASSIGNED").alias("unassigned"),
                      pl.col("overlaps_nucleus").cast(pl.Boolean).alias("nucleus"))
              .group_by(["sample_id", "bx", "by", "gene"])
              .agg(pl.len().alias("n"), pl.col("unassigned").sum().alias("n_unassigned"), pl.col("nucleus").sum().alias("n_nucleus"))
              .collect())
        print(sample, tx["n"].sum(), "transcripts", flush=True)
        parts.append(tx)
    pl.concat(parts).sort(["sample_id", "gene", "bx", "by"]).to_pandas().to_csv(out / "myelin_tx_bins.csv.gz", index=False)
    print("wrote", out / "myelin_tx_bins.csv.gz")


if __name__ == "__main__":
    main()
