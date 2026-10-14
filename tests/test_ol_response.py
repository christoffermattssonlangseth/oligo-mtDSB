"""Tests for utils.ol_response on a synthetic pseudobulk with planted responders.

Run with pytest, or directly:  python tests/test_ol_response.py

Design: 12 sections (2 ages x 2 conditions x 3). Two of the six mtDSB sections
("responders") get 4x ISR gene expression in OLs and half their OLs relabelled
DA-Oligodendrocytes; the other mtDSB sections look like controls.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from utils import ol_response as olr  # noqa: E402
from utils.region_time_de import TISSUE_WIDE, Pseudobulk  # noqa: E402

RESPONDERS = {"RB21m0", "RB60m1"}
REGIONS = ["Thalamus", "Cortex"]


def make_pb(seed: int = 0) -> tuple[Pseudobulk, pd.DataFrame]:
    rng = np.random.default_rng(seed)
    genes = olr.ISR_GENES + [f"g{i}" for i in range(40)]
    rows, meta = [], []
    sections = []
    for age in ("21", "60"):
        for cond in ("control", "mtDSB"):
            for r in range(3):
                sid = f"RB{age}{cond[0]}{r}"
                sections.append((sid, age, cond))
                for reg in REGIONS + [TISSUE_WIDE]:
                    for ct in ("Oligodendrocytes", "DA-Oligodendrocytes", "Astrocytes"):
                        n = 200 if ct != "DA-Oligodendrocytes" else 10
                        if ct == "DA-Oligodendrocytes" and sid in RESPONDERS:
                            n = 200
                        mu = np.full(len(genes), 50.0) * n / 10
                        if ct != "Astrocytes" and sid in RESPONDERS:
                            mu[: len(olr.ISR_GENES)] *= 4
                        rows.append(rng.poisson(mu))
                        meta.append({"cell_type": ct, "region": reg, "sample_id": sid, "age": age, "condition": cond, "n_cells": n})
    meta = pd.DataFrame(meta)
    meta.index = [f"{m.cell_type}|{m.region}|{m.sample_id}" for m in meta.itertuples()]
    counts = pd.DataFrame(np.array(rows), index=meta.index, columns=genes)
    small = meta["n_cells"] < 20
    pb = Pseudobulk(counts.loc[~small], meta.loc[~small], meta.loc[small])
    sec = pd.DataFrame(sections, columns=["sample_id", "age", "condition"]).set_index("sample_id")
    return pb, sec


def test_isr_scores_and_responder_calls():
    pb, sec = make_pb()
    isr = olr.isr_scores(olr.general_pseudobulk(pb))
    assert set(isr.region) == set(REGIONS) | {TISSUE_WIDE}
    calls = olr.responder_calls(isr, sec)
    assert set(calls.index[calls.group == "responder"]) == RESPONDERS
    assert set(calls.index[calls.group == "non-responder"]) == {s for s in sec.index if s.startswith("RB") and "m" in s} - RESPONDERS
    assert (calls.loc[calls.condition == "control", "group"] == "control").all()


def test_cluster_share_uses_dropped_rows():
    pb, sec = make_pb()
    share = olr.cluster_share(pb)
    tw = share[share.region == TISSUE_WIDE].set_index("sample_id")["share"]
    assert np.allclose(tw[list(RESPONDERS)], 0.5)
    non = [s for s in sec.index if s not in RESPONDERS]
    assert np.allclose(tw[non], 10 / 210)  # DA rows had 10 cells -> dropped, still counted


def test_high_score_fraction_and_region_tables():
    pb, sec = make_pb()
    rng = np.random.default_rng(1)
    rows = []
    for sid in sec.index:
        for reg in REGIONS:
            p = rng.uniform(0.85, 1.0, 100) if sid in RESPONDERS and reg == "Thalamus" else rng.uniform(0, 0.6, 100)
            rows.append(pd.DataFrame({"cell_type": "Oligodendrocytes", "region": reg, "sample_id": sid, "p_pos[condition]": p}))
    cs = pd.concat(rows, ignore_index=True)
    hf = olr.high_score_fraction(cs)
    tw = hf[hf.region == TISSUE_WIDE].set_index("sample_id")["frac_high"]
    assert np.allclose(tw[list(RESPONDERS)], 0.5) and (tw.drop(list(RESPONDERS)) == 0).all()
    calls = olr.responder_calls(olr.isr_scores(olr.general_pseudobulk(pb)), sec)
    rk = olr.region_rank_per_animal(hf, calls, "frac_high")
    assert rk.index[0] == "Thalamus"
    t = olr.region_by_group(hf, calls, "frac_high")
    assert t.loc["Thalamus", "responder"] > 0.9 and t.loc["Thalamus", "control"] == 0


def test_gene_score_association_within_section():
    rng = np.random.default_rng(2)
    rows, samples, scores = [], [], []
    for sid in ["A", "B"]:
        n = 400
        sc_ = rng.uniform(0.01, 0.99, n)
        x = rng.normal(size=(n, 5))
        x[:, 0] += 2 * np.log(sc_ / (1 - sc_))       # tracks the score within the section
        x[:, 1] += 5 if sid == "B" else 0            # section offset only, no within-section signal
        rows.append(x); samples += [sid] * n; scores.append(sc_)
    assoc = olr.gene_score_association(np.vstack(rows), np.array(samples), np.concatenate(scores),
                                       [f"g{i}" for i in range(5)], min_cells=100)
    r = assoc.groupby("gene")["r"].mean()
    assert r["g0"] > 0.8 and abs(r["g1"]) < 0.15


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print(f"{name}: ok")
