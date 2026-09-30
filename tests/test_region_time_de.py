"""End-to-end test of utils.region_time_de on a small synthetic Xenium-like data set.

Run with pytest, or directly:  python tests/test_region_time_de.py

The synthetic design mirrors the real study: 12 sections = 2 ages x 2
conditions x 3 sections, 3 regions, 2 cell types, 250 genes. Effects planted
in cell type "OL":
  * genes g0-g4  : up in mtDSB in every region at both ages (global effect)
  * genes g5-g9  : up in mtDSB only at the late age, all regions (age x condition)
  * genes g10-g14: up in mtDSB only in region "Thalamus", both ages (region-specific)
Cell type "Astro" has no planted effects.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from utils import region_time_de as rt  # noqa: E402

GLOBAL = [f"g{i}" for i in range(0, 5)]
LATE_ONLY = [f"g{i}" for i in range(5, 10)]
THAL_ONLY = [f"g{i}" for i in range(10, 15)]
REGIONS = ["Cortex", "Thalamus", "Hypothalamus"]


def make_adata(seed: int = 0, cells_per_group: int = 60, n_genes: int = 250):
    import anndata as ad
    from scipy import sparse

    rng = np.random.default_rng(seed)
    samples = []
    for age in ("21", "60"):
        for cond in ("control", "mtDSB"):
            for r in range(3):
                samples.append((f"RB{age}{cond[0]}{r}", age, cond))
    genes = [f"g{i}" for i in range(n_genes)]
    base = rng.gamma(1.5, 1.2, size=n_genes)  # per-cell mean counts
    base[:15] = np.maximum(base[:15], 0.8)  # planted genes must be expressed enough to be testable
    rows, obs = [], []
    for sid, age, cond in samples:
        for ct in ("OL", "Astro"):
            for reg in REGIONS + ["Unknown"]:
                n = cells_per_group if reg != "Unknown" else 5
                mu = base.copy()
                if ct == "OL" and cond == "mtDSB":
                    mu[:5] *= 4
                    if age == "60":
                        mu[5:10] *= 4
                    if reg == "Thalamus":
                        mu[10:15] *= 4
                x = rng.poisson(mu[None, :] * rng.lognormal(0, 0.3, size=(n, 1)))
                rows.append(sparse.csr_matrix(x))
                obs.append(pd.DataFrame({"sample_id": sid, "age": age, "condition": cond, "cell_class": ct,
                                         "RBD_compartment_simplified": reg}, index=range(n)))
    X = sparse.vstack(rows).tocsr()
    obs = pd.concat(obs, ignore_index=True)
    obs.index = [f"c{i}" for i in range(len(obs))]
    a = ad.AnnData(X=X.astype(np.float32), obs=obs)
    a.var_names = genes
    a.layers["counts"] = X.copy()
    return a


def test_pseudobulk_shapes():
    a = make_adata()
    pb = rt.pseudobulk(a, layer="counts", min_cells_per_sample=20)
    # 2 cell types x (3 regions + tissue-wide) x 12 sections
    assert pb.counts.shape == (2 * 4 * 12, a.n_vars)
    assert set(pb.meta["region"]) == set(REGIONS) | {rt.TISSUE_WIDE}
    assert "Unknown" not in set(pb.meta["region"])
    assert pb.dropped.empty
    assert (pb.counts.loc[pb.meta["region"] == rt.TISSUE_WIDE].sum().sum()
            == pb.counts.loc[pb.meta["region"] != rt.TISSUE_WIDE].sum().sum())
    # rows are integer counts
    assert np.issubdtype(pb.counts.dtypes.iloc[0], np.integer)


def test_min_cells_dropping():
    a = make_adata(cells_per_group=15)
    pb = rt.pseudobulk(a, layer="counts", min_cells_per_sample=20, include_tissue_wide=False)
    assert pb.counts.empty and len(pb.dropped) == 2 * 3 * 12


def test_pipeline_recovers_planted_effects(tmp_path=None):
    tmp_path = Path(tmp_path) if tmp_path else REPO / "tests" / "_tmp_region_time"
    a = make_adata()
    tables = rt.run_all(a, tmp_path, layer="counts", n_cpus=2)

    se = tables["simple_effects"]
    log = tables["simple_effects_log"]
    assert (log["status"] == "ok").all(), log
    assert set(se.columns) >= set(["cell_type", "region", "age", *rt.STANDARD_COLS])

    def hits(ct, reg, age):
        s = se[(se.cell_type == ct) & (se.region == reg) & (se.age == age)]
        return set(s.loc[rt.significant(s, 0.05, 0.5, 5.0), "gene"])

    # global genes everywhere in OL
    for reg in REGIONS:
        for age in ("21", "60"):
            assert set(GLOBAL) <= hits("OL", reg, age), (reg, age)
    # late-only genes: present at 60, absent at 21
    assert set(LATE_ONLY) <= hits("OL", "Cortex", "60")
    assert not (set(LATE_ONLY) & hits("OL", "Cortex", "21"))
    # thalamus-only genes
    assert set(THAL_ONLY) <= hits("OL", "Thalamus", "21")
    assert not (set(THAL_ONLY) & hits("OL", "Cortex", "21"))
    # no planted effect in astrocytes: at most a couple of false positives per stratum
    for reg in REGIONS:
        assert len(hits("Astro", reg, "60")) <= 3

    # age contrast (delta z-test): late-only genes significant, global genes not
    ac = tables["age_contrast"]
    thal = ac[(ac.cell_type == "OL") & (ac.region == "Thalamus")].set_index("gene")
    assert (thal.loc[LATE_ONLY, "padj"] < 0.05).all()
    assert (thal.loc[GLOBAL, "padj"] > 0.05).all()
    assert (thal.loc[LATE_ONLY, "delta_log2FC"] > 1).all()

    # joint interaction model agrees
    it = tables["interaction"]
    assert (tables["interaction_log"]["status"] == "ok").all()
    ix = it[(it.cell_type == "OL") & (it.region == "Cortex") & (it.term == "age_x_condition")].set_index("gene")
    assert (ix.loc[LATE_ONLY, "padj"] < 0.05).all()
    assert (ix.loc[GLOBAL, "padj"] > 0.05).all()
    late = it[(it.cell_type == "OL") & (it.region == "Cortex") & (it.term == "condition_at_60")].set_index("gene")
    assert (late.loc[GLOBAL + LATE_ONLY, "log2FC"] > 1).all()

    # region vs rest: thalamus-specific genes flagged only for Thalamus
    rr = tables["region_vs_rest"]
    th = rr[(rr.cell_type == "OL") & (rr.region == "Thalamus") & (rr.age == "21")].set_index("gene")
    assert (th.loc[THAL_ONLY, "padj"] < 0.05).all()
    assert (th.loc[GLOBAL, "padj"] > 0.05).all()
    cx = rr[(rr.cell_type == "OL") & (rr.region == "Cortex") & (rr.age == "21")].set_index("gene")
    assert (cx.loc[GLOBAL, "padj"] > 0.05).all()

    # summary + persistence round trip
    summ = tables["summary"]
    row = summ[(summ.cell_type == "OL") & (summ.region == "Thalamus") & (summ.age == "60")].iloc[0]
    assert row["n_up"] >= 15
    loaded = rt.load_results(tmp_path)
    assert set(loaded) == set(k for k, v in tables.items() if not v.empty)
    assert loaded["simple_effects"]["age"].dtype == object
    mat = rt.gene_stratum_matrix(loaded["simple_effects"], "OL", GLOBAL + LATE_ONLY, mask_padj=0.05, min_base_mean=5)
    assert mat.shape[0] == 10 and mat.shape[1] == 8
    assert np.isnan(mat.loc[LATE_ONLY[0], "Cortex @ 21"])
    assert mat.loc[GLOBAL[0], "Cortex @ 21"] > 1


def test_plots_render(tmp_path=None):
    import matplotlib

    matplotlib.use("Agg")
    from utils import region_time_plots as rp

    tmp_path = Path(tmp_path) if tmp_path else REPO / "tests" / "_tmp_region_time"
    tables = rt.load_results(tmp_path)
    if not tables:
        tables = rt.run_all(make_adata(), tmp_path, layer="counts", n_cpus=2)
    figdir = tmp_path / "figures"
    rp.save(rp.plot_hits_overview(tables["summary"]), figdir / "overview")
    rp.save(rp.plot_age_scatter_grid(tables["age_contrast"], "OL", forced_genes=GLOBAL[:2]), figdir / "age_scatter_OL")
    rp.save(rp.plot_gene_heatmap(tables["simple_effects"], "OL", min_base_mean=5) if False else
            rp.plot_gene_heatmap(tables["simple_effects"], "OL", genes=GLOBAL + LATE_ONLY + THAL_ONLY), figdir / "heatmap_OL")
    rp.save(rp.plot_interaction_lollipop(tables["interaction"], "OL", "Cortex"), figdir / "interaction_OL_Cortex")
    rp.save(rp.plot_region_specificity(tables["region_vs_rest"], "OL", "21"), figdir / "region_specific_OL_21")
    se = tables["simple_effects"]
    one = se[(se.cell_type == "OL") & (se.region == "Thalamus") & (se.age == "60")]
    rp.save(rp.plot_volcano(one, "OL Thalamus P60", min_base_mean=5), figdir / "volcano")
    assert (figdir / "overview.png").exists()


if __name__ == "__main__":
    import shutil, time

    tmp = REPO / "tests" / "_tmp_region_time"
    shutil.rmtree(tmp, ignore_errors=True)
    for fn in (test_pseudobulk_shapes, test_min_cells_dropping, test_pipeline_recovers_planted_effects, test_plots_render):
        t0 = time.time()
        fn()
        print(f"PASS {fn.__name__} ({time.time() - t0:.1f}s)")
    print(f"figures in {tmp / 'figures'}")
