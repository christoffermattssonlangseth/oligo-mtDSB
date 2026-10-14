"""Tests for utils.sample_classifier on a small synthetic Xenium-like data set.

Run with pytest, or directly:  python tests/test_sample_classifier.py

Synthetic design = the real one (12 sections = 2 ages x 2 conditions x 3),
2 cell types, 3 regions, 200 genes, plus per-section noise. Planted effects:
  * genes g0-g9   : up at age 60 in every cell type (age signal everywhere)
  * genes g10-g19 : up in mtDSB in "OL" only, both ages (condition signal in OL)
"Astro" carries no condition signal, so its condition AUC should be near chance.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from utils import region_time_de as rt  # noqa: E402
from utils import sample_classifier as sc  # noqa: E402

REGIONS = ["Cortex", "Thalamus", "Hypothalamus"]


def make_adata(seed: int = 0, cells_per_group: int = 40, n_genes: int = 200):
    import anndata as ad
    from scipy import sparse

    rng = np.random.default_rng(seed)
    genes = [f"g{i}" for i in range(n_genes)]
    base = rng.gamma(2.0, 1.0, size=n_genes) + 0.3
    rows, obs = [], []
    for age in ("21", "60"):
        for cond in ("control", "mtDSB"):
            for r in range(3):
                sid = f"RB{age}{cond[0]}{r}"
                section_noise = rng.lognormal(0, 0.15, size=n_genes)
                for ct in ("OL", "Astro"):
                    for reg in REGIONS:
                        mu = base * section_noise
                        if age == "60":
                            mu[0:10] *= 2.5
                        if ct == "OL" and cond == "mtDSB":
                            mu[10:20] *= 2.5
                        x = rng.poisson(mu[None, :] * rng.lognormal(0, 0.3, size=(cells_per_group, 1)))
                        rows.append(sparse.csr_matrix(x))
                        obs.append(pd.DataFrame({"sample_id": sid, "age": age, "condition": cond, "cell_class": ct,
                                                 "RBD_compartment_simplified": reg}, index=range(cells_per_group)))
    X = sparse.vstack(rows).tocsr()
    obs = pd.concat(obs, ignore_index=True)
    obs.index = [f"c{i}" for i in range(len(obs))]
    a = ad.AnnData(X=X.astype(np.float32), obs=obs)
    a.var_names = genes
    a.layers["counts"] = X.copy()
    return a


def test_kernel_logistic_matches_sklearn():
    from sklearn.linear_model import LogisticRegression

    rng = np.random.default_rng(1)
    X = rng.normal(size=(11, 300))
    y = np.array([0, 1] * 5 + [0])
    X[y == 1, :5] += 1.0
    Xte = rng.normal(size=(3, 300))
    C = 2.0
    for standardize in (False, True):
        Z, Zt, back, keep, scale = sc._kernel_features(X, Xte, None, standardize)
        theta = sc._ridge_logistic_batch(Z, y[None].astype(float), C)
        p_ours = sc._predict(theta, Zt)[0]

        mu, sd = X.mean(0), (X.std(0) if standardize else 1.0)
        ref = LogisticRegression(C=C / scale, class_weight="balanced", max_iter=10_000, tol=1e-12)
        ref.fit((X - mu) / sd, y)
        p_ref = ref.predict_proba((Xte - mu) / sd)[:, 1]
        assert np.allclose(p_ours, p_ref, atol=1e-4), (standardize, p_ours, p_ref)
        # feature weights match too
        assert np.allclose(back @ theta[0, :-1], ref.coef_[0], atol=1e-4)


def test_permutation_labels():
    meta = pd.DataFrame({"age": ["21"] * 6 + ["60"] * 6, "condition": ["control", "mtDSB"] * 6})
    L, exact = sc.permutation_labels(meta["condition"], meta["age"], n_perm=1000)
    assert exact and L.shape == (400, 12)  # C(6,3)^2
    assert (L[0] == meta["condition"].to_numpy()).all()
    assert len({tuple(r) for r in L}) == 400
    # within-age balance is preserved
    for row in L:
        assert (row[:6] == "mtDSB").sum() == 3 and (row[6:] == "mtDSB").sum() == 3
    L, exact = sc.permutation_labels(meta["age"] + "|" + meta["condition"], None, n_perm=200)
    assert not exact and L.shape == (201, 12)


def _pb():
    return rt.pseudobulk(make_adata(), layer="counts", min_cells_per_sample=20)


def test_pseudobulk_classifier_recovers_planted_signal():
    res = sc.classify_pseudobulk(_pb(), n_perm=400, min_cpm=0)
    m = res["metrics"].set_index(["task", "cell_type", "region"])

    # age is planted everywhere
    for ct in ("OL", "Astro"):
        assert m.loc[("age", ct, rt.TISSUE_WIDE), "auc"] == 1.0
        assert m.loc[("age", ct, rt.TISSUE_WIDE), "perm_p"] <= 0.01
    # condition only in OL
    assert m.loc[("condition", "OL", rt.TISSUE_WIDE), "auc"] >= 0.9
    assert m.loc[("condition", "OL", rt.TISSUE_WIDE), "perm_p"] <= 0.02
    assert m.loc[("condition", "Astro", rt.TISSUE_WIDE), "perm_p"] > 0.05
    # transfer of the OL signature across ages works (AUC on 6 held-out sections)
    assert m.loc[("condition 21->60", "OL", rt.TISSUE_WIDE), "auc"] >= 0.88
    # joint: OL gets both factors -> well above chance (0.25)
    assert m.loc[("joint", "OL", rt.TISSUE_WIDE), "accuracy"] >= 0.75
    # exact permutation floors
    assert m.loc[("condition@21", "OL", rt.TISSUE_WIDE), "n_perm"] == 19
    assert m.loc[("condition", "OL", rt.TISSUE_WIDE), "n_perm"] == 399

    # top coefficients for OL condition are the planted genes
    co = res["coefficients"]
    top = co[(co.task == "condition") & (co.cell_type == "OL") & (co.region == rt.TISSUE_WIDE)].nlargest(10, "coef")
    assert len(set(top.gene) & {f"g{i}" for i in range(10, 20)}) >= 8

    # every section appears once per (task, stratum) in LOSO tasks
    pr = res["predictions"]
    one = pr[(pr.task == "condition") & (pr.cell_type == "OL") & (pr.region == "Cortex")]
    assert sorted(one.sample_id) == sorted(_pb().meta.sample_id.unique())
    assert set(res["sample_summary"].columns) >= {"median_p_true", "frac_correct"}
    # composition features: one log-proportion per (cell type, region); constant here -> stratum skipped
    Xc, mc = sc.composition_features(_pb())
    assert Xc.shape == (12, 2 * len(REGIONS)) and np.allclose(np.exp(Xc).sum(1), 1, atol=1e-2)
    assert sc.COMPOSITION not in set(res["metrics"].cell_type)


def test_cell_classifier_runs_and_separates():
    a = make_adata(cells_per_group=30)
    res = sc.classify_cells(a, cells_per_sample=60, tasks=[t for t in sc.default_tasks() if t.name in ("age", "condition", "joint")])
    cs = res["cell_scores"]
    assert len(cs) == a.n_obs
    assert cs["p_pos[condition]"].notna().all()
    m = res["metrics"].set_index(["task", "cell_type", "region"])
    assert m.loc[("condition", "OL", rt.TISSUE_WIDE), "auc"] >= 0.9
    assert m.loc[("age", "Astro", rt.TISSUE_WIDE), "auc"] >= 0.9
    assert m.loc[("condition", "Astro", rt.TISSUE_WIDE), "auc"] < 0.9
    # relabelling merges cell types before fitting
    merged = sc.classify_cells(a, cells_per_sample=60, tasks=[t for t in sc.default_tasks() if t.name == "age"],
                               cell_type_map={"Astro": "Glia", "OL": "Glia"})
    assert set(merged["cell_scores"].cell_type) == {"Glia"} and len(merged["cell_scores"]) == a.n_obs


def test_cell_classifier_coefficients_find_planted_genes():
    from utils import ol_response as olr

    a = make_adata(cells_per_group=30)
    tasks = [t for t in sc.default_tasks() if t.name in ("condition", "condition 21->60", "joint")]
    res = sc.classify_cells(a, cells_per_sample=60, tasks=tasks, cell_types=["OL"], save_coefficients=True)
    co = res["coefficients"]
    assert set(co.task) == {"condition", "condition 21->60"}  # joint is multi-class: no coefficients
    assert co[co.task == "condition"].fold.nunique() == 12 and set(co[co.task != "condition"].fold) == {"train 21"}
    summ = olr.summarize_coefficients(co)
    top = summ[summ.task == "condition"].nlargest(10, "mean_coef")
    assert len(set(top.gene) & {f"g{i}" for i in range(10, 20)}) >= 8
    assert (top.sign_consistency == 1).all()


def test_cell_classifier_checkpoints_and_resumes(tmp_path=None):
    import shutil

    out = Path(tmp_path) if tmp_path else REPO / "tests" / "_tmp_sample_classifier"
    shutil.rmtree(out, ignore_errors=True)
    a = make_adata(cells_per_group=20)
    tasks = [t for t in sc.default_tasks() if t.name == "condition"]
    # first run only does OL -> partial results readable
    sc.classify_cells(a, cells_per_sample=40, tasks=tasks, cell_types=["OL"], outdir=out)
    part = sc.load_results(out)
    assert set(part["metrics"].cell_type) == {"OL"} and set(part["cell_scores"].cell_type) == {"OL"}
    # second run loads OL from disk, fits Astro; identical to a fresh uninterrupted run
    full = sc.classify_cells(a, cells_per_sample=40, tasks=tasks, cell_types=["Astro", "OL"], outdir=out)
    fresh = sc.classify_cells(a, cells_per_sample=40, tasks=tasks, cell_types=["Astro", "OL"])
    cols = ["p_pos[condition]"]
    assert np.allclose(full["cell_scores"].loc[fresh["cell_scores"].index, cols], fresh["cell_scores"][cols])
    assert set(sc.load_results(out)["metrics"].cell_type) == {"OL", "Astro"}
    shutil.rmtree(out, ignore_errors=True)


if __name__ == "__main__":
    import time

    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            t = time.time()
            fn()
            print(f"{name}: ok ({time.time() - t:.1f} s)")
