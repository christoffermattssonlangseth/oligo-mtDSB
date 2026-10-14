"""Predict a section's age cohort and condition from its expression profile.

Companion to ``utils.region_time_de``: where DE asks *which genes* change, this
module asks *how separable* the groups are, and *in which cell types and
regions* the signal lives.

Study design and what can be predicted
--------------------------------------
12 sections = 2 age cohorts (DOX removal at P21 or P60) x 2 conditions
(control, mtDSB) x 3. Every Xenium section was collected 42 days after
induction, so "time since induction" is constant and cannot be learned;
age at induction and age at sacrifice (~2 vs ~3.5 months) are the same split.
The prediction tasks are therefore (``default_tasks``):

``age``                 P21 vs P60 cohort, all 12 sections (positive control)
``condition``           mtDSB vs control, all 12 sections (ages pooled)
``joint``               4 classes age x condition (one-vs-rest)
``condition@21/@60``    mtDSB vs control within one cohort (6 sections)
``condition 21->60``    train mtDSB vs control on one cohort, test on the other
                        (does the treatment signature transfer across ages?)

Validation
----------
The section is the unit: every task is evaluated with leave-one-section-out
(LOSO) cross-validation, so a section's cells or pseudobulk never inform its
own prediction. With 12 sections the metrics are coarse; significance comes
from a sample-level label-permutation test that keeps the design balanced
(condition labels are permuted *within* age and vice versa). Exact
permutation p-values have a floor: 1/400 for age/condition, 1/20 for the
within-cohort and transfer tasks (3 vs 3 can never reach p < 0.05).

Without real signal, leave-one-out AUC is pessimistic and often far *below*
0.5 (holding a section out shifts the training classes against it), so read
AUC < 0.5 as "no signal", not as a reversed effect. The permutation null is
computed with the same LOSO procedure and carries the same bias, so the
p-values remain valid.

Two levels
----------
1. **Pseudobulk** (``classify_pseudobulk``): one feature vector per
   (cell type, region, section) = centered log-CPM of the summed counts,
   from the pseudobulk that ``region_time_de`` already writes. One model per
   stratum, plus a "Composition" stratum built from cell type x region
   proportions. Fast enough to permute every stratum.
2. **Single cell** (``classify_cells``): one model per cell type trained on a
   per-section subsample of cells, predicting every cell of the held-out
   section. Gives an out-of-fold score per cell (e.g. how mtDSB-like each
   oligodendrocyte looks, for spatial maps) and section-level predictions by
   averaging cell scores per section (and per region). No permutations.

The classifier is L2-regularised logistic regression with balanced class
weights. For the pseudobulk level it is solved in its kernel (dual) form,
which is exact and lets all permutations be fitted in one batched Newton
solve (checked against scikit-learn in ``tests/test_sample_classifier.py``).
Genes are centered but *not* standardized: with 5-11 training sections,
per-gene standardization inflates thousands of noise genes to the weight of
the few informative ones. The kernel is rescaled by its trace so ``C`` does
not depend on panel size.
"""

from __future__ import annotations

import itertools
import json
import logging
import math
import time
import zlib
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.special import expit

from utils.region_time_de import TISSUE_WIDE, Keys, Pseudobulk

logger = logging.getLogger(__name__)

COMPOSITION = "Composition"
"""Pseudo cell-type label of the cell type x region proportion features."""

ALL_CELL_TYPES = "All cell types"


# ----------------------------------------------------------------------------
# Tasks
# ----------------------------------------------------------------------------


@dataclass(frozen=True)
class Task:
    """One prediction problem.

    label
        ``"age"``, ``"condition"`` or ``"joint"`` (age|condition, 4 classes).
    pool
        Restrict to sections of this age (``None`` = all sections).
    train_age
        Transfer task: train on this age, test on the other (no LOSO).
    permute_within
        Sample-level column to permute labels within (keeps the design balanced).
    """

    name: str
    label: str
    pool: str | None = None
    train_age: str | None = None
    permute_within: str | None = None

    @property
    def is_transfer(self) -> bool:
        return self.train_age is not None


def default_tasks(keys: Keys = Keys()) -> list[Task]:
    e, l = keys.early, keys.late
    return [
        Task("age", "age", permute_within="condition"),
        Task("condition", "condition", permute_within="age"),
        Task("joint", "joint"),
        Task(f"condition@{e}", "condition", pool=e),
        Task(f"condition@{l}", "condition", pool=l),
        Task(f"condition {e}->{l}", "condition", train_age=e),
        Task(f"condition {l}->{e}", "condition", train_age=l),
    ]


def task_labels(meta: pd.DataFrame, task: Task, keys: Keys = Keys()) -> tuple[pd.Series, list[str]]:
    """Label per sample and the ordered classes (binary: [negative, positive])."""
    if task.label == "age":
        return meta["age"].astype(str), [keys.early, keys.late]
    if task.label == "condition":
        return meta["condition"].astype(str), [keys.control, keys.case]
    if task.label == "joint":
        classes = [f"{a}|{c}" for a in (keys.early, keys.late) for c in (keys.control, keys.case)]
        return meta["age"].astype(str) + "|" + meta["condition"].astype(str), classes
    raise ValueError(f"unknown label {task.label!r}")


def _folds(sample_meta: pd.DataFrame, task: Task) -> list[tuple[np.ndarray, np.ndarray]]:
    """(train, test) positional index pairs over ``sample_meta`` rows (already pooled)."""
    n = len(sample_meta)
    if task.is_transfer:
        tr = (sample_meta["age"] == task.train_age).to_numpy()
        return [(np.flatnonzero(tr), np.flatnonzero(~tr))]
    return [(np.delete(np.arange(n), i), np.array([i])) for i in range(n)]


def _pool(sample_meta: pd.DataFrame, task: Task) -> np.ndarray:
    if task.pool is None:
        return np.ones(len(sample_meta), dtype=bool)
    return (sample_meta["age"] == task.pool).to_numpy()


# ----------------------------------------------------------------------------
# Kernel ridge-logistic regression, batched over label vectors
# ----------------------------------------------------------------------------


def _kernel_features(Xtr: np.ndarray, Xte: np.ndarray, min_mean: float | None, standardize: bool = False):
    """Center (and optionally standardize) on the training rows and map both sets to the kernel's eigenbasis.

    Features are centered with the training mean, divided by the training SD
    if ``standardize``, and then globally rescaled so that the mean squared
    norm of a training row is 1 (``scale``); this makes ``C`` independent of
    the number of genes. Returns ``(Z_train, Z_test, back, keep, scale)`` where
    ``Z Z^T`` equals the rescaled training kernel, ``back`` (p_kept x r) maps a
    coefficient vector in Z space to per-unit weights on the centered (or
    standardized) features, and ``keep`` is the mask of features used. The fit
    is equivalent to scikit-learn ``LogisticRegression(C=C / scale)`` on those
    centered (standardized) features.
    """
    keep = np.ones(Xtr.shape[1], dtype=bool) if min_mean is None else Xtr.mean(0) >= min_mean
    sd = Xtr.std(0)
    keep &= sd > 1e-8
    if not keep.any():
        raise ValueError("no features left after filtering")
    mu = Xtr[:, keep].mean(0)
    div = sd[keep] if standardize else 1.0
    Xc = (Xtr[:, keep] - mu) / div
    Xt = (Xte[:, keep] - mu) / div
    scale = float((Xc**2).sum(1).mean())
    Xc, Xt = Xc / np.sqrt(scale), Xt / np.sqrt(scale)
    lam, U = np.linalg.eigh(Xc @ Xc.T)
    ok = lam > lam.max() * 1e-10
    lam, U = lam[ok], U[:, ok]
    Z = U * np.sqrt(lam)
    proj = U / np.sqrt(lam)  # Z_test = Xt Xc^T U lam^-1/2
    Zt = Xt @ Xc.T @ proj
    back = Xc.T @ proj / np.sqrt(scale)
    return Z, Zt, back, keep, scale


def _ridge_logistic_batch(Z: np.ndarray, Y: np.ndarray, C: float, n_iter: int = 100, tol: float = 1e-10):
    """Fit P independent L2 logistic regressions sharing the design ``Z``.

    Minimizes ``sum_i s_i * logloss_i + ||beta||^2 / (2C)`` (intercept
    unpenalized) with balanced class weights ``s``. ``Y`` is (P, n) in {0, 1}.
    Returns theta of shape (P, r + 1), the last column being the intercept.
    """
    n, r = Z.shape
    Y = np.asarray(Y, dtype=float)
    Zb = np.hstack([Z, np.ones((n, 1))])
    pen = np.r_[np.full(r, 1.0 / C), 0.0]
    n1 = np.clip(Y.sum(1, keepdims=True), 1, None)
    n0 = np.clip(n - Y.sum(1, keepdims=True), 1, None)
    S = np.where(Y == 1, n / (2 * n1), n / (2 * n0))

    def objective(th):
        F = th @ Zb.T
        return (S * (np.logaddexp(0, F) - Y * F)).sum(1) + 0.5 * (th**2 * pen).sum(1)

    theta = np.zeros((Y.shape[0], r + 1))
    obj = objective(theta)
    for _ in range(n_iter):
        Pr = expit(theta @ Zb.T)
        G = (S * (Pr - Y)) @ Zb + theta * pen
        W = S * Pr * (1 - Pr)
        H = np.einsum("pn,ni,nj->pij", W, Zb, Zb) + np.diag(pen + 1e-10)
        step = np.linalg.solve(H, G[..., None])[..., 0]
        t = np.ones(len(theta))
        for _ in range(30):  # backtracking per row
            new = objective(theta - t[:, None] * step)
            bad = new > obj + 1e-12
            if not bad.any():
                break
            t[bad] *= 0.5
        theta = theta - t[:, None] * step
        obj = objective(theta)
        if np.abs(t[:, None] * step).max() < tol:
            break
    return theta


def _predict(theta: np.ndarray, Z: np.ndarray) -> np.ndarray:
    """(P, n) probabilities of the positive class."""
    return expit(theta[:, :-1] @ Z.T + theta[:, -1:])


# ----------------------------------------------------------------------------
# Permutations
# ----------------------------------------------------------------------------


def _multiset_perms(labels: Sequence[str]):
    """All distinct orderings of a multiset."""
    counts = pd.Series(list(labels)).value_counts().to_dict()
    n = len(labels)
    out: list[str] = []

    def rec():
        if len(out) == n:
            yield tuple(out)
            return
        for lab in sorted(counts):
            if counts[lab]:
                counts[lab] -= 1
                out.append(lab)
                yield from rec()
                out.pop()
                counts[lab] += 1

    yield from rec()


def _n_multiset_perms(labels: Sequence[str]) -> int:
    c = pd.Series(list(labels)).value_counts()
    return math.factorial(len(labels)) // math.prod(math.factorial(int(k)) for k in c)


def permutation_labels(labels: pd.Series, blocks: pd.Series | None, n_perm: int, seed: int = 0) -> tuple[np.ndarray, bool]:
    """Label vectors for a permutation test; row 0 is the observed labelling.

    Permutes within ``blocks``. Enumerates every distinct labelling when there
    are at most ``n_perm`` of them (exact test, observed included once),
    otherwise draws ``n_perm`` random ones. Returns ``(L, exact)`` with ``L`` of
    shape (n_rows, n_samples).
    """
    lab = labels.to_numpy().astype(object)
    groups = [np.arange(len(lab))] if blocks is None else [np.flatnonzero(blocks.to_numpy() == b) for b in pd.unique(blocks)]
    total = math.prod(_n_multiset_perms(lab[g]) for g in groups)
    if total <= n_perm:
        per_block = [list(_multiset_perms(lab[g])) for g in groups]
        rows = []
        for combo in itertools.product(*per_block):
            row = lab.copy()
            for g, vals in zip(groups, combo):
                row[g] = vals
            rows.append(row)
        L = np.array(rows, dtype=object)
        obs_idx = int(np.flatnonzero((L == lab).all(1))[0])
        L[[0, obs_idx]] = L[[obs_idx, 0]]
        return L, True
    rng = np.random.default_rng(seed)
    rows = [lab.copy()]
    for _ in range(n_perm):
        row = lab.copy()
        for g in groups:
            row[g] = rng.permutation(lab[g])
        rows.append(row)
    return np.array(rows, dtype=object), False


def _perm_pvalue(stats: np.ndarray) -> float:
    """``stats[0]`` is observed; larger = better separation.

    The observed labelling is part of the null set in both the exact and the
    random case, so this is ``#(null >= obs) / n`` resp. ``(1 + #) / (1 + n_perm)``.
    """
    return float(np.mean(stats >= stats[0] - 1e-12))


# ----------------------------------------------------------------------------
# Metrics
# ----------------------------------------------------------------------------


def _auc(y: np.ndarray, score: np.ndarray) -> np.ndarray:
    """Rank AUC per row. ``y`` and ``score`` are (P, n); ties count 1/2."""
    y = np.atleast_2d(y).astype(bool)
    score = np.atleast_2d(score)
    out = np.full(len(score), np.nan)
    for i in range(len(score)):
        pos, neg = score[i, y[i]], score[i, ~y[i]]
        if len(pos) and len(neg):
            d = pos[:, None] - neg[None, :]
            out[i] = ((d > 0).sum() + 0.5 * (d == 0).sum()) / d.size
    return out


def score_predictions(pred: pd.DataFrame, group_cols: Sequence[str]) -> pd.DataFrame:
    """Section-level metrics per group from a predictions table.

    Expects columns ``true, pred, p_true`` and ``p_pos`` / ``positive`` for
    binary tasks (``p_pos`` NaN for multi-class).
    """
    rows = []
    for key, g in pred.groupby(list(group_cols), sort=False, observed=True):
        key = key if isinstance(key, tuple) else (key,)
        binary = g["p_pos"].notna().all()
        auc = _auc((g["true"] == g["positive"]).to_numpy()[None], g["p_pos"].to_numpy()[None])[0] if binary else np.nan
        rows.append(
            dict(
                zip(group_cols, key),
                n_samples=len(g),
                auc=auc,
                accuracy=float((g["true"] == g["pred"]).mean()),
                mean_p_true=float(g["p_true"].mean()),
            )
        )
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------
# Core: one feature matrix (samples x features), one task
# ----------------------------------------------------------------------------


def _min_class_ok(labels: pd.Series, classes: list[str], task: Task, min_per_class: int) -> bool:
    counts = labels.value_counts()
    if task.is_transfer:
        return all(counts.get(c, 0) >= 2 for c in classes)  # checked per age below
    return all(counts.get(c, 0) >= min_per_class for c in classes)


def predict_task(
    X: pd.DataFrame,
    sample_meta: pd.DataFrame,
    task: Task,
    keys: Keys = Keys(),
    C: float = 10.0,
    min_mean: float | None = None,
    n_perm: int = 1000,
    seed: int = 0,
    min_per_class: int | None = None,
    standardize: bool = False,
) -> tuple[pd.DataFrame, dict] | None:
    """Out-of-fold predictions and a permutation p-value for one task on one feature matrix.

    ``X`` (samples x features) and ``sample_meta`` (columns ``age, condition``)
    share the index (section ids). Returns ``(predictions, info)`` or ``None``
    if the task is not estimable on these samples.
    """
    pool = _pool(sample_meta, task)
    meta = sample_meta.loc[pool]
    Xv = X.loc[meta.index].to_numpy(dtype=float)
    if not (Xv.std(0) > 1e-8).any():
        return None
    labels, classes = task_labels(meta, task, keys)
    binary = len(classes) == 2
    if min_per_class is None:
        min_per_class = 3 if binary else 2
    if not _min_class_ok(labels, classes, task, min_per_class):
        return None
    if task.is_transfer:
        for age_val, g in labels.groupby(meta["age"]):
            if not all((g == c).sum() >= 2 for c in classes):
                return None
        if meta["age"].nunique() < 2:
            return None

    # label matrix: transfer tasks permute only the training labels
    if task.is_transfer:
        tr_mask = (meta["age"] == task.train_age).to_numpy()
        Ltr, exact = permutation_labels(labels[tr_mask], None, n_perm, seed)
        L = np.empty((len(Ltr), len(labels)), dtype=object)
        L[:, tr_mask] = Ltr
        L[:, ~tr_mask] = labels[~tr_mask].to_numpy()
    else:
        blocks = meta[task.permute_within] if task.permute_within else None
        L, exact = permutation_labels(labels, blocks, n_perm, seed)

    n = len(meta)
    probs = np.zeros((len(L), n, len(classes)))
    for tr, te in _folds(meta, task):
        Z, Zt, *_ = _kernel_features(Xv[tr], Xv[te], min_mean, standardize)
        if binary:
            theta = _ridge_logistic_batch(Z, (L[:, tr] == classes[1]).astype(float), C)
            p = _predict(theta, Zt)
            probs[:, te, 1] = p
            probs[:, te, 0] = 1 - p
        else:
            for k, c in enumerate(classes):
                theta = _ridge_logistic_batch(Z, (L[:, tr] == c).astype(float), C)
                probs[:, te, k] = _predict(theta, Zt)
            s = probs[:, te, :].sum(-1, keepdims=True)
            probs[:, te, :] /= s

    test_idx = np.concatenate([te for _, te in _folds(meta, task)])
    cls_idx = {c: k for k, c in enumerate(classes)}
    Lk = np.vectorize(cls_idx.get)(L)
    p_true = np.take_along_axis(probs, Lk[..., None], axis=-1)[..., 0]
    if binary:
        stat = _auc(Lk[:, test_idx] == 1, probs[:, test_idx, 1])
    else:
        stat = p_true[:, test_idx].mean(1)
    perm_p = _perm_pvalue(stat)

    pred = pd.DataFrame(
        {
            "sample_id": meta.index[test_idx],
            "age": meta["age"].to_numpy()[test_idx],
            "condition": meta["condition"].to_numpy()[test_idx],
            "true": labels.to_numpy()[test_idx],
            "pred": np.array(classes)[probs[0, test_idx].argmax(-1)],
            "p_true": p_true[0, test_idx],
            "p_pos": probs[0, test_idx, 1] if binary else np.nan,
            "positive": classes[1] if binary else None,
        }
    )
    if not binary:
        for k, c in enumerate(classes):
            pred[f"p[{c}]"] = probs[0, test_idx, k]
    info = {"perm_p": perm_p, "n_perm": int(len(L) - 1), "exact": exact, "perm_stat": "auc" if binary else "mean_p_true"}
    return pred, info


def fit_coefficients(
    X: pd.DataFrame,
    sample_meta: pd.DataFrame,
    task: Task,
    keys: Keys = Keys(),
    C: float = 10.0,
    min_mean: float | None = None,
    standardize: bool = False,
) -> pd.Series | None:
    """Feature weights (per unit of centered log-CPM, or per SD if ``standardize``) of a binary task fitted on all its training sections."""
    pool = _pool(sample_meta, task)
    meta = sample_meta.loc[pool]
    if task.is_transfer:
        meta = meta.loc[meta["age"] == task.train_age]
    labels, classes = task_labels(meta, task, keys)
    if len(classes) != 2 or labels.nunique() < 2:
        return None
    Xv = X.loc[meta.index].to_numpy(dtype=float)
    Z, _, back, keep, _ = _kernel_features(Xv, Xv[:1], min_mean, standardize)
    theta = _ridge_logistic_batch(Z, (labels.to_numpy() == classes[1]).astype(float)[None], C)
    return pd.Series(back @ theta[0, :-1], index=X.columns[keep])


# ----------------------------------------------------------------------------
# Pseudobulk level
# ----------------------------------------------------------------------------


def log_cpm(counts: pd.DataFrame) -> pd.DataFrame:
    lib = counts.sum(1).replace(0, 1)
    return np.log1p(counts.div(lib, axis=0) * 1e6)


def composition_features(pb: Pseudobulk, pseudocount: float = 1e-4) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Log proportion of every (cell type, region) among a section's cells.

    Uses kept and dropped pseudobulk rows so small populations still count.
    Tissue-wide rows are ignored (they double-count).
    """
    meta = pd.concat([pb.meta, pb.dropped]) if len(pb.dropped) else pb.meta
    meta = meta.loc[meta["region"] != TISSUE_WIDE].copy()
    meta["n_cells"] = meta["n_cells"].astype(int)
    tab = meta.pivot_table(index="sample_id", columns=["cell_type", "region"], values="n_cells", aggfunc="sum", fill_value=0)
    prop = tab.div(tab.sum(1), axis=0)
    X = np.log(prop + pseudocount)
    X.columns = [f"{c}|{r}" for c, r in X.columns]
    smeta = meta.groupby("sample_id")[["age", "condition"]].first().loc[X.index]
    return X, smeta


def _stratum_inputs(pb: Pseudobulk, min_cells: int):
    """Yield (cell_type, region, X logCPM, sample_meta) per stratum."""
    keep = pb.meta["n_cells"].astype(int) >= min_cells
    meta = pb.meta.loc[keep]
    lc = log_cpm(pb.counts.loc[keep])
    for (ct, reg), m in meta.groupby(["cell_type", "region"], sort=True):
        X = lc.loc[m.index]
        X.index = m["sample_id"].to_numpy()
        smeta = m.set_index("sample_id")[["age", "condition"]]
        yield ct, reg, X, smeta


def classify_pseudobulk(
    pb: Pseudobulk,
    keys: Keys = Keys(),
    tasks: Sequence[Task] | None = None,
    C: float = 10.0,
    min_cpm: float = 10.0,
    min_cells: int = 20,
    n_perm: int = 1000,
    n_top_genes: int = 30,
    include_composition: bool = True,
    seed: int = 0,
) -> dict[str, pd.DataFrame]:
    """Run every task on every (cell type, region) pseudobulk stratum.

    Returns ``{"predictions", "metrics", "coefficients", "sample_summary"}``:

    * predictions: one row per (task, stratum, held-out section)
    * metrics: AUC / accuracy / mean p(true) per (task, stratum), permutation
      p-value and BH ``padj`` across strata within a task
    * coefficients: top ``n_top_genes`` genes per side for binary tasks
      (fitted on all training sections; weight per unit centered log-CPM)
    * sample_summary: per (task, section) median p(true) across strata and the
      fraction of strata that classify it correctly (flags odd sections)
    """
    tasks = list(tasks or default_tasks(keys))
    min_mean = float(np.log1p(min_cpm)) if min_cpm else None
    strata = list(_stratum_inputs(pb, min_cells))
    if include_composition:
        Xc, mc = composition_features(pb)
        strata.append((COMPOSITION, "All regions", Xc, mc))

    preds, infos, coefs = [], [], []
    t0 = time.time()
    for i, (ct, reg, X, smeta) in enumerate(strata):
        mm = None if ct == COMPOSITION else min_mean
        for task in tasks:
            std = ct == COMPOSITION  # proportions span orders of magnitude; genes are only centered
            res = predict_task(X, smeta, task, keys, C=C, min_mean=mm, n_perm=n_perm, seed=seed, standardize=std)
            if res is None:
                continue
            pr, info = res
            pr.insert(0, "region", reg)
            pr.insert(0, "cell_type", ct)
            pr.insert(0, "task", task.name)
            preds.append(pr)
            infos.append({"task": task.name, "cell_type": ct, "region": reg, **info})
            if n_top_genes and task.label != "joint":
                w = fit_coefficients(X, smeta, task, keys, C=C, min_mean=mm, standardize=std)
                if w is not None:
                    top = pd.concat([w.nlargest(n_top_genes), w.nsmallest(n_top_genes)])
                    coefs.append(pd.DataFrame({"task": task.name, "cell_type": ct, "region": reg, "gene": top.index, "coef": top.to_numpy()}))
        if (i + 1) % 20 == 0:
            logger.info("%d/%d strata (%.0f s)", i + 1, len(strata), time.time() - t0)

    predictions = pd.concat(preds, ignore_index=True)
    metrics = score_predictions(predictions, ["task", "cell_type", "region"]).merge(
        pd.DataFrame(infos), on=["task", "cell_type", "region"]
    )
    metrics["padj"] = metrics.groupby("task")["perm_p"].transform(_bh)
    coefficients = pd.concat(coefs, ignore_index=True) if coefs else pd.DataFrame()
    return {
        "predictions": predictions,
        "metrics": metrics,
        "coefficients": coefficients,
        "sample_summary": sample_summary(predictions),
    }


def sample_summary(predictions: pd.DataFrame) -> pd.DataFrame:
    """Per (task, section): how confidently and how often it is classified correctly."""
    g = predictions.assign(correct=predictions["true"] == predictions["pred"]).groupby(["task", "sample_id"], sort=False)
    out = g.agg(
        age=("age", "first"),
        condition=("condition", "first"),
        n_strata=("p_true", "size"),
        median_p_true=("p_true", "median"),
        frac_correct=("correct", "mean"),
    ).reset_index()
    return out


def _bh(p: pd.Series) -> pd.Series:
    from statsmodels.stats.multitest import multipletests

    out = pd.Series(np.nan, index=p.index)
    ok = p.notna()
    if ok.any():
        out[ok] = multipletests(p[ok], method="fdr_bh")[1]
    return out


# ----------------------------------------------------------------------------
# Single-cell level
# ----------------------------------------------------------------------------


def _lognorm_cells(X: sparse.csr_matrix, target: float | None = None) -> tuple[sparse.csr_matrix, float]:
    tot = np.asarray(X.sum(1)).ravel()
    if target is None:
        target = float(np.median(tot[tot > 0]))
    scale = target / np.where(tot > 0, tot, 1)
    Y = sparse.diags(scale) @ X
    Y = sparse.csr_matrix(Y, dtype=np.float32)
    Y.data = np.log1p(Y.data)
    return Y, target


def classify_cells(
    adata,
    keys: Keys = Keys(),
    layer: str | None = "counts",
    tasks: Sequence[Task] | None = None,
    cell_types: Sequence[str] | None = None,
    cells_per_sample: int = 2000,
    min_cells_per_sample: int = 50,
    C: float = 0.1,
    exclude_regions: Sequence[str] = ("Unknown",),
    seed: int = 0,
    chunk: int = 50_000,
    outdir: str | Path | None = None,
    resume: bool = True,
    cell_type_map: dict[str, str] | None = None,
) -> dict[str, pd.DataFrame]:
    """Leave-one-section-out cell classifiers, one per cell type and task.

    Each fold trains ``LogisticRegression`` (L2, balanced classes) on up to
    ``cells_per_sample`` cells per training section (log-normalized,
    standardized genes) and scores *every* cell of that type in the held-out
    section(s). Sections with fewer than ``min_cells_per_sample`` cells of the
    type are left out of that cell type.

    With ``outdir``, results are checkpointed after every cell type: the cell
    scores go to ``outdir/cell_scores/<cell type>.csv.gz`` and
    ``predictions.csv`` / ``metrics.csv`` / ``progress.json`` are rewritten with
    everything finished so far, so a running job can be followed. With
    ``resume``, cell types whose score file already exists are loaded instead
    of refitted. Cell subsampling is seeded per cell type, so resumed and
    uninterrupted runs give identical results.

    ``cell_type_map`` relabels cell types before anything else (e.g.
    ``{"DA-Oligodendrocytes": "Oligodendrocytes"}`` to avoid disease-defined
    clusters); unmapped labels are kept.

    Returns
    -------
    dict with
      * ``cell_scores``: one row per scored cell (index = obs name) with
        ``cell_type, region, sample_id`` and per task ``p_pos[<task>]`` (binary;
        positive = late age / mtDSB) or ``p[<task>][<class>]`` (joint).
      * ``predictions``: section-level predictions (mean cell probability per
        section) per (task, cell type, region incl. ``All regions``).
      * ``metrics``: ``score_predictions`` of the above.
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler

    tasks = list(tasks or default_tasks(keys))
    cols = [keys.cell_type, keys.region, keys.sample, keys.age, keys.condition]
    obs = adata.obs[cols].astype(str).copy()
    obs.columns = ["cell_type", "region", "sample_id", "age", "condition"]
    if cell_type_map:
        obs["cell_type"] = obs["cell_type"].map(lambda c: cell_type_map.get(c, c))
    X = adata.layers[layer] if layer else adata.X
    keep = ~obs["region"].isin([str(r) for r in exclude_regions])
    obs = obs.loc[keep]
    Xn, _ = _lognorm_cells(sparse.csr_matrix(X)[np.flatnonzero(keep.to_numpy())])
    cts = list(cell_types) if cell_types else sorted(obs["cell_type"].unique())
    smeta_sections = obs.groupby("sample_id")[["age", "condition"]].first()
    score_dir = Path(outdir) / "cell_scores" if outdir is not None else None
    if score_dir is not None:
        score_dir.mkdir(parents=True, exist_ok=True)

    score_frames = []
    t0 = time.time()
    for ct in cts:
        ct_file = score_dir / f"{_safe_name(ct)}.csv.gz" if score_dir is not None else None
        if ct_file is not None and resume and ct_file.exists():
            score_frames.append(pd.read_csv(ct_file, index_col=0))
            logger.info("cells: %s loaded from %s", ct, ct_file)
            continue
        rng = np.random.default_rng([seed, zlib.crc32(ct.encode())])
        ct_pos = np.flatnonzero((obs["cell_type"] == ct).to_numpy())
        ct_obs = obs.iloc[ct_pos]
        n_per = ct_obs["sample_id"].value_counts()
        ok_samples = n_per.index[n_per >= min_cells_per_sample]
        smeta_all = ct_obs.groupby("sample_id")[["age", "condition"]].first().loc[sorted(ok_samples)]
        scores = pd.DataFrame(index=ct_obs.index[ct_obs["sample_id"].isin(ok_samples)])
        sample_of = ct_obs["sample_id"].to_numpy()
        for task in tasks:
            pool = _pool(smeta_all, task)
            smeta = smeta_all.loc[pool]
            labels, classes = task_labels(smeta, task, keys)
            binary = len(classes) == 2
            need = 3 if binary and not task.is_transfer else 2
            if not _min_class_ok(labels, classes, task, need) or (task.is_transfer and smeta["age"].nunique() < 2):
                continue
            out_cols = [f"p_pos[{task.name}]"] if binary else [f"p[{task.name}][{c}]" for c in classes]
            for c in out_cols:
                scores[c] = np.nan
            for tr, te in _folds(smeta, task):
                tr_s, te_s = smeta.index[tr], smeta.index[te]
                tr_lab = labels.iloc[tr]
                if tr_lab.nunique() < len(classes):
                    continue
                idx = []
                for s in tr_s:
                    pos = np.flatnonzero(sample_of == s)
                    if len(pos) > cells_per_sample:
                        pos = rng.choice(pos, cells_per_sample, replace=False)
                    idx.append(pos)
                idx = np.concatenate(idx)
                y = pd.Series(sample_of[idx]).map(labels).to_numpy()
                Xtr = Xn[ct_pos[idx]].toarray()
                scaler = StandardScaler().fit(Xtr)
                clf = LogisticRegression(C=C, class_weight="balanced", max_iter=1000)
                clf.fit(scaler.transform(Xtr), y)
                order = [list(clf.classes_).index(c) for c in classes]
                te_pos = np.flatnonzero(np.isin(sample_of, te_s))
                for start in range(0, len(te_pos), chunk):
                    sl = te_pos[start : start + chunk]
                    P = clf.predict_proba(scaler.transform(Xn[ct_pos[sl]].toarray()))[:, order]
                    names = ct_obs.index[sl]
                    if binary:
                        scores.loc[names, out_cols[0]] = P[:, 1]
                    else:
                        scores.loc[names, out_cols] = P
            logger.info("cells: %s / %s done (%.0f s)", ct, task.name, time.time() - t0)
        scores.insert(0, "sample_id", ct_obs.loc[scores.index, "sample_id"])
        scores.insert(0, "region", ct_obs.loc[scores.index, "region"])
        scores.insert(0, "cell_type", ct)
        score_frames.append(scores)
        if ct_file is not None:
            scores.to_csv(ct_file, compression="gzip")
            _checkpoint(outdir, score_frames, smeta_sections, tasks, keys, done=len(score_frames), total=len(cts), t0=t0)
        logger.info("cells: %s done (%.0f s)", ct, time.time() - t0)

    if outdir is not None:
        _checkpoint(outdir, score_frames, smeta_sections, tasks, keys, done=len(score_frames), total=len(cts), t0=t0)
    cell_scores = pd.concat(score_frames)
    predictions = aggregate_cell_scores(cell_scores, smeta_sections, tasks, keys)
    metrics = score_predictions(predictions, ["task", "cell_type", "region"])
    return {"cell_scores": cell_scores, "predictions": predictions, "metrics": metrics}


def _safe_name(s: str) -> str:
    return "".join(c if c.isalnum() or c in "-_" else "_" for c in s)


def _checkpoint(outdir, score_frames, sample_meta, tasks, keys, done, total, t0) -> None:
    """Rewrite predictions / metrics from all cell types finished so far."""
    outdir = Path(outdir)
    preds = aggregate_cell_scores(pd.concat(score_frames), sample_meta, tasks, keys)
    preds.to_csv(outdir / "predictions.csv", index=False)
    score_predictions(preds, ["task", "cell_type", "region"]).to_csv(outdir / "metrics.csv", index=False)
    progress = {"cell_types_done": done, "cell_types_total": total, "seconds": round(time.time() - t0, 1),
                "done": sorted(set(pd.concat([f[["cell_type"]] for f in score_frames])["cell_type"]))}
    (outdir / "progress.json").write_text(json.dumps(progress, indent=2))


def aggregate_cell_scores(
    cell_scores: pd.DataFrame, sample_meta: pd.DataFrame, tasks: Sequence[Task], keys: Keys = Keys(), min_cells: int = 20
) -> pd.DataFrame:
    """Section-level predictions = mean cell probability per (cell type, region, section).

    Also adds region ``All regions`` (all cells of the type in the section).
    Groups with fewer than ``min_cells`` scored cells are skipped.
    """
    cs = pd.concat([cell_scores, cell_scores.assign(region=TISSUE_WIDE)])
    rows = []
    for task in tasks:
        _, classes = task_labels(sample_meta.iloc[:1], task, keys)
        binary = len(classes) == 2
        cols = [f"p_pos[{task.name}]"] if binary else [f"p[{task.name}][{c}]" for c in classes]
        if cols[0] not in cs:
            continue
        sub = cs.dropna(subset=cols)
        agg = sub.groupby(["cell_type", "region", "sample_id"], observed=True)[cols].agg(["mean", "size"])
        n = agg[(cols[0], "size")]
        means = agg.xs("mean", axis=1, level=1).loc[n >= min_cells]
        if means.empty:
            continue
        means = means.reset_index()
        sm = sample_meta.loc[means["sample_id"]]
        labels, _ = task_labels(sm, task, keys)
        df = pd.DataFrame(
            {
                "task": task.name,
                "cell_type": means["cell_type"],
                "region": means["region"],
                "sample_id": means["sample_id"],
                "age": sm["age"].to_numpy(),
                "condition": sm["condition"].to_numpy(),
                "true": labels.to_numpy(),
            }
        )
        P = means[cols].to_numpy()
        if binary:
            P = np.column_stack([1 - P[:, 0], P[:, 0]])
        df["pred"] = np.array(classes)[P.argmax(1)]
        df["p_true"] = P[np.arange(len(P)), [classes.index(t) for t in df["true"]]]
        df["p_pos"] = P[:, 1] if binary else np.nan
        df["positive"] = classes[1] if binary else None
        if not binary:
            for k, c in enumerate(classes):
                df[f"p[{c}]"] = P[:, k]
        rows.append(df)
    return pd.concat(rows, ignore_index=True)


# ----------------------------------------------------------------------------
# I/O
# ----------------------------------------------------------------------------


def save_results(outdir: str | Path, tables: dict[str, pd.DataFrame], run_info: dict | None = None) -> None:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    for name, df in tables.items():
        if name == "cell_scores":
            continue  # written per cell type by classify_cells(outdir=...)
        df.to_csv(outdir / f"{name}.csv", index=False)
    if run_info is not None:
        (outdir / "run_info.json").write_text(json.dumps(run_info, indent=2, default=str))


def load_results(outdir: str | Path) -> dict[str, pd.DataFrame]:
    outdir = Path(outdir)
    out = {}
    for f in outdir.glob("*.csv"):
        out[f.stem] = pd.read_csv(f)
    files = sorted((outdir / "cell_scores").glob("*.csv.gz"))
    if files:
        out["cell_scores"] = pd.concat([pd.read_csv(f, index_col=0) for f in files])
    return out


def tasks_info(tasks: Sequence[Task]) -> list[dict]:
    return [asdict(t) for t in tasks]
