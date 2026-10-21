"""Re-draw the region map with disease-free genes: label transfer from control sections.

The RBD compartments (``RBD_compartment_simplified``) were called by clustering neighbourhood-pooled transcriptomes
of all sections together, using every gene. Where disease changes many genes (stressed OLs, microglia, astrocytes),
it can move a region into another cluster. This module re-assigns regions from genes that do not respond:

* ``select_anatomy_genes``  genes that differ most between regions in control pseudobulks (ANOVA-like F on log2 CPM),
                            excluding a given set of response genes
* ``neighbour_pool``        per section, sum each cell's counts with its ``k`` nearest neighbours (as the RBD pipeline)
* ``RegionTransfer``        log-normalise, PCA and multinomial logistic regression fitted on control cells with their
                            original labels; ``predict`` assigns every cell, ``loco_accuracy`` holds out one control
                            section at a time

Controls are unaffected by disease, so their original labels are the reference; mtDSB sections are then labelled by
the same disease-free model.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial import cKDTree
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression


def select_anatomy_genes(log_cpm: pd.DataFrame, regions: pd.Series, exclude=(), n: int = 400,
                         min_cpm: float = 20.0) -> pd.Series:
    """``log_cpm``: (region, sample) pseudobulks x genes, log2 CPM, controls only; ``regions`` labels its rows.

    Returns the F score (between-region variance / mean within-region variance) of the top ``n`` genes.
    """
    means = log_cpm.groupby(regions.to_numpy()).mean()
    within = log_cpm.groupby(regions.to_numpy()).var().mean()
    F = means.var() / (within + 0.05)
    ok = (np.exp2(log_cpm) - 1).mean() > min_cpm
    return F[ok & ~F.index.isin(set(exclude))].sort_values(ascending=False).head(n)


def neighbour_pool(X: sparse.csr_matrix, sample: np.ndarray, xy: np.ndarray, k: int = 20) -> sparse.csr_matrix:
    """Sum of each cell's counts and those of its ``k`` nearest neighbours in the same section (row order kept)."""
    X = sparse.csr_matrix(X, dtype=np.float32)
    blocks, idx = [], []
    for s in pd.unique(sample):
        ii = np.flatnonzero(sample == s)
        kk = min(k + 1, len(ii))
        _, nb = cKDTree(xy[ii]).query(xy[ii], k=kk)
        nb = nb.reshape(len(ii), kk)
        C = sparse.csr_matrix((np.ones(nb.size, np.float32), (np.repeat(np.arange(len(ii)), kk), nb.ravel())), shape=(len(ii), len(ii)))
        blocks.append(C @ X[ii]); idx.append(ii)
    P = sparse.vstack(blocks).tocsr()
    order = np.concatenate(idx)
    inv = np.empty_like(order); inv[order] = np.arange(len(order))
    return P[inv]


def log_normalise(P: sparse.csr_matrix, scale: float = 1e4) -> sparse.csr_matrix:
    tot = np.asarray(P.sum(1)).ravel()
    N = sparse.diags(scale / np.maximum(tot, 1)) @ P
    N = sparse.csr_matrix(N); N.data = np.log1p(N.data)
    return N


class RegionTransfer:
    """PCA + multinomial logistic regression on neighbourhood-pooled, log-normalised anatomy-gene counts."""

    def __init__(self, n_pcs: int = 30, n_train: int = 200_000, seed: int = 0):
        self.n_pcs, self.n_train, self.rng = n_pcs, n_train, np.random.default_rng(seed)

    def embed_fit(self, N, ref_mask):
        sel = self.rng.choice(np.flatnonzero(ref_mask), min(150_000, ref_mask.sum()), replace=False)
        self.pca = PCA(self.n_pcs, random_state=0).fit(N[sel].toarray())
        return np.vstack([self.pca.transform(N[i:i + 100_000].toarray()) for i in range(0, N.shape[0], 100_000)])

    def fit(self, Z, labels, mask):
        sel = self.rng.choice(np.flatnonzero(mask), min(self.n_train, mask.sum()), replace=False)
        self.clf = LogisticRegression(max_iter=300).fit(Z[sel], labels[sel])
        return self

    def predict(self, Z):
        p = self.clf.predict_proba(Z)
        return self.clf.classes_[p.argmax(1)], p.max(1)

    def loco_accuracy(self, Z, labels, sample, ref_samples, mask) -> pd.Series:
        acc = {}
        for s in ref_samples:
            self.fit(Z, labels, mask & (sample != s))
            te = mask & (sample == s)
            acc[s] = float((self.predict(Z[te])[0] == labels[te]).mean())
        return pd.Series(acc, name="held_out_accuracy")
