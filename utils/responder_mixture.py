"""Probability that each mtDSB mouse responds: a two-component Bayesian mixture, anchored on the controls.

The hard responder call (``utils.ol_response.responder_calls``: whole-section ISR above the highest control) gives
yes / no. This module replaces it with a posterior probability per mouse that accounts for how few controls define
the "normal" range.

Model, for per-mouse readouts ``x`` (d dimensions, already on a roughly normal scale, e.g. log ISR CPM):

* cohort (= Xenium slide) baseline removed by subtracting the mean of that cohort's controls;
* non-responders are drawn like controls. Their predictive density is a multivariate t with the pooled control
  covariance ``S`` (controls minus their cohort mean, ``df = n_controls - n_cohorts``), inflated by
  ``1 + 1 / n_controls_in_cohort`` for the uncertainty in the cohort mean;
* responders are shifted by ``delta`` (a vector, one entry per readout) with covariance ``k^2 S``; by default
  (``one_sided``) they must lie above their cohort's control mean in the first readout (the response raises ISR);
* each mtDSB mouse responds with probability ``pi``.

Priors: ``pi ~ Beta(a, b)`` (default 1, 1), ``delta`` uniform on a grid from ``delta_lo`` to ``delta_hi`` control SDs
in every readout, ``k`` uniform on ``k_grid``. The posterior is computed exactly: all 2^n label vectors are
enumerated, ``pi`` is integrated out (beta-binomial), ``delta`` and ``k`` are summed over the grid. With six mtDSB
mice this is 64 label vectors, so no sampling is needed.

``control_check`` is the null: each control in turn is treated as an unlabelled mouse (and removed from the control
reference) and gets the same posterior; it should be low.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import stats
from scipy.special import betaln, gammaln, logsumexp


@dataclass
class Prior:
    a: float = 1.0
    b: float = 1.0
    delta_lo: float = 0.0
    delta_hi: float = 12.0
    n_delta: int = 49
    k_grid: tuple = (1.0, 1.5, 2.0, 3.0)
    one_sided: bool = True  # responders only above the cohort control mean in the first readout


def centre_by_cohort(x: pd.DataFrame, cohort: pd.Series, is_control: pd.Series) -> pd.DataFrame:
    """Subtract each cohort's control mean from every mouse of that cohort."""
    base = x[is_control].groupby(cohort[is_control]).mean()
    return x - base.loc[cohort].to_numpy()


def _logpdf_t(x: np.ndarray, S: np.ndarray, df: float, c=1.0) -> np.ndarray:
    """Log density of a zero-centred multivariate t with shape ``c * S``; ``x`` is (..., d), ``c`` broadcasts."""
    d = x.shape[-1]
    L = np.linalg.cholesky(S)
    u = np.linalg.solve(L, np.moveaxis(x, -1, 0).reshape(d, -1)).reshape((d,) + x.shape[:-1])
    maha = (u ** 2).sum(0) / c
    return (gammaln((df + d) / 2) - gammaln(df / 2) - d / 2 * np.log(df * np.pi) - np.log(np.diag(L)).sum()
            - d / 2 * np.log(c) - (df + d) / 2 * np.log1p(maha / df))


def _log_resp(x: np.ndarray, G_d: np.ndarray, G_k: np.ndarray, S: np.ndarray, df: float, infl, one_sided: bool):
    """Log density of mice ``x`` (n, d) under the responder component at each grid point -> (G, n).

    With ``one_sided`` the density is truncated to x[:, 0] > 0 (and renormalised), so a mouse *below* the controls
    cannot be explained as a responder with a wide spread.
    """
    c = G_k[:, None] ** 2 * np.broadcast_to(infl, (x.shape[0],))[None, :]
    l1 = _logpdf_t(x[None, :, :] - G_d[:, None, :], S, df, c)
    if one_sided:
        mass = stats.t.cdf(G_d[:, :1] / np.sqrt(S[0, 0] * c), df)
        l1 = np.where(x[None, :, 0] > 0, l1 - np.log(mass), -np.inf)
    return l1


def fit(x: pd.DataFrame, cohort: pd.Series, is_control: pd.Series, prior: Prior = Prior()) -> dict:
    """Posterior P(responder) for every non-control mouse, posterior of pi and of the responder shift.

    ``x``: mice x readouts. ``cohort``: cohort (slide) per mouse. ``is_control``: boolean per mouse.
    """
    x = x.astype(float)
    z = centre_by_cohort(x, cohort, is_control)
    ctrl, test = z[is_control], z[~is_control]
    n_ctrl, n_coh = len(ctrl), cohort[is_control].nunique()
    df = n_ctrl - n_coh
    S = np.atleast_2d(np.cov(ctrl.to_numpy().T, ddof=n_coh))
    sd = np.sqrt(np.diag(S))
    n_in_coh = cohort[is_control].value_counts()
    infl = (1 + 1 / n_in_coh.loc[cohort[~is_control]]).to_numpy()
    T = test.to_numpy()
    d, n = T.shape[1], T.shape[0]

    # log density of each test mouse under the non-responder component
    l0 = _logpdf_t(T, S, df, infl)
    # under the responder component, for every grid point (delta vector, k)
    axis = np.linspace(prior.delta_lo, prior.delta_hi, prior.n_delta)
    deltas = np.array(list(itertools.product(axis, repeat=d))) * sd
    G_d = np.repeat(deltas, len(prior.k_grid), axis=0)  # (G, d)
    G_k = np.tile(np.asarray(prior.k_grid, float), len(deltas))  # (G,)
    l1 = _log_resp(T, G_d, G_k, S, df, infl, prior.one_sided)  # (G, n)
    log_grid_prior = -np.log(len(G_k))

    labels = np.array(list(itertools.product([0, 1], repeat=n)), dtype=bool)
    m = labels.sum(1)
    log_p_labels = betaln(prior.a + m, prior.b + n - m) - betaln(prior.a, prior.b)
    # log p(data | labels, grid point): responders use l1[g], others l0
    ll0 = (l0[None, :] * ~labels).sum(1)  # (L,)
    ll1 = np.where(labels[:, None, :], l1[None, :, :], 0).sum(2)  # (L, G); -inf-safe
    joint = log_p_labels[:, None] + ll0[:, None] + ll1 + log_grid_prior  # (L, G)
    log_ev = logsumexp(joint)
    post_lg = np.exp(joint - log_ev)
    post_labels = post_lg.sum(1)
    p_resp = pd.Series(post_labels @ labels, index=test.index, name="p_responder")

    # posterior of pi on a grid (mixture of betas over label vectors)
    pis = (np.arange(200) + 0.5) / 200  # bin midpoints: Beta(<1, .) densities are infinite at 0 and 1
    pi_dens = sum(w * stats.beta(prior.a + mm, prior.b + n - mm).pdf(pis) for w, mm in zip(post_labels, m))
    # posterior of the shift (in control SDs) per readout
    post_g = post_lg.sum(0)
    delta_sd = deltas / sd
    shift = {c: pd.Series(post_g @ (np.repeat(delta_sd[:, j], len(prior.k_grid))[:, None] == axis[None, :]),
                          index=axis) for j, c in enumerate(x.columns)}
    e_pi = float(sum(w * (prior.a + mm) / (prior.a + prior.b + n) for w, mm in zip(post_labels, m)))
    top = np.argsort(post_labels)[::-1][:5]
    return {"p_responder": p_resp, "z": z, "control_sd": pd.Series(sd, index=x.columns), "df": df,
            "pi_grid": pis, "pi_density": pi_dens, "e_pi": e_pi, "shift_posterior": shift,
            "top_labelings": pd.DataFrame({"labels": ["".join("R" if v else "-" for v in labels[t]) for t in top],
                                           "posterior": post_labels[top]}),
            "log_evidence": float(log_ev), "_grid": (G_d, G_k, post_g, S, prior.one_sided)}


def control_check(x: pd.DataFrame, cohort: pd.Series, is_control: pd.Series, prior: Prior = Prior()) -> pd.Series:
    """Leave-one-control-out: P(responder) for each control when it is treated as an unlabelled mouse."""
    out = {}
    for c in x.index[is_control]:
        ic = is_control.copy()
        ic[c] = False
        out[c] = fit(x, cohort, ic, prior)["p_responder"][c]
    return pd.Series(out, name="p_responder")


def predictive_p(res: dict, new_z: np.ndarray, n_ctrl_in_cohort: int = 3) -> np.ndarray:
    """P(responder) for new mice with cohort-centred readouts ``new_z`` (n x d), given the fitted posterior.

    Uses the posterior mean of pi and the posterior over the responder shift and spread; this is the rule to apply to
    the next cohort (with its own controls for centring).
    """
    G_d, G_k, post_g, S, one_sided = res["_grid"]
    new_z = np.atleast_2d(new_z)
    infl = 1 + 1 / n_ctrl_in_cohort
    df = res["df"]
    l0 = _logpdf_t(new_z, S, df, infl)
    l1 = logsumexp(_log_resp(new_z, G_d, G_k, S, df, infl, one_sided), axis=0, b=post_g[:, None])
    pi = res["e_pi"]
    return 1 / (1 + np.exp(np.log1p(-pi) + l0 - np.log(pi) - l1))
