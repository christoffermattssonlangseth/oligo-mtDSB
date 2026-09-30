"""Treatment effect by brain region and time point, on sample-level pseudobulk.

This module is the single, canonical implementation of the
"condition x region x age" differential-expression comparison for the
oligo-mtDSB Xenium data set. It replaces the inline copies of
``pseudobulk_by`` / ``run_by_region`` / ``combine_delta`` that were pasted into
``notebooks-01/25``, ``notebooks-01/26`` and ``notebooks-02/04``.

Design
------
* The replicate unit is the tissue section (``sample_id``), never the cell.
  Counts are summed per (cell type, region, sample) before any test.
* Four questions are answered, each with its own tidy output table:

  1. **Simple effects** - mtDSB vs control *within* each (cell type, region,
     age) stratum. pyDESeq2 Wald test, one model per stratum.
  2. **Age contrast** - does the treatment effect differ between the two time
     points within a region?  Two variants are provided:
     a. ``run_age_contrasts``: z-test on the difference of the two stratum
        log2 fold changes (fast, uses the tables from step 1);
     b. ``run_interaction_models``: one joint pyDESeq2 model per
        (cell type, region) with design ``~ age + condition + age:condition``
        and a Wald test on the interaction coefficient (shares one dispersion
        estimate across the 12 sections and is the preferred test).
  3. **Region contrast** - is the treatment effect in a region different from
     the effect in the rest of the tissue?  ``run_region_vs_rest`` runs
     mtDSB vs control on the complementary pseudobulk (all other regions) and
     z-tests the difference of log2 fold changes.
  4. **Summaries** - counts of significant genes and an activation index per
     stratum (``summarize_strata``) and a gene x stratum matrix
     (``gene_stratum_matrix``) for heatmaps.

All tables share the same column names (``gene, baseMean, log2FC, lfcSE, stat,
pvalue, padj``) plus the stratification columns ``cell_type, region, age``.

The module only depends on numpy, pandas, scipy, statsmodels and pydeseq2
(>= 0.5, formulaic designs). AnnData is only needed for ``pseudobulk``.
"""

from __future__ import annotations

import contextlib
import io
import json
import logging
import time
import warnings
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd
from scipy import sparse, stats

logger = logging.getLogger(__name__)

# ----------------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------------

TISSUE_WIDE = "All regions"
"""Pseudo-region label used for the whole-section pseudobulk."""

STANDARD_COLS = ["gene", "baseMean", "log2FC", "lfcSE", "stat", "pvalue", "padj"]

_DESEQ_RENAME = {"log2FoldChange": "log2FC"}


@dataclass(frozen=True)
class Keys:
    """Names of the ``adata.obs`` columns and factor levels used throughout."""

    cell_type: str = "cell_class"
    region: str = "RBD_compartment_simplified"
    sample: str = "sample_id"
    age: str = "age"
    condition: str = "condition"
    case: str = "mtDSB"
    control: str = "control"
    ages: tuple[str, str] = ("21", "60")
    """(early, late) time point labels, as strings."""

    @property
    def early(self) -> str:
        return str(self.ages[0])

    @property
    def late(self) -> str:
        return str(self.ages[1])


# ----------------------------------------------------------------------------
# Pseudobulk
# ----------------------------------------------------------------------------


@dataclass
class Pseudobulk:
    """Summed counts per (cell type, region, sample).

    Attributes
    ----------
    counts
        ``DataFrame`` of shape (pseudobulk samples x genes), integer counts.
    meta
        ``DataFrame`` indexed like ``counts`` with columns
        ``cell_type, region, sample_id, age, condition, n_cells``.
    dropped
        Pseudobulk samples that were removed because they had fewer than
        ``min_cells_per_sample`` cells (same columns as ``meta``).
    """

    counts: pd.DataFrame
    meta: pd.DataFrame
    dropped: pd.DataFrame = field(default_factory=pd.DataFrame)

    def subset(self, mask: pd.Series | np.ndarray) -> "Pseudobulk":
        mask = np.asarray(mask, dtype=bool)
        return Pseudobulk(self.counts.loc[mask], self.meta.loc[mask], self.dropped)

    def stratum(self, cell_type: str, region: str, age: str | None = None) -> "Pseudobulk":
        m = (self.meta["cell_type"] == cell_type) & (self.meta["region"] == region)
        if age is not None:
            m &= self.meta["age"] == str(age)
        return self.subset(m)

    @property
    def cell_types(self) -> list[str]:
        return sorted(self.meta["cell_type"].unique())

    @property
    def regions(self) -> list[str]:
        regs = sorted(r for r in self.meta["region"].unique() if r != TISSUE_WIDE)
        if TISSUE_WIDE in set(self.meta["region"]):
            regs.append(TISSUE_WIDE)
        return regs

    def save(self, outdir: str | Path) -> None:
        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        self.counts.to_csv(outdir / "pseudobulk_counts.csv")
        self.meta.to_csv(outdir / "pseudobulk_meta.csv")
        self.dropped.to_csv(outdir / "pseudobulk_dropped.csv")

    @classmethod
    def load(cls, outdir: str | Path) -> "Pseudobulk":
        outdir = Path(outdir)
        counts = pd.read_csv(outdir / "pseudobulk_counts.csv", index_col=0)
        meta = pd.read_csv(outdir / "pseudobulk_meta.csv", index_col=0, dtype=str)
        meta["n_cells"] = meta["n_cells"].astype(int)
        dropped_path = outdir / "pseudobulk_dropped.csv"
        dropped = pd.read_csv(dropped_path, index_col=0, dtype=str) if dropped_path.exists() else pd.DataFrame()
        return cls(counts, meta, dropped)


def read_minimal(path: str | Path, layer: str | None = "counts", keys: Keys = Keys(), extra_obs: Sequence[str] = ()):
    """Read only what the pipeline needs from a (large) h5ad: obs columns + one counts matrix.

    Avoids loading ``X`` and other layers (the annotated Xenium file is > 2 GB).
    Returns an AnnData whose ``layers[layer]`` (or ``X`` if ``layer`` is None)
    holds the raw counts.
    """
    import anndata as ad
    import h5py
    from anndata.io import read_elem

    cols = [keys.cell_type, keys.region, keys.sample, keys.age, keys.condition, *extra_obs]
    with h5py.File(path, "r") as f:
        obs_group = f["obs"]
        idx_key = obs_group.attrs["_index"]
        obs = pd.DataFrame(index=pd.Index(read_elem(obs_group[idx_key]), name=None))
        for c in cols:
            if c not in obs_group:
                raise KeyError(f"obs column {c!r} not in {path}")
            obs[c] = read_elem(obs_group[c])
        var = read_elem(f["var"])
        X = read_elem(f["layers"][layer]) if layer else read_elem(f["X"])
    a = ad.AnnData(X=sparse.csr_matrix(X), obs=obs, var=var[[]])
    if layer:
        a.layers[layer] = a.X
    return a


def merge_cell_types(pb: Pseudobulk, mapping: dict[str, str]) -> Pseudobulk:
    """Collapse cell types in an existing pseudobulk (e.g. Immature OL I-IV -> one class).

    Pseudobulk counts are sums, so merging is exact for the samples that were
    kept. Samples that were dropped for having too few cells in a sub-class are
    not recovered (they are absent from ``pb.counts``); their meta rows stay in
    ``pb.dropped`` for reference. Cell types not in ``mapping`` are kept as is.
    """
    meta = pb.meta.copy()
    meta["cell_type"] = meta["cell_type"].map(lambda c: mapping.get(c, c))
    key = meta["cell_type"] + "|" + meta["region"] + "|" + meta["sample_id"]
    counts = pb.counts.groupby(key.to_numpy()).sum()
    new_meta = meta.groupby(key.to_numpy()).agg(
        cell_type=("cell_type", "first"),
        region=("region", "first"),
        sample_id=("sample_id", "first"),
        age=("age", "first"),
        condition=("condition", "first"),
        n_cells=("n_cells", "sum"),
    )
    counts = counts.loc[new_meta.index]
    counts.index.name = new_meta.index.name = "pb_id"
    return Pseudobulk(counts, new_meta, pb.dropped)


def run_all_from_pseudobulk(
    pb: Pseudobulk,
    outdir: str | Path,
    keys: Keys = Keys(),
    cell_types: Sequence[str] | None = None,
    regions: Sequence[str] | None = None,
    min_per_condition: int = 2,
    region_vs_rest: bool = True,
    interaction: bool = True,
    n_cpus: int | None = None,
    padj_thr: float = 0.05,
    lfc_thr: float = 0.5,
    min_base_mean: float = 20.0,
    run_info_extra: dict | None = None,
) -> dict[str, pd.DataFrame]:
    """All tests + summaries + CSVs starting from an existing pseudobulk."""
    t0 = time.time()
    pb.save(outdir)
    tables: dict[str, pd.DataFrame] = {}
    tables["simple_effects"], tables["simple_effects_log"] = run_simple_effects(
        pb, keys, cell_types, regions, None, min_per_condition, n_cpus
    )
    tables["age_contrast"] = run_age_contrasts(tables["simple_effects"], keys)
    if interaction:
        tables["interaction"], tables["interaction_log"] = run_interaction_models(
            pb, keys, cell_types, regions, min_per_condition, n_cpus
        )
    if region_vs_rest:
        tables["region_vs_rest"] = run_region_vs_rest(pb, tables["simple_effects"], keys, min_per_condition, n_cpus)
    tables["summary"] = summarize_strata(tables["simple_effects"], padj_thr, lfc_thr, min_base_mean)
    run_info = {
        "keys": asdict(keys),
        "min_per_condition": min_per_condition,
        "thresholds": {"padj": padj_thr, "abs_log2FC": lfc_thr, "baseMean": min_base_mean},
        "n_pseudobulk_samples": int(len(pb.meta)),
        "n_pseudobulk_dropped": int(len(pb.dropped)),
        "seconds": round(time.time() - t0, 1),
        "pydeseq2": _pkg_version("pydeseq2"),
        **(run_info_extra or {}),
    }
    save_results(outdir, tables, run_info)
    return tables


def _raw_counts_matrix(adata, layer: str | None) -> sparse.csr_matrix:
    X = adata.layers[layer] if layer else adata.X
    X = sparse.csr_matrix(X)
    sample = X.data[: min(X.nnz, 10_000)]
    if sample.size and not np.allclose(sample, np.round(sample)):
        raise ValueError(
            "Counts matrix does not look like raw integer counts "
            f"(layer={layer!r}). Point `layer` at the raw counts."
        )
    return X


def pseudobulk(
    adata,
    keys: Keys = Keys(),
    layer: str | None = "counts",
    min_cells_per_sample: int = 20,
    exclude_regions: Sequence[str] = ("Unknown",),
    include_tissue_wide: bool = True,
) -> Pseudobulk:
    """Sum raw counts per (cell type, region, section).

    Parameters
    ----------
    adata
        AnnData with raw counts in ``adata.layers[layer]`` (or ``adata.X`` if
        ``layer`` is None) and the ``Keys`` columns in ``obs``.
    keys
        Column names / factor levels.
    layer
        Layer holding raw integer counts.
    min_cells_per_sample
        Pseudobulk samples built from fewer cells are dropped (recorded in
        ``Pseudobulk.dropped``).
    exclude_regions
        Region labels to ignore entirely (e.g. ``"Unknown"``).
    include_tissue_wide
        Also build a whole-section pseudobulk per cell type with region label
        ``TISSUE_WIDE``. Needed for ``run_region_vs_rest``.
    """
    needed = [keys.cell_type, keys.region, keys.sample, keys.age, keys.condition]
    missing = [c for c in needed if c not in adata.obs]
    if missing:
        raise KeyError(f"Missing obs columns: {missing}")

    obs = adata.obs[needed].astype(str).copy()
    obs.columns = ["cell_type", "region", "sample_id", "age", "condition"]
    keep = ~obs["region"].isin([str(r) for r in exclude_regions])
    obs = obs.loc[keep]
    X = _raw_counts_matrix(adata, layer)[np.flatnonzero(keep.to_numpy())]

    # sample-level metadata must be unique per section
    per_sample = obs.groupby("sample_id")[["age", "condition"]].nunique()
    bad = per_sample[(per_sample > 1).any(axis=1)]
    if len(bad):
        raise ValueError(f"age/condition not unique within sample: {bad.index.tolist()}")

    frames = [obs]
    if include_tissue_wide:
        tw = obs.copy()
        tw["region"] = TISSUE_WIDE
        frames.append(tw)

    blocks, metas = [], []
    for frame in frames:
        group_key = frame["cell_type"] + "|" + frame["region"] + "|" + frame["sample_id"]
        codes, uniques = pd.factorize(group_key, sort=True)
        G = sparse.csr_matrix(
            (np.ones(len(codes), dtype=np.float64), (codes, np.arange(len(codes)))),
            shape=(len(uniques), len(codes)),
        )
        pb = G @ X
        blocks.append(sparse.csr_matrix(pb))
        meta = frame.groupby(group_key, sort=True).agg(
            cell_type=("cell_type", "first"),
            region=("region", "first"),
            sample_id=("sample_id", "first"),
            age=("age", "first"),
            condition=("condition", "first"),
            n_cells=("sample_id", "size"),
        )
        meta.index = uniques
        metas.append(meta)

    counts = pd.DataFrame(
        sparse.vstack(blocks).toarray().astype(np.int64),
        index=np.concatenate([m.index.to_numpy() for m in metas]),
        columns=adata.var_names,
    )
    meta = pd.concat(metas)
    counts.index.name = meta.index.name = "pb_id"

    small = meta["n_cells"] < min_cells_per_sample
    dropped = meta.loc[small].copy()
    if len(dropped):
        logger.info("dropping %d pseudobulk samples with < %d cells", len(dropped), min_cells_per_sample)
    return Pseudobulk(counts.loc[~small], meta.loc[~small], dropped)


# ----------------------------------------------------------------------------
# pyDESeq2 wrapper
# ----------------------------------------------------------------------------

_INTERNAL_COND = {"control": "A_control", "case": "B_case"}
_INTERNAL_AGE = {"early": "A_early", "late": "B_late"}


def _internal_metadata(meta: pd.DataFrame, keys: Keys) -> pd.DataFrame:
    """Recode factor levels so that the alphabetical reference is control / early."""
    md = pd.DataFrame(index=meta.index)
    md["cond"] = np.where(meta["condition"] == keys.case, _INTERNAL_COND["case"], _INTERNAL_COND["control"])
    md["age"] = np.where(meta["age"].astype(str) == keys.late, _INTERNAL_AGE["late"], _INTERNAL_AGE["early"])
    return md


def _deseq_fit(counts: pd.DataFrame, meta: pd.DataFrame, design: str, n_cpus: int | None, min_total_count: int):
    from pydeseq2.dds import DeseqDataSet

    counts = counts.loc[:, counts.sum(axis=0) >= min_total_count]
    with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
        # pyDESeq2 warns when the parametric dispersion trend fails and it falls
        # back to the mean trend; that is expected for small strata.
        warnings.simplefilter("ignore")
        dds = DeseqDataSet(counts=counts, metadata=meta, design=design, quiet=True, n_cpus=n_cpus)
        dds.deseq2()
    return dds


def _deseq_test(dds, contrast, n_cpus: int | None) -> pd.DataFrame:
    from pydeseq2.ds import DeseqStats

    with contextlib.redirect_stdout(io.StringIO()):
        st = DeseqStats(dds, contrast=contrast, quiet=True, n_cpus=n_cpus)
        st.summary()
    res = st.results_df.rename(columns=_DESEQ_RENAME)
    res.index.name = "gene"
    return res.reset_index()[STANDARD_COLS]


def _design_vector(dds, terms: Iterable[str]) -> np.ndarray:
    """Numeric contrast: 1 for each design-matrix column whose name contains a term."""
    cols = list(dds.obsm["design_matrix"].columns)
    vec = np.zeros(len(cols))
    for term in terms:
        hits = [i for i, c in enumerate(cols) if term == c]
        if not hits:
            raise KeyError(f"design column {term!r} not in {cols}")
        vec[hits[0]] = 1.0
    return vec


def _cond_term() -> str:
    return f"cond[T.{_INTERNAL_COND['case']}]"


def _age_term() -> str:
    return f"age[T.{_INTERNAL_AGE['late']}]"


def _interaction_term() -> str:
    return f"{_age_term()}:{_cond_term()}"


def de_condition(
    pb: Pseudobulk,
    keys: Keys = Keys(),
    min_per_condition: int = 2,
    n_cpus: int | None = None,
    min_total_count: int = 10,
) -> pd.DataFrame | None:
    """mtDSB vs control on one pseudobulk subset (already restricted to a stratum).

    Returns ``None`` when either condition has fewer than ``min_per_condition``
    sections.
    """
    n = pb.meta["condition"].value_counts()
    if n.get(keys.case, 0) < min_per_condition or n.get(keys.control, 0) < min_per_condition:
        return None
    md = _internal_metadata(pb.meta, keys)
    dds = _deseq_fit(pb.counts, md[["cond"]], "~ cond", n_cpus, min_total_count)
    res = _deseq_test(dds, ["cond", _INTERNAL_COND["case"], _INTERNAL_COND["control"]], n_cpus)
    res["n_case"] = int(n.get(keys.case, 0))
    res["n_control"] = int(n.get(keys.control, 0))
    return res


def run_simple_effects(
    pb: Pseudobulk,
    keys: Keys = Keys(),
    cell_types: Sequence[str] | None = None,
    regions: Sequence[str] | None = None,
    ages: Sequence[str] | None = None,
    min_per_condition: int = 2,
    n_cpus: int | None = None,
    verbose: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """mtDSB vs control within every (cell type, region, age) stratum.

    Returns
    -------
    results
        Long table, one row per gene per stratum, columns
        ``cell_type, region, age`` + ``STANDARD_COLS`` + ``n_case, n_control``.
    log
        One row per stratum with ``status`` (``ok`` or a skip reason) and the
        section counts per condition.
    """
    cell_types = list(cell_types or pb.cell_types)
    regions = list(regions or pb.regions)
    ages = [str(a) for a in (ages or keys.ages)]

    out, log = [], []
    for ct in cell_types:
        for reg in regions:
            for age in ages:
                sub = pb.stratum(ct, reg, age)
                n = sub.meta["condition"].value_counts()
                row = {
                    "cell_type": ct,
                    "region": reg,
                    "age": age,
                    "n_case": int(n.get(keys.case, 0)),
                    "n_control": int(n.get(keys.control, 0)),
                    "n_cells": int(sub.meta["n_cells"].sum()),
                }
                t0 = time.time()
                res = de_condition(sub, keys, min_per_condition, n_cpus) if len(sub.meta) else None
                if res is None:
                    row["status"] = f"skipped: <{min_per_condition} sections per condition"
                    row["n_genes_tested"] = 0
                else:
                    res.insert(0, "age", age)
                    res.insert(0, "region", reg)
                    res.insert(0, "cell_type", ct)
                    out.append(res)
                    row["status"] = "ok"
                    row["n_genes_tested"] = int(res["pvalue"].notna().sum())
                row["seconds"] = round(time.time() - t0, 1)
                log.append(row)
                if verbose:
                    logger.info("%s | %s | %s -> %s", ct, reg, age, row["status"])
    results = pd.concat(out, ignore_index=True) if out else pd.DataFrame(columns=["cell_type", "region", "age", *STANDARD_COLS])
    return results, pd.DataFrame(log)


# ----------------------------------------------------------------------------
# Contrasts between strata
# ----------------------------------------------------------------------------


def delta_lfc_test(res_a: pd.DataFrame, res_b: pd.DataFrame, suffix_a: str, suffix_b: str) -> pd.DataFrame:
    """z-test on the difference of two independent log2 fold changes (B - A).

    Both inputs are standard result tables (``gene, log2FC, lfcSE, padj`` ...).
    Genes must be tested in both (non-missing ``lfcSE``). The Benjamini-Hochberg
    adjustment is over the genes that pass that filter.
    """
    from statsmodels.stats.multitest import multipletests

    a = res_a.set_index("gene")[["baseMean", "log2FC", "lfcSE", "padj"]].add_suffix(f"_{suffix_a}")
    b = res_b.set_index("gene")[["baseMean", "log2FC", "lfcSE", "padj"]].add_suffix(f"_{suffix_b}")
    m = a.join(b, how="inner")
    m = m.dropna(subset=[f"lfcSE_{suffix_a}", f"lfcSE_{suffix_b}"])
    m["delta_log2FC"] = m[f"log2FC_{suffix_b}"] - m[f"log2FC_{suffix_a}"]
    m["se_delta"] = np.sqrt(m[f"lfcSE_{suffix_a}"] ** 2 + m[f"lfcSE_{suffix_b}"] ** 2)
    m["stat"] = m["delta_log2FC"] / m["se_delta"]
    m["pvalue"] = 2 * stats.norm.sf(np.abs(m["stat"]))
    m["padj"] = np.nan
    ok = m["pvalue"].notna()
    if ok.any():
        m.loc[ok, "padj"] = multipletests(m.loc[ok, "pvalue"], method="fdr_bh")[1]
    m.index.name = "gene"
    return m.reset_index()


def run_age_contrasts(results: pd.DataFrame, keys: Keys = Keys()) -> pd.DataFrame:
    """Late-minus-early difference in treatment effect, per (cell type, region).

    Uses the simple-effect tables from :func:`run_simple_effects`.
    """
    out = []
    for (ct, reg), grp in results.groupby(["cell_type", "region"], sort=True):
        a = grp[grp["age"] == keys.early]
        b = grp[grp["age"] == keys.late]
        if a.empty or b.empty:
            continue
        d = delta_lfc_test(a, b, f"age{keys.early}", f"age{keys.late}")
        d.insert(0, "region", reg)
        d.insert(0, "cell_type", ct)
        out.append(d)
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


def run_interaction_models(
    pb: Pseudobulk,
    keys: Keys = Keys(),
    cell_types: Sequence[str] | None = None,
    regions: Sequence[str] | None = None,
    min_per_group: int = 2,
    n_cpus: int | None = None,
    min_total_count: int = 10,
    verbose: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Joint model ``~ age + condition + age:condition`` per (cell type, region).

    For each model three Wald tests are reported, distinguished by ``term``:

    * ``condition_at_<early>`` - treatment effect at the early time point
    * ``condition_at_<late>`` - treatment effect at the late time point
      (main effect + interaction)
    * ``age_x_condition`` - difference of the two (the interaction). This is
      the formal test of "does the treatment effect change with time".

    Every age x condition cell needs at least ``min_per_group`` sections.
    """
    cell_types = list(cell_types or pb.cell_types)
    regions = list(regions or pb.regions)
    out, log = [], []
    for ct in cell_types:
        for reg in regions:
            sub = pb.stratum(ct, reg)
            sub = sub.subset(sub.meta["age"].isin([keys.early, keys.late]))
            cells = sub.meta.groupby(["age", "condition"]).size()
            row = {"cell_type": ct, "region": reg, "n_sections": int(len(sub.meta))}
            for age in (keys.early, keys.late):
                for cond in (keys.control, keys.case):
                    row[f"n_{cond}_{age}"] = int(cells.get((age, cond), 0))
            complete = all(cells.get((a, c), 0) >= min_per_group for a in (keys.early, keys.late) for c in (keys.control, keys.case))
            t0 = time.time()
            if not complete:
                row["status"] = f"skipped: <{min_per_group} sections in some age x condition cell"
                log.append(row)
                continue
            md = _internal_metadata(sub.meta, keys)
            dds = _deseq_fit(sub.counts, md, "~ age + cond + age:cond", n_cpus, min_total_count)
            tests = {
                f"condition_at_{keys.early}": _design_vector(dds, [_cond_term()]),
                f"condition_at_{keys.late}": _design_vector(dds, [_cond_term(), _interaction_term()]),
                "age_x_condition": _design_vector(dds, [_interaction_term()]),
            }
            for term, vec in tests.items():
                res = _deseq_test(dds, vec, n_cpus)
                res.insert(0, "term", term)
                res.insert(0, "region", reg)
                res.insert(0, "cell_type", ct)
                out.append(res)
            row["status"] = "ok"
            row["seconds"] = round(time.time() - t0, 1)
            log.append(row)
            if verbose:
                logger.info("interaction %s | %s -> ok", ct, reg)
    results = pd.concat(out, ignore_index=True) if out else pd.DataFrame()
    return results, pd.DataFrame(log)


def run_region_vs_rest(
    pb: Pseudobulk,
    results: pd.DataFrame,
    keys: Keys = Keys(),
    min_per_condition: int = 2,
    n_cpus: int | None = None,
    verbose: bool = True,
) -> pd.DataFrame:
    """Treatment effect in a region minus the effect in all other regions.

    Needs ``pseudobulk(..., include_tissue_wide=True)``. The "rest" pseudobulk
    for a region is the whole-section pseudobulk minus that region's counts,
    per (cell type, section), so no re-aggregation over cells is needed.
    """
    if TISSUE_WIDE not in set(pb.meta["region"]):
        raise ValueError("run_region_vs_rest needs the tissue-wide pseudobulk (include_tissue_wide=True)")
    out = []
    strata = results.loc[results["region"] != TISSUE_WIDE, ["cell_type", "region", "age"]].drop_duplicates()
    for ct, reg, age in strata.itertuples(index=False):
        region_pb = pb.stratum(ct, reg, age)
        all_pb = pb.stratum(ct, TISSUE_WIDE, age)
        joint = all_pb.meta[["sample_id"]].reset_index().merge(
            region_pb.meta[["sample_id"]].reset_index(), on="sample_id", suffixes=("_all", "_reg")
        )
        if joint.empty:
            continue
        rest_counts = all_pb.counts.loc[joint["pb_id_all"]].to_numpy() - region_pb.counts.loc[joint["pb_id_reg"]].to_numpy()
        rest_counts = pd.DataFrame(np.clip(rest_counts, 0, None), index=joint["pb_id_all"].to_numpy(), columns=pb.counts.columns)
        rest_meta = all_pb.meta.loc[joint["pb_id_all"]].copy()
        rest_meta["region"] = f"rest of tissue ({reg})"
        rest_meta["n_cells"] = rest_meta["n_cells"].to_numpy() - region_pb.meta.loc[joint["pb_id_reg"], "n_cells"].to_numpy()
        rest = Pseudobulk(rest_counts, rest_meta)
        res_rest = de_condition(rest, keys, min_per_condition, n_cpus)
        if res_rest is None:
            continue
        res_reg = results[(results["cell_type"] == ct) & (results["region"] == reg) & (results["age"] == age)]
        d = delta_lfc_test(res_rest, res_reg, "rest", "region")
        d.insert(0, "age", age)
        d.insert(0, "region", reg)
        d.insert(0, "cell_type", ct)
        out.append(d)
        if verbose:
            logger.info("region vs rest %s | %s | %s -> ok", ct, reg, age)
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


# ----------------------------------------------------------------------------
# Summaries
# ----------------------------------------------------------------------------


def significant(df: pd.DataFrame, padj_thr: float = 0.05, lfc_thr: float = 0.5, min_base_mean: float = 20.0) -> pd.Series:
    """Boolean mask: FDR, absolute effect and expression thresholds all met."""
    return (df["padj"] < padj_thr) & (df["log2FC"].abs() >= lfc_thr) & (df["baseMean"] >= min_base_mean)


def summarize_strata(
    results: pd.DataFrame,
    padj_thr: float = 0.05,
    lfc_thr: float = 0.5,
    min_base_mean: float = 20.0,
    top_n: int = 5,
) -> pd.DataFrame:
    """One row per (cell type, region, age) with hit counts and an activation index.

    ``activation_index`` = number of significant up-regulated genes x their mean
    log2FC, the same quantity as the manuscript's "regional activation
    fingerprint" but computed on adjusted p-values and after the expression
    filter.
    """
    df = results.copy()
    df["sig"] = significant(df, padj_thr, lfc_thr, min_base_mean)
    rows = []
    for (ct, reg, age), g in df.groupby(["cell_type", "region", "age"], sort=True):
        s = g[g["sig"]]
        up = s[s["log2FC"] > 0].sort_values("padj")
        down = s[s["log2FC"] < 0].sort_values("padj")
        rows.append(
            {
                "cell_type": ct,
                "region": reg,
                "age": age,
                "n_tested": int(g["padj"].notna().sum()),
                "n_sig": int(len(s)),
                "n_up": int(len(up)),
                "n_down": int(len(down)),
                "mean_abs_log2FC_sig": float(s["log2FC"].abs().mean()) if len(s) else 0.0,
                "activation_index": float(len(up) * up["log2FC"].mean()) if len(up) else 0.0,
                "top_up": ", ".join(up["gene"].head(top_n)),
                "top_down": ", ".join(down["gene"].head(top_n)),
            }
        )
    return pd.DataFrame(rows)


def gene_stratum_matrix(
    results: pd.DataFrame,
    cell_type: str,
    genes: Sequence[str] | None = None,
    value: str = "log2FC",
    mask_padj: float | None = None,
    top_n: int = 30,
    padj_thr: float = 0.05,
    lfc_thr: float = 0.5,
    min_base_mean: float = 20.0,
    regions: Sequence[str] | None = None,
    keys: Keys = Keys(),
) -> pd.DataFrame:
    """Genes x "region @ age" matrix of ``value`` for one cell type.

    If ``genes`` is None, the union of the ``top_n`` most significant genes per
    stratum is used. With ``mask_padj`` set, entries with ``padj >= mask_padj``
    become NaN so only significant cells are coloured.
    """
    df = results[results["cell_type"] == cell_type].copy()
    if regions is not None:
        df = df[df["region"].isin(regions)]
    if genes is None:
        df["sig"] = significant(df, padj_thr, lfc_thr, min_base_mean)
        genes = (
            df[df["sig"]]
            .sort_values("padj")
            .groupby(["region", "age"])
            .head(top_n)["gene"]
            .drop_duplicates()
            .tolist()
        )
    df = df[df["gene"].isin(genes)]
    if mask_padj is not None:
        df.loc[df["padj"] >= mask_padj, value] = np.nan
    df["stratum"] = df["region"] + " @ " + df["age"].astype(str)
    order = [f"{r} @ {a}" for r in (regions or _region_order(df["region"])) for a in keys.ages]
    mat = df.pivot_table(index="gene", columns="stratum", values=value, aggfunc="first")
    mat = mat.reindex(index=[g for g in genes if g in mat.index], columns=[c for c in order if c in mat.columns])
    return mat


def _region_order(regions: Iterable[str]) -> list[str]:
    regs = sorted(set(regions) - {TISSUE_WIDE})
    if TISSUE_WIDE in set(regions):
        regs.append(TISSUE_WIDE)
    return regs


# ----------------------------------------------------------------------------
# Persistence
# ----------------------------------------------------------------------------

FILES = {
    "simple_effects": "de_condition_by_region_age.csv",
    "simple_effects_log": "de_condition_by_region_age_log.csv",
    "age_contrast": "de_age_contrast_by_region.csv",
    "interaction": "de_interaction_by_region.csv",
    "interaction_log": "de_interaction_by_region_log.csv",
    "region_vs_rest": "de_region_vs_rest.csv",
    "summary": "summary_by_region_age.csv",
}


def save_results(outdir: str | Path, tables: dict[str, pd.DataFrame], run_info: dict | None = None) -> None:
    """Write every table in ``tables`` (keys as in ``FILES``) plus ``run_info.json``."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    for key, df in tables.items():
        if df is None or (hasattr(df, "empty") and df.empty):
            continue
        df.to_csv(outdir / FILES[key], index=False)
    if run_info is not None:
        (outdir / "run_info.json").write_text(json.dumps(run_info, indent=2, default=str))


def load_results(outdir: str | Path) -> dict[str, pd.DataFrame]:
    """Read back whatever tables exist in ``outdir``; ``age`` is kept as string."""
    outdir = Path(outdir)
    tables = {}
    for key, name in FILES.items():
        p = outdir / name
        if p.exists():
            tables[key] = pd.read_csv(p, dtype={"age": str})
    return tables


def run_all(
    adata,
    outdir: str | Path,
    keys: Keys = Keys(),
    layer: str | None = "counts",
    min_cells_per_sample: int = 20,
    exclude_regions: Sequence[str] = ("Unknown",),
    cell_types: Sequence[str] | None = None,
    regions: Sequence[str] | None = None,
    min_per_condition: int = 2,
    region_vs_rest: bool = True,
    interaction: bool = True,
    n_cpus: int | None = None,
    padj_thr: float = 0.05,
    lfc_thr: float = 0.5,
    min_base_mean: float = 20.0,
) -> dict[str, pd.DataFrame]:
    """End-to-end: pseudobulk -> all tests -> summaries -> CSVs in ``outdir``."""
    t0 = time.time()
    pb = pseudobulk(adata, keys, layer, min_cells_per_sample, exclude_regions, include_tissue_wide=True)
    pb.save(outdir)
    tables: dict[str, pd.DataFrame] = {}
    tables["simple_effects"], tables["simple_effects_log"] = run_simple_effects(
        pb, keys, cell_types, regions, None, min_per_condition, n_cpus
    )
    tables["age_contrast"] = run_age_contrasts(tables["simple_effects"], keys)
    if interaction:
        tables["interaction"], tables["interaction_log"] = run_interaction_models(
            pb, keys, cell_types, regions, min_per_condition, n_cpus
        )
    if region_vs_rest:
        tables["region_vs_rest"] = run_region_vs_rest(pb, tables["simple_effects"], keys, min_per_condition, n_cpus)
    tables["summary"] = summarize_strata(tables["simple_effects"], padj_thr, lfc_thr, min_base_mean)
    run_info = {
        "keys": asdict(keys),
        "layer": layer,
        "min_cells_per_sample": min_cells_per_sample,
        "exclude_regions": list(exclude_regions),
        "min_per_condition": min_per_condition,
        "thresholds": {"padj": padj_thr, "abs_log2FC": lfc_thr, "baseMean": min_base_mean},
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_pseudobulk_samples": int(len(pb.meta)),
        "n_pseudobulk_dropped": int(len(pb.dropped)),
        "seconds": round(time.time() - t0, 1),
        "pydeseq2": _pkg_version("pydeseq2"),
    }
    save_results(outdir, tables, run_info)
    return tables


def _pkg_version(name: str) -> str | None:
    try:
        from importlib.metadata import version

        return version(name)
    except Exception:  # pragma: no cover
        return None
