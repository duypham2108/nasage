"""Analysis tools on your AnnData (scanpy-style).

Results are stored in ``adata.uns["nasage"]`` and the main table of
each step is also returned. Order: ``score_modules`` -> ``donor_summary`` ->
``age_effects`` (then ``composition``, ``permutation_null``, ``module_correlation``,
``sensitivity`` as needed).
"""
from __future__ import annotations

import json
import warnings

import pandas as pd

from . import aggregate as _agg
from . import modules as _mods
from . import scoring as _scoring
from . import stats as _stats
from ._design import canonical_obs, check_counts, donor_table, store


def _need(adata, key, step):
    st = adata.uns.get("nasage", {})
    if key not in st:
        raise RuntimeError(f"run nas.tl.{step}(adata, ...) first")
    return st[key]


def _feature_obs_cols(score_cols):
    cols = []
    for c in score_cols:
        if c.endswith("_logit"):               # ratio: computed per donor from pooled counts
            base = c[:-len("_logit")]
            cols += [base + "_num", base + "_den"]
            continue
        cols.append(c)
        if c.endswith("_logexpr"):
            cols.append(c[:-len("_logexpr")] + "_cpm")
    return cols


def _clean(df: pd.DataFrame) -> pd.DataFrame:
    """Make a table writable to .h5ad: object columns become plain strings ('' for missing)."""
    df = df.copy()
    for c in df.columns:
        if df[c].dtype == object:
            df[c] = df[c].map(lambda v: "" if v is None or (isinstance(v, float) and v != v)
                              else str(v))
    return df


def _params(adata) -> dict:
    p = dict(_need(adata, "params", "donor_summary"))
    if "age_bins" in p:                        # stored as JSON (mixed int/str lists)
        p["age_bins"] = [tuple(b) for b in json.loads(p["age_bins"])]
    return p


def _obs(adata):
    cols = _need(adata, "score_columns", "score_modules")
    return canonical_obs(adata, _params(adata), extra_cols=_feature_obs_cols(cols))


# ----------------------------------------------------------------------------- scoring
def score_modules(adata, *, layer=None, modules=None, max_rank=_scoring.MAX_RANK, n_jobs=1,
                  resolve_symbols_online=False) -> list[str]:
    """Score every module in every cell; adds columns to ``adata.obs``.

    layer : layer with raw integer counts (None = adata.X). UCell ranks genes
        within a cell, so no normalisation is needed.
    modules : gene modules (default: the built-in 15; e.g.
        ``nas.drop_hla(nas.load_modules())`` for the HLA sensitivity run).

    Returns the score column names. Stores the symbol-mapping report and which
    module genes were found.
    """
    check_counts(adata, layer)
    M = modules if modules is not None else _mods.load_modules()
    mapped, symbol_report = _mods.harmonise_symbols(adata.var_names, M,
                                                    online=resolve_symbols_online)
    cols = _scoring.score_modules(adata, mapped, layer=layer, max_rank=max_rank, n_jobs=n_jobs)
    st = store(adata)
    st["layer"] = "X" if layer is None else layer
    st["modules_mapped"] = mapped
    st["symbol_report"] = _clean(symbol_report)
    st["module_presence"] = _clean(_mods.module_report(mapped, M))
    st["score_columns"] = list(cols)
    return cols


# ------------------------------------------------------------------- donor summaries
def donor_summary(adata, *, donor_key, celltype_key, age_key, age_bins=None, age_order=None,
                  batch_key=None, sex_key=None, subtype_key=None, min_cells=10,
                  statistic="auto", applicability=None) -> pd.DataFrame:
    """One row per donor x cell type: the table every test is based on.

    donor_key, celltype_key, age_key : your ``adata.obs`` columns. Cell-type
        values must match the applicability matrix (basal, ciliated, secretory,
        dendritic, macrophage, Tcell, Bcell, plasma, neutrophil, squamous,
        ionocyte, tuft, mast, NK); others are composition only. Cells with NaN
        cell type are excluded.
    age_bins : ``[(lo, hi, label), ...]`` if age_key is numeric; leave None if it
        already holds age groups (order: `age_order`, else the column's
        categorical order, else number-aware sort).
    batch_key, sex_key : adjusted for in every model.
    subtype_key : your finer labels (goblet/club, cDC1/mregDC, CD4/CD8...). Each
        donor's subtype fractions within a population enter the model, so a
        shift in subtype mix is not read as a change in cell state.
    min_cells : below this a donor x cell type is MISSING, never zero.
    statistic : "auto" (median, or mean where most cells score zero), "median",
        "mean" or "winsorized_mean".

    Also computes module coverage (poorly measured pairs are excluded from
    testing) and the score distribution (ceiling / mostly-zero flags).
    """
    cols = _need(adata, "score_columns", "score_modules")
    params = dict(donor_key=donor_key, celltype_key=celltype_key, age_key=age_key,
                  age_bins=json.dumps([list(b) for b in age_bins]) if age_bins else None,
                  age_order=list(age_order) if age_order else None, batch_key=batch_key,
                  sex_key=sex_key, subtype_key=subtype_key, min_cells=min_cells,
                  statistic=statistic)
    params = {k: v for k, v in params.items() if v is not None}     # h5ad-friendly
    st = store(adata)
    # results computed from an earlier summary would no longer match: drop them
    for k in ("omnibus", "contrasts", "composition", "proportions", "permutation",
              "permutation_rates", "sensitivity", "null_control", "null_tests"):
        st.pop(k, None)
    st["params"] = params
    A = applicability if applicability is not None else _mods.load_applicability()
    obs = _obs(adata)
    unknown = sorted(set(obs.cell_type.astype(str)) - {ct for d in A.values() for ct in d})
    if unknown:
        warnings.warn(f"cell types not in the applicability matrix (composition only): {unknown}")
    layer = None if st["layer"] == "X" else st["layer"]
    st["design"] = _clean(donor_table(obs))
    st["applicability"] = A
    st["coverage"] = _clean(_scoring.module_coverage(adata, st["modules_mapped"], A,
                                                     obs.cell_type, layer=layer))
    st["score_distribution"] = _scoring.score_distribution(obs, cols)
    summary = _agg.donor_celltype_summary(obs, cols, min_cells=min_cells, statistic=statistic)
    st["aggregation"] = _clean(pd.DataFrame(summary.attrs["aggregation"]))
    st["subtype_covariates"] = summary.attrs["subtype_covariates"]
    st["donor_summary"] = _clean(summary)
    return st["donor_summary"]


# ---------------------------------------------------------------------- age models
def age_effects(adata, *, reference, adjust_depth=True, adjust_subtypes=True) -> pd.DataFrame:
    """Age-group differences per (module, cell type), donor as the unit.

    Model: score ~ age group + batch + sex + log_genes (depth) + subtype
    fractions; the optional terms are dropped when there are too few donors
    (recorded in ``covariates`` / ``dropped_covariates``). Primary and
    exploratory pairs are BH-adjusted as separate families; only primary rows
    with ``p_adj_BH < 0.05`` are findings.

    Returns the omnibus table; contrasts vs `reference` (estimate, 95% CI,
    standardised effect) are in ``adata.uns["nasage"]["contrasts"]``.
    """
    summary = _need(adata, "donor_summary", "donor_summary")
    st = store(adata)
    cols, A, cov = st["score_columns"], st["applicability"], st["coverage"]
    kw = dict(adjust_depth=adjust_depth, adjust_subtypes=adjust_subtypes)
    primary = _stats.build_test_family(A, cov, cols, tier=1.0)
    explor = _stats.build_test_family(A, cov, cols, tier=0.5)
    om_p, ct_p = _stats.fit_age_models(summary, primary, reference, **kw)
    om_e, ct_e = _stats.fit_age_models(summary, explor, reference, **kw)
    omnibus = pd.concat([om_p, om_e], ignore_index=True)
    contrasts = pd.concat([ct_p, ct_e], ignore_index=True)
    omnibus = _stats.adjust_family(omnibus, "p_omnibus",
                                   (omnibus.tier == "primary") & omnibus.estimable)
    if len(contrasts):
        contrasts = _stats.adjust_family(contrasts, "p_value", contrasts.tier == "primary")
    st["params"]["reference"] = reference
    st["params"]["adjust_depth"] = adjust_depth
    st["params"]["adjust_subtypes"] = adjust_subtypes
    st["omnibus"], st["contrasts"] = _clean(omnibus), _clean(contrasts)
    return st["omnibus"]


def composition(adata, *, reference=None) -> pd.DataFrame:
    """Cell-type proportions vs age group (CLR, plus arcsine-root as a check).

    Proportions are relative: a rise in one population forces a fall in others.
    """
    st = store(adata)
    ref = reference or _need(adata, "params", "donor_summary").get("reference")
    if ref is None:
        raise ValueError("give reference=, or run nas.tl.age_effects first")
    comp, prop = _stats.composition_analysis(_obs(adata), ref)
    st["composition"], st["proportions"] = _clean(comp), prop
    return st["composition"]


def permutation_null(adata, *, n_perm=200, seed=0, stratify="batch") -> pd.DataFrame:
    """Realised false-positive rate: shuffle age labels between donors (within batch).

    Uses the same primary family and model as ``age_effects``.
    """
    summary = _need(adata, "donor_summary", "donor_summary")
    st = store(adata)
    _need(adata, "omnibus", "age_effects")
    p = st["params"]
    primary = _stats.build_test_family(st["applicability"], st["coverage"], st["score_columns"],
                                       tier=1.0)
    perm = _stats.permutation_null(summary, primary, p["reference"], n_perm=n_perm, seed=seed,
                                   stratify=stratify, adjust_depth=p.get("adjust_depth", True),
                                   adjust_subtypes=p.get("adjust_subtypes", True))
    st["permutation"] = perm
    st["permutation_rates"] = dict(perm.attrs)
    return perm


def module_correlation(adata, cell_type, *, method="spearman", score_cols=None) -> pd.DataFrame:
    """Donor-level correlation between modules within one population.

    If M04, M05 and M06 move together, report them as one interferon-driven
    signal rather than independent findings.
    """
    summary = _need(adata, "donor_summary", "donor_summary")
    return _stats.module_correlation(summary, cell_type, score_cols=score_cols, method=method)


def sensitivity(adata, *, min_cells_grid=(5, 10, 20, 30),
                statistics=("auto", "median", "mean")) -> pd.DataFrame:
    """Re-run the age tests over aggregation choices (min_cells x statistic)."""
    st = store(adata)
    ref = _need(adata, "params", "donor_summary").get("reference")
    if ref is None:
        raise RuntimeError("run nas.tl.age_effects first")
    out = _stats.sensitivity_sweep(_obs(adata), st["coverage"], st["applicability"],
                                   st["score_columns"], ref, min_cells_grid=min_cells_grid,
                                   statistics=statistics)
    st["sensitivity"] = _clean(out)
    return out


def pooled_scores(adata, cell_types, *, min_cells=10, label="pooled") -> pd.DataFrame:
    """Donor summaries that IGNORE cell identity across `cell_types`.

    Compare with the within-cell-type result to see how much of an apparent age
    effect is composition rather than cell state.
    """
    return _agg.pooled_scores(_obs(adata), store(adata)["score_columns"], cell_types,
                              min_cells=min_cells, label=label)


def null_control(adata, *, n_sets=50, seed=0, n_bins=25, n_jobs=1) -> pd.DataFrame:
    """Negative control: expression-matched random gene sets per tested module.

    For every tested (module, cell type) pair, `n_sets` random gene sets of the
    same size are drawn, each gene matched to a module gene on mean expression,
    then scored, aggregated and tested with the same model as ``age_effects``.

    Returned per pair:
      * ``null_frac_p05`` - share of random sets with p < 0.05. About 0.05 is
        expected; much higher means the population carries a generic age signal
        (depth, sample quality, composition), flagged as ``generic_age_signal``
        (> 0.15, a heuristic: about 3 SD above 0.05 for 50 sets).
      * ``empirical_p`` - how often a random set's F statistic is at least the
        module's. A module-specific age effect needs both a small ``p_adj_BH``
        and a small ``empirical_p``.
    """
    from . import controls as _controls
    om = _need(adata, "omnibus", "age_effects")
    st = store(adata)
    p = _params(adata)
    layer = None if st["layer"] == "X" else st["layer"]
    res, null_om = _controls.run_null_control(
        adata, _obs(adata), st["donor_summary"], om, st["modules_mapped"], layer=layer,
        n_sets=n_sets, seed=seed, n_bins=n_bins, min_cells=p.get("min_cells", 10),
        statistic=p.get("statistic", "auto"), reference=p["reference"], n_jobs=n_jobs,
        adjust_depth=p.get("adjust_depth", True), adjust_subtypes=p.get("adjust_subtypes", True))
    st["null_control"] = _clean(res)
    st["null_tests"] = _clean(null_om)
    return st["null_control"]
