"""Negative control: expression-matched random gene sets.

For each multi-gene module, draw random gene sets of the same size, each gene
matched to a module gene on mean expression (expression-bin matching, as used
for control gene sets in module scoring, e.g. Tirosh et al. 2016). Score them
with UCell, aggregate per donor exactly like the module, and test them with the
same model on the same (module, cell type) pairs.

Random genes share no biology, so
  * if many random sets also "differ by age" in a cell type, the age signal there
    is generic (sequencing depth, sample quality, composition), not the module;
  * a module's age effect is specific only if it is stronger than most of its
    matched random sets (empirical p-value).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp

from . import aggregate as _agg
from . import stats as _stats
from .scoring import MAX_RANK

NULL_TAG = "__null"


def gene_expression_bins(X, var_names, cells_mask=None, n_bins=25) -> pd.Series:
    """Mean log1p(CP10K) per gene over the analysed cells, cut into quantile bins."""
    if cells_mask is not None:
        X = X[np.asarray(cells_mask)]
    lib = np.maximum(np.asarray(X.sum(1)).ravel(), 1)
    if sp.issparse(X):
        # Accumulate the per-gene sum of log1p(CP10K) straight off the stored values.
        # Scaling the matrix first (diags @ X) would allocate a second copy of every
        # stored element, which is several GB on a large cohort.
        Xc = X.tocsr() if not sp.isspmatrix_csr(X) else X
        n_rows, n_cols = Xc.shape
        acc = np.zeros(n_cols)
        step = max(1, int(2e7 // max(1, Xc.nnz // max(n_rows, 1) + 1)))
        for s in range(0, n_rows, step):
            e = min(s + step, n_rows)
            lo, hi = Xc.indptr[s], Xc.indptr[e]
            if hi == lo:
                continue
            rows = np.repeat(np.arange(s, e), np.diff(Xc.indptr[s:e + 1]))
            acc += np.bincount(Xc.indices[lo:hi],
                               weights=np.log1p(Xc.data[lo:hi] * (1e4 / lib)[rows]),
                               minlength=n_cols)
        mean = acc / n_rows
    else:
        mean = np.log1p(np.asarray(X) * (1e4 / lib)[:, None]).mean(0)
    s = pd.Series(mean, index=pd.Index(map(str, var_names)))
    expressed = s[s > 0]
    bins = pd.qcut(expressed.rank(method="first"), n_bins, labels=False)
    return bins.reindex(s.index)            # NaN for never-detected genes


def matched_null_sets(module_genes: dict, bins: pd.Series, n_sets=50, seed=0,
                      exclude=()) -> dict:
    """{f"{module}__null{i}": [genes]} with each gene drawn from its partner's bin."""
    rng = np.random.default_rng(seed)
    excluded = set(exclude) | {g for gl in module_genes.values() for g in gl}
    pool = {b: np.array([g for g in idx if g not in excluded])
            for b, idx in bins.dropna().groupby(bins.dropna()).groups.items()}
    out = {}
    for mod, genes in module_genes.items():
        gb = [bins.get(g) for g in genes if pd.notna(bins.get(g))]
        if len(gb) < 3:
            continue
        for i in range(n_sets):
            chosen = []
            for b in gb:
                cand = [g for g in pool.get(b, []) if g not in chosen]
                if cand:
                    chosen.append(str(rng.choice(cand)))
            out[f"{mod}{NULL_TAG}{i}"] = chosen
    return out


def run_null_control(adata, obs: pd.DataFrame, summary: pd.DataFrame, omnibus: pd.DataFrame,
                     modules_mapped: dict, layer=None, n_sets=50, seed=0, n_bins=25,
                     min_cells=10, statistic="auto", reference=None, n_jobs=1,
                     adjust_depth=True, adjust_subtypes=True,
                     tiers=("primary", "exploratory")) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Score, aggregate and test matched random gene sets for every tested UCell pair.

    `tiers` limits the control to pairs of those tiers. Only cells of the tested
    cell types are scored; UCell ranks genes within each cell, so their scores
    do not depend on which other cells are present. Expression bins for gene
    matching still use all analysed cells.

    Returns (per-pair comparison table, all null-set test results).
    """
    import anndata as ad
    from pyucell import compute_ucell_scores

    tested = omnibus[omnibus.estimable.astype(bool)
                     & omnibus.score_column.str.endswith("_UCell")
                     & omnibus.tier.isin(list(tiers))]
    if tested.empty:
        raise ValueError(f"no estimable UCell pairs in tiers {tuple(tiers)}")
    mods = sorted(set(tested.module_id))
    X = adata.X if layer is None else adata.layers[layer]
    in_obs = adata.obs_names.isin(obs.index)
    bins = gene_expression_bins(X, adata.var_names, in_obs, n_bins=n_bins)
    sets = matched_null_sets({m: modules_mapped[m] for m in mods}, bins, n_sets=n_sets,
                             seed=seed)

    # score on a bare AnnData so the user's adata.obs is not filled with null columns
    cts = set(tested.cell_type.astype(str))
    score_ids = obs.index[obs.cell_type.astype(str).isin(cts)]
    sel = adata.obs_names.isin(score_ids)
    tmp = ad.AnnData(X=X[sel], obs=pd.DataFrame(index=adata.obs_names[sel]),
                     var=pd.DataFrame(index=adata.var_names))
    # pyUCell allocates one dense n_cells x n_signatures block, so hundreds of null
    # sets on a large cohort exhaust memory. Score in batches sized to ~0.5 GB and
    # collect the columns; the scores themselves are unchanged by the batching.
    keys = list(sets)
    per_batch = max(1, min(len(keys), int(5e8 / (4 * max(tmp.n_obs, 1)))))
    for b in range(0, len(keys), per_batch):
        compute_ucell_scores(tmp, {k: sets[k] for k in keys[b:b + per_batch]},
                             max_rank=MAX_RANK, suffix="_UCell", n_jobs=n_jobs)
    null_cols = [f"{k}_UCell" for k in sets]
    o = obs.loc[score_ids].join(tmp.obs[null_cols])
    null_summary = _agg.donor_celltype_summary(o, null_cols, min_cells=min_cells,
                                               statistic=statistic)

    pairs = [(f"{m}{NULL_TAG}{i}", ct, "null")
             for m, ct in tested[["module_id", "cell_type"]].itertuples(index=False)
             for i in range(n_sets) if f"{m}{NULL_TAG}{i}" in sets]
    null_om, _ = _stats.fit_age_models(null_summary, pairs, reference,
                                       adjust_depth=adjust_depth, adjust_subtypes=adjust_subtypes)
    null_om["module_id_real"] = null_om.module_id.str.split(NULL_TAG).str[0]

    rows = []
    for r in tested.itertuples(index=False):
        nul = null_om[(null_om.module_id_real == r.module_id) & (null_om.cell_type == r.cell_type)
                      & null_om.estimable.astype(bool)]
        p_null = nul.p_omnibus.astype(float).to_numpy()
        p_obs = float(r.p_omnibus)
        f_obs = float(r.f_stat)
        f_null = nul.f_stat.astype(float).to_numpy()
        rows.append(dict(
            module_id=r.module_id, cell_type=r.cell_type, tier=r.tier, p_omnibus=p_obs,
            n_null=len(p_null),
            null_frac_p05=float((p_null < 0.05).mean()) if len(p_null) else np.nan,
            null_median_p=float(np.median(p_null)) if len(p_null) else np.nan,
            empirical_p=float((1 + (f_null >= f_obs).sum()) / (1 + len(f_null)))
            if len(f_null) else np.nan))
    res = pd.DataFrame(rows)
    # a cell type whose random sets are often "significant" carries a generic age signal
    res["generic_age_signal"] = res.null_frac_p05 > 0.15
    return res, null_om
