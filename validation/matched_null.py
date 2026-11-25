"""Expression- and size-matched random gene sets, the shared specificity control.

Every matched null in this project draws random gene sets from the same expression bins as the
tested module and scores them identically. Set SIZE was previously implicit: a draw took one
random gene per module gene, so size matched automatically within a run, but the null spread
was never reported as a function of size and could not be compared across modules.

`null_spread_by_size` makes that explicit. A 5-gene module and a 20-gene module have different
null distributions, so an empirical p from one is not comparable with the other, and a small
module can clear a null simply by being noisy enough to land in the tail less often.

Nothing here reads an outcome variable.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def expression_bins(expr: pd.Series, n_bins: int = 25) -> pd.Series:
    """Rank genes by mean expression and cut into equal-count bins."""
    return pd.qcut(expr.rank(method="first"), n_bins, labels=False)


def bin_pool(bins: pd.Series) -> dict:
    return {k: np.array(v) for k, v in bins.groupby(bins).groups.items()}


def draw_matched(genes, bins: pd.Series, pool: dict, rng) -> list:
    """One expression-matched random set of the SAME SIZE as `genes`."""
    gb = [bins.get(g) for g in genes]
    return [str(rng.choice(pool[k])) for k in gb if pd.notna(k)]


def null_spread_by_size(bins: pd.Series, pool: dict, score_fn, stat_fn,
                        sizes=(5, 8, 11, 14, 18, 22), n_draws=200, seed=0) -> pd.DataFrame:
    """Null distribution of the test statistic as a function of gene-set size.

    `score_fn(list_of_genes) -> per-sample score`; `stat_fn(score) -> float`. Returns the SD and
    selected quantiles of the null statistic at each size, so an empirical p can be read against
    the right reference rather than a pooled one.
    """
    rng = np.random.default_rng(seed)
    idx = np.array(bins.index)
    rows = []
    for n in sizes:
        vals = []
        for _ in range(n_draws):
            seed_genes = rng.choice(idx, size=min(n, len(idx)), replace=False)
            pick = draw_matched(list(seed_genes), bins, pool, rng)
            if len(pick) < 2:
                continue
            s = score_fn(pick)
            if s is None:
                continue
            v = stat_fn(s)
            if np.isfinite(v):
                vals.append(float(v))
        if vals:
            a = np.asarray(vals)
            rows.append(dict(size=n, n_draws=len(a), sd=float(a.std(ddof=1)),
                             q05=float(np.quantile(a, 0.05)),
                             q50=float(np.quantile(a, 0.50)),
                             q95=float(np.quantile(a, 0.95)),
                             abs_q95=float(np.quantile(np.abs(a), 0.95))))
    return pd.DataFrame(rows)


def empirical_p(observed: float, null_values) -> float:
    """Two-sided empirical p: fraction of null draws at least as extreme."""
    a = np.asarray([v for v in null_values if np.isfinite(v)], dtype=float)
    if a.size == 0:
        return float("nan")
    return float((np.abs(a) >= abs(observed)).mean())
