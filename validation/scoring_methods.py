"""Shared scoring methods for the specification curves.

**Correction, 2026-11-25.** The earlier `sc_rank` in multiverse.py ranked a programme's genes
*against each other* within a sample:

    order = np.argsort(np.argsort(-Xg, axis=1), axis=1) + 1    # Xg = samples x SET genes

Ranking inside the set is scale-free: it captures only the internal ordering of the set's own
genes and discards their expression level entirely, so it is not a depth-robust version of the
score -- it is a different and largely meaningless quantity. UCell ranks the set's genes against
the WHOLE transcriptome. That bug is why mean-z and "rank" correlated at 0.198 on average and
anticorrelated (-0.60 to -0.80) for M06.

`rank_score` below is the corrected version and takes the full expression matrix as context.
Every curve that used the old function must be re-run; see
results/multiverse/PREREG_RESULTS.md.
"""
import numpy as np


def _z(X):
    sd = X.std(0, ddof=1)
    return (X - X.mean(0)) / np.where(sd > 0, sd, 1)


def mean_z(X, full=None, idx=None):
    return _z(X).mean(1)


def median_z(X, full=None, idx=None):
    return np.median(_z(X), axis=1)


def pc1_of_set(X, full=None, idx=None):
    Z = _z(X)
    u, s, _ = np.linalg.svd(Z - Z.mean(0), full_matrices=False)
    v = u[:, 0] * s[0]
    r = np.corrcoef(v, Z.mean(1))[0, 1]
    return v * (np.sign(r) if np.isfinite(r) and r != 0 else 1.0)


def rank_score(X, full=None, idx=None, max_rank=1500):
    """UCell-style: rank the set's genes against ALL genes in the same sample.

    `full` is samples x all-genes; `idx` the column positions of the set's genes within it.
    Returns a score in [0, 1] that rises with the set's relative expression, computed from the
    Mann-Whitney U statistic with a rank cutoff, as UCell does.

    Without `full` this cannot be computed correctly and raises, rather than silently falling
    back to the broken within-set ranking.
    """
    if full is None or idx is None:
        raise ValueError("rank_score needs the full matrix; within-set ranking is invalid")
    n_samples, n_genes = full.shape
    m = len(idx)
    # descending rank of every gene within each sample: 1 = highest expressed
    order = np.argsort(-full, axis=1, kind="stable")
    ranks = np.empty_like(order)
    rows = np.arange(n_samples)[:, None]
    ranks[rows, order] = np.arange(1, n_genes + 1)[None, :]
    r = ranks[:, idx].astype(float)
    r = np.minimum(r, max_rank + 1)          # genes below the cutoff are tied at the cutoff
    u = r.sum(1) - m * (m + 1) / 2.0
    return np.maximum(0.0, 1.0 - u / (m * max_rank))


SCORERS = {"mean_z": mean_z, "median_z": median_z, "rank_ucell": rank_score,
           "pc1_of_set": pc1_of_set}
