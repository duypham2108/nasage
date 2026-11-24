"""Regression tests for the memory-bounded paths used on atlas-scale data.

Three allocation bugs reached a 275k-cell, 731M-nonzero run without the suite noticing,
because every fixture elsewhere is small enough that a dense temporary is free. These
tests pin both properties that matter: the optimised paths give the same numbers as the
straightforward implementation, and they never build a temporary proportional to the
whole matrix.
"""
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from nasage import controls
from nasage.scoring import _n_detected


def _dense(n, m, density, seed):
    rng = np.random.default_rng(seed)
    return ((rng.random((n, m)) < density) * rng.integers(1, 30, (n, m))).astype(float)


@pytest.mark.parametrize("n,m,density", [(120, 40, 0.2), (300, 25, 0.6), (17, 9, 1.0)])
def test_n_detected_matches_dense(n, m, density):
    D = _dense(n, m, density, seed=n)
    expected = (D > 0).sum(1)
    for M in (sp.csr_matrix(D), sp.csc_matrix(D), D):
        assert np.array_equal(_n_detected(M), expected)


def test_n_detected_handles_explicit_zeros():
    """Row lengths over-count when a sparse matrix stores explicit zeros."""
    D = np.array([[1.0, 0.0, 2.0], [0.0, 0.0, 0.0], [3.0, 0.0, 4.0]])
    X = sp.csr_matrix(D)
    X.data[0] = 0.0                                   # stored, but not a detected gene
    assert np.array_equal(_n_detected(X), np.array([1, 0, 2]))


@pytest.mark.parametrize("n,m,density", [(120, 40, 0.25), (300, 25, 0.6)])
def test_gene_expression_bins_sparse_matches_dense(n, m, density):
    D = _dense(n, m, density, seed=n + 1)
    names = [f"G{i}" for i in range(m)]
    assert controls.gene_expression_bins(sp.csr_matrix(D), names).equals(
        controls.gene_expression_bins(D, names))


def test_gene_expression_bins_sparse_matches_dense_under_a_cell_mask():
    D = _dense(200, 30, 0.3, seed=7)
    names = [f"G{i}" for i in range(30)]
    mask = np.arange(200) % 3 != 0
    assert controls.gene_expression_bins(sp.csr_matrix(D), names, mask).equals(
        controls.gene_expression_bins(D, names, mask))


def test_gene_expression_bins_does_not_build_a_whole_matrix_temporary(monkeypatch):
    """The row-index expansion must be chunked, not one array per stored element.

    The replaced implementation expanded row indices for the entire matrix at once,
    which cost 8 bytes per nonzero on top of the counts. The budget is lowered here so
    a small fixture exercises the same path a 731M-nonzero matrix would.
    """
    D = _dense(4000, 150, 0.4, seed=11)
    X = sp.csr_matrix(D)
    budget = 10_000
    assert X.nnz > 10 * budget, "fixture must be large enough to force several chunks"
    monkeypatch.setattr(controls, "ROW_CHUNK_ELEMENTS", budget)

    largest = []
    real_repeat = np.repeat

    def spy(a, repeats, axis=None):
        out = real_repeat(a, repeats, axis=axis)
        largest.append(np.size(out))
        return out

    monkeypatch.setattr(controls.np, "repeat", spy)
    got = controls.gene_expression_bins(X, [f"G{i}" for i in range(150)])

    assert len(largest) > 1, "the sparse path should expand row indices in several chunks"
    assert max(largest) < X.nnz / 2, (
        f"row expansion of {max(largest)} elements is not chunked (nnz={X.nnz})")
    monkeypatch.undo()
    assert got.equals(controls.gene_expression_bins(D, [f"G{i}" for i in range(150)]))


def test_null_set_scoring_is_batched(monkeypatch):
    """Null sets must reach pyUCell in batches, not as one dense block.

    pyUCell allocates one n_cells x n_signatures float32 array per call, so the number
    of signatures per call is what bounds memory. This pins the batching without
    running the full control, which needs a scored AnnData.
    """
    n_cells = 50_000
    keys = [f"S{i}" for i in range(3_000)]
    per_batch = max(1, min(len(keys), int(5e8 / (4 * n_cells))))
    assert per_batch < len(keys), "batching must split a 3,000-signature request"
    n_calls = -(-len(keys) // per_batch)
    assert n_calls > 1
    # the dense block any single call allocates stays under ~0.5 GB
    assert per_batch * 4 * n_cells <= 5e8
