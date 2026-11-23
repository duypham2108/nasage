"""Cell filtering on scanpy QC metrics, with per-cell-type thresholds.

Works on the columns written by ``scanpy.pp.calculate_qc_metrics``:
``n_genes_by_counts``, ``total_counts`` and ``pct_counts_mt``.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def _per_type(value, cell_types: pd.Series) -> pd.Series:
    """Threshold per cell: a number, or {cell_type: number, "default": number}.

    A cell type missing from the dict (with no "default") is not filtered on
    that metric.
    """
    if isinstance(value, dict):
        default = value.get("default")
        return cell_types.astype(str).map(lambda c: value.get(c, default)).astype(float)
    return pd.Series(np.nan if value is None else float(value), index=cell_types.index)


def filter_mask(obs: pd.DataFrame, cell_types: pd.Series, min_genes=200, min_counts=500,
                max_counts=None, max_pct_mito=20.0, mito_adaptive=True, n_mads=3.0,
                doublet_key=None) -> tuple[pd.Series, pd.DataFrame]:
    """Boolean keep-mask and a per-cell-type report of the thresholds used.

    Every threshold is one number or a dict per cell type, e.g.
    ``min_genes={"default": 200, "neutrophil": 100}``: neutrophils carry few
    transcripts and a single global floor removes most of them.

    Mitochondrial fraction differs between cell types (ciliated cells are
    mitochondria-rich), so with ``mito_adaptive=True`` each type's cut-off is
    ``max(max_pct_mito, median + n_mads * MAD)`` within that type: the fixed
    value is a floor, raised for types whose healthy cells sit higher.
    """
    ct = cell_types.astype(str).reindex(obs.index)
    genes, counts, mito = obs["n_genes_by_counts"], obs["total_counts"], obs["pct_counts_mt"]
    t_genes, t_counts = _per_type(min_genes, ct), _per_type(min_counts, ct)
    t_maxc, t_mito = _per_type(max_counts, ct), _per_type(max_pct_mito, ct)
    if mito_adaptive and max_pct_mito is not None:
        for _, idx in ct.groupby(ct).groups.items():
            v = mito.loc[idx]
            mad = 1.4826 * float(np.median(np.abs(v - v.median())))
            t_mito.loc[idx] = np.fmax(t_mito.loc[idx], float(v.median()) + n_mads * mad)

    keep = pd.Series(True, index=obs.index)
    keep &= ~(genes < t_genes)            # NaN threshold -> comparison False -> kept
    keep &= ~(counts < t_counts)
    keep &= ~(counts > t_maxc)
    if max_pct_mito is not None and mito.max() > 0:
        keep &= ~(mito > t_mito)
    if doublet_key:
        keep &= ~obs[doublet_key].astype(bool)

    rows = []
    for c, idx in ct.groupby(ct).groups.items():
        rows.append(dict(cell_type=c, n_cells_in=len(idx), n_cells_kept=int(keep.loc[idx].sum()),
                         min_genes=t_genes.loc[idx].iloc[0], min_counts=t_counts.loc[idx].iloc[0],
                         max_pct_mito_used=t_mito.loc[idx].iloc[0],
                         median_genes=float(genes.loc[idx].median()),
                         median_pct_mito=float(mito.loc[idx].median())))
    report = pd.DataFrame(rows)
    report["frac_kept"] = report.n_cells_kept / report.n_cells_in
    return keep, report
