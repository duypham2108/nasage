"""Preprocessing on your AnnData: design check and cell QC (scanpy-style, in place)."""
from __future__ import annotations

import warnings

import numpy as np

from . import qc as _qc
from ._design import canonical_obs, check_counts, donor_table, store


def check_design(adata, *, donor_key, age_key, celltype_key, batch_key=None, sex_key=None,
                 age_bins=None, age_order=None):
    """One row per donor (age group, batch, sex, cells) and confounding warnings.

    Raises if age, batch or sex vary within a donor. Warns when batch and age
    are completely confounded or an age group has fewer than 3 donors.
    """
    params = dict(donor_key=donor_key, age_key=age_key, celltype_key=celltype_key,
                  batch_key=batch_key, sex_key=sex_key, age_bins=age_bins, age_order=age_order)
    return donor_table(canonical_obs(adata, params))


def filter_cells(adata, *, celltype_key, layer=None, min_genes=200, min_counts=500,
                 max_counts=None, max_pct_mito=20.0, mito_adaptive=True, n_mads=3.0,
                 doublet_key=None, mito_prefix="MT-", inplace=True):
    """Remove low-quality cells, with thresholds per cell type.

    Computes ``scanpy.pp.calculate_qc_metrics`` on `layer` (raw counts; None =
    adata.X) unless its columns are already in ``adata.obs``. Every threshold is
    one number or ``{cell_type: value, "default": value}``. With
    ``mito_adaptive=True`` each cell type's mitochondrial cut-off is
    ``max(max_pct_mito, median + n_mads * MAD)`` within that type, so
    mitochondria-rich ciliated cells are not removed wholesale.

    Returns the per-cell-type report (also in
    ``adata.uns["nasage"]["qc"]``). With ``inplace=False`` returns
    ``(filtered_copy, report)`` and leaves `adata` unchanged.
    """
    import scanpy as sc
    check_counts(adata, layer)
    if celltype_key not in adata.obs:
        raise KeyError(f"column '{celltype_key}' not in adata.obs")
    needed = ("n_genes_by_counts", "total_counts", "pct_counts_mt")
    if not all(c in adata.obs for c in needed):
        if "mt" not in adata.var:
            adata.var["mt"] = np.asarray(adata.var_names.str.upper().str.startswith(mito_prefix))
        if not adata.var["mt"].any():
            warnings.warn("no mitochondrial genes found; the mitochondrial filter is skipped")
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], layer=layer, percent_top=None,
                                   log1p=False, inplace=True)
    keep, report = _qc.filter_mask(adata.obs, adata.obs[celltype_key], min_genes=min_genes,
                                   min_counts=min_counts, max_counts=max_counts,
                                   max_pct_mito=max_pct_mito, mito_adaptive=mito_adaptive,
                                   n_mads=n_mads, doublet_key=doublet_key)
    if not inplace:
        out = adata[keep.values].copy()
        store(out)["qc"] = report
        return out, report
    adata._inplace_subset_obs(keep.values)
    store(adata)["qc"] = report
    return report
