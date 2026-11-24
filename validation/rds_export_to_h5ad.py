"""Rebuild an AnnData (raw counts in .X) from the files written by convert_seurat_rds.R.

Usage:
    python validation/rds_export_to_h5ad.py <export_dir> <out.h5ad>
        [--subset COL=v1,v2] [--recover-lognorm]

--subset keeps only cells whose obs COL is one of the listed values.
--recover-lognorm: the exported matrix is log1p(count / total * scale) (Seurat
  LogNormalize) rather than raw counts. Within a cell, the smallest nonzero value
  corresponds to a count of 1, so count = expm1(x) / expm1(min_x). Recovered values
  are checked to be integers; a cell whose true minimum count is >1 would come out
  scaled down, so cells failing the check are reported and dropped.
"""
import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

ap = argparse.ArgumentParser()
ap.add_argument("export_dir"); ap.add_argument("out")
ap.add_argument("--subset", action="append", default=[])
ap.add_argument("--recover-lognorm", action="store_true")
args = ap.parse_args()
d, out = Path(args.export_dir), Path(args.out)

n_genes, n_cells = (int(x) for x in (d / "counts_dim.txt").read_text().split())
genes = (d / "genes.txt").read_text(encoding="utf8").splitlines()
cells = (d / "cells.txt").read_text(encoding="utf8").splitlines()
obs = pd.read_csv(d / "obs.csv", index_col=0, low_memory=False)
assert list(obs.index) == cells, "obs rows do not match cells"

keep = np.ones(n_cells, dtype=bool)
for s in args.subset:
    col, vals = s.split("=", 1)
    keep &= obs[col].astype(str).isin(vals.split(",")).values
rows = np.flatnonzero(keep)
print(f"{keep.sum()} of {n_cells} cells kept")

p = np.fromfile(d / "counts_p.bin", dtype=np.int32).astype(np.int64)
i = np.fromfile(d / "counts_i.bin", dtype=np.int32)
x = np.fromfile(d / "counts_x.bin", dtype=np.float32)
# genes x cells CSC  ==  cells x genes CSR; select rows before any heavy work
X = sp.csr_matrix((x, i, p), shape=(n_cells, n_genes))[rows]
del i, x
X.sort_indices()

if args.recover_lognorm:
    data = np.expm1(X.data.astype(np.float64))
    nnz = np.diff(X.indptr)
    ok = nnz > 0
    rowmin = np.zeros(X.shape[0])
    rowmin[ok] = np.minimum.reduceat(data, X.indptr[:-1][ok])
    counts = data / np.repeat(rowmin, nnz)
    dev = np.abs(counts - np.round(counts))
    rowdev = np.zeros(X.shape[0])
    rowdev[ok] = np.maximum.reduceat(dev, X.indptr[:-1][ok])
    bad = (~ok) | (rowdev > 0.01)
    print(f"count recovery: max deviation from integer {rowdev[ok].max():.2e}; "
          f"{bad.sum()} cells fail (dropped)")
    X.data = np.round(counts).astype(np.float32)
    X = X[~bad]
    rows = rows[~bad]

o = obs.iloc[rows].copy()
o["n_counts_recovered"] = np.asarray(X.sum(1)).ravel()
o["n_genes_recovered"] = np.diff(X.indptr)
for c in o.columns:
    if o[c].dtype == object:
        o[c] = o[c].astype("category")
a = ad.AnnData(X=X, obs=o, var=pd.DataFrame(index=pd.Index(genes, name="gene")))
a.var_names_make_unique()
a.write_h5ad(out, compression="gzip")
print(a)
print(o[["n_counts_recovered", "n_genes_recovered"]].describe().round(0))
