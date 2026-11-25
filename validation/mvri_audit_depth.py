"""Is the shared axis among the five MVRI components a sequencing-depth artefact?

Dated addition, 2026-11-25, to docs/PREREG_mvri_presentation_audit.md. Not in the locked plan:
it was prompted by the audit's own output, where all five donor-level axis scores correlate
positively (Spearman 0.29-0.68, PC1 58%) even though the deck assigns the mucus axis the
opposite sign to the other four.

A mean-of-z score rises when more of its genes are detected, and detection rises with
sequencing depth, so donors sequenced more deeply score higher on *every* axis. This asks how
much of the shared structure that explains, and whether the age pattern survives adjusting for
it. Reported as a sensitivity analysis, not as a replacement for the locked tests.

Usage: python validation/mvri_audit_depth.py <h5ad> <audit_dir>
"""
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import stats

H5, D = Path(sys.argv[1]), Path(sys.argv[2])
AGE_ORDER = ["Neonate", "Infant", "Young child", "Child", "Adolescent", "Adult", "Eldly"]

S = pd.read_csv(D / "axis_means_by_agebin.csv")       # for the axis names only
AXES = [c for c in S.columns if c != "age"]

A = ad.read_h5ad(H5)
# Depth is computed FROM THE COUNT MATRIX, not from obs: nFeature_RNA / nCount_RNA are
# populated for only 53,150 of 148,978 cells here, and are absent for every donor that
# carries all five axis subsets -- so reading them would have silently produced no data.
import scipy.sparse as sp                                        # noqa: E402
Xc = A.X if sp.issparse(A.X) else sp.csr_matrix(A.X)
ngene_cell = np.asarray((Xc > 0).sum(1)).ravel().astype(float)
ncount_cell = np.asarray(Xc.sum(1)).ravel().astype(float)
meta = (pd.DataFrame({"donor": A.obs["donor"].astype(str).values,
                      "age": A.obs["Age_bin"].astype(str).values,
                      "ngene": ngene_cell, "ncount": ncount_cell})
        .groupby("donor")
        .agg(age=("age", "first"), mean_genes=("ngene", "mean"),
             mean_counts=("ncount", "mean"), n_cells=("ngene", "size")))
meta["log_depth"] = np.log10(meta.mean_counts)

# donor x axis score matrix, rebuilt from the per-axis donor tables the audit saved
P = pd.read_csv(D / "donor_axis_scores.csv", index_col=0)
# The composite exists only for donors that carry EVERY axis's cell subset, so P is a subset
# of the 52 donors. Align to P, not to meta, and report the shortfall.
print(f"composite computable for {len(P)} of {len(meta)} donors "
      f"({len(meta) - len(P)} lack at least one axis's cell subset)")
meta = meta.reindex(P.index)

rows = []
for a in AXES:
    if a not in P.columns:
        continue
    v = P[a]
    ok = v.notna() & meta.mean_genes.notna()
    r_dep = stats.spearmanr(v[ok], meta.mean_genes[ok])
    # age effect before and after removing depth, on donor means
    d = pd.DataFrame({"v": v[ok], "age": meta.age[ok], "g": meta.mean_genes[ok]})
    groups = [d.v[d.age == b].values for b in AGE_ORDER if (d.age == b).sum() > 0]
    s0, p0 = stats.kruskal(*[g for g in groups if len(g)])
    b = np.polyfit(d.g, d.v, 1)
    d["resid"] = d.v - np.polyval(b, d.g)
    groups_r = [d.resid[d.age == b_].values for b_ in AGE_ORDER if (d.age == b_).sum() > 0]
    s1, p1 = stats.kruskal(*[g for g in groups_r if len(g)])
    rows.append(dict(axis=a, rho_with_depth=float(r_dep.statistic),
                     p_depth=float(r_dep.pvalue), KW_p_raw=float(p0),
                     KW_p_depth_adjusted=float(p1)))
    print(f"{a:22s} rho(score, mean genes/cell) = {r_dep.statistic:+.3f} "
          f"(p {r_dep.pvalue:.1e}) | KW p {p0:.4f} -> {p1:.4f} after depth adjustment")

T = pd.DataFrame(rows)
T.to_csv(D / "depth_sensitivity.csv", index=False)

# how much of the between-axis correlation survives partialling out depth
Z = P[AXES].apply(lambda c: (c - c.mean()) / c.std(ddof=1))
resid = Z.apply(lambda c: c - np.polyval(np.polyfit(meta.mean_genes[c.notna()],
                                                    c[c.notna()], 1), meta.mean_genes))
C0 = Z.corr(method="spearman").values
C1 = resid.corr(method="spearman").values
iu = np.triu_indices(len(AXES), 1)
print(f"\nmean pairwise Spearman among axes: {C0[iu].mean():+.3f} raw "
      f"-> {C1[iu].mean():+.3f} after removing depth")
print(f"donor depth range: {meta.mean_genes.min():.0f}-{meta.mean_genes.max():.0f} "
      f"mean genes/cell across {len(meta)} donors")
r_age = stats.spearmanr(meta.mean_genes,
                        [AGE_ORDER.index(a) for a in meta.age])
print(f"depth vs age-bin rank: rho {r_age.statistic:+.3f}, p {r_age.pvalue:.3f}")
pd.DataFrame({"raw_mean_rho": [C0[iu].mean()], "depth_adj_mean_rho": [C1[iu].mean()],
              "rho_depth_vs_agebin": [float(r_age.statistic)],
              "p_depth_vs_agebin": [float(r_age.pvalue)]}).to_csv(
    D / "depth_summary.csv", index=False)
print("\nwrote", D / "depth_sensitivity.csv")
