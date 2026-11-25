"""Corrections and sensitivity analyses for the Figure 1/3 replication.

1. H4 as locked in docs/PREREG_crownlands_figure_replication.md is CIRCULAR: cells are
   selected on CDKN2A and then scored with M11, which contains CDKN2A. The test is
   recomputed on M11 minus CDKN2A (CDKN1A, SERPINE1, IGFBP7, GDF15). This is a
   correction of an invalid test, not a rescue of a failed one.
2. H5 depth sensitivity (post-hoc, labelled as such): CDKN2A+ fraction correlates with
   median detected genes at rho = 0.73, so H5 is re-fitted with depth in the model.
3. H2 depth-matched check (plan section 4.1 item 2, not implemented in the first run):
   CDKN2A counts are downsampled per cell to a common depth before recomputing.

Usage: python validation/crownlands_figrep_corrections.py <subsample.h5ad> <scores.csv.gz> <figrep_dir>
"""
import sys
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats
import statsmodels.formula.api as smf

SUB, SCORES, D = sys.argv[1], sys.argv[2], Path(sys.argv[3])
M11_NO_CDKN2A = ["CDKN1A", "SERPINE1", "IGFBP7", "GDF15"]

# ---------------------------------------------------------------- 1. corrected H4
f = h5py.File(SUB, "r")
sym = np.array([x.decode() if isinstance(x, bytes) else x for x in f["var_symbol"][:]])
ct = np.array([x.decode() for x in f["obs/ct"][:]])
shape = tuple(int(x) for x in f["X"].attrs["shape"])
cil = np.flatnonzero(ct == "ciliated")
ip = f["X/indptr"][:].astype(np.int64)

# pull only the ciliated rows
rows_d, rows_i, lens = [], [], []
for r in cil:
    lo, hi = ip[r], ip[r + 1]
    rows_d.append(f["X/data"][lo:hi]); rows_i.append(f["X/indices"][lo:hi]); lens.append(hi - lo)
X = sp.csr_matrix((np.concatenate(rows_d), np.concatenate(rows_i),
                   np.concatenate([[0], np.cumsum(lens)])), shape=(len(cil), shape[1]))
f.close()
print(f"ciliated subsample: {X.shape}")

A = ad.AnnData(X=X, obs=pd.DataFrame(index=[str(i) for i in cil]),
               var=pd.DataFrame(index=pd.Index(sym)))
from pyucell import compute_ucell_scores
from nasage.scoring import MAX_RANK
compute_ucell_scores(A, {"M11_noCDKN2A": M11_NO_CDKN2A}, max_rank=MAX_RANK,
                     suffix="_UCell", n_jobs=1)

S = pd.read_csv(SCORES).iloc[cil].reset_index(drop=True)
S["M11_noCDKN2A_UCell"] = A.obs["M11_noCDKN2A_UCell"].to_numpy()
cd = np.asarray(X[:, int(np.flatnonzero(sym == "CDKN2A")[0])].todense()).ravel()
S["pos"] = cd >= 1

rows = []
for dn, g in S.groupby("donor"):
    a, b = g[g.pos], g[~g.pos]
    if len(a) < 20 or len(b) < 20:
        continue
    rows.append(dict(donor=dn,
                     m11_full=a.M11_senescence_arrest_UCell.median() - b.M11_senescence_arrest_UCell.median(),
                     m11_nocdkn2a=a.M11_noCDKN2A_UCell.median() - b.M11_noCDKN2A_UCell.median(),
                     m02=a.M02_mucociliary_axoneme_UCell.median() - b.M02_mucociliary_axoneme_UCell.median()))
P = pd.DataFrame(rows)
P.to_csv(D / "within_donor_corrected.csv", index=False)

print(f"\n=== H4 corrected, n = {len(P)} donors ===")
for col, lab in (("m11_full", "M11 as locked (contains CDKN2A)"),
                 ("m11_nocdkn2a", "M11 minus CDKN2A (4 genes)")):
    w = stats.wilcoxon(P[col], alternative="greater")
    print(f"{lab:38s} median delta={P[col].median():+.5f}  "
          f"risen in {int((P[col] > 0).sum())}/{len(P)}  p={w.pvalue:.3g}")
w = stats.wilcoxon(P.m02, alternative="less")
print(f"{'H3 M02 axoneme (locked direction)':38s} median delta={P.m02.median():+.5f}  "
      f"lower in {int((P.m02 < 0).sum())}/{len(P)}  p_one_sided={w.pvalue:.3g}")
print(f"{'H3 opposite direction':38s} p={stats.wilcoxon(P.m02, alternative='greater').pvalue:.3g}")

# ---------------------------------------------------------------- 2. H5 with depth
comp = pd.read_csv(D / "donor_composition.csv")
p16 = pd.read_csv(D / "donor_cdkn2a.csv")
J = comp.dropna(subset=["age"]).merge(
    p16.dropna(subset=["age"])[["donor", "frac_pos", "med_genes"]], on="donor")
print(f"\n=== H5 sensitivity (post-hoc), n = {len(J)} ===")
for formula, lab in (
        ("frac_basal ~ frac_pos + age + sex", "as locked (age + sex)"),
        ("frac_basal ~ frac_pos + age + sex + np.log(med_genes)", "+ log median genes")):
    m = smf.ols(formula, data=J).fit()
    t, p2 = m.tvalues["frac_pos"], m.pvalues["frac_pos"]
    print(f"{lab:28s} beta={m.params['frac_pos']:+.5f}  p_two={p2:.4g}  "
          f"p_one={(p2 / 2 if t < 0 else 1 - p2 / 2):.4g}")
r = stats.spearmanr(J.frac_basal, J.med_genes)
print(f"basal fraction vs median genes: rho={r.statistic:+.3f} p={r.pvalue:.4g}")

# ---------------------------------------------------------------- 3. H2 depth-matched
cilc = pd.read_csv(D / "ciliated_cdkn2a.csv.gz")
target = float(np.percentile(cilc.total, 10))
rng = np.random.default_rng(0)
keep = cilc.total >= target
p = np.clip(target / cilc.total[keep].to_numpy(), 0, 1)
ds = rng.binomial(cilc.cdkn2a[keep].to_numpy().astype(int), p)
sub = cilc[keep].assign(pos_ds=ds >= 1)
g = sub.groupby("donor")
dd = pd.DataFrame({"n": g.size(), "pos": g.pos_ds.sum(), "age": g.age.first()}).reset_index()
dd = dd[(dd.n >= 200)].dropna(subset=["age"])
dd["frac"] = dd.pos / dd.n
r = stats.spearmanr(dd.age, dd.frac)
print(f"\n=== H2 depth-matched (downsampled to {target:,.0f} counts/cell), "
      f"n = {len(dd)} donors ===")
print(f"CDKN2A+ fraction vs age: rho={r.statistic:+.3f}  "
      f"p_one={(r.pvalue / 2 if r.statistic > 0 else 1 - r.pvalue / 2):.4g}")
print(f"cells retained: {int(keep.sum()):,} of {len(cilc):,}")
