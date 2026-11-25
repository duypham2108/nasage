"""Depth check for H3/H4, which the locked plan specified only for H2.

CDKN2A+ status is a detection event, so within a donor the CDKN2A+ cells are likely the
more deeply sequenced ones. Any module whose UCell score rises with detection depth would
then differ between the two groups for purely technical reasons. Tested by (a) comparing
depth between the groups and (b) repeating the contrast on depth-matched cells.
"""
import sys
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats

SUB, SCORES, OUT = sys.argv[1], sys.argv[2], sys.argv[3]

f = h5py.File(SUB, "r")
sym = np.array([x.decode() if isinstance(x, bytes) else x for x in f["var_symbol"][:]])
ct = np.array([x.decode() for x in f["obs/ct"][:]])
cil = np.flatnonzero(ct == "ciliated")
ip = f["X/indptr"][:].astype(np.int64)
j = int(np.flatnonzero(sym == "CDKN2A")[0])
cd = np.zeros(len(cil))
for k, r in enumerate(cil):
    lo, hi = ip[r], ip[r + 1]
    ind = f["X/indices"][lo:hi]
    p = np.searchsorted(ind, j)
    if p < hi - lo and ind[p] == j:
        cd[k] = f["X/data"][lo:hi][p]
f.close()

S = pd.read_csv(SCORES).iloc[cil].reset_index(drop=True)
S["pos"] = cd >= 1
S["depth"] = S.nis_n_genes

print("=== is CDKN2A+ a depth artefact at the cell level? ===")
print(f"CDKN2A+ cells: {int(S.pos.sum()):,} of {len(S):,} ({S.pos.mean():.1%})")
print(f"median detected genes: CDKN2A+ {S.depth[S.pos].median():.0f} vs "
      f"CDKN2A- {S.depth[~S.pos].median():.0f}")
w = stats.mannwhitneyu(S.depth[S.pos], S.depth[~S.pos], alternative="greater")
print(f"CDKN2A+ cells are deeper: p = {w.pvalue:.3g}")

# per-donor depth difference
dd = S.groupby("donor").apply(
    lambda g: g.depth[g.pos].median() - g.depth[~g.pos].median(), include_groups=False)
print(f"per-donor depth difference (CDKN2A+ minus -): median {dd.median():+.0f} genes, "
      f"positive in {int((dd > 0).sum())}/{len(dd)} donors")

# ---- depth-matched contrast: match within donor on decile of detected genes
print("\n=== depth-matched within-donor contrast ===")
rng = np.random.default_rng(0)
rows = []
for dn, g in S.groupby("donor"):
    g = g.copy()
    try:
        g["bin"] = pd.qcut(g.depth, 10, labels=False, duplicates="drop")
    except ValueError:
        continue
    a_idx, b_idx = [], []
    for b, gb in g.groupby("bin"):
        pa, pb = gb.index[gb.pos], gb.index[~gb.pos]
        n = min(len(pa), len(pb))
        if n == 0:
            continue
        a_idx += list(rng.choice(pa, n, replace=False))
        b_idx += list(rng.choice(pb, n, replace=False))
    if len(a_idx) < 20:
        continue
    A, B = g.loc[a_idx], g.loc[b_idx]
    rows.append(dict(donor=dn, n=len(a_idx),
                     depth_a=A.depth.median(), depth_b=B.depth.median(),
                     m02=A.M02_mucociliary_axoneme_UCell.median() - B.M02_mucociliary_axoneme_UCell.median(),
                     m11=A.M11_senescence_arrest_UCell.median() - B.M11_senescence_arrest_UCell.median()))
P = pd.DataFrame(rows)
P.to_csv(OUT, index=False)
print(f"donors with >= 20 matched pairs: {len(P)}; "
      f"residual median depth difference {np.median(P.depth_a - P.depth_b):+.0f} genes")
for col, lab, alt in (("m02", "H3 M02 axoneme lower in CDKN2A+", "less"),
                      ("m11", "H4 M11 higher in CDKN2A+", "greater")):
    w = stats.wilcoxon(P[col], alternative=alt)
    wo = stats.wilcoxon(P[col], alternative="greater" if alt == "less" else "less")
    print(f"{lab:34s} median delta={P[col].median():+.5f}  "
          f"p(locked dir)={w.pvalue:.3g}  p(opposite)={wo.pvalue:.3g}")
