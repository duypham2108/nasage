"""The only valid form of H4: M11 WITHOUT CDKN2A, on depth-matched cells.

The locked H4 is circular (M11 contains CDKN2A, cells are selected on CDKN2A) and the
groups also differ in sequencing depth by ~600 detected genes. Both must be removed
together for the test to mean anything.
"""
import sys, h5py, numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad
from scipy import stats
SUB, SCORES, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
M11_NO = ["CDKN1A", "SERPINE1", "IGFBP7", "GDF15"]

f = h5py.File(SUB, "r")
sym = np.array([x.decode() if isinstance(x, bytes) else x for x in f["var_symbol"][:]])
ct = np.array([x.decode() for x in f["obs/ct"][:]])
cil = np.flatnonzero(ct == "ciliated")
ip = f["X/indptr"][:].astype(np.int64)
d_, i_, l_ = [], [], []
for r in cil:
    lo, hi = ip[r], ip[r + 1]
    d_.append(f["X/data"][lo:hi]); i_.append(f["X/indices"][lo:hi]); l_.append(hi - lo)
X = sp.csr_matrix((np.concatenate(d_), np.concatenate(i_),
                   np.concatenate([[0], np.cumsum(l_)])), shape=(len(cil), len(sym)))
f.close()

A = ad.AnnData(X=X, obs=pd.DataFrame(index=[str(i) for i in cil]),
               var=pd.DataFrame(index=pd.Index(sym)))
from pyucell import compute_ucell_scores
from nasage.scoring import MAX_RANK
compute_ucell_scores(A, {"M11_noCDKN2A": M11_NO}, max_rank=MAX_RANK, suffix="_UCell", n_jobs=1)

S = pd.read_csv(SCORES).iloc[cil].reset_index(drop=True)
S["m11nc"] = A.obs["M11_noCDKN2A_UCell"].to_numpy()
S["pos"] = np.asarray(X[:, int(np.flatnonzero(sym == "CDKN2A")[0])].todense()).ravel() >= 1
S["depth"] = S.nis_n_genes

rng = np.random.default_rng(0)
rows = []
for dn, g in S.groupby("donor"):
    g = g.copy()
    try:
        g["bin"] = pd.qcut(g.depth, 10, labels=False, duplicates="drop")
    except ValueError:
        continue
    a, b = [], []
    for _, gb in g.groupby("bin"):
        pa, pb = gb.index[gb.pos], gb.index[~gb.pos]
        n = min(len(pa), len(pb))
        if n:
            a += list(rng.choice(pa, n, replace=False)); b += list(rng.choice(pb, n, replace=False))
    if len(a) < 20:
        continue
    A_, B_ = g.loc[a], g.loc[b]
    rows.append(dict(donor=dn, n=len(a), dd=A_.depth.median() - B_.depth.median(),
                     m11nc=A_.m11nc.median() - B_.m11nc.median(),
                     m11full=A_.M11_senescence_arrest_UCell.median() - B_.M11_senescence_arrest_UCell.median()))
P = pd.DataFrame(rows); P.to_csv(OUT, index=False)
print(f"donors: {len(P)}; residual depth difference median {P.dd.median():+.0f} genes")
for c, lab in (("m11full", "M11 as locked (contains CDKN2A)"),
               ("m11nc", "M11 minus CDKN2A -- the valid test")):
    w = stats.wilcoxon(P[c], alternative="greater")
    print(f"{lab:34s} median delta={P[c].median():+.5f}  risen in {int((P[c]>0).sum())}/{len(P)}"
          f"  p={w.pvalue:.3g}")
