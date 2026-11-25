"""External test of the manuscript's 28 cilium-assembly genes in Crownlands ciliated cells.

The gene list was derived in the primary cohort (correlation with a senescence score, then
filtered to genes highest in the elderly), so testing it here is an external test of that
list, not a re-derivation. Crownlands has no donor under 18, so only the elderly arm of the
claim -- that cilium assembly rises again in the elderly -- can be tested.

Primary score: mean across the 28 genes of z-scored log1p(CP10K), z across ciliated cells,
averaged per donor. This is what a mean-expression dot plot shows. UCell is a sensitivity.

Controls, all pre-specified here before the result is seen:
  * batch (flex_version x device_guided) and sex in the model;
  * log median detected genes per cell in the model;
  * a depth-matched variant using only cells in a fixed global band of detected genes;
  * 500 expression-matched random 28-gene sets, 25 bins, for an empirical p.

Usage: python validation/crownlands_cilium28.py <subsample.h5ad> <atlas.h5ad> <out_dir>
"""
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats
import statsmodels.formula.api as smf

SUB, ATLAS, OUT = sys.argv[1], sys.argv[2], Path(sys.argv[3])
OUT.mkdir(parents=True, exist_ok=True)
N_SETS, N_BINS = 500, 25
BAND = (2000, 4000)          # detected-gene window for the depth-matched variant
MIN_CELLS = 20

CILIUM28 = ["SEPTIN2", "DNAH3", "SEPTIN7", "ARMC2", "CFAP61", "TTC29", "SEPTIN9", "CFAP206",
            "B9D1", "NEK10", "CC2D2A", "LRGUK", "DNAH6", "TTC12", "KIF27", "BBS9", "CFAP298",
            "DMD", "IFT57", "ARMC9", "DRC1", "ARL3", "CCDC39", "CFAP47", "BBS5", "PIBF1",
            "ROPN1L", "PARVA"]

# ---------------------------------------------------------------- donor-level batch
a = h5py.File(ATLAS, "r")
def acol(n):
    g = a["obs"][n]
    if isinstance(g, h5py.Group):
        c = np.array([x.decode() if isinstance(x, bytes) else x for x in g["categories"][:]])
        return c[g["codes"][:]]
    v = g[:]
    return np.array([x.decode() if isinstance(x, bytes) else x for x in v])
bt = pd.DataFrame({"donor": acol("donor_id"), "flex": acol("flex_version"),
                   "dev": acol("device_guided")}).drop_duplicates("donor")
a.close()
bt["batch"] = bt.flex + "_" + bt.dev

# ---------------------------------------------------------------- ciliated cells
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
obs = pd.DataFrame({c: np.array([x.decode() for x in f[f"obs/{c}"][:]])[cil]
                    for c in ("donor", "age", "age_bin", "sex")})
f.close()
obs["age"] = pd.to_numeric(obs.age, errors="coerce")
obs["n_genes"] = np.diff(X.indptr)
obs = obs.merge(bt[["donor", "batch"]], on="donor", how="left")
print(f"ciliated cells {X.shape[0]:,}; donors {obs.donor.nunique()}")
print("batch:", obs.drop_duplicates('donor').batch.value_counts().to_dict())

present = [g for g in CILIUM28 if g in set(sym)]
print(f"\ncilium genes present: {len(present)}/{len(CILIUM28)}"
      + (f"; missing {sorted(set(CILIUM28) - set(present))}" if len(present) < len(CILIUM28) else ""))

# ---------------------------------------------------------------- scoring
lib = np.maximum(np.asarray(X.sum(1)).ravel(), 1)
Xc = X.tocsr()
rows = np.repeat(np.arange(Xc.shape[0]), np.diff(Xc.indptr))
logv = np.log1p(Xc.data * (1e4 / lib)[rows])
L = sp.csr_matrix((logv, Xc.indices, Xc.indptr), shape=Xc.shape)
mu = np.asarray(L.mean(0)).ravel()
sq = np.asarray(L.multiply(L).mean(0)).ravel()
sd = np.sqrt(np.maximum(sq - mu ** 2, 0))
ok = sd > 0
# expression bins for the matched null (same construction as the package's control)
s = pd.Series(mu, index=sym)
expressed = s[s > 0]
bins = pd.qcut(expressed.rank(method="first"), N_BINS, labels=False).reindex(s.index)

gidx = {g: int(np.flatnonzero(sym == g)[0]) for g in present}


def score(idx):
    """mean z of log1p(CP10K) over the given gene columns, per cell."""
    idx = [i for i in idx if ok[i]]
    if not idx:
        return np.full(X.shape[0], np.nan)
    sub = np.asarray(L[:, idx].todense())
    return ((sub - mu[idx]) / sd[idx]).mean(1)


obs["sc"] = score(list(gidx.values()))

def donor_table(mask=None):
    o = obs if mask is None else obs[mask]
    g = o.groupby("donor")
    t = pd.DataFrame({"sc": g.sc.median(), "n": g.size(), "age": g.age.first(),
                      "age_bin": g.age_bin.first(), "sex": g.sex.first(),
                      "batch": g.batch.first(),
                      "log_genes": np.log(g.n_genes.median())}).reset_index()
    return t[(t.n >= MIN_CELLS) & t.age.notna()]


def test(t, col="sc"):
    """one-sided p for the manuscript's direction (higher in 65+), plus the age slope."""
    t = t[t.age_bin.isin(["40-64", "65+"])].copy()
    t["old"] = (t.age_bin == "65+").astype(int)
    m = smf.ols(f"{col} ~ old + sex + log_genes + C(batch)", data=t).fit()
    tv, p2 = m.tvalues["old"], m.pvalues["old"]
    return dict(n=len(t), beta=m.params["old"], p_two=p2,
                p_one=(p2 / 2 if tv > 0 else 1 - p2 / 2))


T = donor_table()
print("\n=== donors per age bin ===")
print(T.age_bin.value_counts().reindex(["18-39", "40-64", "65+"]).to_dict())
print("\nmean z by age bin:")
print(T.groupby("age_bin").sc.agg(["median", "mean", "size"]).round(4).to_string())

r_all = test(T)
print(f"\n=== 40-64 vs 65+, adjusted for sex, depth and batch ===")
print(f"beta = {r_all['beta']:+.4f}  p(one-sided, higher in 65+) = {r_all['p_one']:.4g}  "
      f"p(two-sided) = {r_all['p_two']:.4g}  n = {r_all['n']}")
rho = stats.spearmanr(T.age, T.sc)
print(f"continuous age: rho = {rho.statistic:+.3f}, p = {rho.pvalue:.4g}")

band = (obs.n_genes >= BAND[0]) & (obs.n_genes <= BAND[1])
Tb = donor_table(band)
r_band = test(Tb)
print(f"\n=== depth-matched (cells with {BAND[0]}-{BAND[1]} detected genes) ===")
print(f"donors {len(Tb)}; beta = {r_band['beta']:+.4f}  "
      f"p(one-sided) = {r_band['p_one']:.4g}  p(two-sided) = {r_band['p_two']:.4g}")

# ---------------------------------------------------------------- matched null
print(f"\n=== {N_SETS} expression-matched random {len(present)}-gene sets ===")
mod_bins = [bins.iloc[i] for i in gidx.values()]
pool = {b: np.flatnonzero((bins == b).to_numpy() & ok) for b in set(pd.Series(mod_bins).dropna())}
rng = np.random.default_rng(0)
exclude = set(gidx.values())
null_b, null_p = [], []
for _ in range(N_SETS):
    idx = []
    for b in mod_bins:
        if pd.isna(b):
            continue
        c = pool[b]
        c = c[~np.isin(c, list(exclude))]
        idx.append(int(rng.choice(c)))
    obs["sc_null"] = score(idx)
    tt = donor_table()
    tt = tt.rename(columns={"sc": "sc_real"})
    tt["sc"] = donor_table().sc  # placeholder, replaced below
    g = obs.groupby("donor").sc_null.median()
    tt["sc"] = tt.donor.map(g)
    r = test(tt)
    null_b.append(r["beta"]); null_p.append(r["p_one"])
null_b, null_p = np.array(null_b), np.array(null_p)
emp = float((null_b >= r_all["beta"]).mean()) if r_all["beta"] > 0 else float((null_b <= r_all["beta"]).mean())
print(f"observed beta {r_all['beta']:+.4f}; null beta median {np.median(null_b):+.4f}, "
      f"95th pct {np.percentile(null_b, 95):+.4f}")
print(f"empirical p (one-sided vs matched null) = {emp:.4f}")
print(f"fraction of random sets with p_one < 0.05 = {(null_p < 0.05).mean():.3f}  "
      f"(nominal 0.05)")

pd.DataFrame(dict(null_beta=null_b, null_p=null_p)).to_csv(OUT / "cilium28_null.csv", index=False)
T.to_csv(OUT / "cilium28_donor_scores.csv", index=False)
print("\nwrote", OUT)
