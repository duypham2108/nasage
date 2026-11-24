"""Regenerate the GSE336399 matched-null distributions, keeping every draw.

`gse336399_age.py` stored only the summary of its 500-set control (median, p95,
empirical p, fraction at p<0.05). Figure 5 panel b2 needs the distribution itself, so
this re-runs the identical loop — same module order, same seed, same bins, same model —
and writes all 500 coefficients and p-values per module.

This is a re-draw of a control that has already been reported, not a new analysis. The
summary statistics it prints must match `results/gse336399/matched_null_control.csv`;
the script asserts that for M02 and CILIUM28 and fails loudly if they drift.

Usage: python validation/figure5_data.py <counts.tsv.gz> <module_scores.csv> <symbol_map.h5ad> <out_dir>
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.formula.api as smf

import nasage as nas

COUNTS, SCORES, SYMMAP, OUT = sys.argv[1], sys.argv[2], sys.argv[3], Path(sys.argv[4])
OUT.mkdir(parents=True, exist_ok=True)
N_SETS, N_BINS = 500, 25
ALIAS = {"CIMAP3": "ODF3B", "STING1": "TMEM173"}
CILIUM28 = ["SEPTIN2", "DNAH3", "SEPTIN7", "ARMC2", "CFAP61", "TTC29", "SEPTIN9", "CFAP206",
            "B9D1", "NEK10", "CC2D2A", "LRGUK", "DNAH6", "TTC12", "KIF27", "BBS9", "CFAP298",
            "DMD", "IFT57", "ARMC9", "DRC1", "ARL3", "CCDC39", "CFAP47", "BBS5", "PIBF1",
            "ROPN1L", "PARVA"]

meta = pd.read_csv(SCORES, index_col=0)
meta = meta[["age", "gender", "cancer", "smoking", "batch", "tin", "log10_depth"]]

counts = pd.read_csv(COUNTS, sep="\t", index_col=0)
counts = counts.loc[:, [c for c in counts.columns if c in meta.index]]
meta = meta.loc[counts.columns]
counts.index = counts.index.astype(str).str.split(".").str[0]
import anndata as ad
vmap = ad.read_h5ad(SYMMAP, backed="r").var["feature_name"].astype(str)
sym = vmap.reindex(counts.index)
counts = counts[sym.notna()].groupby(sym[sym.notna()].values).sum()

cpm = counts / counts.sum(0) * 1e6
logc = np.log2(cpm + 1)
mu, sd = logc.mean(1), logc.std(1, ddof=1)
expressed = (cpm.mean(1) >= 1) & (sd > 0)
Z = logc[expressed].sub(mu[expressed], axis=0).div(sd[expressed], axis=0)
del counts, cpm  # patsy resolves C() in this namespace; keep it clean
print(f"{Z.shape[1]} samples x {Z.shape[0]} expressed symbols")

M = {m: [ALIAS.get(g, g) for g in d["genes"]] for m, d in nas.load_modules().items()}
M["CILIUM28_manuscript"] = CILIUM28
used = {m: [x for x in g if x in Z.index] for m, g in M.items()}
used = {m: g for m, g in used.items() if g}
family = [m for m, g in used.items() if len(g) >= 3 and m != "CILIUM28_manuscript"]
print(f"family size {len(family)}")

S = pd.DataFrame({m: Z.loc[g].mean(0) for m, g in used.items()})
scores = meta.join(S)
FORMULA = "y ~ age + C(gender) + C(cancer) + C(smoking) + tin + log10_depth + C(batch)"


def fit(df, col):
    d = df.copy()
    d["y"] = d[col]
    d = d.dropna(subset=["y", "age", "tin", "log10_depth"])
    m = smf.ols(FORMULA, data=d).fit()
    return m.params["age"], m.pvalues["age"]


s = logc[expressed].mean(1)
bins = pd.qcut(s.rank(method="first"), N_BINS, labels=False)
rng = np.random.default_rng(0)
pool = {b: np.array(idx) for b, idx in bins.groupby(bins).groups.items()}

draws, summary = [], []
for m in family + ["CILIUM28_manuscript"]:
    gb = [bins.get(x) for x in used[m] if x in bins.index]
    obs_b, obs_p = fit(scores, m)
    betas, ps = [], []
    for _ in range(N_SETS):
        pick = [str(rng.choice(pool[b])) for b in gb if pd.notna(b)]
        scores["__null"] = Z.loc[pick].mean(0)
        b, p = fit(scores, "__null")
        betas.append(b); ps.append(p)
        scores.drop(columns=["__null"], inplace=True)
    betas, ps = np.array(betas), np.array(ps)
    draws.append(pd.DataFrame({"module": m, "beta": betas, "p": ps}))
    summary.append(dict(module=m, beta_age=obs_b, p_adj=obs_p,
                        empirical_p=float((np.abs(betas) >= abs(obs_b)).mean()),
                        null_frac_p05=float((ps < 0.05).mean())))
    print(f"{m:42s} beta {obs_b:+.5f}  emp_p {summary[-1]['empirical_p']:.4f}  "
          f"frac_p05 {summary[-1]['null_frac_p05']:.3f}")

D = pd.concat(draws, ignore_index=True)
D.to_csv(OUT / "gse336399_null_draws.csv.gz", index=False)
Sm = pd.DataFrame(summary)
Sm.to_csv(OUT / "gse336399_null_summary.csv", index=False)

ref = pd.read_csv(Path(SCORES).parent / "matched_null_control.csv").set_index("module")
for m in ("M02_mucociliary_axoneme", "CILIUM28_manuscript"):
    r, g = ref.loc[m], Sm.set_index("module").loc[m]
    assert abs(r.beta_age - g.beta_age) < 1e-10, m
    assert abs(r.empirical_p - g.empirical_p) < 1e-12, (m, r.empirical_p, g.empirical_p)
    assert abs(r.null_frac_p05 - g.null_frac_p05) < 1e-12, m
print("\nreproduced the reported control exactly for M02 and CILIUM28")
print("wrote", OUT)
