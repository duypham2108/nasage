"""NasAge axes vs age in older-adult nasal brushings (GSE336399).

Implements docs/PREREG_gse336399_age.md, locked at commit 04f30ef before any expression
value was read. Nothing here deviates from it.

Order of work (plan section 8): expression gate and family size are printed BEFORE any
p-value; then the primary models, the matched-null control, the named cilium secondary,
and the sensitivity analyses.

Usage: python validation/gse336399_age.py <counts.tsv.gz> <gsm_text> <out_dir>
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.formula.api as smf

import nasage as nas

COUNTS, GSM, OUT = sys.argv[1], sys.argv[2], Path(sys.argv[3])
OUT.mkdir(parents=True, exist_ok=True)
N_SETS, N_BINS = 500, 25
ALIAS = {"CIMAP3": "ODF3B", "STING1": "TMEM173"}
CILIUM28 = ["SEPTIN2", "DNAH3", "SEPTIN7", "ARMC2", "CFAP61", "TTC29", "SEPTIN9", "CFAP206",
            "B9D1", "NEK10", "CC2D2A", "LRGUK", "DNAH6", "TTC12", "KIF27", "BBS9", "CFAP298",
            "DMD", "IFT57", "ARMC9", "DRC1", "ARL3", "CCDC39", "CFAP47", "BBS5", "PIBF1",
            "ROPN1L", "PARVA"]

# ---------------------------------------------------------------- metadata
recs, cur = [], None
for line in open(GSM, encoding="utf-8", errors="replace"):
    line = line.rstrip()
    if line.startswith("^SAMPLE"):
        cur = {"gsm": line.split("=")[1].strip()}
        recs.append(cur)
    elif cur is None:
        continue
    elif line.startswith("!Sample_title"):
        cur["sample"] = line.split("=", 1)[1].strip()
    elif line.startswith("!Sample_characteristics_ch1"):
        v = line.split("=", 1)[1].strip()
        if ":" in v:
            k, val = v.split(":", 1)
            cur[k.strip().lower()] = val.strip()
meta = pd.DataFrame(recs).set_index("sample")
meta = meta.rename(columns={"cancer status": "cancer", "smoking status": "smoking",
                            "median tin": "tin", "sequencing batch": "batch"})
meta["age"] = pd.to_numeric(meta.age, errors="coerce")
meta["tin"] = pd.to_numeric(meta.tin, errors="coerce")

C = pd.read_csv(COUNTS, sep="\t", index_col=0)
C = C.loc[:, [c for c in C.columns if c in meta.index]]
meta = meta.loc[C.columns]
C = C.groupby(level=0).sum()
print(f"{C.shape[1]} samples x {C.shape[0]} features; metadata matched for all")
meta["depth"] = C.sum(0)
meta["log10_depth"] = np.log10(meta.depth)
print("age:", meta.age.describe().round(1).to_dict())
print("cancer:", meta.cancer.value_counts().to_dict(),
      "| smoking:", meta.smoking.value_counts().to_dict())

# ---------------------------------------------------------------- scoring (plan 2)
cpm = C / C.sum(0) * 1e6
logc = np.log2(cpm + 1)
mu, sd = logc.mean(1), logc.std(1, ddof=1)
expressed = (cpm.mean(1) >= 1) & (sd > 0)
Z = logc[expressed].sub(mu[expressed], axis=0).div(sd[expressed], axis=0)
print(f"\nfeatures passing CPM >= 1: {int(expressed.sum())} of {len(expressed)}")

M = {m: [ALIAS.get(g, g) for g in d["genes"]] for m, d in nas.load_modules().items()}
M["CILIUM28_manuscript"] = CILIUM28
rows, used = [], {}
for m, genes in M.items():
    g = [x for x in genes if x in Z.index]
    rows.append(dict(module=m, n_genes=len(genes), n_used=len(g),
                     dropped=";".join(sorted(set(genes) - set(g)))))
    if g:
        used[m] = g
cover = pd.DataFrame(rows)
cover.to_csv(OUT / "gene_coverage.csv", index=False)
print("\ngene coverage:")
print(cover.to_string(index=False))

# single-gene and ratio features are kept only where the package defines them as such;
# everything with at least 3 usable genes enters the family
family = [m for m, g in used.items() if len(g) >= 3 and m != "CILIUM28_manuscript"]
print(f"\n=== FAMILY SIZE, fixed before any p-value ===")
print(f"modules in the primary family: {len(family)}")
print(f"smallest attainable BH = p_min x {len(family)}; BH < 0.05 needs "
      f"p < {0.05 / len(family):.2e}")
print("attainable: the model p is not floored, so the plan proceeds.")

S = pd.DataFrame({m: Z.loc[g].mean(0) for m, g in used.items()})
scores = meta.join(S)
scores.to_csv(OUT / "module_scores.csv")

# ---------------------------------------------------------------- models
FORMULA = ("y ~ age + C(gender) + C(cancer) + C(smoking) + tin + log10_depth + C(batch)")


def fit(df, col, formula=FORMULA):
    d = df.copy()
    d["y"] = d[col]
    d = d.dropna(subset=["y", "age", "tin", "log10_depth"])
    m = smf.ols(formula, data=d).fit()
    r = stats.spearmanr(d[col], d.age)
    return dict(score=col, n=len(d), beta_age=m.params["age"], p_adj=m.pvalues["age"],
                rho=r.statistic, p_rho=r.pvalue)


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = np.empty_like(p)
    q[o] = np.minimum.accumulate((p[o] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    return np.minimum(q, 1.0)


P = pd.DataFrame([fit(scores, m) for m in family])
P["BH_adj"] = bh(P.p_adj)
P["BH_rho"] = bh(P.p_rho)
pd.set_option("display.width", 220)
print("\n" + "=" * 78)
print("PRIMARY: age gradient, two-sided, BH within the family")
print("=" * 78)
print(P.sort_values("p_adj").round(4).to_string(index=False))
P.to_csv(OUT / "primary_age_effects.csv", index=False)

# ---------------------------------------------------------------- matched null (plan 2.1)
s = logc[expressed].mean(1)
bins = pd.qcut(s.rank(method="first"), N_BINS, labels=False)
rng = np.random.default_rng(0)
pool = {b: np.array(idx) for b, idx in bins.groupby(bins).groups.items()}
null_rows = []
for m in family + ["CILIUM28_manuscript"]:
    g = used[m]
    gb = [bins.get(x) for x in g if x in bins.index]
    obs = fit(scores, m)
    betas, ps = [], []
    for _ in range(N_SETS):
        pick = [str(rng.choice(pool[b])) for b in gb if pd.notna(b)]
        col = f"__null_{m}"
        scores[col] = Z.loc[pick].mean(0)
        r = fit(scores, col)
        betas.append(r["beta_age"]); ps.append(r["p_adj"])
        scores.drop(columns=[col], inplace=True)
    betas, ps = np.array(betas), np.array(ps)
    emp = float((np.abs(betas) >= abs(obs["beta_age"])).mean())
    null_rows.append(dict(module=m, beta_age=obs["beta_age"], p_adj=obs["p_adj"],
                          null_beta_med=float(np.median(np.abs(betas))),
                          null_beta_p95=float(np.percentile(np.abs(betas), 95)),
                          empirical_p=emp, null_frac_p05=float((ps < 0.05).mean())))
N = pd.DataFrame(null_rows)
N.to_csv(OUT / "matched_null_control.csv", index=False)
print("\n" + "=" * 78)
print(f"SPECIFICITY: {N_SETS} expression-matched random gene sets per module")
print("=" * 78)
print(N.sort_values("empirical_p").round(4).to_string(index=False))

# ---------------------------------------------------------------- cilium secondary
print("\n" + "=" * 78)
print("NAMED SECONDARY: the manuscript's 28 cilium genes")
print("=" * 78)
c = fit(scores, "CILIUM28_manuscript")
cn = N[N.module == "CILIUM28_manuscript"].iloc[0]
print(f"n = {c['n']}, genes used {len(used['CILIUM28_manuscript'])}/28")
print(f"beta_age = {c['beta_age']:+.5f}  p(adjusted, two-sided) = {c['p_adj']:.4g}")
print(f"Spearman rho = {c['rho']:+.3f}  p = {c['p_rho']:.4g}")
print(f"empirical p vs matched null = {cn.empirical_p:.4f}   "
      f"(random sets reaching p<0.05: {cn.null_frac_p05:.3f})")

# ---------------------------------------------------------------- sensitivity (plan 5)
print("\n" + "=" * 78)
print("SENSITIVITY (direction and effect size only where underpowered)")
print("=" * 78)
for label, sel, form in (
        ("benign only (n=40)", scores.cancer == "Benign",
         "y ~ age + C(gender) + C(smoking) + tin + log10_depth"),
        ("never-smokers (n=50)", scores.smoking == "NEVER",
         "y ~ age + C(gender) + C(cancer) + tin + log10_depth")):
    sub = scores[sel]
    rows = [fit(sub, m, form) for m in family + ["CILIUM28_manuscript"]]
    T = pd.DataFrame(rows).sort_values("p_adj")
    print(f"\n-- {label}, n = {len(sub)}")
    print(T.head(6).round(4).to_string(index=False))
    T.to_csv(OUT / f"sensitivity_{label.split()[0]}.csv", index=False)

print("\nwrote", OUT)
