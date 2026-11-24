"""Analysis B of docs/PREREG_celltype_mvri.md, locked at 77701ef.

Does a live attenuated intranasal influenza vaccine engage the seven curated epithelial
programmes? Paired within-participant change, baseline -> day 3, Thwaites GSE230494.

M05 interferon is the POSITIVE CONTROL, not a finding: it already rose in 17/20 participants
(p = 0.008) in results/thwaites2023_laiv/. If it does not rise here, Analysis B is void.

Usage: python validation/laiv_engagement_programs.py <data_dir> <out_dir>
"""
import gzip
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm

sys.path.insert(0, str(Path(__file__).parent))
from mvri_programs import PROGRAMS, PROXY  # noqa: E402

DATA, OUT = Path(sys.argv[1]), Path(sys.argv[2])
OUT.mkdir(parents=True, exist_ok=True)
EXCLUDE = {"P2_junction_integrity"}
N_SETS, N_BINS, N_PERM = 500, 25, 2000
MIN_READS = 5e6

# M05 interferon-stimulated genes, the positive control (same list as the curated module set).
M05 = ["ISG15", "IFI6", "IFI44L", "IFIT1", "IFIT2", "IFIT3", "MX1", "MX2", "OAS1", "OAS2",
       "OAS3", "OASL", "RSAD2", "IFITM1", "IFITM3", "STAT1", "STAT2", "IRF7", "BST2",
       "XAF1", "HERC5", "USP18"]

# ---------------------------------------------------------------- counts
rows = {}
for f in sorted((DATA / "counts").glob("*_counts.txt.gz")):
    name = re.sub(r"^GSM\d+_|_counts\.txt\.gz$", "", f.name)
    with gzip.open(f, "rt") as fh:
        rows[name] = pd.read_csv(fh, sep="\t", header=None, index_col=0, skiprows=1).iloc[:, 0]
C = pd.DataFrame(rows)
C.index = [str(i).split(".")[0] for i in C.index]

# Ensembl -> symbol, using the same reference var table as validation/laiv_thwaites_bulk.py
import anndata as ad  # noqa: E402
sym = ad.read_h5ad(DATA.parent / "ziegler2021_np.h5ad",
                   backed="r").var["feature_name"].astype(str).reindex(C.index)
print(f"{sym.notna().mean():.1%} of genes mapped to symbols")
C = C[sym.notna()].groupby(sym[sym.notna()].values).sum()

reads = C.sum(0)
C = C.loc[:, reads >= MIN_READS]
print(f"{C.shape[1]} samples >= {MIN_READS/1e6:.0f}M reads, {C.shape[0]} genes")

meta = pd.DataFrame(index=C.columns)
meta["participant"] = [re.match(r"(LAIV_\d+)", s).group(1) for s in C.columns]
# same rule as laiv_thwaites_bulk.py: names ending _1 or _Visit1 are baseline
meta["timepoint"] = np.where(meta.index.str.contains(r"_(?:1|Visit1)$"), "baseline", "SD3")

cpm = C / C.sum(0) * 1e6
logc = np.log2(cpm + 1)
base = meta.index[meta.timepoint == "baseline"]
mu, sd = logc.loc[:, base].mean(1), logc.loc[:, base].std(1, ddof=1)
keep = (cpm.loc[:, base].mean(1) >= 1) & (sd > 0)
Z = logc.loc[keep].sub(mu[keep], axis=0).div(sd[keep], axis=0)
print(f"{keep.sum()} expressed genes retained")


def score(genes):
    g = [x for x in genes if x in Z.index]
    return (Z.loc[g].mean(0) if g else None), g


# ---------------------------------------------------------------- paired frame
S = pd.DataFrame(index=Z.columns)
cover = []
for name, prog in PROGRAMS.items():
    v, g = score(prog["genes"])
    cover.append(dict(program=name, n_genes=len(prog["genes"]), n_used=len(g),
                      frac=len(g) / len(prog["genes"])))
    if v is not None:
        S[name] = v
for ct, genes in PROXY.items():
    v, _ = score(genes)
    if v is not None:
        S[f"proxy_{ct}"] = v
v, g05 = score(M05)
S["M05"] = v
S["PTPRC"] = logc.loc["PTPRC"] if "PTPRC" in logc.index else np.nan
S = S.join(meta)
COV = pd.DataFrame(cover)
COV.to_csv(OUT / "coverage.csv", index=False)
print("\n=== PROGRAMME COVERAGE, before any model ===")
print(COV.round(3).to_string(index=False))
print(f"M05 positive control: {len(g05)}/{len(M05)} genes")

piv = S.pivot_table(index="participant", columns="timepoint")
paired = [p for p in piv.index
          if not pd.isna(piv.loc[p, ("M05", "baseline")])
          and not pd.isna(piv.loc[p, ("M05", "SD3")])]
print(f"\npaired participants: {len(paired)}")


def delta(col):
    return pd.Series({p: piv.loc[p, (col, "SD3")] - piv.loc[p, (col, "baseline")]
                      for p in paired}).dropna()


# ---------------------------------------------------------------- positive control FIRST
d05 = delta("M05")
t05 = stats.ttest_rel(np.zeros(len(d05)), -d05)
print("\n" + "=" * 78)
print("POSITIVE CONTROL (plan 7.4): M05 interferon must rise, or Analysis B is void")
print("=" * 78)
print(f"  n = {len(d05)}, mean change = {d05.mean():+.3f} SD, rose in "
      f"{(d05 > 0).sum()}/{len(d05)}, paired t p = {t05.pvalue:.4f}, "
      f"Wilcoxon p = {stats.wilcoxon(d05).pvalue:.4f}")
if not (d05.mean() > 0 and t05.pvalue < 0.05):
    pd.DataFrame([dict(control="M05", mean=d05.mean(), p=t05.pvalue, n=len(d05),
                       verdict="VOID")]).to_csv(OUT / "positive_control.csv", index=False)
    sys.exit("\nPOSITIVE CONTROL FAILED -> Analysis B is void per plan 7.4. Stopped.")
pd.DataFrame([dict(control="M05", mean=d05.mean(), p=t05.pvalue, n=len(d05),
                   n_rose=int((d05 > 0).sum()), verdict="passed")]
             ).to_csv(OUT / "positive_control.csv", index=False)


def paired_test(col, adjust=None):
    d = delta(col)
    if adjust is None:
        t = stats.ttest_rel(np.zeros(len(d)), -d)
        dz = d.mean() / d.std(ddof=1)
        return d.mean(), dz, float(t.pvalue), len(d), float(stats.wilcoxon(d).pvalue)
    X = pd.DataFrame({"const": 1.0, "adj": delta(adjust).reindex(d.index)}).dropna()
    y = d.reindex(X.index)
    r = sm.OLS(y, X).fit()
    return (float(r.params["const"]), np.nan, float(r.pvalues["const"]), len(y), np.nan)


res = []
for name, prog in PROGRAMS.items():
    if name not in S.columns:
        continue
    pcol = f"proxy_{prog['cell_type']}"
    m, dz, p, n, pw = paired_test(name)
    ma, _, pa, na, _ = paired_test(name, adjust=pcol) if pcol in S.columns else (np.nan,) * 5
    mp, _, pp, _, _ = paired_test(name, adjust="PTPRC")
    res.append(dict(program=name, cell_type=prog["cell_type"], n=n, mean_change=m, d_z=dz,
                    p=p, p_wilcoxon=pw, mean_adj=ma, p_adj=pa, p_ptprc_adj=pp))
R = pd.DataFrame(res)
fam = R[~R.program.isin(EXCLUDE)].copy()
pv = fam.p.values
o = np.argsort(pv); q = np.empty_like(pv)
q[o] = np.minimum.accumulate((pv[o] * len(pv) / np.arange(1, len(pv) + 1))[::-1])[::-1]
fam["BH"] = np.minimum(q, 1.0)
R = R.merge(fam[["program", "BH"]], on="program", how="left")

# ---------------------------------------------------------------- matched null
expr = logc.loc[keep, base].mean(1)
bins = pd.qcut(expr.rank(method="first"), N_BINS, labels=False)
pool = {k: np.array(v) for k, v in bins.groupby(bins).groups.items()}
rng = np.random.default_rng(0)
emp = {}
for name, prog in PROGRAMS.items():
    if name in EXCLUDE or name not in S.columns:
        continue
    gb = [bins.get(g) for g in prog["genes"] if g in bins.index]
    obs = abs(R.loc[R.program == name, "mean_change"].iloc[0])
    draws = []
    for _ in range(N_SETS):
        pick = [str(rng.choice(pool[k])) for k in gb if pd.notna(k)]
        v = Z.loc[pick].mean(0)
        tmp = S[[]].copy()
        tmp["__n"] = v
        tmp = tmp.join(meta)
        pv2 = tmp.pivot_table(index="participant", columns="timepoint", values="__n")
        d = (pv2["SD3"] - pv2["baseline"]).dropna()
        if len(d) >= 3:
            draws.append(abs(d.mean()))
    emp[name] = float(np.mean(np.array(draws) >= obs)) if draws else np.nan
R["empirical_p"] = R.program.map(emp)

# ---------------------------------------------------------------- label permutation
rng2 = np.random.default_rng(1)
perm = {}
for name in R.program:
    if name in EXCLUDE or name not in S.columns:
        continue
    obs = abs(R.loc[R.program == name, "mean_change"].iloc[0])
    cnt = 0
    d = delta(name)
    for _ in range(N_PERM):
        sgn = rng2.choice([-1, 1], size=len(d))
        if abs(float((d.values * sgn).mean())) >= obs:
            cnt += 1
    perm[name] = cnt / N_PERM
R["perm_p"] = R.program.map(perm)


def verdict(r):
    if r.program in EXCLUDE:
        return "excluded a priori"
    if not np.isfinite(r.BH) or r.BH >= 0.05:
        return "not engaged (underpowered below d_z 0.5)"
    if np.isfinite(r.empirical_p) and r.empirical_p >= 0.05:
        return "not_specific (fails matched null)"
    if np.isfinite(r.p_adj) and r.p_adj >= 0.05:
        return "composition, not programme"
    return "ENGAGED"


R["verdict"] = R.apply(verdict, axis=1)
R.to_csv(OUT / "engagement.csv", index=False)
S.to_csv(OUT / "scores.csv")
print("\n=== ANALYSIS B: paired change after LAIV (BH over the family of 6) ===")
print(R.round(4).to_string(index=False))
print("\nwrote", OUT)
