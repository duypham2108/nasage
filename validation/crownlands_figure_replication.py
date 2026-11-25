"""Implements docs/PREREG_crownlands_figure_replication.md (locked at commit 2cbb1d2).

H1  basal fraction falls with age                      (full obs table)
H2  CDKN2A+ fraction of ciliated cells rises with age   (full ciliated set)
H3  M02 axoneme lower in CDKN2A+ ciliated, within donor (fixed subsample)
H4  M11 higher in CDKN2A+ ciliated, within donor        (fixed subsample)
H5  basal fraction vs CDKN2A+ fraction, adjusted for age
S1  M05 rises 40-64 -> 65+ (secondary; reported, outside the BH family)

Pass 1 scans the full atlas once for per-cell CDKN2A counts, totals and detected-gene
counts in ciliated cells, and caches them; later runs reuse the cache.

Usage:
  python validation/crownlands_figure_replication.py <atlas.h5ad> <subsample.h5ad> \
      <cell_scores_obs.csv.gz> <out_dir>
"""
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.formula.api as smf

ATLAS, SUB, SCORES, OUT = sys.argv[1], sys.argv[2], sys.argv[3], Path(sys.argv[4])
OUT.mkdir(parents=True, exist_ok=True)

RESP_EPI = {                                   # H1 denominator, fixed by plan section 3
    "respiratory basal cell", "respiratory tract suprabasal cell",
    "multiciliated epithelial cell", "mucus secreting cell", "secretory cell",
    "goblet cell", "club cell", "deuterosomal cell", "ionocyte", "tuft cell",
}
BASAL = {"respiratory basal cell"}
SUPRA = {"respiratory tract suprabasal cell"}
CILIATED = "multiciliated epithelial cell"
MIN_DENOM, MIN_SIDE = 200, 20


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = np.empty_like(p)
    q[o] = np.minimum.accumulate((p[o] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    return np.minimum(q, 1.0)


def one_sided(rho, p_two, want):
    """Convert a two-sided p to the one-sided p for the locked direction."""
    return p_two / 2 if ((rho < 0) if want == "neg" else (rho > 0)) else 1 - p_two / 2


# ---------------------------------------------------------------- obs of the full atlas
f = h5py.File(ATLAS, "r")


def col(name):
    g = f["obs"][name]
    if isinstance(g, h5py.Group):
        cats = np.array([x.decode() if isinstance(x, bytes) else x
                         for x in g["categories"][:]])
        return cats[g["codes"][:]]
    v = g[:]
    return np.array([x.decode() if isinstance(x, bytes) else x for x in v])


obs = pd.DataFrame({"donor": col("donor_id"), "disease": col("disease"),
                    "cell_type": col("cell_type"), "sex": col("sex")})
obs["age"] = pd.Series(col("development_stage")).str.extract(r"(\d+)-year-old")[0].astype(float).values
healthy = (obs.disease == "normal").values
print(f"{len(obs):,} cells; healthy {healthy.sum():,}; donors {obs.donor[healthy].nunique()}")

# ---------------------------------------------------------------- pass 1 (cached)
cache = OUT / "ciliated_cdkn2a.csv.gz"
if cache.exists():
    cil = pd.read_csv(cache)
    print(f"reusing cached ciliated scan: {len(cil):,} cells")
else:
    sym_g = f["var"]["feature_name"]
    symcats = np.array([x.decode() if isinstance(x, bytes) else x
                        for x in sym_g["categories"][:]])
    sym = symcats[sym_g["codes"][:]]
    j = int(np.flatnonzero(sym == "CDKN2A")[0])
    print(f"CDKN2A at column {j}")

    want = healthy & (obs.cell_type == CILIATED).values
    X = f["X"]
    indptr = X["indptr"][:].astype(np.int64)
    rows_out, cd, tot, ng = [], [], [], []
    CH = 50000
    for s in range(0, len(want), CH):
        e = min(s + CH, len(want))
        m = want[s:e]
        if not m.any():
            continue
        lo, hi = indptr[s], indptr[e]
        dat, ind = X["data"][lo:hi], X["indices"][lo:hi]
        lens = np.diff(indptr[s:e + 1])
        starts = np.concatenate([[0], np.cumsum(lens)])
        for r in np.flatnonzero(m):
            a, b = starts[r], starts[r + 1]
            ii, dd = ind[a:b], dat[a:b]
            k = np.searchsorted(ii, j)
            rows_out.append(s + r)
            cd.append(float(dd[k]) if k < len(ii) and ii[k] == j else 0.0)
            tot.append(float(dd.sum()))
            ng.append(int(b - a))
        print(f"  {e:,}/{len(want):,}", flush=True)
    cil = pd.DataFrame({"row": rows_out, "donor": obs.donor.values[rows_out],
                        "age": obs.age.values[rows_out], "sex": obs.sex.values[rows_out],
                        "cdkn2a": cd, "total": tot, "n_genes": ng})
    cil.to_csv(cache, index=False)
    print(f"wrote {cache}: {len(cil):,} ciliated cells")
f.close()

# ---------------------------------------------------------------- donor tables
d = obs[healthy].copy()
d["is_resp"] = d.cell_type.isin(RESP_EPI)
d["is_basal"] = d.cell_type.isin(BASAL)
d["is_basal_supra"] = d.cell_type.isin(BASAL | SUPRA)
g = d[d.is_resp].groupby("donor")
comp = pd.DataFrame({
    "n_resp": g.size(),
    "basal": g.is_basal.sum(),
    "basal_supra": g.is_basal_supra.sum(),
    "age": g.age.first(), "sex": g.sex.first()}).reset_index()
comp["frac_basal"] = comp.basal / comp.n_resp
comp["frac_basal_supra"] = comp.basal_supra / comp.n_resp
comp = comp[comp.n_resp >= MIN_DENOM]
comp_age = comp.dropna(subset=["age"])
print(f"\nH1 donors with >= {MIN_DENOM} respiratory epithelial cells: {len(comp)}; "
      f"with age: {len(comp_age)}")

for thr, tag in ((1, ""), (2, "_t2")):
    cil[f"pos{tag}"] = cil.cdkn2a >= thr
gc = cil.groupby("donor")
p16 = pd.DataFrame({"n_cil": gc.size(), "n_pos": gc.pos.sum(), "n_pos_t2": gc.pos_t2.sum(),
                    "age": gc.age.first(), "sex": gc.sex.first(),
                    "med_genes": gc.n_genes.median()}).reset_index()
p16["frac_pos"] = p16.n_pos / p16.n_cil
p16["frac_pos_t2"] = p16.n_pos_t2 / p16.n_cil
p16 = p16[p16.n_cil >= MIN_DENOM]
p16_age = p16.dropna(subset=["age"])
print(f"H2 donors with >= {MIN_DENOM} ciliated cells: {len(p16)}; with age: {len(p16_age)}")
comp.to_csv(OUT / "donor_composition.csv", index=False)
p16.to_csv(OUT / "donor_cdkn2a.csv", index=False)

res = []

# ---- H1
r = stats.spearmanr(comp_age.age, comp_age.frac_basal)
res.append(dict(id="H1", claim="basal fraction falls with age", n=len(comp_age),
                stat=r.statistic, p=one_sided(r.statistic, r.pvalue, "neg")))
m1 = smf.ols("np.log(frac_basal/(1-frac_basal)) ~ age + sex", data=comp_age).fit()

# ---- H2
r2 = stats.spearmanr(p16_age.age, p16_age.frac_pos)
res.append(dict(id="H2", claim="CDKN2A+ ciliated fraction rises with age", n=len(p16_age),
                stat=r2.statistic, p=one_sided(r2.statistic, r2.pvalue, "pos")))
m2 = smf.ols("np.log((frac_pos+1e-6)/(1-frac_pos+1e-6)) ~ age + sex + np.log(med_genes)",
             data=p16_age).fit()

# ---- H3, H4 (within-donor, on the fixed subsample)
sub = h5py.File(SUB, "r")
sy = np.array([x.decode() if isinstance(x, bytes) else x for x in sub["var_symbol"][:]])
jj = int(np.flatnonzero(sy == "CDKN2A")[0])
sh = tuple(int(x) for x in sub["X"].attrs["shape"])
ip = sub["X/indptr"][:].astype(np.int64)
ct_s = np.array([x.decode() for x in sub["obs/ct"][:]])
cel = np.flatnonzero(ct_s == "ciliated")
cd_s = np.zeros(sh[0])
for s in range(0, sh[0], 50000):
    e = min(s + 50000, sh[0])
    lo, hi = ip[s], ip[e]
    ind, dat = sub["X/indices"][lo:hi], sub["X/data"][lo:hi]
    lens = np.diff(ip[s:e + 1]); st = np.concatenate([[0], np.cumsum(lens)])
    for r_ in range(e - s):
        a, b = st[r_], st[r_ + 1]
        k = np.searchsorted(ind[a:b], jj)
        if k < b - a and ind[a:b][k] == jj:
            cd_s[s + r_] = dat[a:b][k]
sub.close()

S = pd.read_csv(SCORES)
S["cdkn2a"] = cd_s
S["pos"] = S.cdkn2a >= 1
cilS = S[S.ct == "ciliated"]
pairs = []
for dn, grp in cilS.groupby("donor"):
    a, b = grp[grp.pos], grp[~grp.pos]
    if len(a) < MIN_SIDE or len(b) < MIN_SIDE:
        continue
    pairs.append(dict(donor=dn,
                      m02=a.M02_mucociliary_axoneme_UCell.median() - b.M02_mucociliary_axoneme_UCell.median(),
                      m11=a.M11_senescence_arrest_UCell.median() - b.M11_senescence_arrest_UCell.median()))
P = pd.DataFrame(pairs)
P.to_csv(OUT / "within_donor_cdkn2a_contrast.csv", index=False)
print(f"H3/H4 donors with >= {MIN_SIDE} cells on both sides: {len(P)}")
w3 = stats.wilcoxon(P.m02, alternative="less")
w4 = stats.wilcoxon(P.m11, alternative="greater")
res.append(dict(id="H3", claim="M02 axoneme lower in CDKN2A+ ciliated", n=len(P),
                stat=P.m02.median() / P.m02.std(ddof=1), p=w3.pvalue))
res.append(dict(id="H4", claim="M11 higher in CDKN2A+ ciliated", n=len(P),
                stat=P.m11.median() / P.m11.std(ddof=1), p=w4.pvalue))

# ---- H5
J = comp_age.merge(p16_age[["donor", "frac_pos", "med_genes"]], on="donor")
r5 = stats.spearmanr(J.frac_basal, J.frac_pos)
m5 = smf.ols("frac_basal ~ frac_pos + age + sex", data=J).fit()
t5, p5two = m5.tvalues["frac_pos"], m5.pvalues["frac_pos"]
res.append(dict(id="H5", claim="basal vs CDKN2A+ ciliated, adjusted for age", n=len(J),
                stat=m5.params["frac_pos"], p=p5two / 2 if t5 < 0 else 1 - p5two / 2))

R = pd.DataFrame(res)
R["BH"] = bh(R.p)
R.to_csv(OUT / "primary_family.csv", index=False)

pd.set_option("display.width", 220)
print("\n" + "=" * 78)
print("PRIMARY FAMILY (one-sided, directions from the manuscript; BH over 5)")
print("=" * 78)
print(R.round(5).to_string(index=False))

print("\n--- confirmatory models ---")
print(f"H1 logit(basal) ~ age + sex : beta_age={m1.params['age']:+.5f} "
      f"p_two={m1.pvalues['age']:.4g}")
print(f"H2 logit(CDKN2A+) ~ age + sex + log(genes): beta_age={m2.params['age']:+.5f} "
      f"p_two={m2.pvalues['age']:.4g}  beta_genes={m2.params['np.log(med_genes)']:+.4f} "
      f"p={m2.pvalues['np.log(med_genes)']:.3g}")
print(f"H5 unadjusted Spearman rho={r5.statistic:+.3f} p={r5.pvalue:.4g} | "
      f"adjusted beta={m5.params['frac_pos']:+.5f} p_two={p5two:.4g}")

print("\n--- the depth confound on H2 (plan section 4.1) ---")
rr = stats.spearmanr(p16_age.age, p16_age.med_genes)
print(f"age vs median detected genes per cell: rho={rr.statistic:+.3f} p={rr.pvalue:.4g}")
rr2 = stats.spearmanr(p16_age.med_genes, p16_age.frac_pos)
print(f"median genes vs CDKN2A+ fraction:      rho={rr2.statistic:+.3f} p={rr2.pvalue:.4g}")

print("\n--- threshold sensitivity (CDKN2A >= 2) ---")
r2b = stats.spearmanr(p16_age.age, p16_age.frac_pos_t2)
print(f"H2 at >=2 counts: rho={r2b.statistic:+.3f} "
      f"p_one={one_sided(r2b.statistic, r2b.pvalue, 'pos'):.4g}")
print("\n--- H1 sensitivity (basal + suprabasal) ---")
r1b = stats.spearmanr(comp_age.age, comp_age.frac_basal_supra)
print(f"rho={r1b.statistic:+.3f} p_one={one_sided(r1b.statistic, r1b.pvalue, 'neg'):.4g}")
print("\nwrote", OUT)
