"""NasAge modules vs pneumococcal carriage acquisition (GSE261277, Goncalves et al. 2024).

Implements docs/PREREG_vaccine_axes.md section 8.3, locked in commit 192d4c2 before any
expression value from this dataset was read. Sections 8.2 (predicted signs) and 8.3
(model, primary family, multiplicity, power, failure conditions) are the specification;
nothing here deviates from them.

Primary family (one-sided, BH over 3):
    M05 ISG                 higher in carriage POS
    M12 inflammatory SASP   higher in carriage POS
    M01 epithelial barrier  lower  in carriage POS

Usage:
    python validation/pneumo_carriage_gse261277.py <counts.tsv.gz> <gsm_text> <out_dir>
"""
import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import nasage as nas

counts_f, gsm_f, out = sys.argv[1], sys.argv[2], Path(sys.argv[3])
out.mkdir(parents=True, exist_ok=True)

PRIMARY = {                      # module -> locked direction of carriage POS
    "M05_interferon_stimulated_genes": "greater",
    "M12_inflammatory_sasp": "greater",
    "M01_epithelial_barrier": "less",
}

# ---------------------------------------------------------------- metadata
recs, cur = [], None
for line in open(gsm_f, encoding="utf-8", errors="replace"):
    line = line.rstrip()
    if line.startswith("^SAMPLE"):
        cur = {"gsm": line.split("=")[1].strip()}
        recs.append(cur)
    elif cur is None:
        continue
    elif line.startswith("!Sample_title"):
        cur["title"] = line.split("=", 1)[1].strip()
    elif line.startswith("!Sample_characteristics_ch1"):
        v = line.split("=", 1)[1].strip()
        if ":" in v:
            k, val = v.split(":", 1)
            cur[k.strip().lower()] = val.strip()

meta = pd.DataFrame(recs)
meta["barcode"] = meta["title"].str.extract(r"(S\d{8})")
meta = meta.rename(columns={"carriage status": "carriage", "smorkers": "smoker"})
meta["age"] = pd.to_numeric(meta["age"], errors="coerce")
meta = meta.set_index("barcode")[["gsm", "age", "sex", "timepoint", "carriage", "smoker"]]

# ---------------------------------------------------------------- counts
with gzip.open(counts_f, "rt") as fh:
    C = pd.read_csv(fh, sep="\t", index_col=0)
C = C.groupby(level=0).sum()
C = C.loc[:, [c for c in C.columns if c in meta.index]]
meta = meta.loc[C.columns]
print(f"{C.shape[1]} samples x {C.shape[0]} genes; metadata matched for all")

meta["reads"] = C.sum(0)
keep = meta.reads >= 5e6
print(f"depth: min {meta.reads.min():,.0f}  median {meta.reads.median():,.0f}  "
      f"max {meta.reads.max():,.0f}")
print(f"excluded (<5M assigned counts): {int((~keep).sum())} -> "
      f"{sorted(meta.index[~keep])}")
meta, C = meta[keep], C.loc[:, keep]

base = meta.timepoint == "baseline"
print("\nanalysis set after QC:")
print(pd.crosstab(meta.timepoint, meta.carriage).to_string())

# ---------------------------------------------------------------- scoring
# Section 8.3.2: z across the analysed BASELINE samples; mean CPM < 1 at baseline dropped.
cpm = C / C.sum(0) * 1e6
logc = np.log2(cpm + 1)
mu = logc.loc[:, base].mean(1)
sd = logc.loc[:, base].std(1, ddof=1)
expressed = (cpm.loc[:, base].mean(1) >= 1) & (sd > 0)
Z = logc[expressed].sub(mu[expressed], axis=0).div(sd[expressed], axis=0)
print(f"\ngenes passing CPM>=1 at baseline: {int(expressed.sum())} of {len(expressed)}")

M = nas.load_modules()
rows, cover = {}, []
for m, d in M.items():
    used = [g for g in d["genes"] if g in Z.index]
    dropped = sorted(set(d["genes"]) - set(used))
    cover.append(dict(module=m, n_genes=len(d["genes"]), n_used=len(used),
                      dropped=";".join(dropped)))
    if used:
        rows[m] = Z.loc[used].mean(0)
S = pd.DataFrame(rows)
cover = pd.DataFrame(cover)
print("\nmodule gene coverage:")
print(cover.to_string(index=False))

# UCell sensitivity (section 8.3.2): bulk samples ranked like cells
try:
    import anndata as ad
    from pyucell import compute_ucell_scores
    from nasage.scoring import MAX_RANK
    u = ad.AnnData(X=C.T.values.astype(np.float32),
                   obs=pd.DataFrame(index=C.columns), var=pd.DataFrame(index=C.index))
    sigs = {m: [g for g in d["genes"] if g in C.index] for m, d in M.items()
            if len(d["genes"]) > 1}
    compute_ucell_scores(u, sigs, max_rank=MAX_RANK, suffix="_UCell", n_jobs=1)
    S = S.join(u.obs)
    print("UCell sensitivity computed")
except Exception as exc:                                     # pragma: no cover
    print("UCell sensitivity skipped:", exc)

meta["log10_reads"] = np.log10(meta.reads)
scores = meta.join(S)
scores.to_csv(out / "module_scores_all_samples.csv")

# ---------------------------------------------------------------- tests
def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = np.empty_like(p)
    q[o] = np.minimum.accumulate((p[o] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    return np.minimum(q, 1.0)


def cohen_d(a, b):
    na, nb = len(a), len(b)
    s = np.sqrt(((na - 1) * a.var(ddof=1) + (nb - 1) * b.var(ddof=1)) / (na + nb - 2))
    return (a.mean() - b.mean()) / s if s > 0 else np.nan


def contrast(df, col, alt):
    pos = df.loc[df.carriage == "POS", col].dropna()
    neg = df.loc[df.carriage == "NEG", col].dropna()
    u = stats.mannwhitneyu(pos, neg, alternative=alt)
    # adjusted model, section 8.3.3
    import statsmodels.formula.api as smf
    d = df[[col, "carriage", "age", "sex", "smoker", "log10_reads"]].dropna().copy()
    d["y"] = d[col]
    d["pos"] = (d.carriage == "POS").astype(int)
    fit = smf.ols("y ~ pos + age + C(sex) + C(smoker) + log10_reads", data=d).fit()
    t, p2 = fit.tvalues["pos"], fit.pvalues["pos"]
    p1 = p2 / 2 if ((t > 0) == (alt == "greater")) else 1 - p2 / 2
    ci = fit.conf_int().loc["pos"]
    return dict(score=col, n_pos=len(pos), n_neg=len(neg),
                mean_pos=pos.mean(), mean_neg=neg.mean(),
                d=cohen_d(pos, neg), p_rank=u.pvalue,
                beta_adj=fit.params["pos"], ci_lo=ci[0], ci_hi=ci[1], p_adj=p1)


B = scores[base[keep]] if base.index.equals(keep.index) else scores[scores.timepoint == "baseline"]
B = scores[scores.timepoint == "baseline"]

print("\n" + "=" * 78)
print("PRIMARY FAMILY  (one-sided, direction locked in section 8.2; BH over 3)")
print("=" * 78)
prim = pd.DataFrame([contrast(B, m, alt) for m, alt in PRIMARY.items()])
prim["bh_rank"] = bh(prim.p_rank)
prim["bh_adj"] = bh(prim.p_adj)
prim["locked_direction"] = [PRIMARY[m] for m in prim.score]
prim["independence"] = ["replication", "replication", "INDEPENDENT"]
print(prim.round(4).to_string(index=False))
prim.to_csv(out / "primary_family.csv", index=False)

# UCell sensitivity on the same three
uc = [m + "_UCell" for m in PRIMARY if m + "_UCell" in B.columns]
if uc:
    sens = pd.DataFrame([contrast(B, m + "_UCell", alt) for m, alt in PRIMARY.items()
                         if m + "_UCell" in B.columns])
    sens["bh_rank"] = bh(sens.p_rank)
    print("\nUCell sensitivity (same contrasts):")
    print(sens.round(4).to_string(index=False))
    sens.to_csv(out / "primary_family_ucell.csv", index=False)

print("\n" + "=" * 78)
print("EXPLORATORY: all other modules at baseline, two-sided, BH within this family")
print("=" * 78)
other = [c for c in S.columns if not c.endswith("_UCell") and c not in PRIMARY]
ex = pd.DataFrame([contrast(B, m, "two-sided") for m in other])
ex["bh_rank"] = bh(ex.p_rank)
print(ex.sort_values("p_rank").round(4).to_string(index=False))
ex.to_csv(out / "exploratory_baseline.csv", index=False)

print("\n" + "=" * 78)
print("EXPLORATORY: post-exposure carriage contrasts (D2, D9), two-sided")
print("=" * 78)
post = []
for tp in ("D2", "D9"):
    sub = scores[scores.timepoint == tp]
    for m in list(PRIMARY) + other:
        r = contrast(sub, m, "two-sided")
        r["timepoint"] = tp
        post.append(r)
P = pd.DataFrame(post)
P["bh_rank"] = P.groupby("timepoint").p_rank.transform(bh)
print(P[P.score.isin(PRIMARY)].round(4).to_string(index=False))
P.to_csv(out / "exploratory_post_exposure.csv", index=False)

print("\n" + "=" * 78)
print("EXPLORATORY: age gradient within 51-80 at baseline, two-sided")
print("=" * 78)
age = []
for m in list(PRIMARY) + other:
    d = B[[m, "age"]].dropna()
    r = stats.spearmanr(d[m], d.age)
    age.append(dict(score=m, n=len(d), rho=r.statistic, p=r.pvalue))
A = pd.DataFrame(age)
A["bh"] = bh(A.p)
print(A.sort_values("p").round(4).to_string(index=False))
A.to_csv(out / "exploratory_age_gradient.csv", index=False)

print("\nwrote results to", out)
