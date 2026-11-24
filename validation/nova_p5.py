"""Replication of the P5 mucin-barrier result in the NOVA influenza challenge.

Implements docs/PREREG_nova_p5.md, locked at commit 320a5cd before any expression value
from this cohort was linked to its outcome.

Expression source, and why it is the raw arrays: both processed matrices GEO offers are
expressed as a RATIO TO EACH PARTICIPANT'S OWN SD -1 BASELINE, so every baseline value is
exactly 1 (91.7% of baseline cells are 1.0; 0% of SD7 cells are). Baseline is the
denominator and cannot be recovered from them. Per-sample Agilent Feature Extraction files
are therefore used, quantile-normalised across the baseline arrays here.

Usage: python validation/nova_p5.py <gse253686_dir> <out_dir>
"""
import gzip
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

from mvri_programs import PROGRAMS, PROXY

warnings.filterwarnings("ignore")
DATA, OUT = Path(sys.argv[1]), Path(sys.argv[2])
OUT.mkdir(parents=True, exist_ok=True)
N_SETS, N_BINS, N_PERM, MIN_COVER = 500, 25, 2000, 0.60
PRIMARY = "P5_mucin_barrier"
EXCLUDE = {"P2_junction_integrity"}          # plan 4.1: Jaccard 0.824 with M01, dropped a priori


def read_agilent(path):
    """ProbeName -> gProcessedSignal from an Agilent Feature Extraction file."""
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        hdr = None
        for line in fh:
            if line.startswith("FEATURES"):
                hdr = line.rstrip("\n").split("\t")
                break
        if hdr is None:
            raise ValueError(f"no FEATURES block in {path}")
        pi, si = hdr.index("ProbeName"), hdr.index("gProcessedSignal")
        ctrl = hdr.index("ControlType") if "ControlType" in hdr else None
        vals = {}
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) <= max(pi, si):
                continue
            if ctrl is not None and p[ctrl] not in ("0", "0.0"):
                continue
            try:
                v = float(p[si])
            except ValueError:
                continue
            vals.setdefault(p[pi], []).append(v)
    return pd.Series({k: np.mean(v) for k, v in vals.items()})


def quantile_normalise(M):
    """Columns = arrays. Classic quantile normalisation."""
    ranks = M.rank(method="first")
    mean_sorted = pd.Series(np.sort(M.values, axis=0).mean(axis=1), index=range(1, len(M) + 1))
    return ranks.apply(lambda col: col.map(mean_sorted))


# ------------------------------------------------------------------ load arrays
files = pd.read_csv(DATA / "baseline_rawfiles.csv")
cols = {}
for _, r in files.iterrows():
    p = DATA / "raw_baseline" / r["Name"]
    if not p.exists() or p.stat().st_size < 1000:
        print(f"  missing/short, skipped: {r['gsm']}")
        continue
    cols[r["subject"]] = read_agilent(p)
M = pd.DataFrame(cols).dropna(how="any")
print(f"arrays loaded: {M.shape[1]} | probes common to all: {M.shape[0]}")

M = quantile_normalise(M)
L = np.log2(M + 1)

# ------------------------------------------------------------------ probe -> symbol
ann = {}
with gzip.open(DATA / "GPL6480.annot.gz", "rt", encoding="utf-8", errors="replace") as fh:
    started = False
    for line in fh:
        if line.startswith("!platform_table_begin"):
            started = True
            hdr = next(fh).rstrip("\n").split("\t")
            gi = hdr.index("Gene symbol")
            continue
        if not started or line.startswith("!platform_table_end"):
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) > gi and p[gi]:
            ann[p[0]] = p[gi].split("///")[0].strip()
sym = pd.Series(ann).reindex(L.index)
L = L[sym.notna().values]
sym = sym[sym.notna()]
# collapse probes per gene by the probe with the highest mean signal
order = L.mean(1).groupby(sym.values).rank(ascending=False, method="first")
E = L[(order == 1).values]
E.index = sym[(order == 1).values].values
print(f"genes after collapsing probes: {E.shape[0]}")

Z = E.sub(E.mean(1), axis=0).div(E.std(1, ddof=1), axis=0)
Z = Z.dropna(how="any")

# ------------------------------------------------------------------ coverage, before any model
cov = []
for name, p in PROGRAMS.items():
    g = [x for x in p["genes"] if x in Z.index]
    cov.append(dict(program=name, cell_type=p["cell_type"], n_genes=len(p["genes"]),
                    n_used=len(g), frac=len(g) / len(p["genes"]),
                    transferable=len(g) / len(p["genes"]) >= MIN_COVER,
                    dropped=";".join(sorted(set(p["genes"]) - set(g)))))
COV = pd.DataFrame(cov)
COV.to_csv(OUT / "coverage.csv", index=False)
print("\n=== P5 GENE COVERAGE, reported before any model (plan 7.2) ===")
print(COV[["program", "n_used", "n_genes", "frac", "transferable"]].to_string(index=False))
if not COV.set_index("program").loc[PRIMARY, "transferable"]:
    raise SystemExit("P5 below the 60% coverage gate: per plan section 5 the test is NOT run.")

scores = {}
for name, p in PROGRAMS.items():
    g = [x for x in p["genes"] if x in Z.index]
    if len(g) / len(p["genes"]) >= MIN_COVER:
        s = Z.loc[g].mean(0)
        scores[name] = (s - s.mean()) / s.std(ddof=1)
used = {n: [x for x in PROGRAMS[n]["genes"] if x in Z.index] for n in scores}
for ct, markers in PROXY.items():
    g = [x for x in markers if x in Z.index]
    if g:
        s = Z.loc[g].mean(0)
        scores[f"proxy_{ct}"] = (s - s.mean()) / s.std(ddof=1)
S = pd.DataFrame(scores)

# ------------------------------------------------------------------ outcome
s1 = pd.read_csv(DATA / "suppl" / "S1_parsed.csv")
s1["subject"] = s1["Subject ID"].astype(str).str.strip()
s1 = s1.set_index("subject")
D = S.join(s1[["Age", "Sex", "Study Group", "Shedding (BINARY)",
               "Multimodal shedding class"]], how="inner")
D = D.rename(columns={"Shedding (BINARY)": "shed", "Study Group": "arm",
                      "Multimodal shedding class": "mmclass"})
D["shed"] = pd.to_numeric(D.shed, errors="coerce")
D["Age"] = pd.to_numeric(D.Age, errors="coerce")
D["vax"] = (D.arm == "Vax").astype(int)
D["male"] = (D.Sex == "Male").astype(int)
D = D.dropna(subset=["shed", "Age"])
print(f"\nlinked participants: {len(D)}  shedders {int(D.shed.sum())} / "
      f"non-shedders {int((1 - D.shed).sum())}")
print("mmclass:", D.mmclass.value_counts().to_dict())
D.to_csv(OUT / "nova_scores_outcome.csv")


def logit_or(df, col, proxy_col, extra=()):
    cols = [col, "Age", "male", "vax"] + ([proxy_col] if proxy_col else []) + list(extra)
    d = df.dropna(subset=cols + ["shed"])
    X = sm.add_constant(d[cols].astype(float).values, has_constant="add")
    try:
        r = sm.Logit(d["shed"].values.astype(float), X).fit(disp=0)
        b, se, p2 = float(r.params[1]), float(r.bse[1]), float(r.pvalues[1])
        if not np.isfinite(se) or se > 10:
            return (np.nan,) * 5 + (len(d),)
        lo, hi = b - 1.96 * se, b + 1.96 * se
        p1 = stats.norm.cdf(b / se)          # one-sided, protective direction (plan 2)
        return b, se, p2, p1, (np.exp(lo), np.exp(hi)), len(d)
    except Exception:
        return (np.nan,) * 5 + (len(df),)


# ------------------------------------------------------------------ PRIMARY
pcol = f"proxy_{PROGRAMS[PRIMARY]['cell_type']}"
b, se, p2, p1, ci, n = logit_or(D, PRIMARY, pcol)
bu, _, p2u, p1u, ciu, _ = logit_or(D, PRIMARY, None)
print("\n" + "=" * 78)
print("PRIMARY (plan 2): baseline P5 vs shedding, one-sided protective, adjusted")
print("=" * 78)
print(f"  OR per SD = {np.exp(b):.3f}   95% CI {ci[0]:.3f}-{ci[1]:.3f}   n = {n}")
print(f"  one-sided p = {p1:.4f}    (two-sided {p2:.4f})")
print(f"  unadjusted: OR {np.exp(bu):.3f}, one-sided p {p1u:.4f}")

# ------------------------------------------------------------------ matched null
expr_mean = E.mean(1)
expr_mean = expr_mean[~expr_mean.index.duplicated()]
bins = pd.qcut(expr_mean.reindex(Z.index).rank(method="first"), N_BINS, labels=False)
pool = {k: np.array(v) for k, v in bins.groupby(bins).groups.items()}
rng = np.random.default_rng(0)
gb = [bins.get(x) for x in used[PRIMARY]]
nb, npv = [], []
for _ in range(N_SETS):
    pick = [str(rng.choice(pool[k])) for k in gb if pd.notna(k)]
    s = Z.loc[pick].mean(0)
    d2 = D.copy()
    d2["__null"] = ((s - s.mean()) / s.std(ddof=1)).reindex(d2.index).values
    bb, _, pp2, _, _, _ = logit_or(d2, "__null", pcol)
    if np.isfinite(bb):
        nb.append(bb); npv.append(pp2)
nb, npv = np.array(nb), np.array(npv)
emp = float((np.abs(nb) >= abs(b)).mean())
print(f"\nmatched null ({len(nb)} sets): empirical p = {emp:.4f}, "
      f"random sets at p<0.05 = {(npv < 0.05).mean():.3f}")

# ------------------------------------------------------------------ label permutation
rng2 = np.random.default_rng(1)
pb = []
for _ in range(N_PERM):
    d3 = D.copy()
    d3["shed"] = rng2.permutation(d3["shed"].values)
    bb, _, _, _, _, _ = logit_or(d3, PRIMARY, pcol)
    if np.isfinite(bb):
        pb.append(bb)
pb = np.array(pb)
perm_p = float((np.abs(pb) >= abs(b)).mean())
print(f"label permutation ({len(pb)}): fraction |logOR| >= observed = {perm_p:.4f}")

# ------------------------------------------------------------------ secondary
sec = []
for name in scores:
    if name.startswith("proxy_") or name in EXCLUDE or name == PRIMARY:
        continue
    pc = f"proxy_{PROGRAMS[name]['cell_type']}"
    bb, ss, pp2, pp1, cc, nn = logit_or(D, name, pc)
    sec.append(dict(program=name, OR=np.exp(bb), p_two_sided=pp2, n=nn))
SEC = pd.DataFrame(sec)
if len(SEC):
    pv = SEC.p_two_sided.values
    o = np.argsort(pv); q = np.empty_like(pv)
    q[o] = np.minimum.accumulate((pv[o] * len(pv) / np.arange(1, len(pv) + 1))[::-1])[::-1]
    SEC["BH"] = np.minimum(q, 1.0)
print("\nSECONDARY (two-sided, BH over the other programs):")
print(SEC.round(4).to_string(index=False))

# ------------------------------------------------------------------ exploratory
expl = {}
sub = D[D.mmclass.isin(["MDS", "NS"])].copy()
sub["shed"] = (sub.mmclass == "MDS").astype(int)
bb, _, pp2, pp1, cc, nn = logit_or(sub, PRIMARY, pcol)
expl["MDS_vs_NS"] = dict(OR=np.exp(bb), p_one_sided=pp1, p_two_sided=pp2, n=nn)
di = D.copy()
di["__ix"] = di[PRIMARY] * di["vax"]
bb2, _, pp22, _, _, nn2 = logit_or(di, "__ix", pcol, extra=(PRIMARY,))
expl["P5_x_vaccination"] = dict(OR=np.exp(bb2), p_two_sided=pp22, n=nn2)
print("\nEXPLORATORY:", {k: {kk: round(vv, 4) if isinstance(vv, float) else vv
                            for kk, vv in v.items()} for k, v in expl.items()})

# ------------------------------------------------------------------ three-cohort meta
prior = pd.read_csv(OUT.parent / "mvri_outcome" / "per_cohort.csv")
prior = prior[prior.program == PRIMARY]
betas = list(prior.logOR.values) + [b]
ses = list(prior.se.values) + [se]
w = 1 / np.array(ses) ** 2
pooled = float((w * np.array(betas)).sum() / w.sum())
pse = float(np.sqrt(1 / w.sum()))
Q = float((w * (np.array(betas) - pooled) ** 2).sum())
I2 = max(0.0, (Q - (len(betas) - 1)) / Q * 100) if Q > 0 else 0.0
pmeta = float(2 * stats.norm.sf(abs(pooled / pse)))
print(f"\nTHREE-COHORT META: OR {np.exp(pooled):.3f} "
      f"(95% CI {np.exp(pooled - 1.96 * pse):.3f}-{np.exp(pooled + 1.96 * pse):.3f}), "
      f"p = {pmeta:.5f}, Q = {Q:.2f}, I2 = {I2:.1f}%")

pd.DataFrame([dict(program=PRIMARY, OR=np.exp(b), ci_lo=ci[0], ci_hi=ci[1], n=n,
                   p_one_sided=p1, p_two_sided=p2, OR_unadjusted=np.exp(bu),
                   p_one_sided_unadj=p1u, empirical_p=emp,
                   null_frac_p05=float((npv < 0.05).mean()), perm_p=perm_p,
                   meta_OR=np.exp(pooled), meta_p=pmeta, meta_I2=I2)]
             ).to_csv(OUT / "primary.csv", index=False)
SEC.to_csv(OUT / "secondary.csv", index=False)
pd.DataFrame(expl).T.to_csv(OUT / "exploratory.csv")

# ---- source data for the figure panels (added 2026-11-24, deviation 6: output only).
# These persist draws the run already made; no model is recomputed and the RNG seeds are
# fixed, so the saved values are the ones reported above.
pd.DataFrame(dict(logOR=nb, p_two_sided=npv)).to_csv(OUT / "matched_null_draws.csv",
                                                     index=False)
pd.DataFrame(dict(logOR=pb)).to_csv(OUT / "permutation_draws.csv", index=False)

_lab = list(prior.cohort.values) if "cohort" in prior.columns else ["RSV", "pneumococcal"]
forest = [dict(cohort=c, logOR=bb, se=ss, OR=np.exp(bb),
               ci_lo=np.exp(bb - 1.96 * ss), ci_hi=np.exp(bb + 1.96 * ss))
          for c, bb, ss in zip(_lab, prior.logOR.values, prior.se.values)]
forest.append(dict(cohort="NOVA influenza", logOR=b, se=se, OR=np.exp(b),
                   ci_lo=ci[0], ci_hi=ci[1]))
forest.append(dict(cohort="pooled (fixed effect)", logOR=pooled, se=pse,
                   OR=np.exp(pooled), ci_lo=np.exp(pooled - 1.96 * pse),
                   ci_hi=np.exp(pooled + 1.96 * pse)))
pd.DataFrame(forest).to_csv(OUT / "forest_three_cohort.csv", index=False)

# ------------------------------------------------------------------ verdict (plan 2.2)
print("\n" + "=" * 78)
print("VERDICT per plan section 2.2")
print("=" * 78)
if not np.isfinite(b):
    v = "test could not be fitted"
elif np.exp(b) > 1:
    v = "CONTRADICTED: point estimate in the opposite direction"
elif p1 < 0.05:
    v = "REPLICATED: one-sided p < 0.05 with OR < 1"
elif ci[1] < 0.6:
    v = "NOT REPLICATED at the plausible effect size (CI excludes 0.6)"
else:
    v = "UNDERPOWERED, DIRECTION CONSISTENT (OR < 1, p >= 0.05, CI includes 0.5-0.6)"
print(" ", v)
print("\nwrote", OUT)
