"""Do the curated epithelial programs predict a challenge outcome?

Implements docs/PREREG_mvri.md sections 3-6, locked at commit fecab57; the program family is
the seven curated programs of docs/MVRI_PROGRAMS.md, locked at d391a86. Both files were
committed before any outcome value from either cohort was read.

Primary: fixed-effect inverse-variance meta-analysis of the per-SD log odds ratio across the
RSV and pneumococcal challenge cohorts, BH across the seven programs, two-sided.

Usage: python validation/mvri_outcome.py <data_dir> <symbol_map.h5ad> <out_dir>
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
DATA, SYMMAP, OUT = Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3])
OUT.mkdir(parents=True, exist_ok=True)
N_SETS, N_BINS, MIN_COVER = 500, 25, 0.60
MIN_READS = 5e6


def load_rsv():
    s = pd.read_csv(DATA / "gse155237" / "g155_samples.csv")
    s = s[s.timepoint == "D0"].set_index("title")
    M = pd.read_csv(DATA / "gse155237" / "GSE155237_gene_counts.csv.gz", index_col=0)
    M = M.loc[:, [c for c in M.columns if c in s.index]]
    s = s.loc[M.columns]
    import anndata as ad
    vmap = ad.read_h5ad(SYMMAP, backed="r").var["feature_name"].astype(str)
    sym = vmap.reindex(M.index.astype(str).str.split(".").str[0])
    M.index = sym.values
    M = M[sym.notna().values].groupby(level=0).sum()
    y = (s.outcome == "Infected").astype(int)
    return M, pd.DataFrame({"y": y, "batch": s.batch}), "RSV infection"


def load_pneumo():
    recs, cur = [], None
    for line in open(DATA / "gse261277" / "GSE261277_gsm.txt", encoding="utf-8",
                     errors="replace"):
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
    meta = meta.rename(columns={"carriage status": "carriage"}).set_index("barcode")
    with gzip.open(DATA / "gse261277" / "GSE261277_raw_expression_elderly.tsv.gz", "rt") as fh:
        M = pd.read_csv(fh, sep="\t", index_col=0)
    M = M.groupby(level=0).sum()
    M = M.loc[:, [c for c in M.columns if c in meta.index]]
    meta = meta.loc[M.columns]
    sel = (meta.timepoint == "baseline") & (M.sum(0) >= MIN_READS)
    M, meta = M.loc[:, sel.values], meta[sel.values]
    y = (meta.carriage == "POS").astype(int)
    return M, pd.DataFrame({"y": y.values, "batch": "one"}, index=M.columns), \
        "pneumococcal carriage"


def score_cohort(M, info, label):
    depth = M.sum(0)
    cpm = M / depth * 1e6
    logc = np.log2(cpm + 1)
    expressed = (cpm.mean(1) >= 1) & (logc.std(1, ddof=1) > 0)
    Z = logc[expressed].sub(logc[expressed].mean(1), axis=0).div(
        logc[expressed].std(1, ddof=1), axis=0)
    cover, scores, used = [], {}, {}
    for name, p in PROGRAMS.items():
        g = [x for x in p["genes"] if x in Z.index]
        frac = len(g) / len(p["genes"])
        cover.append(dict(cohort=label, program=name, n_used=len(g),
                          n_genes=len(p["genes"]), frac=frac, transferable=frac >= MIN_COVER))
        if frac >= MIN_COVER:
            s = Z.loc[g].mean(0)
            scores[name] = (s - s.mean()) / s.std(ddof=1)   # OR is per SD
            used[name] = g
    for ct, markers in PROXY.items():
        g = [x for x in markers if x in Z.index]
        if g:
            s = Z.loc[g].mean(0)
            scores[f"proxy_{ct}"] = (s - s.mean()) / s.std(ddof=1)
    D = info.copy()
    for k, v in scores.items():
        D[k] = v.reindex(D.index).values
    D["log10_depth"] = np.log10(depth.reindex(D.index).values)
    return D, Z, logc[expressed], used, pd.DataFrame(cover)


def logistic_or(D, col, proxy_col):
    """Per-SD log OR for `col`, adjusted for the composition proxy and library size."""
    cols = [col, "log10_depth"] + ([proxy_col] if proxy_col else [])
    d = D.dropna(subset=cols + ["y"])
    X = sm.add_constant(d[cols].astype(float).values, has_constant="add")
    try:
        r = sm.Logit(d.y.values.astype(float), X).fit(disp=0)
        if not np.isfinite(r.bse[1]) or r.bse[1] > 10:
            return np.nan, np.nan, np.nan, len(d)
        return float(r.params[1]), float(r.bse[1]), float(r.pvalues[1]), len(d)
    except Exception:
        return np.nan, np.nan, np.nan, len(d)


def meta_fixed(betas, ses):
    b, s = np.asarray(betas, float), np.asarray(ses, float)
    ok = np.isfinite(b) & np.isfinite(s) & (s > 0)
    if ok.sum() == 0:
        return np.nan, np.nan, np.nan, np.nan, np.nan
    w = 1 / s[ok] ** 2
    pooled = float((w * b[ok]).sum() / w.sum())
    se = float(np.sqrt(1 / w.sum()))
    p = float(2 * stats.norm.sf(abs(pooled / se)))
    if ok.sum() < 2:
        return pooled, se, p, np.nan, np.nan
    Q = float((w * (b[ok] - pooled) ** 2).sum())
    df = ok.sum() - 1
    I2 = float(max(0.0, (Q - df) / Q * 100)) if Q > 0 else 0.0
    return pooled, se, p, Q, I2


def bh(p):
    p = np.asarray(p, float)
    ok = np.isfinite(p)
    q = np.full(p.shape, np.nan)
    pv = p[ok]
    o = np.argsort(pv)
    qq = np.empty_like(pv)
    qq[o] = np.minimum.accumulate((pv[o] * len(pv) / np.arange(1, len(pv) + 1))[::-1])[::-1]
    q[ok] = np.minimum(qq, 1.0)
    return q


cohorts = {}
for loader in (load_rsv, load_pneumo):
    M, info, label = loader()
    D, Z, logexpr, used, cover = score_cohort(M, info, label)
    cohorts[label] = dict(D=D, Z=Z, logexpr=logexpr, used=used, cover=cover)
    print(f"\n=== {label}: {D.shape[0]} subjects, outcome "
          f"{int(D.y.sum())} positive / {int((1 - D.y).sum())} negative")
    print(cover[["program", "n_used", "n_genes", "frac", "transferable"]].to_string(index=False))

COV = pd.concat([c["cover"] for c in cohorts.values()])
COV.to_csv(OUT / "transferability.csv", index=False)

print("\n" + "=" * 78)
print("PRIMARY: per-SD log OR, composition-adjusted, pooled fixed-effect")
print("=" * 78)
rows, nullrows, percohort = [], [], []
for name in PROGRAMS:
    pcol = f"proxy_{PROGRAMS[name]['cell_type']}"
    bs, ss = [], []
    for label, c in cohorts.items():
        if name not in c["used"]:
            bs.append(np.nan); ss.append(np.nan); continue
        b, se, p, n = logistic_or(c["D"], name, pcol)
        bu, _, pu, _ = logistic_or(c["D"], name, None)
        percohort.append(dict(program=name, cohort=label, logOR=b, se=se, p=p, n=n,
                              OR=np.exp(b) if np.isfinite(b) else np.nan,
                              logOR_unadjusted=bu, p_unadjusted=pu))
        bs.append(b); ss.append(se)
    pooled, se, p, Q, I2 = meta_fixed(bs, ss)
    rows.append(dict(program=name, cell_type=PROGRAMS[name]["cell_type"],
                     pooled_logOR=pooled, pooled_OR=np.exp(pooled), se=se, p=p, Q=Q, I2=I2))

R = pd.DataFrame(rows)
R["BH"] = bh(R.p.values)
PC = pd.DataFrame(percohort)
PC.to_csv(OUT / "per_cohort.csv", index=False)

# ---------------------------------------------------------- matched null (plan 4.2)
print("\nmatched null: 500 expression-matched random sets per program per cohort")
for name in PROGRAMS:
    pcol = f"proxy_{PROGRAMS[name]['cell_type']}"
    obs = R.loc[R.program == name, "pooled_logOR"].iloc[0]
    if not np.isfinite(obs):
        continue
    pools, binsmap, rngs = {}, {}, {}
    for label, c in cohorts.items():
        m = c["logexpr"].mean(1)
        b = pd.qcut(m.rank(method="first"), N_BINS, labels=False)
        binsmap[label] = b
        pools[label] = {k: np.array(v) for k, v in b.groupby(b).groups.items()}
        rngs[label] = np.random.default_rng(abs(hash((label, name))) % (2 ** 32))
    nb, npv = [], []
    for _ in range(N_SETS):
        bs, ss = [], []
        for label, c in cohorts.items():
            if name not in c["used"]:
                bs.append(np.nan); ss.append(np.nan); continue
            gb = [binsmap[label].get(x) for x in c["used"][name]
                  if x in binsmap[label].index]
            pick = [str(rngs[label].choice(pools[label][k])) for k in gb if pd.notna(k)]
            s = c["Z"].loc[pick].mean(0)
            d = c["D"].copy()
            d["__null"] = ((s - s.mean()) / s.std(ddof=1)).reindex(d.index).values
            b2, se2, _, _ = logistic_or(d, "__null", pcol)
            bs.append(b2); ss.append(se2)
        pb, _, pp, _, _ = meta_fixed(bs, ss)
        if np.isfinite(pb):
            nb.append(pb); npv.append(pp)
    nb, npv = np.array(nb), np.array(npv)
    nullrows.append(dict(program=name, pooled_logOR=obs, n_null=len(nb),
                         empirical_p=float((np.abs(nb) >= abs(obs)).mean()),
                         null_frac_p05=float((npv < 0.05).mean())))
    print(f"  {name:32s} pooled logOR {obs:+.3f}  emp_p "
          f"{nullrows[-1]['empirical_p']:.3f}  nullp05 {nullrows[-1]['null_frac_p05']:.3f}")

N = pd.DataFrame(nullrows)
R = R.merge(N[["program", "empirical_p", "null_frac_p05"]], on="program", how="left")
R.to_csv(OUT / "primary_meta.csv", index=False)
N.to_csv(OUT / "matched_null.csv", index=False)

pd.set_option("display.width", 220)
print("\n" + "=" * 78)
print("POOLED RESULT (two-sided, BH over the program family)")
print("=" * 78)
print(R.sort_values("p")[["program", "cell_type", "pooled_OR", "p", "BH", "I2",
                          "empirical_p", "null_frac_p05"]].round(4).to_string(index=False))
print("\nper cohort:")
print(PC[["program", "cohort", "OR", "p", "n", "logOR_unadjusted", "p_unadjusted"]]
      .round(4).to_string(index=False))

# ---------------------------------------------------------- plan section 6
print("\n" + "=" * 78)
print("PLAN SECTION 6: may the word 'index' be used?")
print("=" * 78)
ok = R[(R.BH < 0.05) & (R.empirical_p < 0.05) & (R.I2 <= 75)]
if len(ok) == 0:
    print("NO. No program meets condition 1 (BH < 0.05). The word 'index' is not used.")
else:
    print(ok.to_string(index=False))
print("\nwrote", OUT)
