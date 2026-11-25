"""Do mucosal vaccines engage the curated epithelial programs?

Implements docs/PREREG_mvri_engagement.md, locked at commit d391a86 before any expression
value from either dataset was read. Programs come from validation/mvri_programs.py, whose
gene lists are fixed by docs/MVRI_PROGRAMS.md.

Estimand: the time x arm interaction. Both arms were vaccinated against influenza; only the
LAIV arm received it on the nasal mucosa, so the interaction isolates the mucosal route.

Carriage outcome is NOT read here. This script never touches the 'carriage status' field.

Usage: python validation/mvri_engagement.py <liverpool_dir> <symbol_map.h5ad> <out_dir>
"""
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

from mvri_programs import PROGRAMS, PROXY

warnings.filterwarnings("ignore")
DATA, SYMMAP, OUT = Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3])
OUT.mkdir(parents=True, exist_ok=True)
N_SETS, N_BINS, MIN_COVER = 500, 25, 0.60
rng_global = np.random.default_rng(0)

DATASETS = {
    "GSE117580": dict(counts="GSE117580_counts.csv.gz", sep=",", gsm="GSE117580_gsm.txt",
                      baseline="D-4", posts=["D2", "D9"]),
    "GSE164649": dict(counts="GSE164649_counts.txt.gz", sep="\t", gsm="GSE164649_gsm.txt",
                      baseline="baseline", posts=["D2", "D6", "D9"]),
}


def parse_gsm(path):
    recs, cur = [], None
    for line in open(path, encoding="utf-8", errors="replace"):
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
    m = pd.DataFrame(recs)
    # the expression column name is the S-number inside the title
    import re
    m["sample"] = m.title.str.extract(r"(S\d+A?)")[0]
    m["arm"] = m["vaccine"].str.contains("LAIV").map({True: "LAIV", False: "TIV"})
    m["subject"] = m.get("subject id", pd.Series(index=m.index, dtype=object))
    if m.subject.isna().all():
        m["subject"] = m.title.str.extract(r"(\d+/\d+/\d+)")[0]
    return m.set_index("sample")


import anndata as ad
vmap = ad.read_h5ad(SYMMAP, backed="r").var["feature_name"].astype(str)

summary_rows, cover_rows, null_rows = [], [], []

for ds, cfg in DATASETS.items():
    print("\n" + "=" * 78)
    print(ds)
    print("=" * 78)
    meta = parse_gsm(DATA / cfg["gsm"])
    CNT = pd.read_csv(DATA / cfg["counts"], sep=cfg["sep"], index_col=0)
    CNT.index = CNT.index.astype(str).str.split(".").str[0]
    shared = [c for c in CNT.columns if c in meta.index]
    CNT, meta = CNT[shared], meta.loc[shared]
    sym = vmap.reindex(CNT.index)
    CNT = CNT[sym.notna()].groupby(sym[sym.notna()].values).sum()
    print(f"{CNT.shape[1]} samples x {CNT.shape[0]} symbols; "
          f"arms {meta.arm.value_counts().to_dict()}; "
          f"timepoints {meta.timepoint.value_counts().to_dict()}")

    depth = CNT.sum(0)
    cpm = CNT / depth * 1e6
    logc = np.log2(cpm + 1)
    expressed = (cpm.mean(1) >= 1) & (logc.std(1, ddof=1) > 0)
    Z = logc[expressed].sub(logc[expressed].mean(1), axis=0).div(
        logc[expressed].std(1, ddof=1), axis=0)
    print(f"features at CPM >= 1: {int(expressed.sum())}")

    # ---- transferability, reported before any model (plan section 7.2)
    scores, used = {}, {}
    for name, p in PROGRAMS.items():
        g = [x for x in p["genes"] if x in Z.index]
        frac = len(g) / len(p["genes"])
        cover_rows.append(dict(dataset=ds, program=name, cell_type=p["cell_type"],
                               n_genes=len(p["genes"]), n_used=len(g), frac=frac,
                               transferable=frac >= MIN_COVER,
                               dropped=";".join(sorted(set(p["genes"]) - set(g)))))
        if frac >= MIN_COVER:
            scores[name] = Z.loc[g].mean(0)
            used[name] = g
    for ct, markers in PROXY.items():
        g = [x for x in markers if x in Z.index]
        scores[f"proxy_{ct}"] = Z.loc[g].mean(0) if g else np.nan
    S = pd.DataFrame(scores)
    cov = pd.DataFrame([r for r in cover_rows if r["dataset"] == ds])
    print("\ntransferability (before any model):")
    print(cov[["program", "n_used", "n_genes", "frac", "transferable"]].to_string(index=False))

    D = meta.join(S)
    D["log10_depth"] = np.log10(depth.loc[D.index])
    D["sex"] = D.get("sex", pd.Series("U", index=D.index)).fillna("U")
    D.to_csv(OUT / f"{ds}_scores.csv")

    def fit_delta_ols(d, proxy_col, cfg_baseline):
        """Pre-declared fallback, upgraded: OLS on per-subject paired deltas.

        The plan named a Mann-Whitney comparison of deltas as the fallback. That test
        carries no composition adjustment, which the plan requires of the primary, so the
        fallback is implemented as an OLS on the same deltas with the proxy delta as a
        covariate. Dated deviation; it is only reached if MixedLM fails.
        """
        base = d[d.post == 0].drop_duplicates("subject").set_index("subject")
        post = d[d.post == 1].drop_duplicates("subject").set_index("subject")
        common = base.index.intersection(post.index)
        if len(common) < 6:
            return np.nan, np.nan, len(common), "insufficient_pairs"
        T = pd.DataFrame({
            "dy": post.loc[common, "y"] - base.loc[common, "y"],
            "ddepth": post.loc[common, "log10_depth"] - base.loc[common, "log10_depth"],
            "laiv": (post.loc[common, "arm"] == "LAIV").astype(int),
            "sex": post.loc[common, "sex"].values,
        })
        f = "dy ~ laiv + ddepth"
        if proxy_col is not None:
            T["dproxy"] = (post.loc[common, proxy_col] - base.loc[common, proxy_col]).values
            f += " + dproxy"
        if T.sex.nunique() > 1:
            f += " + C(sex)"
        T = T.dropna()
        if len(T) < 6:
            return np.nan, np.nan, len(T), "insufficient_pairs"
        m = smf.ols(f, data=T).fit()
        return m.params["laiv"], m.pvalues["laiv"], len(T), "delta_ols"

    def fit_interaction(df, col, proxy_col):
        """Time x arm interaction: MixedLM with a subject intercept, as the plan specifies.

        The estimand is post:laiv — the difference between arms in change from baseline.
        """
        d = df.dropna(subset=[col, "log10_depth", "subject"]).copy()
        d["y"] = d[col]
        d["post"] = (d.timepoint != cfg["baseline"]).astype(int)
        d["laiv"] = (d.arm == "LAIV").astype(int)
        if proxy_col is not None:
            d = d.dropna(subset=[proxy_col])
        f = "y ~ post * laiv + log10_depth"
        if proxy_col is not None:
            f += f" + {proxy_col}"
        if d.sex.nunique() > 1:
            f += " + C(sex)"
        try:
            m = smf.mixedlm(f, d, groups=d.subject).fit(reml=False)
            b, p = m.params.get("post:laiv", np.nan), m.pvalues.get("post:laiv", np.nan)
            if np.isfinite(b) and np.isfinite(p):
                return float(b), float(p), len(d), "mixedlm"
        except Exception:
            pass
        return fit_delta_ols(d, proxy_col, cfg["baseline"])

    expr_mean = logc[expressed].mean(1)
    bins = pd.qcut(expr_mean.rank(method="first"), N_BINS, labels=False)
    pool = {b: np.array(idx) for b, idx in bins.groupby(bins).groups.items()}

    for tp in cfg["posts"]:
        sub = D[D.timepoint.isin([cfg["baseline"], tp])]
        for name in used:
            pcol = f"proxy_{PROGRAMS[name]['cell_type']}"
            beta, p, n, how = fit_interaction(sub, name, pcol)
            summary_rows.append(dict(dataset=ds, timepoint=tp, program=name,
                                     cell_type=PROGRAMS[name]["cell_type"],
                                     beta_interaction=beta, p=p, n=n, method=how,
                                     adjusted=True))
            # unadjusted sensitivity: drop the proxy term, do not zero it
            b2, p2, _, _ = fit_interaction(sub, name, None)
            summary_rows.append(dict(dataset=ds, timepoint=tp, program=name,
                                     cell_type=PROGRAMS[name]["cell_type"],
                                     beta_interaction=b2, p=p2, n=n, method=how,
                                     adjusted=False))
            # matched null
            if np.isnan(beta):
                continue
            gb = [bins.get(x) for x in used[name] if x in bins.index]
            rng = np.random.default_rng(abs(hash((ds, tp, name))) % (2**32))
            nb, npv = [], []
            for _ in range(N_SETS):
                pick = [str(rng.choice(pool[b])) for b in gb if pd.notna(b)]
                s3 = sub.copy()
                s3["__null"] = Z.loc[pick].mean(0).reindex(s3.index)
                bb, pp, _, _ = fit_interaction(s3, "__null", pcol)
                if not np.isnan(bb):
                    nb.append(bb); npv.append(pp)
            nb, npv = np.array(nb), np.array(npv)
            null_rows.append(dict(dataset=ds, timepoint=tp, program=name,
                                  beta=beta, p=p, n_null=len(nb),
                                  empirical_p=float((np.abs(nb) >= abs(beta)).mean()),
                                  null_frac_p05=float((npv < 0.05).mean())))
            print(f"  {tp} {name:34s} beta {beta:+.4f} p {p:.4f} "
                  f"emp_p {null_rows[-1]['empirical_p']:.3f} "
                  f"nullp05 {null_rows[-1]['null_frac_p05']:.3f}")


def bh(p):
    p = np.asarray(p, float)
    ok = ~np.isnan(p)
    q = np.full(p.shape, np.nan)
    pv = p[ok]
    o = np.argsort(pv)
    qq = np.empty_like(pv)
    qq[o] = np.minimum.accumulate((pv[o] * len(pv) / np.arange(1, len(pv) + 1))[::-1])[::-1]
    q[ok] = np.minimum(qq, 1.0)
    return q


R = pd.DataFrame(summary_rows)
R["BH"] = np.nan
for ds in DATASETS:
    m = (R.dataset == ds) & R.adjusted
    R.loc[m, "BH"] = bh(R.loc[m, "p"].values)   # family = programs x timepoints, per plan
pd.DataFrame(cover_rows).to_csv(OUT / "transferability.csv", index=False)
R.to_csv(OUT / "engagement_interaction.csv", index=False)
pd.DataFrame(null_rows).to_csv(OUT / "engagement_matched_null.csv", index=False)

print("\n" + "=" * 78)
print("PRIMARY: time x arm interaction, composition-adjusted, BH within dataset")
print("=" * 78)
pd.set_option("display.width", 200)
print(R[R.adjusted].sort_values(["dataset", "p"])
      [["dataset", "timepoint", "program", "beta_interaction", "p", "BH", "n", "method"]]
      .round(4).to_string(index=False))
print("\nwrote", OUT)
