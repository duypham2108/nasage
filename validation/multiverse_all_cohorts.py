"""Specification curve across all three challenge cohorts.

Closes deviation 2 of results/multiverse/PREREG_RESULTS.md: the locked plan
(docs/PREREG_multiverse.md §1) specifies RSV, pneumococcal AND NOVA; only NOVA had been run
because the other two per-gene matrices were not retained. They are now dumped by
mvri_outcome.py under NASAGE_DUMP_EXPR=1.

This tests whether the two cohorts that carry the surviving claim are specification-stable.

Usage: python validation/multiverse_all_cohorts.py <results_dir> <out_dir>
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

sys.path.insert(0, str(Path(__file__).parent))
from mvri_programs import PROGRAMS                 # noqa: E402

R, OUT = Path(sys.argv[1]), Path(sys.argv[2])
OUT.mkdir(parents=True, exist_ok=True)

VARIANTS = {"original": PROGRAMS}
try:
    from mvri_programs_epi import PROGRAMS_EPI
    VARIANTS["epithelial"] = PROGRAMS_EPI
except Exception:
    pass
try:
    sel = pd.read_csv(R / "multiverse" / "orthogonal_selection.csv")
    orth = {}
    for _, r in sel.iterrows():
        genes = [g for g in str(r.kept).split(",") if g]
        if len(genes) >= 8:
            orth[r.program + "_orth"] = dict(cell_type=PROGRAMS[r.program]["cell_type"],
                                             genes=genes)
    if orth:
        VARIANTS["orthogonal"] = orth
except Exception:
    pass


def sc_meanz(X):
    Z = (X - X.mean(0)) / np.where(X.std(0, ddof=1) > 0, X.std(0, ddof=1), 1)
    return Z.mean(1)


def sc_medianz(X):
    Z = (X - X.mean(0)) / np.where(X.std(0, ddof=1) > 0, X.std(0, ddof=1), 1)
    return np.median(Z, axis=1)


def sc_rank(X):
    order = np.argsort(np.argsort(-X, axis=1), axis=1) + 1
    return (1.0 / order).mean(1)


def sc_pc1(X):
    Z = (X - X.mean(0)) / np.where(X.std(0, ddof=1) > 0, X.std(0, ddof=1), 1)
    u, s, _ = np.linalg.svd(Z - Z.mean(0), full_matrices=False)
    v = u[:, 0] * s[0]
    r = np.corrcoef(v, Z.mean(1))[0, 1]
    return v * (np.sign(r) if np.isfinite(r) and r != 0 else 1.0)


SCORERS = {"mean_z": sc_meanz, "median_z": sc_medianz, "rank": sc_rank, "pc1_of_set": sc_pc1}
COVSETS = ["none", "proxy", "proxy+demog", "proxy+demog+PC1"]

COHORTS = {
    "RSV infection": ("expr_rsv.csv.gz", "meta_rsv.csv"),
    "pneumococcal carriage": ("expr_pneumo.csv.gz", "meta_pneumo.csv"),
}

rows = []
for label, (ef, mf) in COHORTS.items():
    ep, mp = R / "mvri_outcome" / ef, R / "mvri_outcome" / mf
    if not ep.exists():
        print(f"{label}: {ep} missing, skipped")
        continue
    E = pd.read_csv(ep, index_col=0)
    meta = pd.read_csv(mp, index_col=0)
    common = [c for c in E.columns if c in meta.index]
    E, meta = E[common], meta.loc[common]
    print(f"\n{label}: {E.shape[0]} genes x {len(common)} subjects, "
          f"{int(meta.y.sum())} positive")
    gpc1 = None
    for vname, VP in VARIANTS.items():
        for pname, p in VP.items():
            genes = [g for g in p["genes"] if g in E.index]
            if len(genes) < 3:
                continue
            Xg = E.loc[genes].T.values
            for sname, fn_ in SCORERS.items():
                v = fn_(Xg)
                sd = v.std(ddof=1)
                v = (v - v.mean()) / (sd if sd > 0 else 1)
                for cname in COVSETS:
                    X = pd.DataFrame({"s": v}, index=meta.index)
                    if "proxy" in cname:
                        pc = f"proxy_{p['cell_type']}"
                        if pc in meta.columns:
                            X[pc] = pd.to_numeric(meta[pc], errors="coerce").values
                    if "demog" in cname and "log10_depth" in meta.columns:
                        X["depth"] = pd.to_numeric(meta["log10_depth"],
                                                   errors="coerce").values
                    if "PC1" in cname:
                        if gpc1 is None:
                            Z = E.T.values
                            Z = (Z - Z.mean(0)) / np.where(Z.std(0, ddof=1) > 0,
                                                           Z.std(0, ddof=1), 1)
                            uu, ss, _ = np.linalg.svd(Z - Z.mean(0), full_matrices=False)
                            gpc1 = uu[:, 0] * ss[0]
                        X["gpc1"] = gpc1
                    try:
                        Xf = sm.add_constant(X.astype(float)).dropna()
                        res = sm.Logit(meta.y.astype(int).loc[Xf.index], Xf).fit(disp=0)
                        rows.append(dict(cohort=label, variant=vname, program=pname,
                                         scoring=sname, covars=cname,
                                         est=float(res.params["s"]),
                                         p=float(res.pvalues["s"])))
                    except Exception:
                        pass

SC = pd.DataFrame(rows)
if not len(SC):
    sys.exit("no specifications ran")
SC.to_csv(OUT / "specification_curve_rsv_pneumo.csv", index=False)
print(f"\n{len(SC)} specifications run")

print("\n" + "=" * 78)
print("P5 mucin barrier - the programme that carries the surviving claim")
print("=" * 78)
for coh, g in SC[SC.program == "P5_mucin_barrier"].groupby("cohort"):
    med = g.est.median()
    print(f"\n{coh}: {len(g)} specs | median OR {np.exp(med):.3f} | "
          f"{100 * (g.p < 0.05).mean():.0f}% p<0.05 | "
          f"{100 * (np.sign(g.est) == np.sign(med)).mean():.0f}% same sign")
    t = g.pivot_table(index="scoring", columns="covars", values="est").apply(np.exp)
    print(t.round(3).to_string())

print("\n" + "=" * 78)
print("ALL PROGRAMMES, per cohort (locked criterion: >=90% sign AND median p<0.05)")
print("=" * 78)
summ = []
for (coh, pr), g in SC.groupby(["cohort", "program"]):
    med = float(g.est.median())
    sign = float((np.sign(g.est) == np.sign(med)).mean())
    # median specification's p, taken as the p of the spec nearest the median estimate
    pmed = float(g.iloc[(g.est - med).abs().argsort().iloc[0]].p)
    summ.append(dict(cohort=coh, program=pr, n=len(g), median_OR=float(np.exp(med)),
                     frac_p05=float((g.p < 0.05).mean()), sign_stability=sign,
                     p_at_median=pmed,
                     robust=bool(sign >= 0.90 and pmed < 0.05)))
S = pd.DataFrame(summ).sort_values(["cohort", "median_OR"])
S.to_csv(OUT / "specification_summary.csv", index=False)
print(S.round(3).to_string(index=False))
print(f"\nrobust by the locked criterion: {int(S.robust.sum())} of {len(S)}")
print("\nwrote", OUT)
