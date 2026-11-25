"""Re-run all three specification curves with the CORRECTED rank scorer.

See validation/scoring_methods.py for the bug: the old rank statistic ranked a programme's
genes against each other rather than against the transcriptome. Every curve that used it is
re-run here, and the old and new results are compared directly.

Usage: python validation/multiverse_rerun.py <data_dir> <results_dir> <out_dir>
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp
import statsmodels.api as sm
import statsmodels.formula.api as smf

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from nasage.modules import load_modules                       # noqa: E402
from mvri_programs import PROGRAMS                            # noqa: E402
from scoring_methods import SCORERS                           # noqa: E402

DATA, R, OUT = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
OUT.mkdir(parents=True, exist_ok=True)
MOD = load_modules()
FAMS = ["Basal", "Secretory", "Ciliated"]
rows = []


def run_grid(E, outcome, label, programs, covbuilder, covsets, logistic=True):
    """E: samples x genes DataFrame. outcome: Series aligned to E.index."""
    full = E.values
    colpos = {g: i for i, g in enumerate(E.columns)}
    gpc1 = None
    for pname, p in programs.items():
        idx = [colpos[g] for g in p["genes"] if g in colpos]
        if len(idx) < 3:
            continue
        Xg = full[:, idx]
        for sname, fn_ in SCORERS.items():
            try:
                v = fn_(Xg, full=full, idx=idx)
            except Exception:
                continue
            sd = np.std(v, ddof=1)
            v = (v - np.mean(v)) / (sd if sd > 0 else 1)
            for cname in covsets:
                X = pd.DataFrame({"s": v}, index=E.index)
                X, gpc1 = covbuilder(X, cname, p, full, gpc1)
                try:
                    Xf = sm.add_constant(X.astype(float)).dropna()
                    y = outcome.loc[Xf.index]
                    if logistic:
                        res = sm.Logit(y.astype(int), Xf).fit(disp=0)
                    else:
                        res = sm.OLS(y.astype(float), Xf).fit()
                    rows.append(dict(target=label, program=pname, scoring=sname,
                                     covars=cname, est=float(res.params["s"]),
                                     p=float(res.pvalues["s"])))
                except Exception:
                    pass


# ======================================================== challenge cohorts
def cov_bulk(X, cname, p, full, gpc1):
    meta = cov_bulk.meta
    if "proxy" in cname:
        pc = f"proxy_{p['cell_type']}"
        if pc in meta.columns:
            X[pc] = pd.to_numeric(meta[pc], errors="coerce").values
    if "demog" in cname:
        for cc in ("Age", "male", "vax", "log10_depth"):
            if cc in meta.columns:
                X[cc] = pd.to_numeric(meta[cc], errors="coerce").values
    if "PC1" in cname:
        if gpc1 is None:
            Z = (full - full.mean(0)) / np.where(full.std(0, ddof=1) > 0,
                                                 full.std(0, ddof=1), 1)
            u, s, _ = np.linalg.svd(Z - Z.mean(0), full_matrices=False)
            gpc1 = u[:, 0] * s[0]
        X["gpc1"] = gpc1
    return X, gpc1


COVSETS_B = ["none", "proxy", "proxy+demog", "proxy+demog+PC1"]
SOURCES = [
    ("nova", "expr_symbols.csv.gz", "nova_scores_outcome.csv", "shed", "NOVA influenza"),
    ("mvri_outcome", "expr_rsv.csv.gz", "meta_rsv.csv", "y", "RSV infection"),
    ("mvri_outcome", "expr_pneumo.csv.gz", "meta_pneumo.csv", "y", "pneumococcal carriage"),
]
for d, ef, mf, ycol, label in SOURCES:
    ep, mp = R / d / ef, R / d / mf
    if not ep.exists() or not mp.exists():
        print(f"  skip {label}: missing {ep if not ep.exists() else mp}")
        continue
    E = pd.read_csv(ep, index_col=0).T            # samples x genes
    meta = pd.read_csv(mp, index_col=0)
    common = [i for i in E.index if i in meta.index]
    E, meta = E.loc[common], meta.loc[common]
    cov_bulk.meta = meta
    print(f"{label}: {E.shape[0]} samples x {E.shape[1]} genes")
    run_grid(E, meta[ycol], label, PROGRAMS, cov_bulk, COVSETS_B, logistic=True)

# ======================================================== lifespan (Loske pseudobulk)
print("\nlifespan: building donor x cell-type pseudobulk")
A = ad.read_h5ad(DATA / "loske2022_controls.h5ad")
BASAL_K, CIL_K = ("basal",), ("ciliated", "cil_diff", "cil_acute", "multiciliated")
SEC_K = ("secretory", "goblet", "club", "glandular", "sec_diff")
DROP_K = ("foxn4", "squamous", "ionocyte", "lowrna", "low_rna", "deutero")


def to_family(lbl):
    s = str(lbl).lower().replace(" ", "_")
    if any(k in s for k in DROP_K):
        return None
    if any(k in s for k in BASAL_K):
        return "Basal"
    if any(k in s for k in CIL_K):
        return "Ciliated"
    if any(k in s for k in SEC_K):
        return "Secretory"
    return None


A.obs["__fam"] = A.obs["detailed_clusters"].map(to_family)
A = A[A.obs.__fam.notna()].copy()
A.obs["__donor"] = A.obs["Patient_ID"].astype(str)
A.obs["__age"] = pd.to_numeric(A.obs["paper_age"], errors="coerce")
A = A[A.obs.__age.notna()].copy()
X = A.X
if float(X.max()) > 50:
    tot = np.asarray(X.sum(1)).ravel(); tot[tot == 0] = 1
    A.X = sp.diags(1e4 / tot) @ X
    A.X.data = np.log1p(A.X.data)
Xd = np.asarray(A.X.todense()) if sp.issparse(A.X) else np.asarray(A.X)
pb = pd.DataFrame(Xd, columns=A.var_names)
pb["__donor"] = A.obs.__donor.values
pb["__fam"] = A.obs.__fam.values
cnt = pb.groupby(["__donor", "__fam"]).size().rename("n_cells")
agg = pb.groupby(["__donor", "__fam"]).mean(numeric_only=True).join(cnt)
agg = agg[agg.n_cells >= 20]
meta_l = (pd.DataFrame({"__donor": A.obs.__donor.values, "age": A.obs.__age.values,
                        "sex": A.obs["gender"].astype(str).values,
                        "lg": np.log((Xd > 0).sum(1) + 1)})
          .groupby("__donor").agg({"age": "mean", "sex": "first", "lg": "mean"}))
del A, Xd, pb

AXES = {
    "Mucus barrier (P5)": PROGRAMS["P5_mucin_barrier"]["genes"],
    "Antigen presentation (M06)": MOD["M06_mhc_I_processing"]["genes"],
    "IgA support (P4)": PROGRAMS["P4_iga_transport"]["genes"],
}   # M04 and M12 fail the measurability gate and are NOT tested (see deviation 1)
GENECOLS = [c for c in agg.columns if c not in ("n_cells",)]
for f in FAMS:
    sub = agg.xs(f, level="__fam") if f in agg.index.get_level_values("__fam") else None
    if sub is None or len(sub) < 8:
        continue
    full = sub[GENECOLS].values
    colpos = {g: i for i, g in enumerate(GENECOLS)}
    mm = meta_l.reindex(sub.index)
    for nm, genes in AXES.items():
        idx = [colpos[g] for g in genes if g in colpos]
        if len(idx) < 3:
            continue
        Xg = full[:, idx]
        for sname, fn_ in SCORERS.items():
            try:
                v = fn_(Xg, full=full, idx=idx)
            except Exception:
                continue
            sd = np.std(v, ddof=1)
            d = pd.DataFrame({"y": (v - np.mean(v)) / (sd if sd > 0 else 1),
                              "age": mm.age.values, "sex": mm.sex.values,
                              "lg": mm.lg.values}).dropna()
            for cname in ["none", "depth", "depth+sex"]:
                terms = ["age"] + (["lg"] if "depth" in cname else []) + \
                        (["sex"] if "sex" in cname else [])
                try:
                    m = smf.ols("y ~ " + " + ".join(terms), data=d).fit()
                    rows.append(dict(target=f"lifespan|{f}", program=nm, scoring=sname,
                                     covars=cname, est=float(m.params["age"]),
                                     p=float(m.pvalues["age"])))
                except Exception:
                    pass

SC = pd.DataFrame(rows)
SC.to_csv(OUT / "specification_curve_corrected.csv", index=False)
print(f"\n{len(SC)} specifications run with the corrected scorer")

print("\n" + "=" * 78)
print("CORRECTED CURVES (>=90% sign AND median p<0.05 = robust)")
print("=" * 78)
summ = []
for (tg, pr), g in SC.groupby(["target", "program"]):
    med = float(g.est.median())
    sign = float((np.sign(g.est) == np.sign(med)).mean())
    pmed = float(g.iloc[(g.est - med).abs().argsort().iloc[0]].p)
    summ.append(dict(target=tg, program=pr, n=len(g), median_est=med,
                     frac_p05=float((g.p < 0.05).mean()), sign_stability=sign,
                     p_at_median=pmed, robust=bool(sign >= 0.90 and pmed < 0.05)))
S = pd.DataFrame(summ).sort_values(["target", "median_est"])
S.to_csv(OUT / "specification_summary_corrected.csv", index=False)
key = S[S.program.isin(["P5_mucin_barrier", "Mucus barrier (P5)",
                        "Antigen presentation (M06)"])]
print(key.round(4).to_string(index=False))
print(f"\nrobust overall: {int(S.robust.sum())} of {len(S)}")
print("\nwrote", OUT)
