"""Specification curve on the lifespan trajectories, plus scoring-method concordance.

Closes deviation 3 of results/multiverse/PREREG_RESULTS.md (the lifespan arm of the locked grid,
docs/PREREG_multiverse.md §1) and quantifies WEAKNESSES.md D4 -- bulk scoring uses mean-z while
single cell uses a rank statistic, so results from the two paths have never been commensurable.

Usage: python validation/multiverse_lifespan.py <data_dir> <results_dir> <out_dir>
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp
import statsmodels.formula.api as smf

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from nasage.modules import load_modules      # noqa: E402
from mvri_programs import PROGRAMS           # noqa: E402

DATA, R, OUT = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
OUT.mkdir(parents=True, exist_ok=True)
MOD = load_modules()
AXES = {
    "Mucus barrier (P5)": PROGRAMS["P5_mucin_barrier"]["genes"],
    "Antigen presentation (M06)": MOD["M06_mhc_I_processing"]["genes"],
    "Innate sensing (M04)": MOD["M04_innate_sensing"]["genes"],
    "IgA support (P4)": PROGRAMS["P4_iga_transport"]["genes"],
    "Cytokine tone (M12)": MOD["M12_inflammatory_sasp"]["genes"],
}
FAMS = ["Basal", "Secretory", "Ciliated"]
COVSETS = ["none", "depth", "depth+sex", "depth+sex+PC1"]

BASAL_K = ("basal",)
CIL_K = ("ciliated", "cil_diff", "cil_acute", "multiciliated")
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


A = ad.read_h5ad(DATA / "loske2022_controls.h5ad")
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
ngl = np.asarray((A.X > 0).sum(1)).ravel()
OBS = pd.DataFrame({"__donor": A.obs.__donor.values, "__fam": A.obs.__fam.values,
                    "age": A.obs.__age.values,
                    "sex": A.obs["gender"].astype(str).values,
                    "log_genes": np.log(ngl + 1)})
print(f"{A.n_obs} cells, {OBS.__donor.nunique()} donors, age "
      f"{OBS.age.min():.0f}-{OBS.age.max():.0f}")


def getX(genes):
    present = [g for g in genes if g in A.var_names]
    if len(present) < 3:
        return None
    Z = A[:, present].X
    return np.asarray(Z.todense()) if sp.issparse(Z) else np.asarray(Z)


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

# donor x cell-type scores for every (axis, scoring) -- also the input to the D4 concordance
DON = {}
for nm, genes in AXES.items():
    Xg = getX(genes)
    if Xg is None:
        continue
    for sname, fn_ in SCORERS.items():
        v = fn_(Xg)
        sd = v.std(ddof=1)
        t = OBS.copy()
        t["y"] = (v - v.mean()) / (sd if sd > 0 else 1)
        g = (t.groupby(["__donor", "__fam"])
             .agg({"y": "mean", "age": "mean", "log_genes": "mean", "sex": "first"})
             .join(t.groupby(["__donor", "__fam"]).size().rename("n_cells"))
             .reset_index())
        DON[(nm, sname)] = g[g.n_cells >= 20]

# global PC1 per cell type, from the mean_z score matrix (same covariate as before)
PC1 = {}
for f in FAMS:
    cols = {}
    for nm in AXES:
        k = (nm, "mean_z")
        if k in DON:
            d = DON[k]
            cols[nm] = d[d.__fam == f].set_index("__donor").y
    if len(cols) >= 3:
        Mx = pd.DataFrame(cols).dropna()
        if len(Mx) >= 10:
            Z = (Mx - Mx.mean()) / Mx.std(ddof=1)
            u, s, _ = np.linalg.svd(Z.values - Z.values.mean(0), full_matrices=False)
            PC1[f] = pd.Series(u[:, 0] * s[0], index=Mx.index)

rows = []
for (nm, sname), g in DON.items():
    for f in FAMS:
        d = g[g.__fam == f].copy()
        if len(d) < 8:
            continue
        if f in PC1:
            d["pc1"] = PC1[f].reindex(d.__donor).values
        for cname in COVSETS:
            terms = ["age"]
            if "depth" in cname:
                terms.append("log_genes")
            if "sex" in cname:
                terms.append("sex")
            if "PC1" in cname:
                if "pc1" not in d.columns or d.pc1.isna().all():
                    continue
                terms.append("pc1")
            dd = d.dropna(subset=["y", "age"] + [t for t in terms if t in d.columns])
            if len(dd) < 8:
                continue
            try:
                m = smf.ols("y ~ " + " + ".join(terms), data=dd).fit()
                rows.append(dict(axis=nm, family=f, scoring=sname, covars=cname,
                                 n=len(dd), est=float(m.params["age"]),
                                 p=float(m.pvalues["age"])))
            except Exception:
                pass

SC = pd.DataFrame(rows)
SC.to_csv(OUT / "specification_curve_lifespan.csv", index=False)
print(f"\n{len(SC)} specifications run")

print("\n" + "=" * 78)
print("LIFESPAN SPECIFICATION CURVE (locked criterion: >=90% sign AND median p<0.05)")
print("=" * 78)
summ = []
for (nm, f), g in SC.groupby(["axis", "family"]):
    med = float(g.est.median())
    sign = float((np.sign(g.est) == np.sign(med)).mean())
    pmed = float(g.iloc[(g.est - med).abs().argsort().iloc[0]].p)
    summ.append(dict(axis=nm, family=f, n=len(g), median_beta=med,
                     frac_p05=float((g.p < 0.05).mean()), sign_stability=sign,
                     p_at_median=pmed, robust=bool(sign >= 0.90 and pmed < 0.05)))
S = pd.DataFrame(summ).sort_values("median_beta")
S.to_csv(OUT / "specification_summary_lifespan.csv", index=False)
print(S.round(4).to_string(index=False))
print(f"\nrobust: {int(S.robust.sum())} of {len(S)}")

# ---------------------------------------------------------------- D4: scoring concordance
print("\n" + "=" * 78)
print("D4: do mean-z and rank scoring agree? (bulk path vs single-cell path)")
print("=" * 78)
conc = []
for nm in AXES:
    a, b = DON.get((nm, "mean_z")), DON.get((nm, "rank"))
    if a is None or b is None:
        continue
    j = a.merge(b, on=["__donor", "__fam"], suffixes=("_meanz", "_rank"))
    for f in FAMS:
        s = j[j.__fam == f]
        if len(s) < 8:
            continue
        conc.append(dict(axis=nm, family=f, n=len(s),
                         pearson=float(np.corrcoef(s.y_meanz, s.y_rank)[0, 1]),
                         spearman=float(pd.Series(s.y_meanz).corr(
                             pd.Series(s.y_rank).reset_index(drop=True), method="spearman"))))
CD = pd.DataFrame(conc)
CD.to_csv(OUT / "scoring_concordance.csv", index=False)
print(CD.round(3).to_string(index=False))
if len(CD):
    print(f"\nmean Pearson between the two scoring paths: {CD.pearson.mean():.3f}")
print("\nwrote", OUT)
