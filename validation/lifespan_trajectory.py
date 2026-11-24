"""Lifespan trajectory of mucosal defence programmes, by cell type.

Plan: docs/PREREG_lifespan_trajectory.md, locked at 74cbff8.
Primary cohort: Loske 2022, ages 4-77, 41 donors, cell types annotated.

Usage: python validation/lifespan_trajectory.py <data_dir> <out_dir>
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
from nasage.modules import load_modules          # noqa: E402
from nasage.stats import effective_tests         # noqa: E402
from mvri_programs import PROGRAMS               # noqa: E402

DATA, OUT = Path(sys.argv[1]), Path(sys.argv[2])
OUT.mkdir(parents=True, exist_ok=True)
MIN_CELLS, N_SETS, N_BINS, N_PERM = 20, 500, 25, 2000
FAMS = ("Basal", "Secretory", "Ciliated")

M = load_modules()
# the five axes of the user's design, mapped to referenced sets (plan section 1)
AXES = {
    "Mucus barrier (P5)": PROGRAMS["P5_mucin_barrier"]["genes"],
    "Antigen presentation (M06)": M["M06_mhc_I_processing"]["genes"],
    "Innate sensing (M04)": M["M04_innate_sensing"]["genes"],
    "IgA support (P4)": PROGRAMS["P4_iga_transport"]["genes"],
    "Cytokine tone (M12)": M["M12_inflammatory_sasp"]["genes"],
}

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
A.obs["__sex"] = A.obs["gender"].astype(str)

X = A.X
if float(X.max()) > 50:
    tot = np.asarray(X.sum(1)).ravel(); tot[tot == 0] = 1
    A.X = sp.diags(1e4 / tot) @ X
    A.X.data = np.log1p(A.X.data)

print(f"{A.n_obs} epithelial cells, {A.obs.__donor.nunique()} donors, "
      f"age {A.obs.__age.min():.0f}-{A.obs.__age.max():.0f}")


def detection(sub, genes):
    present = [g for g in genes if g in A.var_names]
    if not present:
        return 0.0
    Z = sub[:, present].X
    frac = np.asarray((Z > 0).mean(0)).ravel()
    return float((frac >= 0.05).mean())


def score(genes):
    present = [g for g in genes if g in A.var_names]
    if not present:
        return None
    Z = A[:, present].X
    Z = np.asarray(Z.todense()) if sp.issparse(Z) else np.asarray(Z)
    m = Z.mean(1)
    sd = m.std(ddof=1)
    return (m - m.mean()) / sd if sd > 0 else m * 0.0


# ---------------------------------------------------------------- measurability FIRST
cov = []
for name, genes in AXES.items():
    for f in FAMS:
        sub = A[A.obs.__fam == f]
        if sub.n_obs:
            cov.append(dict(axis=name, family=f, n_cells=sub.n_obs,
                            frac_detected=detection(sub, genes)))
COV = pd.DataFrame(cov)
COV["passes_gate"] = COV.frac_detected >= 0.60
COV.to_csv(OUT / "measurability.csv", index=False)
print("\n=== MEASURABILITY GATE, before any model (plan 3) ===")
print(COV.pivot_table(index="axis", columns="family", values="frac_detected").round(2).to_string())
failed = sorted(COV.loc[~COV.passes_gate, "axis"].unique())
print("fails the 60% gate somewhere:", failed if failed else "none")

# ---------------------------------------------------------------- donor x cell-type frame
S = pd.DataFrame(index=A.obs_names)
for name, genes in AXES.items():
    v = score(genes)
    if v is not None:
        S[name] = v
ng = np.asarray((A.X > 0).sum(1)).ravel()
S["log_genes"] = np.log(ng + 1)
S["__donor"] = A.obs.__donor.values
S["__fam"] = A.obs.__fam.values
S["age"] = A.obs.__age.values
S["sex"] = A.obs.__sex.values
cnt = S.groupby(["__donor", "__fam"]).size().rename("n_cells")
D = (S.groupby(["__donor", "__fam"]).agg(
        {**{c: "mean" for c in list(AXES) + ["log_genes", "age"]}, "sex": "first"})
     .join(cnt).reset_index())
D = D[D.n_cells >= MIN_CELLS]
D.to_csv(OUT / "donor_celltype_scores.csv", index=False)
print(f"\n{len(D)} donor x cell-type units after the {MIN_CELLS}-cell floor")
print(D.groupby("__fam").size().to_dict())


def fit(d, col, quad=False):
    d = d.rename(columns={col: "y"})
    f = "y ~ age + sex + log_genes" + (" + I(age**2)" if quad else "")
    try:
        return smf.ols(f, data=d).fit()
    except Exception:
        return None


rows = []
for name in AXES:
    for f in FAMS:
        d = D[D.__fam == f].dropna(subset=[name, "age", "log_genes"])
        gate = COV[(COV.axis == name) & (COV.family == f)]
        if len(d) < 8 or (len(gate) and not bool(gate.passes_gate.iloc[0])):
            rows.append(dict(axis=name, family=f, n=len(d), tested=False,
                             reason="fails measurability gate" if len(gate) and
                             not bool(gate.passes_gate.iloc[0]) else "too few donors"))
            continue
        r1, r2 = fit(d, name), fit(d, name, quad=True)
        if r1 is None:
            continue
        rows.append(dict(
            axis=name, family=f, n=len(d), tested=True,
            beta_age=float(r1.params["age"]), p_age=float(r1.pvalues["age"]),
            beta_quad=float(r2.params["I(age ** 2)"]) if r2 is not None else np.nan,
            p_quad=float(r2.pvalues["I(age ** 2)"]) if r2 is not None else np.nan,
            dAIC=float(r1.aic - r2.aic) if r2 is not None else np.nan))
R = pd.DataFrame(rows)
tested = R[R.tested == True]                                      # noqa: E712
pv = tested.p_age.values
o = np.argsort(pv); q = np.empty_like(pv)
q[o] = np.minimum.accumulate((pv[o] * len(pv) / np.arange(1, len(pv) + 1))[::-1])[::-1]
R.loc[tested.index, "BH_age"] = np.minimum(q, 1.0)
pq = tested.p_quad.values
o = np.argsort(pq); q2 = np.empty_like(pq)
q2[o] = np.minimum.accumulate((pq[o] * len(pq) / np.arange(1, len(pq) + 1))[::-1])[::-1]
R.loc[tested.index, "BH_quad"] = np.minimum(q2, 1.0)

# ---------------------------------------------------------------- matched null
expr = pd.Series(np.asarray(A.X.mean(0)).ravel(), index=A.var_names)
bins = pd.qcut(expr.rank(method="first"), N_BINS, labels=False)
pool = {k: np.array(v) for k, v in bins.groupby(bins).groups.items()}
rng = np.random.default_rng(0)
emp = {}
for name, genes in AXES.items():
    gb = [bins.get(g) for g in genes if g in bins.index]
    if not gb:
        continue
    draws = {f: [] for f in FAMS}
    for _ in range(N_SETS):
        pick = [str(rng.choice(pool[k])) for k in gb if pd.notna(k)]
        v = score(pick)
        if v is None:
            continue
        t = pd.DataFrame({"y": v, "d": A.obs.__donor.values, "f": A.obs.__fam.values,
                          "age": A.obs.__age.values, "lg": S.log_genes.values,
                          "sex": A.obs.__sex.values})
        g = t.groupby(["d", "f"]).agg({"y": "mean", "age": "mean", "lg": "mean",
                                       "sex": "first"}).reset_index()
        for f in FAMS:
            dd = g[g.f == f]
            if len(dd) < 8:
                continue
            try:
                m = smf.ols("y ~ age + sex + lg", data=dd).fit()
                draws[f].append(float(m.params["age"]))
            except Exception:
                pass
    for f in FAMS:
        row = R[(R.axis == name) & (R.family == f) & (R.tested == True)]   # noqa: E712
        if len(row) and draws[f]:
            emp[(name, f)] = float(np.mean(np.abs(draws[f]) >= abs(row.beta_age.iloc[0])))
R["empirical_p"] = [emp.get((a, f), np.nan) for a, f in zip(R.axis, R.family)]

# ---------------------------------------------------------------- effective tests
eff = {}
for f in FAMS:
    s = D[D.__fam == f][list(AXES)].dropna(axis=1, how="all").dropna()
    if len(s) >= 10 and s.shape[1] >= 2:
        eff[f] = effective_tests(s.corr(method="spearman"))
EFF = pd.DataFrame(eff).T
EFF.to_csv(OUT / "effective_tests.csv")

# ---------------------------------------------------------------- composition vs age
comp = (S.groupby(["__donor", "__fam"]).size().rename("n").reset_index())
tot = comp.groupby("__donor").n.sum().rename("tot")
comp = comp.join(tot, on="__donor")
comp["frac"] = comp.n / comp.tot
comp = comp.merge(D[["__donor", "age"]].drop_duplicates(), on="__donor", how="left")
crow = []
for f in FAMS:
    d = comp[comp.__fam == f].dropna(subset=["age"])
    if len(d) >= 8:
        m = smf.ols("frac ~ age", data=d).fit()
        crow.append(dict(family=f, beta=float(m.params["age"]), p=float(m.pvalues["age"]),
                         n=len(d)))
CMP = pd.DataFrame(crow)
CMP.to_csv(OUT / "composition_vs_age.csv", index=False)


def verdict(r):
    if not r.tested:
        return r.get("reason", "not tested")
    if np.isfinite(r.empirical_p) and r.empirical_p >= 0.05:
        return "not_specific (fails matched null)"
    if r.BH_age >= 0.05:
        return "no age effect"
    if np.isfinite(r.dAIC) and r.dAIC > 2 and r.BH_quad < 0.05:
        return "NON-MONOTONIC (peak may be described)"
    return "monotonic age effect"


R["verdict"] = R.apply(verdict, axis=1)
R.to_csv(OUT / "trajectory.csv", index=False)

print("\n=== EFFECTIVE TESTS among the five axes (plan 3) ===")
print(EFF.round(2).to_string())
print("\n=== CELL-TYPE COMPOSITION vs AGE (plan 4.3) ===")
print(CMP.round(4).to_string(index=False))
print("\n=== TRAJECTORIES (BH over all tested cells) ===")
print(R.round(4).to_string(index=False))
print("\nwrote", OUT)
