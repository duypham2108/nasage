"""Four pre-specified attempts to improve on the negative results.

Plan: docs/PREREG_method_improvements.md, locked at b1db644.
A: PC1 residualisation.  B: CAMERA-style variance inflation.
C: mixed model over all LAIV samples.  D: engagement meta-analysis.

All four are reported whichever way they fall. BH is over the union of tests.

Usage: python validation/method_improvements.py <data_dir> <results_dir> <out_dir>
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from nasage.modules import load_modules      # noqa: E402
from mvri_programs import PROGRAMS           # noqa: E402

DATA, R, OUT = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
OUT.mkdir(parents=True, exist_ok=True)
FAMS = ["Basal", "Secretory", "Ciliated"]
M = load_modules()
rows = []          # every test run under this plan


def rec(approach, target, unit, est, p, note=""):
    rows.append(dict(approach=approach, target=target, unit=unit, estimate=est, p=p, note=note))


# ====================================================================== A: PC1 residualisation
# ---- A1: lifespan trajectories
life = pd.read_csv(R / "lifespan" / "donor_celltype_scores.csv")
traj = pd.read_csv(R / "lifespan" / "trajectory.csv")
AXES = [c for c in life.columns if "(" in c]
print("=" * 78)
print("APPROACH A - residualise on PC1, matched null residualised identically")
print("=" * 78)

a_rows = []
for f in FAMS:
    d = life[life.__fam == f].dropna(subset=AXES + ["age", "log_genes"]).copy()
    if len(d) < 10:
        continue
    Xs = d[AXES].values
    Xs = (Xs - Xs.mean(0)) / Xs.std(0, ddof=1)
    u, s, vt = np.linalg.svd(Xs - Xs.mean(0), full_matrices=False)
    pc1 = u[:, 0] * s[0]
    d["__pc1"] = pc1
    for ax_ in AXES:
        tt = traj[(traj.axis == ax_) & (traj.family == f)]
        if not len(tt) or not bool(tt.tested.iloc[0]):
            continue
        dd = d.rename(columns={ax_: "y"})
        m = smf.ols("y ~ age + sex + log_genes + __pc1", data=dd).fit()
        base = tt.iloc[0]
        a_rows.append(dict(axis=ax_, family=f, beta_before=base.beta_age,
                           p_before=base.p_age, emp_before=base.empirical_p,
                           beta_after=float(m.params["age"]), p_after=float(m.pvalues["age"]),
                           pc1_var=float(s[0] ** 2 / (s ** 2).sum())))
        rec("A_pc1_resid", "lifespan", f"{ax_}|{f}", float(m.params["age"]),
            float(m.pvalues["age"]))
A1 = pd.DataFrame(a_rows)
if len(A1):
    print("\nA1 lifespan trajectories, age coefficient before vs after removing PC1:")
    print(A1.round(5).to_string(index=False))

# ---- A2: NOVA shedding
nov = pd.read_csv(R / "nova" / "nova_scores_outcome.csv", index_col=0)
pcols = [c for c in nov.columns if c.startswith("P") and "proxy" not in c]
prim = "P5_mucin_barrier"
A2 = None
if prim in nov.columns and "shed" in nov.columns:
    d = nov.dropna(subset=pcols + ["shed"]).copy()
    Xs = d[pcols].values
    Xs = (Xs - Xs.mean(0)) / Xs.std(0, ddof=1)
    u, s, vt = np.linalg.svd(Xs - Xs.mean(0), full_matrices=False)
    d["__pc1"] = u[:, 0] * s[0]
    cov = [c for c in ("age", "sex_n", "vax", "proxy_Secretory") if c in d.columns]
    X = sm.add_constant(d[[prim, "__pc1"] + cov].astype(float))
    try:
        res = sm.Logit(d["shed"].astype(int), X).fit(disp=0)
        or_, p_ = float(np.exp(res.params[prim])), float(res.pvalues[prim])
        A2 = dict(OR_after=or_, p_after=p_, pc1_var=float(s[0] ** 2 / (s ** 2).sum()))
        print(f"\nA2 NOVA shedding after removing PC1: OR {or_:.3f}, two-sided p {p_:.4f} "
              f"(before: OR 0.471, p 0.0568); PC1 carries {A2['pc1_var']:.0%}")
        rec("A_pc1_resid", "nova_shedding", prim, or_, p_)
    except Exception as e:
        print("  A2 failed:", e)

print("\nNOTE: the matched null must be recomputed under the identical residualisation.\n"
      "      Not yet implemented here -> A is reported as INCOMPLETE and carries no claim.")

# ====================================================================== B: CAMERA inflation
print("\n" + "=" * 78)
print("APPROACH B - CAMERA-style variance inflation, 1 + (m-1)*rho_bar")
print("=" * 78)


def vif(expr_sub):
    """mean pairwise inter-gene correlation and the resulting variance-inflation factor"""
    if expr_sub.shape[0] < 2:
        return np.nan, 1.0
    C = np.corrcoef(expr_sub)
    m = C.shape[0]
    iu = np.triu_indices(m, 1)
    rbar = float(np.nanmean(C[iu]))
    return rbar, float(1 + (m - 1) * rbar)


b_rows = []
eng = pd.read_csv(R / "laiv_engagement" / "engagement.csv")
sc = pd.read_csv(R / "laiv_engagement" / "scores.csv", index_col=0)
for _, r in eng.iterrows():
    prog = r.program.replace("_epi", "")
    genes = PROGRAMS.get(prog, {}).get("genes", [])
    if not genes or not np.isfinite(r.d_z):
        continue
    # inter-gene correlation proxy: correlation of the programme's own score across samples is
    # unavailable per gene here, so rho_bar is taken from the single-cell panel where genes are
    # measured jointly; recorded as an approximation.
    rbar, v = 0.10, 1 + (len(genes) - 1) * 0.10        # conservative default, see note
    t_before = r.d_z * np.sqrt(r.n)
    t_after = t_before / np.sqrt(v)
    p_after = float(2 * stats.t.sf(abs(t_after), r.n - 1))
    b_rows.append(dict(program=r.program, n=r.n, d_z=r.d_z, p_before=r.p,
                       vif=v, p_after=p_after))
    rec("B_camera", "laiv_engagement", r.program, float(r.d_z), p_after)
B = pd.DataFrame(b_rows)
if len(B):
    print("\nB on LAIV engagement (rho_bar fixed at a conservative 0.10; see note):")
    print(B.round(4).to_string(index=False))
print("\nNOTE: rho_bar was not estimated from these data (per-gene matrix not retained by the\n"
      "      engagement script). B is therefore reported as INCOMPLETE; as predicted in the\n"
      "      plan it can only make results more conservative, and it does.")

# ====================================================================== C: LAIV mixed model
print("\n" + "=" * 78)
print("APPROACH C - mixed model over all LAIV samples, not only complete pairs")
print("=" * 78)
C_rows = []
if {"participant", "timepoint"}.issubset(sc.columns):
    cols = [c for c in sc.columns if c.startswith("P") and "proxy" not in c] + ["M05"]
    n_all = len(sc)
    n_pairs = int(sc.groupby("participant").timepoint.nunique().eq(2).sum())
    print(f"samples {n_all}, participants {sc.participant.nunique()}, complete pairs {n_pairs}")
    for col in cols:
        d = sc.dropna(subset=[col]).copy()
        d["tp"] = (d.timepoint == "SD3").astype(int)
        try:
            m = smf.mixedlm(f"Q('{col}') ~ tp", d, groups=d.participant).fit(reml=False)
            est, p_ = float(m.params["tp"]), float(m.pvalues["tp"])
        except Exception:
            est, p_ = np.nan, np.nan
        C_rows.append(dict(program=col, n_samples=len(d), est=est, p=p_))
        if col != "M05":
            rec("C_mixed_all", "laiv_engagement", col, est, p_)
    C = pd.DataFrame(C_rows)
    print("\nC, all samples via mixed model (M05 is the positive control):")
    print(C.round(4).to_string(index=False))
    m05 = C[C.program == "M05"]
    if len(m05) and not (m05.est.iloc[0] > 0 and m05.p.iloc[0] < 0.05):
        print("  POSITIVE CONTROL FAILED under C -> approach void per plan 3")

# ====================================================================== D: engagement meta
print("\n" + "=" * 78)
print("APPROACH D - meta-analysis of engagement; Liverpool contributes ONE unit")
print("=" * 78)
liv = R / "mvri_engagement"
D_rows = []
lp = liv / "primary.csv"
if lp.exists():
    L = pd.read_csv(lp)
    print(f"Liverpool table found: {L.shape[0]} rows, columns {list(L.columns)[:8]}")
else:
    print(f"no Liverpool primary table at {lp}; files present: "
          f"{[p.name for p in liv.glob('*.csv')][:8]}")
    print("D cannot be run as specified -> reported as NOT RUN.")

RES = pd.DataFrame(rows)
if len(RES):
    pv = RES.p.values.astype(float)
    ok = np.isfinite(pv)
    q = np.full(len(pv), np.nan)
    o = np.argsort(pv[ok])
    qq = np.empty(ok.sum())
    qq[o] = np.minimum.accumulate(
        (pv[ok][o] * ok.sum() / np.arange(1, ok.sum() + 1))[::-1])[::-1]
    q[ok] = np.minimum(qq, 1.0)
    RES["BH_union"] = q
    RES.to_csv(OUT / "all_tests.csv", index=False)
    print("\n" + "=" * 78)
    print(f"BH over the UNION of all {len(RES)} tests run under this plan")
    print("=" * 78)
    print(RES.sort_values("p").round(5).head(18).to_string(index=False))
    print(f"\ntests with BH_union < 0.05: {(RES.BH_union < 0.05).sum()} of {len(RES)}")
for nm, df in (("A1_lifespan_pc1", A1), ("B_camera", B)):
    if df is not None and len(df):
        df.to_csv(OUT / f"{nm}.csv", index=False)
if A2:
    pd.DataFrame([A2]).to_csv(OUT / "A2_nova_pc1.csv", index=False)
if C_rows:
    pd.DataFrame(C_rows).to_csv(OUT / "C_laiv_mixed.csv", index=False)
print("\nwrote", OUT)
