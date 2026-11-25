"""Complete approaches A, B and D of docs/PREREG_method_improvements.md (locked b1db644).

The first pass left A without its residualised matched null (which the plan says voids it),
B with a placeholder inter-gene correlation, and D unrun. This finishes all three.

Usage: python validation/method_improvements_complete.py <data_dir> <results_dir> <out_dir>
"""
import gzip
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp
import statsmodels.formula.api as smf
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from nasage.modules import load_modules      # noqa: E402
from mvri_programs import PROGRAMS           # noqa: E402

DATA, R, OUT = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
OUT.mkdir(parents=True, exist_ok=True)
FAMS = ["Basal", "Secretory", "Ciliated"]
N_SETS, N_BINS = 500, 25
M = load_modules()
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


# ================================================== A: residualised matched null (lifespan)
print("=" * 78)
print("A (completed): matched null recomputed WITH the identical PC1 residualisation")
print("=" * 78)
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


def score(genes):
    present = [g for g in genes if g in A.var_names]
    if not present:
        return None
    Z = A[:, present].X
    Z = np.asarray(Z.todense()) if sp.issparse(Z) else np.asarray(Z)
    m = Z.mean(1)
    sd = m.std(ddof=1)
    return (m - m.mean()) / sd if sd > 0 else m * 0.0


ng = np.asarray((A.X > 0).sum(1)).ravel()
base = pd.DataFrame({"__donor": A.obs.__donor.values, "__fam": A.obs.__fam.values,
                     "age": A.obs.__age.values, "sex": A.obs.__sex.values,
                     "log_genes": np.log(ng + 1)})
for nm, g in AXES.items():
    v = score(g)
    if v is not None:
        base[nm] = v
cnt = base.groupby(["__donor", "__fam"]).size().rename("n_cells")
D = (base.groupby(["__donor", "__fam"])
     .agg({**{c: "mean" for c in list(AXES) + ["log_genes", "age"]}, "sex": "first"})
     .join(cnt).reset_index())
D = D[D.n_cells >= 20]

# PC1 per family, from the real axis scores (the same covariate the real test used)
pc1 = {}
for f in FAMS:
    d = D[D.__fam == f].dropna(subset=list(AXES))
    if len(d) < 10:
        continue
    Xs = d[list(AXES)].values
    Xs = (Xs - Xs.mean(0)) / Xs.std(0, ddof=1)
    u, s, _ = np.linalg.svd(Xs - Xs.mean(0), full_matrices=False)
    pc1[f] = pd.Series(u[:, 0] * s[0], index=d.index)

traj = pd.read_csv(R / "lifespan" / "trajectory.csv")
A1 = pd.read_csv(R / "method_improve" / "A1_lifespan_pc1.csv")

expr = pd.Series(np.asarray(A.X.mean(0)).ravel(), index=A.var_names)
bins = pd.qcut(expr.rank(method="first"), N_BINS, labels=False)
pool = {k: np.array(v) for k, v in bins.groupby(bins).groups.items()}
rng = np.random.default_rng(0)

null_beta = {f: {nm: [] for nm in AXES} for f in FAMS}
for nm, genes in AXES.items():
    gb = [bins.get(g) for g in genes if g in bins.index]
    if not gb:
        continue
    for _ in range(N_SETS):
        pick = [str(rng.choice(pool[k])) for k in gb if pd.notna(k)]
        v = score(pick)
        if v is None:
            continue
        t = base[["__donor", "__fam", "age", "sex", "log_genes"]].copy()
        t["y"] = v
        g = (t.groupby(["__donor", "__fam"])
             .agg({"y": "mean", "age": "mean", "log_genes": "mean", "sex": "first"})
             .reset_index())
        for f in FAMS:
            if f not in pc1:
                continue
            dd = D[D.__fam == f].dropna(subset=list(AXES))
            gg = g[g.__fam == f].set_index("__donor").reindex(
                dd.__donor.values).reset_index()
            gg["__pc1"] = pc1[f].values
            gg = gg.dropna(subset=["y", "age", "log_genes"])
            if len(gg) < 8:
                continue
            try:
                m = smf.ols("y ~ age + sex + log_genes + __pc1", data=gg).fit()
                null_beta[f][nm].append(float(m.params["age"]))
            except Exception:
                pass

arows = []
for _, r in A1.iterrows():
    nb = null_beta.get(r.family, {}).get(r.axis, [])
    emp = float(np.mean(np.abs(nb) >= abs(r.beta_after))) if nb else np.nan
    arows.append(dict(axis=r.axis, family=r.family,
                      beta_before=r.beta_before, emp_before=r.emp_before,
                      beta_after=r.beta_after, p_after=r.p_after,
                      emp_after=emp, n_null=len(nb),
                      survives=bool(np.isfinite(emp) and emp < 0.05 and r.p_after < 0.05)))
AA = pd.DataFrame(arows)
AA.to_csv(OUT / "A_completed_residualised_null.csv", index=False)
print(AA.round(5).to_string(index=False))

# ================================================== B: rho_bar estimated from the data
print("\n" + "=" * 78)
print("B (completed): inter-gene correlation estimated from the Thwaites matrix")
print("=" * 78)
cols = {}
for fp in sorted((DATA / "thwaites2023" / "counts").glob("*_counts.txt.gz")):
    name = re.sub(r"^GSM\d+_|_counts\.txt\.gz$", "", fp.name)
    with gzip.open(fp, "rt") as fh:
        cols[name] = pd.read_csv(fh, sep="\t", header=None, index_col=0,
                                 skiprows=1).iloc[:, 0]
Cm = pd.DataFrame(cols)
Cm.index = [str(i).split(".")[0] for i in Cm.index]
sym = ad.read_h5ad(DATA / "ziegler2021_np.h5ad",
                   backed="r").var["feature_name"].astype(str).reindex(Cm.index)
Cm = Cm[sym.notna()].groupby(sym[sym.notna()].values).sum()
Cm = Cm.loc[:, Cm.sum(0) >= 5e6]
cpm = Cm / Cm.sum(0) * 1e6
logc = np.log2(cpm + 1)
logc = logc[(cpm.mean(1) >= 1)]

eng = pd.read_csv(R / "laiv_engagement" / "engagement.csv")
brows = []
for _, r in eng.iterrows():
    prog = r.program.replace("_epi", "")
    genes = [g for g in PROGRAMS.get(prog, {}).get("genes", []) if g in logc.index]
    if len(genes) < 3 or not np.isfinite(r.d_z):
        continue
    Cc = np.corrcoef(logc.loc[genes].values)
    iu = np.triu_indices(len(genes), 1)
    rbar = float(np.nanmean(Cc[iu]))
    v = max(1.0, 1 + (len(genes) - 1) * rbar)
    t_b = r.d_z * np.sqrt(r.n)
    p_a = float(2 * stats.t.sf(abs(t_b / np.sqrt(v)), r.n - 1))
    brows.append(dict(program=r.program, m_genes=len(genes), rho_bar=rbar, vif=v,
                      d_z=r.d_z, p_before=r.p, p_after=p_a))
B = pd.DataFrame(brows)
B.to_csv(OUT / "B_completed_camera.csv", index=False)
print(B.round(4).to_string(index=False))

# ================================================== D: engagement meta-analysis
print("\n" + "=" * 78)
print("D: meta-analysis, Liverpool as ONE unit (plan 4) + Thwaites")
print("=" * 78)
LI = pd.read_csv(R / "mvri_engagement" / "engagement_interaction.csv")
drows = []
for prog in sorted(LI.program.unique()):
    sub = LI[(LI.program == prog) & LI.beta_interaction.notna() & LI.p.notna()].copy()
    if not len(sub):
        continue
    # SE back-calculated from beta and its two-sided p
    z = stats.norm.isf(sub.p.values / 2)
    se = np.abs(sub.beta_interaction.values) / np.where(z > 0, z, np.nan)
    ok = np.isfinite(se) & (se > 0)
    if ok.sum() == 0:
        continue
    w = 1 / se[ok] ** 2
    liv_b = float((w * sub.beta_interaction.values[ok]).sum() / w.sum())
    liv_se = float(np.sqrt(1 / w.sum()))
    e = eng[eng.program.str.replace("_epi", "", regex=False) == prog]
    if not len(e):
        continue
    er = e.iloc[0]
    if not np.isfinite(er.d_z) or er.n < 3:
        continue
    thw_b = float(er.d_z)
    thw_se = float(np.sqrt(1 / er.n + er.d_z ** 2 / (2 * er.n)))
    bs = np.array([liv_b, thw_b]); ses = np.array([liv_se, thw_se])
    w2 = 1 / ses ** 2
    pooled = float((w2 * bs).sum() / w2.sum())
    pse = float(np.sqrt(1 / w2.sum()))
    Q = float((w2 * (bs - pooled) ** 2).sum())
    I2 = max(0.0, (Q - 1) / Q * 100) if Q > 0 else 0.0
    drows.append(dict(program=prog, liverpool=liv_b, thwaites=thw_b, pooled=pooled,
                      p=float(2 * stats.norm.sf(abs(pooled / pse))), Q=Q, I2=I2))
Dm = pd.DataFrame(drows)
if len(Dm):
    pv = Dm.p.values
    o = np.argsort(pv); q = np.empty_like(pv)
    q[o] = np.minimum.accumulate((pv[o] * len(pv) / np.arange(1, len(pv) + 1))[::-1])[::-1]
    Dm["BH"] = np.minimum(q, 1.0)
Dm.to_csv(OUT / "D_engagement_meta.csv", index=False)
print(Dm.round(4).to_string(index=False))
print("\nNOTE: Liverpool is a LAIV-vs-TIV between-arm contrast and Thwaites an uncontrolled\n"
      "      within-participant change. The pooled value mixes two estimands and is a\n"
      "      direction-only summary, not an effect size.")
print("\nwrote", OUT)
