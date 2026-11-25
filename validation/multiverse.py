"""Specification curve + axis-orthogonal redesign.

Plan: docs/PREREG_multiverse.md, locked at 51bb606.

Grid: 4 scoring x 4 covariate sets x 3 gene-set variants = 48 specifications per programme per
outcome. The whole distribution is reported, never the best specification.

Usage: python validation/multiverse.py <data_dir> <results_dir> <out_dir>
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
from nasage.modules import load_modules                      # noqa: E402
from mvri_programs import PROGRAMS, PROXY                     # noqa: E402

DATA, R, OUT = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
OUT.mkdir(parents=True, exist_ok=True)
FAMS = ["Basal", "Secretory", "Ciliated"]
MOD = load_modules()
EPI_FRAC = pd.read_csv(R / "module_audit" / "gene_epithelial_fraction.csv"
                       ).set_index("gene")["epi_frac"].to_dict()

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


# ============================================================ scoring methods (plan 1)
def sc_meanz(Xg):
    Z = (Xg - Xg.mean(0)) / np.where(Xg.std(0, ddof=1) > 0, Xg.std(0, ddof=1), 1)
    return Z.mean(1)


def sc_medianz(Xg):
    Z = (Xg - Xg.mean(0)) / np.where(Xg.std(0, ddof=1) > 0, Xg.std(0, ddof=1), 1)
    return np.median(Z, axis=1)


def sc_rank(Xg):
    """mean reciprocal rank within sample -- depth-robust, UCell-style"""
    order = np.argsort(np.argsort(-Xg, axis=1), axis=1) + 1
    return (1.0 / order).mean(1)


def sc_pc1(Xg):
    Z = (Xg - Xg.mean(0)) / np.where(Xg.std(0, ddof=1) > 0, Xg.std(0, ddof=1), 1)
    u, s, vt = np.linalg.svd(Z - Z.mean(0), full_matrices=False)
    v = u[:, 0] * s[0]
    return v * np.sign(np.corrcoef(v, Z.mean(1))[0, 1] or 1.0)


SCORERS = {"mean_z": sc_meanz, "median_z": sc_medianz, "rank": sc_rank, "pc1_of_set": sc_pc1}
# Deviation, 2026-11-25: the locked grid names "depth" as the third axis. NOVA is a microarray
# cohort with no sequencing depth, so for it the slot holds the demographic adjustment set the
# locked NOVA plan used (age, sex, vaccination arm). Recorded rather than silently substituted.
COVSETS = ["none", "proxy", "proxy+demog", "proxy+demog+PC1"]

# ============================================================ axis-orthogonal redesign (plan 2)
print("=" * 78)
print("REDESIGN: axis-orthogonal programmes, built from Ziegler + Loske structure only")
print("=" * 78)


def load_sc(fn, dcol, ccol):
    A = ad.read_h5ad(DATA / fn)
    for c in ("feature_name", "gene_symbols"):
        if c in A.var.columns:
            symv = A.var[c].astype(str)
            if symv.str.startswith("ENSG").mean() < 0.5:
                A.var_names = pd.Index(symv.values)
                A.var_names_make_unique()
                break
    A.obs["__fam"] = A.obs[ccol].map(to_family)
    A = A[A.obs.__fam.notna()].copy()
    A.obs["__donor"] = A.obs[dcol].astype(str)
    X = A.X
    if float(X.max()) > 50:
        tot = np.asarray(X.sum(1)).ravel(); tot[tot == 0] = 1
        A.X = sp.diags(1e4 / tot) @ X
        A.X.data = np.log1p(A.X.data)
    return A


gene_pc1corr = {}
for fn, dcol, ccol in (("ziegler2021_np.h5ad", "donor_id", "Coarse_Cell_Annotations"),
                       ("loske2022_controls.h5ad", "Patient_ID", "detailed_clusters")):
    A = load_sc(fn, dcol, ccol)
    allg = sorted({g for p in PROGRAMS.values() for g in p["genes"] if g in A.var_names})
    Xg = A[:, allg].X
    Xg = np.asarray(Xg.todense()) if sp.issparse(Xg) else np.asarray(Xg)
    df = pd.DataFrame(Xg, columns=allg)
    df["__d"] = A.obs.__donor.values
    df["__f"] = A.obs.__fam.values
    agg = df.groupby(["__d", "__f"]).mean(numeric_only=True)
    for f in FAMS:
        sub = agg.xs(f, level="__f") if f in agg.index.get_level_values("__f") else None
        if sub is None or len(sub) < 10:
            continue
        # PC1 of the PROGRAMME SCORE matrix in this cell type
        S = np.column_stack([sc_meanz(sub[[g for g in p["genes"] if g in sub.columns]].values)
                             for p in PROGRAMS.values()
                             if len([g for g in p["genes"] if g in sub.columns]) >= 3])
        Sz = (S - S.mean(0)) / S.std(0, ddof=1)
        u, s, _ = np.linalg.svd(Sz - Sz.mean(0), full_matrices=False)
        pc = u[:, 0] * s[0]
        for g in sub.columns:
            v = sub[g].values
            if np.std(v) > 0:
                gene_pc1corr.setdefault((f, g), []).append(
                    abs(float(np.corrcoef(v, pc)[0, 1])))
    del A

PC1C = {k: float(np.mean(v)) for k, v in gene_pc1corr.items()}
ORTH = {}
rows_sel = []
for name, p in PROGRAMS.items():
    f = p["cell_type"]
    cand = [(g, PC1C.get((f, g), np.nan)) for g in p["genes"]]
    scored = [(g, c) for g, c in cand if np.isfinite(c)]
    if len(scored) < 3:
        continue
    cs = np.array([c for _, c in scored])
    cut = np.quantile(cs, 2 / 3)                       # lowest two tertiles of |corr|
    keep = [g for g, c in scored
            if c <= cut and EPI_FRAC.get(g, 1.0) >= 0.5]
    rows_sel.append(dict(program=name, n_in=len(p["genes"]), n_scored=len(scored),
                         n_kept=len(keep), cutoff=float(cut),
                         kept=",".join(sorted(keep)),
                         dropped=",".join(sorted(set(p["genes"]) - set(keep)))))
    if len(keep) >= 8:
        ORTH[name + "_orth"] = dict(cell_type=f, genes=keep)
SEL = pd.DataFrame(rows_sel)
SEL.to_csv(OUT / "orthogonal_selection.csv", index=False)
print(SEL[["program", "n_in", "n_scored", "n_kept"]].to_string(index=False))
print("retained orthogonal programmes:", sorted(ORTH) or "NONE")

VARIANTS = {"original": PROGRAMS}
try:
    from mvri_programs_epi import PROGRAMS_EPI
    VARIANTS["epithelial"] = PROGRAMS_EPI
except Exception:
    pass
if ORTH:
    VARIANTS["orthogonal"] = ORTH

# ============================================================ NOVA specification curve
print("\n" + "=" * 78)
print("SPECIFICATION CURVE on NOVA shedding")
print("=" * 78)
nd = Path(DATA / "gse253686")
EXPR = R / "nova" / "expr_symbols.csv.gz"
rows = []
try:
    nov = pd.read_csv(R / "nova" / "nova_scores_outcome.csv", index_col=0)
    have_nova = "shed" in nov.columns
except Exception:
    have_nova = False

if have_nova and EXPR.exists():
    E = pd.read_csv(EXPR, index_col=0)
    meta = nov
    common = [c for c in E.columns if c in meta.index]
    E = E[common]
    meta = meta.loc[common]
    glob_pc1 = None
    for vname, VP in VARIANTS.items():
        for pname, p in VP.items():
            genes = [g for g in p["genes"] if g in E.index]
            if len(genes) < 3:
                continue
            Xg = E.loc[genes].T.values
            for sname, fn_ in SCORERS.items():
                v = fn_(Xg)
                v = (v - v.mean()) / (v.std(ddof=1) or 1)
                d = pd.DataFrame({"y": meta["shed"].astype(int).values, "s": v},
                                 index=meta.index)
                for cname in COVSETS:
                    X = pd.DataFrame({"s": d.s})
                    if "proxy" in cname:
                        pc = f"proxy_{p['cell_type']}"
                        if pc in meta.columns:
                            X[pc] = meta[pc].values
                    if "demog" in cname:
                        for cc in ("Age", "male", "vax"):
                            if cc in meta.columns:
                                X[cc] = pd.to_numeric(meta[cc], errors="coerce").values
                    if "PC1" in cname:
                        if glob_pc1 is None:
                            Z = E.T.values
                            Z = (Z - Z.mean(0)) / np.where(Z.std(0, ddof=1) > 0,
                                                           Z.std(0, ddof=1), 1)
                            uu, ss, _ = np.linalg.svd(Z - Z.mean(0), full_matrices=False)
                            glob_pc1 = uu[:, 0] * ss[0]
                        X["gpc1"] = glob_pc1
                    try:
                        Xf = sm.add_constant(X.astype(float)).dropna()
                        res = sm.Logit(d.y.loc[Xf.index], Xf).fit(disp=0)
                        rows.append(dict(target="nova_shedding", variant=vname,
                                         program=pname, scoring=sname, covars=cname,
                                         est=float(res.params["s"]),
                                         p=float(res.pvalues["s"])))
                    except Exception:
                        pass
else:
    print(f"NOVA per-gene matrix not available at {EXPR}; "
          "specification curve for NOVA reported as NOT RUN.")

SC = pd.DataFrame(rows)
if len(SC):
    SC.to_csv(OUT / "specification_curve.csv", index=False)
    print(f"\n{len(SC)} specifications run")
    for (tg, pr), g in SC.groupby(["target", "program"]):
        print(f"\n--- {pr} ({tg}), {len(g)} specs")
        print(f"    median effect {g.est.median():+.3f} | "
              f"{100 * (g.p < 0.05).mean():.0f}% reach p<0.05 | "
              f"{100 * (np.sign(g.est) == np.sign(g.est.median())).mean():.0f}% same sign")
print("\nwrote", OUT)
