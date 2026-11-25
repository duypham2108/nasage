"""Adult-portion replication of the lifespan trajectories.

Closes deviation 1 of results/lifespan/PREREG_RESULTS.md: plan
docs/PREREG_lifespan_trajectory.md §5 declared a replication in Ziegler 2021 (decade bins) and
GSE336399 (47-89, continuous) that was never run.

Criterion, as locked: agreement in DIRECTION. Magnitudes are not expected to match across
platforms, so they are reported but not compared.

Usage: python validation/lifespan_replication.py <data_dir> <results_dir> <out_dir>
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
from nasage.modules import load_modules       # noqa: E402
from mvri_programs import PROGRAMS            # noqa: E402

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
ADULT_MIN = 18
rows = []

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


# ============================================================ 1. Loske adult portion (reference)
print("=" * 78)
print("REFERENCE: Loske adult portion (>= 18 y) -- the estimate being replicated")
print("=" * 78)
life = pd.read_csv(R / "lifespan" / "donor_celltype_scores.csv")
ad_l = life[life.age >= ADULT_MIN]
print(f"{ad_l.__donor.nunique()} adult donors, {len(ad_l)} donor x cell-type units")
for nm in AXES:
    if nm not in ad_l.columns:
        continue
    for f in FAMS:
        d = ad_l[ad_l.__fam == f].dropna(subset=[nm, "age", "log_genes"])
        if len(d) < 8:
            continue
        m = smf.ols(f"Q('{nm}') ~ age + sex + log_genes", data=d).fit()
        rows.append(dict(cohort="Loske_adult", axis=nm, family=f, n=len(d),
                         beta=float(m.params["age"]), p=float(m.pvalues["age"])))

# ============================================================ 2. Ziegler, decade bins
print("\n" + "=" * 78)
print("REPLICATION 1: Ziegler 2021, decade bins (ordinal), donor x cell type")
print("=" * 78)
A = ad.read_h5ad(DATA / "ziegler2021_np.h5ad")
symv = A.var["feature_name"].astype(str)
A.var_names = pd.Index(symv.values)
A.var_names_make_unique()
A.obs["__fam"] = A.obs["Coarse_Cell_Annotations"].map(to_family)
A = A[A.obs.__fam.notna()].copy()
A.obs["__donor"] = A.obs["donor_id"].astype(str)
dec = A.obs["development_stage"].astype(str).str.extract(
    r"(third|fourth|fifth|sixth|seventh|eighth|ninth)", expand=False)
DECMAP = {"third": 3, "fourth": 4, "fifth": 5, "sixth": 6, "seventh": 7, "eighth": 8,
          "ninth": 9}
A.obs["__dec"] = dec.map(DECMAP).astype(float)
A = A[A.obs.__dec.notna()].copy()
X = A.X
if float(X.max()) > 50:
    tot = np.asarray(X.sum(1)).ravel(); tot[tot == 0] = 1
    A.X = sp.diags(1e4 / tot) @ X
    A.X.data = np.log1p(A.X.data)


def cellscore(Aobj, genes):
    present = [g for g in genes if g in Aobj.var_names]
    if len(present) < 3:
        return None, 0.0
    Z = Aobj[:, present].X
    Z = np.asarray(Z.todense()) if sp.issparse(Z) else np.asarray(Z)
    frac = float(((Z > 0).mean(0) >= 0.05).mean())
    m = Z.mean(1)
    sd = m.std(ddof=1)
    return ((m - m.mean()) / sd if sd > 0 else m * 0.0), frac


ngz = np.asarray((A.X > 0).sum(1)).ravel()
base = pd.DataFrame({"__donor": A.obs.__donor.values, "__fam": A.obs.__fam.values,
                     "dec": A.obs.__dec.values,
                     "sex": A.obs["sex"].astype(str).values if "sex" in A.obs else "u",
                     "log_genes": np.log(ngz + 1)})
gate_z = {}
for nm, g in AXES.items():
    v, frac = cellscore(A, g)
    gate_z[nm] = frac
    if v is not None:
        base[nm] = v
cnt = base.groupby(["__donor", "__fam"]).size().rename("n_cells")
Dz = (base.groupby(["__donor", "__fam"])
      .agg({**{c: "mean" for c in list(AXES) if c in base.columns},
            "dec": "mean", "log_genes": "mean", "sex": "first"})
      .join(cnt).reset_index())
Dz = Dz[Dz.n_cells >= 20]
print(f"{Dz.__donor.nunique()} donors, {len(Dz)} donor x cell-type units; "
      f"decades {int(Dz.dec.min())}-{int(Dz.dec.max())}")
print("measurability:", {k: round(v, 2) for k, v in gate_z.items()})
for nm in AXES:
    if nm not in Dz.columns or gate_z.get(nm, 0) < 0.60:
        continue
    for f in FAMS:
        d = Dz[Dz.__fam == f].dropna(subset=[nm, "dec", "log_genes"])
        if len(d) < 8:
            continue
        m = smf.ols(f"Q('{nm}') ~ dec + sex + log_genes", data=d).fit()
        rows.append(dict(cohort="Ziegler", axis=nm, family=f, n=len(d),
                         beta=float(m.params["dec"]), p=float(m.pvalues["dec"])))
del A

# ============================================================ 3. GSE336399, continuous 47-89
print("\n" + "=" * 78)
print("REPLICATION 2: GSE336399, bulk, continuous age 47-89")
print("=" * 78)
meta = pd.read_csv(R / "gse336399" / "module_scores.csv")
cnts = pd.read_csv(DATA / "gse336399" / "counts.tsv.gz", sep="\t", index_col=0)
# counts ship Ensembl IDs; map to symbols with the same reference the other scripts use
cnts.index = [str(i).split(".")[0] for i in cnts.index]
_sym = ad.read_h5ad(DATA / "ziegler2021_np.h5ad",
                    backed="r").var["feature_name"].astype(str).reindex(cnts.index)
print(f"  {_sym.notna().mean():.1%} of Ensembl IDs mapped to symbols")
cnts = cnts[_sym.notna()].groupby(_sym[_sym.notna()].values).sum()
cpm = cnts / cnts.sum(0) * 1e6
logc = np.log2(cpm + 1)
keep = (cpm.mean(1) >= 1) & (logc.std(1, ddof=1) > 0)
Zb = logc[keep].sub(logc[keep].mean(1), axis=0).div(logc[keep].std(1, ddof=1), axis=0)
# counts columns are the library IDs held in the unnamed first column of module_scores.csv
idcol = "Unnamed: 0" if "Unnamed: 0" in meta.columns else meta.columns[0]
meta = meta.set_index(meta[idcol].astype(str))
common = [c for c in Zb.columns if str(c) in meta.index]
print(f"{Zb.shape[0]} genes, {len(common)} samples matched to metadata")
if common:
    mm = meta.loc[[str(c) for c in common]]
    for nm, g in AXES.items():
        gg = [x for x in g if x in Zb.index]
        frac = len(gg) / len(g)
        if frac < 0.60:
            print(f"  {nm}: {frac:.2f} coverage - fails gate, not tested")
            continue
        s = Zb.loc[gg, common].mean(0).values
        d = pd.DataFrame({"y": (s - s.mean()) / s.std(ddof=1),
                          "age": pd.to_numeric(mm["age"].values, errors="coerce"),
                          "sex": mm["gender"].astype(str).values,
                          "depth": pd.to_numeric(mm["log10_depth"].values,
                                                 errors="coerce")}).dropna()
        if len(d) < 10:
            continue
        m = smf.ols("y ~ age + sex + depth", data=d).fit()
        rows.append(dict(cohort="GSE336399", axis=nm, family="bulk", n=len(d),
                         beta=float(m.params["age"]), p=float(m.pvalues["age"])))

# ============================================================ verdicts
RES = pd.DataFrame(rows)
RES.to_csv(OUT / "replication_estimates.csv", index=False)
print("\n" + "=" * 78)
print("ALL ESTIMATES")
print("=" * 78)
print(RES.round(5).to_string(index=False))

print("\n" + "=" * 78)
print("DIRECTION AGREEMENT (the locked criterion)")
print("=" * 78)
ver = []
for nm in AXES:
    ref = RES[(RES.cohort == "Loske_adult") & (RES.axis == nm)]
    if not len(ref):
        continue
    ref_sign = np.sign(ref.beta.mean())
    for coh in ("Ziegler", "GSE336399"):
        sub = RES[(RES.cohort == coh) & (RES.axis == nm)]
        if not len(sub):
            ver.append(dict(axis=nm, cohort=coh, status="not tested"))
            continue
        agree = float((np.sign(sub.beta) == ref_sign).mean())
        ver.append(dict(axis=nm, cohort=coh, n_cells=len(sub),
                        ref_sign="neg" if ref_sign < 0 else "pos",
                        frac_same_sign=agree,
                        min_p=float(sub.p.min()),
                        replicates=bool(agree >= 0.5)))
V = pd.DataFrame(ver)
V.to_csv(OUT / "direction_agreement.csv", index=False)
print(V.round(3).to_string(index=False))
print("\nwrote", OUT)
