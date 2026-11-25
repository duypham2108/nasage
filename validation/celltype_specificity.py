"""Analysis A of docs/PREREG_celltype_mvri.md, locked at 77701ef.

Is each of the seven curated programmes highest in the cell type it was curated for?
Unit of analysis is the donor x cell type, never the cell.

Usage: python validation/celltype_specificity.py <data_dir> <out_dir>
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from mvri_programs import PROGRAMS, PROXY  # noqa: E402

DATA, OUT = Path(sys.argv[1]), Path(sys.argv[2])
OUT.mkdir(parents=True, exist_ok=True)

EXCLUDE = {"P2_junction_integrity"}          # locked: Jaccard 0.824 with M01
MIN_CELLS = 20                               # plan 2
FAMILIES = ("Basal", "Secretory", "Ciliated")
N_SETS, N_BINS, N_PERM = 500, 25, 2000

# Cohort -> (file, donor column, cell-type column). Crownlands is excluded by the plan.
COHORTS = {
    "Ziegler2021": ("ziegler2021_np.h5ad", "donor_id", "Coarse_Cell_Annotations"),
    "Loske2022": ("loske2022_controls.h5ad", "Patient_ID", "detailed_clusters"),
    "Xu2022": ("xu2022_turbinate.h5ad", "donor_id", "cell_type"),
}

# Annotation strings -> the three epithelial families. Fixed before any score was computed.
# Deuterosomal/precursor (FOXN4), squamous, ionocyte, low-RNA and immune labels belong to no
# family and are dropped rather than forced into one.
BASAL_K = ("basal",)
CILIATED_K = ("ciliated", "cil_diff", "cil_acute", "multiciliated")
SECRETORY_K = ("secretory", "goblet", "club", "glandular", "sec_diff")
DROP_K = ("foxn4", "squamous", "ionocyte", "lowrna", "low_rna", "deutero")

# Map each cohort's annotation strings onto the three epithelial families.
# Goblet cells are secretory by lineage and are grouped there; squamous, immune and
# ionocytes are not one of the three families and are dropped.
def to_family(lbl):
    s = str(lbl).lower().replace(" ", "_")
    if any(k in s for k in DROP_K):
        return None
    if any(k in s for k in BASAL_K):
        return "Basal"
    if any(k in s for k in CILIATED_K):
        return "Ciliated"
    if any(k in s for k in SECRETORY_K):
        return "Secretory"
    return None


def to_symbols(A):
    """Return A with var_names as HGNC symbols; these objects ship Ensembl IDs + a symbol col."""
    for c in ("feature_name", "gene_symbols", "features"):
        if c in A.var.columns:
            sym = A.var[c].astype(str)
            if sym.str.startswith("ENSG").mean() < 0.5:
                A.var_names = pd.Index(sym.values)
                A.var_names_make_unique()
                return A
    return A


def ucell_scores(A, genes, bg_ranks=None):
    """Mean reciprocal-rank score per cell (UCell-style), computed on ranks within cell."""
    present = [g for g in genes if g in A.var_names]
    if not present:
        return None, present
    X = A[:, present].X
    X = np.asarray(X.todense()) if sp.issparse(X) else np.asarray(X)
    return X, present


def score_cells(A, genes):
    """z-of-mean-logexpr per cell, z taken across cells within the cohort."""
    present = [g for g in genes if g in A.var_names]
    if not present:
        return None, present
    X = A[:, present].X
    X = np.asarray(X.todense()) if sp.issparse(X) else np.asarray(X)
    m = X.mean(1)
    sd = m.std(ddof=1)
    return ((m - m.mean()) / sd if sd > 0 else m * 0.0), present


def detection(A, genes):
    present = [g for g in genes if g in A.var_names]
    if not present:
        return 0.0
    X = A[:, present].X
    frac = np.asarray((X > 0).mean(0)).ravel() if sp.issparse(X) else (np.asarray(X) > 0).mean(0)
    return float((frac >= 0.05).mean())


rows, cover_rows, null_rows = [], [], []
for coh, (fn, dcol, ccol) in COHORTS.items():
    p = DATA / fn
    if not p.exists():
        print(f"  {coh}: file missing, skipped ({p})")
        continue
    A = ad.read_h5ad(p)
    A = to_symbols(A)
    # pick donor / cell-type columns if not declared
    if dcol is None:
        for c in A.obs.columns:
            if c.lower() in ("donor_id", "donor", "patient", "sample", "subject"):
                dcol = c
                break
    if ccol is None:
        for c in ("Coarse_Cell_Annotations", "cell_type", "celltype", "annotation"):
            if c in A.obs.columns:
                ccol = c
                break
    if dcol is None or ccol is None:
        print(f"  {coh}: no donor/celltype column, skipped")
        continue

    fam = A.obs[ccol].map(to_family)
    keep = fam.notna().values
    A = A[keep].copy()
    A.obs["__fam"] = fam[keep].values
    A.obs["__donor"] = A.obs[dcol].astype(str).values

    # normalise to CP10K + log1p if the matrix looks like raw counts
    X = A.X
    mx = float(X.max())
    if mx > 50:
        tot = np.asarray(X.sum(1)).ravel()
        tot[tot == 0] = 1
        if sp.issparse(X):
            A.X = sp.diags(1e4 / tot) @ X
            A.X.data = np.log1p(A.X.data)
        else:
            A.X = np.log1p(X * (1e4 / tot)[:, None])

    print(f"  {coh}: {A.n_obs} epithelial cells, {A.obs.__donor.nunique()} donors, "
          f"families {sorted(A.obs.__fam.unique())}")

    # ---- measurability gate, reported before any model (plan 2 / 7.2)
    for name, prog in PROGRAMS.items():
        for f in FAMILIES:
            sub = A[A.obs.__fam == f]
            if sub.n_obs == 0:
                continue
            cover_rows.append(dict(cohort=coh, program=name, family=f,
                                   n_cells=sub.n_obs, frac_detected=detection(sub, prog["genes"])))

    # ---- per-cell scores, then collapse to donor x family
    sc = {}
    for name, prog in PROGRAMS.items():
        v, present = score_cells(A, prog["genes"])
        if v is not None:
            sc[name] = v
    ngene = np.asarray((A.X > 0).sum(1)).ravel()
    S = pd.DataFrame(sc, index=A.obs_names)
    S["__donor"] = A.obs.__donor.values
    S["__fam"] = A.obs.__fam.values
    S["log_genes"] = np.log(ngene + 1)
    cnt = S.groupby(["__donor", "__fam"]).size().rename("n_cells")
    agg = S.groupby(["__donor", "__fam"]).mean(numeric_only=True).join(cnt)
    agg = agg[agg.n_cells >= MIN_CELLS].reset_index()
    agg["cohort"] = coh
    rows.append(agg)

    # ---- matched null: 500 expression-matched random sets, same contrast
    expr = np.asarray(A.X.mean(0)).ravel()
    order = pd.Series(expr, index=A.var_names).rank(method="first")
    bins = pd.qcut(order, N_BINS, labels=False)
    pool = {k: np.array(v) for k, v in bins.groupby(bins).groups.items()}
    rng = np.random.default_rng(0)
    for name, prog in PROGRAMS.items():
        if name in EXCLUDE:
            continue
        own = prog["cell_type"]
        gb = [bins.get(g) for g in prog["genes"] if g in bins.index]
        if not gb:
            continue
        for _ in range(N_SETS):
            pick = [str(rng.choice(pool[k])) for k in gb if pd.notna(k)]
            v, _ = score_cells(A, pick)
            if v is None:
                continue
            t = pd.DataFrame({"v": v, "d": A.obs.__donor.values, "f": A.obs.__fam.values})
            g = t.groupby(["d", "f"]).v.mean().unstack()
            if own not in g.columns:
                continue
            oth = [c for c in FAMILIES if c in g.columns and c != own]
            if not oth:
                continue
            d = (g[own] - g[oth].mean(1)).dropna()
            if len(d) >= 3:
                null_rows.append(dict(cohort=coh, program=name, delta=float(d.mean())))
    del A

if not rows:
    sys.exit("no cohort loaded")

AG = pd.concat(rows, ignore_index=True)
AG.to_csv(OUT / "donor_celltype_scores.csv", index=False)
COV = pd.DataFrame(cover_rows)
COV.to_csv(OUT / "coverage.csv", index=False)

print("\n=== MEASURABILITY, before any model (plan 2) ===")
piv = COV.pivot_table(index="program", columns=["cohort", "family"], values="frac_detected")
print(piv.round(2).to_string())

# ---------------------------------------------------------------- Analysis A tests
res = []
for name, prog in PROGRAMS.items():
    own = prog["cell_type"]
    for coh, sub in AG.groupby("cohort"):
        g = sub.pivot_table(index="__donor", columns="__fam", values=name)
        if own not in g.columns:
            continue
        oth = [c for c in FAMILIES if c in g.columns and c != own]
        if not oth:
            continue
        d = (g[own] - g[oth].mean(1)).dropna()
        if len(d) < 3:
            continue
        w = stats.wilcoxon(d)[1] if len(d) >= 6 else np.nan
        smd = float(d.mean() / d.std(ddof=1)) if d.std(ddof=1) > 0 else np.nan
        se = float(np.sqrt(1 / len(d) + smd ** 2 / (2 * len(d)))) if np.isfinite(smd) else np.nan
        res.append(dict(program=name, cohort=coh, own=own, n_donors=len(d),
                        delta=float(d.mean()), smd=smd, se=se, p_wilcoxon=w))
R = pd.DataFrame(res)

meta = []
for name, sub in R.groupby("program"):
    s = sub.dropna(subset=["smd", "se"])
    if not len(s):
        continue
    w = 1 / s.se.values ** 2
    pooled = float((w * s.smd.values).sum() / w.sum())
    pse = float(np.sqrt(1 / w.sum()))
    Q = float((w * (s.smd.values - pooled) ** 2).sum())
    I2 = max(0.0, (Q - (len(s) - 1)) / Q * 100) if Q > 0 and len(s) > 1 else 0.0
    pm = float(2 * stats.norm.sf(abs(pooled / pse)))
    nl = [r["delta"] for r in null_rows if r["program"] == name]
    emp = float(np.mean(np.abs(nl) >= abs(s.delta.mean()))) if nl else np.nan
    meta.append(dict(program=name, smd=pooled, ci_lo=pooled - 1.96 * pse,
                     ci_hi=pooled + 1.96 * pse, p=pm, I2=I2, n_cohorts=len(s),
                     all_same_sign=bool((np.sign(s.smd) == np.sign(pooled)).all()),
                     empirical_p=emp, n_null=len(nl)))
M = pd.DataFrame(meta)
fam_m = M[~M.program.isin(EXCLUDE)].copy()
pv = fam_m.p.values
o = np.argsort(pv); q = np.empty_like(pv)
q[o] = np.minimum.accumulate((pv[o] * len(pv) / np.arange(1, len(pv) + 1))[::-1])[::-1]
fam_m["BH"] = np.minimum(q, 1.0)
M = M.merge(fam_m[["program", "BH"]], on="program", how="left")


def verdict(r):
    if not np.isfinite(r.smd):
        return "not_tested"
    if r.smd < 0 and r.BH < 0.05:
        return "MISASSIGNED"
    if r.BH >= 0.05 or not np.isfinite(r.BH):
        return "not demonstrated"
    if np.isfinite(r.empirical_p) and r.empirical_p >= 0.05:
        return "not_specific (fails matched null)"
    return "SPECIFIC" if r.all_same_sign else "specific but heterogeneous"


M["verdict"] = M.apply(verdict, axis=1)
R.to_csv(OUT / "per_cohort.csv", index=False)
M.to_csv(OUT / "meta.csv", index=False)
pd.DataFrame(null_rows).to_csv(OUT / "matched_null_draws.csv", index=False)

print("\n=== ANALYSIS A: own cell type vs the other two ===")
print(R.round(3).to_string(index=False))
print("\n=== META (BH over the family of 6; P2 excluded a priori) ===")
print(M.round(4).to_string(index=False))
print("\nwrote", OUT)
