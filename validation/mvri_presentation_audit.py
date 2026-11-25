"""Audit of the MVRI slide deck's analysis, in the primary cohort.

Implements docs/PREREG_mvri_presentation_audit.md. The five axes and their cell subsets are
transcribed verbatim from slides 9-13; the weights from slide 14. Nothing here is a new
analysis of the biology -- it re-runs the deck's own analysis with the unit of analysis
corrected to the donor, a measurability gate, expression-matched nulls, and a weight
sensitivity sweep.

Usage: python validation/mvri_presentation_audit.py <h5ad> <out_dir> [n_null]
"""
import sys
import zlib
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats

H5, OUT = Path(sys.argv[1]), Path(sys.argv[2])
NNULL = int(sys.argv[3]) if len(sys.argv) > 3 else 500
OUT.mkdir(parents=True, exist_ok=True)

AGE_ORDER = ["Neonate", "Infant", "Young child", "Child", "Adolescent", "Adult", "Eldly"]

# ---------------------------------------------------------------- the deck, transcribed
AXES = {
    "Mucus permeability": dict(
        genes=["MUC1", "MUC4", "MUC5AC", "MUC5B", "MUC16",
               "GALNT4", "PLAU", "BPIFA2", "LYZ", "SCGB1A1",
               "CLDN1", "CLDN4", "OCLN", "TJP1", "DSG3"],
        subset=["Secretory", "Goblet", "Club", "Basal"], weight=-0.20),
    "Antigen presentation": dict(
        genes=["FCGR1A", "FCGR2A", "FCGR2B", "FCGR3A", "SCARB1", "LRP1", "LY6E",
               "HLA-A", "HLA-B", "HLA-C", "HLA-DRA", "HLA-DRB1", "CD74",
               "CD40", "CD80", "IFNGR1", "IFNGR2"],
        subset=["Secretory", "Goblet", "DC", "Mac", "Mono"], weight=+0.25),
    "Innate sensing": dict(
        genes=["TLR1", "TLR2", "TLR3", "TLR4", "TLR5", "TLR6", "TLR7", "TLR8", "TLR9",
               "TLR10", "DDX58", "IFIH1", "CGAS", "TMEM173", "IRF3", "IRF7", "IRF9",
               "NFKB1", "NFKB2", "RELA"],
        subset=["__EPI__", "DC", "Mac", "Mono"], weight=+0.25),
    "IgA support": dict(
        genes=["PIGR", "TGFB1", "TNFSF13", "TNFSF13B"],
        subset=["Secretory", "Goblet", "Plasma", "B mem", "DC"], weight=+0.20),
    "Cytokine tone": dict(
        genes=["IL33", "TSLP", "IL6", "CXCL8", "IL8", "CXCL1", "CXCL10", "CCL20",
               "ISG15", "IFIT1", "IFIT2", "IFIT3", "MX1", "OAS1"],
        subset=["__EPI__", "__MYELOID__"], weight=+0.10),
}
# symbol fallbacks permitted by the plan
ALIASES = {"CGAS": ["CGAS", "MB21D1", "C6orf150"], "TMEM173": ["TMEM173", "STING1"],
           "CXCL8": ["CXCL8", "IL8"], "IL8": ["CXCL8", "IL8"]}


def seed_of(*p):
    return zlib.crc32("|".join(map(str, p)).encode()) % (2 ** 32)


print("reading", H5.name)
A = ad.read_h5ad(H5)
print(f"  {A.n_obs} cells x {A.n_vars} genes, {A.obs.donor.nunique()} donors")

# normalise to log CPM-10k once
X = A.X
if not sp.issparse(X):
    X = sp.csr_matrix(X)
tot = np.asarray(X.sum(1)).ravel()
tot[tot == 0] = 1
X = sp.diags(1e4 / tot) @ X
X = X.tocsr()
X.data = np.log1p(X.data)

var = pd.Index(A.var_names)
pos = {g: i for i, g in enumerate(var)}
L2 = A.obs["Cell_type_annotation_level2"].astype(str).values
L1 = A.obs["Cell_type_annotation_level1"].astype(str).values
donor = A.obs["donor"].astype(str).values
agebin = A.obs["Age_bin"].astype(str).values
print("  level2 labels:", sorted(set(L2))[:40])


def subset_mask(spec):
    m = np.zeros(A.n_obs, bool)
    for s in spec:
        if s == "__EPI__":
            m |= (L1 == "Epi")
        elif s == "__MYELOID__":
            m |= (L1 == "Myeloid")
        else:
            m |= np.char.startswith(L2.astype(str), s)
    return m


def resolve(genes):
    out = []
    for g in genes:
        for cand in ALIASES.get(g, [g]):
            if cand in pos:
                if pos[cand] not in out:
                    out.append(pos[cand])
                break
    return out


def score_cells(idx, Xs):
    """Mean of per-gene z, computed within the subset -- the deck's style of score.

    Takes the already-sliced subset matrix; re-slicing X per call exhausted memory on the
    largest subset (Epi + Myeloid, 166M nonzeros).
    """
    M = np.asarray(Xs[:, idx].todense())
    sd = M.std(0, ddof=1)
    return ((M - M.mean(0)) / np.where(sd > 0, sd, 1)).mean(1)


def kw(groups):
    gs = [g for g in groups if len(g) > 0]
    if len(gs) < 3:
        return np.nan, np.nan
    try:
        s, p = stats.kruskal(*gs)
        return float(s), float(p)
    except Exception:
        return np.nan, np.nan


def by_bin(values, keys):
    return [values[keys == b] for b in AGE_ORDER if (keys == b).sum() > 0]


# ---------------------------------------------------------------- A, B, C per axis
rows_out, donor_tbl = [], {}
for name, ax in AXES.items():
    rows = np.where(subset_mask(ax["subset"]))[0]
    idx = resolve(ax["genes"])
    n_req = len({g for g in ax["genes"] if g not in ("IL8",)})
    Xs = X[rows].tocsc()                      # slice ONCE; reused for score, gate and null
    det = np.asarray((Xs[:, idx] > 0).sum(0)).ravel() / len(rows)
    frac_meas = float((det >= 0.05).mean()) if len(idx) else 0.0
    gate = frac_meas >= 0.60

    v = score_cells(idx, Xs)
    dkey, akey = donor[rows], agebin[rows]
    # per-cell test (what the deck's boxplots imply)
    s_cell, p_cell = kw(by_bin(v, akey))
    # per-donor test (the valid unit)
    d = (pd.DataFrame({"donor": dkey, "age": akey, "v": v})
         .groupby(["donor", "age"], observed=True).v.mean().reset_index())
    s_don, p_don = kw(by_bin(d.v.values, d.age.values))
    donor_tbl[name] = d.set_index("donor")

    # C: expression-matched null at donor level
    mean_all = np.asarray(Xs.sum(0)).ravel() / len(rows)
    bins = pd.qcut(pd.Series(mean_all).rank(method="first"), 25, labels=False)
    pool = {k: np.array(g) for k, g in bins.groupby(bins).groups.items()}
    gb = [bins.iloc[i] for i in idx]
    rng = np.random.default_rng(seed_of("mvri_audit", name))
    nstat, npv = [], []
    for _ in range(NNULL):
        pick = [int(rng.choice(pool[k])) for k in gb if pd.notna(k)]
        if len(pick) < 3:
            continue
        vn = score_cells(pick, Xs)
        dn = (pd.DataFrame({"donor": dkey, "age": akey, "v": vn})
              .groupby(["donor", "age"], observed=True).v.mean().reset_index())
        sn, pn = kw(by_bin(dn.v.values, dn.age.values))
        if np.isfinite(sn):
            nstat.append(sn); npv.append(pn)
    nstat, npv = np.asarray(nstat), np.asarray(npv)
    emp = float((nstat >= s_don).mean()) if nstat.size else np.nan

    rows_out.append(dict(axis=name, n_genes_listed=n_req, n_genes_found=len(idx),
                         n_cells=len(rows), n_donors=d.donor.nunique(),
                         frac_measurable=frac_meas, gate_pass=bool(gate),
                         KW_p_percell=p_cell, KW_p_perdonor=p_don,
                         p_ratio=(p_don / p_cell) if p_cell and p_cell > 0 else np.nan,
                         KW_stat_perdonor=s_don, empirical_p=emp,
                         null_frac_p05=float((npv < 0.05).mean()) if npv.size else np.nan))
    r = rows_out[-1]
    print(f"{name:22s} genes {len(idx):2d}/{n_req:2d} meas {frac_meas:.2f} "
          f"{'PASS' if gate else 'FAIL'} | p_cell {p_cell:.2e} p_donor {p_don:.4f} "
          f"| emp_p {emp:.3f} null05 {r['null_frac_p05']:.3f}")
    del Xs                                    # free the subset before the next axis

R = pd.DataFrame(rows_out)
ok = R.KW_p_perdonor.notna()
pv = R.loc[ok, "KW_p_perdonor"].values
order = np.argsort(pv)
q = np.empty_like(pv)
q[order] = np.minimum.accumulate(
    (pv[order] * len(pv) / np.arange(1, len(pv) + 1))[::-1])[::-1]
R.loc[ok, "BH_perdonor"] = np.minimum(q, 1.0)
R["age_varying"] = R.BH_perdonor < 0.05
R["specific"] = R.empirical_p < 0.05
R.to_csv(OUT / "axis_audit.csv", index=False)

# ---------------------------------------------------------------- D: weight sensitivity
S = pd.DataFrame({k: v.groupby("donor").v.mean() for k, v in donor_tbl.items()})
age_of = (pd.DataFrame({"donor": donor, "age": agebin})
          .drop_duplicates().set_index("donor").age)
S["age"] = age_of.reindex(S.index).values
S = S.dropna()
Z = S[list(AXES)].apply(lambda c: (c - c.mean()) / c.std(ddof=1))
binmeans = Z.assign(age=S.age.values).groupby("age").mean()
binmeans = binmeans.reindex([b for b in AGE_ORDER if b in binmeans.index])
binmeans.to_csv(OUT / "axis_means_by_agebin.csv")
S.to_csv(OUT / "donor_axis_scores.csv")       # donor x axis, for the depth sensitivity

w_deck = np.array([AXES[k]["weight"] for k in AXES])
deck_curve = binmeans.values @ w_deck
deck_peak = binmeans.index[int(np.argmax(deck_curve))]

rng = np.random.default_rng(seed_of("mvri_audit", "weights"))
signs = np.array([-1 if AXES[k]["weight"] < 0 else 1 for k in AXES])
peaks = []
for _ in range(10000):
    w = rng.uniform(0.05, 1.0, size=len(AXES)) * signs
    w = w / np.abs(w).sum()
    peaks.append(binmeans.index[int(np.argmax(binmeans.values @ w))])
pk = pd.Series(peaks).value_counts(normalize=True).sort_values(ascending=False)
child = float(pk.get("Young child", 0.0) + pk.get("Child", 0.0))
pd.DataFrame({"age_bin": pk.index, "frac_peak": pk.values}).to_csv(
    OUT / "weight_sensitivity.csv", index=False)

# ---------------------------------------------------------------- E: dimensionality
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from nasage.stats import effective_tests          # noqa: E402
C = Z.corr(method="spearman")
eff = effective_tests(C)
C.to_csv(OUT / "axis_correlation.csv")

print("\n" + "=" * 78)
print("PRE-SPECIFIED VERDICTS (docs/PREREG_mvri_presentation_audit.md section 3)")
print("=" * 78)
print(R[["axis", "frac_measurable", "gate_pass", "KW_p_percell", "KW_p_perdonor",
         "BH_perdonor", "empirical_p", "age_varying", "specific"]].round(4).to_string(index=False))
print(f"\nD  deck weights peak at: {deck_peak}")
print(f"   childhood peaks in {child*100:.1f}% of 10,000 sign-respecting weightings "
      f"-> {'ROBUST' if child >= 0.90 else 'NOT robust'}")
print("   peak distribution:", dict(pk.round(3)))
print(f"\nE  5 axes -> {eff['eff_li_ji']:.1f} effective tests "
      f"(PC1 {eff['pc1_var_frac']*100:.0f}% of variance)")
print("\nwrote", OUT)
