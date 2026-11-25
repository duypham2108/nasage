"""Dated amendment, 2026-11-25: re-run the MVRI audit on the complete cell labels.

Amends `docs/PREREG_mvri_presentation_audit.md`. The locked run
(`validation/mvri_presentation_audit.py`, results in `results/mvri_audit/`) scored cells using
`Cell_type_annotation_level1/2`. Those two columns are NaN for 53,150 of 148,978 cells -- the
whole of 16 donors, including four of the five Elderly. Every axis subset is defined by a cell
label, so all five axes silently dropped the same 16 donors, and the composite's "36 of 52
donors, Elderly n = 1" was a property of the label column rather than of the index.

`Cell_type_prediction_level1/2/3` are populated for every cell, and level3 is the column the
manuscript's own Supplementary Figure 1a is drawn from. This script re-runs the complete audit
(tests A-E of the plan) on those labels and writes the per-panel detail Figure 5 needs.

The locked results are not touched. This writes to a separate directory so the two can be
compared, and the figure legend reports this run with the original as the sensitivity.

Nothing else changes: the same gene lists, the same subsets, the same scorer, the same seeds,
the same verdict rules.

Usage: python validation/mvri_audit_amendment.py <h5ad> <out_dir> [n_null] [annotation|prediction]
"""
import gc
import sys
import zlib
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from mvri_axes import AGE_ORDER, AXES, resolve, score_cells   # noqa: E402
from nasage.stats import effective_tests                      # noqa: E402

H5, OUT = Path(sys.argv[1]), Path(sys.argv[2])
NNULL = int(sys.argv[3]) if len(sys.argv) > 3 else 500
WHICH = sys.argv[4] if len(sys.argv) > 4 else "prediction"
OUT.mkdir(parents=True, exist_ok=True)

COL1 = f"Cell_type_{WHICH}_level1"
COL2 = f"Cell_type_{WHICH}_level2"

PANEL_D_TYPES = ["Basal", "Club", "Secretory", "Goblet", "Goblet inflammatory", "Ciliated",
                 "Macrophage", "Monocyte", "cDC", "pDC", "B mem", "Plasma"]
MIN_CELLS_PANEL_D = 200


def seed_of(*p):
    return zlib.crc32("|".join(map(str, p)).encode()) % (2 ** 32)


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


print("reading", H5.name, "| labels:", COL1, "/", COL2)
A = ad.read_h5ad(H5)
var = pd.Index(A.var_names)
pos = {g: i for i, g in enumerate(var)}
L1 = A.obs[COL1].astype(str).values
L2 = A.obs[COL2].astype(str).values
donor = A.obs["donor"].astype(str).values
agebin = A.obs["Age_bin"].astype(str).values
N_OBS = A.n_obs
n_unlab = int((L1 == "nan").sum())
print(f"  {N_OBS} cells, {len(set(donor))} donors, {n_unlab} cells without a level1 label")

X = A.X if sp.issparse(A.X) else sp.csr_matrix(A.X)
A.X = None
del A
gc.collect()
tot = np.asarray(X.sum(1)).ravel()
tot[tot == 0] = 1
X = (sp.diags(1e4 / tot) @ X).tocsr()
np.log1p(X.data, out=X.data)
gc.collect()

# level1 category names differ between the two label sets; map to the audit's vocabulary
EPI = {"Epi", "Epithelial"}
MYE = {"Myeloid"}


def subset_mask(spec):
    m = np.zeros(N_OBS, bool)
    for s in spec:
        if s == "__EPI__":
            m |= np.isin(L1, list(EPI))
        elif s == "__MYELOID__":
            m |= np.isin(L1, list(MYE))
        else:
            m |= np.char.startswith(L2.astype(str), s)
    return m


rows_out, donor_tbl = [], {}
det_rows, donor_rows, cell_rows, ct_rows, null_rows, win_rows = [], [], [], [], [], []
rank = {b: i for i, b in enumerate(AGE_ORDER)}

for name, ax in AXES.items():
    rows = np.where(subset_mask(ax["subset"]))[0]
    idx = resolve(ax["genes"], pos)
    n_req = len({g for g in ax["genes"] if g != "IL8"})
    Xs = X[rows].tocsc()
    dkey, akey, tkey = donor[rows], agebin[rows], L2[rows]

    det = np.asarray((Xs[:, idx] > 0).sum(0)).ravel() / len(rows)
    frac_meas = float((det >= 0.05).mean()) if len(idx) else 0.0
    gate = frac_meas >= 0.60
    for g, d_ in zip([var[j] for j in idx], det):
        det_rows.append(dict(axis=name, gene=g, detection=float(d_),
                             measurable=bool(d_ >= 0.05)))

    v = score_cells(idx, Xs)
    cell_rows.append(pd.DataFrame({"axis": name, "age": akey, "score": v}))
    s_cell, p_cell = kw(by_bin(v, akey))

    d = (pd.DataFrame({"donor": dkey, "age": akey, "v": v})
         .groupby(["donor", "age"], observed=True)
         .agg(score=("v", "mean"), n_cells=("v", "size")).reset_index())
    s_don, p_don = kw(by_bin(d.score.values, d.age.values))
    donor_tbl[name] = d.set_index("donor")
    donor_rows.append(d.assign(axis=name))

    ct = (pd.DataFrame({"t": tkey, "v": v}).groupby("t", observed=True)
          .agg(mean_score=("v", "mean"), n_cells=("v", "size")).reset_index())
    ct = ct[(ct.n_cells >= MIN_CELLS_PANEL_D) & ct.t.isin(PANEL_D_TYPES)]
    ct_rows.append(ct.assign(axis=name))

    for wname, bins_ in [("Neonate to Young child", AGE_ORDER[:3]),
                         ("Child to Elderly", AGE_ORDER[3:])]:
        sub = d[d.age.isin(bins_)]
        if sub.age.nunique() < 3:
            continue
        r = stats.spearmanr(sub.score.values, [rank[a] for a in sub.age])
        rng_b = np.random.default_rng(seed_of("fig5_window", name, wname))
        boot = []
        for _ in range(2000):
            k = rng_b.integers(0, len(sub), len(sub))
            s2 = sub.iloc[k]
            if s2.age.nunique() < 2:
                continue
            b_ = stats.spearmanr(s2.score.values, [rank[a] for a in s2.age]).statistic
            if np.isfinite(b_):
                boot.append(b_)
        boot = np.asarray(boot)
        win_rows.append(dict(axis=name, window=wname, n_donors=len(sub),
                             rho=float(r.statistic), p=float(r.pvalue),
                             lo=float(np.percentile(boot, 2.5)),
                             hi=float(np.percentile(boot, 97.5))))

    mean_all = np.asarray(Xs.sum(0)).ravel() / len(rows)
    qbins = pd.qcut(pd.Series(mean_all).rank(method="first"), 25, labels=False)
    pool = {k: np.array(g) for k, g in qbins.groupby(qbins).groups.items()}
    gb = [qbins.iloc[i] for i in idx]
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
            nstat.append(sn)
            npv.append(pn)
    nstat, npv = np.asarray(nstat), np.asarray(npv)
    emp = float((nstat >= s_don).mean()) if nstat.size else np.nan
    null_rows.append(pd.DataFrame({"axis": name, "draw": np.arange(len(nstat)),
                                   "kw_stat": nstat}))

    rows_out.append(dict(axis=name, n_genes_listed=n_req, n_genes_found=len(idx),
                         n_cells=len(rows), n_donors=d.donor.nunique(),
                         frac_measurable=frac_meas, gate_pass=bool(gate),
                         KW_p_percell=p_cell, KW_p_perdonor=p_don,
                         KW_stat_perdonor=s_don, empirical_p=emp,
                         null_frac_p05=float((npv < 0.05).mean()) if npv.size else np.nan))
    r = rows_out[-1]
    print(f"{name:22s} genes {len(idx):2d}/{n_req:2d} meas {frac_meas:.2f} "
          f"{'PASS' if gate else 'FAIL'} | donors {r['n_donors']:2d} | "
          f"p_cell {p_cell:.2e} p_donor {p_don:.4f} | emp_p {emp:.3f} "
          f"null05 {r['null_frac_p05']:.3f}")
    del Xs
    gc.collect()

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

S = pd.DataFrame({k: v.groupby("donor").score.mean() for k, v in donor_tbl.items()})
age_of = (pd.DataFrame({"donor": donor, "age": agebin})
          .drop_duplicates().set_index("donor").age)
S["age"] = age_of.reindex(S.index).values
S = S.dropna()
Z = S[list(AXES)].apply(lambda c: (c - c.mean()) / c.std(ddof=1))
binmeans = Z.assign(age=S.age.values).groupby("age").mean()
binmeans = binmeans.reindex([b for b in AGE_ORDER if b in binmeans.index])
binmeans.to_csv(OUT / "axis_means_by_agebin.csv")
S.to_csv(OUT / "donor_axis_scores.csv")

w_deck = np.array([AXES[k]["weight"] for k in AXES])
deck_peak = binmeans.index[int(np.argmax(binmeans.values @ w_deck))]
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

C = Z.corr(method="spearman")
C.to_csv(OUT / "axis_correlation.csv")
eff = effective_tests(C)

F5 = OUT / "figure5"
F5.mkdir(exist_ok=True)
pd.DataFrame(det_rows).to_csv(F5 / "per_gene_detection.csv", index=False)
pd.concat(donor_rows).to_csv(F5 / "donor_axis_long.csv", index=False)
pd.concat(cell_rows).to_csv(F5 / "cell_scores.csv.gz", index=False)
pd.concat(ct_rows).to_csv(F5 / "celltype_axis_means.csv", index=False)
pd.concat(null_rows).to_csv(F5 / "null_stats.csv.gz", index=False)
pd.DataFrame(win_rows).to_csv(F5 / "window_trends.csv", index=False)

print("\n" + "=" * 78)
print(f"AMENDED RUN on {WHICH} labels ({n_unlab} unlabelled cells)")
print("=" * 78)
print(R[["axis", "n_donors", "frac_measurable", "gate_pass", "KW_p_perdonor",
         "BH_perdonor", "empirical_p", "age_varying", "specific"]].round(4).to_string(index=False))
print(f"\nD  deck weights peak at: {deck_peak}")
print(f"   childhood peaks in {child*100:.1f}% of 10,000 weightings "
      f"-> {'ROBUST' if child >= 0.90 else 'NOT robust'}")
print("   peak distribution:", dict(pk.round(3)))
print(f"\nE  5 axes -> {eff['eff_li_ji']:.1f} effective tests "
      f"(PC1 {eff['pc1_var_frac']*100:.0f}% of variance)")
print(f"\ncomposite computable for {len(S)} donors")
print("wrote", OUT)
