"""Dated amendment, 2026-11-25: the axes restricted to epithelium, undetected genes removed.

Second amendment to `docs/PREREG_mvri_presentation_audit.md`, after
`validation/mvri_audit_amendment.py`. Three declared deviations from the locked plan, all in
the same direction -- scoring less, and scoring it in one compartment:

1. **Epithelium only.** Every axis is scored in the epithelial members of its own declared
   subset. The antigen-presentation axis loses DC/Mac/Mono and the IgA axis loses Plasma and
   B mem. This matches the manuscript's primary family (basal / secretory / ciliated) and
   removes the composition artefact in which an axis's donor score moves with how many
   myeloid cells that donor happened to yield.
2. **Undetected genes are dropped from the score, not merely reported.** A gene detected in
   under 5% of the cells its axis is scored in contributes noise to a mean-of-z, so the score
   is recomputed on the surviving genes. The gate is still evaluated on the full declared
   list, so the measurability verdict is unchanged by this.
3. **Axes that fail the gate are reported but not carried forward.** All five are tested and
   BH-corrected across five, so the multiple-testing correction is unaffected; only the
   surviving axes reach the figure.

Matched nulls are drawn at the *reduced* set size and matched on the expression decile of the
reduced set, so the null answers the question actually being asked.

`scope` selects deviation 1: "all" keeps each axis's declared subset as written (the state
the figure was in before the epithelial restriction), "epi" intersects it with epithelium.
Deviations 2 and 3 apply either way.

Usage: python validation/mvri_audit_epithelial.py <h5ad> <out_dir> [n_null] [all|epi]
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
SCOPE = sys.argv[4] if len(sys.argv) > 4 else "epi"
OUT.mkdir(parents=True, exist_ok=True)
F5 = OUT / "figure5"
F5.mkdir(exist_ok=True)

COL1, COL2 = "Cell_type_prediction_level1", "Cell_type_prediction_level2"
EPI_L1 = "Epi"
GATE, FLOOR = 0.60, 0.05
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


def by_bin(v, k):
    return [v[k == b] for b in AGE_ORDER if (k == b).sum() > 0]


print("reading", H5.name)
A = ad.read_h5ad(H5)
var = pd.Index(A.var_names)
pos = {g: i for i, g in enumerate(var)}
L1 = A.obs[COL1].astype(str).values
L2 = A.obs[COL2].astype(str).values
donor = A.obs["donor"].astype(str).values
agebin = A.obs["Age_bin"].astype(str).values
N_OBS = A.n_obs
print(f"  {N_OBS} cells, {len(set(donor))} donors, "
      f"{int((L1 == EPI_L1).sum())} epithelial cells | scope = {SCOPE}")

X = A.X if sp.issparse(A.X) else sp.csr_matrix(A.X)
A.X = None
del A
gc.collect()
tot = np.asarray(X.sum(1)).ravel()
tot[tot == 0] = 1
X = (sp.diags(1e4 / tot) @ X).tocsr()
np.log1p(X.data, out=X.data)
gc.collect()


def subset_mask(spec):
    """The declared subset, intersected with epithelium when SCOPE == 'epi'."""
    m = np.zeros(N_OBS, bool)
    for s in spec:
        if s == "__EPI__":
            m |= (L1 == EPI_L1)
        elif s == "__MYELOID__":
            m |= (L1 == EPI_L1) if SCOPE == "epi" else (L1 == "Myeloid")
        else:
            m |= np.char.startswith(L2.astype(str), s)
    return (m & (L1 == EPI_L1)) if SCOPE == "epi" else m


rows_out, donor_tbl = [], {}
det_rows, donor_rows, cell_rows, ct_rows, null_rows, win_rows, mem_rows = (
    [], [], [], [], [], [], [])
rank = {b: i for i, b in enumerate(AGE_ORDER)}

for name, ax in AXES.items():
    rows = np.where(subset_mask(ax["subset"]))[0]
    idx_all = resolve(ax["genes"], pos)
    n_req = len({g for g in ax["genes"] if g != "IL8"})
    Xs = X[rows].tocsc()
    dkey, akey, tkey = donor[rows], agebin[rows], L2[rows]

    det = np.asarray((Xs[:, idx_all] > 0).sum(0)).ravel() / len(rows)
    frac_meas = float((det >= FLOOR).mean()) if len(idx_all) else 0.0
    gate = frac_meas >= GATE
    for j, g, d_ in zip(idx_all, [var[j] for j in idx_all], det):
        det_rows.append(dict(axis=name, gene=g, detection=float(d_),
                             measurable=bool(d_ >= FLOOR)))

    # deviation 2: score only the genes that clear the floor
    idx = [j for j, d_ in zip(idx_all, det) if d_ >= FLOOR]
    if len(idx) < 3:
        print(f"{name:22s} only {len(idx)} genes above the floor -- not scored")
        del Xs
        gc.collect()
        continue
    for g in [var[j] for j in idx]:
        mem_rows.append(dict(axis=name, gene=g))

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
    ct = ct[(ct.n_cells >= MIN_CELLS_PANEL_D) & (ct.t != "nan")]
    ct_rows.append(ct.assign(axis=name))

    for wname, bins_ in [("Neonate to Young child", AGE_ORDER[:3]),
                         ("Child to Elderly", AGE_ORDER[3:])]:
        sub = d[d.age.isin(bins_)]
        if sub.age.nunique() < 3:
            continue
        r = stats.spearmanr(sub.score.values, [rank[a] for a in sub.age])
        rb = np.random.default_rng(seed_of("fig5_window", name, wname))
        boot = []
        for _ in range(2000):
            s2 = sub.iloc[rb.integers(0, len(sub), len(sub))]
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
    qb = pd.qcut(pd.Series(mean_all).rank(method="first"), 25, labels=False)
    pool = {k: np.array(g) for k, g in qb.groupby(qb).groups.items()}
    gb = [qb.iloc[i] for i in idx]
    rng = np.random.default_rng(seed_of("mvri_epi", name))
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

    rows_out.append(dict(axis=name, n_genes_listed=n_req, n_genes_found=len(idx_all),
                         n_genes_scored=len(idx), n_cells=len(rows),
                         n_donors=d.donor.nunique(), frac_measurable=frac_meas,
                         gate_pass=bool(gate), KW_p_percell=p_cell, KW_p_perdonor=p_don,
                         KW_stat_perdonor=s_don, empirical_p=emp,
                         null_frac_p05=float((npv < 0.05).mean()) if npv.size else np.nan))
    r = rows_out[-1]
    print(f"{name:22s} genes {len(idx):2d}/{len(idx_all):2d} scored  meas {frac_meas:.2f} "
          f"{'PASS' if gate else 'FAIL'} | cells {len(rows):6d} donors {r['n_donors']:2d} "
          f"| p_donor {p_don:.4f} | emp_p {emp:.3f} null05 {r['null_frac_p05']:.3f}")
    del Xs
    gc.collect()

R = pd.DataFrame(rows_out)
ok = R.KW_p_perdonor.notna()
pv = R.loc[ok, "KW_p_perdonor"].values
o = np.argsort(pv)
q = np.empty_like(pv)
q[o] = np.minimum.accumulate((pv[o] * len(pv) / np.arange(1, len(pv) + 1))[::-1])[::-1]
R.loc[ok, "BH_perdonor"] = np.minimum(q, 1.0)
R["age_varying"] = R.BH_perdonor < 0.05
R["specific"] = R.empirical_p < 0.05
R["carried_forward"] = R.gate_pass & R.specific
R.to_csv(OUT / "axis_audit.csv", index=False)

S = pd.DataFrame({k: v.groupby("donor").score.mean() for k, v in donor_tbl.items()})
age_of = (pd.DataFrame({"donor": donor, "age": agebin})
          .drop_duplicates().set_index("donor").age)
S["age"] = age_of.reindex(S.index).values
S = S.dropna()
S.to_csv(OUT / "donor_axis_scores.csv")
cols = [c for c in S.columns if c != "age"]
Z = S[cols].apply(lambda c: (c - c.mean()) / c.std(ddof=1))
C = Z.corr(method="spearman")
C.to_csv(OUT / "axis_correlation.csv")
eff = effective_tests(C)

pd.DataFrame(det_rows).to_csv(F5 / "per_gene_detection.csv", index=False)
pd.DataFrame(mem_rows).to_csv(F5 / "gene_membership.csv", index=False)
pd.concat(donor_rows).to_csv(F5 / "donor_axis_long.csv", index=False)
pd.concat(cell_rows).to_csv(F5 / "cell_scores.csv.gz", index=False)
pd.concat(ct_rows).to_csv(F5 / "celltype_axis_means.csv", index=False)
pd.concat(null_rows).to_csv(F5 / "null_stats.csv.gz", index=False)
pd.DataFrame(win_rows).to_csv(F5 / "window_trends.csv", index=False)

print("\n" + "=" * 78)
print(f"scope = {SCOPE}, undetected genes dropped from the score")
print("=" * 78)
print(R[["axis", "n_genes_scored", "n_cells", "n_donors", "frac_measurable", "gate_pass",
         "KW_p_perdonor", "BH_perdonor", "empirical_p",
         "carried_forward"]].round(4).to_string(index=False))
print(f"\n{len(cols)} axes -> {eff['eff_li_ji']:.1f} effective tests "
      f"(PC1 {eff['pc1_var_frac']*100:.0f}%)")
print("carried into Figure 5:", list(R.loc[R.carried_forward, "axis"]))
print("wrote", OUT)
