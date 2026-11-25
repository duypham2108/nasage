"""Everything manuscript Figure 5 and Supplementary Figure 5 need, in one pass over the object.

Plan: `docs/FIGURE5_build_plan.md`. This computes no new test -- the verdicts stay those of
`docs/PREREG_mvri_presentation_audit.md`, already written to `axis_audit.csv`. What it adds is
the per-panel detail the audit summarised away: per-gene detection, per-cell scores, the
donor table for every donor carrying an axis's subset (the audit saved only the 36 donors with
*all five*), per-cell-type means, the 500 matched-null statistics, and the two-window trends.

The null loop reproduces the audit's draws exactly -- same `seed_of`, same rng sequence, same
scorer -- and asserts the recomputed empirical p matches `axis_audit.csv` before writing
anything, so a drift in either copy of the axis definitions is caught here rather than in a
figure.

Usage: python validation/figure5_panel_data.py <h5ad> <audit_dir> [n_null]
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
from mvri_axes import (AGE_ORDER, AXES, check_against_audit,  # noqa: E402
                       resolve, score_cells)

H5, AUD = Path(sys.argv[1]), Path(sys.argv[2])
NNULL = int(sys.argv[3]) if len(sys.argv) > 3 else 500
OUT = AUD / "figure5"
OUT.mkdir(parents=True, exist_ok=True)

AUDIT = check_against_audit(AUD / "axis_audit.csv")

# cell types shown in panel d, in epithelium-then-immune order; others pooled out
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


print("reading", H5.name)
A = ad.read_h5ad(H5)
print(f"  {A.n_obs} cells x {A.n_vars} genes, {A.obs.donor.nunique()} donors")

var = pd.Index(A.var_names)
pos = {g: i for i, g in enumerate(var)}
L1 = A.obs["Cell_type_annotation_level1"].astype(str).values
L2 = A.obs["Cell_type_annotation_level2"].astype(str).values
donor = A.obs["donor"].astype(str).values
agebin = A.obs["Age_bin"].astype(str).values
N_OBS = A.n_obs

X = A.X if sp.issparse(A.X) else sp.csr_matrix(A.X)
# Drop the AnnData before normalising. The raw matrix is ~2 GB on its own, and holding it
# alongside the normalised copy and a per-axis CSC slice is what exhausted memory here.
A.X = None
del A
gc.collect()

tot = np.asarray(X.sum(1)).ravel()
tot[tot == 0] = 1
X = (sp.diags(1e4 / tot) @ X).tocsr()
np.log1p(X.data, out=X.data)      # in place: a copy of the 278M-element data array is 2.1 GB
gc.collect()


def subset_mask(spec):
    m = np.zeros(N_OBS, bool)
    for s in spec:
        if s == "__EPI__":
            m |= (L1 == "Epi")
        elif s == "__MYELOID__":
            m |= (L1 == "Myeloid")
        else:
            m |= np.char.startswith(L2.astype(str), s)
    return m


det_rows, donor_rows, cell_rows, ct_rows, null_rows, win_rows = [], [], [], [], [], []

for name, ax in AXES.items():
    rows = np.where(subset_mask(ax["subset"]))[0]
    idx = resolve(ax["genes"], pos)
    Xs = X[rows].tocsc()
    dkey, akey, tkey = donor[rows], agebin[rows], L2[rows]

    # ---- per-gene detection, for the measurability panel
    det = np.asarray((Xs[:, idx] > 0).sum(0)).ravel() / len(rows)
    for g, d_ in zip([var[j] for j in idx], det):
        det_rows.append(dict(axis=name, gene=g, detection=float(d_),
                             measurable=bool(d_ >= 0.05)))

    # ---- scores
    v = score_cells(idx, Xs)
    cell_rows.append(pd.DataFrame({"axis": name, "age": akey, "score": v}))

    d = (pd.DataFrame({"donor": dkey, "age": akey, "v": v})
         .groupby(["donor", "age"], observed=True)
         .agg(score=("v", "mean"), n_cells=("v", "size")).reset_index())
    d["axis"] = name
    donor_rows.append(d)

    # ---- per-cell-type mean (panel d)
    ct = (pd.DataFrame({"t": tkey, "v": v}).groupby("t", observed=True)
          .agg(mean_score=("v", "mean"), n_cells=("v", "size")).reset_index())
    ct = ct[(ct.n_cells >= MIN_CELLS_PANEL_D) & ct.t.isin(PANEL_D_TYPES)]
    ct["axis"] = name
    ct_rows.append(ct)

    # ---- two age windows, Spearman of donor score against bin rank
    rank = {b: i for i, b in enumerate(AGE_ORDER)}
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
            boot.append(stats.spearmanr(s2.score.values,
                                        [rank[a] for a in s2.age]).statistic)
        boot = np.asarray([b for b in boot if np.isfinite(b)])
        win_rows.append(dict(axis=name, window=wname, n_donors=len(sub),
                             rho=float(r.statistic), p=float(r.pvalue),
                             lo=float(np.percentile(boot, 2.5)),
                             hi=float(np.percentile(boot, 97.5))))

    # ---- matched nulls, reproducing the audit's draws exactly
    s_don, _ = kw(by_bin(d.score.values, d.age.values))
    mean_all = np.asarray(Xs.sum(0)).ravel() / len(rows)
    qbins = pd.qcut(pd.Series(mean_all).rank(method="first"), 25, labels=False)
    pool = {k: np.array(g) for k, g in qbins.groupby(qbins).groups.items()}
    gb = [qbins.iloc[i] for i in idx]
    rng = np.random.default_rng(seed_of("mvri_audit", name))
    nstat = []
    for _ in range(NNULL):
        pick = [int(rng.choice(pool[k])) for k in gb if pd.notna(k)]
        if len(pick) < 3:
            continue
        vn = score_cells(pick, Xs)
        dn = (pd.DataFrame({"donor": dkey, "age": akey, "v": vn})
              .groupby(["donor", "age"], observed=True).v.mean().reset_index())
        sn, _ = kw(by_bin(dn.v.values, dn.age.values))
        if np.isfinite(sn):
            nstat.append(sn)
    nstat = np.asarray(nstat)
    emp = float((nstat >= s_don).mean())
    null_rows.append(pd.DataFrame({"axis": name, "draw": np.arange(len(nstat)),
                                   "kw_stat": nstat}))

    want = float(AUDIT.loc[name, "empirical_p"])
    print(f"{name:22s} obs KW {s_don:7.3f} | empirical p {emp:.3f} "
          f"(audit {want:.3f}) {'OK' if NNULL == 500 and abs(emp - want) < 1e-9 else ''}")
    if NNULL == 500:
        assert abs(emp - want) < 1e-9, (
            f"{name}: empirical p {emp} does not reproduce the audit's {want}; "
            "the axis definitions or the scorer have drifted")
    del Xs
    gc.collect()

pd.DataFrame(det_rows).to_csv(OUT / "per_gene_detection.csv", index=False)
pd.concat(donor_rows).to_csv(OUT / "donor_axis_long.csv", index=False)
pd.concat(cell_rows).to_csv(OUT / "cell_scores.csv.gz", index=False)
pd.concat(ct_rows).to_csv(OUT / "celltype_axis_means.csv", index=False)
pd.concat(null_rows).to_csv(OUT / "null_stats.csv.gz", index=False)
pd.DataFrame(win_rows).to_csv(OUT / "window_trends.csv", index=False)
print("\nwrote", OUT)
