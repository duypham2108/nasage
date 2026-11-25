"""Per-gene age trends inside each surviving axis, epithelium only.

Companion to `validation/mvri_audit_epithelial.py`, for panel e of Figure 5: it answers
whether an axis's age trend is carried by the whole gene set or by one or two genes. No new
test of the axis is involved -- these are descriptive per-gene trends, and they are not
corrected or used for any verdict.

For every gene that survives the detection floor in a gate-passing axis: donor-level mean
expression within that axis's epithelial subset, then Spearman rho against age-bin rank over
all seven bins and within the two age windows.

Usage: python validation/figure5_gene_trends.py <h5ad> <audit_dir>
"""
import gc
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from mvri_axes import AGE_ORDER, AXES, resolve      # noqa: E402

H5, AUD = Path(sys.argv[1]), Path(sys.argv[2])
F5 = AUD / "figure5"
COL1, COL2 = "Cell_type_prediction_level1", "Cell_type_prediction_level2"
EPI_L1 = "Epi"

R = pd.read_csv(AUD / "axis_audit.csv").set_index("axis")
KEEP = [a for a in AXES if a in R.index and bool(R.loc[a, "gate_pass"])]
DET = pd.read_csv(F5 / "per_gene_detection.csv")
print("axes:", KEEP)

A = ad.read_h5ad(H5)
var = pd.Index(A.var_names)
pos = {g: i for i, g in enumerate(var)}
L1 = A.obs[COL1].astype(str).values
L2 = A.obs[COL2].astype(str).values
donor = A.obs["donor"].astype(str).values
agebin = A.obs["Age_bin"].astype(str).values
N_OBS = A.n_obs
X = A.X if sp.issparse(A.X) else sp.csr_matrix(A.X)
A.X = None
del A
gc.collect()
tot = np.asarray(X.sum(1)).ravel()
tot[tot == 0] = 1
X = (sp.diags(1e4 / tot) @ X).tocsr()
np.log1p(X.data, out=X.data)
gc.collect()

rank = {b: i for i, b in enumerate(AGE_ORDER)}
rows = []
for name in KEEP:
    spec = AXES[name]["subset"]
    m = np.zeros(N_OBS, bool)
    for s in spec:
        if s in ("__EPI__", "__MYELOID__"):
            m |= (L1 == EPI_L1)
        else:
            m |= np.char.startswith(L2.astype(str), s)
    m &= (L1 == EPI_L1)
    ridx = np.where(m)[0]
    keep_g = set(DET[(DET.axis == name) & DET.measurable].gene)
    idx = [j for j in resolve(AXES[name]["genes"], pos) if var[j] in keep_g]
    Xs = X[ridx].tocsc()
    M = np.asarray(Xs[:, idx].todense())
    dk, ak = donor[ridx], agebin[ridx]
    df = pd.DataFrame(M, columns=[var[j] for j in idx])
    df["donor"], df["age"] = dk, ak
    dm = df.groupby(["donor", "age"], observed=True).mean().reset_index()
    ages = np.array([rank[a] for a in dm.age])
    for g in [var[j] for j in idx]:
        v = dm[g].values
        r_all = stats.spearmanr(v, ages)
        sel1 = np.isin(dm.age, AGE_ORDER[:3])
        sel2 = np.isin(dm.age, AGE_ORDER[3:])
        r1 = stats.spearmanr(v[sel1], ages[sel1])
        r2 = stats.spearmanr(v[sel2], ages[sel2])
        det = float(DET[(DET.axis == name) & (DET.gene == g)].detection.iloc[0])
        rows.append(dict(axis=name, gene=g, detection=det,
                         rho_all=float(r_all.statistic), p_all=float(r_all.pvalue),
                         rho_early=float(r1.statistic), p_early=float(r1.pvalue),
                         rho_late=float(r2.statistic), p_late=float(r2.pvalue)))
        print(f"{name:22s} {g:10s} det {det:.2f}  rho_all {r_all.statistic:+.2f} "
              f"(P {r_all.pvalue:.3f})  early {r1.statistic:+.2f}  late {r2.statistic:+.2f}")
    del Xs, M, df
    gc.collect()

pd.DataFrame(rows).to_csv(F5 / "gene_trends.csv", index=False)
print("\nwrote", F5 / "gene_trends.csv")
