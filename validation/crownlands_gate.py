"""Step 3 of docs/PREREG_crownlands_age.md: the measurability gate.

Reports, per module and primary cell type, the fraction of the module's available genes that
reach 5% detection. A module is tested in a cell type only if that fraction is >= 60% (plan
section 3.1). The family size and the smallest attainable BH are printed BEFORE any p-value is
computed anywhere, per plan section 8 step 3.

M07 is excluded outright (all HLA class II absent from the probe panel, plan section 3).
M06 enters only as the HLA-free variant.

Usage: python validation/crownlands_gate.py <subsample.h5ad> <out_dir>
"""
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

import nasage as nas

IN, OUT = sys.argv[1], Path(sys.argv[2])
OUT.mkdir(parents=True, exist_ok=True)

GATE_DETECT = 0.05      # a gene counts as measurable at >= 5% of cells
GATE_FRACTION = 0.60    # a module passes if >= 60% of its available genes count
PRIMARY = ("basal", "secretory", "ciliated")
ALIAS = {"CIMAP3": "ODF3B", "STING1": "TMEM173"}

f = h5py.File(IN, "r")
var = np.array([x.decode() if isinstance(x, bytes) else x for x in f["var_symbol"][:]])
g2i = {}
for i, g in enumerate(var):          # first occurrence wins; symbols are unique in this file
    g2i.setdefault(g, i)
ct = np.array([x.decode() for x in f["obs/ct"][:]])
cats, ccode = np.unique(ct, return_inverse=True)
print(f"{len(ct):,} cells, {len(var):,} genes, cell types {dict(zip(cats, np.bincount(ccode)))}")

# ---- module gene sets as they will actually be scored
M = nas.load_modules()
sets = {}
for m, d in M.items():
    if m == "M07_mhc_II_presentation":
        continue                                   # excluded by plan section 3
    genes = [ALIAS.get(g, g) for g in d["genes"]]
    if m == "M06_mhc_I_processing":
        m, genes = "M06_noHLA", [g for g in genes if not g.startswith("HLA-")]
    sets[m] = sorted({g for g in genes if g in g2i})

genes = sorted({g for v in sets.values() for g in v})
gidx = np.array([g2i[g] for g in genes])
lut = np.full(len(var), -1, np.int32)
lut[gidx] = np.arange(len(genes))

# ---- detection scan
X = f["X"]
indptr = X["indptr"][:].astype(np.int64)
n = len(indptr) - 1
det = np.zeros((len(cats), len(genes)), np.int64)
ncell = np.bincount(ccode, minlength=len(cats))

CH = 20000
for s in range(0, n, CH):
    e = min(s + CH, n)
    lo, hi = indptr[s], indptr[e]
    ind = X["indices"][lo:hi]
    rows = np.repeat(np.arange(s, e), np.diff(indptr[s:e + 1]))
    m = lut[ind] >= 0
    if m.any():
        np.add.at(det.reshape(-1), ccode[rows[m]] * len(genes) + lut[ind[m]], 1)

frac = pd.DataFrame(det / ncell[:, None], index=cats, columns=genes)
frac.T.to_csv(OUT / "gene_detection.csv")

# ---- the gate
rows = []
for m, gl in sets.items():
    for c in PRIMARY:
        v = frac.loc[c, gl]
        n_ok = int((v >= GATE_DETECT).sum())
        rows.append(dict(module=m, cell_type=c, n_genes_in_module=len(M.get(
            "M06_mhc_I_processing" if m == "M06_noHLA" else m, {"genes": []})["genes"]),
            n_available=len(gl), n_detected=n_ok,
            frac_detected=round(n_ok / len(gl), 3) if gl else 0.0,
            median_detection=round(float(v.median()), 4) if gl else 0.0,
            passes=bool(gl and n_ok / len(gl) >= GATE_FRACTION)))
G = pd.DataFrame(rows)
G.to_csv(OUT / "measurability_gate.csv", index=False)

pd.set_option("display.width", 200)
print("\n=== measurability gate (>= 60% of available genes at >= 5% detection) ===")
print(G.pivot(index="module", columns="cell_type", values="frac_detected")
      .reindex(columns=list(PRIMARY)).to_string())
print("\npass / fail:")
print(G.pivot(index="module", columns="cell_type", values="passes")
      .reindex(columns=list(PRIMARY)).to_string())

n_tests = int(G.passes.sum())
print(f"\n=== FAMILY SIZE, fixed before any p-value ===")
print(f"module x cell-type combinations passing the gate: {n_tests}")
print(f"modules passing in at least one primary cell type: "
      f"{G[G.passes].module.nunique()} of {len(sets)}")
print(f"smallest attainable BH on the primary family = p_min x {n_tests}")
print(f"  -> BH < 0.05 requires p < {0.05 / n_tests:.2e} for the top-ranked test")
print("\nThe primary p-values come from the parametric age model, which is not floored,")
print("so BH < 0.05 is attainable. The >= 500-set random and permutation controls of plan")
print("section 5.1 have an empirical floor of p = 1/501 = 0.002; they are a separate")
print("specificity gate, not part of this BH family.")

failed = G[~G.passes]
if len(failed):
    print("\n=== not tested (poorly_measurable), reported with detection and excluded "
          "from the family ===")
    print(failed[["module", "cell_type", "n_available", "n_detected", "median_detection"]]
          .to_string(index=False))
