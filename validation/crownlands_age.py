"""Steps 4-5 of docs/PREREG_crownlands_age.md: the age model and its specificity controls.

Runs the package's normal path on the fixed subsample built by crownlands_subsample.py:
UCell per cell -> donor-level means per cell type -> age_effects -> random-set and
permutation controls.

Everything about the design is fixed by the plan, not here:
  * primary family = the 24 module x cell-type combinations passing the gate (step 3);
  * age_bin with 18-39 as reference; two-sided throughout, no directional prior;
  * BH within the primary family; >= 500 expression-matched random sets; 500 permutations.

Usage: python validation/crownlands_age.py <subsample.h5ad> <gate.csv> <out_dir> [n_sets] [n_perm]
"""
import sys
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

import nasage as nas

IN, GATE, OUT = sys.argv[1], sys.argv[2], Path(sys.argv[3])
N_SETS = int(sys.argv[4]) if len(sys.argv) > 4 else 500
N_PERM = int(sys.argv[5]) if len(sys.argv) > 5 else 500
OUT.mkdir(parents=True, exist_ok=True)
AGE_ORDER = ["18-39", "40-64", "65+"]

# ---------------------------------------------------------------- load
f = h5py.File(IN, "r")
shape = tuple(int(x) for x in f["X"].attrs["shape"])
X = sp.csr_matrix((f["X/data"][:], f["X/indices"][:], f["X/indptr"][:]), shape=shape)
sym = np.array([x.decode() if isinstance(x, bytes) else x for x in f["var_symbol"][:]])
obs = pd.DataFrame({c: np.array([x.decode() for x in f[f"obs/{c}"][:]])
                    for c in ("donor", "ct", "age", "age_bin", "sex", "cell_type_raw")})
f.close()

# collapse duplicate symbols so the var index is unique
if pd.Index(sym).has_duplicates:
    keep = ~pd.Index(sym).duplicated()
    X, sym = X[:, keep], sym[keep]
    print(f"dropped {int((~keep).sum())} duplicate gene symbols")

A = ad.AnnData(X=X, obs=obs.reset_index(drop=True),
               var=pd.DataFrame(index=pd.Index(sym, name=None)))
A.obs["age_bin"] = pd.Categorical(A.obs.age_bin, categories=AGE_ORDER, ordered=True)
print(f"AnnData {A.shape}, donors {A.obs.donor.nunique()}")
print(A.obs.groupby("ct", observed=True).donor.nunique().to_string())

# ---------------------------------------------------------------- modules + applicability
gate = pd.read_csv(GATE)
passing = gate[gate.passes]
family = sorted(zip(passing.module, passing.cell_type))
print(f"\nprimary family: {len(family)} module x cell-type combinations")

M = nas.load_modules()
mods = {m: d for m, d in M.items() if m != "M07_mhc_II_presentation"}
# M06 enters only as the HLA-free variant (plan section 3)
hla_free = nas.drop_hla({"M06_mhc_I_processing": M["M06_mhc_I_processing"]})
mods["M06_noHLA"] = hla_free["M06_mhc_I_processing"]
del mods["M06_mhc_I_processing"]

# applicability: 1.0 exactly for the gate-passing combinations, 0.0 elsewhere.
# This makes the tested primary family identical to the family declared in step 3.
base = nas.load_applicability()
cts = ["basal", "secretory", "ciliated"]
appl = {m: {c: (1.0 if (m, c) in set(family) else 0.0) for c in cts} for m in mods}

nas.tl.score_modules(A, modules=mods, n_jobs=1)
print("scored:", len([c for c in A.obs.columns if c.endswith("_UCell")]), "modules")

# ---------------------------------------------------------------- donor summary + age model
nas.tl.donor_summary(
    A, donor_key="donor", celltype_key="ct", age_key="age_bin", age_order=AGE_ORDER,
    sex_key="sex", subtype_key="cell_type_raw", batch_key=None,
    min_cells=20, applicability=appl, primary_celltypes=("basal", "secretory", "ciliated"),
)
eff = nas.tl.age_effects(A, reference="18-39", adjust_depth=True, adjust_subtypes=True)
eff.to_csv(OUT / "age_effects_raw.csv", index=False)
print("\nage_effects columns:", list(eff.columns))
print(eff.head(3).to_string())

# ---------------------------------------------------------------- controls
print(f"\nrandom-set control, {N_SETS} expression-matched sets ...", flush=True)
nullc = nas.tl.null_control(A, n_sets=N_SETS, seed=0, n_bins=25, n_jobs=1, tiers=("primary",))
nullc.to_csv(OUT / "random_set_control.csv", index=False)

print(f"permutation null, {N_PERM} permutations ...", flush=True)
perm = nas.tl.permutation_null(A, n_perm=N_PERM, seed=0, stratify=None)
perm.to_csv(OUT / "permutation_null.csv", index=False)

import json
(OUT / "uns_nasage.json").write_text(json.dumps(
    {k: str(v)[:2000] for k, v in A.uns.get("nasage", {}).items()}, indent=2))
A.obs.to_csv(OUT / "cell_scores_obs.csv.gz", compression="gzip")
print("\nwrote", OUT)
