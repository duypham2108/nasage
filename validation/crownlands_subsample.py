"""Step 2 of docs/PREREG_crownlands_age.md: build and save the fixed cell subsample.

Single sequential pass over the 27.65 GB CSR matrix, keeping only the pre-selected rows and
writing them incrementally to a small h5ad. Nothing here computes a score or a p-value.

Fixed by the plan, not by this script:
  * healthy donors only (disease == "normal");
  * primary cell types basal / secretory / ciliated, per the plan's mapping;
  * at most CAP cells per donor x cell type, drawn without replacement with rng(SEED);
  * a donor contributes to a cell type only with >= 20 cells of it.

Usage: python validation/crownlands_subsample.py <in.h5ad> <out.h5ad> [cap] [seed]
"""
import sys
import h5py
import numpy as np
import pandas as pd

IN, OUT = sys.argv[1], sys.argv[2]
CAP = int(sys.argv[3]) if len(sys.argv) > 3 else 500
SEED = int(sys.argv[4]) if len(sys.argv) > 4 else 0
MIN_CELLS = 20

# Author annotation (obs["fine_celltype"]), not the CELLxGENE ontology label.
# Deviation recorded 2026-11-24: the ontology mapping sent "mucus secreting cell" to
# secretory, which is the submucosal Mucous_gland population, and left Goblet, Club and
# Serous out. That is not the population the primary cohort calls Secretory. Basal and
# ciliated are unchanged by the switch (exact same cells).
CELLTYPE_KEY = "fine_celltype"
MAPPING = {
    "Quiescent_HBC": "basal",
    "Cycling_HBC": "basal",
    "Club": "secretory",
    "Serous": "secretory",
    "Proliferating_Secretory": "secretory",
    "Multiciliated": "ciliated",
    # exploratory, carried but not in the primary family
    "Goblet": "goblet",
    "Mucous_gland": "mucous_gland",
    "Suprabasal": "suprabasal",
}
PRIMARY = ("basal", "secretory", "ciliated")
BATCH_COLS = ("flex_version", "device_guided", "cohort", "sample_id")

f = h5py.File(IN, "r")


def col(name):
    g = f["obs"][name]
    if isinstance(g, h5py.Group):
        cats = np.array([x.decode() if isinstance(x, bytes) else x
                         for x in g["categories"][:]])
        return cats[g["codes"][:]]
    v = g[:]
    return np.array([x.decode() if isinstance(x, bytes) else x for x in v])


donor = col("donor_id")
disease = col("disease")
ctype = col(CELLTYPE_KEY)
stage = col("development_stage")
sex = col("sex")

obs = pd.DataFrame({"donor": donor, "disease": disease, "cell_type_raw": ctype,
                    "stage": stage, "sex": sex})
for b in BATCH_COLS:
    obs[b] = col(b)
# one batch factor for the model: assay version x guided device
obs["batch"] = obs.flex_version.astype(str) + "_" + obs.device_guided.astype(str)
obs["ct"] = pd.Series(ctype).map(MAPPING).values
obs["age"] = pd.Series(stage).str.extract(r"(\d+)-year-old")[0].astype(float).values

elig = (obs.disease == "normal") & obs.ct.notna()
print(f"total cells {len(obs):,}; healthy {int((obs.disease=='normal').sum()):,}; "
      f"healthy in mapped cell types {int(elig.sum()):,}")
print("batch levels:", obs.loc[elig, 'batch'].value_counts().to_dict())
print(f"donors: all {obs.donor.nunique()}, healthy {obs.loc[obs.disease=='normal','donor'].nunique()}")

rng = np.random.default_rng(SEED)
pick, dropped = [], []
for (d, c), idx in obs[elig].groupby(["donor", "ct"], observed=True).groups.items():
    idx = np.asarray(idx)
    if len(idx) < MIN_CELLS:
        dropped.append((d, c, len(idx)))
        continue
    pick.append(idx if len(idx) <= CAP else rng.choice(idx, CAP, replace=False))
sel = np.sort(np.concatenate(pick))
print(f"\nselected {len(sel):,} cells (cap {CAP}, seed {SEED})")
print(f"donor x cell-type groups dropped for < {MIN_CELLS} cells: {len(dropped)}")

sub = obs.iloc[sel].copy()
bins = pd.cut(sub.age, [17, 39, 64, 200], labels=["18-39", "40-64", "65+"])
sub["age_bin"] = bins.values
print("\ncells per cell type:\n", sub.ct.value_counts().to_string())
print("\ndonors per cell type:\n", sub.groupby("ct").donor.nunique().to_string())
print("\ndonors per cell type x age bin:")
print(sub.drop_duplicates(["donor", "ct"]).pivot_table(
    index="ct", columns="age_bin", values="donor", aggfunc="count", observed=False).to_string())

# ---------------------------------------------------------------- extract rows
X = f["X"]
indptr = X["indptr"][:].astype(np.int64)
n_genes = int(X.attrs["shape"][1])
keep = np.zeros(len(indptr) - 1, bool)
keep[sel] = True
new_pos = np.full(len(keep), -1, np.int64)
new_pos[sel] = np.arange(len(sel))

import os

# Resume: if X is already complete from an earlier run, only the obs block is rewritten.
have_X = False
if os.path.exists(OUT):
    try:
        with h5py.File(OUT, "r") as chk:
            have_X = ("X/indptr" in chk
                      and tuple(chk["X"].attrs["shape"]) == (len(sel), n_genes))
    except OSError:
        have_X = False

if have_X:
    print("\nX already written and matches the selection; rewriting obs only")
    out = h5py.File(OUT, "a")
    written = out["X/data"].shape[0]
    for c in list(out["obs"].keys()) if "obs" in out else []:
        del out[f"obs/{c}"]
else:
    out = h5py.File(OUT, "w")
    gd = out.create_dataset("X/data", (0,), maxshape=(None,), dtype="float32",
                            chunks=(1 << 20,), compression="lzf")
    gi = out.create_dataset("X/indices", (0,), maxshape=(None,), dtype="int32",
                            chunks=(1 << 20,), compression="lzf")
    rowlen = np.zeros(len(sel), np.int64)

    CH = 50000
    written = 0
    for s in range(0, len(keep), CH):
        e = min(s + CH, len(keep))
        m = keep[s:e]
        if not m.any():
            continue
        lo, hi = indptr[s], indptr[e]
        dat = X["data"][lo:hi]
        ind = X["indices"][lo:hi]
        lens = np.diff(indptr[s:e + 1])
        starts = np.concatenate([[0], np.cumsum(lens)])
        rows = np.flatnonzero(m)
        take = np.concatenate([np.arange(starts[r], starts[r + 1]) for r in rows])
        gd.resize((written + len(take),)); gd[written:] = dat[take].astype("float32")
        gi.resize((written + len(take),)); gi[written:] = ind[take].astype("int32")
        rowlen[new_pos[s + rows]] = lens[rows]
        written += len(take)
        print(f"  {e:,}/{len(keep):,}  nnz {written:,}", flush=True)

    out.create_dataset("X/indptr",
                       data=np.concatenate([[0], np.cumsum(rowlen)]).astype("int64"))
    out["X"].attrs["encoding-type"] = "csr_matrix"
    out["X"].attrs["encoding-version"] = "0.1.0"
    out["X"].attrs["shape"] = np.array([len(sel), n_genes], "int64")

out.attrs["cap"], out.attrs["seed"] = CAP, SEED
S = h5py.string_dtype()
for k in ("var_index", "var_symbol"):
    if k in out:
        del out[k]
out.create_dataset("var_index", data=f["var"]["_index"][:])
# symbols live in var/feature_name (categorical); the index is Ensembl
fn = f["var"]["feature_name"]
symcats = np.array([x.decode() if isinstance(x, bytes) else x for x in fn["categories"][:]])
sym = np.asarray(symcats[fn["codes"][:]], dtype=object)
out.create_dataset("var_symbol", shape=sym.shape, dtype=S, data=sym)
for c in ("donor", "ct", "age", "age_bin", "sex", "cell_type_raw", "batch",
          "flex_version", "device_guided", "cohort", "sample_id"):
    v = np.asarray(sub[c].astype(str).values, dtype=object)
    out.create_dataset(f"obs/{c}", shape=v.shape, dtype=S, data=v)
out.close()
print(f"\nwrote {OUT}: {len(sel):,} cells x {n_genes:,} genes, nnz {written:,}")
sub.to_csv(OUT.replace(".h5ad", "_obs.csv"), index=False)
