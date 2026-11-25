"""Measurability check of the curated modules (esp. M12) in the user's own nasal AnnData.

Chunked CSR scan: per cell type, fraction of cells with >0 counts for each module gene,
and mean counts-per-10k. No full load.
"""
import sys, h5py, numpy as np, pandas as pd

F = sys.argv[1]
OUT = sys.argv[2]
MODS = r"D:\research\nasage_validation\nasage_module_genes_long.csv"
ALIAS = {"CIMAP3": "ODF3B", "STING1": "TMEM173"}

mods = pd.read_csv(MODS)
mods["gene"] = mods["gene"].replace(ALIAS)

f = h5py.File(F, "r")
var = np.array([x.decode() if isinstance(x, bytes) else x for x in f["var"]["_index"][:]])
g2i = {g: i for i, g in enumerate(var)}

genes = sorted(set(mods["gene"]) & set(g2i))
gidx = np.array([g2i[g] for g in genes])
lut = np.full(len(var), -1, dtype=np.int32)
lut[gidx] = np.arange(len(gidx))

o = f["obs"]
def col(n):
    g = o[n]
    if isinstance(g, h5py.Group):
        c = np.array([x.decode() if isinstance(x, bytes) else x for x in g["categories"][:]])
        return c[g["codes"][:]]
    v = g[:]
    return np.array([x.decode() if isinstance(x, bytes) else x for x in v]) if v.dtype.kind in "SO" else v

ct = col("Cell_type_annotation_level2")
cats, ccode = np.unique(ct, return_inverse=True)
nct = len(cats)

X = f["X"]
indptr = X["indptr"][:].astype(np.int64)
n_cells = len(indptr) - 1

det = np.zeros((nct, len(genes)), dtype=np.int64)     # cells with >0
tot = np.zeros((nct, len(genes)), dtype=np.float64)   # summed counts
depth = np.zeros(nct, dtype=np.float64)               # summed total counts (module-gene universe excluded; full)
ncell = np.bincount(ccode, minlength=nct).astype(np.int64)

CH = 20000  # cells per chunk
for s in range(0, n_cells, CH):
    e = min(s + CH, n_cells)
    lo, hi = indptr[s], indptr[e]
    ind = X["indices"][lo:hi]
    dat = X["data"][lo:hi]
    rows = np.repeat(np.arange(s, e), np.diff(indptr[s:e + 1]))
    np.add.at(depth, ccode[rows], dat)
    m = lut[ind] >= 0
    if m.any():
        gg = lut[ind[m]]
        cc = ccode[rows[m]]
        flat = cc * len(genes) + gg
        np.add.at(det.reshape(-1), flat, 1)
        np.add.at(tot.reshape(-1), flat, dat[m])
    print(f"  {e}/{n_cells}", flush=True)

frac = det / ncell[:, None]
cp10k = tot / depth[:, None] * 1e4

rec = []
gm = mods.drop_duplicates("gene").set_index("gene")["module"]
for i, c in enumerate(cats):
    for j, g in enumerate(genes):
        rec.append((c, int(ncell[i]), gm[g], g, frac[i, j], cp10k[i, j]))
df = pd.DataFrame(rec, columns=["cell_type", "n_cells", "module", "gene", "frac_detected", "cp10k"])
df.to_csv(OUT, index=False)
print("wrote", OUT)

PRIM = ["Basal", "Secretory", "Ciliated"]
print("\n=== per-module summary in primary epithelial types ===")
for ctn in PRIM:
    sub = df[df.cell_type == ctn]
    print(f"\n-- {ctn} (n={sub.n_cells.iloc[0]})")
    s = sub.groupby("module").agg(median_frac=("frac_detected", "median"),
                                  min_frac=("frac_detected", "min"),
                                  genes_ge_5pct=("frac_detected", lambda x: int((x >= 0.05).sum())),
                                  n=("gene", "size"))
    print(s.round(3).to_string())

print("\n=== M12 gene detail ===")
m12 = df[df.module == "M12_inflammatory_sasp"]
print(m12.pivot_table(index="gene", columns="cell_type", values="frac_detected")
      .reindex(columns=[c for c in PRIM if c in set(df.cell_type)]).round(4).to_string())
print("\n=== M11 gene detail ===")
m11 = df[df.module == "M11_senescence_arrest"]
print(m11.pivot_table(index="gene", columns="cell_type", values="frac_detected")
      .reindex(columns=[c for c in PRIM if c in set(df.cell_type)]).round(4).to_string())
