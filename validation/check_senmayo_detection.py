"""Detection of the SenMayo (SAUL_SEN_MAYO, 125 genes) senescence/SASP set in the
primary cohort, per cell type, compared with the curated 7-gene M12.

Decides whether Figure 3 can say "SASP" or only "senescence arrest".
"""
import sys, h5py, numpy as np, pandas as pd

F, GRP, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
M12 = ["CXCL1", "CXCL8", "IL1A", "IL1B", "IL6", "TNF", "CCL2"]
M11 = ["CDKN1A", "CDKN2A", "GDF15", "IGFBP7", "SERPINE1"]

sen = [l.strip() for l in open(GRP) if l.strip() and not l.startswith("#")][1:]
print(f"SenMayo: {len(sen)} genes")

f = h5py.File(F, "r")
var = np.array([x.decode() if isinstance(x, bytes) else x for x in f["var"]["_index"][:]])
g2i = {g: i for i, g in enumerate(var)}

sets = {"SenMayo": sen, "M12": M12, "M11": M11}
genes = sorted({g for s in sets.values() for g in s if g in g2i})
print("present in data:", {k: sum(g in g2i for g in v) for k, v in sets.items()},
      "of", {k: len(v) for k, v in sets.items()})

gidx = np.array([g2i[g] for g in genes])
lut = np.full(len(var), -1, np.int32)
lut[gidx] = np.arange(len(genes))

o = f["obs"]["Cell_type_annotation_level2"]
raw = np.array([x.decode() if isinstance(x, bytes) else x for x in o["categories"][:]])
ct = raw[o["codes"][:]]                      # same path as the measurability script
cats, ccode = np.unique(ct, return_inverse=True)
nct = len(cats)

X = f["X"]
indptr = X["indptr"][:].astype(np.int64)
n = len(indptr) - 1
det = np.zeros((nct, len(genes)), np.int64)
ncell = np.bincount(ccode, minlength=nct)

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
frac.to_csv(OUT)

PRIM = ["Basal", "Secretory", "Ciliated", "Goblet"]
print("\n=== genes reaching >=5% detection ===")
rows = []
for name, gl in sets.items():
    gl = [g for g in gl if g in genes]
    for ct in PRIM:
        if ct not in frac.index:
            continue
        v = frac.loc[ct, gl]
        rows.append(dict(set=name, cell_type=ct, n_avail=len(gl),
                         n_ge5pct=int((v >= 0.05).sum()),
                         pct_ge5pct=round(100 * (v >= 0.05).mean(), 1),
                         median_frac=round(float(v.median()), 4)))
print(pd.DataFrame(rows).to_string(index=False))

print("\n=== SenMayo genes clearing 5% in Ciliated (the Fig 3 population) ===")
if "Ciliated" in frac.index:
    v = frac.loc["Ciliated", [g for g in sen if g in genes]].sort_values(ascending=False)
    print(v[v >= 0.05].round(3).to_string())
    print(f"\n{int((v >= 0.05).sum())} of {len(v)} available SenMayo genes clear 5% in ciliated")
