"""Does the CDKN2A depth confound survive an HVG-based pipeline?

The manuscript's figures are built on top highly variable genes, not on the curated
modules. This repeats the Crownlands CDKN2A+ contrasts the HVG way, so the question
"is p16-positivity a detection event" is answered in the pipeline actually used.

Steps, on the fixed ciliated subsample:
  1. HVGs by scanpy seurat_v3 on raw counts (n_top configurable).
  2. How much of each curated module is HVG at all.
  3. Depth of CDKN2A+ vs CDKN2A- cells (a property of the grouping, not of the genes).
  4. CDKN2A+ vs CDKN2A- contrast, within donor, scored on HVGs only:
       a. a whole-transcriptome-free readout: mean z of HVG expression (the thing an
          HVG-based heatmap or PCA shows);
       b. module genes restricted to HVGs;
     each with and without within-donor depth matching.

Usage: python validation/crownlands_hvg_check.py <subsample.h5ad> <scores.csv.gz> <out_dir> [n_top]
"""
import sys
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy import stats

SUB, SCORES, OUT = sys.argv[1], sys.argv[2], Path(sys.argv[3])
N_TOP = int(sys.argv[4]) if len(sys.argv) > 4 else 2000
OUT.mkdir(parents=True, exist_ok=True)

import nasage as nas
ALIAS = {"CIMAP3": "ODF3B", "STING1": "TMEM173"}

f = h5py.File(SUB, "r")
sym = np.array([x.decode() if isinstance(x, bytes) else x for x in f["var_symbol"][:]])
ct = np.array([x.decode() for x in f["obs/ct"][:]])
cil = np.flatnonzero(ct == "ciliated")
ip = f["X/indptr"][:].astype(np.int64)
d_, i_, l_ = [], [], []
for r in cil:
    lo, hi = ip[r], ip[r + 1]
    d_.append(f["X/data"][lo:hi]); i_.append(f["X/indices"][lo:hi]); l_.append(hi - lo)
X = sp.csr_matrix((np.concatenate(d_), np.concatenate(i_),
                   np.concatenate([[0], np.cumsum(l_)])), shape=(len(cil), len(sym)))
f.close()

S = pd.read_csv(SCORES).iloc[cil].reset_index(drop=True)
A = ad.AnnData(X=X, obs=S.reset_index(drop=True),
               var=pd.DataFrame(index=pd.Index(sym)))
A.var_names_make_unique()
j = int(np.flatnonzero(A.var_names == "CDKN2A")[0])
A.obs["pos"] = np.asarray(A.X[:, j].todense()).ravel() >= 1
A.obs["depth"] = A.obs.nis_n_genes
print(f"ciliated cells {A.n_obs:,}; CDKN2A+ {int(A.obs.pos.sum()):,} ({A.obs.pos.mean():.1%})")

# ---- 1. HVGs, the way an HVG pipeline picks them
# scanpy's default "seurat" flavour on log-normalised counts; seurat_v3 would need
# scikit-misc, which is not installed. Both are standard; this one needs no new package.
norm = A.copy()
sc.pp.normalize_total(norm, target_sum=1e4)
sc.pp.log1p(norm)
sc.pp.highly_variable_genes(norm, flavor="seurat", n_top_genes=N_TOP, batch_key="donor")
hvg = norm.var_names[norm.var.highly_variable].tolist()
print(f"HVGs (seurat flavour, n_top={N_TOP}, batch=donor): {len(hvg)}")
print(f"CDKN2A is an HVG: {'CDKN2A' in set(hvg)}")

# ---- 2. module overlap with the HVG set
M = nas.load_modules()
rows = []
for m, dd in M.items():
    g = [ALIAS.get(x, x) for x in dd["genes"]]
    present = [x for x in g if x in set(A.var_names)]
    inhvg = [x for x in present if x in set(hvg)]
    rows.append(dict(module=m, n=len(g), present=len(present), in_hvg=len(inhvg),
                     hvg_genes=";".join(inhvg)))
MO = pd.DataFrame(rows)
MO.to_csv(OUT / "module_hvg_overlap.csv", index=False)
print("\nmodule genes that are HVGs:")
print(MO[["module", "present", "in_hvg"]].to_string(index=False))

# ---- 3. the grouping itself
print("\n=== depth of the CDKN2A+ grouping (independent of gene selection) ===")
print(f"median detected genes: CDKN2A+ {A.obs.depth[A.obs.pos].median():.0f} vs "
      f"CDKN2A- {A.obs.depth[~A.obs.pos].median():.0f}")
pd_ = A.obs.groupby("donor").apply(
    lambda g: g.depth[g.pos].median() - g.depth[~g.pos].median(), include_groups=False)
print(f"per-donor difference: median {pd_.median():+.0f} genes, "
      f"positive in {int((pd_ > 0).sum())}/{len(pd_)} donors")

# ---- 4. HVG-based readouts, matched and unmatched
sub = norm[:, hvg].copy()
sc.pp.scale(sub, max_value=10)
del norm
A.obs["hvg_mean_z"] = np.asarray(sub.X).mean(1)          # what an HVG heatmap/PCA shows

for mod, tag in (("M02_mucociliary_axoneme", "m02_hvg"),
                 ("M11_senescence_arrest", "m11_hvg")):
    g = [x for x in MO.loc[MO.module == mod, "hvg_genes"].iloc[0].split(";") if x]
    g = [x for x in g if x != "CDKN2A"]                   # never score the defining gene
    A.obs[tag] = (np.asarray(sub[:, g].X).mean(1) if g else np.nan)
    print(f"{mod}: {len(g)} HVG genes used ({','.join(g) if g else 'none'})")

cols = ["hvg_mean_z", "m02_hvg", "m11_hvg"]
rng = np.random.default_rng(0)
raw, matched = [], []
for dn, g in A.obs.groupby("donor"):
    a, b = g[g.pos], g[~g.pos]
    if len(a) >= 20 and len(b) >= 20:
        raw.append(dict(donor=dn, **{c: a[c].median() - b[c].median() for c in cols}))
    gg = g.copy()
    try:
        gg["bin"] = pd.qcut(gg.depth, 10, labels=False, duplicates="drop")
    except ValueError:
        continue
    ai, bi = [], []
    for _, gb in gg.groupby("bin"):
        pa, pb = gb.index[gb.pos], gb.index[~gb.pos]
        n = min(len(pa), len(pb))
        if n:
            ai += list(rng.choice(pa, n, replace=False)); bi += list(rng.choice(pb, n, replace=False))
    if len(ai) >= 20:
        A_, B_ = gg.loc[ai], gg.loc[bi]
        matched.append(dict(donor=dn, dd=A_.depth.median() - B_.depth.median(),
                            **{c: A_[c].median() - B_[c].median() for c in cols}))
R, Mt = pd.DataFrame(raw), pd.DataFrame(matched)
R.to_csv(OUT / "hvg_contrast_raw.csv", index=False)
Mt.to_csv(OUT / "hvg_contrast_matched.csv", index=False)

print(f"\n=== CDKN2A+ vs CDKN2A- on HVGs (n={len(R)} raw, {len(Mt)} matched; "
      f"residual depth {Mt.dd.median():+.0f} genes) ===")
print(f"{'readout':14s} {'raw delta':>11s} {'raw p':>11s} {'matched delta':>14s} {'matched p':>11s}")
for c in cols:
    pr = stats.wilcoxon(R[c].dropna()).pvalue if R[c].notna().any() else np.nan
    pm = stats.wilcoxon(Mt[c].dropna()).pvalue if Mt[c].notna().any() else np.nan
    print(f"{c:14s} {R[c].median():+11.4f} {pr:11.3g} {Mt[c].median():+14.4f} {pm:11.3g}")
print("\nwrote", OUT)
