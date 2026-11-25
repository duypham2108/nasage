"""Derive epithelial cell-type programs for the mucosal vaccine readiness screen.

Implements docs/PREREG_mvri.md section 2, locked at commit fecab57 before any program
was computed. No outcome data is read here, by construction: this script opens exactly
one input file, the primary cohort h5ad.

Usage: python validation/mvri_derive.py <primary.h5ad> <out_dir>
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import igraph as ig
import leidenalg

import nasage as nas

H5AD, OUT = sys.argv[1], Path(sys.argv[2])
OUT.mkdir(parents=True, exist_ok=True)

CELL_TYPES = ["Basal", "Secretory", "Ciliated"]
MIN_CELLS_PER_DONOR = 50
DETECT_FRAC = 0.10          # plan 2.2
DONOR_FRAC = 0.80           # plan 2.2
KNN, RESOLUTION, SEED = 15, 1.0, 0
SIZE_MIN, SIZE_MAX = 15, 150
MAX_FAMILY = 8
JACCARD_REDUNDANT = 0.5

PROXY = {  # plan 3: forbidden from entering any program
    "Ciliated": ["FOXJ1", "TPPP3", "PIFO", "SNTN", "CAPS"],
    "Basal": ["KRT5", "KRT14", "TP63", "DLK2", "MIR205HG"],
    "Secretory": ["SCGB1A1", "BPIFA1", "MSMB", "WFDC2", "SLPI"],
}
FORBIDDEN = {g for v in PROXY.values() for g in v}

CURATED = {m: set(d["genes"]) for m, d in nas.load_modules().items()}


def rankdata_cols(X):
    """Column-wise ranks of a dense 2-D array, average ties."""
    order = np.argsort(X, axis=0, kind="mergesort")
    ranks = np.empty_like(X, dtype=np.float64)
    n = X.shape[0]
    idx = np.arange(n)
    for j in range(X.shape[1]):
        col = X[order[:, j], j]
        r = idx.astype(np.float64)
        # average ties
        i = 0
        while i < n:
            k = i
            while k + 1 < n and col[k + 1] == col[i]:
                k += 1
            if k > i:
                r[i:k + 1] = (i + k) / 2.0
            i = k + 1
        ranks[order[:, j], j] = r
    return ranks


def corr_within_donor(Xd):
    """Spearman gene-gene correlation within one donor. Xd: cells x genes dense."""
    R = rankdata_cols(Xd)
    R -= R.mean(0)
    sd = np.sqrt((R ** 2).sum(0))
    sd[sd == 0] = np.inf
    R /= sd
    return R.T @ R


print(f"reading {H5AD}")
A = ad.read_h5ad(H5AD)
ct_col = "Cell_type_annotation_level2"
print(f"{A.shape[0]} cells, {A.shape[1]} genes, {A.obs.donor.nunique()} donors")

rows = []
programs = {}
for ct in CELL_TYPES:
    sub = A[A.obs[ct_col] == ct]
    counts = sub.obs.donor.value_counts()
    donors = counts[counts >= MIN_CELLS_PER_DONOR].index.tolist()
    sub = sub[sub.obs.donor.isin(donors)]
    print(f"\n=== {ct}: {sub.shape[0]} cells, {len(donors)} donors "
          f"(>= {MIN_CELLS_PER_DONOR} cells)")

    X = sub.X.tocsr() if sp.issparse(sub.X) else sp.csr_matrix(sub.X)
    # per-donor detection fraction
    det_per_donor = []
    for d in donors:
        m = (sub.obs.donor == d).values
        Xi = X[m]
        det = np.asarray((Xi > 0).sum(0)).ravel() / Xi.shape[0]
        det_per_donor.append(det >= DETECT_FRAC)
    det_per_donor = np.vstack(det_per_donor)
    keep = det_per_donor.mean(0) >= DONOR_FRAC
    genes = sub.var_names[keep]
    genes = pd.Index([g for g in genes if g not in FORBIDDEN])
    print(f"gene universe: {len(genes)} genes detected in >= {DETECT_FRAC:.0%} of cells "
          f"in >= {DONOR_FRAC:.0%} of donors ({len(FORBIDDEN & set(sub.var_names[keep]))} "
          f"proxy markers removed)")

    gi = sub.var_names.get_indexer(genes)
    # log1p CP10K
    C = sp.csr_matrix(X[:, gi], dtype=np.float32)
    tot = np.asarray(X.sum(1)).ravel()
    tot[tot == 0] = 1
    C = sp.diags(1e4 / tot) @ C
    C = C.tocsr()
    C.data = np.log1p(C.data)

    acc = np.zeros((len(genes), len(genes)), dtype=np.float64)
    for d in donors:
        m = (sub.obs.donor == d).values
        acc += corr_within_donor(np.asarray(C[m].todense(), dtype=np.float64))
    Corr = acc / len(donors)
    np.fill_diagonal(Corr, 0.0)

    # kNN gene graph on the averaged within-donor correlation
    nn = np.argsort(-Corr, axis=1)[:, :KNN]
    edges, weights = [], []
    for i in range(len(genes)):
        for j in nn[i]:
            if Corr[i, j] > 0:
                edges.append((i, int(j)))
                weights.append(float(Corr[i, j]))
    g = ig.Graph(n=len(genes), edges=edges, directed=False)
    g.es["weight"] = weights
    g.simplify(combine_edges="max")
    part = leidenalg.find_partition(g, leidenalg.RBConfigurationVertexPartition,
                                    weights="weight", resolution_parameter=RESOLUTION,
                                    seed=SEED)
    lab = np.array(part.membership)
    print(f"leiden: {lab.max() + 1} raw clusters, sizes "
          f"{sorted(np.bincount(lab))[::-1][:8]}")

    for c in range(lab.max() + 1):
        gs = list(genes[lab == c])
        if not (SIZE_MIN <= len(gs) <= SIZE_MAX):
            continue
        idx = np.where(lab == c)[0]
        coh = float(Corr[np.ix_(idx, idx)][np.triu_indices(len(idx), 1)].mean())
        best_j, best_m = 0.0, ""
        for m, mg in CURATED.items():
            inter = len(mg & set(gs))
            jac = inter / len(mg | set(gs)) if (mg or gs) else 0.0
            if jac > best_j:
                best_j, best_m = jac, m
        name = f"P_{ct[:3].lower()}_{c:02d}"
        rows.append(dict(program=name, cell_type=ct, n_genes=len(gs), coherence=coh,
                         max_jaccard=best_j, nearest_module=best_m,
                         redundant=best_j > JACCARD_REDUNDANT))
        programs[name] = gs

P = pd.DataFrame(rows)
if P.empty:
    raise SystemExit("no program passed the size filter; plan section 7 applies")
P = P.sort_values(["cell_type", "coherence"], ascending=[True, False])
print("\n=== all programs passing the size filter ===")
print(P.to_string(index=False))

kept = P[~P.redundant].copy()
print(f"\ndropped as redundant with a curated module (Jaccard > {JACCARD_REDUNDANT}): "
      f"{int(P.redundant.sum())}")

# plan 2.1: round-robin across cell types, highest coherence first
queues = {ct: kept[kept.cell_type == ct].program.tolist() for ct in CELL_TYPES}
family, i = [], 0
while len(family) < MAX_FAMILY and any(queues.values()):
    ct = CELL_TYPES[i % len(CELL_TYPES)]
    if queues[ct]:
        family.append(queues[ct].pop(0))
    i += 1
P["in_family"] = P.program.isin(family)

print(f"\n=== FAMILY, fixed before any outcome is touched: {len(family)} programs ===")
print(P[P.in_family].to_string(index=False))
if len(family) < 3:
    print("\nFEWER THAN 3 PROGRAMS: plan section 7 says this is a derivation failure.")

P.to_csv(OUT / "programs_summary.csv", index=False)
with open(OUT / "programs_genes.csv", "w", encoding="utf-8") as f:
    f.write("program,cell_type,in_family,genes\n")
    for _, r in P.iterrows():
        f.write(f"{r.program},{r.cell_type},{r.in_family},"
                f"{';'.join(sorted(programs[r.program]))}\n")
print("\nwrote", OUT)
