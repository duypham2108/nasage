"""Audit the current module design: redundancy, effective dimensionality, size, measurability.

This is a DIAGNOSTIC on the module panel itself, not a test of any hypothesis. It reads module
definitions and already-computed donor x cell-type scores; it fits no outcome model and touches
no outcome variable, so it cannot leak outcome information into module design.

Usage: python validation/module_design_audit.py <results_dir> <out_dir>
"""
import itertools
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from nasage.modules import load_modules          # noqa: E402
from mvri_programs import PROGRAMS               # noqa: E402

R, OUT = Path(sys.argv[1]), Path(sys.argv[2])
OUT.mkdir(parents=True, exist_ok=True)

SC_COHORTS = ["ziegler2021_v2", "loske2022_v2"]
CELL_TYPES = ["ciliated", "secretory", "basal"]

# ---------------------------------------------------------------- 1. size and gene overlap
sets = {}
for k, v in load_modules().items():
    sets[k] = set(v["genes"] if isinstance(v, dict) else v)
for k, v in PROGRAMS.items():
    sets["*" + k] = set(v["genes"])

size = pd.DataFrame([dict(set_id=k, n_genes=len(v)) for k, v in sets.items()]
                    ).sort_values("n_genes")
size.to_csv(OUT / "module_sizes.csv", index=False)

jac = []
for a, b in itertools.combinations(sorted(sets), 2):
    A, B = sets[a], sets[b]
    u = A | B
    if u:
        jac.append(dict(a=a, b=b, jaccard=len(A & B) / len(u), shared=len(A & B)))
J = pd.DataFrame(jac).sort_values("jaccard", ascending=False)
J.to_csv(OUT / "gene_overlap.csv", index=False)

# ---------------------------------------------------------------- 2. score redundancy
frames = []
for coh in SC_COHORTS:
    f = R / coh / "donor_celltype_scores.csv"
    if f.exists():
        d = pd.read_csv(f)
        d["cohort"] = coh
        frames.append(d)
D = pd.concat(frames, ignore_index=True)
mcols = [c for c in D.columns if c.startswith("M")]

corr_rows, dim_rows = [], []
for ct in CELL_TYPES:
    s = D[(D.cell_type == ct) & (D.eligible == True)][mcols]   # noqa: E712
    s = s.loc[:, s.notna().mean() > 0.8].dropna()
    if len(s) < 10:
        continue
    C = s.corr(method="spearman")
    for a, b in itertools.combinations(C.columns, 2):
        corr_rows.append(dict(cell_type=ct, a=a, b=b, rho=C.loc[a, b], n=len(s)))
    ev = np.sort(np.linalg.eigvalsh(C.values))[::-1]
    ev[ev < 0] = 0
    k = len(ev)
    li_ji = float(sum((1 if e >= 1 else 0) + (e - np.floor(e)) for e in ev))
    galwey = float(np.sqrt(ev).sum() ** 2 / ev.sum())
    dim_rows.append(dict(cell_type=ct, n_donors=len(s), k_modules=k,
                         eff_tests_li_ji=li_ji, eff_tests_galwey=galwey,
                         pc1_var_frac=float(ev[0] / ev.sum())))
CORR = pd.DataFrame(corr_rows)
DIM = pd.DataFrame(dim_rows)
CORR.to_csv(OUT / "score_correlation.csv", index=False)
DIM.to_csv(OUT / "effective_dimensionality.csv", index=False)

# redundancy that gene overlap misses: high |rho|, low Jaccard
jmap = {(r.a, r.b): r.jaccard for r in J.itertuples()}
hidden = []
for r in CORR.itertuples():
    a, b = r.a.replace("_UCell", "").replace("_logexpr", "").replace("_logit", ""), \
           r.b.replace("_UCell", "").replace("_logexpr", "").replace("_logit", "")
    j = jmap.get((a, b), jmap.get((b, a), 0.0))
    if abs(r.rho) > 0.70:
        hidden.append(dict(cell_type=r.cell_type, a=a, b=b, rho=r.rho, jaccard=j,
                           missed_by_jaccard=bool(j <= 0.5)))
H = pd.DataFrame(hidden).sort_values("rho", key=abs, ascending=False)
H.to_csv(OUT / "redundancy_missed_by_jaccard.csv", index=False)

print("=== SIZE: sets with fewer than 8 genes ===")
print(size[size.n_genes < 8].to_string(index=False))
print("\n=== GENE OVERLAP: Jaccard > 0.3 ===")
print(J[J.jaccard > 0.3].round(3).to_string(index=False))
print("\n=== EFFECTIVE DIMENSIONALITY of the module panel ===")
print(DIM.round(2).to_string(index=False))
print("\n=== REDUNDANCY MISSED BY THE JACCARD RULE (|rho| > 0.7, Jaccard <= 0.5) ===")
print(H[H.missed_by_jaccard].round(3).to_string(index=False))
print("\nwrote", OUT)
