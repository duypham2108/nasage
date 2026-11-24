"""Label-permutation null for the MVRI outcome result.

The matched random-gene-set control in mvri_outcome.py asks whether the *gene set* is
special. This asks a different question: whether the *pipeline* -- scoring, composition
adjustment, per-cohort logistic fits, fixed-effect pooling -- produces a p as small as the
observed one when the outcome labels carry no information.

Outcome labels are permuted within cohort, so cohort size, outcome prevalence and all
expression data are preserved exactly.

Usage: python validation/mvri_outcome_permute.py <data_dir> <symbol_map.h5ad> <out_dir> [n_perm]
"""
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

import mvri_outcome as mo   # importing runs the primary analysis, then we reuse its objects

warnings.filterwarnings("ignore")
N_PERM = int(sys.argv[4]) if len(sys.argv) > 4 else 1000
OUT = Path(sys.argv[3])
PROGRAM = "P5_mucin_barrier"

obs = mo.R.loc[mo.R.program == PROGRAM].iloc[0]
pcol = f"proxy_{mo.PROGRAMS[PROGRAM]['cell_type']}"
print(f"\nobserved pooled logOR {obs.pooled_logOR:+.4f}  p {obs.p:.5f}")

rng = np.random.default_rng(0)
pool_b, pool_p = [], []
for _ in range(N_PERM):
    bs, ss = [], []
    for label, c in mo.cohorts.items():
        d = c["D"].copy()
        d["y"] = rng.permutation(d["y"].values)
        b, se, _, _ = mo.logistic_or(d, PROGRAM, pcol)
        bs.append(b); ss.append(se)
    pb, _, pp, _, _ = mo.meta_fixed(bs, ss)
    if np.isfinite(pb):
        pool_b.append(pb); pool_p.append(pp)

pool_b, pool_p = np.array(pool_b), np.array(pool_p)
res = dict(program=PROGRAM, n_perm=len(pool_b),
           observed_logOR=float(obs.pooled_logOR), observed_p=float(obs.p),
           perm_p_by_effect=float((np.abs(pool_b) >= abs(obs.pooled_logOR)).mean()),
           perm_p_by_pvalue=float((pool_p <= obs.p).mean()),
           perm_frac_p05=float((pool_p < 0.05).mean()))
print(f"permutations: {res['n_perm']}")
print(f"  fraction with |pooled logOR| >= observed : {res['perm_p_by_effect']:.4f}")
print(f"  fraction with pooled p <= observed       : {res['perm_p_by_pvalue']:.4f}")
print(f"  fraction reaching p < 0.05 (calibration) : {res['perm_frac_p05']:.4f}  (nominal 0.05)")
pd.DataFrame([res]).to_csv(OUT / "label_permutation.csv", index=False)
pd.DataFrame({"logOR": pool_b, "p": pool_p}).to_csv(
    OUT / "label_permutation_draws.csv.gz", index=False)
print("wrote", OUT / "label_permutation.csv")
