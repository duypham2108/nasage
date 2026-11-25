"""Pre-specified SECONDARY endpoint and SENSITIVITY analyses for GSE155237.

Written after seeing the primary results, but implements ONLY what was already
specified in the docstring of validation/rsv_challenge_gse155237.py (commit
5d0b924): the `symptom` endpoint, and sensitivity variants (a) no PTPRC term,
(b) no batch term, (c) new-batch donors only, (d) UCell scores. Nothing is added
and nothing is dropped; every variant is reported.

Usage: python validation/rsv_secondary_sensitivity.py <work_dir> <out_dir>
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from rsv_challenge_gse155237 import (DESCRIPTIVE, POSCTRL, PRIMARY, SECONDARY,  # noqa: E402
                                     fit, load_counts, score_modules)

import nasage as nas  # noqa: E402

work, out = Path(sys.argv[1]), Path(sys.argv[2])
meta = pd.read_csv(work / "g155_samples.csv")
counts = load_counts(work)
cols = pd.Index(counts.columns.astype(str))
meta = meta.set_index(meta.title.astype(str).values).loc[cols[cols.isin(meta.title.astype(str))]]
counts = counts[cols[cols.isin(meta.index)]]

keep = (meta.timepoint == "D0").values
md, C = meta[keep], counts.loc[:, keep]
cpm = C / C.sum(0) * 1e6
logc = np.log2(cpm[cpm.mean(1) >= 1] + 1)
M = {m: d["genes"] for m, d in nas.load_modules().items()}
S = score_modules(logc, M)
S.index = md.index

print("symptom x outcome:\n", pd.crosstab(md.symptom, md.outcome), "\n")

cov_full = pd.DataFrame({"batch": md.batch.values,
                         "log10_counts": np.log10(C.sum(0).values),
                         "PTPRC": logc.loc["PTPRC"].values}, index=S.index)

# UCell sensitivity
try:
    import anndata as ad
    from pyucell import compute_ucell_scores
    from nasage.scoring import MAX_RANK
    u = ad.AnnData(X=C.T.values.astype(np.float32), obs=pd.DataFrame(index=C.columns),
                   var=pd.DataFrame(index=C.index))
    compute_ucell_scores(u, {m: [g for g in g_ if g in C.index] for m, g_ in M.items()
                             if len(g_) > 1}, max_rank=MAX_RANK, suffix="", n_jobs=1)
    U = u.obs; U.index = S.index
except Exception as exc:
    print("UCell skipped:", exc); U = None

tested = {**PRIMARY, **SECONDARY, **POSCTRL}
endpoints = {"infected": (md.outcome.values == "Infected").astype(int),
             "symptom_cold": (md.symptom.values == "cold").astype(int)}
rows = []
for ep, yv in endpoints.items():
    y = pd.Series(yv, index=S.index)
    for m, d in {**tested, **{k: "two" for k in DESCRIPTIVE}}.items():
        if m not in S:
            continue
        variants = {
            "main": (S[m], cov_full, slice(None)),
            "no_PTPRC": (S[m], cov_full.drop(columns="PTPRC"), slice(None)),
            "no_batch": (S[m], cov_full.drop(columns="batch"), slice(None)),
            "new_batch_only": (S[m], cov_full.drop(columns="batch"), md.batch.values == "new"),
        }
        if U is not None and m in U:
            variants["ucell"] = (U[m], cov_full, slice(None))
        for vname, (sc, cv, sel) in variants.items():
            r = fit(y[sel], sc[sel], cv[sel], d)
            r.update(endpoint=ep, module=m, variant=vname, predicted=d,
                     n=int(np.sum(sel)) if isinstance(sel, np.ndarray) else int(len(y)),
                     tier=("primary" if m in PRIMARY else "secondary" if m in SECONDARY
                           else "positive_control" if m in POSCTRL else "descriptive"))
            rows.append(r)

R = pd.DataFrame(rows)
R.to_csv(out / "rsv_secondary_sensitivity.csv", index=False)
pd.set_option("display.width", 230)
cols_show = ["endpoint", "module", "tier", "variant", "n", "or_per_sd", "ci_low", "ci_high",
             "p_one_sided", "p_two_sided"]
print("=== SECONDARY ENDPOINT: symptom (matches the source paper's claim)")
print(R[(R.endpoint == "symptom_cold") & (R.variant == "main")][cols_show].round(4)
      .to_string(index=False))
print("\n=== SENSITIVITY, primary modules + positive control, both endpoints")
sub = R[R.module.isin(list(PRIMARY) + list(POSCTRL))]
print(sub[cols_show].round(4).to_string(index=False))
