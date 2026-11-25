"""Project-wide test registry: how many tests has this project actually run?

Addresses WEAKNESSES.md B8. BH has always been applied inside declared families; nobody has
counted the tests across the whole programme of work. This walks every results table, collects
every p-value with its plan and family, and reports the global picture.

It does NOT retroactively re-correct published results -- family-wise correction inside a
pre-registered family is the right procedure and stays. The registry answers a different
question: how much multiplicity is there in total, and how many findings survive if the whole
project is treated as one family.

Usage: python validation/test_registry.py <results_dir> <out_dir>
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

R, OUT = Path(sys.argv[1]), Path(sys.argv[2])
OUT.mkdir(parents=True, exist_ok=True)

# (directory, file, p-column, label-columns, plan)
SOURCES = [
    ("mvri_outcome", "per_cohort.csv", "p", ["program", "cohort"], "PREREG_mvri"),
    ("mvri_engagement", "engagement_interaction.csv", "p",
     ["dataset", "program", "timepoint"], "PREREG_mvri_engagement"),
    ("nova", "primary.csv", "p_two_sided", ["program"], "PREREG_nova_p5"),
    ("nova", "secondary.csv", "p_two_sided", ["program"], "PREREG_nova_p5"),
    ("celltype_spec", "meta.csv", "p", ["program"], "PREREG_celltype_mvri_A"),
    ("laiv_engagement", "engagement.csv", "p", ["program"], "PREREG_celltype_mvri_B"),
    ("lifespan", "trajectory.csv", "p_age", ["axis", "family"], "PREREG_lifespan"),
    ("lifespan_repl", "replication_estimates.csv", "p", ["cohort", "axis", "family"],
     "PREREG_lifespan_repl"),
    ("method_improve", "union_BH.csv", "p", ["approach", "unit"], "PREREG_method_improvements"),
    ("multiverse", "specification_curve.csv", "p",
     ["program", "scoring", "covars"], "PREREG_multiverse"),
    ("multiverse", "specification_curve_rsv_pneumo.csv", "p",
     ["cohort", "program", "scoring", "covars"], "PREREG_multiverse"),
]

rows = []
for d, f, pcol, labs, plan in SOURCES:
    fp = R / d / f
    if not fp.exists():
        print(f"  missing: {fp}")
        continue
    t = pd.read_csv(fp)
    if pcol not in t.columns:
        print(f"  no column {pcol} in {fp}")
        continue
    have = [c for c in labs if c in t.columns]
    for _, r in t.iterrows():
        p = r[pcol]
        if not np.isfinite(p):
            continue
        rows.append(dict(plan=plan, source=f"{d}/{f}",
                         label=" | ".join(str(r[c]) for c in have), p=float(p)))

T = pd.DataFrame(rows)
if not len(T):
    sys.exit("no tests collected")

# Exploratory sweeps are counted separately: a specification curve is 16 views of ONE question,
# not 16 questions. Counting them as independent tests would be as wrong as ignoring them.
T["kind"] = np.where(T.plan.isin(["PREREG_multiverse", "PREREG_method_improvements"]),
                     "sensitivity_sweep", "primary_or_secondary")

print("=" * 78)
print("PROJECT-WIDE TEST REGISTRY")
print("=" * 78)
print(T.groupby(["kind", "plan"]).size().rename("n_tests").to_string())
print(f"\nTOTAL p-values recorded: {len(T)}")
print(T.groupby('kind').size().rename('n').to_string())


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = np.empty_like(p)
    q[o] = np.minimum.accumulate((p[o] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    return np.minimum(q, 1.0)


H = T[T.kind == "primary_or_secondary"].copy()
H["BH_project"] = bh(H.p.values)
H = H.sort_values("p")
H.to_csv(OUT / "registry_primary.csv", index=False)
T.to_csv(OUT / "registry_all.csv", index=False)

print("\n" + "=" * 78)
print(f"BH across ALL {len(H)} primary/secondary tests treated as ONE family")
print("=" * 78)
print(H.head(20).round(5).to_string(index=False))
print(f"\nsurvive project-wide BH < 0.05: {int((H.BH_project < 0.05).sum())} of {len(H)}")
print(f"nominally p < 0.05:              {int((H.p < 0.05).sum())} of {len(H)}")
print(f"expected false positives at 5%:  {0.05 * len(H):.1f}")
print("\nwrote", OUT)
