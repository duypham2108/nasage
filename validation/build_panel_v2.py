"""Derive module panel v2 from the v1 panel using STRUCTURE ONLY.

Inputs used: v1 gene lists, donor x cell-type module scores (expression structure), and
measurability coverage. **No outcome variable of any kind is read.** That is the property that
makes v2 usable for a pre-registered test later; see docs/MODULE_DESIGN_AUDIT.md section 5.

v1 is not modified. v2 is written to a separate resource file, and every v1 result stays
attached to v1.

Changes applied (audit section 3, items 2, 4, 5, 6 and the correlation rule):
  1. Collapse the interferon / antigen-presentation axis into one declared dimension.
  2. Retire M12 from primary families (fails the 60% measurability gate in all four cohorts).
  3. Mark single-gene sets as gene-level readouts, not modules.
  4. Apply the score-correlation redundancy rule in place of Jaccard.

Usage: python validation/build_panel_v2.py <results_dir> <out_json> [<audit_out_dir>]
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from nasage.modules import load_modules                       # noqa: E402
from nasage.stats import effective_tests, redundancy_report    # noqa: E402

R = Path(sys.argv[1])
OUT_JSON = Path(sys.argv[2])
AUDIT = Path(sys.argv[3]) if len(sys.argv) > 3 else OUT_JSON.parent

RHO = 0.75
SC_COHORTS = ["ziegler2021_v2", "loske2022_v2"]
CELL_TYPES = ["ciliated", "secretory", "basal"]

# The axis identified in docs/MODULE_DESIGN_AUDIT.md section 1.1. Members correlate 0.72-0.93
# across all three epithelial cell types while sharing no genes, because all are
# interferon-inducible. They become declared COMPONENTS of one dimension, kept and reportable
# individually, but contributing ONE test to the primary family.
IFN_AXIS = ["M05_interferon_stimulated_genes", "M06_mhc_I_processing",
            "M07_mhc_II_presentation", "M09_plasma_secretory_machinery",
            "M03_pigr_antibody_transport"]
RETIRE_FROM_PRIMARY = ["M12_inflammatory_sasp"]   # fails the measurability gate everywhere
MIN_MODULE_SIZE = 8

v1 = load_modules()

# ---------------------------------------------------------------- evidence: score structure
frames = []
for coh in SC_COHORTS:
    f = R / coh / "donor_celltype_scores.csv"
    if f.exists():
        frames.append(pd.read_csv(f))
if not frames:
    sys.exit("no donor x cell-type score tables found; cannot derive v2 from structure")
D = pd.concat(frames, ignore_index=True)
mcols = [c for c in D.columns if c.startswith("M")]

ev_rows, red_rows = [], []
for ct in CELL_TYPES:
    s = D[(D.cell_type == ct) & (D.eligible == True)][mcols]      # noqa: E712
    s = s.loc[:, s.notna().mean() > 0.8].dropna()
    if len(s) < 10:
        continue
    e = effective_tests(s.corr(method="spearman"))
    e["cell_type"] = ct
    e["n_donors"] = len(s)
    ev_rows.append(e)
    r = redundancy_report(s, v1, rho_threshold=RHO)
    r["cell_type"] = ct
    red_rows.append(r)
EV = pd.DataFrame(ev_rows)
RED = pd.concat(red_rows, ignore_index=True) if red_rows else pd.DataFrame()
AUDIT.mkdir(parents=True, exist_ok=True)
EV.to_csv(AUDIT / "v2_effective_dimensionality_v1panel.csv", index=False)
RED.to_csv(AUDIT / "v2_redundancy_v1panel.csv", index=False)

print("=== v1 panel structure (the evidence v2 is built from) ===")
print(EV.round(2).to_string(index=False))
if len(RED):
    rr = RED[RED.redundant]
    print(f"\nredundant pairs at |rho| >= {RHO}: {len(rr)}; "
          f"missed by the Jaccard rule: {int(rr.missed_by_jaccard.sum())}")

# ---------------------------------------------------------------- build v2
v2 = {}
for mid, m in v1.items():
    e = dict(m)
    n = len(m["genes"])
    e["panel_version"] = "v2"
    e["derived_from"] = mid
    e["kind"] = "gene_readout" if n == 1 else ("small_module" if n < MIN_MODULE_SIZE
                                               else "module")
    e["axis"] = "IFN_antigen_presentation" if mid in IFN_AXIS else "independent"
    e["axis_role"] = "component" if mid in IFN_AXIS else "primary_test"
    e["primary_family"] = mid not in RETIRE_FROM_PRIMARY and mid not in IFN_AXIS
    if mid in RETIRE_FROM_PRIMARY:
        e["retired_reason"] = ("fails the 60% measurability gate in all four single-cell "
                               "cohorts (0-14% of genes detected); retained for the "
                               "measurement demonstration only")
    v2[mid] = e

# the collapsed axis, scored as the mean of its components' z-scores
members = [m for m in IFN_AXIS if m in v1]
axis_genes = sorted({g for m in members for g in v1[m]["genes"]})
v2["A01_ifn_antigen_presentation"] = dict(
    genes=axis_genes,
    name="Interferon / antigen-presentation axis (collapsed)",
    definition=("One declared dimension replacing " + ", ".join(members) + " as independent "
                "primary tests. Members correlate 0.72-0.93 across ciliated, secretory and "
                "basal cells while sharing no genes: all are interferon-inducible. Components "
                "remain individually reportable but contribute one test to the primary family."),
    anchor_go="GO:0034340 (response to type I interferon)",
    anchor_reactome="R-HSA-913531 (Interferon Signaling)",
    pmids=[],
    interpretation=("High score indicates an interferon-activated epithelial state with "
                    "coordinated antigen-presentation and secretory-transport machinery."),
    limitations=("A composite of five previously separate modules. Agreement among its "
                 "components is not independent corroboration -- that is why it exists."),
    panel_version="v2",
    derived_from=members,
    kind="axis",
    axis="IFN_antigen_presentation",
    axis_role="axis",
    primary_family=True,
)

OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
with open(OUT_JSON, "w") as fh:
    json.dump(v2, fh, indent=1)

prim = [k for k, v in v2.items() if v.get("primary_family")]
print(f"\nv2 written: {len(v2)} sets, {len(prim)} in the primary family")
print(f"  collapsed axis A01 from {len(members)} modules, {len(axis_genes)} genes")
print(f"  retired from primary: {RETIRE_FROM_PRIMARY}")
print(f"  gene-level readouts: {[k for k, v in v2.items() if v.get('kind') == 'gene_readout']}")
print(f"  small modules (<{MIN_MODULE_SIZE} genes): "
      f"{[k for k, v in v2.items() if v.get('kind') == 'small_module']}")
print("\nwrote", OUT_JSON)
print("v1 is unmodified. No outcome variable was read in building v2.")
