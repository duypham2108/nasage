# Pre-registered plan: does panel v2 hold outside the cohorts it was built from?

Locked 2026-11-25, before any held-out cohort is scored under v2. Closes the last item of
`docs/WEAKNESSES.md` "Still open": *the v2 panel exists but has never been tested on a cohort not
used to build it.*

## 1. The problem

`validation/build_panel_v2.py` derived v2 from **expression structure in Ziegler 2021 and Loske
2022 only**. No outcome variable was read, which is what makes a later test admissible — but the
structure itself was fitted to two cohorts, and the two central claims of v2 are structural:

1. **M05, M06, M07, M09 and M03 are one dimension**, not five, because they correlate 0.72–0.93
   despite sharing no genes (`docs/MODULE_DESIGN_AUDIT.md` §1.1).
2. **The declared family of 14 modules behaves as 7–10 effective tests**, so BH over 14 is wrong.

Both could be properties of those two cohorts rather than of nasal epithelium. Until they
replicate somewhere else, v2 is a reorganisation fitted to its own training data.

## 2. Held-out cohorts

None of these contributed to the derivation.

| Cohort | Assay | Scores file |
|---|---|---|
| GSE336399 | bulk nasal brushing | `results/gse336399/module_scores.csv` |
| Thwaites 2023 LAIV | bulk nasal curettage | `results/thwaites2023_laiv/module_scores_all_samples.csv` |
| Crownlands | single-cell, donor × cell type | `results/crownlands/cell_scores_obs.csv.gz` |

Crownlands is used for **structure only**. It remains withdrawn for biological claims over the
CDKN2A depth confound (`WEAKNESSES.md` C5); that withdrawal concerns age effects, not whether two
modules covary, so a correlation structure read from it is admissible and is labelled as such.

## 3. Pre-specified predictions

| # | Prediction | Pass criterion |
|---|---|---|
| **V1** | The five A01 members form a tight axis | median pairwise Spearman ρ ≥ 0.60 among members, in ≥ 2 of 3 held-out cohorts |
| **V2** | That axis is tighter than chance | observed median ρ exceeds the 95th percentile of 1,000 random 5-module subsets of the same panel |
| **V3** | v2 reduces effective dimensionality | Li & Ji effective tests under v2 < under v1 in every held-out cohort |
| **V4** | M12 remains unmeasurable | < 60% of M12 genes detected, in every held-out cohort where gene-level detection is available |

ρ ≥ 0.60 rather than the 0.75 used in derivation: a threshold chosen on training data must be
relaxed for a held-out test, or the test is asking the new cohorts to reproduce a value that was
selected to be extreme. This relaxation is declared here, before any held-out number is seen.

## 4. Verdicts, fixed in advance

- **V1 and V2 both pass** → the IFN/antigen-presentation axis is a property of nasal tissue, not
  of two cohorts. v2's collapse is justified and may be used for correction going forward.
- **V1 passes, V2 fails** → the members covary, but no more than any five modules do. That is the
  dominant-axis problem (B1) again, not a specific axis. v2's collapse is then **not** a special
  claim about interferon biology and must be described as a global-axis correction.
- **V1 fails** → the collapse does not replicate. v2 is reported as unsupported out of sample and
  v1 remains the panel of record.
- **V3 fails** → the effective-test reduction is cohort-specific; report both counts.

## 5. What this cannot establish

v2 has never been used for a primary test, and nothing here tests whether it improves any
biological result. A pass means only that the structural rearrangement generalises.

## 6. Outputs

`D:\research\nasage_validation\results\panel_v2_heldout\`. Script:
`validation/panel_v2_heldout.py`.
