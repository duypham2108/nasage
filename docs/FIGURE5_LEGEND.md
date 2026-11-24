# Figure 5 and Supplementary Figure S5 — draft legends and source data

| | File |
|---|---|
| **Figure 5** (main) | `D:\research\nasage_validation\results\figure5\figure5.pdf` — 180 × 128 mm |
| **Supplementary Figure S5** | `D:\research\nasage_validation\results\figure5\figure5_supp.pdf` — 180 × 140 mm |

Both vector, Type-42 fonts (text stays editable), with 400 dpi PNGs beside them. Built by
`validation/figure5_plot.py` from `validation/figure5_data.py`.

**Revision history.** 2026-11-24 (second revision): panels S5e and S5f added, carrying the NOVA
influenza replication of the P5 mucin-barrier result (`results/nova/PREREG_RESULTS.md`). The
result met the §2.2 verdict table but **failed the §3 specificity control**, so it enters the
figure as a bounded negative, not as a third replication; the GSE253686 row in S5d records it as
NOT SPECIFIC. S5 is now 180 × 210 mm. Main Figure 5 is unchanged.

**Revision history.** 2026-11-24: the ALI panel was removed (no source data exists); Figure 4 was
narrowed to the childhood contrast; the pre-registration table's PENDING row was deleted; and the
nulls, measurement caveats and pre-registration table were moved from the main figure to S5 at
the author's direction.

---

## Figure 5 — draft legend

**Figure 5. The nasal ageing axes replicate in independent cohorts.**

**(a)** Mucociliary axoneme score (module M02, mean z of log₂(CPM+1) over 18 expressed genes)
against age in an independent cohort of 120 nasal brushings, ages 47–89 (GSE336399, bulk
RNA-seq). Line, ordinary least squares; shaded band, 95% confidence interval. Inset: Spearman ρ
with its unadjusted p, and the age coefficient from a model adjusted for sex, cancer status,
smoking, RNA integrity (median TIN), log₁₀ library size and sequencing batch. M02 is
literature-curated and was fixed before this cohort was identified. The 13-module
Benjamini–Hochberg family is not significant (BH = 0.22); the analysis plan was committed before
the data were accessed.

**(b)** The 28 cilium-assembly genes derived in the primary cohort (CILIUM28), scored and
modelled identically in the same 120 brushings. This set was selected by correlation with a
senescence score and then filtered to genes highest in the elderly, so it is not independent of
the discovery cohort; (a) is.

**(c)** Specificity of (a) and (b). Grey, age coefficients from 500 expression-matched random
gene sets (25 expression bins, identical model); vertical lines, the observed coefficients.
Neither observed coefficient is approached by any random set, and **0 of 500 random sets reached
p < 0.05** against a nominal 5%, indicating no generic age-related inflation in this cohort.

**(d)** Interferon-stimulated gene score (module M05, UCell, donor means) in basal, secretory and
ciliated epithelium of SARS-CoV-2-negative children (n = 18) and adults (n = 23) (Loske et al.
2022). Horizontal bars, group medians. Below each cell type: Benjamini–Hochberg-adjusted p from
the primary family and the empirical p against 500 expression-matched random gene sets
(*, empirical p < 0.05).

**(e)** As (d), for MHC-I antigen processing (module M06). A second, independently curated module
separates children from adults in the same direction, and clears the matched-null control in all
three epithelial types.

---

**(f)** Mean programme score per epithelial family, pooled over 224 donor × cell-type units
from three single-cell nasal cohorts (Ziegler 2021 n = 105, Loske 2022 n = 103, Xu 2022 n = 16;
minimum 20 cells per donor × cell type). Solid box, the cell type the programme was curated for;
dashed orange, the one that is misassigned.

**(g)** Own cell type minus the mean of the other two, standardised mean difference, fixed-effect
meta-analysis across the three cohorts; bars, 95% CI. Five of six programmes are specific to
their curated cell type (BH < 0.05, sign unanimous in every cohort, all clearing a 500-draw
expression-matched null). **P6 entry-factor availability is significantly LOWER in ciliated
cells** (SMD −1.66, I² = 0) and also fails its matched null, so it is misassigned and loses its
cell-type label. Effect sizes of 2–4 are partly tautological — curated cell-type genes should be
enriched in that cell type — so this is a validation that the labels are correct, not a
discovery; its value is that the test could fail, and for P6 it did.

---

## Supplementary Figure S5 — draft legend

**Supplementary Figure S5. Nulls, measurement limits, and the pre-registered outcome tests.**

**(a)** M05 interferon score against age *within* older adults, in two independent cohorts:
GSE336399 (n = 120, ages 47–89) and GSE261277 baseline samples (n = 39, ages 51–80). Dashed line,
ordinary least squares; Spearman ρ and its unadjusted p inset. Neither cohort shows an age
gradient in interferon tone among adults, which is why the claim in Figure 4 is restricted to the
child–adult contrast. Neither cohort is well powered for a small effect, so this is an absence of
support, not evidence of absence.

**(b)** Detected genes per cell in ciliated cells, split by CDKN2A detection (CDKN2A⁻, blue;
CDKN2A⁺, orange); dashed vertical lines, medians. CDKN2A⁺ cells are sequenced more deeply (median
4,265 vs 2,888 detected genes; same sign in 181 of 187 donors), so p16 positivity called from
single-cell detection is substantially a sequencing-depth event. Data from an olfactory
brush-biopsy cohort, used here for measurement properties only and supporting no biological
claim.

**(c)** Fraction of each module's genes reaching ≥ 5% detection in ciliated cells across four
independent single-cell cohorts spanning two chemistries. Dashed line, the 60% measurability gate
applied throughout. M02, M05 and M11 clear the gate in every cohort; **M12 inflammatory SASP
fails in all four** (0–14% of genes detected), so secreted-cytokine SASP is below the detection
floor of droplet scRNA-seq in airway epithelium and carries no claim in this study.

**(d)** Outcome-anchored tests. Predictions, direction and multiplicity were fixed by a committed
plan before the data were accessed, and results are reported whichever way they fell.

**(e)** Baseline P5 mucin-barrier score and the odds of acquiring infection, per SD, in three
independent human challenge cohorts: RSV (GSE155237, n = 39), pneumococcal carriage (GSE261277,
n = 39) and A/H1N1 influenza (GSE253686 NOVA, n = 74). Horizontal bars, 95% confidence intervals;
diamond, fixed-effect pooled estimate; dashed line, no effect. All three point estimates are
protective and the cohorts are statistically homogeneous (I² = 4%). **The pooled estimate is shown
for completeness and is not claimed**: the NOVA cohort, which carries the most weight, fails its
specificity control (panel f, orange).

**(f)** Specificity control in the NOVA cohort. Distribution of the log odds ratio for 500
expression-matched random gene sets fitted with the identical model; orange line, the observed P5
estimate. 11.4% of random gene sets are at least as extreme (empirical p = 0.114), and 16.4% reach
p < 0.05 against a nominal 5%, so baseline bulk nasal expression in this cohort predicts shedding
non-specifically. The same control returned 0 of 500 in both earlier cohorts. Under the locked
plan (`docs/PREREG_nova_p5.md` §3) this voids the NOVA test irrespective of its p-value.

---

## What each panel is allowed to claim

| Panel | Claim it supports | Claim it does **not** support |
|---|---|---|
| 5a | axoneme transcription rises with age across 47–89, by an unbiased module | significance within the pre-specified family (BH = 0.22); anything about childhood |
| 5b | the manuscript's own gene set behaves the same way in a new cohort | independence — this set was selected in the discovery cohort |
| 5c | the rise is specific to these gene sets, not generic | that the effect is large (ρ ≈ 0.17) |
| 5d, 5e | higher interferon tone and MHC-I in children's epithelium, module-specific | a lifespan trajectory; the contrast tested is child vs adult |
| S5a | no interferon gradient within older adults, in two cohorts | that interferon tone is flat across the whole lifespan |
| S5b | p16⁺ calls track sequencing depth | any biological conclusion — the cohort is withdrawn |
| S5c | SASP cytokines cannot be measured in droplet scRNA-seq of nasal epithelium | that SASP is biologically absent |
| S5d | the axes were tested against real outcomes and reported honestly | that any axis predicts an outcome |
| S5e | baseline mucin-barrier score trends protective in three challenge cohorts and never reverses | **that P5 predicts infection** — the pooled OR is shown but voided by S5f |
| S5f | in the largest cohort the association is not specific to P5 | that the two smaller cohorts are also non-specific — there the control returned 0/500 |

**The sentence the discussion turns on.** CILIUM28 (panel b) is biased toward this result by
construction. **M02 (panel a) is not** — literature-curated, fixed before this cohort was
identified, never filtered on age or outcome — and it reproduces CILIUM28 in direction, magnitude
and specificity. The replication claim rests on panel a, with panel b as corroboration. Writing
it the other way round inverts the evidence.

**Panel 5a's BH = 0.22 is not optional.** It is the multiplicity of the test the figure leads
with, and the annotation stays. Presenting +0.0234 at p = 0.017 without it would misreport the
pre-specified primary analysis.

---

## Source data per panel

| Panel | File |
|---|---|
| 5a, 5b | `results/gse336399/module_scores.csv`, `results/gse336399/primary_age_effects.csv`, `results/gse336399/matched_null_control.csv` |
| 5c | `results/figure5/gse336399_null_draws.csv.gz` |
| 5d, 5e | `results/loske2022_v2/donor_celltype_scores.csv`, `results/loske2022_v2/age_effects_with_null.csv` |
| S5a | `results/gse336399/module_scores.csv`, `results/pneumo_gse261277/module_scores_all_samples.csv` |
| S5b | `results/crownlands_figrep/ciliated_cdkn2a.csv.gz` |
| S5c | `results/own_cohort/module_measurability.csv`, `results/ziegler2021_v2/module_coverage.csv`, `results/loske2022_v2/module_coverage.csv`, `results/crownlands/measurability_gate.csv` |
| S5d | `results/pneumo_gse261277/PREREG_RESULTS.md`, `results/laiv_gambia/go_vs_outcomes.csv`, `results/nova/primary.csv` |
| S5e | `results/nova/forest_three_cohort.csv` |
| S5f | `results/nova/matched_null_draws.csv`, `results/nova/primary.csv` |
| 5f, 5g | `results/celltype_spec/donor_celltype_scores.csv`, `results/celltype_spec/meta.csv` |

`validation/figure5_data.py` regenerates the 500 matched-null draws behind panel 5c. It re-runs
the control exactly as originally reported — same module order, seed, bins and model — and asserts
that the observed coefficient, empirical p and fraction of random sets at p < 0.05 for M02 and
CILIUM28 match `matched_null_control.csv` to machine precision. They do.

---

## Known gaps

1. **No functional or protein-level corroboration of the ciliary claim.** Every panel in Figure 5
   is transcript abundance measured by sequencing — the same measurement whose confounds S5b is
   about. Nothing here distinguishes "more ciliary programme per cell" from "more ciliated
   cells". Until that experiment exists the text must say *axoneme transcription increases with
   age*, never *ciliated cells increase* or anything about ciliary function. See
   `FIGURE5_PLAN.md` for the three experiments that would close it.
2. **GSE253686 (NOVA) was dropped from the table**, not completed. 75 participants with baseline
   SD −1 nasal arrays are downloaded; age and shedding outcome sit behind a PMC download gate. If
   obtained, the plan must be locked *before* they are opened, and the row returns.
3. **The "Chua 2020" citation for the SASP floor has no result file.** `RESULTS.md` and
   `FIGURE_WORDING.md` both say M12 fails in "4 of 4" cohorts including Chua 2020, but no Chua
   output exists in either tree. Panel S5c uses the four cohorts that **do** have verifiable
   detection tables: the primary cohort, Ziegler 2021, Loske 2022 and the olfactory Flex cohort.
   The text must be reconciled with that before submission.
4. **Figure 5 now shows only results that replicated.** That is the author's decision, and S5
   carries the rest. The main text must still state the two nulls and point to S5, or the figure
   set reads as selective reporting.
