# Replication of the P5 mucin-barrier result in the NOVA influenza challenge

Status: **LOCKED** by the commit that adds this file. Later changes are dated deviations in the
results file, never edits here.

This is a **replication of a specific prior result with a pre-specified direction**, not a
screen. The result being replicated is in `results/mvri_outcome/PREREG_RESULTS.md`: P5 mucin
barrier, pooled OR 0.190 per SD, BH 0.009, across RSV challenge and pneumococcal carriage.

**Nothing in this plan was written with knowledge of any expression–outcome relationship in
NOVA.** Only the outcome marginals and the ID linkage below were inspected, which is what the
power table requires.

---

## 1. Cohort and linkage, seen before locking

GEO GSE253686, Walters et al. 2025 *Sci Transl Med*
[doi:10.1126/scitranslmed.adt1452](https://doi.org/10.1126/scitranslmed.adt1452). Nasal
epithelium microarray (GPL6480), 297 samples, 75 participants, four timepoints:
**SD −1 (pre-challenge baseline)**, SD1, SD3, SD7. Adults challenged with A/H1N1; half had
received licensed quadrivalent inactivated vaccine a mean of 64 days earlier.

Outcome recovered from supplementary data file S1 (`S1_parsed.csv`):

| | |
|---|---|
| participants in S1 | 74 |
| **linkage to GEO nasal subjects** | **74 of 75 matched by Subject ID**; NOVA06 is in GEO only and is excluded |
| shedding, binary | **53 shedders / 21 non-shedders** |
| multimodal class | NS 21, 1DS 19, MDS 34 |
| study group | Vax 37 / Unvax 37 |
| sex | 40 F / 34 M |
| age | 19–55, median 32 |

No identifiers were reconstructed: S1 carries the NOVA subject IDs openly and GEO sample titles
carry the same IDs.

**S1 resolves a GEO defect.** Three participants (NOVA01, NOVA03, NOVA10) carry conflicting
vaccination labels across their own timepoints in GEO, flagged in `FIGURE5_PLAN.md`. S1's
`Study Group` is taken as authoritative for all 74, and the conflict is recorded rather than the
participants dropped.

## 2. Primary test — one test, one direction

**Hypothesis.** Higher baseline (SD −1) P5 mucin-barrier score is associated with **lower odds of
viral shedding** after A/H1N1 challenge.

- **Direction is locked as protective (OR < 1), one-sided**, because it was established in two
  prior cohorts before these data were opened. The two-sided p is reported alongside for
  transparency and the verdict rests on the one-sided test.
- **Outcome:** `Shedding (BINARY)` from S1, 53 / 21.
- **Model:** logistic regression of shedding on the per-SD P5 score, adjusted for **age, sex,
  vaccination arm, and the secretory composition proxy**.
- **One test. No multiplicity correction is applied to it**, which is the point of
  pre-specifying a single directional replication.
- **Scoring:** mean z of the P5 genes present on the array, z across the 74 analysed baseline
  samples. Programs and the composition proxy are fixed in `validation/mvri_programs.py`.
  If fewer than 60% of P5 genes are on GPL6480, the test is **not run** and that is the result.

### 2.1 Power, declared before any p-value

Logistic regression, 4,000 simulations, α = 0.05:

| true OR per SD | NOVA alone (53/21) | MDS vs NS (34/21) | all three cohorts pooled |
|---|---|---|---|
| 0.40 | 0.93 | 0.88 | 1.00 |
| 0.50 | 0.74 | 0.66 | 0.99 |
| 0.60 | 0.48 | 0.42 | 0.85 |
| 0.70 | 0.26 | 0.22 | 0.55 |

**Interpretation rule.** The prior estimate (OR 0.19) is almost certainly inflated by winner's
curse; the honest expectation is OR 0.4–0.6. NOVA alone is well powered at 0.4, adequately at
0.5, and **underpowered at 0.6 and above**. A null is therefore reported with this table and is
**not** read as refutation unless the confidence interval excludes OR 0.6.

### 2.2 What counts as replication

| Observed | Verdict |
|---|---|
| one-sided p < 0.05 **and** point estimate OR < 1 | **replicated** |
| OR < 1, p ≥ 0.05, CI includes 0.5 | **underpowered, direction consistent** |
| OR < 1, CI excludes 0.6 but p ≥ 0.05 | **not replicated at the plausible effect size** |
| OR > 1 | **contradicted**; reported as prominently as a success would be |

## 3. Specificity control

500 expression-matched random gene sets, 25 expression bins, identical model, as everywhere else
in this project. Empirical p is the fraction with a log OR at least as extreme. P5 failing this
control voids the test whatever the p-value.

A **label permutation** of shedding status (2,000 draws) is also run, as it was for the original
result, to confirm the pipeline is not anticonservative on this cohort.

## 4. Secondary and exploratory, clearly subordinate

1. **The other five programs** (P1, P3, P4, P6, P7), two-sided, BH across the five. P2 is
   excluded a priori for Jaccard 0.824 with M01 — the redundancy rule that was applied late last
   time and is applied up front here.
2. **Three-cohort pooled meta-analysis** of P5: RSV + pneumococcal + NOVA, fixed effect, with
   Cochran's Q and I². Reported as the summary estimate, with the caveat that NOVA is the only
   one of the three that is influenza.
3. **Dose–response:** the three-level shedding class (NS / 1DS / MDS) by ordinal logistic, and
   MDS vs NS directly. Exploratory; the binary test above is primary.
4. **Effect modification by vaccination arm:** a P5 × vaccination-status interaction.
   **Exploratory and underpowered at 37 per arm** — it is reported for direction only and is the
   nearest thing in this project to a vaccine-relevant readout. It cannot support a claim.
5. **Post-baseline timepoints are not used.** Only SD −1 enters any test here.

## 5. Failure conditions

| Condition | Consequence |
|---|---|
| < 60% of P5 genes on the array | test not run; reported as untransferable |
| P5 passes but fails the matched null | `not_specific`; no claim |
| Significant unadjusted, null after composition adjustment | reported as composition, not programme |
| Primary null | reported with the power table; **no switching to the ordinal or MDS-vs-NS test to rescue it** |
| Fewer than 70 participants link | reported; the power table is recomputed and stated |

## 6. What this still cannot claim

- **Not a vaccine outcome.** The endpoint is viral shedding after challenge. Half the cohort was
  vaccinated, which permits the exploratory interaction in §4.4, but no test here measures
  response to a vaccine. The phrase "mucosal vaccine readiness index" remains unlicensed
  whatever the result.
- **Not cell-type specific**, since nasal tissue is profiled in bulk.
- **Not causal.**
- Nothing about children: NOVA is 19–55.

## 7. Order of work

1. Lock this plan by commit.
2. Map probes to symbols; report P5 gene coverage on GPL6480 **before** any model is fitted.
3. Primary one-sided test, with the matched null and the label permutation.
4. Secondary and exploratory analyses, in the order listed.
5. Apply §2.2 and state the verdict.
