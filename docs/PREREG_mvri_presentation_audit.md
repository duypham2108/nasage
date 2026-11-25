# Pre-registered plan: re-analysis of the MVRI slide deck in the primary cohort

Locked 2026-11-25, before any axis score is compared across age bins. Data:
`integrated_Healthy_covid_airways_ethnicity_senescent_SoupX_April252024_raw_minimal.h5ad`
(148,978 cells, 38,269 genes, **52 donors**, nose only, healthy only).

This re-analyses the five axes of the "Mucosal vaccine readiness" deck (slides 9–14) and the
weighted composite index built from them. It is an audit of an existing analysis, not a new
claim, so every test below is declared here with its verdict rule.

## 1. What is being checked, and why

The deck reports five axis scores as **per-cell boxplots across seven age bins** and combines
them into a weighted index that peaks in childhood. Four properties of that analysis are
testable without any new data:

| # | Concern | Test |
|---|---|---|
| **A** | Are the axes measurable at all in droplet data? | Detection gate, before any model |
| **B** | Per-cell boxplots treat cells as independent samples of age; the unit is the donor (n = 52 across 7 bins) | Same test on cells and on donor means, side by side |
| **C** | Is any age difference specific to the curated genes? | 500 expression- and size-matched random gene sets |
| **D** | Does the childhood peak follow from the data or from the chosen weights? | Weight sensitivity over the full sign-respecting simplex |

## 2. Axes, transcribed verbatim from the deck

Gene lists and cell subsets are taken from slides 9–13 without modification, so that this tests
the deck's analysis rather than a variant of it. `CGAS` and `TMEM173` are resolved to the
symbols present in this object; no other substitution is permitted. Weights are from slide 14:
mucus −0.20, antigen +0.25, innate +0.25, IgA +0.20, cytokine +0.10.

**Each axis is scored only in the cell subset the deck assigns to it.**

## 3. Pre-specified tests and verdicts

**A — Measurability.** An axis passes if ≥ 60% of its genes are detected in ≥ 5% of cells of its
own subset. *Declared in advance: cytokine tone is expected to fail, because the SASP/cytokine
family has failed this gate in four other cohorts. An axis that fails the gate is reported as
unmeasured, and any flat trajectory for it carries no information.*

**B — Unit of analysis.** Kruskal–Wallis across the seven age bins, computed twice per axis:
once over cells, once over donor means. Report both p-values and the ratio.
*Declared in advance: the per-cell p will be many orders of magnitude smaller. The donor-level
p is the valid one.* An axis is called **age-varying** only if the donor-level BH-adjusted p
across the five axes is < 0.05.

**C — Specificity.** For each axis, 500 random gene sets matched on expression decile and set
size, scored and tested identically at donor level. Empirical p is the fraction of null sets
whose Kruskal–Wallis statistic is at least as large as observed. An axis is **specific** only if
empirical p < 0.05. The fraction of random sets reaching donor-level p < 0.05 is reported as the
calibration of the test in this cohort.

**D — Weight sensitivity.** The composite is recomputed over 10,000 random weight vectors drawn
so that the mucus term is negative and the other four positive — that is, respecting every sign
the deck asserts — with magnitudes drawn uniformly and normalised. Report the distribution of
*which age bin attains the maximum*. The childhood peak is called **robust to weighting** only
if a childhood bin (Young child or Child) peaks in ≥ 90% of draws.

**E — Dimensionality.** Li & Ji effective number of tests among the five donor-level axis
scores, with PC1 variance fraction.

## 4. What this cannot do

It cannot test readiness. No vaccine outcome exists in this cohort, so the index is unvalidated
regardless of how A–E resolve. A favourable result here means only that the age pattern is real
and specific, not that it predicts response to vaccination.

It also cannot test a continuous age trend: this object carries `Age_bin` only, with no numeric
age, so every test is across seven ordered bins.

## 5. Declared expectations

Recorded so the result cannot be presented as having confirmed a prior that was never stated:

1. Antigen presentation and IgA will remain age-varying at donor level.
2. Innate sensing (PRR) and cytokine tone will not.
3. The childhood peak will **not** be robust to weighting.
4. The five axes will behave as fewer than five effective tests.

If 1–4 all hold, the recommendation is to report the surviving axes individually and withdraw
the composite index.

## 6. Outputs

`D:\research\nasage_validation\results\mvri_audit\`. Script:
`validation/mvri_presentation_audit.py`.
