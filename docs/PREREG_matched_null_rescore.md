# Pre-registered plan: matched nulls under the corrected rank scorer

Locked 2026-11-25, **before** any null is drawn under the corrected scorer. Closes the item left
explicitly open in commit `7de3a52`:

> "What would make this more convincing: re-running the matched nulls under the corrected scorer,
> which has **not** been done. Until that is done, the corrected curves show stability of the
> estimate, not specificity."

## 1. Why this is needed

Every matched null in this project was computed with **mean-z** scoring
(`validation/mvri_outcome.py` §4.2, `validation/nova_p5.py`, `validation/lifespan_trajectory.py`).
The corrected specification curves report results from **four** scorers, including the repaired
UCell rank statistic. A specificity control computed under one scorer does not license a claim
made under another: the null distribution of a rank statistic has a different shape from the null
of a mean-z statistic, so the empirical p is not transferable.

The curves therefore currently demonstrate **stability of the estimate across analytic choices**,
which is not the same as **specificity against expression-matched gene sets**. This plan supplies
the missing half.

## 2. Design

For each of three outcome cohorts (NOVA influenza shedding, RSV infection, pneumococcal carriage)
and the lifespan analysis (Loske, three epithelial families):

- **Observed statistic**: the same model as the locked plan, with the programme scored by scorer
  *s* ∈ {mean_z, median_z, rank_ucell, pc1_of_set}.
- **Null**: 500 expression-matched random gene sets per programme per cohort (25 bins, same size
  as the programme), **scored with the same scorer *s*** and run through the identical model.
- **Empirical p**: the two-sided fraction of null draws whose statistic is at least as extreme as
  the observed one, exactly as `matched_null.empirical_p` defines it.
- **Null calibration**: the fraction of null draws reaching model p < 0.05. Under a well-behaved
  test this is 0.05; the mean-z runs gave 0.164 in NOVA, which is the B1 dominant-axis problem.

The primary programme is **P5 (mucin barrier)**. All other programmes are reported, uncorrected
for the extra scorers, as context.

## 3. Verdicts, fixed before any number is seen

| Outcome | Reading |
|---|---|
| P5 empirical p < 0.05 under **all four** scorers in a cohort | Specificity holds in that cohort and is not a scoring artefact |
| P5 empirical p < 0.05 under the **corrected rank scorer only** | **Not** a rescue. Report as scorer-dependent; the mean-z null remains the locked result |
| P5 empirical p ≥ 0.05 under the corrected scorer where mean-z passed | The locked specificity claim weakens; report the weaker figure |
| Null fraction p<0.05 stays ≫ 0.05 under every scorer | The dominant axis is scorer-invariant — a stronger statement of B1 than currently made |

**No scorer may be selected after the fact.** The headline for each cohort is the **worst**
empirical p across the four scorers, not the best. This is declared here precisely because the
rank fix already moved results in the favourable direction once.

## 4. What this cannot do

It cannot create a vaccine outcome (A1), cannot increase n (C1), and cannot make P5 orthogonal to
the dominant axis (B1/B2). A pass here means only: the P5–infection association is not reproduced
by random gene sets matched on expression and size, under the scorer used to report it.

## 5. Outputs

`D:\research\nasage_validation\results\matched_null_rescore\` — `rescore_nulls.csv`,
`PREREG_RESULTS.md`. Script: `validation/matched_null_rescore.py`.
