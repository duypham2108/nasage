# Four pre-specified attempts to improve on the negative results

Status: **LOCKED** by the commit that adds this file.

---

## 0. What this is, and the limit that cannot be argued away

Several results in this project are null or fail their specificity control. Each has a
**diagnosed cause**, and each cause has a principled fix. This plan declares the **complete list
of fixes to be tried, before any of them is run**, so that the exercise is a methods comparison
and not a search for a significant result.

**Three rules bind everything below.**

1. **All four approaches are reported, whichever way each falls.** The list is fixed here. No
   approach may be added after seeing a result, and none may be dropped.
2. **Nothing here can confirm anything.** Every approach re-analyses data already analysed. A
   positive result from a better method on the same cohort is **hypothesis-generating only**, and
   must be labelled so in every sentence that reports it. Confirmation requires a cohort not used
   here.
3. **Multiplicity across approaches is counted.** Four approaches × the declared families is the
   real test count. Any result is reported with a BH correction over the **union** of tests run
   under this plan, not within one approach.

If all four fail, that is the answer, and the nulls stand as they are.

## 1. Approach A — remove the dominant axis before testing

**Diagnosis.** `docs/MODULE_DESIGN_AUDIT.md` §1.2: PC1 carries 41–47% of module-score variance
(55–69% among the five lifespan axes). `results/nova/` §3: 16.4% of expression-matched random
gene sets reach p < 0.05 for shedding. `results/lifespan/` §3: 7 of 11 trajectory cells reach
BH < 0.05 but fail the matched null. These are one phenomenon — a global axis that almost any
gene set partially loads on.

**Fix.** Compute PC1 of the full score matrix (all programmes, within cell type / cohort) and
**residualise each programme score on it** before testing. A programme that still predicts after
the shared axis is removed is carrying something of its own.

**Applied to:** the NOVA shedding test, and the 11 lifespan trajectory cells.

**Control.** The matched null is recomputed **with the identical residualisation applied to each
random gene set**. This is the critical detail: residualising only the real programme and not the
null would manufacture specificity. If the matched null is not recomputed the same way, the
approach is void.

**Risk, stated now.** PC1 may itself be the biology. Removing it can destroy a real effect as
easily as a confound, and a *loss* of signal under Approach A is not evidence the signal was
spurious.

## 2. Approach B — a correlation-aware competitive test

**Diagnosis.** The project's primary statistics are self-contained (a model p-value per
programme) while its control is competitive (versus random gene sets). The two ask different
questions, which is why they disagree so often. Gene-set tests are known to be anticonservative
when genes within a set are correlated.

**Fix.** Use a **CAMERA-style variance-inflation correction** (Wu & Smyth 2012, *Nucleic Acids
Res* 40:e133, [doi:10.1093/nar/gks461](https://doi.org/10.1093/nar/gks461)): estimate the mean
pairwise inter-gene correlation within each programme and inflate the test variance by
`1 + (m − 1)·ρ̄`. This makes the primary statistic competitive by construction, aligning it with
the control rather than fighting it.

**Applied to:** NOVA shedding, the lifespan trajectories, and the LAIV engagement test.

**Expected direction, declared:** this should make results **more** conservative, not less. If an
existing null becomes significant under B, that is a red flag about the implementation, not a
discovery, and will be investigated before being reported as anything.

## 3. Approach C — use all LAIV samples, not only complete pairs

**Diagnosis.** `results/laiv_engagement/` used 20 complete baseline/day-3 pairs out of 59 QC-passing
samples; the locked power table gives 0.27 power at the observed effect sizes. The null may be
power, not absence.

**Fix.** A **linear mixed model across all 59 samples**, `score ~ timepoint + (1 | participant)`,
which uses the unpaired samples instead of discarding them.

**Positive control unchanged:** M05 interferon must still rise, or the approach is void for the
same reason as before.

**Honest ceiling, declared now:** this recovers at most the information in 19 extra singleton
samples. It is a modest power gain, not a new experiment, and will not rescue an effect of
d_z ≈ 0.1.

## 4. Approach D — meta-analyse the three vaccine-engagement tests

**Diagnosis.** Three separate LAIV engagement nulls (Thwaites paired, GSE117580, GSE164649), each
underpowered alone.

**Fix.** Fixed-effect inverse-variance meta-analysis of the per-programme standardised change
across all three, with Cochran's Q and I².

**The non-independence rule still applies.** `PREREG_mvri_engagement.md` §0.1 declares GSE117580
and GSE164649 to be **one Liverpool cohort family, not two independent replications**. They
therefore contribute **one** pooled estimate, not two, and the meta-analysis has **two
independent units**: Liverpool and Thwaites. Treating them as three would inflate the apparent
evidence, and that is not done.

## 5. What counts as success, fixed now

| Outcome | Reading |
|---|---|
| a null becomes significant **and** clears its matched null **and** survives BH over the union | **hypothesis for a new cohort**; never reported as a replication or a confirmed finding |
| a null stays null | the null stands, now with a better method behind it — a stronger negative than before |
| an existing positive **weakens** under A or B | reported as prominently as any gain; the earlier estimate was inflated |
| an approach disagrees with its own control | the approach is void, not the control |

## 6. What will not be done

- **No new gene sets, no gene added or removed, no reweighting.** The panels are fixed.
- **No change to any outcome definition, cohort inclusion, or covariate set.**
- **No approach invented after seeing results.** The four above are the complete list.
- **No dropping of a cohort or a cell type** that fails under a new method.
- **"Mucosal vaccine readiness" stays forbidden** regardless of outcome, per
  `PREREG_celltype_mvri.md` §5 — none of these approaches introduces a vaccine outcome.

## 7. Order of work

1. Lock this plan by commit.
2. Run A, B, C, D in that order. Record each result before starting the next.
3. BH over the union of all tests run under this plan.
4. Write every result into one file, including the ones that changed nothing.
