# Every known weak point in the MVRI line of work

Date: 2026-11-25. Compiled from all results and audits in this project. Ordered by how much
each one threatens a conclusion, not by how awkward it is to state.

**Severity key.** **FATAL** — blocks a claim outright. **MAJOR** — changes what may be said.
**MINOR** — bookkeeping or presentation.

---

## A. The one that governs everything

| # | Weakness | Sev. | Fixable here? |
|---|---|---|---|
| **A1** | **No vaccine outcome exists in any reachable cohort.** Four attempts: Liverpool LAIV-vs-TIV (null), Gambia IgA (participant link severed), NOVA vaccination interaction (p = 0.41), Thwaites engagement (null, positive control passed). | **FATAL** | **No.** Not an analysis problem — the data do not exist publicly (`docs/COHORT_SEARCH_mvri.md`). |

Nothing below changes A1. Any amount of method work still leaves "mucosal vaccine readiness"
unsupported, because readiness requires a vaccine outcome and there is none.

## B. Statistical weaknesses — the real battleground

| # | Weakness | Sev. | Evidence |
|---|---|---|---|
| **B1** | **A dominant axis makes gene-set tests non-specific.** PC1 carries 41–47% of module variance, 55–69% of the five lifespan axes. In NOVA, **16.4% of expression-matched random gene sets reach p < 0.05** against a nominal 5%. In the lifespan analysis, **7 of 11** tested cells pass BH but fail the matched null. | **FATAL** for specificity claims | `MODULE_DESIGN_AUDIT.md` §1.2, `results/nova/` §3, `results/lifespan/` §3 |
| **B2** | **Removing that axis destroys most of the signal.** After PC1 residualisation, NOVA P5 goes OR 0.471 → 0.663 (p 0.057 → 0.37); M06 dies in basal; several effects **flip sign**. | **MAJOR** | `results/method_improve/` |
| **B3** | **The declared family is not the effective family.** 14 modules behave as 7–10 independent tests; 5 axes as 3. BH over the declared count is wrong in both directions. | **MAJOR** | `MODULE_DESIGN_AUDIT.md` §1.2 |
| **B4** | **The redundancy rule catches none of the redundancy.** 24 module pairs exceed \|ρ\| = 0.70; **all 24 have Jaccard 0.0**. | **MAJOR** | `MODULE_DESIGN_AUDIT.md` §1.1 |
| **B5** | **Matched null and model p disagree systematically.** P7 clears 0/500 yet p = 0.34; P1/P7 clear in LAIV yet BH ≥ 0.53. The primary statistic is self-contained, the control is competitive: they ask different questions. | **MAJOR** | `mvri_outcome/` §5.5, `laiv_engagement/` §3 |
| **B6** | **Winner's curse.** The original OR 0.19 shrank to 0.471 on replication, exactly as predicted. | **MAJOR**, now handled | `results/nova/` §2 |
| **B7** | **P5 in NOVA is significant only after composition adjustment** (unadjusted OR 0.848, p = 0.28) — a different estimand from the two cohorts where it held both ways. | **MAJOR** | `results/nova/` §4 |
| **B8** | **Project-wide multiplicity is untracked.** Dozens of tests across many plans; BH is applied within families, never across the programme of work. | **MAJOR** | — |
| **B9** | **Composition shifts with age** (secretory −0.0077/yr, p < 0.0001), so any pooled-tissue score inherits a false trend. | MAJOR, controlled | `results/lifespan/` §5 |

## C. Data and cohort weaknesses

| # | Weakness | Sev. |
|---|---|---|
| **C1** | **Small n everywhere**: RSV 39, pneumococcal 39, NOVA 74, LAIV 20 complete pairs (39 participants). Power 0.27 at the observed LAIV effects. | **MAJOR** |
| **C2** | **No neonates or infants.** Youngest donor anywhere is 4 y, so the childhood peak in the original design is untestable and the span is 4–77, not "lifespan". | **MAJOR** |
| **C3** | **Liverpool non-independence** — GSE117580 and GSE164649 are one cohort family, not two replications. | **MAJOR** |
| **C4** | **Yoshida 2022 is excluded** (the user's own primary data), removing the largest nasal atlas from validation. | **MAJOR**, unavoidable |
| **C5** | **Crownlands withdrawn** for biological claims over the CDKN2A depth confound. | MAJOR |
| **C6** | **Specificity rests on two cohorts.** Xu 2022 contributes 4–7 donors and almost no weight. | MAJOR |
| **C7** | **The lifespan replication declared in the plan was never run** (Ziegler, GSE336399). | **MAJOR**, outstanding |
| **C8** | **NOVA's GEO matrices are ratio-to-own-baseline**, so baseline is unrecoverable without the raw arrays. | MINOR, solved |

## D. Measurement weaknesses

| # | Weakness | Sev. |
|---|---|---|
| **D1** | **M12 SASP is unmeasurable** in droplet scRNA-seq — 0–17% of genes detected, failing in all four cohorts and all three cell types. | **FATAL** for that axis |
| **D2** | **M04 innate sensing fails the gate in 2 of 3 cell types** (0.50 basal, 0.57 ciliated). | **MAJOR** |
| **D3** | **P4 IgA transport is only 49.7% epithelial** — JCHAIN (plasma cell), TNFSF13B (BAFF), IL6. Valid in single cell, reports immune infiltration in bulk. | **MAJOR** |
| **D4** | **Scoring is inconsistent**: rank-based UCell in single cell, mean-z of log2 CPM in bulk. The depth confound lives in the mean-z path. | **MAJOR** |
| **D5** | **M03 and M03b are single genes** called modules. | MINOR |
| **D6** | **Gene-set sizes range 1–22** with no size-matched null, so empirical p values are not comparable across sets. | MAJOR |

## E. Design and inference weaknesses

| # | Weakness | Sev. |
|---|---|---|
| **E1** | **P6 is misassigned** — significantly *lower* in ciliated cells (SMD −1.66), and fails its matched null. | **MAJOR**, now reported |
| **E2** | **Specificity is only within epithelium.** A programme can be secretory-specific among epithelial cells and higher still in an immune type (exactly P4). | MAJOR |
| **E3** | **Cell-type specificity is partly tautological** — curated cell-type genes should be enriched in that cell type. SMDs of 2–4 overstate the achievement. | MAJOR |
| **E4** | **The M06 trajectory is not independent** of the existing child-vs-adult panels in the same cohort. | MAJOR |
| **E5** | **Circularity risk**: re-testing the same cohorts with improved methods cannot confirm anything. | **MAJOR**, structural |
| **E6** | **No causal evidence**, and no functional or protein-level corroboration of any claim. | MAJOR |
| **E7** | **Programmes were curated, not learned.** The unsupervised derivation failed (MT/rRNA artefacts) and was abandoned. | MINOR |

---

## What is actually attackable

Of the 30 weaknesses, most are data limits that no method fixes. **Three are genuinely
addressable by better analysis**, and they are the only ones worth more compute:

1. **B1/B2 — the dominant axis.** The deepest problem. Everything that fails the matched null
   fails because of it. Addressing it properly means either removing it (tried; destroys signal)
   or **designing gene sets that are orthogonal to it by construction**.
2. **B3/B4/B5 — the testing framework.** Use a competitive test as the primary statistic rather
   than bolting one on as a control, with the effective family size.
3. **D4/D6 — scoring and null calibration.** Make scoring consistent and size-matched.

**C1 (power) and A1 (no vaccine outcome) are not analysis problems and will not yield.** Honest
accounting: even if 1–3 all succeed, the result is a better-specified association with infection
in ≤ 152 challenge participants, and still no vaccine claim.

## The one redesign worth attempting

**Axis-orthogonal programme redesign.** Select, within each programme, the genes that are
(a) specific to the owning cell type and (b) **least loaded on the global PC1**, then re-test.
This attacks B1 at its root instead of correcting for it afterwards.

It is admissible only under strict conditions, because it is a redesign:

- **Gene selection must use no outcome variable** — only single-cell expression structure, as
  the epithelial filter did.
- **Selection and testing must use different cohorts.** Selecting on Ziegler/Loske structure and
  testing on RSV/pneumococcal/NOVA outcomes keeps them separate.
- **The original programmes stay frozen**; the redesign is a new, separately named set.
- **It may fail**, and a redesigned set that loses the association is evidence the association
  was the global axis all along — which is a result worth having.


---

# Progress, 2026-11-25

| # | Weakness | Status |
|---|---|---|
| **B1/B2** | dominant axis | **characterised, not solved.** Residualisation destroys signal; P5 cannot be made axis-orthogonal (7 genes of 11, below the floor). Reported as a limit of the claim. |
| **B3** | effective vs declared family | **fixed in code.** `effective_tests()` reports Li & Ji and Galwey beside every BH. |
| **B4** | Jaccard rule | **fixed.** `redundancy_report()` applies the score-correlation rule; Jaccard kept as a secondary check. |
| **B5** | matched null vs model p | **characterised.** Approach B shows the primary statistic is anticonservative under inter-gene correlation; P7's VIF is 12.5. |
| **B6** | winner's curse | **handled.** Predicted 0.4–0.6, observed 0.471. |
| **B8** | project-wide multiplicity | **FIXED.** `test_registry.py`: 832 p-values, 140 primary. P5's outcome results do not survive project-wide BH; both figures now reported. |
| **C7** | lifespan replication unrun | **FIXED.** Run in Ziegler and GSE336399. Adult trajectories are flat in all three cohorts; the decline is a childhood-to-adult contrast. |
| **D1/D2** | M12, M04 unmeasurable | **confirmed independently** in Ziegler (0.29, 0.40), a fourth cohort. |
| **D6** | size-matched nulls | **fixed in code** (`validation/matched_null.py`). |
| multiverse dev. 2 | RSV/pneumo curve unrun | **FIXED.** 448 specifications; P5 is 100% sign-stable in both, robust in pneumococcal. |
| **A1** | no vaccine outcome | **unchanged and unfixable here.** |
| **C1** | n = 39/39/74 | **unchanged and unfixable here.** |

**Still open:** the lifespan arm of the multiverse grid (`PREREG_multiverse` §1); D4 scoring
consistency between bulk and single cell; E6 functional corroboration; and the v2 panel exists
but has never been tested on a cohort not used to build it.
