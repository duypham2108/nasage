# External replication of the manuscript's composition and senescence claims

Status: **LOCKED** by the commit that adds this file. Later changes are dated deviations in the
results file, never edits here.

Companion to `docs/PREREG_crownlands_age.md` (commit `0d897a4`), which covers the module age
gradients. This plan covers the claims of Figures 1 and 3, and one secondary claim from
Figure 4, in the same cohort.

## 0. Why this plan needs to exist, and its one structural weakness

The claims below were derived from the primary cohort and written into the manuscript outline
before the Crownlands atlas was found. Testing them here is therefore a genuine external
replication of pre-existing claims, not a search.

**The weakness, stated plainly:** I chose *which* of the manuscript's claims to test after
seeing this cohort's cell-type composition and gene coverage. I did not see any outcome — no
basal fraction by age, no CDKN2A+ fraction, no axoneme contrast — but the selection of testable
claims was informed by what the cohort contains. The honest description is **a pre-specified
test of a selected subset of claims**, and any write-up says which claims were left untested and
why (section 6).

## 1. Claims under test, with their locked directions

Directions are taken from the manuscript's own figures. They are predictions, so the primary
tests are one-sided.

| ID | Source | Claim | Direction |
|---|---|---|---|
| **H1** | Fig 1c | Basal cells decrease with age | basal fraction **falls** with age |
| **H2** | Fig 3 | CDKN2A+ ciliated cells accumulate with age | CDKN2A+ fraction of ciliated cells **rises** with age |
| **H3** | Fig 3 | CDKN2A+ ciliated cells show reduced cilium movement | M02 axoneme **lower** in CDKN2A+ than CDKN2A- ciliated cells, within donor |
| **H4** | Fig 3 | CDKN2A+ ciliated cells show increased senescence programmes | M11 senescence arrest **higher** in CDKN2A+ than CDKN2A- ciliated cells, within donor |
| **H5** | Fig 3 | Basal loss and p16+ ciliated rise are linked | across donors, basal fraction and CDKN2A+ ciliated fraction are **negatively** correlated after adjusting for age |
| **S1** | Fig 4 | A distinct IFN-signalling programme is high again in the elderly | M05 **rises** from 40-64 to 65+ |

S1 is **secondary**, for two reasons. The childhood half of Figure 4 cannot be tested here at
all (no donor under 18), and GSE261277 gave M05 rho = -0.16, p = 0.32 across ages 51-80, which
points the other way. S1 is reported with that conflicting evidence beside it, one-sided for the
predicted rise, with the two-sided p also given.

## 2. Data, and which table each test may use

- **Composition tests (H1, H2, H5 numerator/denominator)** use the **full obs table** of all
  4,028,275 cells, healthy donors only. The capped subsample must not be used for any fraction
  that compares *across* cell types: capping at 500 cells per donor x cell type destroys the
  between-type proportions by construction.
- **Within-cell-type fractions** (CDKN2A+ among ciliated) may use the subsample, because the cap
  draws at random *within* each donor x cell type, so a within-type proportion stays unbiased.
  They are nonetheless computed on the full ciliated set where feasible, and both are reported.
- **Expression contrasts (H3, H4, S1)** use the fixed subsample
  `subsample_cap500_seed0.h5ad` from `PREREG_crownlands_age.md` section 4, unchanged.

Donors: the 192 with `disease == "normal"`. Note that 5 of them have no parseable age in
`development_stage`, so age-based tests use **187**; this is already a recorded deviation of the
companion plan and is restated in the results.

## 3. Definitions, fixed before analysis

**CDKN2A+** = at least **1** CDKN2A count in that cell. Sensitivity analysis at **>= 2** counts;
if the verdict differs between the two the claim is reported as `threshold_dependent` and
carries no weight. CDKN2A is detected in 20.6% of ciliated cells here (step-3 gate table), so
both thresholds are usable.

**Basal fraction** = respiratory basal cells / all **respiratory epithelial** cells of that
donor, where the denominator is the fixed list in `PREREG_crownlands_age.md` section 2 plus
goblet, suprabasal, club, deuterosomal, ionocyte and tuft cells. Olfactory populations and
immune cells are **excluded from both numerator and denominator**, because this cohort samples
the olfactory cleft and its olfactory and immune content is a property of the sampling site, not
of the donor's airway. Suprabasal cells are **not** counted as basal; a sensitivity analysis
adds them.

**Per-donor minimum:** a donor enters a test only with >= 200 cells in the relevant denominator,
and >= 20 cells on each side of a within-donor contrast. Counts are reported.

## 4. Models

Unit of analysis is the **donor** throughout. Never the cell: a cell-level test would set n to
the cell count and every p-value would be meaningless.

- **H1, H2:** one-sided Spearman of the donor fraction against age, plus a confirmatory
  logit-transformed linear model `logit(fraction) ~ age + sex + log(median genes per cell)`.
  Both reported.
- **H3, H4:** per donor, the difference in donor-level module score between CDKN2A+ and CDKN2A-
  ciliated cells, tested across donors with a one-sided Wilcoxon signed-rank test. Scoring is
  the package's normal path (UCell on raw counts).
- **H5:** one-sided Spearman of the two donor fractions, plus partial correlation adjusting for
  age, and a linear model with age and sex. The adjusted estimate is the primary one: without
  adjusting for age, H5 is guaranteed by H1 and H2 and tests nothing.
- **S1:** the age model of `PREREG_crownlands_age.md` section 5, restricted to the 40-64 versus
  65+ contrast.

**Multiplicity:** BH across the **five** primary hypotheses H1-H5. Smallest attainable BH is
p_min x 5, so BH < 0.05 needs p < 0.01, which 187 donors can reach (section 5). S1 is outside
this family and is corrected within the companion plan's family.

### 4.1 The depth confound, which is the main threat to H2

Whether a cell is called CDKN2A+ depends on how deeply it was sequenced. If older donors'
samples were sequenced deeper, H2 would be confirmed by an artefact. Three measures, all fixed
now:

1. Median detected genes per cell enters the confirmatory model as a covariate.
2. **Depth-matched sensitivity:** counts are downsampled per cell to a common total (the 10th
   percentile of per-cell totals among ciliated cells) and H2 is recomputed. If the verdict
   changes, H2 is reported as `depth_dependent` and carries no claim.
3. The correlation between donor age and median genes per cell is reported whatever it is.

## 5. Power, declared in advance

Simulated, 4,000 replicates.

Spearman at n = 192:

| true rho | alpha 0.05 | alpha 0.01 |
|---|---|---|
| 0.15 | 0.51 | 0.28 |
| 0.20 | 0.76 | 0.54 |
| 0.25 | 0.92 | 0.78 |
| 0.30 | 0.98 | 0.93 |

Within-donor paired contrast at n = 187, one-sided:

| true dz | alpha 0.05 | alpha 0.01 |
|---|---|---|
| 0.2 | 0.85 | 0.63 |
| 0.3 | 0.99 | 0.96 |
| 0.4 | 1.00 | 1.00 |

**Interpretation rule, fixed now:** the correlation tests are well powered for rho >= 0.25 and
underpowered below rho ~ 0.15; the paired tests are well powered for dz >= 0.3. A null is
reported with the relevant row of these tables, never as evidence of absence.

## 6. Claims deliberately not tested here, and why

Recorded so the selection in section 0 is visible.

| Claim | Why not tested |
|---|---|
| Fig 1d, the Pasta transcriptomic clock | trained on other platforms; this cohort is probe-based 10x Flex with 18,127 genes, so a poor correlation could not be distinguished from failure to transfer |
| Fig 2, Basal 1 / Basal 2 and pseudotime | no basal subsets annotated here; mapping them would be a derivation, and cross-platform pseudotime is not comparable |
| Fig 4, the childhood interferon peak | no donor under 18 |
| Fig 3, SASP | M12 fails the measurability gate in all three cell types, as it does in 4 other cohorts |

## 7. Failure conditions, stated in advance

| Condition | Consequence |
|---|---|
| A hypothesis fails to reach BH < 0.05 | reported as a failure to replicate, with the power table; no re-specification, no switch of denominator or threshold |
| H2 verdict changes under depth matching | `depth_dependent`; no claim |
| A verdict changes between the 1-count and 2-count CDKN2A definitions | `threshold_dependent`; no claim |
| H5 holds unadjusted but not after adjusting for age | reported as **not** supporting a link beyond shared ageing; this is the expected null and must be stated as such |
| Basal fraction differs from the primary cohort in absolute terms | **not** reported as a discrepancy: dissociation protocol and sampling site set absolute fractions, so only the age gradient within each cohort is comparable |

## 8. Order of work

1. Lock this plan by commit.
2. Composition tables from the full obs (H1, H2, H5), with per-donor counts reported first.
3. Within-donor expression contrasts (H3, H4).
4. Depth and threshold sensitivity analyses.
5. S1, last, after the primary family is complete.
