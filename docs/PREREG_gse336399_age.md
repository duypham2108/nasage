# Age gradient of the NasAge axes in older-adult nasal brushings (GSE336399)

Status: **LOCKED** by the commit that adds this file. Later changes are dated deviations in the
results file, never edits here.

## 0. Why this cohort, and why its companion was dropped

Crownlands is withdrawn (olfactory brush biopsy; age effects inseparable from per-cell detection
depth). This plan replaces it with a **true nasal brushing** cohort, accepting a narrower age
range in exchange for the right tissue and a different failure mode (bulk composition rather
than single-cell depth).

**GSE309130 was examined and is not used.** Its metadata, seen before any expression value:
every participant is SARS-CoV-2 infected, with no uninfected control; only 42 subjects; and the
asymptomatic group contains **no participant over 39**, so symptom status, strain
(omicron/delta/pre-delta) and vaccination are all collinear with age. An interferon-axis age
gradient cannot be separated from acute infection there. Recorded as examined-and-rejected so
the choice is visible.

## 1. Cohort

GEO GSE336399, "Characterization of the cancer-associated field of injury in the nasal
epithelium in never-smokers". **120 nasal brushings, one per subject**, bulk RNA-seq
(GPL30173), 57,736 features.

Seen before locking (GEO metadata only, no expression values):

| | |
|---|---|
| age | 47-89; bins 47-54: 12, 55-64: 30, 65-74: 49, 75+: 29 |
| sex | 64 female, 56 male |
| cancer status | Benign 40, Cancer 52, Metastatic 27, not collected 1 |
| smoking | NEVER 50, FORMER 60, CURRENT 10 |
| sequencing batch | 10 batches of exactly 12 - a balanced design |
| median TIN (RNA integrity) | mean 52.3, range 24.9-57.5; **uncorrelated with age** (rho = 0.065, p = 0.48) |

Linkage: expression columns match `!Sample_title` for **120 of 120** samples. No severed
identifiers are involved.

### 1.1 Declared weaknesses

1. **No young reference.** This tests the gradient *within* older adults, like GSE261277. It
   cannot address childhood, and therefore cannot test Figure 4's bimodality.
2. **Cancer status is confounded with age** (mean age 70.1 in Cancer vs 65.2 in Benign). It
   enters every model as a covariate, and a benign-only sensitivity analysis is pre-specified.
3. **Bulk, so composition is confounded.** An apparent state change may be a shift in cell-type
   proportions. This cannot be resolved here and is stated, not worked around.
4. Smoking history is present in all three cancer groups and enters as a covariate.

## 2. Primary analysis

- **Scoring.** Per module, mean over module genes of z-scored log2(CPM+1), z across the 120
  analysed samples. Genes with mean CPM < 1 across the cohort are dropped and reported. UCell on
  raw counts is a **sensitivity**, not a second primary test; a module whose verdict differs
  between the two is `scoring_dependent` and carries no claim.
- **Model.** `score ~ age + sex + cancer_status + smoking + median_TIN + log10(total counts) +
  sequencing_batch`, OLS, with `age` continuous. Reported alongside a Spearman correlation of
  score against age, unadjusted.
- **Direction: none predicted. Two-sided**, for every module. The gene sets were fixed long
  before this cohort was found.
- **Family and multiplicity.** All modules passing the expression gate, BH across them. The
  family size and the smallest attainable BH are reported **before** any p-value, as in the
  Crownlands plan. If BH < 0.05 is unattainable, no verdict is issued.

### 2.1 Specificity control

**500 expression-matched random gene sets** per module, 25 bins, same model. The empirical p is
the fraction of random sets with an age effect at least as extreme. A module that does not clear
the matched null is `not_specific` and carries no claim, whatever its BH.

This control is not optional: in the withdrawn Crownlands analysis 11.8-33.4% of matched random
sets tracked age, and it disqualified five BH-significant results.

## 3. Named secondary: the manuscript's 28 cilium genes

The 28 cilium-assembly genes shown in the manuscript's Figure 3 dot plot (SEPTIN2, DNAH3,
SEPTIN7, ARMC2, CFAP61, TTC29, SEPTIN9, CFAP206, B9D1, NEK10, CC2D2A, LRGUK, DNAH6, TTC12,
KIF27, BBS9, CFAP298, DMD, IFT57, ARMC9, DRC1, ARL3, CCDC39, CFAP47, BBS5, PIBF1, ROPN1L,
PARVA) are tested as a named gene set, same model, same matched-null control, **two-sided**.

Independence, fixed now: that list was derived in the primary cohort by correlation with a
senescence score **and then filtered to genes highest in the elderly**. The filter cannot
transfer to an independent cohort, so a positive result here is informative; a null is ambiguous
between "not real" and "selection effect". This is stated whichever way it falls.

## 4. Power, declared in advance

Spearman, 4,000 simulations:

| true rho | n = 120, alpha 0.05 | n = 120, alpha 0.0036 (BH over 14) | n = 40 (benign only) |
|---|---|---|---|
| 0.15 | 0.36 | - | 0.15 |
| 0.20 | 0.56 | - | 0.22 |
| 0.25 | 0.76 | - | 0.32 |
| 0.30 | 0.89 | 0.60 | 0.43 |
| 0.35 | - | 0.80 | - |
| 0.40 | 0.99 | 0.92 | 0.70 |

**Interpretation rule:** the full cohort is well powered for rho >= 0.3 after multiplicity and
underpowered below rho ~ 0.2. The **benign-only sensitivity is underpowered throughout** and is
read for direction and effect size only, never for significance. A null is reported with the
relevant row, never as absence of effect.

## 5. Sensitivity analyses, fixed now

1. **Benign only** (n = 40), same model without the cancer term. Direction and effect size only.
2. **Never-smokers only** (n = 50), same model without the smoking term.
3. UCell in place of mean-z.
4. Age as a 4-level bin in place of continuous, to expose non-monotonic patterns.

## 6. Failure conditions

| Condition | Consequence |
|---|---|
| Smallest attainable BH > 0.05 | no verdict; reported as unattainable |
| No module reaches BH < 0.05 | reported as a null with the power table; no re-specification |
| A module passes BH but not the matched null | `not_specific`; no claim |
| A verdict differs between mean-z and UCell | `scoring_dependent`; no claim |
| The age effect disappears when cancer status is added | reported as confounded by cancer status, not as an age effect |
| A module's genes are mostly below CPM 1 | reported via the coverage flag, not re-specified |

## 7. What this can and cannot claim

| Can claim | Cannot claim |
|---|---|
| "Axis X changes across ages 47-89 in nasal brushings, beyond expression-matched random sets" | Anything about children, adolescents or young adults |
| "The manuscript's cilium gene set does / does not change with age in an independent nasal cohort" | That a change is epithelial cell state rather than composition |
| Direction and effect size for a future, better-powered test | That any axis predicts a vaccine or infection outcome - no outcome exists here |

## 8. Order of work

1. Lock this plan by commit.
2. Expression gate and family size, reported before any p-value.
3. Primary models and the matched-null control.
4. The named cilium secondary.
5. Sensitivity analyses.
