# Mucosal vaccine readiness: epithelial cell-type programs and their outcome test

Status: **LOCKED** by the commit that adds this file. Later changes are dated deviations in the
results file, never edits here.

This is the **second** attempt at a derived mucosal module. The first
(`results/setA_derivation/RESULT_attempt1_FAILED.md`) failed at the discovery stage because the
locked plan paired a per-gene t-test with BH < 0.05 over ~15,000 genes at n = 20 pairs, where the
best attainable BH was ≈ 0.06. That failure is the reason for three choices below: a restricted
gene universe, a hard cap on family size, and error control placed in the replication step rather
than in discovery.

---

## 0. What this can and cannot be

**It cannot be a fitted index.** Fitting weights needs a training cohort with both a nasal
transcriptome and a vaccine outcome. No such cohort is available: the one mucosal-vaccine dataset
with open outcomes (LAIV Gambia, GSE230494) has its participant link severed, and this project
will not reconstruct it or contact the custodians.

**So this is a pre-registered screen, not an index.** Programs are defined with no access to any
outcome, then tested. An *index* is only justified if programs survive, and its definition is
fixed in advance in §6 so that the decision is not made after seeing which ones did.

Naming discipline for the manuscript: until §6 is met, these are **candidate epithelial programs**,
never "the mucosal vaccine readiness index".

---

## 1. Derivation cohort and cell types

Primary cohort: `integrated_Healthy_covid_airways_ethnicity_senescent_SoupX_April252024_raw_minimal.h5ad`
— 148,978 cells, 38,269 genes, 52 healthy donors, nose only, raw counts, no COVID cases.

Cell types, from `Cell_type_annotation_level2`: **Basal, Secretory, Ciliated** only. This matches
the primary family used throughout the project; immune populations are not used, and goblet,
suprabasal and the rest are not used.

**No outcome data of any kind is loaded during derivation.** The derivation script does not import
or read any file under `results/rsv_gse155237/`, `results/pneumo_gse261277/` or
`results/thwaites2023_laiv/`.

## 2. Program derivation (unsupervised)

Per cell type, in this order:

1. **Donor filter.** Donors with >= 50 cells of that type.
2. **Gene universe, declared a priori.** Genes detected in >= 10% of that cell type's cells, and
   present in >= 80% of retained donors at that threshold. A program whose genes cannot be
   detected cannot transfer to bulk, and the first attempt's lesson is that the universe must be
   cut on grounds independent of any result.
3. **Normalisation.** log1p of CP10K, per cell.
4. **Within-donor gene–gene correlation.** Spearman correlation computed *within each donor* and
   averaged across donors, weighted equally per donor. Correlation computed across all donors at
   once would be dominated by donor and kit-version structure; this is the step that makes the
   programs about within-tissue co-variation rather than between-donor batch.
5. **Program detection.** kNN gene graph (k = 15) on the averaged correlation, Leiden clustering,
   resolution 1.0, seed 0.
6. **Program filter.** Keep programs with **15 to 150 genes**.
7. **Coherence.** For each surviving program, coherence = mean within-donor correlation among its
   genes, averaged over donors.

### 2.1 Family size cap, fixed now at 8

Power (§4) is the binding constraint, not biology. **At most 8 programs enter the test**,
selected by the pre-declared rule: rank all surviving programs by coherence within each cell
type, then take programs round-robin across Basal, Secretory, Ciliated (highest-coherence first
in each) until 8 are held, so each cell type is represented. If fewer than 8 survive §2.6, all of
them enter and the family is smaller.

This rule uses no outcome information and is fixed before any program is seen.

### 2.2 Redundancy with the existing curated modules

A derived program whose gene overlap with any curated module M01-M13 exceeds **Jaccard 0.5** is
**dropped before the cap is applied** and recorded as redundant. The point of this exercise is to
add axes, not to rediscover M05 and claim it as new.

## 3. Transfer to bulk, and the composition problem

Both validation cohorts are **bulk**. A program derived inside one cell type, scored in bulk, is
a mixture of that program's activity and that cell type's abundance. This is the single largest
threat to the interpretation and it is handled by design, not in discussion.

- **Score.** Mean z of log2(CPM+1) across program genes present at CPM >= 1 in that cohort. A
  program with fewer than 60% of its genes available is **not tested** in that cohort and is
  reported as untransferable.
- **Composition proxies**, declared now, built from marker genes that are **forbidden** from
  entering any program (enforced in code — a program containing one of these genes has it
  removed before scoring):
  - ciliated: FOXJ1, TPPP3, PIFO, SNTN, CAPS
  - basal: KRT5, KRT14, TP63, DLK2, MIR205HG
  - secretory: SCGB1A1, BPIFA1, MSMB, WFDC2, SLPI
- **Primary model adjusts for the program's own cell-type proxy.** An association that survives
  this is a state effect; one that does not is abundance. The unadjusted model is a
  **sensitivity**, reported alongside, never as the primary.

## 4. Validation design and power

Two independent human challenge cohorts, both baseline (pre-challenge) nasal transcriptomes with
a binary outcome:

| cohort | n | outcome | platform |
|---|---|---|---|
| GSE155237 (RSV challenge, Habibi 2020) | 39 (23 infected / 16 uninfected) | infection | nasal curettage RNA-seq |
| GSE261277 (pneumococcal challenge) | 39 (17 carriage+ / 22 carriage−) | carriage acquisition | nasal RNA-seq |

**Primary test: fixed-effect inverse-variance meta-analysis** of the per-SD log odds ratio across
the two cohorts, per program, with BH across the program family. Per-cohort logistic models
adjust for the composition proxy and for log10 library size.

Power, logistic regression, 3,000 simulations, two-sided:

| true OR per SD | RSV alone | pneumococcal alone | pooled (n = 78) | pooled, BH over 8 |
|---|---|---|---|---|
| 1.5 | 0.19 | 0.19 | 0.41 | 0.11 |
| 1.8 | 0.37 | 0.39 | 0.71 | 0.36 |
| 2.0 | 0.50 | 0.51 | 0.85 | 0.53 |
| 2.5 | 0.75 | 0.76 | 0.98 | 0.83 |
| 3.0 | 0.89 | 0.88 | 1.00 | 0.96 |

**Interpretation rule, fixed now.** The pooled test is adequately powered only for OR >= 2.5 per
SD after multiplicity, and is underpowered below OR ≈ 2. A null result is reported with this
table and is **not** evidence that the programs are uninformative. Single-cohort results are read
for direction and effect size only.

**Direction: none predicted, two-sided**, for every program. These are unsupervised programs; no
prior exists.

### 4.1 Heterogeneity

RSV is a virus and pneumococcus a bacterium; a program may be relevant to one and not the other.
Cochran's Q and I² are reported per program. A program with I² > 75% is reported as
**pathogen-specific**, with both cohort estimates shown, and does **not** count towards §6.

### 4.2 Specificity control

**500 expression-matched random gene sets** per program per cohort, 25 expression bins, same
model, meta-analysed identically. The empirical p is the fraction of random sets with a pooled
|log OR| at least as large. A program that does not clear this is `not_specific` and carries no
claim whatever its BH. This control is not optional: in one cohort analysed earlier in this
project 11.8-33.4% of matched random sets tracked the real variable.

## 5. Exploratory endpoint: response to live vaccine

GSE230494 (Thwaites 2023, LAIV) has 20 participants with paired baseline and 72 h samples linked
by study ID. The outcome table is unusable (severed link), but the **within-participant response
magnitude** is not.

- Endpoint: change in M05 ISG score from baseline to 72 h, within participant.
- Test: Spearman correlation of baseline program score against that change, n = 20.
- **Circularity guard, fixed now:** any program overlapping M05 at Jaccard > 0.1 is excluded from
  this endpoint.
- Power: r = 0.5 gives 0.62, r = 0.3 gives 0.25 at n = 20. **Exploratory only, no BH, no verdict**,
  reported for direction and effect size.

Also recorded now so it is not re-attempted: influenza reads in these libraries are **zero** in
all 13 samples screened (`results/thwaites2023_laiv/viral_read_precheck.csv`), so a shedding proxy
cannot be recovered from the RNA-seq itself.

## 6. When an "index" may be declared

All of the following, or the word index is not used:

1. At least one program reaches **BH < 0.05** in the pooled primary test;
2. **and** clears the matched-null control (§4.2);
3. **and** survives composition adjustment (§3), i.e. the primary model, not only the sensitivity;
4. **and** has I² <= 75% (§4.1), i.e. is not pathogen-specific.

If two or more programs meet all four in the same direction, the index is their **unweighted mean
z**. No weights are fitted — there is no cohort to fit them on, and fitting on the validation
cohorts would destroy the only independent evidence available.

## 7. Failure conditions

| Condition | Consequence |
|---|---|
| Fewer than 3 programs survive §2.6 | reported as a derivation failure; no test is run |
| A program is untransferable in a cohort (§3) | excluded from that cohort, reported, not re-specified |
| No program reaches BH < 0.05 | reported as a null with the power table; **no re-specification, no threshold relaxation** |
| A program passes BH but fails the matched null | `not_specific`; no claim |
| A program passes unadjusted but fails composition-adjusted | reported as **composition, not state** |
| I² > 75% | `pathogen_specific`; shown but excluded from §6 |
| Leiden returns a single giant program | reported; resolution is **not** tuned after the fact |

## 8. What this cannot claim, whatever the result

- That any program predicts response to an actual **vaccine**. Neither validation endpoint is a
  vaccine outcome; both are pathogen challenge. The LAIV endpoint in §5 is a host response
  magnitude, not protection, and is exploratory.
- That a program is epithelial-cell-type-specific **in the validation cohorts**. It was derived in
  one cell type; in bulk it is scored across whole tissue.
- Anything causal.

## 9. Order of work

1. Lock this plan by commit.
2. Derive programs (§2), with no outcome file accessible. Report program count, sizes, coherence,
   redundancy drops and the final family of <= 8 **before** any outcome is touched.
3. Transfer and score in both cohorts; report transferability.
4. Primary pooled test, matched-null control, heterogeneity.
5. Exploratory LAIV endpoint.
6. Apply §6 and state whether the word "index" is earned.
