# Pre-specified vaccine-type directions for NasAge axes, and their first test (LAIV, GSE230494)

Status: **locked before analysis.** The git commit that adds this file is its timestamp. Any later
change is a labelled deviation (section 7); this text is not edited in place.

## 0. What was known when this was written (2026-11-24)

**Seen:**
- GEO metadata of GSE230494: 40 young adults given LAIV; bulk RNA-seq of inferior-turbinate
  curettage at baseline (up to 7 days before) and 72 h after vaccination; Salmon counts. GEO holds
  no outcome, age or sex fields.
- The first sentences of the GEO series summary.
- Our descriptive cross-check of 26 published GO-term scores in the Gambian paediatric LAIV
  cohort: a different cohort, not NasAge modules. There, baseline innate/defence regulation was
  associated with less vaccine shedding.

**Not seen:**
- any GSE230494 expression values;
- any outcome data;
- the results sections or figures of the linked paper (Thwaites et al. 2023, PMID 38052824).

Mechanistic reasoning per module is from `docs/MODULE_EVIDENCE.md` ("How to read it" entries).

## 1. Directions by vaccine type (hypotheses)

Signs mean: `+` = a higher baseline module score predicts a stronger vaccine take or response;
`−` = predicts a weaker one; `0` = no directional prediction; `±` = opposing mechanisms, so no
direction. Only the LAIV column is tested here (section 2); the other columns are recorded now so
that later tests are pre-specified too.

| Axis | Live attenuated, replicating (LAIV) | Replication-deficient viral vector, intranasal | Protein/subunit + mucosal adjuvant, intranasal | Mechanism |
|---|---|---|---|---|
| M05 ISG (interferon tone) | **− replication, − mucosal IgA** | − (weak) | + (weak; endogenous adjuvant) | antiviral state restricts replicating vaccine; IFN is also adjuvant-like |
| M06 MHC-I | same as M05 (not independent; IFN-inducible) | same as M05 | 0 | see M05 |
| M04 innate sensing | − replication (weak) | ± | + (TLR/STING-adjuvant responsiveness) | sensing triggers antiviral response; also adjuvant responsiveness |
| M03 PIGR | **+ mucosal IgA**; 0 replication | + mucosal IgA | + mucosal IgA | transport of secretory IgA to the lumen |
| M02b mucins | **− replication** | − | − (traps antigen) | mucus and sialylated decoys trap virus/antigen |
| M02 mucociliary | ± | − (clearance) | − (shorter residence) | clears vaccine; ciliated cells are also influenza targets |
| M01 barrier | ± | − (access) | − (antigen access to APCs) | junctions restrict subepithelial antigen access |
| M11 senescence | 0 | 0 | 0 | no mechanistic basis for a direction |
| M12 SASP, M07/M08/M09/M10/M13 | not tested here (bulk tissue cannot isolate immune cell states; M12 poorly measured in swab epithelium) | | | |

## 2. Primary hypotheses tested in GSE230494 (LAIV, young adults)

All use **baseline** (pre-vaccination) samples only.

| ID | Module | Endpoint | Predicted sign |
|---|---|---|---|
| H1 | M05 ISG | E1 vaccine replication | − |
| H2 | M05 ISG | E2 mucosal IgA response | − |
| H3 | M03 PIGR | E2 mucosal IgA response | + |
| H4 | M02b mucins | E1 vaccine replication | − |

**Secondary (reported, never promoted to primary):**
- M06 → E1/E2, predicted −. This is the M05 signal again if r(M05, M06) > 0.7, so it is not
  counted as independent support.
- M04 → E1, predicted −.
- All modules → E3, serum response: same signs as E2, weaker expectation.
- M01 and M02 → E1/E2, two-sided with no prediction (descriptive only).

**Manipulation check (must pass before H1–H4 are interpreted):**
- M05 increases from baseline to 72 h (paired), and increases more in participants with E1 = 1
  than E1 = 0.
- If M05 does not rise after a live influenza vaccine, the scoring or data are not fit for
  purpose; H1–H4 are then reported but marked "uninterpretable".

## 3. Endpoints

Taken from the published study data. The first available definition in each list is used:
- **E1 vaccine replication:**
  1. any LAIV strain detected by PCR in nasal samples after vaccination (binary);
  2. otherwise peak or area-under-curve viral load (continuous, log).
- **E2 mucosal IgA response:**
  1. fold-rise in nasal (nasosorption or nasal-wash) vaccine-strain HA-specific IgA, baseline to
     the post-vaccination peak, log2 (continuous);
  2. otherwise the authors' IgA responder status (binary).

  If several strains are reported, use the mean of log2 fold-rises across strains.
- **E3 serum response:** HAI or serum HA-IgG fold-rise, log2, same rules as E2.

If a pre-specified endpoint is not available, the hypotheses that use it are recorded as **not
testable**. No other endpoint is substituted, except as a labelled deviation.

## 4. Scoring (bulk tissue)

- **Input:** Salmon counts summed to gene symbols, restricted to baseline samples with any read
  depth. Samples are excluded only if fewer than 5 million assigned reads, or (for the timepoint
  pair) a participant ID cannot be matched.
- **Primary score:** for each module, the mean across module genes of the z-scored
  log2(CPM + 1). Genes with mean CPM < 1 across baseline samples are dropped. Gene lists are
  NasAge 0.2.0 `gene_modules.json` as of this commit, with HLA genes kept.
- **Sensitivity scores:**
  1. UCell on counts (as the package does for single cells);
  2. M05/M06 scored without HLA genes (`nas.drop_hla`).
- **Module-specific check:** for each primary hypothesis, 1,000 random gene sets
  expression-matched to the module genes (25 expression bins, as in `nasage.controls`), scored
  and modelled identically. Empirical p = (1 + number of random sets whose test statistic is in
  the predicted direction and at least as extreme) / 1,001.

## 5. Models

One row per participant (baseline sample).
- **Binary endpoint:** logistic regression. **Continuous endpoint:** linear regression.
- Predictor: the module score, standardised (per SD).
- **Covariates:**
  - log10 assigned reads;
  - immune-content proxy, log2(CPM + 1) of PTPRC (bulk samples differ in immune content);
  - sex and age if the published data provide them;
  - baseline (pre-vaccination) antibody level on the same log scale, for E2 and E3 only (prior
    immunity both blocks replication and changes fold-rise).
- **Sensitivity:**
  1. no PTPRC covariate;
  2. no baseline-antibody covariate;
  3. UCell scores.

## 6. Decision rules

Multiple testing: Holm across H1–H4 (4 tests), using one-sided p in the predicted direction.

| Verdict | Rule |
|---|---|
| **supported** | predicted sign, Holm p < 0.05, **and** empirical p < 0.05 against random gene sets |
| **consistent** | predicted sign, one-sided nominal p < 0.05, but fails Holm or the random-set check |
| **not supported** | anything else in the predicted direction or null |
| **contradicted** | opposite sign with two-sided p < 0.05 |

Power (stated in advance): with about 35 participants, a linear model detects a standardised
slope of about 0.45 (80% power, one-sided α 0.0125). A logistic model detects an odds ratio per
SD of about 3 if E1 is roughly balanced. Moderate true effects can therefore be missed, so a
"not supported" verdict is not evidence of no effect.

**Independence caveat:**
- If the linked paper already reported a baseline transcript–outcome association that overlaps
  a hypothesis, that hypothesis is marked "not independent of the source publication". It is
  then a replication with a prespecified gene set, not a new discovery.
- This is checked only after the analysis.

## 7. Deviations

Any change after this commit (endpoint substitution, covariate unavailable, sample exclusions not
covered above, coding errors found later) goes into a dated "Deviations" list appended to the
results file, with the reason. Hypothesis signs and the decision rules are never changed.

## 8. Extensions (dated; sections 0-7 above are unchanged)

### 8.1 — 2026-11-24, RSV challenge (GSE155237)

The live-replicating arm of section 1 is tested in a challenge cohort, because LAIV outcomes
cannot be linked to the public transcriptomes. Signs are taken unchanged from section 1; the
full specification is the docstring of `validation/rsv_challenge_gse155237.py`, committed before
the expression data were downloaded. It is a **proxy** for the live-replicating arm: a challenge
virus, not a vaccine.

### 8.2 — 2026-11-24, pneumococcal carriage after mucosal vaccination

Datasets: GSE164649, GSE117580, GSE261277, GSE124949 (nasal transcriptome around experimental
pneumococcal challenge; the first two randomise LAIV against inactivated vaccine). Only GEO
sample metadata has been seen: no expression values. Endpoint: carriage acquisition (POS/NEG).

Predicted signs, recorded now so they cannot be chosen later:

| Axis | Carriage acquisition | Mechanism |
|---|---|---|
| M05 ISG (interferon tone) | **+** (higher tone, higher odds of carriage) | type I/III interferon impairs innate control of pneumococcal colonisation; this is the basis of the influenza/LAIV-pneumococcus interaction |
| M12 inflammatory SASP | **+** | neutrophilic mucosal inflammation accompanies impaired control |
| M01 barrier | **−** | an intact epithelial barrier resists colonisation |
| M04, M06, M02b, M03, M11 | 0 | no directional prediction |

Note this is the OPPOSITE sign to M05 in section 1: interferon tone is predicted to restrict a
replicating *virus* while impairing control of a *bacterium*. That contrast is itself a
prediction, and it is the clearest statement of why no single "mucosal readiness" score can
exist.

Independence must be assessed per dataset before any verdict: these source papers studied
exactly this contrast, so several of these are likely replications rather than new tests.

### 8.3 — 2026-11-24, GSE261277 analysis specification (LOCKED before expression values seen)

Section 8.2 fixed the predicted signs for pneumococcal carriage. This section fixes the model,
the test and the interpretation rules for one of the datasets it named, GSE261277.

**State of knowledge at lock time.** Seen: the GEO series record, the 115 GSM characteristics
blocks, the column *names* of `GSE261277_raw_expression_elderly.tsv.gz`, and the published
supplementary tables mmc1-mmc6. **Not seen: any expression value from this dataset.** The signs
in 8.2 were locked before even the metadata were available.

#### 8.3.1 Cohort and endpoint

Experimental human pneumococcal challenge (EHPC), serotype 6B strain BHN418, 80,000 CFU per
nostril, healthy adults aged 50-84. Nasal curettage bulk RNA-seq at baseline (pre-challenge),
D2 and D9. Endpoint: **carriage acquisition, POS vs NEG**, as deposited in the GSM
characteristics.

Primary analysis set: the **41 baseline samples** (24 NEG, 17 POS). Baseline is pre-inoculation,
so this is a prospective test of whether mucosal state predicts acquisition, not a response
measurement.

#### 8.3.2 Scoring

As in 8.1 and the module confirmation plan: mean over module genes of z-scored log2(CPM+1),
z across the 41 analysed baseline samples; genes with mean CPM < 1 across those samples are
dropped and the dropped genes are reported. UCell on raw counts is computed as a **sensitivity
analysis**, not a second primary test; a module whose verdict differs between the two is
reported as `scoring_dependent` and carries no claim (the rule already used for M12 in
GSE155237).

Library-size QC: samples below 5 M assigned counts are excluded, as in 8.1. If that removes any
sample the counts above are restated in the results file.

#### 8.3.3 Primary family: the three locked signs

One-sided tests, direction taken from 8.2 and not revisited:

| Module | Locked direction in carriage POS | Test |
|---|---|---|
| M05 ISG | higher (+) | one-sided Mann-Whitney, POS > NEG |
| M12 inflammatory SASP | higher (+) | one-sided Mann-Whitney, POS > NEG |
| M01 epithelial barrier | lower (-) | one-sided Mann-Whitney, POS < NEG |

BH across these **three** tests. Smallest attainable BH is p_min x 3, so BH < 0.05 is reachable
whenever p < 0.0167. This attainability was checked before locking; the failure of the two
previous plans (`setA_derivation`, `module_confirmation`) was caused by omitting exactly this
check.

Confirmatory model, reported alongside the rank test and not replacing it:
`score ~ carriage + age + sex + smoker + log10(total counts)`, OLS, one-sided p on the carriage
coefficient. Charlson index is **not** included: it is a composite of the comorbidity fields and
would be collinear with them. Both the unadjusted rank test and the adjusted model are reported
whatever they show.

#### 8.3.4 Independence: which of these is a new test

Determined from the source publication's own supplementary tables (mmc1-mmc4), which list the
22 Reactome pathways it tested on the Susceptible-vs-Protected contrast at -5, 2 and 9 days.

| Status | Modules | Basis |
|---|---|---|
| **Replication** (the paper tested an overlapping programme on the same contrast) | M05, M06, M07, M08, M12 | their list contains *Interferon Signaling*, *DDX58/IFIH1-mediated induction of interferon-alpha/beta*, *Class I MHC mediated antigen processing & presentation*, *MHC class II antigen presentation*, *Costimulation by the CD28 family*, *Signaling by Interleukins*, *Neutrophil degranulation* |
| **Independent** (absent from their list) | M01, M02, M02b, M03, M03b, M09, M10, M11, M13 | no epithelial barrier, ciliary, mucin, immunoglobulin-transport or senescence pathway appears |

Consequence, fixed now: **M01 is the only novel directional test in the primary family.** A
positive result for M05 or M12 is reported as *replication of the source paper's finding using
the NasAge modules*, never as independent support for the axes. This wording is fixed here so
it cannot be softened after seeing the result.

#### 8.3.5 Secondary and exploratory, reported but not claimed

- All remaining modules at baseline, **two-sided**, BH within that family separately. Exploratory.
- D2 and D9 carriage contrasts, two-sided, exploratory. These are post-exposure, so they measure
  response to colonisation, not readiness, and cannot support a predictive claim.
- Age association of each module across 51-80 at baseline, two-sided, exploratory. This cohort
  has no one under 50, so it tests the gradient *within* older adults only and is not a
  substitute for the cross-lifespan age analysis.

#### 8.3.6 Power, declared in advance

Simulated one-sided Mann-Whitney, 17 vs 24, 4,000 replicates:

| True effect | power at alpha 0.05 | power at alpha 0.0167 |
|---|---|---|
| d = 0.4 | 0.32 | 0.17 |
| d = 0.6 | 0.56 | 0.37 |
| d = 0.8 | 0.78 | 0.60 |
| d = 1.0 | 0.91 | 0.78 |
| d = 1.2 | 0.98 | 0.93 |

**Interpretation rule, fixed now:** a null result in this cohort is reported as *underpowered
for effects below d ~ 0.8*, never as evidence that an axis does not predict carriage. A
significant result is reported with its confidence interval, and the effect size is expected to
be inflated by selection at this sample size (winner's curse), which is stated in any write-up.

#### 8.3.7 Known limitations, fixed now so they are not discovered later

1. **Bulk, so composition is confounded.** A baseline inflammation signal may be immune-cell
   presence rather than epithelial state. This cannot be resolved in this dataset. The
   cell-type-specificity check of the module confirmation plan section 4 is run against Ziegler
   2021 and reported; a module whose genes concentrate > 50% in one population is flagged
   `composition_dominated` and not carried forward as a cell-state claim.
2. **No donor identifier.** Donors can be reconstructed only from identical covariate profiles;
   5 such groups contain two donors each. Paired within-donor analyses across timepoints are
   therefore **not performed**. The primary baseline test does not require pairing.
3. **Elderly arm only.** Only the 51-80 cohort is deposited; the younger comparison cohort shown
   in the paper's mmc2 is not in GEO, and the submitters state the raw data were lost. No
   age-group interaction can be tested.
4. **M12 is at the detection floor in single-cell data** (4 of 4 datasets; see
   `results/own_cohort/MEASURABILITY.md`). It is measurable in bulk, so it is tested here, but a
   positive M12 result in this bulk cohort **cannot be transferred** to the single-cell analysis.

#### 8.3.8 Failure conditions

| Condition | Consequence |
|---|---|
| No module in the primary family reaches BH < 0.05 | reported as a null, with the power table; no re-specification, no switch to D2/D9 as primary |
| M05 or M12 significant but M01 not | reported as replication only; the axes gain no independent support |
| A module's verdict differs between mean-z and UCell | `scoring_dependent`; no claim |
| > 50% of a module's genes concentrate in one cell type | `composition_dominated`; not transferred |
| Fewer than 41 baseline samples survive QC | counts restated; analysis proceeds; power table recomputed and reported |
