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
