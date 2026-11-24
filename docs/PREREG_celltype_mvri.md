# Are the seven epithelial programmes actually cell-type-specific, and does a live
# intranasal vaccine engage them?

Status: **LOCKED** by the commit that adds this file. Later changes are dated deviations in the
results file, never edits here.

Programmes and composition proxies are fixed in `validation/mvri_programs.py`, locked at `d391a86`
with per-gene references in `docs/MVRI_PROGRAMS.md`. **No gene is added, removed or reweighted
here.** P2 remains excluded a priori (Jaccard 0.824 with M01).

---

## 0. Why this exists, and what it is not

The seven programmes are *defined* as cell-type-specific — P1/P2 basal, P3/P4/P5 secretory,
P6/P7 ciliated. **That assignment has never been tested.** Every outcome test so far (RSV,
pneumococcal, NOVA) used **bulk** tissue, where cell-type specificity is unmeasurable by
construction. A programme could carry its whole signal in the wrong cell type and every bulk
result would look identical.

**This plan does not test vaccine readiness and cannot.** Per-participant mucosal IgA and
shedding outcomes are not public for any cohort in reach (`docs/COHORT_SEARCH_mvri.md`), the
Liverpool between-arm engagement test was null (`results/mvri_engagement/`), and the NOVA
vaccination interaction was flat (OR 0.586, p = 0.41). The phrase **"mucosal vaccine readiness
index" remains unlicensed and is forbidden in any figure title or results sentence produced under
this plan**, whatever the outcome. §5 sets what may be claimed.

Two questions are asked, in order:

- **A. Specificity.** Is each programme highest in the cell type it was curated for?
- **B. Engagement.** Does a live attenuated intranasal influenza vaccine move these programmes
  within participant?

## 1. Prior literature fixing the directions

Wilson et al. 2023, *Viruses* [doi:10.3390/v15091956](https://doi.org/10.3390/v15091956) — in
primary human nasal epithelial cultures, **ciliated cells are the dominant cell type infected by
influenza B**, with a Victoria preference for MUC5AC⁺ cells and a Yamagata preference for basal
cells, and a strong interferon response by 48 h. LAIV is a live, replicating influenza virus
delivered to this epithelium.

This fixes the engagement directions in §3 **before** the data are opened. It is not a
cell-type-resolved vaccine study in humans and is used only to set direction.

## 2. Analysis A — cell-type specificity

**Cohorts** (single-cell, nasal, cell types annotated, none containing Yoshida 2022 donors):

| Cohort | object | donors |
|---|---|---|
| Ziegler 2021 nasopharyngeal | `ziegler2021_np.h5ad` | 58 |
| Loske 2022 controls | `loske2022_controls.h5ad` | — |
| Xu 2022 turbinate | `xu2022_turbinate.h5ad` | 7 |

Crownlands is **excluded**: the cohort is withdrawn for biological claims over the CDKN2A
sequencing-depth confound (Supp. Fig. S5b) and is used in this project for measurement properties
only.

**Prediction, locked.** For each programme, the mean score in its **own** cell type exceeds the
mean across the other two epithelial families. Own-cell-type assignment is taken from
`mvri_programs.py` and is not revisited:

| Programme | own cell type |
|---|---|
| P1 basal progenitor/repair, P2* | Basal |
| P3 antimicrobial, P4 IgA transport, P5 mucin barrier | Secretory |
| P6 entry-factor availability, P7 ciliary beat machinery | Ciliated |

\* P2 is scored and reported but **excluded from the family** and from BH, as locked previously.

**Unit of analysis is the donor, not the cell.** Scores are computed per cell (UCell), averaged
to donor × cell type, and only donor × cell-type cells with **≥ 20 cells** are kept. Cells, not
donors, are the pseudo-replicates that would otherwise inflate significance.

**Model.** Per cohort, paired across cell types within donor: Δ = own − mean(other two). Tested
with a two-sided Wilcoxon signed-rank across donors, then combined across the three cohorts by
fixed-effect inverse-variance meta-analysis of the standardised mean difference, with Q and I².
**BH across the family of 6.**

**Depth.** `log_genes` per donor × cell type is reported alongside. If a programme's specificity
tracks depth rather than identity, that is stated and the programme is marked `depth_confounded`.

**Measurability gate** as everywhere: a programme needs ≥ 60% of genes at ≥ 5% detection in that
cell type, or it is not tested there and that is the result.

### 2.1 What counts as a pass

| Observed | Verdict |
|---|---|
| Δ > 0, BH < 0.05, consistent sign in all three cohorts | **specific** |
| Δ > 0, BH < 0.05, sign flips in a cohort | **specific but heterogeneous**; I² reported |
| BH ≥ 0.05 | **not demonstrated**; the cell-type label on that programme is withdrawn |
| Δ < 0 | **misassigned**; reported as prominently as a pass, and the programme is renamed |

## 3. Analysis B — engagement by a live intranasal vaccine

**Cohort.** Thwaites et al. 2023, GEO GSE230494, *Nat Commun*
[doi:10.1038/s41467-023-43842-7](https://doi.org/10.1038/s41467-023-43842-7). Nasal curettage
bulk RNA-seq, **20 participants with paired baseline and day-3 samples** after LAIV.

This cohort already has a **passed manipulation check** in this project
(`results/thwaites2023_laiv/PREREG_RESULTS.md`): M05 interferon rose in 17/20 participants,
paired t p = 0.008, d_z = 0.66. **M05 is therefore the positive control for Analysis B**, not a
finding. If M05 does not rise in the re-run, Analysis B is void and reported as such.

**Design.** Paired within-participant change, baseline → day 3, for each programme. Two-sided,
**BH across the family of 6**. n = 20 pairs.

**Locked directions** from §1, stated for interpretation but **not used to make the tests
one-sided** — the family is screened two-sided because the literature concerns wild-type
influenza B in culture, not LAIV in humans:

- P6 entry-factor availability and P7 ciliary beat machinery (ciliated, the infected compartment):
  expected to **change**, direction not fixed.
- P5 mucin barrier (secretory, the barrier compartment): expected to **rise** if the epithelium
  mounts a barrier response.
- P1 basal repair: expected to **rise** if the vaccine causes epithelial damage.

**Composition adjustment is the cell-type-specific part of this analysis.** Thwaites is bulk, so a
programme can appear to move purely because the cell mixture shifted. Every paired test is
repeated with the change in that programme's composition proxy (`PROXY` in `mvri_programs.py`)
as a covariate. A programme that is significant unadjusted and null adjusted is reported as
**composition, not programme**, exactly as in the earlier plans.

**PTPRC (CD45) change** is carried as a covariate in a sensitivity model, because immune-cell
infiltration after vaccination would move every epithelial programme downward by dilution alone.

### 3.1 Power, declared before any p-value

Paired t, n = 20, α = 0.05 two-sided, 4,000 simulations:

| true d_z | power |
|---|---|
| 0.80 | 0.93 |
| 0.66 (the observed M05 effect) | 0.77 |
| 0.50 | 0.50 |
| 0.35 | 0.27 |

**n = 20 is adequate only for a large paired effect.** A null for any programme is reported with
this table and is **not** evidence of no engagement.

## 4. Controls

1. **Expression-matched random gene sets**, 500 draws, 25 bins, identical model, both analyses.
   A programme failing its matched null is `not_specific` whatever its p-value — the rule that
   voided the NOVA result and is applied here from the start.
2. **Label permutation**, 2,000 draws (cell-type labels in A, timepoint labels in B).
3. **Leave-one-gene-out** for any programme that passes, to confirm it is not one gene.

## 5. What may be claimed, and the naming rule

A programme may be described as **cell-type-specific** only if it passes A. It may be described as
**vaccine-engaged** only if it passes B *after* composition adjustment, with M05 having risen.

**The composite may be given a name only if ≥ 3 programmes pass A and ≥ 2 of those also pass B.**
Even then the permitted name is a **cell-type-resolved epithelial barrier score**, never
"vaccine readiness index": engagement is not readiness, and no outcome is measured anywhere in
this plan. If fewer pass, the figure reports programmes individually and no composite is formed.

**Forbidden in any output of this plan:** "readiness", "predicts vaccine response", "vaccine
efficacy", and any statement that these programmes forecast an immune outcome.

## 6. What this cannot claim regardless of result

- **No outcome of any kind** is measured — not IgA, not shedding, not seroconversion.
- Analysis B is **bulk**, so its cell-type reading is inferred from Analysis A, not measured.
  This limitation is stated in the figure legend, not only here.
- Analysis A is **cross-sectional** and says nothing about ageing. Age is not tested in this plan.
- LAIV in Thwaites is given to adults; nothing here extends to children or the elderly.

## 7. Order of work

1. Lock this plan by commit.
2. Analysis A: measurability gate and cell counts reported **before** any model.
3. Analysis A tests, matched null, permutation.
4. Analysis B: M05 positive control **first**; if it fails, stop and report.
5. Analysis B tests, adjusted and unadjusted, matched null, permutation.
6. Apply §2.1 and §5; state the verdicts; build panels only from what passed.
