# Figure 5: external validation — panel plan

Status: proposal, 2026-11-24. Written against what has actually been produced, not against what
we hoped to produce. Every number cited here exists in `D:\research\nasage_validation\results\`
and every test was run from a plan committed before the data were touched.

Figures 1–4 are the discovery story in the primary cohort. **Figure 5 is the paper's defence.**
Its job is not to add a new claim — it is to show that the central claims survive outside the
cohort that generated them, and to show the reader exactly where they do not.

---

## The claim set Figure 5 has to answer

| From | Claim | External status |
|---|---|---|
| Fig 3 | ciliary/axoneme programme is elevated in elderly nasal epithelium | **replicated** (GSE336399 + ALI) |
| Fig 4 | interferon tone is high in childhood | **replicated** (Loske 2022) |
| Fig 4 | interferon tone rises again in the elderly (bimodality) | **not replicated**, two cohorts |
| Fig 3 | inflammatory SASP accompanies senescence | **untestable** — unmeasurable in 4/4 cohorts |
| Fig 1–2 | the axes predict mucosal vaccine / infection outcome | **one null, one descriptive** |

Four of these have a real answer. Figure 5 should present all five, including the two that went
against us, because the two negatives are what make the three positives credible.

---

## Panel a — functional validation in air–liquid interface culture

**Your own ALI data.** Proportion of α-tubulin+ (ciliated) cells, elderly vs young donors.

This panel carries the ciliary claim at the protein level. It is the only evidence in the paper
that is immune to every confound that dominates the rest of this work — sequencing depth, cell
composition inferred from transcript abundance, gene-set selection bias. It should lead.

Layout: donor-level points with group medians, not bars. n donors per group stated on the axis.

## Panel b — the ciliary age gradient in an independent nasal cohort

**GSE336399**, 120 nasal brushings, ages 47–89, bulk RNA-seq, 10 balanced batches.
Plan locked in `docs/PREREG_gse336399_age.md` (commit `04f30ef`) before any expression value was
read; results in `results/gse336399/PREREG_RESULTS.md`.

| gene set | β_age | p (adjusted) | 500 matched random sets at p < 0.05 |
|---|---|---|---|
| **M02 mucociliary axoneme** (literature-curated) | **+0.0234** | 0.0167 | **0.000** |
| **CILIUM28** (this manuscript's own genes) | **+0.0180** | 0.0277 | **0.000** |

Two sub-panels:

- **b1** — scatter of M02 score against age, 120 points, fitted line with CI.
- **b2** — the specificity control: histogram of the 500 matched random-set age coefficients with
  M02 and CILIUM28 marked. Not one random set reached p < 0.05 against a nominal 5%, so this
  cohort has *zero* generic age inflation. That makes the control informative rather than merely
  conservative, and it is the single most persuasive object in the figure.

**The argument the legend must make.** CILIUM28 was derived in the primary cohort by correlating
with a senescence score and then filtering to genes highest in the elderly, so it is biased
toward this result by construction — stated in the locked plan before the run. **M02 carries no
such bias**: literature-curated, fixed long before this cohort was found, never filtered on age
or outcome. The two sets barely overlap, and the unbiased module reproduces the biased one in
direction, magnitude and specificity. That is the claim. Not "the 28 genes replicate" — "an
independent, unselected module reproduces them."

**What the legend must also say**, in the same breath: the pre-specified 13-module primary family
is **null** (smallest BH = 0.218); ρ ≈ 0.17; unadjusted Spearman p = 0.056–0.059; both
pre-declared sensitivity subsets show the direction without significance; the cohort is bulk, so
an increase in axoneme transcription may be a higher ciliated-cell fraction rather than more
transcription per cell; and the cohort starts at 47, so childhood is untestable here.

## Panel c — the childhood arm of the bimodal claim

**Loske 2022**, SARS-CoV-2-negative children (n = 18) vs adults (n = 23), single-cell.
M05 ISG and M06 MHC-I are higher in children in basal, secretory and ciliated cells
(standardised effect 1.6–2.4, BH ≤ 0.002), and the basal and secretory M05 effects beat
expression-matched random sets (empirical p 0.004–0.012 in the final design).

This independently reproduces the published paediatric interferon-tone finding and, with it, the
young half of Figure 4.

**Paired with it, in the same panel, the negative result.** Two independent cohorts show **no**
interferon rise within older adults:

| cohort | range | result |
|---|---|---|
| GSE261277 | 51–80 | ρ = −0.16, p = 0.32 |
| GSE336399 | 47–89 | β = −0.0071, p = 0.306 |

So the bimodality in Figure 4 is supported on its left arm and unsupported on its right. That has
to be stated in the main text, not buried. Either Figure 4's claim is narrowed to the childhood
contrast, or the elderly arm is explicitly marked as cohort-specific and not replicated. **This is
a decision for you, and it should be made before the figure is drawn, not after review.**

## Panel d — outcome anchoring, reported as pre-registered

A small table-figure, one row per outcome-anchored attempt, each with its locked direction, its
result, and its power.

| cohort | outcome | locked prediction | result |
|---|---|---|---|
| GSE261277 pneumococcal challenge (17 POS / 22 NEG) | carriage acquisition | M05 +, M12 +, M01 − | **null**, all p ≥ 0.20, stable across scoring methods; underpowered below d ≈ 0.8 |
| LAIV Gambia | vaccine virus shedding | innate/defence axis − | OR 0.15 and 0.48 per SD, p = 0.007, BH 0.064 — **descriptive only**, those pathways are the published model's own pre-selected features |
| GSE253686 NOVA influenza challenge (74 adults) | viral shedding days | *to be locked* | **pending** — see below |

Reviewers at this journal will look for exactly this panel and will not find it in most papers.
A clean null with a declared power rule is worth more than a borderline positive with an
undeclared one. Present it as a strength.

## Panel e — the measurement caveats (candidate for supplementary)

Three findings that are about measurement, not nasal biology, and that hold regardless of which
cohorts survive:

1. **CDKN2A+ status is largely a sequencing-depth event.** 4,265 vs 2,888 median detected genes;
   +606 in the same direction in 181 of 187 donors. Depth matching within donor reverses an
   axoneme contrast from p = 4.5e-27 one way to p = 0.020 the other.
2. **Never score a senescence module containing the gene that defines the groups.** Removing
   CDKN2A and matching depth shrank the contrast sixfold.
3. **Expression-matched random sets are not optional.** In one cohort 11.8–33.4% of matched random
   sets tracked real age, and 41.6% of random 4-gene sets reached p < 0.05 in the CDKN2A
   contrast. In GSE336399, 0%. The control tells you which regime you are in; nothing else does.

Point 1 is a genuine methods contribution to the p16 single-cell senescence literature, where
groups are routinely defined by CDKN2A detection without depth matching. It may be worth a short
paragraph in the discussion as well as a supplementary panel.

**Also for supplementary:** M12 inflammatory SASP is unmeasurable in airway epithelium on 10x 3′
chemistry — **4 datasets out of 4** (primary cohort, Ziegler 2021, Loske 2022, Chua 2020). In our
own ciliated cells 0 of 7 genes reach 5% detection. Agreed wording is in `docs/FIGURE_WORDING.md`.

---

## The pending piece: GSE253686 (NOVA influenza challenge)

This is the one addition that would change panel d from "pre-registered hypotheses with their
nulls honestly reported" to "a pre-registered prediction tested against a real clinical endpoint."

**What is in hand.** The series matrix is downloaded
(`D:\research\nasage_validation\data\gse253686\`): 297 nasal-epithelium microarrays, **75
participants**, four timepoints each (**SD −1 baseline**, SD1, SD3, SD7), vaccinated vs
unvaccinated, GPL6480. Baseline samples are exactly what NasAge needs: pre-challenge nasal state.

**What is missing.** GEO carries no age and no outcome. Per-participant shedding group
(NS / 1DS / MDS) and age are in the paper's *data file S1* and *S16*, which are hosted on PMC
behind a download gate that blocks scripted retrieval. One browser click gets them.

**One metadata defect already found**, recorded now rather than after the analysis: three
participants (NOVA01, NOVA03, NOVA10) carry **conflicting vaccination labels across their own
timepoints** in GEO. They will be excluded, or their status resolved from the paper, and this is
pre-declared.

**If the file can be obtained**, the locked plan would be: baseline (SD −1) NasAge axes →
multi-level shedding outcome, Poisson or ordinal model adjusted for age, sex and vaccination
status, the usual 500 expression-matched random sets, direction fixed in advance per module,
BH within the measurable family, power table declared before any p-value. Same discipline as
GSE261277 and GSE336399.

---

## What Figure 5 should not try to do

- It should not re-run or re-draw Figures 1–4. Those stand on the primary cohort.
- It should not claim bimodality is externally validated. Half of it is.
- It should not present CILIUM28 alone as the replication. M02 is the load-bearing result.
- It should not rest any claim on Crownlands, which is withdrawn (olfactory brush biopsy, and
  module age effects inseparable from per-cell detection depth).
- It should not assert that any axis predicts a vaccine or infection outcome. On current evidence
  one such test is null and one is descriptive.
