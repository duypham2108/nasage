# Lifespan trajectory of mucosal defence programmes, by cell type

Status: **LOCKED** by the commit that adds this file.

Origin: the design is the user's — five mucosal defence axes profiled across age bins, with a
childhood peak as the proposed shape. This plan implements that idea with **referenced gene
sets instead of AI-generated ones**, per cell type, with the project's existing controls.

---

## 1. What changes from the original sketch, and why

The user's five axes were assembled manually with an AI assistant and carry no references. Every
gene set in this project carries a per-gene reason and a citation (`docs/MVRI_PROGRAMS.md`,
`gene_set_provenance_genes.csv`). Four of the five axes already have curated equivalents, so the
idea is testable without introducing an unreferenced set:

| User's axis | Curated equivalent used here | genes |
|---|---|---|
| Mucus permeability | **P5 mucin barrier** | 11 |
| Antigen presentation | **M06 MHC-I processing** | 15 |
| Innate sensing | **M04 innate sensing** | 15 |
| IgA support | **P4 IgA transport** | 10 |
| Cytokine tone | **M12 inflammatory SASP** | 7 |

**M12 is predicted to fail the measurability gate** (0–14% of genes detected in all four
single-cell cohorts tested, Supp. Fig. S5c). That prediction is recorded here, before the run:
if it fails again, the Cytokine tone axis cannot be plotted and the compressed-near-zero
appearance of that panel in the original sketch is a detection-floor artefact, not flat biology.

**The title changes.** No panel under this plan measures a vaccine outcome, so the words
"vaccine readiness" are forbidden here exactly as in `PREREG_celltype_mvri.md` §5 and
`PREREG_epithelial_only.md` §4.

## 2. Cohorts, and the honest span

| Cohort | age | donors | role |
|---|---|---|---|
| **Loske 2022** | **4–77, continuous** | 41 | **primary** — the only cohort spanning childhood to old age with cell types annotated |
| Ziegler 2021 | 3rd–9th decade (bins) | 58 | adult/older replication |
| GSE336399 | 47–89 | 120 | older-adult replication, bulk-style donor scores |

**There are no neonates or infants in any available cohort.** The claim is therefore
**childhood to old age (4–77)**, never "across the lifespan", and the figure must not imply
coverage the data do not have. Crownlands is excluded (withdrawn for biological claims).

## 3. Primary analysis

**Unit: donor × cell type**, ≥ 20 cells, as in `PREREG_celltype_mvri.md`. Cell types: basal,
secretory, ciliated.

For each (programme, cell type) in Loske:

- **Model:** score ~ age + sex + log_genes, donor-level OLS. `log_genes` is in the model because
  the CDKN2A depth confound (Supp. Fig. S5b) showed depth can manufacture an age-like effect.
- **Primary estimate:** the linear age coefficient per year.
- **Non-monotonic shape, pre-specified:** a quadratic age term is fitted **as a declared second
  model**, not chosen after inspecting the curve. The childhood peak in the original sketch is a
  non-monotonic claim, so it needs a test that can express one. Reported as ΔAIC between linear
  and quadratic, with the quadratic coefficient and its p.
- **Family:** 5 programmes × 3 cell types = 15 tests. **BH across all 15**, with the effective
  number of tests reported beside it (`nasage.stats.effective_tests`).
- **Measurability gate** applied per cell type before any model; a programme failing in a cell
  type is not tested there.

## 4. Controls

1. **500 expression-matched random gene sets**, same model, per programme × cell type. A
   programme failing its matched null is `not_specific`, whatever its p.
2. **Age-label permutation**, 2,000 draws, to confirm the pipeline is not anticonservative on
   this cohort.
3. **Composition check.** A trajectory can arise because the cell mixture changes with age
   rather than the programme changing within a cell. Because the unit is donor × cell type, this
   is already largely controlled; in addition, **cell-type proportion versus age is reported
   alongside**, so a reader can see whether composition moves.

## 5. Replication, declared in advance

The adult portion (≥ 18 y) is re-tested in Ziegler (decade bins, ordinal) and GSE336399
(47–89, continuous). Agreement in **direction** is the replication criterion; magnitudes are not
expected to match across platforms.

## 6. What is not independent, stated now

**M06 and M04 trajectories in Loske are not independent of current Figure 5 panels d and e**,
which already report a child-versus-adult contrast for M05 and M06 in this same cohort. A
childhood peak for M06 here is the **same observation in a finer-grained form**, not
corroboration, and must be described that way.

Further, the module audit found M05/M06/M07/M09/M03 correlate 0.72–0.93 while sharing no genes
(`docs/MODULE_DESIGN_AUDIT.md` §1.1). **M06 "antigen presentation" and M04 "innate sensing" may
be one interferon axis**, in which case two of the five panels are one finding. The
effective-tests figure in §3 will show this, and the legend must state it.

**Genuinely new under this plan:** the age trajectories of **P5 mucin barrier, P4 IgA transport
and M12** per cell type across 4–77. None has been tested against age in this cohort.

## 7. Verdicts

| Observed | Reading |
|---|---|
| linear age effect, BH < 0.05, clears matched null | programme changes with age in that cell type |
| quadratic better by ΔAIC > 2 **and** quadratic term BH < 0.05 | **non-monotonic**; a peak may be described |
| quadratic not supported | the trajectory is described as monotonic; **no peak is claimed** even if the plotted line looks curved |
| fails matched null | `not_specific`; plotted only as grey, with no claim |
| fails measurability gate | not plotted |

A peak is a claim and requires the quadratic test to pass. An eye-catching curve in a line plot
is not evidence of one.
