# Figure 5 — draft legend and source data

Figure files: `D:\research\nasage_validation\results\figure5\figure5.pdf` (vector, 180 mm wide,
Type-42 fonts, editable text) and `figure5.png` (400 dpi).
Built by `validation/figure5_plot.py` from `validation/figure5_data.py`.

Revised 2026-11-24: the ALI panel was removed (no source data exists), and Figure 4 is now
narrowed to the childhood contrast, so panel d documents the basis for that narrowing rather
than reporting a failed replication.

---

## Draft legend

**Figure 5. External validation of the nasal ageing axes, with the tests that failed.**

**(a)** Mucociliary axoneme score (module M02, mean z of log₂(CPM+1) over 18 expressed genes)
against age in an independent cohort of 120 nasal brushings, ages 47–89 (GSE336399, bulk
RNA-seq). Line, ordinary least squares; shaded band, 95% confidence interval. Inset: Spearman ρ
with its unadjusted p, and the age coefficient from a model adjusted for sex, cancer status,
smoking, RNA integrity (median TIN), log₁₀ library size and sequencing batch. The 13-module
Benjamini–Hochberg family is not significant (BH = 0.22); the analysis plan was committed before
the data were accessed.

**(b)** Specificity of the result in (a). Grey, the distribution of age coefficients from 500
expression-matched random gene sets (25 expression bins, identical model); vertical lines, the
observed coefficients for M02 (orange) and for the 28 cilium-assembly genes derived in the
primary cohort (CILIUM28, green). Neither gene set's coefficient is approached by any random set,
and **0 of 500 random sets reached p < 0.05** against a nominal 5%, indicating no generic
age-related inflation in this cohort.

**(c)** Interferon-stimulated gene score (module M05, UCell, donor means) in basal, secretory and
ciliated epithelium of SARS-CoV-2-negative children (n = 18) and adults (n = 23) (Loske et al.
2022). Horizontal bars, group medians. Below each cell type: Benjamini–Hochberg-adjusted p from
the primary family and the empirical p against 500 expression-matched random gene sets
(*, empirical p < 0.05).

**(d)** The same module against age *within* older adults, in two independent cohorts: GSE336399
(n = 120, ages 47–89) and GSE261277 baseline samples (n = 39, ages 51–80). Dashed line, ordinary
least squares; Spearman ρ and its unadjusted p inset. Neither cohort shows an age gradient in
interferon tone among adults, which is why the claim in Figure 4 is restricted to the
child–adult contrast.

**(e)** Detected genes per cell in ciliated cells, split by CDKN2A detection (CDKN2A⁻, blue;
CDKN2A⁺, orange); dashed vertical lines, medians. CDKN2A⁺ cells are sequenced more deeply
(median 4,265 vs 2,888 detected genes; the difference has the same sign in 181 of 187 donors), so
p16 positivity called from single-cell detection is substantially a sequencing-depth event. Data
from an olfactory brush-biopsy cohort used here for measurement properties only.

**(f)** Fraction of each module's genes reaching ≥ 5% detection in ciliated cells across four
independent single-cell cohorts spanning two chemistries. Dashed line, the 60% measurability
gate applied throughout. M02, M05 and M11 clear the gate in every cohort; **M12 inflammatory
SASP fails in all four** (0–14% of genes detected), so secreted-cytokine SASP is below the
detection floor of droplet scRNA-seq in airway epithelium and is reported with no claim attached.

**(g)** Outcome-anchored tests, each reported as locked. Predictions, direction and multiplicity
were fixed by a committed plan before the data were accessed, and results are reported whichever
way they fell.

---

## What each panel is allowed to claim

| Panel | Claim it supports | Claim it does **not** support |
|---|---|---|
| a | axoneme transcription rises with age across 47–89 | significance within the pre-specified family (BH = 0.22); anything about childhood |
| b | the rise is specific to these gene sets, not generic | that the effect is large (ρ ≈ 0.17) |
| c | higher interferon tone in children's epithelium, module-specific in basal and secretory | ciliated cells (empirical p = 0.076) |
| d | no interferon gradient within older adults, in two cohorts | that interferon tone is flat across the whole lifespan |
| e | p16⁺ calls track sequencing depth | any biological conclusion — the cohort is withdrawn |
| f | SASP cytokines cannot be measured in droplet scRNA-seq of nasal epithelium | that SASP is biologically absent |
| g | the axes were tested against real outcomes and reported honestly | that any axis predicts an outcome |

**The single most important sentence for the discussion:** the 28-gene cilium set was derived in
the primary cohort by correlation with a senescence score and then filtered to genes highest in
the elderly, so it is biased toward the result in panels a/b by construction. **M02 carries no
such bias** — literature-curated, fixed before this cohort was identified, never filtered on age
or outcome — and reproduces it in direction, magnitude and specificity. The replication claim
should rest on M02, with CILIUM28 as the corroborating secondary.

**Panel f is a positive contribution, not an apology.** "This programme cannot be measured on
this platform, in four cohorts across two chemistries" is a useful result for anyone scoring SASP
modules in droplet data, and it is the reason M12 carries no claim anywhere in the paper.

---

## Source data per panel

| Panel | File |
|---|---|
| a | `results/gse336399/module_scores.csv`, `results/gse336399/primary_age_effects.csv` |
| b | `results/figure5/gse336399_null_draws.csv.gz`, `results/gse336399/matched_null_control.csv` |
| c | `results/loske2022_v2/donor_celltype_scores.csv`, `results/loske2022_v2/age_effects_with_null.csv` |
| d | `results/gse336399/module_scores.csv`, `results/pneumo_gse261277/module_scores_all_samples.csv` |
| e | `results/crownlands_figrep/ciliated_cdkn2a.csv.gz` |
| f | `results/own_cohort/module_measurability.csv`, `results/ziegler2021_v2/module_coverage.csv`, `results/loske2022_v2/module_coverage.csv`, `results/crownlands/measurability_gate.csv` |
| g | `results/pneumo_gse261277/PREREG_RESULTS.md`, `results/laiv_gambia/go_vs_outcomes.csv` |

`validation/figure5_data.py` regenerates the 500 matched-null draws behind panel b. It re-runs the
control exactly as originally reported — same module order, seed, bins and model — and asserts
that the observed coefficient, empirical p and fraction of random sets at p < 0.05 for M02 and
CILIUM28 match `matched_null_control.csv` to machine precision. They do.

---

## Known gaps in this draft

1. **Panel g row 3 is pending.** GSE253686 (NOVA, A/H1N1 challenge) has 75 participants with
   baseline SD −1 nasal arrays downloaded, but age and shedding outcome live in the paper's
   supplementary data files S1/S16, behind a PMC download gate. The analysis plan must be locked
   *before* those are opened.
2. **Panel e comes from a withdrawn cohort.** The finding is methodological and holds, but the
   legend must say the cohort is olfactory brush biopsy and supports no biological claim.
3. **The "Chua 2020" citation for the SASP floor has no result file.** `RESULTS.md` and
   `FIGURE_WORDING.md` both say M12 fails in "4 of 4" cohorts including Chua 2020, but no Chua
   output exists in the repository or in `nasage_validation/results/`. Panel f therefore uses the
   four cohorts that **do** have verifiable detection tables: the primary cohort, Ziegler 2021,
   Loske 2022 and the olfactory Flex cohort. Either locate the Chua analysis or change the text
   to cite these four — the current wording cannot be reproduced from the files on disk.
4. **No functional/protein-level panel.** The transcriptomic ciliary result in panels a/b is not
   corroborated at protein level anywhere in the figure. See `FIGURE5_PLAN.md` for what would
   fill that gap.
