# Build plan: manuscript Figure 5

Drafted 2026-11-25, after `docs/PREREG_mvri_presentation_audit.md` and its results
(`D:\research\nasage_validation\results\mvri_audit\PREREG_RESULTS.md`).

Scope: the atlas manuscript's Figure 5 (currently "Clinical translation: the mucosal vaccine
readiness index across the lifespan"). Not the NasAge validation figures, which keep their own
numbering in `manuscript/figures_main.tex`.

Figures 1–4 are untouched by this plan.

---

## 0. What the audit forces

| Draft panel (page 6 of the plan PDF) | Status after the audit |
|---|---|
| a — five axis cards | keep, but two axes must be marked unmeasured |
| b — five per-cell boxplots across age bins | **cannot be published as drawn**: cells are not independent samples of donor age; every p underflows to 0.0 |
| c — axis × age-bin heat map (z) | rebuild at donor level |
| d — weighted composite curve | **demote**: computable for 36/52 donors, Elderly n = 1 |
| e — empty | replaced below |

Donor-level numbers already computed, for reuse rather than recomputation:

| Axis | measurable | gate | KW p (cell) | KW p (donor) | BH | empirical p |
|---|---|---|---|---|---|---|
| Mucus permeability | 0.93 | pass | 0.0 | 0.090 | 0.150 | 0.340 |
| Antigen presentation | 0.71 | pass | 0.0 | 0.069 | 0.150 | 0.252 |
| Innate sensing (PRR) | 0.55 | **FAIL** | 0.0 | 0.322 | 0.403 | 0.838 |
| IgA support | 0.75 | pass | 0.0 | 0.074 | 0.150 | 0.146 |
| Cytokine tone | 0.54 | **FAIL** | 0.0 | 0.467 | 0.467 | 0.870 |

So the figure has **three measurable axes, none significant at donor level, none specific
against matched random gene sets**. The honest headline is *maturation of the mucosal
interface*, with the trends reported as suggestive and underpowered. "Readiness" cannot be
claimed: no vaccine outcome exists in this cohort.

## 1. Proposed title

> **Figure 5 | Maturation of the mucosal interface across the lifespan**

Dropping "readiness" and "clinical translation". See §5 for the alternative that would earn
those words back.

## 2. Main-figure panels

**a — The interface, and what is measurable in droplet data.**
Five axis cards as in the draft, each annotated with the fraction of its genes detected in ≥ 5%
of the cells it is scored in (0.93 / 0.71 / 0.55 / 0.75 / 0.54). The two below the 60% gate are
greyed with the label *below detection floor*. Cell subset and gene count printed on each card.
This panel does the work of saying, on the figure itself, why only three axes continue.

**b — Donor-level axis trajectories.** Three measurable axes, seven age bins, **one point per
donor**, median bar, n per bin on the axis. Kruskal–Wallis p printed per panel (0.090 / 0.069 /
0.074) with the BH-adjusted value. This is the direct replacement for draft panel b and is the
single most important change in the figure.

**c — Two age windows.** Per-axis slope with CI for Neonate→Young child and Child→Elderly,
drawn in the same idiom as Figure 4a so the reader can place it. Separates "rises through
childhood" from "declines after", which the seven-bin KW test cannot.

**d — Cellular source of each axis.** Dot plot, axis × cell type (Basal, Secretory, Goblet,
Ciliated, plus Myeloid and Plasma/B where the deck assigns them), mean score and fraction of
cells. Ties Figure 5 back to the basal/secretory/ciliated primary family and shows the reader
that the antigen-presentation axis is not an epithelial-only quantity.

**e — Composite, honestly drawn.** Donor-level MVRI across bins with per-bin n printed
**as computed, not as declared**: 36 of 52 donors, Elderly = 1. Grey band over Elderly marked
*n = 1*. If you would rather not carry a composite at all, this panel moves to supplementary
and the main figure ends at d — my preference, and it is the cleaner paper.

## 3. Supplementary Figure 5

a. **Per cell versus per donor, same axis, side by side** (p = 0.0 vs 0.090). One panel that
   justifies the whole rebuild.
b. **Matched-null specificity**: observed KW statistic against 500 expression- and
   size-matched random gene sets, per axis. Shows no axis clears its null.
c. **Per-gene detection** within each axis, so the two gate failures are legible gene by gene
   (TLR family, CGAS/STING, IL33/TSLP/IL6 at the droplet floor).
d. **Weight sensitivity**: distribution of which bin peaks over 10,000 sign-respecting weight
   vectors (childhood 94.6%; Child 84.5%, Young child 10.1%). Note this *passed* — the weights
   are not what undermines the index.
e. **Depth sensitivity**: axis score vs mean genes per cell (rho +0.38 to +0.81), mean
   inter-axis rho 0.468 → 0.367 after removing depth, depth vs age-bin rho +0.085 (p = 0.62).
f. **Composite coverage**: which donors lack which axis subset, by age bin.
g. **Dimensionality**: pairwise axis correlation matrix, Li & Ji 4.0 effective tests, PC1 58%.

## 4. The overlap to decide before building

Figure 4 is titled "Programmes that peak in childhood: interferon signalling and **antigen
presentation**", and Figure 4a/4c already carry antigen presentation across Basal, Secretory,
Goblet and Ciliated with LTSR support. Figure 5b's antigen-presentation axis is the same
biology at lower resolution and lower significance.

Three ways out, in my order of preference:

1. **Fold the surviving content into Figure 4** and drop Figure 5. The paper becomes four
   figures, and nothing is lost except an index that does not survive its own audit.
2. **Keep Figure 5 but make mucus permeability and IgA support its subject**, since neither
   appears in Figure 4, and cite Figure 4 for antigen presentation rather than redrawing it.
3. Keep as planned in §2 and accept the redundancy.

## 5. What would earn "clinical translation" back

The index has no outcome attached to it, and no amount of re-drawing fixes that. The only data
in reach with a real mucosal outcome is the **pneumococcal carriage** cohort used in the NasAge
validation work, where a nasal programme reaches empirical p = 0.000 under all four scorers.
A Figure 5 built on that would carry the translational claim the current one asserts.

That is a scope decision, not a plotting one: it brings a second paper's cohorts into the atlas
manuscript. Flagging it, not assuming it.

## 6. Build order and inputs

Everything needed for §2b–e and §3a–g is already on disk; nothing requires re-scoring.

- `D:\research\nasage_validation\results\mvri_audit\donor_axis_scores.csv` — donor × axis
- `.../axis_means_by_agebin.csv`, `.../verdicts.csv`, `.../null_draws*.csv`
- `.../depth_sensitivity.csv`, `.../depth_summary.csv`
- axis definitions and subsets: `validation/mvri_presentation_audit.py` (verbatim from slides 9–14)

New script: `validation/figure5_maturation_plot.py`, house style as the other figure scripts
(Okabe–Ito, Type-42 fonts, 180 mm width, two-row grid for the main figure).

Panel d is the one piece not yet computed — it needs a per-cell-type mean of each axis score,
a single pass over the primary object, a few minutes.

---

# Built, 2026-11-25

`validation/figure5_maturation_plot.py` →
`D:\research\nasage_validation\results\mvri_audit_pred\figures\figure5_maturation.pdf`
and `figure5_maturation_supp.pdf` (180 mm wide, Type-42 fonts, the age ramp and title style
of Figures 1–4).

**What changed from the plan above, and why.**

The plan was written against the locked audit, where nothing was significant and the composite
covered 36 of 52 donors. Building it uncovered the reason: `Cell_type_annotation_level1/2` are
NaN for 53,150 cells and 16 whole donors, so the audit — and the deck — had been scoring 69% of
the object. Re-run on `Cell_type_prediction_level*`, which covers every cell
(`validation/mvri_audit_amendment.py`, results in `results/mvri_audit_pred/`):

* **antigen presentation** (empirical P 0.010) and **IgA support** (empirical P 0.016) are
  age-varying and clear their matched nulls, so Figure 5 now carries a positive result;
* mucus permeability is age-varying but **not** specific (empirical P 0.128);
* innate sensing and cytokine tone still fail the measurability gate (0.45, 0.54);
* the composite is computable for **all 52 donors**, so panel e no longer carries an n = 1 bin.

Two consequences for the figure as drawn:

1. **Panel b headlines the empirical P, not the BH P.** With all three source datasets pooled,
   53–80% of random matched gene sets reach a nominal donor-level P < 0.05, so nominal and BH
   values are uninterpretable. Supplementary panel b states the calibration under each null.
2. **Supplementary panel f was replaced.** The planned donor-coverage heat map became
   uninformative once coverage was complete; it now shows source dataset by age bin, which is
   the dominant caveat — the Elderly bin is four-fifths the `senescence` dataset and the Adult
   bin is entirely `ucl`, so the post-childhood decline cannot be separated from batch.

The §4 overlap decision (Figure 4 already owns antigen presentation) and the §5 question
(whether to bring in a cohort with a real outcome) are both still open and still yours.

---

# Revised, 2026-11-25 (second pass)

Panel a became a Sankey (genes -> axes -> cell types) on request, and the figure was then
trimmed: **axes that fail the measurability gate and genes below the 5% detection floor are
removed from the score as well as from the figure**, so panel a is a true statement of the
inputs rather than a filtered view of them.

Run: `validation/mvri_audit_epithelial.py <h5ad> <out> 500 all` ->
`results/mvri_audit_trim/`. The same script with `epi` restricts every axis to epithelium;
that variant was explored and is recorded below but is not what the figure uses.

Dropping the undetected genes **strengthens** both surviving results, because a mean-of-z
over a set padded with undetected genes is diluted:

| Axis | genes scored | empirical P before | after |
|---|---|---|---|
| Antigen presentation | 11 of 17 | 0.010 | **< 0.002** (0 of 500) |
| IgA support | 3 of 4 | 0.016 | **0.002** |
| Mucus permeability | 13 of 15 | 0.128 | 0.066 (still not specific) |
| Innate sensing | — | — | gate fail, 45% |
| Cytokine tone | — | — | gate fail, 54% |

**Epithelium-only variant** (`results/mvri_audit_epi/`, not used for the figure): antigen
presentation strengthens further to empirical P < 0.002 on 66,925 epithelial cells, but
**IgA support becomes unscorable** — only 2 of its 4 genes clear the floor once Plasma and
B mem are excluded, which shows its signal is partly carried by the B/plasma compartment
rather than by epithelium. Worth stating in the text even though the figure keeps all cell
types.

Document: `manuscript/figure5.tex`, two pages (figure, then legend), compiled with tectonic
from the figure directory. The §4 overlap question (Figure 4 already owns antigen
presentation) and the §5 question (a cohort with a real outcome) remain open.
