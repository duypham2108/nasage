# Figure wording, with the evidence each claim rests on

Working file. Each claim below is paired with the analysis that supports it and, where it
matters, the thing it must **not** be allowed to say. Update as results land; do not delete
superseded wording, strike it and date the replacement.

Last updated 2026-11-24.

---

## Figure 1 — basal cell loss and a transcriptomic clock

*Outline:* basal cells decrease with age (c); a pre-trained transcriptomic clock (Pasta,
Salignon et al., Adv Sci, doi:10.1002/advs.76740) tracks donor age in the epithelium (d).

**Status: not yet checked by me.** No wording fixed. Two things to settle before drafting:

1. Composition claims need the donor-level denominator stated. The primary cohort has
   Neonate 6, Infant 7, Young child 9, Child 8, Adolescent 7, Adult 6, **Elderly 5** donors.
   A basal-fraction trend across the lifespan rests on 5 elderly donors at the far end.
2. Pasta is pre-trained on other tissues; whether it was trained on any nasal or airway data
   needs checking, because a clock trained on data overlapping this cohort would be circular.

---

## Figure 2 — basal subsets diverge

*Outline:* Basal 1 does not change with age; Basal 2 advances along pseudotime with age.

**Status: not yet checked by me.** No wording fixed. Note for drafting: a pseudotime-versus-age
correlation needs a donor-level unit of analysis (mean pseudotime per donor), not a cell-level
one, or the n is the cell count and every p-value is meaningless.

---

## Figure 3 — senescence differs by cell type

*Outline:* CDKN2A+ ciliated cells accumulate with age and show reduced cilium movement and
increased senescence/SASP programmes. Basal cell loss and basal pseudotime progression both
correlate with the rise in p16+ ciliated cells.

### Agreed wording (2026-11-24)

> CDKN2A+ ciliated cells accumulate with age and show increased senescence-arrest programmes
> and a senescence-associated secretory phenotype dominated by IGFBP and matrix components.
> Classical pro-inflammatory SASP cytokines (IL6, IL1B, TNF, CXCL8) are below the detection
> floor of droplet scRNA-seq in this tissue.

### What supports each clause

Evidence: `validation/check_module_measurability.py` and
`validation/check_senmayo_detection.py`, run on the primary cohort
(148,978 cells, 52 healthy donors). Tables in
`D:\research\nasage_validation\results\own_cohort\` (`module_measurability.csv`,
`senmayo_detection.csv`, `MEASURABILITY.md`).

**"senescence-arrest programmes" — supported.** Fraction of ciliated cells with any count:
CDKN2A 0.206, CDKN1A 0.297, GDF15 0.284, IGFBP7 0.764. SERPINE1 (0.022) is below the floor and
contributes noise. 4 of 5 M11 genes clear 5%.

**"reduced cilium movement" — supported.** M02 mucociliary axoneme is 20/20 genes detected at
0.689 median in ciliated cells, the best-measured module in the dataset. A decrease in this
score is a real measurement, not a detection artefact.

**"dominated by IGFBP and matrix components" — supported, and the hedge is required.**
SenMayo (SAUL_SEN_MAYO, 124 genes) has 37 genes clearing 5% detection in ciliated cells, versus
**0 of 7 for the curated M12 module**. The detected genes are led by IGFBP2 0.846, IGFBP7 0.764,
CD9 0.692, MIF 0.676, HMGB1 0.574, CTSB 0.551, IGFBP5 0.441, CTNNB1 0.338, EGFR 0.213,
JUN 0.306. These are matrix, growth-factor-binding and intracellular genes. The only secreted
inflammatory mediators clearing the floor are C3 0.186, IL18 0.173, CXCL16 0.118, ANG 0.114,
IL7 0.061.

**"below the detection floor" — supported in 4 of 4 datasets.** In ciliated cells: CXCL1 0.030,
CXCL8 0.019, IL1B 0.004, IL1A 0.001, TNF 0.001, IL6 0.000, CCL2 0.000. The same seven genes fail
in Ziegler 2021, Loske 2022 and Chua 2020. This is an assay property, not a cohort property, and
it should be stated as a positive contribution: nasal SASP requires protein measurement
(nasosorption fluid, Luminex) or a targeted amplified assay.

### Two constraints on how Figure 3 may be built

**1. SenMayo and M11 are not independent.** SenMayo contains all 7 M12 genes and 3 of the 5 M11
genes (GDF15, IGFBP7, SERPINE1). IGFBP7 is simultaneously M11's best-detected gene in ciliated
cells (0.764) and one of SenMayo's top drivers. Reporting "increased senescence **and**
increased SASP" as two findings would be reporting one measurement twice.

*Rule:* score SenMayo with the M11-overlapping genes removed. 35 of the 37 detected genes
remain, so nothing material is lost. State the removal in the legend.

**2. The SenMayo panel needs the random-set control.** Thirty-seven broadly expressed genes will
track sequencing depth and library complexity. Without expression-matched random gene sets
(>= 500, 25 expression bins, as the package implements), a reviewer can reasonably claim the
panel plots detection rate rather than biology. This control is not optional for this panel.

### What Figure 3 must not say

- Not "increased SASP" unqualified. The cytokine arm is unmeasured, not measured-and-high.
- Not "M12 rises with age in ciliated cells." M12 is `poorly_measurable` in every cohort and is
  excluded from the primary family. Any M12 number in this cohort is a CXCL1 number.
- Not a transfer of the bulk LAIV M12 result. M12 rose after live intranasal vaccine in bulk
  RNA-seq, but that rise was shown not to be live-vaccine-specific
  (`results/module_confirmation/RESULT_confirmation_FAILED.md`), and it cannot be measured at
  single-cell level at all.

---

## Figure 4 — interferon programmes peak in childhood, bimodally

*Outline:* in line with Yoshida et al. Nature 2022, IFN programmes peak in childhood; the
response is bimodal, with classical ISGs peaking only in childhood and a distinct IFN-signalling
cluster high again in the elderly.

**Status: not yet checked by me.** No wording fixed. Three notes for drafting:

1. M05 is well measured in this cohort: 22/22 genes present, 14-20 genes clearing 5% detection
   across basal, secretory and ciliated. The measurement is sound.
2. This figure cannot be replicated in the Crownlands cohort, which has no donor under 18. It
   stays a single-cohort result, and the external replication (Fig 5b) covers the adult-to-
   elderly half only. Say so rather than letting a reader assume Fig 5b covers Fig 4.
3. The bimodality is the interesting claim and will attract the most scrutiny. A second
   IFN-signalling cluster high in the elderly needs the random-set control and an explicit
   statement of which genes separate the two clusters.

---

## Figure 5 — external replication and the vaccine framing

Planned panels (analysis locked in `docs/PREREG_crownlands_age.md`, commit `0d897a4`, and
`docs/PREREG_vaccine_axes.md` sections 8.1-8.3):

- **a.** Module x age-bin heatmap across basal / secretory / ciliated in the primary cohort,
  excluding the axes Figures 3 and 4 already own (so: barrier, mucins, PIGR/IgA transport,
  MHC-I, plasma secretory machinery; not ISG, not senescence).
- **b.** The same age effects in the Crownlands Gateway atlas: 192 healthy donors aged 20-90,
  independent cohort, different platform and sampling site. Readout is a concordance scatter of
  age-effect sizes with a rank correlation.
- **c.** Expression-matched random-gene-set null, >= 500 sets.
- **d.** The pre-registered direction table: predicted sign of each axis per vaccine type, with
  each outcome-anchored test marked by what it returned, nulls included.

### Wording fixed in advance for panel b

Whatever it shows:

> Age effects were compared across cohorts by rank concordance of effect sizes, never pooled.
> The two cohorts use different platforms (10x 3' versus 10x Flex probe-based) and different
> sampling sites (inferior turbinate versus olfactory cleft), so absolute module scores are not
> comparable and only within-cohort age gradients were tested.

And the limitation that must appear in the legend, because it is the weak side of the
comparison:

> The 65+ versus young-adult contrast rests on 6 adult and 5 elderly donors in the primary
> cohort against 39 and 90 in the replication cohort; a module significant in the replication
> cohort and not in the primary cohort is reported as concordant if the directions agree, not
> as a failure.

### Wording fixed in advance for panel d

> The interferon axis is predicted to act in opposite directions for different mucosal vaccine
> types: higher interferon tone restricts a replicating viral vaccine while impairing control of
> pneumococcal colonisation. No single mucosal readiness score can therefore exist, and
> readiness must be defined per pathogen.

This prediction was registered in `docs/PREREG_vaccine_axes.md` section 8.2 (commit `5d0b924`)
before the relevant data were accessible.

### What Figure 5 must not say

**No claim that NasAge predicts mucosal vaccine response.** Six outcome-anchored attempts have
been made and none produced a positive result:

| Attempt | Endpoint | Outcome |
|---|---|---|
| GSE230494 H1-H4 | LAIV shedding, nasal IgA | not testable, no participant key |
| GSE230494 viral reads | vaccine replication | endpoint absent, 0/9 day-3 samples |
| GSE155237 | RSV challenge | positive control failed |
| Set A derivation | per-gene | failed on power |
| Module confirmation | LAIV vs TIV arms | failed on power |
| GSE261277 | pneumococcal carriage acquisition | **null**, correctly specified, underpowered below d ~ 0.8 |

These datasets are public and a reviewer will find this. The vaccine axes enter the manuscript
as **pre-registered, direction-specific hypotheses with their tests reported, nulls included**,
never as validated predictors. That framing is a strength: it states exactly what the public
data can and cannot establish.

---

## Standing hypothesis, not a result

Recorded so the direction is on record before any future test
(`results/pneumo_gse261277/PREREG_RESULTS.md`):

> Baseline M03 (PIGR) and M09 (plasma secretory machinery) are lower in individuals who go on to
> acquire pneumococcal carriage.

In GSE261277 this was the largest exploratory signal (M09 d = -0.91, M03 d = -0.58, both lower
in future carriers) and it is biologically coherent, PIGR being the receptor that transports
secretory IgA across the epithelium. It did **not** survive BH correction, it was one of twelve
two-sided exploratory tests, and the rank and adjusted p-values disagree by an order of
magnitude. It must not be reported as a finding from that dataset.
