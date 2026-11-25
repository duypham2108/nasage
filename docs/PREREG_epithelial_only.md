# Epithelial-only restriction of the programmes, and the re-test

Status: **LOCKED** by the commit that adds this file.

This is a **declared sensitivity analysis on cohorts already analysed**, not a fresh
pre-registered discovery test. That distinction is the whole point of this section and is
carried into every sentence of the results.

---

## 1. Why a restriction, and why it is legitimate

The seven programmes are labelled epithelial — basal, secretory, ciliated. Measured against
Ziegler 2021 nasopharyngeal single-cell data (28,948 epithelial vs 3,640 immune/other cells),
several are not:

| Programme | mean epithelial fraction | genes below 0.5 |
|---|---|---|
| **P4 IgA transport** | **0.497** | **4** — TNFSF13B 0.10, JCHAIN 0.22, IL6 0.24, IL17RA 0.30 |
| P6 entry-factor availability | 0.638 | 6 — NEU1, CTSL, NEU3, ANPEP, FURIN, ST3GAL6 |
| P3 secretory antimicrobial | 0.660 | 2 — LYZ 0.21, S100A8 0.46 |
| P1 basal progenitor/repair | 0.689 | 2 — MKI67 0.44, TOP2A 0.49 |
| P7 ciliary beat machinery | 0.744 | 0 |
| P5 mucin barrier | 0.752 | 0 |
| P2 junction integrity (excluded) | 0.854 | 0 |

Epithelial fraction = mean log-normalised expression in epithelial cells ÷ (epithelial +
immune), computed per gene across all cells. **JCHAIN is a plasma-cell transcript and TNFSF13B
is BAFF**; in bulk nasal tissue P4 has been substantially reporting immune infiltration under
the name "epithelial IgA transport".

**The filter is outcome-blind.** It is computed from single-cell expression only, in a cohort
(Ziegler) that contributes no outcome to any test in this project. No outcome variable, in any
cohort, was consulted in defining it. That is what makes it admissible at all.

## 2. The filter, fixed now

- A gene is **epithelial** if epithelial fraction **≥ 0.5**.
- A programme is retained if **≥ 60% of its genes are epithelial** (the same 60% threshold used
  everywhere else in this project) **and** it retains ≥ 8 genes after filtering.
- Filtered programmes are named `<P>_epi` and are **distinct objects** from the originals. The
  originals are not edited; v1 results stay attached to v1 definitions.
- Composition-proxy genes remain forbidden, as asserted in `mvri_programs.py`.

Applying this: **P4 fails the 60% rule and is dropped entirely** (6/10 epithelial = 60% but only
6 genes remain, below the size floor of 8). P6 retains 8 of 14. All others retain ≥ 10.

## 3. What is re-tested, and how it is reported

1. **Outcome meta-analysis** across RSV (GSE155237), pneumococcal (GSE261277) and NOVA
   (GSE253686), per-SD logistic, composition-adjusted, fixed effect.
2. **LAIV engagement**, paired baseline → day 3, Thwaites GSE230494, with M05 as positive
   control.

**Both are re-analyses of data already analysed.** Therefore:

- Every result is reported **alongside the unfiltered result**, never instead of it.
- The matched null and label permutation are re-run; the NOVA specificity failure is expected to
  persist and **if it does, the NOVA cohort stays voided** — the filter is not a route to
  rescuing it.
- **No p-value from this analysis may be presented as a pre-registered finding.** The honest
  label is: a sensitivity analysis showing whether removing immune contamination strengthens,
  weakens or leaves unchanged a result obtained earlier.
- Effective number of tests is reported beside BH, per `docs/MODULE_DESIGN_AUDIT.md`.

## 4. The naming rule still stands, unchanged

`PREREG_mvri.md` §6 and `PREREG_celltype_mvri.md` §5 are **not** relaxed by this plan. For the
phrase "mucosal vaccine readiness index" to be used, a programme must predict a **vaccine
outcome**. No cohort in reach provides one:

| Attempt | Endpoint | Result |
|---|---|---|
| Liverpool LAIV vs TIV | carriage, between-arm | null, all 7 |
| Gambia LAIV | IgA, participant link severed | not testable |
| NOVA | shedding × vaccination interaction | null (p = 0.41) |
| Thwaites LAIV | engagement, paired, control passed | null, all 7 |

**Nothing in this plan can change that, and a positive result here would not license the phrase.**
The most this analysis can establish is that a *cleaner epithelial* programme associates with
*infection susceptibility* — which is the claim already on the table, not a vaccine claim.

## 5. Pre-specified verdicts

| Observed | Reading |
|---|---|
| filtered effect stronger, NOVA still fails matched null | contamination was diluting a real epithelial signal; NOVA still voided |
| filtered ≈ unfiltered | the immune genes were never driving it; the result is robust |
| filtered weaker | the earlier result depended partly on immune infiltration; stated plainly |
| P5 unchanged | expected — P5 has **zero** genes below 0.5 and should be unaffected. If P5 moves materially, the filter or the scoring is wrong and that is investigated before anything is claimed |

P5 is the built-in negative control for the filter itself.
