# Addendum: do mucosal vaccines engage the epithelial programs?

Status: **LOCKED** by the commit that adds this file. Addendum to `docs/PREREG_mvri.md`
(commit `fecab57`), which it does not override: that plan's outcome test is unchanged and
remains unrun.

Programs are fixed in `docs/MVRI_PROGRAMS.md` / `validation/mvri_programs.py`, with a reason per
gene and references per program.

---

## 0. Why this experiment exists, and what it is not

The locked plan tests whether a baseline program predicts an **infection** outcome. It never
asks the prior question: **is this program touched by a mucosal vaccine at all?** A readiness
axis that a live intranasal vaccine does not perturb may still be a valid baseline predictor,
but one that *is* perturbed, specifically by the mucosal route, has a mechanism behind it.

Two datasets make this testable, and they are unusually good for it because participants were
**randomised between a mucosal and a systemic influenza vaccine**:

| | GSE117580 | GSE164649 |
|---|---|---|
| samples | 96 nasal | 151 nasal |
| arms | LAIV 48 / TIV 48 | LAIV 73 / TIV 78 |
| timepoints | D−4 (pre-vaccine), D2, D9 | baseline, D2, D6, D9 |
| carriage | POS 50 / NEG 46 | POS 84 / NEG 67 |
| subject IDs | present | present |
| source | Liverpool EHPC pneumococcal challenge programme | Liverpool EHPC pneumococcal challenge programme |

**The TIV arm is the control that makes this worth doing.** Both arms receive an influenza
vaccine; only one receives it on the nasal mucosa. A change seen in LAIV but not TIV is
attributable to the mucosal route rather than to vaccination, season, sampling or the
challenge procedure.

This experiment uses **post-vaccination expression but no challenge outcome**. Running it does
not consume the outcome test in `PREREG_mvri.md`.

## 0.1 A declared non-independence

GSE117580, GSE164649 and GSE261277 all come from the **Liverpool experimental human
pneumococcal carriage programme**, and participant identifiers share a format
(`1060/016/045`, `1480/019/118`). Public metadata does not let me exclude participant overlap
between them.

Consequence, fixed now: the pneumococcal carriage endpoint is treated as **one cohort family,
not three independent replications**. Any carriage result from GSE117580 or GSE164649 is
reported as a *consistency check* on GSE261277, never as independent replication, and is not
pooled with it. The RSV cohort (GSE155237, Imperial) is the only genuinely independent
challenge cohort.

---

## 1. Primary question and model

**Does a mucosal vaccine change each program's score more than a systemic vaccine does?**

Per dataset, per program:

```
score ~ time * arm + sex + log10_depth,  with subject as a random intercept
```

`time` = post timepoint vs baseline; `arm` = LAIV vs TIV. The estimand is the
**time × arm interaction**: the difference between arms in change from baseline.

- Subjects with a baseline and at least one post sample enter. Participants contribute paired
  observations, so a mixed model is used; where it fails to converge, the paired
  baseline-to-post delta is compared between arms by Mann–Whitney, declared now as the fallback.
- **Two-sided for every program.** `MVRI_PROGRAMS.md` records which directions would be
  unsurprising; none is a locked prediction.
- Timepoints are tested separately (D2 and D9 in GSE117580; D2, D6, D9 in GSE164649), and the
  family for multiplicity is **programs × timepoints within a dataset**.

### 1.1 Secondary: does the mucosal vaccine move the program at all?

Within the LAIV arm only, the paired baseline → post change, Wilcoxon signed rank. Reported
alongside the interaction so that "LAIV moves it, TIV does not" can be distinguished from "TIV
moves it downwards".

### 1.2 Replication requirement

A program is called **engaged by the mucosal route** only if the interaction is same-signed in
**both** datasets and reaches BH < 0.05 in at least one. One dataset alone is reported as
`single-cohort`, with both estimates shown.

## 2. Scoring and the composition problem

Identical to `PREREG_mvri.md` §3 and not re-specified here:

- mean z of log2(CPM+1) over program genes at CPM ≥ 1; a program with < 60% of genes available
  in a dataset is **not tested** there and is reported as untransferable;
- the **primary model adjusts for the program's own cell-type composition proxy**, built from
  the marker genes that no program may contain. A vaccine that changes ciliated-cell abundance
  would otherwise move every ciliated program;
- the unadjusted model is a **sensitivity**, reported beside the primary, never instead of it.

## 3. Specificity control

**500 expression-matched random gene sets** per program per dataset, 25 expression bins, same
model, same interaction term. Empirical p is the fraction of random sets with an interaction at
least as extreme. A program that does not clear this is `not_specific` and carries no claim.

Required here for a specific reason: LAIV is a live virus and provokes a broad interferon
response, so a great many gene sets will move after LAIV and not after TIV. Without the matched
null, "engaged by the mucosal route" would be close to meaningless.

## 4. Power, declared before any p-value

GSE117580: 35 baseline, 34 D2, 27 D9 across two arms — roughly 15–18 subjects per arm per
contrast. GSE164649: 44 baseline, 35–36 per post timepoint — roughly 18–22 per arm.

For a between-arm comparison of paired deltas at n ≈ 17 per arm, two-sided α = 0.05:

| Cohen's d of the interaction | power |
|---|---|
| 0.5 | 0.29 |
| 0.8 | 0.60 |
| 1.0 | 0.78 |
| 1.2 | 0.90 |

After BH over a family of ~14–21 (7 programs × 2–3 timepoints), **only large interactions
(d ≳ 1.0) are reliably detectable**. A null is reported with this table. The requirement in §1.2
that both datasets agree in sign is the main protection against a chance finding, not the
p-value.

## 5. Failure conditions

| Condition | Consequence |
|---|---|
| A program is untransferable in a dataset | excluded there, reported, not re-specified |
| No program reaches BH < 0.05 | reported as a null with the power table; no re-specification |
| BH passed but matched null failed | `not_specific`; no claim |
| Significant unadjusted, null after composition adjustment | reported as **composition change, not program activity** |
| Signs disagree between the two datasets | `inconsistent`; both shown, no claim |
| Mixed model fails to converge | the §1 fallback is used and the deviation is dated in the results file |

## 6. What this cannot claim

- **Not protection.** Engagement by LAIV is not evidence that a program predicts a vaccine
  outcome. That question belongs to `PREREG_mvri.md` and is still unrun.
- **Not cell-type specificity in these data.** Programs were curated as cell-type programs and
  are scored in bulk nasal tissue.
- **Not generalisation beyond influenza vaccines.** Both datasets use LAIV against TIV.
- Nothing about children: GSE164649 is 17–27 years old, GSE117580 adult.

## 7. Order of work

1. Lock this addendum and the program definitions by commit.
2. Report transferability (gene availability) per dataset **before** any model is fitted.
3. Primary interaction models, composition-adjusted, with the matched-null control.
4. Secondary within-LAIV change, and the unadjusted sensitivity.
5. Apply §1.2 and state which programs, if any, are engaged by the mucosal route.
