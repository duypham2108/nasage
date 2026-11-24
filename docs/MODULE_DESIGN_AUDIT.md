# Audit of the module panel, and how to optimise it

Date: 2026-11-24. Script: `validation/module_design_audit.py`. Source data:
`results/module_audit/`.

This is a diagnostic on the **module panel itself**. It fits no outcome model and reads no
outcome variable, so running it cannot leak outcome information into module design. That
property is what makes it safe to act on; see §5, which is the part that constrains everything
else.

---

## 1. Four defects, quantified

### 1.1 The redundancy rule in use catches none of the actual redundancy

The project's rule is **Jaccard > 0.5 on gene membership → drop**. Applied to the current panel
it flags exactly two pairs: P2/M01 (0.824, already acted on) and P7/M02 (0.462, retained).

Measured on donor × cell-type scores in Ziegler 2021 + Loske 2022, **24 module pairs exceed
|ρ| = 0.70 — and all 24 have Jaccard 0.0.**

| pair | ciliated | secretory | basal |
|---|---|---|---|
| M06 MHC-I ↔ M07 MHC-II | **+0.93** | **+0.89** | **+0.81** |
| M05 ISG ↔ M07 MHC-II | +0.81 | **+0.91** | — |
| M05 ISG ↔ M06 MHC-I | +0.78 | +0.84 | — |
| M06 MHC-I ↔ M09 plasma secretory | +0.88 | +0.80 | +0.76 |
| M07 MHC-II ↔ M03 PIGR | +0.74 | +0.83 | +0.75 |
| M01 barrier ↔ M06 MHC-I | — | +0.80 | +0.72 |
| M08 APC costim ↔ M09 plasma secretory | — | — | **−0.83** |

These modules share **no genes**. They are redundant because MHC-I, MHC-II, PIGR and the plasma
secretory machinery are all interferon-inducible: one biological axis measured five ways.

**The rule is not merely weak, it is pointed at the wrong space.** Gene overlap measures how a
module was written; score correlation measures what it detects.

### 1.2 The panel has roughly half the dimensions it claims

| cell type | donors | modules | effective tests (Li & Ji) | (Galwey) | PC1 % of variance |
|---|---|---|---|---|---|
| ciliated | 55 | 14 | **7.0** | 9.0 | **43%** |
| secretory | 54 | 14 | **8.0** | 8.9 | **47%** |
| basal | 39 | 14 | **10.0** | 9.8 | **41%** |

Fourteen modules behave as **7–10 independent tests**, and a single component carries
**41–47%** of all variance across the panel.

Two consequences, pulling in opposite directions:

- **BH over 14 modules over-corrects**, because the tests are positively correlated. Declared
  families are more conservative than intended.
- **Agreement across modules is not replication.** "M05, M06 and M07 all move with age" is close
  to one observation reported three times. Any passage in `RESULTS.md` that reads as mutual
  corroboration between modules on this axis is weaker than it looks.

### 1.3 Two "modules" are single genes, and five more are small

| set | genes |
|---|---|
| M03 PIGR antibody transport | **1** |
| M03b FCGRT IgG transport | **1** |
| M02b airway mucin panel | 5 |
| M11 senescence arrest | 5 |
| M10 IgA fraction, M12 SASP, M13 TRM core | 7 |

A one-gene module is a gene. Calling it a module implies averaging-away of measurement noise
that is not happening, and it is the reason M03 shows up correlated with everything interferon-
driven in §1.1 — a single transcript has no internal averaging to stabilise it.

### 1.4 Measurability and depth, already known, restated as design faults

- **M12 inflammatory SASP fails the 60% gate in all four single-cell cohorts** (0–14% of genes
  detected). It is a design error to carry a module that droplet scRNA-seq cannot measure.
- **M11 senescence arrest** is 5 genes and is entangled with sequencing depth (Supp. Fig. S5b:
  CDKN2A⁺ status is largely a depth event).
- Scoring is **inconsistent**: UCell (rank-based, depth-robust) in single cell, mean-z of
  log2 CPM in bulk. The depth confound lives in the mean-z path.

## 2. What this costs the existing results

Nothing here invalidates a computed number. It changes what the numbers mean:

1. **Module-to-module agreement is weak evidence**, not strong. Must be restated wherever used.
2. **Family-wise corrections were conservative.** No claim becomes significant that was not;
   some nulls were tested more harshly than declared.
3. **It partly explains the NOVA matched-null failure.** If 41–47% of module-score variance is
   one axis, then in bulk tissue a random expression-matched gene set will load on that axis
   too — which is precisely the 16.4% of random sets reaching p < 0.05
   (`results/nova/` §3). The specificity failure is not only a property of that cohort; it is
   partly a property of scoring any modest gene set against a strongly structured matrix.

## 3. Proposed changes, in priority order

| # | Change | Why | Cost |
|---|---|---|---|
| 1 | **Replace the Jaccard rule with a score-correlation rule**: drop or merge at \|ρ\| > 0.75 measured in held-out single-cell data, keeping Jaccard as a secondary check | catches 24 real redundancies that Jaccard misses entirely | low — the diagnostic already exists |
| 2 | **Declare one `IFN_response` axis** covering M05/M06/M07/M09/M03 and test it **once**, with the members as pre-declared components of it rather than as independent modules | they are one dimension; §1.1 | medium — rewrites a family definition |
| 3 | **Report the effective number of tests** alongside every BH correction | makes the correction honest in both directions | low |
| 4 | **Retire M12 SASP** from all primary families; keep it only as the measurability demonstration it already serves as in Supp. Fig. S5c | cannot be measured on the platform | low — already carries no claim |
| 5 | **Rename M03 and M03b as genes, not modules** (`PIGR_expression`, `FCGRT_expression`) | one gene is not a module | low — naming only |
| 6 | **Move bulk scoring to a rank-based statistic** matching the single-cell path | removes the depth confound at the source; makes bulk and single-cell results commensurable | medium — re-runs bulk scoring |
| 7 | **Size-matched nulls**: the matched null should match set size as well as expression bin | a 5-gene set and a 20-gene set have different null spreads | low |
| 8 | **Minimum module size of 8**, with smaller sets either merged or explicitly labelled single-gene readouts | §1.3 | medium |

## 4. What is deliberately **not** proposed

- **No gene is added or removed on the basis of how it performed in an outcome test.** That is
  the one change that would be genuinely fatal, and it is not on the list.
- **No re-derivation from data.** The unsupervised attempt failed and is recorded in
  `DERIVATION_v1_FAILED.md`; curation with references is the current design and stays.
- **No retuning of the 60% measurability gate** to let M12 through.

## 5. The constraint that governs all of this

**Module optimisation using the cohorts that produced the results invalidates those results.**
Every outcome test in this project was pre-registered against a module panel fixed in advance.
If the panel is re-fitted — even using only expression structure — and the outcome tests are
re-run on the same cohorts, the pre-registration no longer holds and the p-values become
uninterpretable.

The only safe sequence:

1. Freeze the current panel as **v1**. Every result already computed stays attached to v1 and is
   reported as such. Nothing is recomputed to look better.
2. Build **v2** using only the structural diagnostics above — expression correlation, size,
   measurability — and **never an outcome variable**.
3. Lock v2 by commit, with the derivation script.
4. Test v2 **only on a cohort not used to build it**, under a new pre-registration.
5. Report v1 and v2 side by side. If v2 is better, that is a finding; if it is not, that is also
   a finding and it gets reported.

The audit in §1 satisfies step 2's constraint: it read scores and gene lists, never an outcome.

## 6. Recommendation

Do **1, 3, 4, 5 and 7 now** — they are bookkeeping and honesty fixes that cost nothing and
change no computed result. Specifically: report effective tests next to every BH, rename the two
single-gene modules, formally retire M12 from primary families, add size matching to the null,
and swap the redundancy rule.

Treat **2 and 6** as a v2 panel under §5 — they change what gets tested and must not be applied
retrospectively.

The single highest-value item is **2**. Collapsing the interferon/antigen-presentation axis into
one declared dimension would remove the largest source of false reassurance in the current
design, and it is the change a reviewer is most likely to demand once they see the correlation
matrix.
