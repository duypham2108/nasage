# Specification curve and axis-orthogonal redesign

Status: **LOCKED** by the commit that adds this file.
Follows `docs/WEAKNESSES.md`, which identifies the only three attackable weaknesses.

---

## 0. Why a specification curve and not "trying approaches"

Running many analyses and reporting the ones that work is the failure mode this project exists
to avoid. Running **every** reasonable analysis and reporting the **whole distribution** is a
recognised method — specification curve / multiverse analysis (Simonsohn, Simmons & Nelson 2020,
*Nat Hum Behav* 4:1208, [doi:10.1038/s41562-020-0912-z](https://doi.org/10.1038/s41562-020-0912-z);
Steegen et al. 2016, *Perspect Psychol Sci* 11:702).

The difference is entirely in what is reported:

- **The headline is the distribution**: the median effect across all specifications, the
  fraction of specifications reaching p < 0.05, and the sign stability.
- **The maximum is never the headline.** The best specification is reported only alongside how
  many were run to find it.
- **No specification may be added after seeing results**, and none removed. The grid below is
  complete and fixed.

## 1. The grid, fixed now

Varied jointly, giving **4 × 4 × 3 = 48 specifications** per programme per outcome:

**Scoring (4):** mean of z-scored genes · median of z-scored genes · mean reciprocal rank
(UCell-style, depth-robust) · first principal component of the gene set.

**Covariates (4):** none · composition proxy · composition proxy + depth · composition proxy +
depth + global PC1.

**Gene set (3):** original programme · epithelial-filtered (`PREREG_epithelial_only.md`) ·
**axis-orthogonal redesign** (§2).

**Applied to:** the three challenge outcomes (RSV, pneumococcal, NOVA) for P5 and the other
epithelial programmes, and the lifespan age trajectories.

**Unit is always the donor.** Cells are never the unit; that is not a researcher degree of
freedom, it is a correctness requirement.

## 2. The axis-orthogonal redesign

**Rationale.** `WEAKNESSES.md` B1: everything that fails the matched null fails because a single
dominant axis is shared by almost all gene sets. Correcting afterwards (Approach A) destroyed the
signal. This instead builds sets that do not load on the axis in the first place.

**Procedure, fixed:**

1. In **Ziegler 2021 + Loske 2022 single-cell data only** (neither contributes an outcome),
   compute the donor × cell-type score matrix for all programmes and extract **PC1 per cell
   type**.
2. For every gene in a programme, compute its correlation with that PC1 across donor × cell-type
   units.
3. **Retain genes in the lowest two tertiles of |correlation with PC1|** that are also
   **epithelial** (fraction ≥ 0.5) — the filter already locked.
4. A redesigned programme is kept only if **≥ 8 genes** remain, the same size floor as before.
5. The redesigned set is named `<P>_orth` and is a **separate object**. Originals are untouched.

**No outcome variable is read at any step.** Selection uses Ziegler/Loske expression structure;
testing uses RSV/pneumococcal/NOVA outcomes. Selection and test cohorts are disjoint.

**Pre-declared failure modes:**

| Observed | Reading |
|---|---|
| orthogonal set keeps the association **and** now clears the matched null | the association is real and was previously confounded by the shared axis — **hypothesis for a new cohort**, not a confirmation |
| orthogonal set loses the association | **the association was the shared axis.** This is the most likely outcome and is reported as the main result if it occurs |
| orthogonal set keeps the association but still fails the matched null | the axis was not the cause; the non-specificity is something else |
| fewer than 8 genes survive | that programme cannot be made orthogonal and is reported as such |

## 3. Reporting rules

1. **The full specification curve is plotted and tabulated**, all 48 points per test, ordered by
   effect size, with significance marked.
2. **Three numbers lead**: median effect, % of specifications with p < 0.05, and % with the same
   sign.
3. **A result is called robust only if ≥ 90% of specifications agree in sign and the median
   reaches p < 0.05.** Fixed now.
4. **Multiplicity:** BH over the union of all tests run under this plan, as in
   `PREREG_method_improvements.md` §0.3.
5. **Everything is hypothesis-generating.** These are the same cohorts already analysed. Nothing
   here can confirm anything (`WEAKNESSES.md` E5).

## 4. What is still not attempted, and why

- **No new power.** C1 stands: n = 39/39/74. No specification creates participants.
- **No vaccine outcome.** A1 stands and is untouched by all of this.
- **No gene chosen for how it performed against an outcome**, in any specification.
- **"Mucosal vaccine readiness" remains forbidden**, whatever the curve shows.

## 5. Honest prior, recorded before running

Approach A already showed that removing the shared axis removes most of the signal. The most
likely result of §2 is therefore that **the orthogonal programmes lose the association**. If that
happens it is the single most informative finding of this whole line of work, and it will be
written as the headline rather than buried: *the baseline mucin-barrier association with
infection is substantially a general expression axis, not a specific programme.*
