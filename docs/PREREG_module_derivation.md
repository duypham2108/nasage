# Derivation plan for outcome-anchored NasAge modules (Set A)

Status: **LOCKED** by this commit. Later changes are dated deviations recorded in the results file, never edits here. Nothing in
the derivation datasets has been examined beyond GEO sample metadata and published abstracts.

## 0. Objective and what it can claim

Derive a gene module that captures **the nasal mucosal programme engaged by a live attenuated
intranasal vaccine (LAIV)**, defined entirely in public cohorts, then test how that programme
varies across age bins in our own cohort.

| Can claim | Cannot claim |
|---|---|
| "This is the mucosal programme a live intranasal vaccine engages, replicated across independent cohorts" | "A high score means a better vaccine response" |
| "Its baseline level differs / does not differ by age" | Anything about IgA or protection: no immunogenicity endpoint enters the derivation |
| Direction-specific statements per pathogen | A single 'mucosal readiness' score |

The anchor is **vaccine-induced mucosal response**, not vaccine immunogenicity. Set B
(outcome-associated modules, anchored on carriage or infection susceptibility) is deferred; its
predicted signs are already recorded in §8.2 of `PREREG_vaccine_axes.md`.

## 1. Derivation data (public only; our cohort is NEVER used to derive)

| Role | Dataset | Design | Usable contrast |
|---|---|---|---|
| **Discovery** | GSE230494 (Thwaites 2023) | 40 young adults, LAIV, **no bacterial challenge**; nasal curettage at baseline and 72 h | **Paired within participant**, baseline → 72 h (n = 20–21 pairs) |
| **Replication 1** | GSE164649 (Carniel 2021) | Randomised LAIV vs TIV; pneumococcus day 0, **vaccine day +3**; baseline, D2, D6, D9 | **LAIV vs TIV at D6 and D9** (= 3 and 6 days post-vaccine) |
| **Replication 2** | GSE117580 (Jochems 2018) | Randomised LAIV vs TIV; **vaccine 3 days before** pneumococcus; D−4, D2, D9 | **LAIV vs TIV at D2 and D9** (= 5 and 12 days post-vaccine) |

Facts established from metadata before locking:
- Day numbering in both Liverpool studies is relative to **pneumococcal inoculation**, not
  vaccination. Every post-vaccine sample there is also post-challenge, which is why the
  arm contrast, not a time contrast, is used: randomisation balances bacterial exposure.
- GSE164649: 44 subjects, 22 LAIV / 22 TIV at baseline; D6 18/17, D9 16/20. Carries age and sex.
- GSE117580: 35 subjects; D2 18/16, D9 12/15.
- **No participant overlap** between the two Liverpool series (subject IDs differ in the study
  segment, `/019/` vs `/016/`; intersection empty). Replication between them is therefore real,
  though both come from one group and one challenge model.
- All three are bulk RNA-seq on GPL20301/GPL16791-class platforms.

## 2. Gene entry rule

A gene enters Set A only if **all** hold:

1. **Expressed**: mean CPM ≥ 1 in the discovery dataset's samples.
2. **Discovery**: paired differential expression baseline → 72 h in GSE230494, BH-adjusted
   p < 0.05 **and** |log2 fold change| ≥ 0.5, by a paired test on log2(CPM+1) with participant as
   a blocking factor.
3. **Replication**: same direction **and** nominal p < 0.05 in at least one of the four LAIV-vs-TIV
   contrasts (GSE164649 D6, D9; GSE117580 D2, D9), each adjusted for the study's own carriage
   status and sequencing depth.
4. **Size cap**: if more than 100 genes qualify, keep the 100 with the largest discovery effect
   size. Modules larger than this are unwieldy to interpret and to score.

If fewer than 10 genes survive, **Set A fails** and is reported as a failure, not rescued by
relaxing thresholds.

## 3. The composition check (the main safeguard)

Bulk-derived modules often encode **cell-type proportions** rather than cell state. Our package
scores within cell types, which removes exactly that — so such a module would mean something
different, or nothing, after transfer. Today's RSV result makes this a live risk, not a
theoretical one.

For each candidate gene, compute cell-type specificity in an independent reference nasal
single-cell dataset (**Ziegler 2021**, already held locally; not Yoshida):

- Report, for the final module, the fraction of genes whose expression is concentrated in a
  single population (top-decile specificity).
- If **> 50%** of genes concentrate in one population, the module is flagged
  `composition_dominated` and is **reported but not transferred** as a cell-state module.
- The per-cell-type expression profile is published alongside the module either way.

## 4. Relationship to the existing 15 curated modules

The LAIV response is likely to be interferon-dominated, so Set A may largely recapitulate M05.
That must be measured and stated, not glossed:

- Report Jaccard overlap of Set A with each of the 15 curated modules.
- **If overlap with M05 ≥ 0.5**, Set A is reported as *empirical confirmation of M05*, not as a
  new module, and the paper says so plainly.
- Genes in Set A that are **not** in any curated module are listed separately; these are the
  genuinely new content and get their own literature annotation, as in `MODULE_EVIDENCE.md`.

## 5. Transfer to our cohort (fixed now, before any derivation)

- Scoring: the package's normal path — UCell on raw counts, per cell type.
- Primary populations: **basal, secretory, ciliated** (`primary_celltypes`); immune populations
  exploratory, per the 2026-11-24 decision.
- Model: the standard `nas.tl.age_effects` design — age bin + batch + sex + depth + subtype
  fractions; BH within the primary family.
- Negative controls: ≥ 500 expression-matched random gene sets, plus the permutation null.
- **Direction: none is predicted.** Whether the vaccine-response programme rises or falls with
  age is the open question, so the age test is **two-sided and exploratory in direction**, while
  remaining confirmatory in the sense that the gene set is fixed before our data are touched.
  Claiming a predicted direction here would be false: we have no prior strong enough.
- An applicability entry for Set A is assigned from the §3 cell-type profile, not by hand.

## 6. What would make this fail, stated in advance

| Condition | Consequence |
|---|---|
| < 10 genes survive replication | Set A fails; reported as such |
| > 50% of genes concentrate in one cell type | flagged `composition_dominated`, not transferred |
| Jaccard with M05 ≥ 0.5 | reported as confirmation of M05, not a new module |
| Module not measurable in our data (coverage flag) | reported, not re-specified to fit |

## 7. Order of work

1. Lock this plan by commit.
2. Download the two Liverpool count matrices; build sample tables.
3. Discovery in GSE230494 → candidate list.
4. Replication in the four arm contrasts → Set A.
5. Composition check and curated-module overlap.
6. Only then: score Set A on our cohort and run the pre-specified age model.

Deviations are appended, dated, to the results file. Thresholds in §2 and rules in §6 are not
changed after step 3 begins.
