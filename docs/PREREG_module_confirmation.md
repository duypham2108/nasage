# Confirmation and refinement of the curated modules against live mucosal vaccination

Status: **LOCKED** by the commit that adds this file. Later changes are dated deviations in the
results file, never edits here.

Context: this follows the documented failure of the Set A derivation
(`results/setA_derivation/RESULT_attempt1_FAILED.md`, plan `b99961e`). That attempt tried to
*derive* modules per gene from a 20-pair cohort and had no power to do so. The present plan asks
the better-powered question on the same data: **which of the 15 existing curated modules does a
live attenuated intranasal vaccine engage, and which genes inside them carry that response?**

## 0. What this can and cannot claim

| Can claim | Cannot claim |
|---|---|
| "These curated programmes are engaged by a live intranasal vaccine, replicated across independent randomised cohorts" | "A high score predicts a better vaccine response" — no immunogenicity endpoint enters |
| "These genes within a module carry the vaccine response; these do not" | That a refined module is a new discovery: the gene universe is the curated module |
| "The engaged programmes do / do not differ across age bins" | A single mucosal 'readiness' score |

## 1. Data and contrasts (unchanged from `PREREG_module_derivation.md` §1)

| Role | Dataset | Contrast |
|---|---|---|
| Discovery | GSE230494 | paired baseline → 72 h post-LAIV, no bacterial challenge, 20 pairs |
| Replication 1 | GSE164649 | LAIV vs TIV at D6, D9 (3 and 6 days post-vaccine) |
| Replication 2 | GSE117580 | LAIV vs TIV at D2, D9 (5 and 12 days post-vaccine) |

Day numbering in both Liverpool studies is relative to pneumococcal inoculation, so every
post-vaccine sample there is also post-challenge; the randomised **arm** contrast is used because
it balances bacterial exposure. The two Liverpool series share no participants. In both, TIV is a
systemically administered comparator, so LAIV-vs-TIV at the nasal mucosa isolates the live
intranasal response. Our own cohort is never used here.

## 2. Primary analysis: which modules are engaged

Unit of analysis is the **module score**, not the gene. All 15 curated modules are tested;
multi-gene modules and the single-gene features (M03, M03b) and the M10 ratio all enter as the
package defines them.

- **Scoring.** Two methods are computed: (a) mean over module genes of z-scored log2(CPM+1),
  z across that dataset's analysed samples; (b) UCell on raw counts. Genes with mean CPM < 1 in
  the dataset are dropped. Samples with < 5 M assigned counts are excluded.
- **Discovery test.** Paired test of the module score, baseline vs 72 h, across the 20 paired
  participants. BH across the 15 modules.
- **Replication test.** Module score ~ arm + carriage + log10 depth in each of the four arm
  contrasts.
- **A module is ENGAGED only if all hold:**
  1. discovery BH < 0.05,
  2. same direction and nominal p < 0.05 in **at least one** arm contrast,
  3. **the same verdict under both scoring methods.**

Condition 3 is deliberate. We found earlier that mean-z and UCell can give materially different
answers on bulk data (M12: OR 1.93 vs 4.60 in GSE155237), so a module whose status depends on the
scoring choice is not robust enough to carry a claim, and is reported as `scoring_dependent`.

Modules that are **not** engaged are reported too: that a programme is untouched by a live
mucosal vaccine at these timepoints is itself informative (for example, if M03 PIGR does not
move, the IgA transport machinery is not transcriptionally engaged by 72 h).

## 3. Secondary analysis: which genes carry the response

Only for modules confirmed as ENGAGED, and only within that module's own curated gene list — the
universe is fixed by the literature, so this is refinement, not discovery.

- For each gene: paired log2FC and p in discovery; same model as §2 in each arm contrast.
- A gene is a **core responder** if: discovery BH < 0.05 (BH computed across the genes of all
  engaged modules only, which is a universe of order 100, not 15,000) **and** same direction with
  nominal p < 0.05 in at least one arm contrast.
- For each engaged module with **>= 5** core responders, a refined variant `<module>_core` is
  defined. Fewer than 5 means no variant is created for that module.
- Report, per engaged module, the proportion of its genes that are core responders. A module
  where most genes respond is a coherent programme; one where two genes carry everything is not,
  and is labelled `driven_by_few`.

## 4. Composition check (unchanged in substance from `PREREG_module_derivation.md` §3)

All three datasets are bulk, so an apparent response can reflect immune-cell infiltration rather
than a change of cell state — and our package scores within cell types, which removes exactly
that. For every engaged module and every `_core` variant, compute per-gene cell-type specificity
in an independent reference nasal single-cell dataset (Ziegler 2021; not Yoshida):

- Report the fraction of genes concentrated in a single population.
- If **> 50%** concentrate in one population, flag `composition_dominated`: reported, and **not**
  transferred as a cell-state module.
- Publish the per-cell-type profile either way.

This limitation is stated explicitly in any write-up: the discovery dataset has no cell-type
resolution, so "engaged" at bulk level cannot distinguish infiltration from state change.

## 5. Transfer to our cohort

Fixed before our data are touched.

- Score with the package's normal path: UCell on raw counts, per cell type.
- Primary populations **basal, secretory, ciliated**; immune populations exploratory.
- Both the **full curated module** and its `_core` variant are tested, reported side by side.
- Model: standard `nas.tl.age_effects` (age bin + batch + sex + depth + subtype fractions), BH
  within the primary family; >= 500 expression-matched random gene sets and the permutation null.
- **Direction: none predicted.** Whether a vaccine-engaged programme rises or falls with age is
  the open question. The age test is two-sided; the gene sets are fixed beforehand, which is what
  removes circularity.

## 6. Failure conditions, stated in advance

| Condition | Consequence |
|---|---|
| No module meets all three ENGAGED criteria | reported as a failure; no transfer |
| A module passes under one scoring method only | `scoring_dependent`; reported, not claimed |
| < 5 core responders in an engaged module | no `_core` variant for it |
| > 50% of genes in one cell type | `composition_dominated`; not transferred |
| A module is not measurable in our cohort | reported via the coverage flag, not re-specified |

## 7. Order of work

1. Lock this plan by commit.
2. Module-level discovery and replication (§2).
3. Gene-level refinement for engaged modules only (§3).
4. Composition check and write-up (§4).
5. Only then: transfer to our cohort and run the pre-specified age model (§5).
