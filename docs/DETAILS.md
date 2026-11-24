# nasage — full description

Age-associated **nasal mucosal immune-state profiling** from single-cell RNA-seq
of nasal brushings: pre-registered gene modules, a module-cell-type
applicability matrix, UCell scoring, donor-level inference, and a simulator for
checking that your design can detect what you are looking for.

It answers: *do donors in different age bins differ in these mucosal programmes,
within the same annotated cell population?*

It does **not** answer: *which age group should be vaccinated?* — see
[docs/INTERPRETATION.md](INTERPRETATION.md). No composite "vaccine
readiness" score is produced, by design.

---

## Install

```bash
conda env create -f environment.yml && conda activate nasage
pip install -e .
pytest -q                 # end-to-end smoke test on simulated data
```

or, in an existing environment: `pip install -r requirements.txt && pip install -e .`

## Quickstart

Use the Python functions shown in the [README](../README.md#3-run-it-in-python).

## What your data needs

| Requirement | Why |
|---|---|
| `.h5ad` with **raw integer counts** (in `.X` or a named layer) | QC metrics, single-gene pseudobulk and the IgA ratio are computed from counts |
| `donor` column | the unit of replication; cells are not replicates |
| `cell_type` annotation | all primary comparisons are within a population |
| `age` (years) or `age_bin` | the exposure |
| `batch`, `sex` | adjusted for; a warning is raised if missing |
| ≥ 3 donors per age bin, ideally ≥ 8 | below this the age term is not estimable |

The population column (`celltype_key`) must use the canonical names (basal,
ciliated, secretory, dendritic, macrophage, Tcell, Bcell, plasma, and, mostly for
composition, neutrophil, squamous, ionocyte, tuft, mast, NK); cells with `NaN`
are left out. Pass your finer annotation as `subtype_key`: when several labels
share a population (goblet and club cells -> secretory), each donor's subtype
fractions enter the model as covariates, so a change in subtype mix is not read
as a change in cell state.

## What the pipeline does

1. **Check** — refuses non-integer counts; `nas.pp.check_design` stops if age,
   batch or sex vary within a donor, and flags batch-age confounding and age
   groups with too few donors.
2. **QC** (`nas.pp.filter_cells`) — scanpy QC metrics on the counts layer; prespecified
   thresholds, each settable per cell type; adaptive mitochondrial cut-off
   (`max(fixed, median + 3 MAD)` within each cell type) so mitochondria-rich
   ciliated cells are not removed preferentially; optional doublet column.
3. **Symbol harmonisation** — exact match, then a verified alias table
   (`DDX58→RIGI`, `PIFO→CIMAP3`, `MB21D1→CGAS`, `TMEM173→STING1`, each checked
   against mygene.info), optionally mygene.info for the remainder. Reports what
   is missing per module instead of silently shrinking a module.
4. **Coverage and detectability** — per module x cell type; poorly measured
   combinations are **flagged and excluded from testing**, never interpreted as
   biological absence.
5. **Scoring** — [UCell](https://doi.org/10.1016/j.csbj.2021.06.043) via the
   published `pyUCell` implementation (`max_rank=1500`). `PIGR` and `FCGRT` are
   **single-gene features**, reported as donor-level pseudobulk
   `log1p(mean CP10K)` — they are not padded into artificial signatures. The
   IgA fraction (M10) is a **ratio feature**: logit of pooled IgA over IgA+IgG+IgM
   heavy-chain counts. `score_distribution` (in `adata.uns`) flags modules whose UCell
   score sits near its ceiling or is zero in most cells.
6. **Donor x cell-type aggregation** — minimum 10 cells by default; below the
   threshold the value is **missing, not zero**. `statistic="auto"` uses the
   median, or the mean where more than half of a population's cells score zero
   (decided per module x cell type, blind to age). The summary also carries the
   donor's median detected genes (`log_genes`) and subtype fractions.
7. **Inference** — per (module, cell type): OLS of the donor-level summary on
   age bin + batch + sex + `log_genes` + subtype fractions (the last two dropped
   when the donor count cannot support them, and recorded); omnibus Wald test for age plus contrasts against a
   reference bin, with 95% CIs, donor counts and donors-per-bin. BH-FDR within
   one declared family (the primary combinations of the applicability matrix, restricted to
   basal, ciliated and secretory cells by `primary_celltypes`; immune populations are
   exploratory because they are sparse in many donors, especially older ones);
   exploratory combinations are adjusted separately and never promoted.
8. **Composition** — CLR regression of donor proportions (arcsine-root as a
   robustness check), BH across populations.
9. **Calibration** — `nas.tl.permutation_null` permutes age labels across donors
   (within batch) to check the realised false-positive rate of *your* family on
   *your* data.
10. **Biological checks** — `nas.tl.module_correlation` (do M04, M05 and M06 move
    together across donors?) and a sensitivity rerun without HLA genes
    (`nas.drop_hla`).
11. **Report** — `nas.write_report` writes `report.md` plus the main tables; all
    results are also kept in `adata.uns["nasage"]`.

## Gene modules

Fifteen features, kept deliberately separate (`src/nasage/resources/gene_modules.json`,
provenance in `src/nasage/resources/gene_set_provenance_modules.csv`, references with
PMIDs/DOIs in `src/nasage/resources/references.csv`; the evidence for
each module is summarised in [MODULE_EVIDENCE.md](MODULE_EVIDENCE.md)):

| Module | Content | Reading for mucosal vaccination |
|---|---|---|
| M01 epithelial barrier | claudins, occludin, ZO scaffolds, junctional adhesion | context-dependent |
| M02 mucociliary axoneme | dynein arms, radial spokes, FOXJ1 programme | context-dependent |
| M02b airway mucin panel | MUC5AC/MUC5B/SCGB1A1/BPIFA1/BPIFB1 (small panel) | context-dependent |
| M03 PIGR transport | **single gene** | supportive (capacity only) |
| M03b FCGRT transport | **single gene** | context-dependent |
| M04 innate sensing | constitutive PRRs and adaptors; sensors that are themselves ISGs (RIGI, IFIH1, DHX58, TLR3, AIM2, CGAS) removed | context-dependent |
| M05 interferon-stimulated genes | ISG effectors; type I and III interferon induce the same genes | context-dependent, **not** automatically beneficial |
| M06 MHC-I processing | immunoproteasome, TAP, peptide loading (interferon-inducible) | context-dependent |
| M07 MHC-II presentation | class II, CD74, CIITA, CTSS | supportive in professional APC; exploratory in epithelium |
| M08 APC costimulation | B7/TNFSF ligands, mature-DC axis | supportive |
| M09 plasma secretory machinery | PRDM1/XBP1/IRF4 and ER secretory genes; no Ig genes | supportive |
| M10 IgA fraction | **ratio**: IgA share of IgA+IgG+IgM heavy-chain counts | supportive |
| M11 senescence arrest | CDKN1A, CDKN2A, SERPINE1, IGFBP7, GDF15 | context-dependent |
| M12 inflammatory SASP | IL6, IL1A/B, CXCL8, CXCL1, CCL2, TNF | context-dependent |
| M13 TRM core | CD69, ITGAE, ITGA1, CXCR6, PDCD1, CRTAM, DUSP6 | supportive |

Sensing is separate from interferon response, presentation from costimulation,
and secretory capacity from isotype, because merging them hides which part of
the pathway differs. Interferon still drives M04, M05 and M06 partly together;
check `nas.tl.module_correlation` before reporting them as separate findings.

## Designing a study with this

`nasage.simulate.simulate` generates synthetic data at your
intended donor count (`n_per_bin`), cells per donor and sequencing depth
(`libsize_scale`), with planted effects (`"planted"`) or none (`"null"`). Run the
pipeline on repeated simulations to estimate power and the realised
false-positive rate.
Published power for 8 donors per bin on these modules is **not** a general
claim — rerun it for your design.

## Outputs

`nas.write_report(adata, outdir)` writes:

```
outdir/
├── report.md                      # plain-language summary with the caveats attached
├── module_coverage.csv            # detectability per module x cell type
├── donor_celltype_eligibility.csv # donors and cells per cell type
├── age_effects_omnibus.csv        # F, p, BH-FDR, donors per bin, estimability
├── age_effects_contrasts.csv      # estimates, 95% CI, standardised effects
├── composition_results.csv        # CLR and arcsine-root
├── permutation_null.csv           # realised false-positive behaviour (if run)
├── donor_celltype_scores.csv      # the donor-level table every test is based on
├── donor_design.csv, qc_report.csv, score_distribution.csv, symbol_mapping.csv
└── pooled_vs_within.csv           # composition-confounding check (if given)
```

The same tables are in `adata.uns["nasage"]`; figures come from `nas.pl.*`.

## Provenance and honesty notes

* Every gene symbol was checked against mygene.info; Reactome (release 97) and
  GO identifiers in the provenance table were retrieved from the live services,
  and literature references carry PMIDs resolved through Europe PMC.
* Reactome covers the axonemal module poorly (2/20 genes) — recorded as an
  evidence limitation rather than papered over.
* There is no Reactome pathway for pIgR transcytosis containing `PIGR`; the
  anchor is GO:0002415 (IDA, PMID 10229845).
* The package was developed and validated on **synthetic data only**. Recovery
  of planted effects demonstrates behaviour under the simulator's assumptions;
  it is not biological validation.

## Citing the methods used

UCell (Andreatta & Carmona 2021, PMID 34285779); pyUCell (PMID 41665985);
Scanpy (PMID 29409532); donor-level aggregation rationale (Squair et al. 2021,
PMID 34584091; Crowell et al. 2020, PMID 33257685); compositional analysis
(Büttner et al. 2021, PMID 34824236; Phipson et al. 2022, PMID 36005887);
Benjamini & Hochberg (1995). Full table: `src/nasage/resources/references.csv`.

## Licence

MIT (see `../LICENSE`).
