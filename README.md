# NasAge

**Nas**al mucosal immune state across **age**, from single-cell RNA-seq.
Install as `nasage`; import as `import nasage as nas`.

Scores nasal mucosal immune programmes (gene modules) in single-cell RNA-seq of
nasal brushings and tests whether they differ between age groups **within each
cell type**. Each donor counts once in the statistics, however many cells it has.

It works on your own AnnData object, scanpy-style: you keep your column names
and pass them as `*_key` arguments, the functions live in `nas.pp`, `nas.tl` and
`nas.pl`, and results are stored in `adata.uns["nasage"]`.

- Full description of the method, modules and outputs: [docs/DETAILS.md](docs/DETAILS.md)
- How to read the results, and how not to: [docs/INTERPRETATION.md](docs/INTERPRETATION.md)
- The literature behind each module: [docs/MODULE_EVIDENCE.md](docs/MODULE_EVIDENCE.md)

---

## 1. Install

```bash
conda create -n nasage python=3.11 pip -c conda-forge
conda activate nasage
pip install -e ".[dev]"
pytest -q            # should report all tests passed
```

## 2. What your AnnData needs

| What | Where | Notes |
|---|---|---|
| **Raw integer counts** | `adata.X` or a layer, e.g. `adata.layers["counts"]` | Normalised values are refused. No normalisation is needed: UCell ranks genes within each cell |
| **Gene symbols** | `adata.var_names` | HGNC symbols (`PIGR`, `HLA-DRA`), not Ensembl IDs |
| **Donor ID** | `obs` column | One ID per donor and sample |
| **Population** | `obs` column | Values from the list below; `NaN` for cells to leave out |
| **Fine annotation** | `obs` column (optional, recommended) | e.g. Goblet, Club, cDC1, mregDC, CD4 T, CD8 T |
| **Age** | `obs` column | Age in years, **or** age groups such as `"18-39"` |
| **Batch, sex** | `obs` columns | Optional but recommended; adjusted for in every model |

Age, batch and sex must be the same for every cell of a donor. If a donor was
sequenced in two batches or sampled twice, give each sample its own donor ID;
otherwise the analysis stops with an error.

**Populations.** Modules are tested only in the populations where they make
biological sense, so the population column must use these names:
`basal`, `ciliated`, `secretory`, `dendritic`, `macrophage`, `Tcell`, `Bcell`,
`plasma`, and, mostly for composition, `neutrophil`, `squamous`, `ionocyte`,
`tuft`, `mast`, `NK`. Build it from your own annotation with pandas:

```python
import scanpy as sc
import nasage as nas

adata = sc.read_h5ad("nasal_brush.h5ad")

adata.obs["population"] = adata.obs["celltype_fine"].map({
    "Basal": "basal", "Suprabasal": "basal",
    "Multiciliated": "ciliated", "Deuterosomal": "ciliated",
    "Goblet": "secretory", "Club": "secretory",
    "cDC1": "dendritic", "cDC2": "dendritic", "mregDC": "dendritic",
    "Macrophage": "macrophage", "Monocyte": "macrophage",
    "CD4 T": "Tcell", "CD8 T": "Tcell", "B cell": "Bcell", "Plasma cell": "plasma",
    "Neutrophil": "neutrophil", "Squamous": "squamous", "NK": "NK",
})                                   # labels not in the dict become NaN and are left out
```

Keep the fine annotation (`celltype_fine` here) and pass it as `subtype_key`:
the model then adjusts for each donor's subtype mix within a population, so a
shift from club to goblet cells is not read as a change in cell state.

## 3. Run the analysis

Choose every setting **before** you look at the age results.

```python
# 1. design check: one row per donor; warns about batch-age confounding and small groups
design = nas.pp.check_design(adata, donor_key="donor_id", celltype_key="population",
                             age_key="age_group", batch_key="chemistry", sex_key="sex")

# 2. cell QC, in place; thresholds can differ per population
qc_report = nas.pp.filter_cells(
    adata, celltype_key="population", layer="counts",
    min_genes={"default": 200, "neutrophil": 100}, min_counts=500,
    max_pct_mito=20, mito_adaptive=True)            # one row per population

# 3. score the 15 modules in every cell (adds columns to adata.obs)
nas.tl.score_modules(adata, layer="counts")

# 4. one value per donor x population: the table every test is based on
summary = nas.tl.donor_summary(
    adata, donor_key="donor_id", celltype_key="population", age_key="age_group",
    batch_key="chemistry", sex_key="sex", subtype_key="celltype_fine", min_cells=10)

# 5. age effects, compositions and a false-positive check
omnibus = nas.tl.age_effects(adata, reference="18-39")
composition = nas.tl.composition(adata)
perm = nas.tl.permutation_null(adata, n_perm=200)

# 6. save
nas.write_report(adata, "results/cohortA", label="cohort A nasal brushings")
adata.write_h5ad("results/cohortA/adata_scored.h5ad")   # results are kept in adata.uns
```

### Settings

| Argument | Meaning |
|---|---|
| `layer` | Layer with raw counts (`None` = `adata.X`) |
| `min_genes`, `min_counts`, `max_pct_mito` | One number, or `{population: value, "default": value}`. Neutrophils carry few transcripts, so a single gene floor removes most of them |
| `mito_adaptive` | Each population's mitochondrial cut-off becomes `max(max_pct_mito, median + 3 MAD)`, so mitochondria-rich ciliated cells are not removed wholesale |
| `min_cells` | A donor with fewer cells of a population is **missing** for it, never zero |
| `reference` | The age group every other group is compared against |
| `statistic` (`donor_summary`) | `"auto"` (default): median per donor, or mean where most cells score zero |

### Age: groups or years

If `age_key` already holds age groups, nothing else is needed. Their order is
`age_order=[...]` if given, otherwise the column's categorical order, otherwise
a number-aware sort (`"5-11"` before `"12-17"`).

If `age_key` holds age in years, add `age_bins` to `check_design` and `donor_summary`:

```python
AGE_BINS = [(5, 11, "5-11"), (12, 17, "12-17"), (18, 39, "18-39"),
            (40, 64, "40-64"), (65, 120, "65+")]       # (min, max, label), inclusive

nas.tl.donor_summary(adata, donor_key="donor_id", celltype_key="population",
                     age_key="age_years", age_bins=AGE_BINS, batch_key="chemistry",
                     sex_key="sex", subtype_key="celltype_fine")
```

Donors outside every bin are left out, with a warning naming them.

### What the model adjusts for

Each test is `score ~ age group + batch + sex + log_genes + subtype fractions`,
fitted on one row per donor. `log_genes` (median detected genes per cell) is
there because UCell scores rise with detected genes, and sample quality can
differ between age groups. The last two terms are dropped automatically when
there are too few donors; `omnibus["covariates"]` and
`omnibus["dropped_covariates"]` record what was used. Switch them off with
`nas.tl.age_effects(adata, reference=..., adjust_depth=False, adjust_subtypes=False)`.

## 4. Read the results

```python
st = adata.uns["nasage"]

findings = omnibus[(omnibus.tier == "primary") & (omnibus.p_adj_BH.astype(float) < 0.05)]
contrasts = st["contrasts"]          # each group vs reference: estimate, 95% CI, std. effect
st["permutation_rates"]              # realised false-positive rates
```

Only **primary** rows with `p_adj_BH < 0.05` are findings. Exploratory rows are
adjusted separately and are never promoted.

| `adata.uns["nasage"]` key | Content |
|---|---|
| `omnibus` | One row per module x population: F test for age, p, `p_adj_BH`, donors per group, `estimable`, covariates used |
| `contrasts` | Each age group vs `reference`: estimate, 95% CI, standardised effect |
| `donor_summary` | The donor-level values every test is based on |
| `coverage` | How well each module is measured in each population; poorly measured pairs are not tested |
| `score_distribution` | Flags modules near the UCell ceiling or zero in most cells |
| `composition`, `proportions` | Population proportions by age group |
| `qc`, `design`, `symbol_report` | QC per population, donor table, gene-symbol matching |

`nas.write_report` writes `report.md` and these tables as CSV.

### Checks before interpreting

```python
# interferon drives M04, M05 and M06 together: if they correlate across donors,
# report them as one signal, not three findings
nas.tl.module_correlation(adata, "ciliated")

# HLA alleles align to the reference with allele-dependent efficiency:
# rerun without HLA genes and check M06/M07 age effects survive
hla_free = adata.copy()
nas.tl.score_modules(hla_free, layer="counts", modules=nas.drop_hla(nas.load_modules()))
nas.tl.donor_summary(hla_free, donor_key="donor_id", celltype_key="population",
                     age_key="age_group", batch_key="chemistry", sex_key="sex",
                     subtype_key="celltype_fine")
nas.tl.age_effects(hla_free, reference="18-39")

# aggregation choices
nas.tl.sensitivity(adata, min_cells_grid=(5, 10, 20), statistics=("auto", "median", "mean"))

# how much of an apparent effect is composition rather than cell state
nas.tl.pooled_scores(adata, ["basal", "ciliated", "secretory"])
```

### Figures

```python
nas.pl.effect_heatmap(adata).savefig("results/cohortA/age_effects.png", dpi=300)
nas.pl.composition(adata).savefig("results/cohortA/composition.png", dpi=300)
combos = [("M05_interferon_stimulated_genes", "ciliated"), ("M03_pigr_antibody_transport", "secretory")]
nas.pl.donor_scores(adata, combos)       # donor values by age group
nas.pl.forest(adata, combos)             # estimates with 95% CIs
```

## 5. Modules

| Module | What it reads out | Tested (primary) in |
|---|---|---|
| M01 epithelial barrier | Tight, adherens and desmosomal junction genes | basal, ciliated, secretory |
| M02 mucociliary axoneme | Motile-cilia structure and FOXJ1 programme | ciliated |
| M02b airway mucins | MUC5AC, MUC5B, SCGB1A1, BPIFA1, BPIFB1 (small panel) | secretory |
| M03 PIGR | Single gene: secretory IgA transport | secretory |
| M03b FCGRT | Single gene: IgG transport (FcRn) | exploratory only |
| M04 innate sensing | Constitutive PRRs and adaptors; sensors that are themselves ISGs removed | epithelium, DC, macrophage, B cell |
| M05 interferon-stimulated genes | ISG response to type I **and** III interferon | all populations |
| M06 MHC-I processing | Immunoproteasome, TAP, HLA class I (interferon-inducible) | all populations |
| M07 MHC-II presentation | HLA class II, CD74, CIITA | DC, macrophage, B cell |
| M08 APC costimulation | Mature/regulatory DC programme | DC, macrophage, B cell |
| M09 plasma secretory machinery | PRDM1, XBP1, IRF4, ER secretory genes (no Ig genes) | plasma |
| M10 IgA fraction | IgA share of Ig heavy-chain counts (logit) | plasma |
| M11 senescence arrest | CDKN1A, CDKN2A and secreted senescence markers (SERPINE1, IGFBP7, GDF15) | basal, ciliated, secretory |
| M12 inflammatory SASP | IL6, IL1A/B, CXCL8, CXCL1, CCL2, TNF (inflammaging core) | epithelium, macrophage |
| M13 resident memory T cells | TRM core programme (CD69, ITGAE, ITGA1, CXCR6, PDCD1, CRTAM, DUSP6) | T cell |

## 6. Function reference

| Function | What it does |
|---|---|
| `nas.pp.check_design` | Donor table; stops if age/batch/sex vary within a donor; warns about confounding |
| `nas.pp.filter_cells` | Scanpy QC metrics, then filtering with per-population thresholds, in place |
| `nas.tl.score_modules` | UCell scores; PIGR and FCGRT as single genes; IgA fraction as a ratio of counts |
| `nas.tl.donor_summary` | One value per donor x population, plus coverage, score distribution and covariates |
| `nas.tl.age_effects` | Age-group tests per module x population, BH within primary and exploratory families |
| `nas.tl.composition` | Population proportions vs age (CLR and arcsine-root) |
| `nas.tl.permutation_null` | False-positive check: age labels shuffled between donors, within batch |
| `nas.tl.module_correlation` | Donor-level correlation between modules within one population |
| `nas.tl.sensitivity` | Re-runs the tests over `min_cells` and summary statistics |
| `nas.tl.pooled_scores` | Scores ignoring cell type, to show how composition can fake an age effect |
| `nas.pl.*` | `effect_heatmap`, `composition`, `donor_scores`, `forest` |
| `nas.write_report` | `report.md` plus all result tables as CSV |
| `nas.load_modules`, `nas.drop_hla` | The module definitions, and a copy without HLA genes |

## Important limits

- The pipeline has not yet been checked against a real cohort with a known age
  effect, so treat first results as exploratory.
- A higher score is not "better for vaccination", and no combined readiness
  score is produced. See [docs/INTERPRETATION.md](docs/INTERPRETATION.md).
- Interferon scores largely reflect recent infection, and ambient RNA can
  inflate the mucin, PIGR and antibody modules.

## Licence

MIT, see [LICENSE](LICENSE). Methods to cite are listed in
[docs/DETAILS.md](docs/DETAILS.md#citing-the-methods-used).
