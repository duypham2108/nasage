"""
Synthetic nasal-brush scRNA-seq simulator.

ALL DATA PRODUCED BY THIS MODULE ARE SYNTHETIC. No real donor, sample or
individual-level dataset is used anywhere in this project.

Generative model (per cell i, of cell type c, from donor d in batch b, age bin a):

    log mu[g] = base[g]
              + celltype_offset[g, c]
              + donor_gene_effect[g, d]
              + donor_module_tone[module(g), d]
              + batch_effect[g, b]
              + age_effect[g, c, a]              (scenario "planted" only)

    rel[g]    = exp(log mu[g]) / sum_g exp(log mu[g])     (per cell type/donor/bin)
    m[g]      = L_i * rel[g]                              (L_i = cell library size)
    counts[g] ~ NegativeBinomial(mean = m[g], var = m[g] + phi[g] * m[g]^2)

Cell-type composition per donor is Dirichlet-distributed around age-bin-specific
mean proportions; rare immune populations are dropped entirely in a subset of
donors to mimic a brush sample that simply failed to recover them.

The ground-truth parameter object returned alongside the data is written to a
separate directory and is NEVER read by the analysis pipeline.
"""
from __future__ import annotations

import json
import numpy as np
import pandas as pd
import scipy.sparse as sp

# ----------------------------------------------------------------------------- config
AGE_BINS = ["5-11", "12-17", "18-39", "40-65"]
CELLTYPES = ["basal", "ciliated", "secretory", "dendritic", "macrophage",
             "Tcell", "Bcell", "plasma"]

# nasal brush: epithelium dominant, immune populations sparse
BASE_PROPS = {"basal": 0.20, "ciliated": 0.33, "secretory": 0.32,
              "dendritic": 0.028, "macrophage": 0.035, "Tcell": 0.052,
              "Bcell": 0.020, "plasma": 0.015}   # epithelium 85%, immune 15%

MEDIAN_LIBSIZE = {"basal": 6000, "ciliated": 7000, "secretory": 6500,
                  "dendritic": 4000, "macrophage": 5500, "Tcell": 2500,
                  "Bcell": 2500, "plasma": 9000}

# curated marker genes used only to make simulated clusters identifiable
MARKERS = {
    "basal": ["KRT5", "KRT14", "TP63", "DLK2", "MIR205HG"],
    "ciliated": ["FOXJ1", "TPPP3", "DNAH5", "CAPS", "C20orf85"],
    "secretory": ["SCGB1A1", "MUC5AC", "MUC5B", "BPIFA1", "SERPINB3"],
    "dendritic": ["ITGAX", "CD1C", "CLEC9A", "FCER1A", "LAMP3"],
    "macrophage": ["CD68", "LYZ", "MRC1", "MSR1", "C1QA"],
    "Tcell": ["CD3D", "CD3E", "CD2", "IL7R", "TRAC"],
    "Bcell": ["CD79A", "MS4A1", "CD19", "BANK1", "TNFRSF13C"],
    "plasma": ["MZB1", "JCHAIN", "DERL3", "TNFRSF17", "SSR4"],
}

# --------------------------------------------------------------------- planted effects
# Prespecified before any analysis. Values are natural-log fold changes applied
# to the per-gene expression rate, by age bin, within the named cell types.
PLANTED = {
    "P1_pigr_secretory": dict(
        module="M03_pigr_antibody_transport", cell_types=["secretory"],
        logfc_by_bin=[0.00, 0.20, 0.38, 0.45], gene_heterogeneity=False,
        note="single-gene feature (PIGR) increasing with age bin"),
    "P2_ifn_epithelial": dict(
        module="M05_interferon_stimulated_genes", cell_types=["ciliated", "secretory"],
        logfc_by_bin=[0.45, 0.28, 0.00, -0.08], gene_heterogeneity=True,
        frac_nonresponding=0.25, het_sd=0.45,
        note="context-dependent ISG tone, highest in the youngest bin; NOT to be "
             "labelled beneficial"),
    "P4_costim_dc": dict(
        module="M08_apc_costimulation", cell_types=["dendritic"],
        logfc_by_bin=[0.00, 0.12, 0.22, 0.26], gene_heterogeneity=True,
        frac_nonresponding=0.30, het_sd=0.40,
        note="costimulation increasing with age bin, so that no age bin is "
             "universally superior"),
}

# P3: composition shift with NO change in cell state
COMPOSITION_MULTIPLIER = {
    "ciliated":  [0.80, 0.92, 1.06, 1.16],
    "secretory": [1.18, 1.07, 0.96, 0.88],
    "Tcell":     [0.70, 0.85, 1.10, 1.30],
}


# ----------------------------------------------------------------------------- helpers
def _build_gene_universe(modules: dict, n_genes: int, rng) -> list:
    module_genes = sorted({g for v in modules.values() for g in v})
    marker_genes = sorted({g for v in MARKERS.values() for g in v})
    named = list(dict.fromkeys(module_genes + marker_genes))
    n_bg = n_genes - len(named)
    if n_bg < 0:
        raise ValueError("n_genes too small for the module+marker gene set")
    bg = [f"BG{str(i).zfill(5)}" for i in range(n_bg)]
    genes = named + bg
    rng.shuffle(genes)
    return genes


def _celltype_offsets(genes, modules, applicability, rng):
    """log-scale cell-type specific offsets, genes x celltypes."""
    G, C = len(genes), len(CELLTYPES)
    gidx = {g: i for i, g in enumerate(genes)}
    off = np.zeros((G, C), dtype=np.float32)

    # generic cell-type programmes on background/unassigned genes
    for ci, ct in enumerate(CELLTYPES):
        up = rng.random(G) < 0.06
        dn = rng.random(G) < 0.06
        off[up, ci] += rng.normal(1.8, 0.6, up.sum()).astype(np.float32)
        off[dn, ci] -= rng.normal(1.5, 0.5, dn.sum()).astype(np.float32)

    # module genes follow the applicability matrix (defined before any age effect)
    tier_offset = {1.0: 1.30, 0.5: 0.45, 0.0: -1.10}
    for mod, genelist in modules.items():
        for ci, ct in enumerate(CELLTYPES):
            tier = applicability[mod].get(ct, 0.0)
            base = tier_offset[float(tier)]
            for g in genelist:
                if g in gidx:
                    off[gidx[g], ci] = base + rng.normal(0, 0.25)

    # curated markers: strong, specific
    for ct, mk in MARKERS.items():
        ci = CELLTYPES.index(ct)
        for g in mk:
            if g in gidx:
                off[gidx[g], :] -= 1.2
                off[gidx[g], ci] = 2.6 + rng.normal(0, 0.3)
    return off


def _donor_table(n_per_bin, rng):
    rows = []
    for bi, ab in enumerate(AGE_BINS):
        sexes = (["F", "M"] * ((n_per_bin + 1) // 2))[:n_per_bin]
        batches = [1, 2] * (n_per_bin // 2)
        rng.shuffle(sexes)
        rng.shuffle(batches)
        for k in range(n_per_bin):
            rows.append(dict(donor=f"D{bi+1}{str(k+1).zfill(2)}", age_bin=ab,
                             age_bin_index=bi, sex=sexes[k], batch=f"B{batches[k]}"))
    return pd.DataFrame(rows)


def _donor_proportions(donors, scenario, rng, conc=55.0, p_drop=0.28):
    """Dirichlet donor proportions, optional age-driven composition shift,
    plus complete dropout of a rare immune population in some donors."""
    rare = ["dendritic", "Bcell", "plasma", "macrophage"]
    props, dropped = {}, {}
    for _, r in donors.iterrows():
        mean = np.array([BASE_PROPS[c] for c in CELLTYPES], dtype=float)
        if scenario == "planted":
            for ct, mult in COMPOSITION_MULTIPLIER.items():
                mean[CELLTYPES.index(ct)] *= mult[r.age_bin_index]
        mean = mean / mean.sum()
        p = rng.dirichlet(mean * conc)
        drop = []
        if rng.random() < p_drop:
            ct = rare[rng.integers(len(rare))]
            p[CELLTYPES.index(ct)] = 0.0
            drop.append(ct)
            p = p / p.sum()
        props[r.donor] = p
        dropped[r.donor] = drop
    return props, dropped


# ------------------------------------------------------------------------------ main
def simulate(scenario: str,
             modules: dict,
             applicability: dict,
             seed: int = 20261017,
             n_per_bin: int = 8,
             cells_per_donor: int = 500,
             n_genes: int = 4000,
             libsize_scale: float = 1.0,
             capture_efficiency: float = 1.0):
    """Simulate one synthetic nasal-brush scRNA-seq dataset.

    scenario : "null" (no true age effect) or "planted" (prespecified effects)
    Returns (AnnData, truth_dict).
    """
    import anndata as ad

    if scenario not in ("null", "planted"):
        raise ValueError("scenario must be 'null' or 'planted'")
    rng = np.random.default_rng(seed)

    genes = _build_gene_universe(modules, n_genes, rng)
    gidx = {g: i for i, g in enumerate(genes)}
    G = len(genes)

    base = rng.normal(0.0, 1.7, G).astype(np.float32)          # log relative abundance
    phi = np.exp(rng.normal(np.log(0.35), 0.30, G)).astype(np.float32)   # NB dispersion
    ct_off = _celltype_offsets(genes, modules, applicability, rng)

    donors = _donor_table(n_per_bin, rng)
    batch_eff = {}
    for b in ["B1", "B2"]:
        e = np.zeros(G, dtype=np.float32)
        hit = rng.random(G) < 0.5
        e[hit] = rng.normal(0, 0.28, hit.sum())
        batch_eff[b] = e
    batch_lib_mult = {"B1": 1.0, "B2": 0.78}

    donor_gene_eff = {d: rng.normal(0, 0.25, G).astype(np.float32) for d in donors.donor}
    donor_module_tone = {d: {m: rng.normal(0, 0.15) for m in modules} for d in donors.donor}

    # ---- age effects (planted scenario only), with per-gene heterogeneity
    age_eff = np.zeros((G, len(CELLTYPES), len(AGE_BINS)), dtype=np.float32)
    truth_genes = []
    if scenario == "planted":
        for pid, spec in PLANTED.items():
            genelist = [g for g in modules[spec["module"]] if g in gidx]
            if spec["gene_heterogeneity"]:
                w = np.clip(rng.normal(1.0, spec["het_sd"], len(genelist)), 0.0, 2.5)
                nz = rng.random(len(genelist)) < spec["frac_nonresponding"]
                w[nz] = 0.0
            else:
                w = np.ones(len(genelist))
            for ct in spec["cell_types"]:
                ci = CELLTYPES.index(ct)
                for gi, g in enumerate(genelist):
                    for ai in range(len(AGE_BINS)):
                        fc = spec["logfc_by_bin"][ai] * w[gi]
                        age_eff[gidx[g], ci, ai] = fc
                        truth_genes.append(dict(planted_id=pid, module=spec["module"],
                                                gene=g, cell_type=ct,
                                                age_bin=AGE_BINS[ai],
                                                true_logfc=float(fc),
                                                gene_weight=float(w[gi])))

    props, dropped = _donor_proportions(donors, scenario, rng)

    # ---- generate counts donor by donor
    X_blocks, obs_rows = [], []
    for _, r in donors.iterrows():
        n_cells = int(np.clip(rng.lognormal(np.log(cells_per_donor), 0.28), 150, 1200))
        p = props[r.donor]
        counts_by_ct = rng.multinomial(n_cells, p)
        for ci, ct in enumerate(CELLTYPES):
            n_c = int(counts_by_ct[ci])
            if n_c == 0:
                continue
            tone = np.zeros(G, dtype=np.float32)
            for m, gl in modules.items():
                t = donor_module_tone[r.donor][m]
                for g in gl:
                    if g in gidx:
                        tone[gidx[g]] += t
            logmu = (base + ct_off[:, ci] + donor_gene_eff[r.donor]
                     + tone + batch_eff[r.batch] + age_eff[:, ci, r.age_bin_index])
            rel = np.exp(logmu.astype(np.float64))
            rel /= rel.sum()
            L = rng.lognormal(np.log(MEDIAN_LIBSIZE[ct] * batch_lib_mult[r.batch]
                                     * libsize_scale * capture_efficiency), 0.40, n_c)
            m_mat = np.outer(L, rel)                           # cells x genes
            nb_n = (1.0 / phi)[None, :]
            pvec = nb_n / (nb_n + m_mat)
            blk = rng.negative_binomial(np.broadcast_to(nb_n, m_mat.shape), pvec)
            X_blocks.append(sp.csr_matrix(blk.astype(np.int32)))
            for k in range(n_c):
                obs_rows.append((f"{r.donor}_{ct}_{k}", r.donor, r.age_bin,
                                 r.age_bin_index, r.sex, r.batch, ct))

    X = sp.vstack(X_blocks, format="csr")
    obs = pd.DataFrame(obs_rows, columns=["cell_id", "donor", "age_bin",
                                          "age_bin_index", "sex", "batch", "cell_type"]
                       ).set_index("cell_id")
    for c in ["donor", "age_bin", "sex", "batch", "cell_type"]:
        obs[c] = obs[c].astype("category")
    adata = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=pd.Index(genes, name="gene")))
    # X holds the RAW COUNTS; no duplicate layer is stored (file size).
    adata.uns["X_is_raw_counts"] = True
    adata.uns["synthetic"] = True
    adata.uns["scenario"] = scenario
    adata.uns["seed"] = seed
    adata.uns["disclaimer"] = ("SYNTHETIC DATA generated by nasal_sim.simulate(); "
                               "no real individual-level data were used.")

    truth = dict(
        scenario=scenario, seed=seed, n_genes=G, n_cells=int(adata.n_obs),
        age_bins=AGE_BINS, cell_types=CELLTYPES,
        planted_specs=PLANTED if scenario == "planted" else {},
        composition_multiplier=COMPOSITION_MULTIPLIER if scenario == "planted" else {},
        donor_proportions={d: dict(zip(CELLTYPES, map(float, v))) for d, v in props.items()},
        dropped_populations=dropped,
        true_gene_effects=truth_genes,
        base_props=BASE_PROPS, median_libsize=MEDIAN_LIBSIZE,
        batch_libsize_multiplier=batch_lib_mult,
        note=("Ground truth. The analysis pipeline must not read this object: it is "
              "used only to evaluate recovery after the fact."),
    )
    return adata, truth


def module_true_effect_table(truth: dict) -> pd.DataFrame:
    """Module-level mean planted log fold change per cell type and age bin."""
    if not truth["true_gene_effects"]:
        return pd.DataFrame(columns=["module", "cell_type", "age_bin",
                                     "mean_true_logfc", "n_genes_affected"])
    df = pd.DataFrame(truth["true_gene_effects"])
    out = (df.groupby(["module", "cell_type", "age_bin"], as_index=False)
             .agg(mean_true_logfc=("true_logfc", "mean"),
                  n_genes_affected=("true_logfc", lambda s: int((s != 0).sum()))))
    return out


if __name__ == "__main__":
    from .modules import load_applicability, load_modules
    modules = {k: v["genes"] for k, v in load_modules().items()}
    appl = load_applicability()
    a, t = simulate("planted", modules, appl)
    print(a)
    print(module_true_effect_table(t).head())
