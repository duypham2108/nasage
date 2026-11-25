"""Module scoring (UCell), ratio/single-gene features and measurability reporting."""
from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp

#: Features anchored on a single defensible gene. These are NOT expanded into
#: artificial signatures and are NOT scored with UCell; they are reported as
#: donor-level pseudobulk expression.
SINGLE_GENE_FEATURES = {"M03_pigr_antibody_transport": "PIGR",
                        "M03b_fcgrt_igg_transport": "FCGRT"}

#: Features that are a ratio of pooled counts, not a UCell score. The donor
#: value is logit(sum(numerator) / sum(denominator)) over the donor's cells of a
#: population. M10: IgA share of Ig heavy-chain transcripts, so that a shift in
#: isotype is not read as a change in secretory capacity (M09).
RATIO_FEATURES = {
    "M10_iga_fraction": dict(
        numerator=["IGHA1", "IGHA2"],
        denominator=["IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHM"]),
}
MAX_RANK = 1500


def feature_type(module_id: str) -> str:
    if module_id in SINGLE_GENE_FEATURES:
        return "single_gene_feature"
    if module_id in RATIO_FEATURES:
        return "ratio_feature"
    return "multi_gene_module"


def module_coverage(adata, modules_mapped: dict, applicability: dict, cell_types: pd.Series,
                    layer=None, detection_floor=0.05) -> pd.DataFrame:
    """Per module x cell type: genes present, detection rate, measurability flag.

    `cell_types` is indexed by obs_names (cells not in it are ignored); `layer`
    holds counts (None = adata.X). A module that is poorly detected is flagged,
    NOT interpreted as biologically absent. Flagged pairs are excluded from testing.
    """
    rows = []
    L = adata.X if layer is None else adata.layers[layer]
    gpos = {g: i for i, g in enumerate(adata.var_names)}
    ct_all = pd.Series(cell_types).reindex(adata.obs_names).astype(object)
    for mod, genes in modules_mapped.items():
        idx = [gpos[g] for g in genes if g in gpos]
        for ct in pd.unique(ct_all.dropna()):
            m = (ct_all == ct).values
            if m.sum() == 0 or not idx:
                continue
            sub = L[m][:, idx]
            pos = sub > 0
            det = np.asarray(pos.mean(0)).ravel()
            per_cell = np.asarray(pos.sum(1)).ravel()
            rows.append(dict(module_id=mod, cell_type=ct, n_genes_in_data=len(idx),
                             n_cells=int(m.sum()), mean_detection_rate=float(det.mean()),
                             frac_genes_detected=float((det > detection_floor).mean()),
                             frac_cells_any_detected=float((per_cell > 0).mean()),
                             median_genes_detected_per_cell=float(np.median(per_cell)),
                             applicability=float(applicability.get(mod, {}).get(ct, 0.0))))
    cov = pd.DataFrame(rows)
    cov["feature_type"] = cov.module_id.map(feature_type)
    cov["poorly_measurable"] = np.select(
        [cov.feature_type == "single_gene_feature", cov.feature_type == "ratio_feature"],
        [cov.mean_detection_rate < 0.05,
         # a ratio needs Ig transcript in most cells of the population
         cov.frac_cells_any_detected < 0.30],
        (cov.mean_detection_rate < 0.10) | (cov.n_genes_in_data < 3)).astype(bool)
    return cov


def _n_detected(M) -> np.ndarray:
    """Per-cell count of detected genes, without materialising a boolean copy.

    ``(M > 0).sum(1)`` builds a second sparse matrix and reduces it in int64, which
    costs ~12 bytes per stored element -- several GB on a large cohort. For CSR with
    no explicit zeros the row lengths already give the answer; otherwise fall back to
    counting the nonzero values per row.
    """
    if sp.issparse(M):
        M = M.tocsr() if not sp.isspmatrix_csr(M) else M
        lens = np.diff(M.indptr)
        if M.nnz and not (M.data == 0).any():
            return lens.astype(np.int64)
        rows = np.repeat(np.arange(M.shape[0]), lens)
        return np.bincount(rows[M.data != 0], minlength=M.shape[0]).astype(np.int64)
    return np.asarray((np.asarray(M) > 0).sum(1)).ravel()


def score_modules(adata, modules_mapped: dict, layer=None, max_rank=MAX_RANK,
                  n_jobs=1) -> list[str]:
    """UCell scores for multi-gene modules, plus single-gene and ratio features.

    `layer` holds raw counts (None = adata.X). UCell ranks genes within each
    cell, so raw counts and normalised values give the same scores. Uses the
    published pyUCell implementation (Andreatta & Carmona 2021). Scores are
    written into adata.obs and the column names returned. Each cell's detected
    gene count is stored as ``obs["nis_n_genes"]`` (the depth covariate).
    """
    from pyucell import compute_ucell_scores
    present = set(adata.var_names)
    sigs = {m: [g for g in gl if g in present] for m, gl in modules_mapped.items()
            if feature_type(m) == "multi_gene_module"}
    sigs = {m: gl for m, gl in sigs.items() if len(gl) >= 3}
    compute_ucell_scores(adata, sigs, layer=layer, max_rank=max_rank,
                         suffix="_UCell", n_jobs=n_jobs)
    cols = [f"{m}_UCell" for m in sigs]
    gpos = {g: i for i, g in enumerate(adata.var_names)}
    counts = adata.X if layer is None else adata.layers[layer]
    adata.obs["nis_n_genes"] = _n_detected(counts)
    lib = np.maximum(np.asarray(counts.sum(1)).ravel(), 1)
    for mod, gene in SINGLE_GENE_FEATURES.items():
        if gene not in gpos:
            continue
        x = _column(counts, gpos[gene])
        adata.obs[f"{mod}_cpm"] = x / lib * 1e4                     # CP10K
        adata.obs[f"{mod}_logexpr"] = np.log1p(adata.obs[f"{mod}_cpm"].to_numpy())
        cols.append(f"{mod}_logexpr")
    for mod, spec in RATIO_FEATURES.items():
        num = [gpos[g] for g in spec["numerator"] if g in gpos]
        den = [gpos[g] for g in spec["denominator"] if g in gpos]
        if not num or len(den) <= len(num):        # need at least one non-numerator gene
            continue
        adata.obs[f"{mod}_num"] = np.asarray(counts[:, num].sum(1)).ravel().astype(float)
        adata.obs[f"{mod}_den"] = np.asarray(counts[:, den].sum(1)).ravel().astype(float)
        cols.append(f"{mod}_logit")
    return cols


def score_distribution(obs: pd.DataFrame, score_cols, ceiling=0.8) -> pd.DataFrame:
    """Per UCell module x cell type: is the score near its ceiling or mostly zero?

    UCell scores lie in [0, 1]. A module made of the most highly expressed genes
    of a cell (mucins in secretory cells, Ig machinery in plasma cells,
    HLA-A/B/C and B2M everywhere) sits near the top and has little room to
    show an increase; a module of rarely detected genes (TLRs, cytokines) is
    zero in most cells, and a median then discards the signal. Both are
    flagged here. The thresholds are heuristics for inspection, not tests.
    """
    rows = []
    for col in [c for c in score_cols if c.endswith("_UCell")]:
        for ct, v in obs.groupby("cell_type", observed=True)[col]:
            v = v.to_numpy()
            rows.append(dict(module_id=col[:-len("_UCell")], cell_type=ct, n_cells=len(v),
                             median=float(np.median(v)), q90=float(np.quantile(v, 0.9)),
                             frac_zero=float((v == 0).mean()),
                             frac_above_ceiling=float((v > ceiling).mean())))
    df = pd.DataFrame(rows)
    df["near_ceiling"] = df["median"] > ceiling
    df["mostly_zero"] = df.frac_zero > 0.5
    return df


def _column(M, j) -> np.ndarray:
    """Column j as a flat float array, for sparse or dense matrices."""
    col = M[:, j]
    if sp.issparse(col):
        col = col.toarray()
    return np.asarray(col, dtype=float).ravel()


def score_column(module_id: str) -> str:
    return {"single_gene_feature": f"{module_id}_logexpr",
            "ratio_feature": f"{module_id}_logit"}.get(feature_type(module_id),
                                                       f"{module_id}_UCell")
