"""Donor-level inference: age-bin models, composition, permutation null."""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

from .scoring import score_column

MIN_DONORS_PER_BIN = 3


def build_test_family(applicability: dict, coverage: pd.DataFrame, score_cols,
                      tier: float = 1.0, modules_present=None) -> list[tuple]:
    """The family of (module, cell type) pairs to be tested, fixed in advance.

    Only pairs marked `tier` in the applicability matrix and not flagged as
    poorly measurable enter. Exploratory pairs (tier 0.5) are analysed but are
    adjusted in their own family and never promoted to primary.
    """
    bad = set(map(tuple, coverage.loc[coverage.poorly_measurable,
                                      ["module_id", "cell_type"]].values))
    name = "primary" if tier == 1.0 else "exploratory"
    pairs = []
    for mod, d in applicability.items():
        if modules_present is not None and mod not in modules_present:
            continue
        if score_column(mod) not in score_cols:
            continue
        for ct, v in d.items():
            if float(v) == tier and (mod, ct) not in bad:
                pairs.append((mod, ct, name))
    return pairs


def _formula(score_col, covariates, ref_bin):
    rhs = [f'C(age_bin, Treatment(reference="{ref_bin}"))'] + list(covariates)
    return f'Q("{score_col}") ~ ' + " + ".join(rhs)


MIN_RESID_DF = 3
#: Optional numeric covariates from `aggregate.donor_celltype_summary`, dropped
#: (last first) when the donor count cannot support them.
DEPTH_COVARIATES = ("log_genes",)
SUBTYPE_COVARIATES = ("subfrac_1", "subfrac_2")


def _usable(d: pd.DataFrame, c: str) -> bool:
    return c in d and d[c].notna().all() and d[c].astype(str).nunique() > 1


def _n_params(d, levels, covs) -> int:
    p = len(levels)                                   # intercept + (levels - 1)
    for c in covs:
        p += 1 if pd.api.types.is_numeric_dtype(d[c]) else d[c].astype(str).nunique() - 1
    return p


def fit_one(df, score_col, ref_bin, covariates=("batch", "sex"),
            min_donors_per_bin=MIN_DONORS_PER_BIN, adjust_depth=True, adjust_subtypes=True,
            min_resid_df=MIN_RESID_DF):
    """OLS of one donor-level score on age group + covariates, within one population.

    Besides `covariates` (categorical: batch, sex), the model adds by default:
      * ``log_genes`` (adjust_depth) - UCell scores rise with detected genes per
        cell, and sequencing depth/sample quality can differ between age groups;
      * ``subfrac_1/2`` (adjust_subtypes) - subtype fractions within the
        population, so a shift in subtype mix is not read as cell state.
    Optional covariates are dropped (subtype first, then depth) if they leave
    fewer than `min_resid_df` residual degrees of freedom or make the design
    rank-deficient; what was used and dropped is reported.
    """
    d = df.dropna(subset=[score_col]).copy()
    counts = d.age_bin.astype(str).value_counts()
    levels = [l for l in pd.Categorical(df.age_bin).categories if l in set(d.age_bin.astype(str))]
    if len(counts) < 2 or counts.min() < min_donors_per_bin or ref_bin not in counts:
        return None, dict(reason="insufficient donors per age bin", n_donors=int(len(d)),
                          donors_per_bin=counts.to_dict())
    d["age_bin"] = pd.Categorical(d.age_bin.astype(str), categories=levels)
    base = [c for c in covariates if _usable(d, c)]
    for c in base:
        if not pd.api.types.is_numeric_dtype(d[c]):
            d[c] = d[c].astype(str)
    optional = ([c for c in DEPTH_COVARIATES if adjust_depth] +
                [c for c in SUBTYPE_COVARIATES if adjust_subtypes])
    optional = [c for c in optional if _usable(d, c)]
    dropped = []
    while True:
        covs = base + optional
        if len(d) - _n_params(d, levels, covs) >= min_resid_df or not optional:
            try:
                m = smf.ols(_formula(score_col, covs, ref_bin), data=d).fit()
            except Exception as exc:                                 # pragma: no cover
                return None, dict(reason=f"fit failed: {exc}", n_donors=int(len(d)))
            if _full_rank(m) or not optional:
                break
        dropped.append(optional.pop())
    if not _full_rank(m):
        # statsmodels fits a rank-deficient design via pinv without error, which
        # would report arbitrary age estimates and p-values.
        return None, dict(reason="design not full rank (age bin aliased with "
                                 "batch/sex in this cell type)",
                          n_donors=int(len(d)), donors_per_bin=counts.to_dict())
    if m.df_resid < 1:
        return None, dict(reason="no residual degrees of freedom",
                          n_donors=int(len(d)), donors_per_bin=counts.to_dict())
    return m, dict(n_donors=int(len(d)), donors_per_bin=counts.to_dict(),
                   covariates=",".join(covs), dropped_covariates=",".join(dropped))


def _full_rank(m) -> bool:
    return int(np.linalg.matrix_rank(m.model.exog)) == m.model.exog.shape[1]


def _age_wald(m):
    terms = [t for t in m.model.exog_names if t.startswith("C(age_bin")]
    R = np.zeros((len(terms), len(m.model.exog_names)))
    for i, t in enumerate(terms):
        R[i, m.model.exog_names.index(t)] = 1
    return terms, m.f_test(R)


def fit_age_models(summary: pd.DataFrame, pairs, ref_bin: str,
                   covariates=("batch", "sex"), adjust_depth=True,
                   adjust_subtypes=True) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Omnibus age test and contrasts for every (module, cell type) pair.

    Donor is the unit of replication; cells are not. Effect sizes come with
    95% confidence intervals, the donor count and the donors per age bin.
    See `fit_one` for the depth and subtype covariates.
    """
    groups = sorted(map(str, pd.unique(summary.age_bin.dropna())))
    if str(ref_bin) not in groups:
        raise ValueError(f"ref_bin '{ref_bin}' is not one of the age groups in the data: {groups}")
    # Speed: pandas deep-copies .attrs on every slice, and the summary can hold
    # thousands of score columns (random-gene-set control). Drop attrs once, split
    # rows by cell type once, and hand each fit only the columns fit_one reads.
    summary = summary.copy(deep=False)
    summary.attrs = {}
    meta = [c for c in ("age_bin", *covariates, *DEPTH_COVARIATES, *SUBTYPE_COVARIATES)
            if c in summary.columns]
    rows = summary.groupby("cell_type", observed=True).indices      # positions, not labels
    omni, contr = [], []
    for mod, ct, tier in pairs:
        col = score_column(mod)
        cols = meta + ([col] if col in summary.columns else [])
        sub = (summary.iloc[rows[ct], summary.columns.get_indexer(cols)]
               if ct in rows else summary.iloc[:0])
        m, info = fit_one(sub, col, ref_bin, covariates,
                          adjust_depth=adjust_depth, adjust_subtypes=adjust_subtypes)
        base = dict(module_id=mod, cell_type=ct, tier=tier, score_column=col,
                    n_donors=info.get("n_donors"), reason=info.get("reason"),
                    donors_per_bin=str(info.get("donors_per_bin", {})),
                    covariates=info.get("covariates", ""),
                    dropped_covariates=info.get("dropped_covariates", ""))
        if m is None:
            omni.append(dict(**base, estimable=False, f_stat=np.nan, p_omnibus=np.nan,
                             residual_sd=np.nan))
            continue
        terms, ft = _age_wald(m)
        sd = float(np.sqrt(m.mse_resid))
        omni.append(dict(**base, estimable=True, f_stat=float(np.squeeze(ft.fvalue)),
                         p_omnibus=float(np.squeeze(ft.pvalue)), residual_sd=sd))
        ci = m.conf_int(alpha=0.05)
        for t in terms:
            lvl = t.split("T.")[-1].rstrip("]")
            contr.append(dict(module_id=mod, cell_type=ct, tier=tier, score_column=col,
                              contrast=f"{lvl} vs {ref_bin}", estimate=float(m.params[t]),
                              ci_low=float(ci.loc[t, 0]), ci_high=float(ci.loc[t, 1]),
                              se=float(m.bse[t]), p_value=float(m.pvalues[t]),
                              std_effect=float(m.params[t] / sd) if sd > 0 else np.nan,
                              n_donors=info["n_donors"]))
    return pd.DataFrame(omni), pd.DataFrame(contr)


def adjust_family(df, pcol="p_omnibus", mask=None, newcol="p_adj_BH"):
    """Benjamini-Hochberg inside one explicitly declared family of tests."""
    df = df.copy()
    if mask is None:
        mask = pd.Series(True, index=df.index)
    sel = mask & df[pcol].notna()
    df[newcol] = np.nan
    if sel.sum():
        df.loc[sel, newcol] = multipletests(df.loc[sel, pcol], method="fdr_bh")[1]
    df["family_size"] = int(sel.sum())
    return df


def composition_analysis(obs: pd.DataFrame, ref_bin: str, covariates=("batch", "sex"),
                         pseudocount=0.5) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Donor-level compositional regression (CLR, plus arcsine-root as a check).

    Proportions are RELATIVE. A rise in one population forces a fall in the
    others, and none of this measures absolute abundance in the tissue: a
    nasal brush samples what detaches, not what is there.
    """
    tab = obs.groupby(["donor", "cell_type"], observed=False).size().unstack(fill_value=0)
    meta = obs[["donor", "age_bin", "sex", "batch"]].drop_duplicates("donor").set_index("donor").loc[tab.index]
    prop = tab.div(tab.sum(1), axis=0)
    lg = np.log(tab + pseudocount)
    clr = lg.sub(lg.mean(1), axis=0)
    asr = np.arcsin(np.sqrt(prop.clip(0, 1)))
    present = set(meta.age_bin.dropna().astype(str))      # empty age groups would alias
    levels = [l for l in pd.Categorical(meta.age_bin).categories if str(l) in present]
    out = []
    for name, M in (("CLR", clr), ("asin_sqrt", asr)):
        for ct in tab.columns:
            d = meta.copy(); d["y"] = M[ct].values
            d["age_bin"] = pd.Categorical(d.age_bin.astype(str), categories=levels)
            covs = [c for c in covariates if c in d and d[c].astype(str).nunique() > 1]
            for c in covs:
                d[c] = d[c].astype(str)
            m = smf.ols(_formula("y", covs, ref_bin).replace('Q("y")', "y"), data=d).fit()
            rec = dict(transform=name, cell_type=ct,
                       mean_proportion=float(prop[ct].mean()),
                       n_donors_with_zero=int((tab[ct] == 0).sum()))
            if not _full_rank(m):
                out.append(dict(rec, f_stat=np.nan, p_omnibus=np.nan,
                                reason="design not full rank"))
                continue
            terms, ft = _age_wald(m)
            ci = m.conf_int()
            rec.update(f_stat=float(np.squeeze(ft.fvalue)),
                       p_omnibus=float(np.squeeze(ft.pvalue)))
            for t in terms:
                lvl = t.split("T.")[-1].rstrip("]")
                rec[f"beta_{lvl}"] = float(m.params[t])
                rec[f"ci_low_{lvl}"] = float(ci.loc[t, 0])
                rec[f"ci_high_{lvl}"] = float(ci.loc[t, 1])
            out.append(rec)
    res = pd.DataFrame(out)
    res["p_adj_BH"] = np.nan
    for tr in res["transform"].unique():
        s = (res["transform"] == tr) & res["p_omnibus"].notna()
        if s.any():
            res.loc[s, "p_adj_BH"] = multipletests(res.loc[s, "p_omnibus"], method="fdr_bh")[1]
    return res, prop


def permutation_null(summary: pd.DataFrame, pairs, ref_bin: str, n_perm: int = 200,
                     covariates=("batch", "sex"), seed: int = 0,
                     stratify: str | None = "batch", adjust_depth=True,
                     adjust_subtypes=True) -> pd.DataFrame:
    """Calibrate the family by permuting age labels ACROSS DONORS.

    With real data you cannot re-simulate, so this is how you check that the
    chosen family and model are not generating significant results on their own.
    The permutation keeps each donor's cells together and permutes the donor's
    age label, so the donor-level correlation structure is preserved.

    By default labels are permuted WITHIN batch, so the observed batch-age
    association is kept; an unrestricted permutation breaks it and understates
    the false-positive rate when batch and age are partly confounded. Pass
    ``stratify=None`` for an unrestricted permutation.
    """
    rng = np.random.default_rng(seed)
    keep = ["donor", "age_bin"] + ([stratify] if stratify and stratify in summary else [])
    donor_meta = summary[keep].drop_duplicates("donor").reset_index(drop=True)
    strata = (donor_meta.groupby(stratify, observed=True).indices.values()
              if len(keep) == 3 else [np.arange(len(donor_meta))])
    strata = list(strata)
    rows = []
    for k in range(n_perm):
        perm = donor_meta[["donor", "age_bin"]].copy()
        lab = perm.age_bin.to_numpy().copy()
        for idx in strata:
            lab[idx] = rng.permutation(lab[idx])
        perm["age_bin"] = lab
        s = summary.drop(columns="age_bin").merge(perm, on="donor")
        om, _ = fit_age_models(s, pairs, ref_bin, covariates, adjust_depth=adjust_depth,
                               adjust_subtypes=adjust_subtypes)
        om = adjust_family(om, "p_omnibus", om.estimable)
        rows.append(dict(perm=k, n_tests=int(om.estimable.sum()),
                         n_p_lt_05=int((om.p_omnibus < 0.05).sum()),
                         n_bh_lt_05=int((om.p_adj_BH < 0.05).sum()),
                         min_p=float(om.p_omnibus.min()) if om.estimable.any() else np.nan))
    df = pd.DataFrame(rows)
    df.attrs["per_test_false_positive_rate"] = float(df.n_p_lt_05.sum() / max(df.n_tests.sum(), 1))
    df.attrs["family_wise_rate_BH05"] = float((df.n_bh_lt_05 > 0).mean())
    return df


def sensitivity_sweep(obs: pd.DataFrame, coverage: pd.DataFrame, applicability, score_cols, ref_bin,
                      min_cells_grid=(5, 10, 20, 30),
                      statistics=("auto", "median", "mean")) -> pd.DataFrame:
    """Re-run the donor-level analysis over aggregation choices."""
    from .aggregate import donor_celltype_summary
    pairs = build_test_family(applicability, coverage, score_cols)
    rows = []
    for mc in min_cells_grid:
        for st in statistics:
            s = donor_celltype_summary(obs, score_cols, min_cells=mc, statistic=st)
            om, _ = fit_age_models(s, pairs, ref_bin)
            om = adjust_family(om, "p_omnibus", om.estimable)
            om["min_cells"] = mc; om["statistic"] = st
            rows.append(om)
    return pd.concat(rows, ignore_index=True)


def module_correlation(summary: pd.DataFrame, cell_type: str, score_cols=None,
                       method: str = "spearman") -> pd.DataFrame:
    """Donor-level correlation between modules within one population.

    Interferon drives ISGs (M05), MHC-I processing (M06) and part of innate
    sensing (M04) together. If two modules correlate strongly across donors,
    report them as one signal rather than as independent findings.
    """
    d = summary[(summary.cell_type == cell_type) & summary.eligible]
    cols = [c for c in (score_cols or [c for c in d.columns
                                       if c.endswith(("_UCell", "_logexpr", "_logit"))])
            if c in d and d[c].notna().sum() >= 3]
    return d[cols].corr(method=method)
