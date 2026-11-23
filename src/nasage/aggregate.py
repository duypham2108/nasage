"""Donor x cell-type aggregation with explicit missingness."""
from __future__ import annotations

import numpy as np
import pandas as pd

MIN_CELLS_DEFAULT = 10
MAX_SUBTYPE_COVARIATES = 2
_STATS = {"median": np.median, "mean": np.mean,
          # winsorised at the 10th/90th percentiles
          "winsorized_mean": lambda v: float(pd.Series(v).clip(*np.percentile(v, [10, 90])).mean())}
_STATS["trimmed_mean"] = _STATS["winsorized_mean"]   # backwards-compatible alias


def _cell_value(g: pd.DataFrame, c: str, fn) -> float:
    """Donor x cell-type value of one score column."""
    if c.endswith("_logexpr"):
        # pseudobulk: a per-cell median of a sparsely detected gene is usually zero
        return float(np.log1p(np.mean(g[c[:-len("_logexpr")] + "_cpm"].values)))
    if c.endswith("_logit"):
        base = c[:-len("_logit")]
        num, den = float(g[base + "_num"].sum()), float(g[base + "_den"].sum())
        if den == 0:
            return np.nan
        p = (num + 0.5) / (den + 1.0)
        return float(np.log(p / (1 - p)))
    return float(fn(g[c].values))


def choose_statistic(obs: pd.DataFrame, score_cols, zero_frac=0.5) -> pd.DataFrame:
    """Median or mean per (UCell module, cell type), decided without looking at age.

    Where more than `zero_frac` of a population's cells score exactly zero (rarely
    detected genes: TLRs, cytokines), the median is zero for most donors and
    throws the signal away, so the mean is used; otherwise the median.
    """
    rows = []
    for c in [c for c in score_cols if c.endswith("_UCell")]:
        for ct, v in obs.groupby("cell_type", observed=True)[c]:
            fz = float((v.to_numpy() == 0).mean())
            rows.append(dict(score_column=c, cell_type=ct, frac_zero=fz,
                             statistic="mean" if fz > zero_frac else "median"))
    return pd.DataFrame(rows, columns=["score_column", "cell_type", "frac_zero", "statistic"])


def _subtype_plan(obs: pd.DataFrame, max_cov=MAX_SUBTYPE_COVARIATES) -> dict:
    """{cell_type: [subtype for subfrac_1, subfrac_2, ...]}; most abundant = reference."""
    if "cell_subtype" not in obs:
        return {}
    plan = {}
    for ct, s in obs.groupby("cell_type", observed=True)["cell_subtype"]:
        order = s.astype(str).value_counts().index.tolist()
        if len(order) > 1:
            plan[str(ct)] = order[1:1 + max_cov]
    return plan


def donor_celltype_summary(obs: pd.DataFrame, score_cols, min_cells=MIN_CELLS_DEFAULT,
                           statistic="auto") -> pd.DataFrame:
    """One row per donor x cell type.

    Donor x cell-type cells with fewer than `min_cells` cells are MISSING, never
    zero: a population that was not recovered carries no information about its
    state.

    statistic : "auto" (default; see `choose_statistic`), "median", "mean" or
        "winsorized_mean", applied to UCell scores. Single-gene features are
        pseudobulk log1p(mean CP10K); ratio features are logit of pooled counts.

    Also returned, for use as covariates:
      * ``log_genes`` - log of the median detected genes per cell. UCell scores
        depend on how many genes a cell detects, and depth can differ between
        age groups (sample quality).
      * ``subfrac_1``, ``subfrac_2`` - each donor's fraction of the 2nd and 3rd
        most common subtype within the population (from ``obs["cell_subtype"]``),
        so a shift in subtype mix (e.g. goblet vs club) is not read as cell state.
        Which subtype each column holds is in ``summary.attrs["subtype_covariates"]``.

    ``pd.DataFrame(summary.attrs["aggregation"])`` lists the statistic used per
    module x cell type when statistic="auto".
    """
    if statistic == "auto":
        plan = choose_statistic(obs, score_cols)
        stat_for = {(r.score_column, r.cell_type): r.statistic for r in plan.itertuples()}
    else:
        if statistic not in _STATS:
            raise ValueError(f"statistic must be 'auto' or one of {sorted(_STATS)}")
        plan = None
        stat_for = {}
    sub_plan = _subtype_plan(obs)
    has_genes = "n_genes" in obs
    rows = []
    for (d, ct), g in obs.groupby(["donor", "cell_type"], observed=True):
        n = len(g)
        rec = dict(donor=d, cell_type=ct, n_cells=n, eligible=bool(n >= min_cells),
                   age_bin=g.age_bin.iloc[0], sex=g.sex.iloc[0], batch=g.batch.iloc[0])
        if has_genes:
            rec["log_genes"] = float(np.log(max(np.median(g.n_genes.values), 1)))
        subs = sub_plan.get(str(ct), [])
        for k in range(MAX_SUBTYPE_COVARIATES):
            rec[f"subfrac_{k + 1}"] = (float((g.cell_subtype.astype(str) == subs[k]).mean())
                                       if k < len(subs) else np.nan)
        for c in score_cols:
            if n < min_cells:
                rec[c] = np.nan
            else:
                fn = _STATS[stat_for.get((c, ct), "median") if statistic == "auto" else statistic]
                rec[c] = _cell_value(g, c, fn)
        rows.append(rec)
    donors = obs[["donor", "age_bin", "sex", "batch"]].drop_duplicates("donor")
    have = {(r["donor"], r["cell_type"]) for r in rows}
    for _, d in donors.iterrows():
        for ct in obs.cell_type.cat.categories:
            if (d.donor, ct) not in have:
                rows.append(dict(donor=d.donor, cell_type=ct, n_cells=0, eligible=False,
                                 age_bin=d.age_bin, sex=d.sex, batch=d.batch,
                                 **{c: np.nan for c in score_cols}))
    out = pd.DataFrame(rows)
    if isinstance(obs.age_bin.dtype, pd.CategoricalDtype):
        # keep the user's age-group order for contrasts and figures
        out["age_bin"] = pd.Categorical(out.age_bin.astype(str),
                                        categories=[str(c) for c in obs.age_bin.cat.categories])
    # attrs must hold plain Python objects (a DataFrame here breaks pandas concat/repr);
    # pd.DataFrame(summary.attrs["aggregation"]) gives the table back.
    out.attrs["aggregation"] = [] if plan is None else plan.to_dict("records")
    out.attrs["subtype_covariates"] = sub_plan
    return out


def eligibility_table(summary: pd.DataFrame) -> pd.DataFrame:
    return (summary.groupby("cell_type", observed=True)
            .agg(n_donors=("donor", "nunique"), n_eligible=("eligible", "sum"),
                 median_cells=("n_cells", "median"), min_cells=("n_cells", "min"),
                 max_cells=("n_cells", "max")).reset_index())


def pooled_scores(obs: pd.DataFrame, score_cols, cell_types, min_cells=MIN_CELLS_DEFAULT,
                  label="pooled") -> pd.DataFrame:
    """Donor summaries that IGNORE cell identity, for the confounding check only.

    Comparing this with the within-cell-type result shows how much of an
    apparent age effect is composition rather than cell state.
    """
    sub = obs[obs.cell_type.isin(cell_types)]
    rows = []
    for d, g in sub.groupby("donor", observed=True):
        if len(g) < min_cells:
            continue
        rec = dict(donor=d, cell_type=label, n_cells=len(g), eligible=True,
                   age_bin=g.age_bin.iloc[0], sex=g.sex.iloc[0], batch=g.batch.iloc[0])
        for c in score_cols:
            rec[c] = _cell_value(g, c, np.median)
        rows.append(rec)
    return pd.DataFrame(rows)
