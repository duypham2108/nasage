"""Internal: turn the user's AnnData columns into the canonical frame the backend uses.

Users keep their own column names and pass them as ``*_key`` arguments; nothing
is renamed in ``adata.obs``.
"""
from __future__ import annotations

import re
import warnings

import numpy as np
import pandas as pd
import scipy.sparse as sp

UNS_KEY = "nasage"


def store(adata) -> dict:
    """The package's results dict in ``adata.uns``."""
    return adata.uns.setdefault(UNS_KEY, {})


def get_matrix(adata, layer):
    if layer is None:
        return adata.X
    if layer not in adata.layers:
        raise KeyError(f"layer '{layer}' not in adata.layers {list(adata.layers)}; "
                       "pass the layer that holds raw counts (layer=None for adata.X)")
    return adata.layers[layer]


def check_counts(adata, layer) -> None:
    """Refuse non-integer input: QC metrics, CP10K and the IgA ratio need raw counts."""
    X = get_matrix(adata, layer)
    head = X[:200]
    vals = head.data if sp.issparse(head) else np.asarray(head).ravel()
    if vals.size and not np.allclose(vals, np.round(vals)):
        where = "adata.X" if layer is None else f"adata.layers['{layer}']"
        raise ValueError(f"{where} does not hold integer counts. Pass layer= the layer with "
                         "raw counts (e.g. layer='counts').")


def natural_key(s):
    """Sort '5-11' before '12-17' (numbers compared as numbers)."""
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", str(s))]


def _age_groups(obs, age_key, age_bins, age_order, donor):
    col = obs[age_key]
    if age_bins is not None:
        age = pd.to_numeric(col, errors="coerce")
        lab = pd.Series(None, index=obs.index, dtype="object")
        for lo, hi, name in age_bins:
            lab[(age >= lo) & (age <= hi)] = str(name)
        order = [str(b[2]) for b in age_bins]
        what = "an age outside age_bins (or missing age)"
    else:
        if pd.api.types.is_numeric_dtype(col) and col.nunique() > 12:
            raise ValueError(f"'{age_key}' looks like numeric age; pass age_bins=[(lo, hi, label), ...] "
                             "or point age_key at an age-group column")
        lab = col.astype("object").where(col.notna(), None).map(
            lambda x: None if x is None else str(x))
        present = [str(x) for x in pd.unique(lab.dropna())]
        if age_order is not None:
            order = [str(x) for x in age_order]
            missing = sorted(set(present) - set(order))
            if missing:
                raise ValueError(f"age groups {missing} are not listed in age_order {order}")
        elif isinstance(col.dtype, pd.CategoricalDtype):
            order = [str(x) for x in col.cat.categories if str(x) in present]
        else:
            order = sorted(present, key=natural_key)
        what = "no age group"
    out = lab.isna().values
    if out.any():
        donors = sorted(map(str, pd.unique(donor[out])))
        warnings.warn(f"{len(donors)} donors ({int(out.sum())} cells) have {what} and are "
                      f"excluded: {donors[:20]}{' ...' if len(donors) > 20 else ''}")
    return pd.Categorical(lab, categories=order)


def _check_donor_constant(obs: pd.DataFrame, cols) -> None:
    """Donor-level attributes must be constant within a donor."""
    bad = {}
    for c in cols:
        n = obs.groupby("donor", observed=True)[c].nunique(dropna=False)
        if (n > 1).any():
            bad[c] = sorted(map(str, n.index[n > 1]))
    if bad:
        detail = "; ".join(f"{c}: {v[:10]}{' ...' if len(v) > 10 else ''}" for c, v in bad.items())
        raise ValueError(
            "donor-level variables vary within donor (" + detail + "). The model has one row "
            "per donor x cell type; give each sample its own donor ID or harmonise the metadata.")


def canonical_obs(adata, params: dict, extra_cols=()) -> pd.DataFrame:
    """Cells usable for analysis, with columns donor, cell_type, age_bin, batch, sex,
    optional cell_subtype and n_genes, plus `extra_cols` copied from adata.obs.

    Cells with no cell type (NaN in celltype_key) or no age group are excluded.
    """
    o = adata.obs
    for k in ("donor_key", "celltype_key", "age_key"):
        if params.get(k) is None:
            raise KeyError(f"{k} is required")
        if params[k] not in o:
            raise KeyError(f"column '{params[k]}' ({k}) not in adata.obs")
    df = pd.DataFrame(index=o.index)
    df["donor"] = o[params["donor_key"]].astype(str).values
    df["cell_type"] = o[params["celltype_key"]].values
    df["age_bin"] = _age_groups(o, params["age_key"], params.get("age_bins"),
                                params.get("age_order"), df["donor"].values)
    for name in ("batch", "sex"):
        key = params.get(f"{name}_key")
        if key is None:
            warnings.warn(f"no {name}_key given; {name} will not be adjusted for. Unmodelled "
                          "technical or biological structure inflates false positives.")
            df[name] = "unknown"
        else:
            if key not in o:
                raise KeyError(f"column '{key}' ({name}_key) not in adata.obs")
            df[name] = o[key].astype(str).values
    if params.get("subtype_key"):
        df["cell_subtype"] = o[params["subtype_key"]].astype(str).values
    for c in ("nis_n_genes", "n_genes_by_counts"):
        if c in o:
            df["n_genes"] = o[c].values
            break
    age_raw = o[params["age_key"]].values
    _check_donor_constant(df.assign(age_raw=age_raw),
                          ["age_raw", "batch", "sex"])
    keep = df.cell_type.notna().values & df.age_bin.notna().values
    n_ct = int(df.cell_type.isna().sum())
    if n_ct:
        warnings.warn(f"{n_ct} cells have no value in '{params['celltype_key']}' and are excluded")
    df = df[keep].copy()
    df["cell_type"] = pd.Categorical(df.cell_type.astype(str))
    for c in ("donor", "batch", "sex"):
        df[c] = df[c].astype("category")
    for c in extra_cols:
        df[c] = o.loc[df.index, c].values
    return df


def donor_table(obs: pd.DataFrame) -> pd.DataFrame:
    """Donor-level design table plus confounding warnings."""
    d = (obs.groupby("donor", observed=True)
         .agg(age_bin=("age_bin", "first"), batch=("batch", "first"),
              sex=("sex", "first"), n_cells=("age_bin", "size")))
    tab = pd.crosstab(d.age_bin, d.batch)
    if (tab > 0).sum(1).max() == 1 and tab.shape[1] > 1:
        warnings.warn("every age group sits in a single batch: batch and age are completely "
                      "confounded and the age term is not identifiable.")
    small = d.age_bin.value_counts()
    if (small < 3).any():
        warnings.warn(f"age groups with fewer than 3 donors: {small[small < 3].to_dict()}; "
                      "these are not estimable.")
    return d.reset_index()
