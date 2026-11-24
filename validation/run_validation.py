"""Run NasAge on one external dataset (healthy donors only) and save every validation output.

Usage:
    python validation/run_validation.py --config validation/configs/<dataset>.json

The config (JSON) says where the h5ad is and which columns to use:
{
  "name": "deprez2020",
  "h5ad": "D:/research/nasage_validation/data/deprez2020.h5ad",
  "counts": "raw" | "X" | "<layer name>",     # where raw integer counts are
  "symbol_col": "feature_name",               # var column with gene symbols (null if var_names are symbols)
  "subset": {"disease": ["normal"], "tissue": ["nose"]},   # obs filters, all must hold
  "exclude_studies": {"<obs col>": ["<value>"]},           # e.g. remove Yoshida 2022 donors
  "donor_key": "donor_id", "age_key": "age_group", "age_bins": null, "age_order": null,
  "age_from_development_stage": false,        # CELLxGENE: parse "NN-year-old human stage"
  "fine_key": "cell_type",                    # fine annotation (kept as subtype_key)
  "population_map": {"basal cell": "basal", ...},
  "batch_key": null, "sex_key": "sex",
  "reference": "18-39", "min_cells": 10, "n_null": 50, "n_perm": 200,
  "out": "D:/research/nasage_validation/results/deprez2020"
}
"""
from __future__ import annotations

import argparse
import json
import re
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd


def log(msg):
    print(time.strftime("%H:%M:%S"), msg, flush=True)


def load(cfg):
    import anndata as ad
    log(f"reading {cfg['h5ad']}")
    a = ad.read_h5ad(cfg["h5ad"])
    for col, keep in (cfg.get("subset") or {}).items():
        a = a[a.obs[col].astype(str).isin([str(k) for k in keep]).values]
    for col, drop in (cfg.get("exclude_studies") or {}).items():
        a = a[~a.obs[col].astype(str).isin([str(k) for k in drop]).values]
    a = a.copy()
    counts = cfg.get("counts", "X")
    if counts == "raw":
        a.layers["counts"] = a.raw[:, a.var_names].X.copy() if a.raw.n_vars != a.n_vars \
            else a.raw.X.copy()
    elif counts == "X":
        a.layers["counts"] = a.X.copy()
    else:
        a.layers["counts"] = a.layers[counts].copy()
    a.raw = None
    if cfg.get("symbol_col"):
        a.var_names = a.var[cfg["symbol_col"]].astype(str).values
        a.var_names_make_unique()
    if cfg.get("age_map"):                       # e.g. decade labels -> analysis groups
        src_col, mapping = cfg["age_map"]["column"], cfg["age_map"]["map"]
        a.obs[cfg["age_key"]] = pd.Categorical(
            a.obs[src_col].astype(str).map(mapping), categories=cfg.get("age_order"))
    if cfg.get("age_numeric_from"):              # e.g. a string column of ages in years
        a.obs[cfg["age_key"]] = pd.to_numeric(a.obs[cfg["age_numeric_from"]].astype(str),
                                              errors="coerce")
    if cfg.get("age_from_development_stage"):
        m = a.obs["development_stage"].astype(str).str.extract(r"(\d+)-year-old")[0]
        a.obs["age_years"] = pd.to_numeric(m, errors="coerce")
    a.obs["population"] = a.obs[cfg["fine_key"]].astype(str).map(cfg["population_map"])
    log(f"{a.n_obs} cells, {a.obs[cfg['donor_key']].nunique()} donors after subsetting")
    return a


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    cfg = json.loads(Path(ap.parse_args().config).read_text())
    out = Path(cfg["out"]); out.mkdir(parents=True, exist_ok=True)
    warnings.simplefilter("default")
    import nasage as nas

    a = load(cfg)
    keys = dict(donor_key=cfg["donor_key"], celltype_key="population", age_key=cfg["age_key"],
                age_bins=cfg.get("age_bins"), age_order=cfg.get("age_order"),
                batch_key=cfg.get("batch_key"), sex_key=cfg.get("sex_key"))
    design = nas.pp.check_design(a, **keys)
    design.to_csv(out / "donor_design.csv", index=False)
    log("donors per age group: " + str(design.age_bin.value_counts().to_dict()))

    log("QC"); nas.pp.filter_cells(a, celltype_key="population", layer="counts")
    log("scoring"); nas.tl.score_modules(a, layer="counts", n_jobs=8)
    log("donor summary")
    nas.tl.donor_summary(a, subtype_key=cfg["fine_key"], min_cells=cfg.get("min_cells", 10), **keys)
    log("age effects"); om = nas.tl.age_effects(a, reference=cfg["reference"])
    log("composition"); nas.tl.composition(a)
    log("permutation null"); nas.tl.permutation_null(a, n_perm=cfg.get("n_perm", 200))
    log("random-gene-set null")
    nas.tl.null_control(a, n_sets=cfg.get("n_null", 50), n_jobs=8,
                        tiers=tuple(cfg.get("null_tiers", ("primary", "exploratory"))))
    nas.write_report(a, out, label=cfg["name"])

    st = a.uns["nasage"]
    corr = {ct: nas.tl.module_correlation(a, ct) for ct in ("ciliated", "secretory", "basal")
            if ct in set(st["donor_summary"].cell_type.astype(str))}
    for ct, c in corr.items():
        c.to_csv(out / f"module_correlation_{ct}.csv")
    st["contrasts"].to_csv(out / "age_effects_contrasts.csv", index=False)
    om2 = om.merge(st["null_control"][["module_id", "cell_type", "empirical_p", "null_frac_p05"]],
                   on=["module_id", "cell_type"], how="left")
    om2.to_csv(out / "age_effects_with_null.csv", index=False)
    summary = dict(
        name=cfg["name"], n_cells=int(a.n_obs), n_donors=int(design.shape[0]),
        donors_per_group=design.age_bin.value_counts().to_dict(),
        primary_tests=int(((om.tier == "primary") & om.estimable.astype(bool)).sum()),
        primary_significant=om[(om.tier == "primary") & (om.p_adj_BH.replace("", np.nan)
                                                         .astype(float) < 0.05)]
        [["module_id", "cell_type"]].values.tolist(),
        permutation=dict(st.get("permutation_rates", {})),
        null_generic_signal=st["null_control"].loc[
            st["null_control"].generic_age_signal.astype(str) == "True",
            ["module_id", "cell_type"]].values.tolist())
    (out / "summary.json").write_text(json.dumps(summary, indent=1, default=str))
    log("done: " + json.dumps(summary, default=str)[:600])


if __name__ == "__main__":
    main()
