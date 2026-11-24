"""Derive Set A: the nasal mucosal programme engaged by a live attenuated intranasal vaccine.

Implements docs/PREREG_module_derivation.md, locked in commit b99961e. No threshold or rule
here may be changed after discovery begins; departures go in the results file as dated
deviations.

  Discovery   GSE230494  paired baseline -> 72 h post-LAIV, no bacterial challenge
  Replication GSE164649  LAIV vs TIV at D6, D9     (randomised; day numbers are relative to
              GSE117580  LAIV vs TIV at D2, D9      pneumococcal inoculation, hence arm contrast)

Entry: mean CPM >= 1; discovery BH < 0.05 AND |log2FC| >= 0.5; same direction and nominal
p < 0.05 in >= 1 of the four arm contrasts; cap 100 by |log2FC|; fail if < 10 survive.

Usage: python validation/derive_setA_laiv_module.py <out_dir>
"""
from __future__ import annotations

import gzip
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

import nasage as nas

THW = Path("D:/research/nasage_validation/data/thwaites2023")
LIV = Path("D:/research/nasage_validation/data/liverpool")
SCRATCH = Path("C:/Users/duy/AppData/Local/Temp/claude/a802a2006c0920b5/"
               "dfa6adbf-4875-4870-80ea-436b025fe68a/scratchpad")
ZIEGLER = "D:/research/nasage_validation/data/ziegler2021_np.h5ad"
MIN_CPM, BH_MAX, LFC_MIN, CAP, MIN_GENES = 1.0, 0.05, 0.5, 100, 10


def symbol_map(index: pd.Index) -> pd.Series:
    import anndata as ad
    var = ad.read_h5ad(ZIEGLER, backed="r").var
    return var["feature_name"].astype(str).reindex(
        pd.Index(index.astype(str).str.split(".").str[0]))


def to_symbols(C: pd.DataFrame) -> pd.DataFrame:
    sym = symbol_map(C.index)
    C = C[sym.notna().values]
    return C.groupby(sym[sym.notna()].values).sum()


def logcpm(C: pd.DataFrame):
    cpm = C / C.sum(0) * 1e6
    return np.log2(cpm + 1), cpm.mean(1)


def vectorised_ols(Y: np.ndarray, X: np.ndarray, col: int):
    """Per-gene OLS; return (beta, t, p) for design column `col`. Y is genes x samples."""
    n, p = X.shape
    pinv = np.linalg.pinv(X)
    B = Y @ pinv.T                                   # genes x p
    resid = Y - B @ X.T
    dof = n - np.linalg.matrix_rank(X)
    s2 = (resid ** 2).sum(1) / dof
    xtx_inv = np.linalg.pinv(X.T @ X)
    se = np.sqrt(np.maximum(s2 * xtx_inv[col, col], 1e-300))
    beta = B[:, col]
    t = beta / se
    return beta, t, 2 * stats.t.sf(np.abs(t), dof)


# ----------------------------------------------------------------- discovery
def discovery():
    cols = {}
    for f in sorted((THW / "counts").glob("*_counts.txt.gz")):
        name = re.sub(r"^GSM\d+_|_counts\.txt\.gz$", "", f.name)
        with gzip.open(f, "rt") as fh:
            cols[name] = pd.read_csv(fh, sep="\t", header=None, index_col=0,
                                     skiprows=1).iloc[:, 0]
    C = to_symbols(pd.DataFrame(cols))
    meta = pd.DataFrame(index=C.columns)
    meta["participant"] = meta.index.str.extract(r"(LAIV_\d+)", expand=False).values
    meta["tp"] = np.where(meta.index.str.contains(r"_(?:1|Visit1)$"), "baseline", "SD3")
    C = C.loc[:, C.sum(0) >= 5e6]                    # same depth rule as the earlier prereg
    meta = meta.loc[C.columns]
    pairs = [p for p, g in meta.groupby("participant")
             if set(g.tp) == {"baseline", "SD3"}]
    print(f"discovery: {C.shape[1]} samples, {len(pairs)} paired participants")

    lg, meancpm = logcpm(C)
    keep = meancpm >= MIN_CPM
    lg = lg[keep]
    b = lg[[meta.index[(meta.participant == p) & (meta.tp == "baseline")][0] for p in pairs]]
    s = lg[[meta.index[(meta.participant == p) & (meta.tp == "SD3")][0] for p in pairs]]
    d = s.to_numpy() - b.to_numpy()
    t, p = stats.ttest_rel(s.to_numpy(), b.to_numpy(), axis=1)
    out = pd.DataFrame({"gene": lg.index, "log2FC": d.mean(1), "t": t, "p": p})
    out["p_adj_BH"] = multipletests(np.nan_to_num(out.p, nan=1.0), method="fdr_bh")[1]
    out["n_pairs"] = len(pairs)
    print(f"  genes tested: {len(out)}; BH<0.05: {(out.p_adj_BH < BH_MAX).sum()}; "
          f"and |log2FC|>={LFC_MIN}: {((out.p_adj_BH < BH_MAX) & (out.log2FC.abs() >= LFC_MIN)).sum()}")
    return out


# --------------------------------------------------------------- replication
def liverpool_tables():
    s164 = pd.read_csv(SCRATCH / "GSE164649_samples.csv")
    s164["sid"] = s164.title.str.extract(r"\(([^)]+)\)")
    s164["arm"] = s164.vaccine
    s164["carr"] = s164.carriage
    s117 = pd.read_csv(SCRATCH / "GSE117580_samples.csv")
    s117["sid"] = s117.title.str.extract(r"_([A-Z]?\d+)$")
    s117["arm"] = np.where(s117.vaccine.str.contains("LAIV"), "LAIV", "TIV")
    s117["carr"] = s117.carriage_status

    with gzip.open(LIV / "GSE164649_counts.txt.gz", "rt") as fh:
        C164 = pd.read_csv(fh, sep="\t", index_col=0)
    with gzip.open(LIV / "GSE117580_counts.csv.gz", "rt") as fh:
        C117 = pd.read_csv(fh, index_col=0)
    return (to_symbols(C164), s164.set_index("sid")), (to_symbols(C117), s117.set_index("sid"))


def arm_contrast(C, meta, timepoint, candidates, label):
    sel = meta.index[(meta.timepoint == timepoint)]
    sel = [s for s in sel if s in C.columns]
    md, sub = meta.loc[sel], C[sel]
    lg, _ = logcpm(sub)
    lg = lg.reindex(candidates).dropna(how="all")
    X = pd.DataFrame({"intercept": 1.0,
                      "arm_LAIV": (md.arm.values == "LAIV").astype(float),
                      "carriage_POS": pd.Series(md.carr.values).astype(str)
                                        .str.upper().str.startswith("POS").astype(float).values,
                      "log10_depth": np.log10(sub.sum(0).values)}, index=md.index)
    if X.arm_LAIV.nunique() < 2:
        return None
    beta, t, p = vectorised_ols(lg.to_numpy(), X.to_numpy(), col=1)
    n_laiv, n_tiv = int(X.arm_LAIV.sum()), int((1 - X.arm_LAIV).sum())
    print(f"  {label}: LAIV {n_laiv} vs TIV {n_tiv}, {lg.shape[0]} candidate genes")
    return pd.DataFrame({f"{label}_beta": beta, f"{label}_p": p}, index=lg.index)


def main():
    out = Path(sys.argv[1]); out.mkdir(parents=True, exist_ok=True)
    disc = discovery()
    cand = disc[(disc.p_adj_BH < BH_MAX) & (disc.log2FC.abs() >= LFC_MIN)].copy()
    cand = cand.set_index("gene")
    if cand.empty:
        raise SystemExit("no discovery candidates")

    (C164, s164), (C117, s117) = liverpool_tables()
    reps = []
    for C, md, tps, tag in ((C164, s164, ["D6", "D9"], "GSE164649"),
                            (C117, s117, ["D2", "D9"], "GSE117580")):
        for tp in tps:
            r = arm_contrast(C, md, tp, cand.index, f"{tag}_{tp}")
            if r is not None:
                reps.append(r)
    R = cand.join(pd.concat(reps, axis=1))

    # replication: same direction AND nominal p < 0.05 in >= 1 contrast
    labels = [c[:-5] for c in R.columns if c.endswith("_beta")]
    ok = pd.DataFrame({lab: (np.sign(R[f"{lab}_beta"]) == np.sign(R.log2FC))
                            & (R[f"{lab}_p"] < 0.05) for lab in labels})
    R["n_contrasts_replicated"] = ok.sum(1)
    R["replicated"] = R.n_contrasts_replicated >= 1
    print(f"candidates: {len(R)}; replicated in >=1 contrast: {int(R.replicated.sum())}")
    print(ok.sum().to_dict())

    setA = R[R.replicated].reindex(R[R.replicated].log2FC.abs().sort_values(ascending=False).index)
    capped = len(setA) > CAP
    setA = setA.head(CAP)
    if len(setA) < MIN_GENES:
        R.to_csv(out / "setA_candidates_FAILED.csv")
        raise SystemExit(f"SET A FAILS: only {len(setA)} genes survived (< {MIN_GENES})")
    print(f"Set A: {len(setA)} genes" + (f" (capped from {int(R.replicated.sum())})" if capped else ""))

    # ---- composition check against an independent reference atlas
    import anndata as ad
    a = ad.read_h5ad(ZIEGLER)
    a.var_names = a.var["feature_name"].astype(str).values
    a.var_names_make_unique()
    ct = a.obs["Coarse_Cell_Annotations"].astype(str)
    genes = [g for g in setA.index if g in a.var_names]
    X = a[:, genes].X
    X = np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X)
    prof = pd.DataFrame(X, columns=genes).groupby(ct.values).mean()
    prof = prof.div(prof.sum(0).replace(0, np.nan), axis=1)          # fraction per cell type
    top = prof.idxmax(0)
    conc = prof.max(0)
    comp = pd.DataFrame({"top_celltype": top, "frac_in_top": conc})
    dominant = comp.top_celltype.value_counts()
    frac_dom = float(dominant.iloc[0] / len(comp)) if len(comp) else np.nan
    flagged = frac_dom > 0.5
    print(f"composition: {len(genes)}/{len(setA)} genes in reference; "
          f"most common top cell type = {dominant.index[0]} for {frac_dom:.0%} of genes "
          f"-> {'FLAGGED composition_dominated' if flagged else 'not flagged'}")

    # ---- overlap with the 15 curated modules
    M = {m: set(d["genes"]) for m, d in nas.load_modules().items()}
    s = set(setA.index)
    jac = pd.Series({m: len(s & g) / len(s | g) for m, g in M.items()}).sort_values(ascending=False)
    m05 = jac.get("M05_interferon_stimulated_genes", 0.0)
    print("\ntop Jaccard overlaps with curated modules:\n", jac.head(5).round(3).to_string())
    print(f"M05 Jaccard = {m05:.3f} -> "
          f"{'CONFIRMATION OF M05 (not a new module)' if m05 >= 0.5 else 'treated as a distinct module'}")
    novel = sorted(s - set().union(*M.values()))
    print(f"genes not in any curated module: {len(novel)}")

    setA.to_csv(out / "setA_genes.csv"); R.to_csv(out / "setA_all_candidates.csv")
    comp.to_csv(out / "setA_composition.csv"); prof.to_csv(out / "setA_celltype_profile.csv")
    json.dump({"n_genes": len(setA), "capped": bool(capped),
               "n_replicated": int(R.replicated.sum()),
               "contrasts": {k: int(v) for k, v in ok.sum().items()},
               "composition_dominant_celltype": str(dominant.index[0]),
               "composition_frac": frac_dom, "composition_flagged": bool(flagged),
               "jaccard_M05": float(m05), "jaccard_top": jac.head(5).round(4).to_dict(),
               "n_novel_genes": len(novel), "novel_genes": novel},
              open(out / "setA_summary.json", "w"), indent=1)
    print("\ntop 25 by |log2FC|:\n",
          setA.head(25)[["log2FC", "p_adj_BH", "n_contrasts_replicated"]].round(3).to_string())


if __name__ == "__main__":
    main()
