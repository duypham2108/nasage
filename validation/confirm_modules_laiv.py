"""Confirm and refine the 15 curated modules against live mucosal vaccination.

Implements docs/PREREG_module_confirmation.md, locked in commit abbde23.

  Discovery   GSE230494  paired baseline -> 72 h post-LAIV (20 pairs), no bacterial challenge
  Replication GSE164649  LAIV vs TIV at D6, D9      (randomised; arm contrast, because day
              GSE117580  LAIV vs TIV at D2, D9       numbering is relative to pneumococcus)

ENGAGED = discovery BH < 0.05, replicated (same direction, p < 0.05) in >= 1 arm contrast, and
the same verdict under BOTH scoring methods.

Implementation notes for "as the package defines them" (§2):
  * multi-gene modules: scoring (a) mean of z-scored log2(CPM+1), (b) UCell on counts;
  * single-gene features (M03 PIGR, M03b FCGRT): both methods reduce to the same quantity, so
    the second scoring is log2(CPM+1); condition 3 is then automatically met if (a) passes;
  * M10 is the logit IgA share of immunoglobulin heavy-chain counts, computed directly.
  * mean-z is computed across all analysed samples of a dataset, not per timepoint subset, so
    the scale is stable across contrasts within that dataset.

Usage: python validation/confirm_modules_laiv.py <out_dir>
"""
from __future__ import annotations

import gzip
import json
import re
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from statsmodels.stats.multitest import multipletests

import nasage as nas

warnings.simplefilter("ignore")
THW = Path("D:/research/nasage_validation/data/thwaites2023")
LIV = Path("D:/research/nasage_validation/data/liverpool")
SCRATCH = Path("C:/Users/duy/AppData/Local/Temp/claude/a802a2006c0920b5/"
               "dfa6adbf-4875-4870-80ea-436b025fe68a/scratchpad")
ZIEGLER = "D:/research/nasage_validation/data/ziegler2021_np.h5ad"
MIN_CPM, BH_MAX, MIN_CORE = 1.0, 0.05, 5


def symbols(C):
    import anndata as ad
    var = ad.read_h5ad(ZIEGLER, backed="r").var
    s = var["feature_name"].astype(str).reindex(
        pd.Index(C.index.astype(str).str.split(".").str[0]))
    return C[s.notna().values].groupby(s[s.notna()].values).sum()


def prep(C):
    """-> (log2CPM of expressed genes, z-scores, raw counts) for one dataset."""
    C = C.loc[:, C.sum(0) >= 5e6]
    cpm = C / C.sum(0) * 1e6
    lg = np.log2(cpm[cpm.mean(1) >= MIN_CPM] + 1)
    sd = lg.std(1, ddof=1)
    Z = lg[sd > 0].sub(lg[sd > 0].mean(1), axis=0).div(sd[sd > 0], axis=0)
    return lg, Z, C


def module_scores(lg, Z, C, M):
    """Return {'meanz': df, 'alt': df} of samples x modules ('alt' = UCell / log2CPM / ratio)."""
    mz, alt = {}, {}
    for m, d in M.items():
        g = [x for x in d["genes"] if x in Z.index]
        if not g:
            continue
        mz[m] = Z.loc[g].mean(0)
    # single-gene features: second method is log2CPM (monotone in mean-z)
    for m, d in M.items():
        if len(d["genes"]) == 1 and d["genes"][0] in lg.index:
            alt[m] = lg.loc[d["genes"][0]]
    # IgA ratio feature
    iga = [g for g in M["M10_iga_fraction"]["genes"] if g in lg.index]
    a = [g for g in iga if g.startswith("IGHA")]
    if a and len(iga) > len(a):
        cpm = 2 ** lg - 1
        frac = cpm.loc[a].sum(0) / cpm.loc[iga].sum(0).replace(0, np.nan)
        alt["M10_iga_fraction"] = np.log((frac + .01) / (1 - frac + .01))
        mz["M10_iga_fraction"] = alt["M10_iga_fraction"]
    try:
        import anndata as ad
        from pyucell import compute_ucell_scores
        from nasage.scoring import MAX_RANK
        u = ad.AnnData(X=C.T.values.astype(np.float32), obs=pd.DataFrame(index=C.columns),
                       var=pd.DataFrame(index=C.index))
        sigs = {m: [g for g in d["genes"] if g in C.index] for m, d in M.items()
                if len(d["genes"]) > 1 and m != "M10_iga_fraction"}
        compute_ucell_scores(u, {k: v for k, v in sigs.items() if v}, max_rank=MAX_RANK,
                             suffix="", n_jobs=1)
        for m in sigs:
            if m in u.obs:
                alt[m] = u.obs[m]
    except Exception as exc:
        print("  UCell failed:", exc)
    return pd.DataFrame(mz), pd.DataFrame(alt)


def discovery_tables():
    cols = {}
    for f in sorted((THW / "counts").glob("*_counts.txt.gz")):
        n = re.sub(r"^GSM\d+_|_counts\.txt\.gz$", "", f.name)
        with gzip.open(f, "rt") as fh:
            cols[n] = pd.read_csv(fh, sep="\t", header=None, index_col=0, skiprows=1).iloc[:, 0]
    C = symbols(pd.DataFrame(cols))
    lg, Z, C = prep(C)
    md = pd.DataFrame(index=C.columns)
    md["participant"] = md.index.str.extract(r"(LAIV_\d+)", expand=False).values
    md["tp"] = np.where(md.index.str.contains(r"_(?:1|Visit1)$"), "baseline", "SD3")
    pairs = [p for p, g in md.groupby("participant") if set(g.tp) == {"baseline", "SD3"}]
    return lg, Z, C, md, pairs


def paired_test(S, md, pairs):
    bi = [md.index[(md.participant == p) & (md.tp == "baseline")][0] for p in pairs]
    si = [md.index[(md.participant == p) & (md.tp == "SD3")][0] for p in pairs]
    rows = []
    for m in S.columns:
        b, s = S.loc[bi, m].to_numpy(), S.loc[si, m].to_numpy()
        ok = np.isfinite(b) & np.isfinite(s)
        if ok.sum() < 5:
            continue
        d = s[ok] - b[ok]
        rows.append(dict(module=m, n_pairs=int(ok.sum()), delta=float(d.mean()),
                         frac_up=float((d > 0).mean()),
                         p=float(stats.ttest_rel(s[ok], b[ok]).pvalue),
                         p_wilcoxon=float(stats.wilcoxon(d).pvalue) if len(set(d)) > 1 else np.nan))
    R = pd.DataFrame(rows)
    R["p_adj_BH"] = multipletests(R.p.fillna(1), method="fdr_bh")[1]
    return R.sort_values("p")


def liverpool():
    out = {}
    s164 = pd.read_csv(SCRATCH / "GSE164649_samples.csv")
    s164["sid"] = s164.title.str.extract(r"\(([^)]+)\)")
    s164["arm"], s164["carr"] = s164.vaccine, s164.carriage
    with gzip.open(LIV / "GSE164649_counts.txt.gz", "rt") as fh:
        out["GSE164649"] = (symbols(pd.read_csv(fh, sep="\t", index_col=0)),
                            s164.set_index("sid"), ["D6", "D9"])
    s117 = pd.read_csv(SCRATCH / "GSE117580_samples.csv")
    s117["sid"] = s117.title.str.extract(r"_([A-Z]?\d+)$")
    s117["arm"] = np.where(s117.vaccine.str.contains("LAIV"), "LAIV", "TIV")
    s117["carr"] = s117.carriage_status
    with gzip.open(LIV / "GSE117580_counts.csv.gz", "rt") as fh:
        out["GSE117580"] = (symbols(pd.read_csv(fh, index_col=0)),
                            s117.set_index("sid"), ["D2", "D9"])
    return out


def arm_model(y, md, depth):
    X = pd.DataFrame({"arm_LAIV": (md.arm.values == "LAIV").astype(float),
                      "carr_POS": md.carr.astype(str).str.upper().str.startswith("POS")
                                    .astype(float).values,
                      "depth": np.log10(depth.values)}, index=md.index)
    X = sm.add_constant(X)
    ok = np.isfinite(y.values)
    if ok.sum() < 8 or X.arm_LAIV[ok].nunique() < 2:
        return np.nan, np.nan
    m = sm.OLS(y.values[ok], X[ok]).fit()
    return float(m.params["arm_LAIV"]), float(m.pvalues["arm_LAIV"])


def main():
    out = Path(sys.argv[1]); out.mkdir(parents=True, exist_ok=True)
    M = nas.load_modules()
    lg, Z, C, md, pairs = discovery_tables()
    print(f"discovery: {C.shape[1]} samples, {len(pairs)} pairs, {lg.shape[0]} genes")
    Smz, Salt = module_scores(lg, Z, C, M)
    disc = {"meanz": paired_test(Smz, md, pairs), "alt": paired_test(Salt, md, pairs)}
    for k, v in disc.items():
        print(f"\n--- discovery ({k}): BH<0.05 -> "
              f"{sorted(v.loc[v.p_adj_BH < BH_MAX, 'module'])}")

    # ---- replication: arm contrasts
    rep = {}
    for tag, (Cl, mdl, tps) in liverpool().items():
        lgl, Zl, Cl2 = prep(Cl)
        mdl = mdl.loc[[s for s in mdl.index if s in Cl2.columns]]
        Cl2 = Cl2[mdl.index]
        Zl, lgl = Zl[mdl.index], lgl[mdl.index]
        mzl, altl = module_scores(lgl, Zl, Cl2, M)
        depth = Cl2.sum(0)
        for tp in tps:
            sel = mdl.index[mdl.timepoint == tp]
            for scoring, S in (("meanz", mzl), ("alt", altl)):
                for m in S.columns:
                    b, p = arm_model(S.loc[sel, m], mdl.loc[sel], depth[sel])
                    rep[(scoring, m, f"{tag}_{tp}")] = (b, p)
            print(f"  {tag} {tp}: n={len(sel)}, LAIV={(mdl.loc[sel].arm=='LAIV').sum()}")

    # ---- ENGAGED determination
    rows = []
    for scoring in ("meanz", "alt"):
        d = disc[scoring].set_index("module")
        for m in d.index:
            contrasts = {k[2]: v for k, v in rep.items() if k[0] == scoring and k[1] == m}
            agree = {c: (b, p) for c, (b, p) in contrasts.items()
                     if np.isfinite(p) and p < 0.05 and np.sign(b) == np.sign(d.delta[m])}
            rows.append(dict(scoring=scoring, module=m, delta=d.delta[m], frac_up=d.frac_up[m],
                             p_disc=d.p[m], bh_disc=d.p_adj_BH[m], n_rep=len(agree),
                             rep_contrasts=",".join(sorted(agree)),
                             passes=bool(d.p_adj_BH[m] < BH_MAX and agree)))
    V = pd.DataFrame(rows)
    piv = V.pivot_table(index="module", columns="scoring", values="passes", aggfunc="first")
    piv = piv.reindex(columns=["meanz", "alt"]).fillna(False)
    engaged = sorted(piv.index[piv.meanz & piv.alt])
    scoring_dep = sorted(piv.index[(piv.meanz | piv.alt) & ~(piv.meanz & piv.alt)])
    print(f"\n=== ENGAGED (both scorings): {engaged}")
    print(f"=== scoring_dependent: {scoring_dep}")
    V.to_csv(out / "module_confirmation.csv", index=False)

    if not engaged:
        json.dump({"engaged": [], "scoring_dependent": scoring_dep, "status": "FAILED"},
                  open(out / "confirmation_summary.json", "w"), indent=1)
        print("\nNO MODULE ENGAGED -> failure per plan §6")
        return

    # ---- gene refinement within engaged modules
    universe = sorted({g for m in engaged for g in M[m]["genes"] if g in lg.index})
    bi = [md.index[(md.participant == p) & (md.tp == "baseline")][0] for p in pairs]
    si = [md.index[(md.participant == p) & (md.tp == "SD3")][0] for p in pairs]
    B, S_ = lg.loc[universe, bi].to_numpy(), lg.loc[universe, si].to_numpy()
    t, p = stats.ttest_rel(S_, B, axis=1)
    G = pd.DataFrame({"gene": universe, "log2FC": (S_ - B).mean(1), "p": p})
    G["p_adj_BH"] = multipletests(np.nan_to_num(G.p, nan=1.0), method="fdr_bh")[1]
    G = G.set_index("gene")
    print(f"\ngene refinement: universe {len(G)} genes from {len(engaged)} engaged modules; "
          f"BH<0.05: {(G.p_adj_BH < BH_MAX).sum()}")

    for tag, (Cl, mdl, tps) in liverpool().items():
        lgl, _, Cl2 = prep(Cl)
        mdl = mdl.loc[[s for s in mdl.index if s in Cl2.columns]]
        Cl2, lgl = Cl2[mdl.index], lgl[mdl.index]
        depth = Cl2.sum(0)
        for tp in tps:
            sel = mdl.index[mdl.timepoint == tp]
            bs, ps = [], []
            for g in G.index:
                if g in lgl.index:
                    b, pp = arm_model(lgl.loc[g, sel], mdl.loc[sel], depth[sel])
                else:
                    b, pp = np.nan, np.nan
                bs.append(b); ps.append(pp)
            G[f"{tag}_{tp}_beta"], G[f"{tag}_{tp}_p"] = bs, ps
    labs = [c[:-5] for c in G.columns if c.endswith("_beta")]
    agree = pd.DataFrame({l: (np.sign(G[f"{l}_beta"]) == np.sign(G.log2FC))
                             & (G[f"{l}_p"] < 0.05) for l in labs})
    G["n_rep"] = agree.sum(1)
    G["core"] = (G.p_adj_BH < BH_MAX) & (G.n_rep >= 1)
    G["module"] = [",".join(m for m in engaged if g in M[m]["genes"]) for g in G.index]
    G.to_csv(out / "gene_refinement.csv")

    cores = {}
    for m in engaged:
        gs = [g for g in M[m]["genes"] if g in G.index and G.core[g]]
        frac = len(gs) / len([g for g in M[m]["genes"] if g in G.index])
        print(f"  {m}: {len(gs)} core / {len([g for g in M[m]['genes'] if g in G.index])} "
              f"measured ({frac:.0%})" + ("  -> _core variant" if len(gs) >= MIN_CORE else
                                          "  -> too few, no variant"))
        if len(gs) >= MIN_CORE:
            cores[f"{m}_core"] = gs

    # ---- composition check
    import anndata as ad
    a = ad.read_h5ad(ZIEGLER)
    a.var_names = a.var["feature_name"].astype(str).values
    a.var_names_make_unique()
    ct = a.obs["Coarse_Cell_Annotations"].astype(str).values
    comp = {}
    for name, gs in list(cores.items()) + [(m, [g for g in M[m]["genes"] if g in a.var_names])
                                           for m in engaged]:
        gs = [g for g in gs if g in a.var_names]
        if not gs:
            continue
        X = a[:, gs].X
        X = np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X)
        prof = pd.DataFrame(X, columns=gs).groupby(ct).mean()
        prof = prof.div(prof.sum(0).replace(0, np.nan), axis=1)
        top = prof.idxmax(0).value_counts()
        frac = float(top.iloc[0] / len(gs))
        comp[name] = dict(n_genes=len(gs), dominant=str(top.index[0]), frac=frac,
                          flagged=bool(frac > 0.5))
        print(f"  composition {name}: {top.index[0]} for {frac:.0%} of {len(gs)} genes"
              + ("  [FLAGGED]" if frac > 0.5 else ""))

    json.dump({"engaged": engaged, "scoring_dependent": scoring_dep,
               "core_modules": {k: v for k, v in cores.items()}, "composition": comp},
              open(out / "confirmation_summary.json", "w"), indent=1)
    print("\nwrote", out)


if __name__ == "__main__":
    main()
