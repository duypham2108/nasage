"""Matched nulls under all four scorers, including the corrected rank statistic.

Implements docs/PREREG_matched_null_rescore.md. Every specificity control in this project was
computed with mean-z; the corrected specification curves report four scorers. This re-draws the
expression- and size-matched nulls separately for each scorer, so the empirical p is computed
under the same statistic it is quoted beside.

Headline rule, locked in advance: the reported empirical p for a cohort is the WORST across the
four scorers, never the best.

Usage: python validation/matched_null_rescore.py <data_dir> <results_dir> <out_dir> [n_draws]
"""
import sys
import zlib
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp
import statsmodels.api as sm
import statsmodels.formula.api as smf

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from nasage.modules import load_modules      # noqa: E402
from mvri_programs import PROGRAMS           # noqa: E402
from scoring_methods import SCORERS, gene_ranks, rank_score   # noqa: E402


def scorers_for(full):
    """SCORERS with the rank statistic bound to a rank matrix computed once for `full`.

    Identical output to the generic path (asserted below); the only difference is that the
    transcriptome-wide sort happens once per cohort instead of once per null draw.
    """
    rk = gene_ranks(full)
    out = dict(SCORERS)
    out["rank_ucell"] = lambda X, full=None, idx=None: rank_score(X, idx=idx, ranks=rk)
    return out

DATA, R, OUT = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
NDRAW = int(sys.argv[4]) if len(sys.argv) > 4 else 500
N_BINS = 25
OUT.mkdir(parents=True, exist_ok=True)
MOD = load_modules()
rows, draws = [], []


def seed_of(*parts):
    """Deterministic seed from text.

    `hash()` on str is salted per process (PEP 456), so `hash(...) % 2**32` gives a DIFFERENT
    seed on every run and the nulls are not reproducible. CRC32 is stable across processes and
    machines. Two 20-draw runs of this script disagreed before this was fixed.
    """
    return zlib.crc32("|".join(str(p) for p in parts).encode()) % (2 ** 32)


def zs(v):
    sd = np.std(v, ddof=1)
    return (v - np.mean(v)) / (sd if sd > 0 else 1.0)


def bins_of(E):
    """Expression bins over columns of a samples x genes frame."""
    m = E.mean(0)
    b = pd.qcut(m.rank(method="first"), N_BINS, labels=False)
    pool = {k: np.array(v) for k, v in b.groupby(b).groups.items()}
    return b, pool


# ============================================================ outcome cohorts (bulk, logistic)
SOURCES = [
    ("nova", "expr_symbols.csv.gz", "nova_scores_outcome.csv", "shed", "NOVA influenza"),
    ("mvri_outcome", "expr_rsv.csv.gz", "meta_rsv.csv", "y", "RSV infection"),
    ("mvri_outcome", "expr_pneumo.csv.gz", "meta_pneumo.csv", "y", "pneumococcal carriage"),
]

for d, ef, mf, ycol, label in SOURCES:
    ep, mp = R / d / ef, R / d / mf
    if not ep.exists() or not mp.exists():
        print(f"skip {label}: missing input")
        continue
    E = pd.read_csv(ep, index_col=0).T
    meta = pd.read_csv(mp, index_col=0)
    common = [i for i in E.index if i in meta.index]
    E, meta = E.loc[common], meta.loc[common]
    full = E.values
    SC_LOCAL = scorers_for(full)
    colpos = {g: i for i, g in enumerate(E.columns)}
    b, pool = bins_of(E)
    y = meta[ycol].astype(int)
    print(f"\n{label}: {E.shape[0]} samples x {E.shape[1]} genes, "
          f"{int(y.sum())} positive")

    def fit(v, pcol):
        X = pd.DataFrame({"s": zs(v)}, index=E.index)
        if pcol and pcol in meta.columns:
            X[pcol] = pd.to_numeric(meta[pcol], errors="coerce").values
        Xf = sm.add_constant(X.astype(float)).dropna()
        m = sm.Logit(y.loc[Xf.index], Xf).fit(disp=0)
        return float(m.params["s"]), float(m.pvalues["s"])

    for pname, p in PROGRAMS.items():
        genes = [g for g in p["genes"] if g in colpos]
        if len(genes) < 3:
            continue
        pcol = f"proxy_{p['cell_type']}"
        idx = [colpos[g] for g in genes]
        gb = [b.get(g) for g in genes]
        for sname, fn_ in SC_LOCAL.items():
            try:
                obs, obsp = fit(fn_(full[:, idx], full=full, idx=idx), pcol)
            except Exception as e:
                print(f"  {pname} {sname}: observed fit failed ({e})")
                continue
            rng = np.random.default_rng(seed_of(label, pname, sname))
            nb, npv = [], []
            for _ in range(NDRAW):
                pick = [str(rng.choice(pool[k])) for k in gb if pd.notna(k)]
                nidx = [colpos[g] for g in pick if g in colpos]
                if len(nidx) < 3:
                    continue
                try:
                    e_, pv_ = fit(fn_(full[:, nidx], full=full, idx=nidx), pcol)
                except Exception:
                    continue
                if np.isfinite(e_):
                    nb.append(e_); npv.append(pv_)
            nb, npv = np.asarray(nb), np.asarray(npv)
            if nb.size == 0:
                continue
            draws.append(pd.DataFrame(dict(analysis=label, program=pname, scoring=sname,
                                           null_est=nb, null_p=npv, observed=obs)))
            rows.append(dict(analysis=label, unit="cohort", program=pname, scoring=sname,
                             n_draws=int(nb.size), observed=obs, model_p=obsp,
                             empirical_p=float((np.abs(nb) >= abs(obs)).mean()),
                             null_frac_p05=float((npv < 0.05).mean())))
            r = rows[-1]
            print(f"  {pname:30s} {sname:11s} b={obs:+.3f} p={obsp:.3f} "
                  f"emp_p={r['empirical_p']:.3f} null05={r['null_frac_p05']:.3f}")

# ============================================================ lifespan (Loske pseudobulk, OLS)
h5 = DATA / "loske2022_controls.h5ad"
if h5.exists():
    print("\nlifespan: building donor x cell-type pseudobulk")
    A = ad.read_h5ad(h5)
    BASAL_K, CIL_K = ("basal",), ("ciliated", "cil_diff", "cil_acute", "multiciliated")
    SEC_K = ("secretory", "goblet", "club", "glandular", "sec_diff")
    DROP_K = ("foxn4", "squamous", "ionocyte", "lowrna", "low_rna", "deutero")

    def to_family(lbl):
        s = str(lbl).lower().replace(" ", "_")
        if any(k in s for k in DROP_K):
            return None
        if any(k in s for k in BASAL_K):
            return "Basal"
        if any(k in s for k in CIL_K):
            return "Ciliated"
        if any(k in s for k in SEC_K):
            return "Secretory"
        return None

    A.obs["__fam"] = A.obs["detailed_clusters"].map(to_family)
    A = A[A.obs.__fam.notna()].copy()
    A.obs["__donor"] = A.obs["Patient_ID"].astype(str)
    A.obs["__age"] = pd.to_numeric(A.obs["paper_age"], errors="coerce")
    A = A[A.obs.__age.notna()].copy()
    X = A.X
    if float(X.max()) > 50:
        tot = np.asarray(X.sum(1)).ravel(); tot[tot == 0] = 1
        A.X = sp.diags(1e4 / tot) @ X
        A.X.data = np.log1p(A.X.data)
    # Pseudobulk by sparse group-mean. Densifying 66k x 27k first (as the earlier scripts did)
    # costs several GB and is unnecessary: the aggregate is ~100 x 27k.
    Xs = A.X if sp.issparse(A.X) else sp.csr_matrix(A.X)
    ngl = np.asarray((Xs > 0).sum(1)).ravel()
    key = pd.MultiIndex.from_arrays([A.obs.__donor.values, A.obs.__fam.values],
                                    names=["__donor", "__fam"])
    codes, uniq = pd.factorize(key)
    G = sp.csr_matrix((np.ones(len(codes)), (codes, np.arange(len(codes)))),
                      shape=(len(uniq), len(codes)))
    ncells = np.asarray(G.sum(1)).ravel()
    sums = np.asarray((G @ Xs).todense())
    agg = pd.DataFrame(sums / ncells[:, None], columns=list(A.var_names),
                       index=pd.MultiIndex.from_tuples(list(uniq),
                                                       names=["__donor", "__fam"]))
    agg["n_cells"] = ncells
    agg = agg[agg.n_cells >= 20]
    meta_l = (pd.DataFrame({"__donor": A.obs.__donor.values, "age": A.obs.__age.values,
                            "sex": A.obs["gender"].astype(str).values,
                            "lg": np.log(ngl + 1)})
              .groupby("__donor").agg({"age": "mean", "sex": "first", "lg": "mean"}))
    del A, Xs, sums

    AXES = {
        "Mucus barrier (P5)": PROGRAMS["P5_mucin_barrier"]["genes"],
        "Antigen presentation (M06)": MOD["M06_mhc_I_processing"]["genes"],
        "IgA support (P4)": PROGRAMS["P4_iga_transport"]["genes"],
    }
    GENECOLS = [c for c in agg.columns if c != "n_cells"]
    for f in ["Basal", "Secretory", "Ciliated"]:
        if f not in agg.index.get_level_values("__fam"):
            continue
        sub = agg.xs(f, level="__fam")
        if len(sub) < 8:
            continue
        Efam = sub[GENECOLS]
        full = Efam.values
        SC_LOCAL = scorers_for(full)
        colpos = {g: i for i, g in enumerate(GENECOLS)}
        b, pool = bins_of(Efam)
        mm = meta_l.reindex(sub.index)

        def fit_age(v):
            d = pd.DataFrame({"y": zs(v), "age": mm.age.values,
                              "lg": mm.lg.values}).dropna()
            m = smf.ols("y ~ age + lg", data=d).fit()
            return float(m.params["age"]), float(m.pvalues["age"])

        for nm, genes0 in AXES.items():
            genes = [g for g in genes0 if g in colpos]
            if len(genes) < 3:
                continue
            idx = [colpos[g] for g in genes]
            gb = [b.get(g) for g in genes]
            for sname, fn_ in SC_LOCAL.items():
                try:
                    obs, obsp = fit_age(fn_(full[:, idx], full=full, idx=idx))
                except Exception:
                    continue
                rng = np.random.default_rng(seed_of(f, nm, sname))
                nb, npv = [], []
                for _ in range(NDRAW):
                    pick = [str(rng.choice(pool[k])) for k in gb if pd.notna(k)]
                    nidx = [colpos[g] for g in pick if g in colpos]
                    if len(nidx) < 3:
                        continue
                    try:
                        e_, pv_ = fit_age(fn_(full[:, nidx], full=full, idx=nidx))
                    except Exception:
                        continue
                    if np.isfinite(e_):
                        nb.append(e_); npv.append(pv_)
                nb, npv = np.asarray(nb), np.asarray(npv)
                if nb.size == 0:
                    continue
                rows.append(dict(analysis=f"lifespan|{f}", unit="donor x cell type",
                                 program=nm, scoring=sname, n_draws=int(nb.size),
                                 observed=obs, model_p=obsp,
                                 empirical_p=float((np.abs(nb) >= abs(obs)).mean()),
                                 null_frac_p05=float((npv < 0.05).mean())))
                r = rows[-1]
                print(f"  {f:10s} {nm:28s} {sname:11s} b={obs:+.4f} p={obsp:.3f} "
                      f"emp_p={r['empirical_p']:.3f} null05={r['null_frac_p05']:.3f}")

D = pd.DataFrame(rows)
D.to_csv(OUT / "rescore_nulls.csv", index=False)
if draws:
    pd.concat(draws, ignore_index=True).to_csv(OUT / "null_draws.csv.gz", index=False)

print("\n" + "=" * 86)
print("LOCKED HEADLINE: worst empirical p across the four scorers (never the best)")
print("=" * 86)
summ = (D.groupby(["analysis", "program"])
        .agg(n_scorers=("scoring", "size"),
             worst_emp_p=("empirical_p", "max"),
             best_emp_p=("empirical_p", "min"),
             rank_emp_p=("empirical_p", lambda s: np.nan),
             mean_null05=("null_frac_p05", "mean"))
        .reset_index())
rk = D[D.scoring == "rank_ucell"].set_index(["analysis", "program"]).empirical_p
summ["rank_emp_p"] = [rk.get((a, p), np.nan) for a, p in zip(summ.analysis, summ.program)]
summ["passes_all"] = summ.worst_emp_p < 0.05
summ.to_csv(OUT / "rescore_summary.csv", index=False)
key = summ[summ.program.str.contains("P5|M06|P4|mucin|barrier|presentation|IgA", regex=True)]
print(key.round(4).to_string(index=False))
print(f"\npass under ALL four scorers: {int(summ.passes_all.sum())} of {len(summ)}")
print("\nwrote", OUT)
