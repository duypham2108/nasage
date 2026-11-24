"""Tests: the scanpy-style API end to end (on simulated counts), plus unit tests."""
import numpy as np
import pandas as pd
import pytest

import nasage as nas
from nasage import _design, aggregate, qc, scoring, stats
from nasage import modules as mods
from nasage.simulate import simulate

KEYS = dict(donor_key="donor", celltype_key="cell_type", age_key="age_bin",
            batch_key="batch", sex_key="sex")


@pytest.fixture(scope="module")
def adata():
    M, A = nas.load_modules(), nas.load_applicability()
    a, _ = simulate("planted", {k: v["genes"] for k, v in M.items()}, A, seed=1, n_per_bin=4,
                    cells_per_donor=250, n_genes=1200)
    a.layers["counts"] = a.X.copy()
    # a finer annotation: a third of secretory cells are "Club"
    a.obs["fine"] = np.where((a.obs.cell_type == "secretory").values
                             & (np.arange(a.n_obs) % 3 == 0), "Club",
                             a.obs.cell_type.astype(str))
    with pytest.warns(UserWarning, match="mitochondrial"):
        nas.pp.filter_cells(a, celltype_key="cell_type", layer="counts", max_pct_mito=None)
    nas.tl.score_modules(a, layer="counts")
    nas.tl.donor_summary(a, subtype_key="fine", **KEYS)
    nas.tl.age_effects(a, reference="18-39")
    nas.tl.composition(a)
    return a


def _st(a):
    return a.uns["nasage"]


# ------------------------------------------------------------------- pipeline (API)
def test_modules_load_and_map(adata):
    M, A = nas.load_modules(), nas.load_applicability()
    assert len(M) == 15 and set(M) == set(A)
    st = _st(adata)
    assert st["module_presence"].frac_found.min() > 0.9
    cols = st["score_columns"]
    assert any(c.endswith("_UCell") for c in cols)
    assert "M03_pigr_antibody_transport_logexpr" in cols       # single gene, not a signature
    assert "M10_iga_fraction_logit" in cols                     # ratio feature


def test_scores_are_cell_type_specific(adata):
    med = adata.obs.groupby("cell_type", observed=True)[
        "M09_plasma_secretory_machinery_UCell"].median()
    assert med.idxmax() == "plasma"


def test_primary_family_is_epithelial_by_default(adata):
    om = _st(adata)["omnibus"]
    prim = set(om.loc[om.tier == "primary", "cell_type"].astype(str))
    assert prim <= set(nas.tl.PRIMARY_CELLTYPES)
    # immune pairs that the matrix marks primary are still tested, as exploratory
    ex = om[(om.module_id == "M05_interferon_stimulated_genes") & (om.cell_type == "Tcell")]
    assert len(ex) == 1 and ex.tier.iloc[0] == "exploratory"
    # the matrix itself is unchanged
    assert nas.load_applicability()["M05_interferon_stimulated_genes"]["Tcell"] == 1.0


def test_primary_celltypes_none_keeps_matrix(adata):
    a = adata.copy()
    nas.tl.donor_summary(a, subtype_key="fine", primary_celltypes=None, **KEYS)
    nas.tl.age_effects(a, reference="18-39")
    om = _st(a)["omnibus"]
    assert "Tcell" in set(om.loc[om.tier == "primary", "cell_type"].astype(str))
    with pytest.raises(ValueError, match="not in the applicability matrix"):
        nas.tl.donor_summary(a, primary_celltypes=("epithelium",), **KEYS)


def test_missing_is_missing_not_zero(adata):
    s = _st(adata)["donor_summary"]
    inel = s[~s.eligible]
    assert len(inel) > 0 and inel[_st(adata)["score_columns"]].isna().all().all()


def test_age_effects_run_and_adjust(adata):
    om, ct = _st(adata)["omnibus"], _st(adata)["contrasts"]
    assert om.estimable.any()
    sel = om.p_adj_BH.notna()
    assert (om.loc[sel, "p_adj_BH"] >= om.loc[sel, "p_omnibus"] - 1e-9).all()
    assert {"estimate", "ci_low", "ci_high", "n_donors"} <= set(ct.columns)
    assert om.covariates.str.contains("log_genes").any()        # depth covariate used


def test_composition_is_compositional(adata):
    assert np.allclose(_st(adata)["proportions"].sum(1), 1.0)
    assert {"CLR", "asin_sqrt"} == set(_st(adata)["composition"]["transform"])


def test_iga_ratio_finite(adata):
    s = _st(adata)["donor_summary"]
    v = s.loc[(s.cell_type == "plasma") & s.eligible, "M10_iga_fraction_logit"]
    assert v.notna().all() and np.isfinite(v).all()


def test_subtype_fraction_and_correlation(adata):
    st = _st(adata)
    assert st["subtype_covariates"]["secretory"] == ["Club"]
    s = st["donor_summary"]
    assert s[(s.cell_type == "secretory") & s.eligible].subfrac_1.between(0.2, 0.5).all()
    r = nas.tl.module_correlation(adata, "secretory")
    assert r.shape[0] == r.shape[1] and np.allclose(np.diag(r), 1)


def test_permutation_and_report(adata, tmp_path):
    perm = nas.tl.permutation_null(adata, n_perm=3)
    assert len(perm) == 3 and "per_test_false_positive_rate" in _st(adata)["permutation_rates"]
    path = nas.write_report(adata, tmp_path, label="test")
    assert path.exists()
    for f in ("age_effects_omnibus.csv", "donor_celltype_scores.csv", "qc_report.csv"):
        assert (tmp_path / f).exists()
    for fig in (nas.pl.effect_heatmap(adata), nas.pl.composition(adata)):
        assert fig is not None


def test_results_survive_h5ad_roundtrip(adata, tmp_path):
    import anndata as ad
    adata.write_h5ad(tmp_path / "x.h5ad")
    back = ad.read_h5ad(tmp_path / "x.h5ad")
    assert "omnibus" in back.uns["nasage"]


def test_qc_report_and_counts_check(adata):
    rep = _st(adata)["qc"]
    assert {"cell_type", "frac_kept", "max_pct_mito_used"} <= set(rep.columns)
    b = adata[:50].copy()
    b.layers["lognorm"] = np.log1p(b.layers["counts"].toarray() / 3.0)
    with pytest.raises(ValueError, match="integer counts"):
        nas.tl.score_modules(b, layer="lognorm")


# --------------------------------------------------------------------- unit tests
def test_packaged_resources_resolve():
    for name in ("gene_modules.json", "module_celltype_applicability.json",
                 "gene_set_provenance_modules.csv", "module_evidence.csv"):
        assert mods._resource(name).exists()


def test_column_helper_sparse_and_dense():
    import scipy.sparse as sp
    M = np.array([[0, 2], [3, 0], [0, 5]], dtype=float)
    for X in (M, sp.csr_matrix(M), sp.csc_matrix(M)):
        assert np.array_equal(scoring._column(X, 1), [2, 0, 5])


def _toy_summary(batch_equals_age: bool):
    rows = []
    for b, ab in enumerate(["5-11", "18-39"]):
        for k in range(4):
            rows.append(dict(donor=f"{ab}_{k}", cell_type="basal", age_bin=ab,
                             sex="F" if k < 2 else "M",
                             batch=f"B{b}" if batch_equals_age else f"B{k % 2}",
                             M_UCell=float(k) + b))
    return pd.DataFrame(rows)


def test_rank_deficient_design_not_estimable():
    m, info = stats.fit_one(_toy_summary(True), "M_UCell", "18-39")
    assert m is None and "not full rank" in info["reason"]
    m, _ = stats.fit_one(_toy_summary(False), "M_UCell", "18-39")
    assert m is not None


def test_donor_level_variables_must_be_constant():
    obs = pd.DataFrame(dict(donor=["d1", "d1", "d2"], batch=["B1", "B2", "B1"]))
    with pytest.raises(ValueError, match="d1"):
        _design._check_donor_constant(obs, ["batch"])
    _design._check_donor_constant(obs.assign(batch="B1"), ["batch"])


def test_permutation_is_within_batch(monkeypatch):
    s = _toy_summary(False)
    seen = []
    real = stats.fit_age_models

    def spy(summary, *a, **k):
        seen.append(pd.crosstab(summary.batch, summary.age_bin))
        return real(summary, *a, **k)

    monkeypatch.setattr(stats, "fit_age_models", spy)
    stats.permutation_null(s, [("M", "basal", "primary")], "18-39", n_perm=5)
    observed = pd.crosstab(s.batch, s.age_bin)
    assert all(t.equals(observed) for t in seen)


def test_age_groups_and_bins():
    import anndata as ad
    obs = pd.DataFrame(dict(pid=["a", "a", "b", "c"], ct=["basal"] * 4,
                            grp=["18-39", "18-39", "5-11", None], yrs=[30, 30, 7, 80]),
                       index=[f"c{i}" for i in range(4)])
    a = ad.AnnData(X=np.ones((4, 2), dtype=np.float32), obs=obs)
    p = dict(donor_key="pid", celltype_key="ct", age_key="grp", batch_key="pid", sex_key="pid")
    with pytest.warns(UserWarning, match="no age group"):
        out = _design.canonical_obs(a, p)
    assert list(out.age_bin.cat.categories) == ["5-11", "18-39"] and len(out) == 3
    with pytest.raises(ValueError, match="not listed"):
        _design.canonical_obs(a, dict(p, age_order=["18-39"]))
    with pytest.warns(UserWarning, match="outside age_bins"):
        out = _design.canonical_obs(a, dict(p, age_key="yrs",
                                            age_bins=[(5, 11, "5-11"), (18, 39, "18-39")]))
    assert set(out.donor) == {"a", "b"}


def test_module_content_changes():
    M = nas.load_modules()
    for g in ("RIGI", "IFIH1", "DHX58", "TLR3", "AIM2", "CGAS"):  # ISG sensors out of M04
        assert g not in M["M04_innate_sensing"]["genes"]
    m9 = M["M09_plasma_secretory_machinery"]["genes"]
    assert not any(g.startswith("IGH") or g == "JCHAIN" for g in m9)
    no_hla = nas.drop_hla(M)
    assert not any(g.startswith("HLA-") for m in no_hla.values() for g in m["genes"])
    assert "B2M" in no_hla["M06_mhc_I_processing"]["genes"]


def test_auto_statistic_uses_mean_for_sparse():
    obs = pd.DataFrame(dict(cell_type=pd.Categorical(["a"] * 4 + ["b"] * 4),
                            dense_UCell=[.2, .3, .4, .5] * 2,
                            sparse_UCell=[0, 0, 0, .4, .1, .2, .3, .4]))
    plan = aggregate.choose_statistic(obs, ["dense_UCell", "sparse_UCell"])
    got = {(r.score_column, r.cell_type): r.statistic for r in plan.itertuples()}
    assert got[("sparse_UCell", "a")] == "mean" and got[("sparse_UCell", "b")] == "median"
    assert got[("dense_UCell", "a")] == "median"


def test_adaptive_mito_keeps_high_mito_cell_type():
    rng = np.random.default_rng(0)
    n = 200
    obs = pd.DataFrame(dict(n_genes_by_counts=1000, total_counts=3000,
                            pct_counts_mt=np.r_[rng.normal(25, 2, n), rng.normal(5, 1, n)]),
                       index=[f"c{i}" for i in range(2 * n)])
    ct = pd.Series(["ciliated"] * n + ["basal"] * n, index=obs.index)
    keep, _ = qc.filter_mask(obs, ct, max_pct_mito=20, mito_adaptive=False)
    assert keep[ct == "ciliated"].mean() < 0.05                  # fixed cut-off wipes them out
    keep, rep = qc.filter_mask(obs, ct, max_pct_mito=20, mito_adaptive=True)
    rep = rep.set_index("cell_type")
    assert rep.loc["ciliated", "frac_kept"] > 0.9 and rep.loc["basal", "max_pct_mito_used"] == 20
    keep, _ = qc.filter_mask(obs, ct, min_genes={"default": 2000, "ciliated": 500},
                             max_pct_mito=None)
    assert set(ct[keep]) == {"ciliated"}


def test_subtype_and_depth_covariates_enter_model():
    rng = np.random.default_rng(1)
    rows = []
    for ab in ["5-11", "18-39"]:
        for k in range(8):
            rows.append(dict(donor=f"{ab}_{k}", cell_type="secretory", age_bin=ab,
                             sex="F" if k < 4 else "M", batch=f"B{k % 2}",
                             log_genes=rng.normal(7, .2), subfrac_1=rng.uniform(.1, .6),
                             subfrac_2=np.nan, M_UCell=rng.normal()))
    s = pd.DataFrame(rows)
    m, info = stats.fit_one(s, "M_UCell", "18-39")
    assert m is not None
    assert "log_genes" in info["covariates"] and "subfrac_1" in info["covariates"]
    assert "subfrac_2" not in info["covariates"]                 # all-NaN -> unused
    # 8 donors, balanced sex/batch: 6 parameters would leave 2 residual df < 3
    m, info = stats.fit_one(s.iloc[[0, 1, 4, 5, 8, 9, 12, 13]], "M_UCell", "18-39")
    assert m is not None and info["dropped_covariates"] == "subfrac_1"
    assert m.df_resid >= stats.MIN_RESID_DF


def test_summary_attrs_do_not_break_pandas(adata):
    s = _st(adata)["donor_summary"]
    repr(pd.concat([s, s], ignore_index=True))
    assert len(_st(adata)["aggregation"]) > 0


def test_composition_ignores_empty_age_groups(adata):
    obs = nas.tl._obs(adata)
    obs["age_bin"] = obs.age_bin.cat.add_categories(["65+"])       # defined, but no donors
    res, _ = stats.composition_analysis(obs, "18-39")
    assert res.p_omnibus.notna().all()


def test_null_control(adata):
    res = nas.tl.null_control(adata, n_sets=8, seed=1)
    assert {"null_frac_p05", "empirical_p", "generic_age_signal"} <= set(res.columns)
    assert (res.n_null > 0).all() and res.empirical_p.between(0, 1).all()
    assert res.empirical_p.min() >= 1 / 9 - 1e-9                # floor is 1 / (n_sets + 1)
    # no generic age signal is planted in the simulator: random sets rarely reach p < 0.05
    assert res.null_frac_p05.mean() < 0.15


def test_matched_null_sets_match_expression():
    from nasage import controls
    bins = pd.Series([0, 0, 1, 1, 2, 2, 2], index=list("abcdefg"), dtype=float)
    sets = controls.matched_null_sets({"M": ["a", "c", "e"]}, bins, n_sets=5, seed=0)
    for genes in sets.values():
        assert sorted(bins[genes]) == [0, 1, 2]                 # one gene per partner bin
        assert not {"a", "c", "e"} & set(genes)                 # module genes excluded
