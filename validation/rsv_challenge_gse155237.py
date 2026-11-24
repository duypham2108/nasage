"""PRE-SPECIFIED analysis: NasAge axes vs RSV challenge outcome (GSE155237).

================================ SPECIFICATION ================================
This docstring is the pre-specification. It is committed BEFORE the expression
data are downloaded or inspected; the commit timestamp is the registration. It
extends docs/PREREG_vaccine_axes.md (locked in commit 1988050) to a second
dataset. The predicted SIGNS are not new: they are taken unchanged from section 1
of that file, which was locked before any data were seen.

--- What had been seen when this was written (2026-11-24) -----------------------
SEEN: GEO sample metadata only (donor_id, batch, timepoint, outcome, symptom) and
the series summary; the PubMed abstract of the source paper (Habibi et al. 2020,
Science, doi:10.1126/science.aba9301). NOT SEEN: any expression value from this
series, and any relation between expression and outcome.

--- Rationale for the dataset ---------------------------------------------------
Our locked table predicts that baseline interferon tone OPPOSES replication of a
live, replicating respiratory virus. LAIV outcome data cannot be linked to the
public transcriptomes, so this tests the same mechanistic arm with a challenge
virus instead of a live vaccine. It is therefore a PROXY for the live-replicating
arm, never a vaccine-response result, and must be reported as such.

--- Design (from GEO metadata) --------------------------------------------------
40 donors, nasal curettage, timepoints D0 (pre-challenge) and D3, RSV challenge.
Outcome per donor: Infected (24) / Uninfected (16); at D0, 23 vs 16.
Only D0 (pre-challenge) samples enter the primary analysis.

--- Independence from the source publication ------------------------------------
The authors ran a genome-wide pre-challenge analysis (80 DEGs vs susceptibility).
Their REPORTED baseline predictor is mucosal neutrophil activation, with IL-17 and
TNF pathways in those who resisted infection; interferon is not their finding.
Therefore:
  * H1/H4 below are NOT restatements of their result, but they are also not fully
    independent, because their screen covered all genes. Verdicts are labelled
    "pre-specified, source-paper screen covered the same contrast".
  * M12 overlaps their reported finding directly (it contains CXCL1, IL6, IL1A/B,
    TNF). It is therefore NOT a test: it is a POSITIVE CONTROL. If NasAge scoring
    on this dataset cannot recover their published neutrophil-inflammation signal,
    the other results here are uninterpretable.

--- Hypotheses (signs fixed in commit 1988050; one-sided) -----------------------
Endpoint E: Infected vs Uninfected (binary, per donor), D0 samples.
  H1  M05 interferon-stimulated genes -> OR per SD < 1  (higher ISG tone resists)
  H4  M02b airway mucins             -> OR per SD < 1
Primary family = {H1, H4}; Holm across 2 tests, one-sided in the predicted
direction.

Secondary (reported, never promoted):
  M06 MHC-I  -> OR < 1 (not independent of M05 if r(M05,M06) > 0.7)
  M04 innate sensing -> OR < 1
  the same tests using `symptom` (cold / no_cold) as the endpoint
  M01, M02, M03, M11 -> two-sided, no prediction, descriptive only

Positive control (NOT a hypothesis, direction follows the source paper):
  M12 inflammatory SASP -> OR per SD > 1

--- Scoring (as docs/PREREG_vaccine_axes.md section 4) --------------------------
Module score = mean over module genes of z-scored log2(CPM+1); z computed across
the analysed (D0) samples; genes with mean CPM < 1 across D0 dropped. Samples with
< 5 M assigned counts excluded. Sensitivity: UCell on counts.

--- Model -----------------------------------------------------------------------
Logistic regression, one row per donor:
    outcome ~ module_score_z + batch + log10_counts + PTPRC_log2cpm
Age and sex are NOT in this series' metadata and so cannot be adjusted for; this
is recorded as a limitation, not a deviation.
Sensitivity: (a) no PTPRC term; (b) no batch term; (c) new-batch donors only;
(d) UCell scores.

--- Module-specificity check ----------------------------------------------------
For every tested module, 1000 expression-matched random gene sets (25 expression
bins, as in nasage.controls) are scored and modelled identically.
empirical_p = (1 + #{random sets at least as extreme IN THE PREDICTED DIRECTION})
              / 1001.

--- Decision rules (as section 6 of the locked file) ----------------------------
  supported    : predicted sign, Holm p < 0.05, AND empirical p < 0.05
  consistent   : predicted sign, one-sided nominal p < 0.05, fails Holm or empirical
  not supported: anything else in the predicted direction, or null
  contradicted : opposite sign, two-sided p < 0.05
Power stated in advance: with 23 vs 16 donors, ~80% power needs OR per SD of about
2.6. Moderate true effects will be missed; "not supported" is not evidence of no
effect.

--- Deviations ------------------------------------------------------------------
Any departure from the above is appended, dated, to the results file. Signs and
decision rules are never changed.
================================ END SPECIFICATION ============================

Usage: python validation/rsv_challenge_gse155237.py <work_dir> <out_dir>
"""
from __future__ import annotations

import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

import nasage as nas

COUNTS_URL = ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE155nnn/GSE155237/suppl/"
              "GSE155237_gene_counts.csv.gz")
PRIMARY = {"M05_interferon_stimulated_genes": "neg", "M02b_airway_mucin_panel": "neg"}
SECONDARY = {"M06_mhc_I_processing": "neg", "M04_innate_sensing": "neg"}
POSCTRL = {"M12_inflammatory_sasp": "pos"}
DESCRIPTIVE = ["M01_epithelial_barrier", "M02_mucociliary_axoneme",
               "M03_pigr_antibody_transport", "M11_senescence_arrest"]
N_NULL = 1000


def load_counts(work: Path) -> pd.DataFrame:
    f = work / "GSE155237_gene_counts.csv.gz"
    if not f.exists():
        print("downloading counts")
        urllib.request.urlretrieve(COUNTS_URL, f)
    c = pd.read_csv(f, index_col=0)
    c = c.loc[:, [x for x in c.columns if pd.api.types.is_numeric_dtype(c[x])]]
    if c.index.astype(str).str.startswith("ENSG").mean() > 0.5:
        import anndata as ad
        var = ad.read_h5ad("D:/research/nasage_validation/data/ziegler2021_np.h5ad",
                           backed="r").var
        sym = var["feature_name"].astype(str).reindex(
            c.index.astype(str).str.split(".").str[0])
        c = c[sym.notna().values].groupby(sym[sym.notna()].values).sum()
    return c


def score_modules(cpm_log: pd.DataFrame, modules: dict) -> pd.DataFrame:
    mu, sd = cpm_log.mean(1), cpm_log.std(1, ddof=1)
    ok = sd > 0
    Z = cpm_log[ok].sub(mu[ok], axis=0).div(sd[ok], axis=0)
    return pd.DataFrame({m: Z.loc[[g for g in d if g in Z.index]].mean(0)
                         for m, d in modules.items()
                         if any(g in Z.index for g in d)})


def fit(y: pd.Series, score: pd.Series, cov: pd.DataFrame, direction: str):
    x = pd.concat([(score - score.mean()) / score.std(ddof=1), cov], axis=1)
    x = sm.add_constant(pd.get_dummies(x, drop_first=True).astype(float))
    try:
        m = sm.Logit(y.astype(float), x).fit(disp=0)
    except Exception as exc:
        return dict(reason=str(exc)[:70])
    b, p2 = float(m.params.iloc[1]), float(m.pvalues.iloc[1])
    p1 = p2 / 2 if ((b < 0) == (direction == "neg")) else 1 - p2 / 2
    return dict(log_or=b, or_per_sd=float(np.exp(b)), z=float(m.tvalues.iloc[1]),
                p_two_sided=p2, p_one_sided=p1,
                ci_low=float(np.exp(m.conf_int().iloc[1, 0])),
                ci_high=float(np.exp(m.conf_int().iloc[1, 1])))


def main():
    work, out = Path(sys.argv[1]), Path(sys.argv[2])
    work.mkdir(parents=True, exist_ok=True); out.mkdir(parents=True, exist_ok=True)
    meta = pd.read_csv(work / "g155_samples.csv")
    counts = load_counts(work)
    print("counts:", counts.shape, "| first columns:", list(counts.columns[:3]))

    # ---- match samples to metadata (GSM id or title)
    cols = pd.Index(counts.columns.astype(str))
    for key in ("title", "gsm"):
        hit = meta[key].astype(str)
        if cols.isin(hit).mean() > 0.5:
            meta = meta.set_index(hit.values).loc[cols[cols.isin(hit)]]
            counts = counts[cols[cols.isin(hit)]]
            break
    else:
        raise SystemExit(f"cannot match counts columns to metadata: {list(cols[:4])}")
    print("matched", counts.shape[1], "samples")

    tot = counts.sum(0)
    keep = (meta.timepoint == "D0").values & (tot >= 5e6).values
    print("D0 samples:", int((meta.timepoint == "D0").sum()),
          "| excluded <5M counts:", int(((meta.timepoint == "D0").values & ~(tot >= 5e6).values).sum()))
    md, C = meta[keep], counts.loc[:, keep]
    cpm = C / C.sum(0) * 1e6
    expressed = cpm.mean(1) >= 1
    logc = np.log2(cpm[expressed] + 1)
    print("genes kept:", int(expressed.sum()))
    print("outcome:", md.outcome.value_counts().to_dict())

    M = {m: d["genes"] for m, d in nas.load_modules().items()}
    S = score_modules(logc, M)
    cov = pd.DataFrame({"batch": md.batch.values,
                        "log10_counts": np.log10(C.sum(0).values),
                        "PTPRC": logc.loc["PTPRC"].values if "PTPRC" in logc.index else 0.0},
                       index=S.index)
    y = pd.Series((md.outcome.values == "Infected").astype(int), index=S.index)

    # ---- null gene sets, expression matched
    from nasage.controls import gene_expression_bins, matched_null_sets
    bins = gene_expression_bins(C.T.values, list(C.index))
    tested = {**PRIMARY, **SECONDARY, **POSCTRL}
    nulls = matched_null_sets({m: [g for g in M[m] if g in logc.index] for m in tested},
                              bins, n_sets=N_NULL, seed=0)
    SN = score_modules(logc, nulls)

    rows = []
    for m, d in {**tested, **{k: "two" for k in DESCRIPTIVE}}.items():
        if m not in S:
            rows.append(dict(module=m, note="not scoreable")); continue
        r = fit(y, S[m], cov, d)
        tier = ("primary" if m in PRIMARY else "secondary" if m in SECONDARY
                else "positive_control" if m in POSCTRL else "descriptive")
        r.update(module=m, tier=tier, predicted=d, n=len(y))
        if m in tested and "z" in r:
            zs = [fit(y, SN[c], cov, d).get("z", np.nan)
                  for c in SN.columns if c.startswith(m + "__null")]
            zs = np.array([v for v in zs if np.isfinite(v)])
            more = (zs <= r["z"]).sum() if d == "neg" else (zs >= r["z"]).sum()
            r["n_null"] = len(zs)
            r["empirical_p"] = float((1 + more) / (1 + len(zs)))
            r["null_frac_p05"] = float(np.mean(np.abs(zs) > 1.96))
        rows.append(r)
    R = pd.DataFrame(rows)

    prim = R.tier == "primary"
    from statsmodels.stats.multitest import multipletests
    R.loc[prim, "p_holm"] = multipletests(R.loc[prim, "p_one_sided"], method="holm")[1]

    def verdict(r):
        if r.tier not in ("primary",) or not np.isfinite(r.get("p_one_sided", np.nan)):
            return ""
        wrong = (r.log_or > 0) if r.predicted == "neg" else (r.log_or < 0)
        if wrong and r.p_two_sided < 0.05:
            return "contradicted"
        if not wrong and r.p_holm < 0.05 and r.empirical_p < 0.05:
            return "supported"
        if not wrong and r.p_one_sided < 0.05:
            return "consistent"
        return "not supported"
    R["verdict"] = R.apply(verdict, axis=1)

    corr = S[[c for c in ("M05_interferon_stimulated_genes", "M06_mhc_I_processing",
                          "M12_inflammatory_sasp") if c in S]].corr(method="spearman")
    R.to_csv(out / "rsv_gse155237_results.csv", index=False)
    S.join(md.reset_index(drop=True).set_index(S.index)).to_csv(out / "rsv_module_scores.csv")
    pd.set_option("display.width", 220)
    print("\n", R[[c for c in ("module", "tier", "predicted", "n", "or_per_sd", "ci_low",
                               "ci_high", "p_one_sided", "p_holm", "empirical_p",
                               "null_frac_p05", "verdict") if c in R]].round(4).to_string(index=False))
    print("\nmodule correlation (Spearman):\n", corr.round(2).to_string())


if __name__ == "__main__":
    main()
