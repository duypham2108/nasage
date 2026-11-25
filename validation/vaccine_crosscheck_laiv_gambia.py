"""Cross-omics check on the open Gambian LAIV dataset (Zenodo 10.5281/zenodo.14719593).

The open table links 244 children's vaccine outcomes to 26 baseline nasal GO
pathway scores (121 children with nasal data). The baseline nasal RNA-seq on
Zenodo 10.5281/zenodo.15224712 uses hashed IDs with no public key (linking needs
a data agreement with LSHTM), so NasAge modules cannot be linked to outcomes here.

This script asks: do the baseline nasal programmes closest to NasAge themes
(innate immune / defence regulation) relate to LAIV outcomes?
Logistic regression per GO term and outcome, adjusted for sex, study year and
baseline respiratory-virus positivity; BH within each outcome.

Caveat: the 26 GO terms are the published model's selected features, likely
chosen on these same children, so associations are descriptive, not validation.

Usage: python validation/vaccine_crosscheck_laiv_gambia.py <data_dir> <out_dir>
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

data, out = Path(sys.argv[1]), Path(sys.argv[2])
out.mkdir(parents=True, exist_ok=True)
c = pd.read_csv(data / "laiv_gambia_integrated.csv")
names = pd.read_csv(data / "laiv_gambia_go_names.csv", index_col=0)["name"]
go = [x for x in c.columns if x.startswith("nasal_baseline_go.")]
c["shed_any_v2"] = (c[["h1_v2_shed", "h3_v2_shed", "b_v2_shed"]].sum(1) > 0).astype(float)
c["n_strains_shed_v2"] = c[["h1_v2_shed", "h3_v2_shed", "b_v2_shed"]].sum(1)
outcomes = ["max_iga_responder", "max_HAI_responder", "shed_any_v2"]

d = c[c[go[0]].notna()].copy()
d["year"] = d.year.astype(str)
rows = []
for y in outcomes:
    for g in go:
        sub = d[[y, g, "sex", "year", "v0_resp_virus_positive"]].dropna().copy()
        sub["z"] = (sub[g] - sub[g].mean()) / sub[g].std()
        try:
            m = smf.logit(f"{y} ~ z + sex + year + v0_resp_virus_positive", data=sub).fit(disp=0)
            ci = m.conf_int().loc["z"]
            rows.append(dict(outcome=y, go=g.split("go.")[1], term=names.get(g, ""), n=len(sub),
                             odds_ratio_per_sd=float(np.exp(m.params["z"])),
                             ci_low=float(np.exp(ci[0])), ci_high=float(np.exp(ci[1])),
                             p=float(m.pvalues["z"])))
        except Exception as e:                                       # separation etc.
            rows.append(dict(outcome=y, go=g.split("go.")[1], term=names.get(g, ""),
                             n=len(sub), p=np.nan, note=str(e)[:80]))
res = pd.DataFrame(rows)
res["p_adj_BH"] = np.nan
for y in outcomes:
    s = (res.outcome == y) & res.p.notna()
    res.loc[s, "p_adj_BH"] = multipletests(res.loc[s, "p"], method="fdr_bh")[1]
res["nasage_theme"] = res.go.map({"0045088": "M04/M05 innate sensing & IFN",
                                  "0031347": "M04/M05 innate defence",
                                  "0050671": "M08/M13 lymphocyte activation",
                                  "0008630": "M11 DNA-damage / senescence-adjacent",
                                  "0060562": "M01 epithelial"}).fillna("")
res.sort_values(["outcome", "p"]).to_csv(out / "go_vs_outcomes.csv", index=False)
print(res[res.nasage_theme != ""].sort_values(["outcome", "p"])
      [["outcome", "term", "n", "odds_ratio_per_sd", "ci_low", "ci_high", "p", "p_adj_BH"]]
      .round(3).to_string(index=False))
print("\nBH<0.05 overall:", int((res.p_adj_BH < 0.05).sum()), "of", int(res.p.notna().sum()))
