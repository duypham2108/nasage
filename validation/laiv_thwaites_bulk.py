"""NasAge modules in bulk nasal curettage RNA-seq after LAIV (GSE230494; Thwaites et al. 2023).

Implements docs/PREREG_vaccine_axes.md (locked in commit 1988050):
  * scoring: per module, mean over genes of z-scored log2(CPM+1) (z relative to
    baseline samples); genes with mean CPM < 1 across baseline dropped; UCell on
    counts as sensitivity; samples with < 5 M assigned reads excluded;
  * manipulation check: M05 rises from baseline to 72 h within participant;
  * baseline scores + covariates (log10 reads, PTPRC) saved for H1-H4, which need a
    participant-level link to outcomes that the public data do not provide.

Usage: python validation/laiv_thwaites_bulk.py <data_dir> <symbol_map.h5ad> <out_dir>
"""
import gzip
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import nasage as nas

data, map_h5ad, out = Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3])
out.mkdir(parents=True, exist_ok=True)

# ---- counts: Ensembl gene x sample -> symbol x sample
cols = {}
for f in sorted((data / "counts").glob("*_counts.txt.gz")):
    name = re.sub(r"^GSM\d+_|_counts\.txt\.gz$", "", f.name)
    with gzip.open(f, "rt") as fh:
        s = pd.read_csv(fh, sep="\t", header=None, index_col=0, skiprows=1).iloc[:, 0]
    cols[name] = s
C = pd.DataFrame(cols)
C.index = C.index.str.split(".").str[0]
import anndata as ad
var = ad.read_h5ad(map_h5ad, backed="r").var
sym = var["feature_name"].astype(str).reindex(C.index)
print(f"{C.shape[1]} samples, {C.shape[0]} genes; {sym.notna().mean():.1%} mapped to symbols")
C = C[sym.notna()].groupby(sym[sym.notna()].values).sum()

meta = pd.DataFrame(index=C.columns)
meta["participant"] = meta.index.str.extract(r"(LAIV_\d+)", expand=False).values
meta["timepoint"] = np.where(meta.index.str.contains(r"_(?:1|Visit1)$"), "baseline", "SD3")
meta["reads"] = C.sum(0)
keep = meta.reads >= 5e6
print("excluded (<5M reads):", list(meta.index[~keep]))
meta, C = meta[keep], C.loc[:, keep]

cpm = C / C.sum(0) * 1e6
base = meta.timepoint == "baseline"
logc = np.log2(cpm + 1)
mu, sd = logc.loc[:, base].mean(1), logc.loc[:, base].std(1, ddof=1)
expressed = (cpm.loc[:, base].mean(1) >= 1) & (sd > 0)
Z = logc[expressed].sub(mu[expressed], axis=0).div(sd[expressed], axis=0)

M = nas.load_modules()
M_nohla = nas.drop_hla(M)
rows, cover = {}, []
for tag, mods in (("", M), ("_noHLA", M_nohla)):
    for m, d in mods.items():
        if tag and m not in ("M05_interferon_stimulated_genes", "M06_mhc_I_processing"):
            continue
        g = [x for x in d["genes"] if x in Z.index]
        cover.append(dict(module=m + tag, n_genes=len(d["genes"]), n_used=len(g)))
        if g:
            rows[m + tag] = Z.loc[g].mean(0)
S = pd.DataFrame(rows)
cover = pd.DataFrame(cover)

# UCell sensitivity (bulk samples ranked like cells)
try:
    from pyucell import compute_ucell_scores
    from nasage.scoring import MAX_RANK
    u = ad.AnnData(X=C.T.values.astype(np.float32), obs=pd.DataFrame(index=C.columns),
                   var=pd.DataFrame(index=C.index))
    sigs = {m: [g for g in d["genes"] if g in C.index] for m, d in M.items()
            if len(d["genes"]) > 1}
    compute_ucell_scores(u, sigs, max_rank=MAX_RANK, suffix="_UCell", n_jobs=1)
    S = S.join(u.obs)
except Exception as exc:                                     # pragma: no cover
    print("UCell sensitivity skipped:", exc)

meta["log10_reads"] = np.log10(meta.reads)
meta["PTPRC_log2cpm"] = logc.loc["PTPRC"] if "PTPRC" in logc.index else np.nan
scores = meta.join(S)
scores.to_csv(out / "module_scores_all_samples.csv")
scores[base[keep]].to_csv(out / "baseline_scores_for_H1_H4.csv")
cover.to_csv(out / "module_gene_coverage.csv", index=False)

# ---- manipulation check (part 1): within-participant rise of M05, baseline -> SD3
res = []
for col in [c for c in S.columns]:
    w = scores.pivot_table(index="participant", columns="timepoint", values=col).dropna()
    if len(w) < 5:
        continue
    d = w["SD3"] - w["baseline"]
    t = stats.ttest_rel(w["SD3"], w["baseline"])
    wx = stats.wilcoxon(d)
    res.append(dict(score=col, n_pairs=len(w), mean_change=d.mean(), median_change=d.median(),
                    frac_increased=(d > 0).mean(), t_p=t.pvalue, wilcoxon_p=wx.pvalue,
                    dz=d.mean() / d.std(ddof=1)))
R = pd.DataFrame(res).sort_values("t_p")
R.to_csv(out / "baseline_to_SD3_paired_change.csv", index=False)
pd.set_option("display.width", 200)
print(meta.groupby("timepoint").size().to_dict(), "| pairs:", int(R.n_pairs.max()))
print(cover.to_string(index=False))
print(R.round(4).to_string(index=False))
