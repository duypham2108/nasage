"""Does panel v2 replicate outside the two cohorts it was derived from?

Implements docs/PREREG_panel_v2_heldout.md. v2 collapsed five modules into one declared axis
(A01) on the basis of Ziegler/Loske structure. Here the same structure is read in three cohorts
that contributed nothing to that derivation: GSE336399, Thwaites LAIV, and Crownlands.

No outcome variable is read anywhere in this script.

Usage: python validation/panel_v2_heldout.py <results_dir> <out_dir>
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from nasage.modules import load_modules        # noqa: E402
from nasage.stats import effective_tests       # noqa: E402

R, OUT = Path(sys.argv[1]), Path(sys.argv[2])
OUT.mkdir(parents=True, exist_ok=True)

A01 = ["M05_interferon_stimulated_genes", "M06_mhc_I_processing", "M07_mhc_II_presentation",
       "M09_plasma_secretory_machinery", "M03_pigr_antibody_transport"]
V1_PRIMARY = [m for m in load_modules() if m.startswith("M")]
RHO_HELDOUT = 0.60
N_PERM = 1000


def module_columns(df):
    """Module score columns, tolerating the _UCell / _logexpr suffixes used in single cell."""
    out = {}
    for c in df.columns:
        base = c
        for suf in ("_UCell", "_logexpr", "_cpm", "_num", "_den"):
            if base.endswith(suf):
                base = base[: -len(suf)]
        if base.startswith("M") and pd.api.types.is_numeric_dtype(df[c]):
            if base not in out:          # first matching representation wins
                out[base] = c
    return out


def load_heldout():
    coh = {}
    f = R / "gse336399" / "module_scores.csv"
    if f.exists():
        coh["GSE336399"] = pd.read_csv(f)
    f = R / "thwaites2023_laiv" / "module_scores_all_samples.csv"
    if f.exists():
        d = pd.read_csv(f)
        # one row per participant (baseline only) so repeated timepoints do not inflate rho
        if "timepoint" in d.columns:
            base = d[d.timepoint.astype(str).str.upper().isin(["BASELINE", "D0", "SD0", "V1"])]
            d = base if len(base) >= 10 else d.groupby("participant", as_index=False).first()
        coh["Thwaites LAIV"] = d
    f = R / "crownlands" / "cell_scores_obs.csv.gz"
    if f.exists():
        d = pd.read_csv(f)
        mc = module_columns(d)
        keys = [k for k in ("donor", "ct") if k in d.columns]
        if keys:
            d = d.groupby(keys, as_index=False)[list(mc.values())].mean()
        coh["Crownlands"] = d
    return coh


rows, pred = [], []
for label, D in load_heldout().items():
    mc = module_columns(D)
    present = [m for m in A01 if m in mc]
    panel = [m for m in mc if m in V1_PRIMARY]
    if len(present) < 3 or len(panel) < 6:
        print(f"{label}: too few modules present ({len(present)} of A01, {len(panel)} panel)")
        continue
    S = D[[mc[m] for m in panel]].apply(pd.to_numeric, errors="coerce")
    S.columns = panel
    S = S.dropna(axis=1, how="all")
    panel = list(S.columns)
    present = [m for m in present if m in panel]
    C = S.corr(method="spearman")

    def med_rho(members):
        v = [C.loc[a, b] for i, a in enumerate(members) for b in members[i + 1:]]
        v = [x for x in v if np.isfinite(x)]
        return float(np.median(v)) if v else np.nan

    obs = med_rho(present)
    rng = np.random.default_rng(0)
    nullv = []
    for _ in range(N_PERM):
        pick = list(rng.choice(panel, size=len(present), replace=False))
        m = med_rho(pick)
        if np.isfinite(m):
            nullv.append(m)
    nullv = np.asarray(nullv)
    q95 = float(np.quantile(nullv, 0.95)) if nullv.size else np.nan
    perm_p = float((nullv >= obs).mean()) if nullv.size else np.nan

    # V3: effective tests, v1 panel vs v2 (A01 members contribute one column, their mean)
    Sv1 = S[panel]
    v2cols = {c: Sv1[c] for c in panel if c not in present}
    v2cols["A01_ifn_antigen_presentation"] = Sv1[present].mean(1)
    Sv2 = pd.DataFrame(v2cols)
    try:
        e1 = effective_tests(Sv1.corr(method="spearman"))
        e2 = effective_tests(Sv2.corr(method="spearman"))
        eff1, eff2 = float(e1["eff_li_ji"]), float(e2["eff_li_ji"])
    except Exception as exc:
        print(f"  {label}: effective_tests failed ({exc})")
        eff1 = eff2 = np.nan

    rows.append(dict(cohort=label, n_units=len(S), n_A01=len(present), n_panel=len(panel),
                     median_rho_A01=obs, null_q95=q95, perm_p=perm_p,
                     V1_pass=bool(obs >= RHO_HELDOUT), V2_pass=bool(perm_p < 0.05),
                     declared_v1=len(panel), eff_v1=eff1,
                     declared_v2=Sv2.shape[1], eff_v2=eff2,
                     V3_pass=bool(np.isfinite(eff1) and np.isfinite(eff2) and eff2 < eff1)))
    r = rows[-1]
    print(f"{label:16s} n={len(S):4d}  median rho(A01)={obs:+.3f} "
          f"(null q95 {q95:+.3f}, perm p={perm_p:.3f})  "
          f"eff v1={eff1:.1f} -> v2={eff2:.1f}")

T = pd.DataFrame(rows)
T.to_csv(OUT / "v2_heldout.csv", index=False)

print("\n" + "=" * 80)
print("PRE-SPECIFIED VERDICTS (docs/PREREG_panel_v2_heldout.md section 4)")
print("=" * 80)
if len(T):
    v1 = int(T.V1_pass.sum()); v2 = int(T.V2_pass.sum()); v3 = int(T.V3_pass.sum())
    n = len(T)
    print(f"V1  median rho >= {RHO_HELDOUT}:          {v1} of {n} cohorts "
          f"-> {'PASS' if v1 >= 2 else 'FAIL'}")
    print(f"V2  tighter than random 5-module sets: {v2} of {n} cohorts "
          f"-> {'PASS' if v2 >= 2 else 'FAIL'}")
    print(f"V3  effective tests fall under v2:     {v3} of {n} cohorts "
          f"-> {'PASS' if v3 == n else 'FAIL'}")
    print()
    if v1 >= 2 and v2 >= 2:
        print("VERDICT: the axis replicates AND is specific. v2's collapse is justified.")
    elif v1 >= 2:
        print("VERDICT: the members covary out of sample, but no more than any five modules do.")
        print("         This is the global-axis problem (B1), not an interferon-specific axis.")
        print("         v2 must be described as a global-axis correction, not a biological one.")
    else:
        print("VERDICT: the collapse does not replicate out of sample. v1 remains the panel")
        print("         of record and v2 is reported as unsupported.")
    print(T.round(3).to_string(index=False))
print("\nwrote", OUT)
