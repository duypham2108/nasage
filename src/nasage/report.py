"""Assemble result tables and a plain-language markdown report."""
from __future__ import annotations

from pathlib import Path

import pandas as pd

HEADER = """# Age-associated nasal mucosal immune-state profile

**What this is.** Module scores for prespecified nasal mucosal programmes,
compared across age bins *within* annotated cell populations, with the donor as
the unit of replication.

**What this is not.** It is not a vaccine-readiness score, not a ranking of age
groups, and not evidence that any age group would respond better to a mucosal
vaccine. Baseline transcription does not measure vaccine response; that needs
vaccination and an immunogenicity or protection endpoint.
"""


def _write(outdir, *, design, coverage, eligibility, omnibus, contrasts,
                 composition, permutation=None, pooled_vs_within=None,
                 scenario_label="real data", extra_notes=()):
    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    tables = dict(module_coverage=coverage, donor_celltype_eligibility=eligibility,
                  age_effects_omnibus=omnibus, age_effects_contrasts=contrasts,
                  composition_results=composition)
    if permutation is not None:
        tables["permutation_null"] = permutation
    if pooled_vs_within is not None:
        tables["pooled_vs_within"] = pooled_vs_within
    for name, df in tables.items():
        df.to_csv(outdir / f"{name}.csv", index=False)

    prim = omnibus[(omnibus.tier == "primary") & omnibus.estimable]
    sig = prim[prim.p_adj_BH < 0.05].sort_values("p_adj_BH")
    nonest = omnibus[(omnibus.tier == "primary") & (~omnibus.estimable)]
    poor = coverage[(coverage.poorly_measurable) & (coverage.applicability > 0)]

    lines = [HEADER, f"\nDataset: **{scenario_label}**\n",
             f"- donors: {design.shape[0]}",
             f"- primary family: {len(prim)} estimable (module, cell type) combinations",
             f"- not estimable (too few donors with enough cells): {len(nonest)}",
             f"- flagged poorly measurable and excluded: {len(poor)}", ""]
    lines.append("## Age-associated differences (BH-FDR < 0.05 in the primary family)\n")
    if len(sig):
        lines.append(sig[["module_id", "cell_type", "n_donors", "f_stat",
                          "p_omnibus", "p_adj_BH"]].to_markdown(index=False))
    else:
        lines.append("No combination reached BH-FDR < 0.05. With this donor count that is "
                     "an absence of evidence, not evidence of no difference.")
    lines.append("\n## Combinations that could not be assessed\n")
    lines.append("These are reported as missing, not as zero and not as absence of the "
                 "biology. A population a brush fails to recover tells you about sampling.\n")
    if len(nonest):
        lines.append(nonest[["module_id", "cell_type", "n_donors", "donors_per_bin"]].to_markdown(index=False))
    lines.append("\n## Composition\n")
    clr = composition[composition["transform"] == "CLR"].sort_values("p_adj_BH")
    lines.append(clr[["cell_type", "mean_proportion", "p_omnibus", "p_adj_BH",
                      "n_donors_with_zero"]].round(4).to_markdown(index=False))
    lines.append("\nProportions are compositional: they sum to one, so a change in one "
                 "population forces changes in the others, and none of them measures "
                 "absolute abundance in the mucosa.\n")
    if permutation is not None:
        lines.append("## Permutation calibration\n")
        lines.append(f"- per-test false-positive rate at p < 0.05: "
                     f"{permutation.attrs.get('per_test_false_positive_rate', float('nan')):.3f}")
        lines.append(f"- permutations with at least one BH-significant result: "
                     f"{permutation.attrs.get('family_wise_rate_BH05', float('nan')):.3f}\n")
    lines.append("## Interpretation limits\n")
    lines += ["- Module scores are relative within a cell; they are not absolute "
              "expression and are not comparable across modules.",
              "- A higher barrier, mucociliary or interferon score is **not** automatically "
              "favourable for mucosal vaccination; each has a plausible opposing effect.",
              "- Interferon tone largely reflects recent respiratory exposure; without "
              "pathogen testing it cannot be attributed to age.",
              "- Epithelial MHC class II is inducible and is not professional "
              "antigen-presenting function.",
              "- No composite readiness score is produced, by design."]
    lines += list(extra_notes)
    (outdir / "report.md").write_text("\n".join(lines), encoding="utf-8")
    return outdir / "report.md"


def write_report(adata, outdir, *, label="", pooled_vs_within=None):
    """Write report.md and the result tables stored by ``nas.tl`` to `outdir`.

    Needs ``tl.donor_summary``, ``tl.age_effects`` and ``tl.composition``; adds
    the permutation check if ``tl.permutation_null`` was run. Also writes the
    donor summary, QC report, score distribution and symbol mapping.
    """
    from . import aggregate
    from .tl import _need
    st = adata.uns["nasage"]
    summary = _need(adata, "donor_summary", "donor_summary")
    perm = st.get("permutation")
    if perm is not None:
        perm = perm.copy()
        perm.attrs.update(dict(st.get("permutation_rates", {})))
    path = _write(outdir, design=st["design"], coverage=st["coverage"],
                  eligibility=aggregate.eligibility_table(summary),
                  omnibus=_need(adata, "omnibus", "age_effects"), contrasts=st["contrasts"],
                  composition=_need(adata, "composition", "composition"), permutation=perm,
                  pooled_vs_within=pooled_vs_within, scenario_label=label or "real data")
    outdir = Path(outdir)
    summary.to_csv(outdir / "donor_celltype_scores.csv", index=False)
    for key, name in (("qc", "qc_report"), ("score_distribution", "score_distribution"),
                      ("null_control", "null_control"),
                      ("symbol_report", "symbol_mapping"), ("design", "donor_design")):
        if key in st:
            st[key].to_csv(outdir / f"{name}.csv", index=False)
    return path
