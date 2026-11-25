"""Epithelial-only restriction of the curated programmes.

Plan: docs/PREREG_epithelial_only.md, locked at 198e01b.

The filter is derived from single-cell expression in Ziegler 2021 (a cohort contributing no
outcome to any test in this project) and is outcome-blind. Originals in mvri_programs.py are
NOT modified; these are distinct objects.

Rule: keep genes with epithelial fraction >= 0.5; keep a programme if >= 60% of its genes are
epithelial AND >= 8 genes remain.
"""
from pathlib import Path

import pandas as pd

from mvri_programs import PROGRAMS, PROXY

EPI_FRAC_MIN = 0.5
PROG_FRAC_MIN = 0.60
MIN_GENES = 8
_FRAC = Path("D:/research/nasage_validation/results/module_audit/gene_epithelial_fraction.csv")


def _load_fraction():
    if not _FRAC.exists():
        raise FileNotFoundError(
            f"{_FRAC} missing; run the epithelial-fraction step first (see plan section 1)")
    return pd.read_csv(_FRAC).set_index("gene")["epi_frac"].to_dict()


def build():
    frac = _load_fraction()
    kept, dropped, report = {}, {}, []
    for name, prog in PROGRAMS.items():
        genes = prog["genes"]
        # a gene with no measurement in the reference is kept: absence of evidence is not
        # evidence that it is immune, and dropping it would be a silent, untracked filter
        epi = [g for g in genes if frac.get(g, 1.0) >= EPI_FRAC_MIN]
        f = len(epi) / len(genes)
        ok = f >= PROG_FRAC_MIN and len(epi) >= MIN_GENES
        report.append(dict(program=name, n_genes=len(genes), n_epithelial=len(epi),
                           frac_epithelial=f, retained=ok,
                           removed=",".join(sorted(set(genes) - set(epi)))))
        if ok:
            kept[name + "_epi"] = dict(cell_type=prog["cell_type"], genes=epi,
                                       derived_from=name)
        else:
            dropped[name] = dict(frac_epithelial=f, n_remaining=len(epi))
    return kept, dropped, pd.DataFrame(report).sort_values("frac_epithelial")


PROGRAMS_EPI, DROPPED, FILTER_REPORT = build()

# the composition-proxy rule survives the filter
_forbidden = {g for v in PROXY.values() for g in v}
for _n, _p in PROGRAMS_EPI.items():
    _bad = sorted(set(_p["genes"]) & _forbidden)
    assert not _bad, f"{_n} contains composition-proxy genes: {_bad}"

if __name__ == "__main__":
    print(FILTER_REPORT.round(3).to_string(index=False))
    print("\nretained:", sorted(PROGRAMS_EPI))
    print("dropped :", DROPPED)
