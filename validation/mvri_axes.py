"""The five MVRI axes, transcribed verbatim from slides 9-14 of the deck.

Lifted unchanged out of `validation/mvri_presentation_audit.py` so that the figure scripts and
the audit score the same thing. The audit script keeps its own copy deliberately: it is the
locked, committed record of what was tested, and editing it would change the provenance of
results already written. Anything importing from here asserts against `axis_audit.csv` that
the two copies still agree -- see `check_against_audit`.
"""
import numpy as np
import pandas as pd

AGE_ORDER = ["Neonate", "Infant", "Young child", "Child", "Adolescent", "Adult", "Eldly"]

# Display names. "Eldly" is the spelling in the object's Age_bin column and must stay that way
# for every lookup; only the figure labels use the corrected spelling.
AGE_LABEL = {b: b for b in AGE_ORDER} | {"Eldly": "Elderly"}

AXES = {
    "Mucus permeability": dict(
        genes=["MUC1", "MUC4", "MUC5AC", "MUC5B", "MUC16",
               "GALNT4", "PLAU", "BPIFA2", "LYZ", "SCGB1A1",
               "CLDN1", "CLDN4", "OCLN", "TJP1", "DSG3"],
        subset=["Secretory", "Goblet", "Club", "Basal"], weight=-0.20),
    "Antigen presentation": dict(
        genes=["FCGR1A", "FCGR2A", "FCGR2B", "FCGR3A", "SCARB1", "LRP1", "LY6E",
               "HLA-A", "HLA-B", "HLA-C", "HLA-DRA", "HLA-DRB1", "CD74",
               "CD40", "CD80", "IFNGR1", "IFNGR2"],
        subset=["Secretory", "Goblet", "DC", "Mac", "Mono"], weight=+0.25),
    "Innate sensing": dict(
        genes=["TLR1", "TLR2", "TLR3", "TLR4", "TLR5", "TLR6", "TLR7", "TLR8", "TLR9",
               "TLR10", "DDX58", "IFIH1", "CGAS", "TMEM173", "IRF3", "IRF7", "IRF9",
               "NFKB1", "NFKB2", "RELA"],
        subset=["__EPI__", "DC", "Mac", "Mono"], weight=+0.25),
    "IgA support": dict(
        genes=["PIGR", "TGFB1", "TNFSF13", "TNFSF13B"],
        subset=["Secretory", "Goblet", "Plasma", "B mem", "DC"], weight=+0.20),
    "Cytokine tone": dict(
        genes=["IL33", "TSLP", "IL6", "CXCL8", "IL8", "CXCL1", "CXCL10", "CCL20",
               "ISG15", "IFIT1", "IFIT2", "IFIT3", "MX1", "OAS1"],
        subset=["__EPI__", "__MYELOID__"], weight=+0.10),
}

ALIASES = {"CGAS": ["CGAS", "MB21D1", "C6orf150"], "TMEM173": ["TMEM173", "STING1"],
           "CXCL8": ["CXCL8", "IL8"], "IL8": ["CXCL8", "IL8"]}

# The three axes that clear the 60% measurability gate, in the order the figure draws them.
MEASURABLE = ["Mucus permeability", "Antigen presentation", "IgA support"]
UNMEASURED = ["Innate sensing", "Cytokine tone"]


def check_against_audit(audit_csv):
    """Fail loudly if this copy of the axes has drifted from the locked audit's."""
    R = pd.read_csv(audit_csv).set_index("axis")
    assert list(R.index) == list(AXES), (
        f"axis names differ from the audit: {list(R.index)} vs {list(AXES)}")
    for name, ax in AXES.items():
        n_listed = len({g for g in ax["genes"] if g != "IL8"})
        assert int(R.loc[name, "n_genes_listed"]) == n_listed, (
            f"{name}: {n_listed} genes here, {R.loc[name, 'n_genes_listed']} in the audit")
    return R


def resolve(genes, pos):
    """Gene symbols -> column indices, applying the alias table, order preserved."""
    out = []
    for g in genes:
        for cand in ALIASES.get(g, [g]):
            if cand in pos and pos[cand] not in out:
                out.append(pos[cand])
                break
    return out


def score_cells(idx, Xs):
    """Mean of per-gene z within the subset -- the deck's own score.

    Takes an already-sliced CSC subset matrix: re-slicing the full matrix per call exhausts
    memory on the largest subset (Epi + Myeloid, 166M nonzeros).
    """
    M = np.asarray(Xs[:, idx].todense())
    sd = M.std(0, ddof=1)
    return ((M - M.mean(0)) / np.where(sd > 0, sd, 1)).mean(1)
