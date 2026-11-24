"""Literature-curated epithelial cell-type programs for the mucosal vaccine readiness screen.

Gene lists are fixed by docs/MVRI_PROGRAMS.md, which carries the reason for every gene and at
least one reference per program. This file is the machine-readable copy; the two must not
diverge, and neither is edited after an outcome has been seen.

Composition-proxy genes (PROXY below) are forbidden from every program by construction, so a
program score can be adjusted for its cell type's abundance without adjusting away itself.
"""

PROXY = {
    "Ciliated": ["FOXJ1", "TPPP3", "PIFO", "SNTN", "CAPS"],
    "Basal": ["KRT5", "KRT14", "TP63", "DLK2", "MIR205HG"],
    "Secretory": ["SCGB1A1", "BPIFA1", "MSMB", "WFDC2", "SLPI"],
}

PROGRAMS = {
    "P1_basal_progenitor_repair": dict(cell_type="Basal", genes=[
        "ITGA6", "ITGB4", "NGFR", "COL17A1", "LAMB3", "LAMA3", "LAMC2",
        "DST", "PLEC", "CDH3", "BCAM", "SOX2", "MKI67", "TOP2A"]),
    "P2_junction_integrity": dict(cell_type="Basal", genes=[
        "CLDN1", "CLDN3", "CLDN4", "CLDN7", "OCLN", "TJP1", "TJP2", "TJP3",
        "CGN", "F11R", "CDH1", "DSP", "PKP3", "JUP"]),
    "P3_secretory_antimicrobial": dict(cell_type="Secretory", genes=[
        "LTF", "LYZ", "BPIFB1", "LCN2", "DEFB1", "DEFB4A", "CST3", "PI3",
        "CP", "AZGP1", "S100A8", "S100A9"]),
    "P4_iga_transport": dict(cell_type="Secretory", genes=[
        "PIGR", "FCGRT", "CCL28", "TNFSF13", "TNFSF13B", "IL6", "CXCL17",
        "JCHAIN", "IL17RA", "IL17RC"]),
    "P5_mucin_barrier": dict(cell_type="Secretory", genes=[
        "MUC5AC", "MUC5B", "MUC4", "MUC1", "MUC16", "TFF1", "TFF3", "AGR2",
        "CLCA1", "FCGBP", "ZG16B"]),
    "P6_entry_factor_availability": dict(cell_type="Ciliated", genes=[
        "ST6GAL1", "ST3GAL4", "ST3GAL6", "ST6GALNAC1", "NEU1", "NEU3",
        "TMPRSS2", "TMPRSS4", "TMPRSS11D", "TMPRSS11E", "CTSL", "FURIN",
        "ANPEP", "BSG"]),
    "P7_ciliary_beat_machinery": dict(cell_type="Ciliated", genes=[
        "DNAH5", "DNAH9", "DNAI1", "DNAI2", "DNAL1", "DNAAF1", "DNAAF2",
        "CCDC39", "CCDC40", "RSPH1", "RSPH4A", "RSPH9", "SPAG6", "SPAG1",
        "HYDIN", "ARMC4", "CFAP43", "CFAP44"]),
}

_forbidden = {g for v in PROXY.values() for g in v}
for _n, _p in PROGRAMS.items():
    _bad = sorted(set(_p["genes"]) & _forbidden)
    assert not _bad, f"{_n} contains composition-proxy genes: {_bad}"
    assert len(_p["genes"]) == len(set(_p["genes"])), f"{_n} has duplicate genes"
