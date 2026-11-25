"""Supplementary Table 1: the gene-level evidence behind Figure 5a.

Joins three things for every gene drawn in Figure 5a:
  * which programme it belongs to and whether it cleared the 5% detection floor
    (`results/.../figure5/per_gene_detection.csv`, written by the audit);
  * the cell type that expresses it most widely **in this atlas**, computed per gene over
    the cells its own programme is scored in (`figure5/gene_celltype.csv`);
  * a one-line function and a literature reference, curated below.

The literature column is curated by hand and deliberately short: one indicative primary
paper or review per claim, not an exhaustive citation. Every reference here was retrieved
from PubMed and its DOI checked against the record; none is reconstructed from memory.

Usage: python validation/supp_table_genes.py <audit_dir> <out_dir>
"""
import sys
from pathlib import Path

import pandas as pd

AUD, OUT = Path(sys.argv[1]), Path(sys.argv[2])
F5 = AUD / "figure5"
OUT.mkdir(parents=True, exist_ok=True)

REFS = {
    "R1": ("Ma J, Rubin BK, Voynow JA. Mucins, mucus, and goblet cells. "
           "Chest 154:169-176 (2018).", "10.1016/j.chest.2017.11.008"),
    "R2": ("Kim KC. Role of epithelial mucins during airway infection. "
           "Pulm Pharmacol Ther 25:415-419 (2012).", "10.1016/j.pupt.2011.12.003"),
    "R3": ("Okuda K et al. Localization of secretory mucins MUC5AC and MUC5B in "
           "normal/healthy human airways. Am J Respir Crit Care Med 199:715-727 (2019).",
           "10.1164/rccm.201804-0734OC"),
    "R4": ("Gon Y, Hashimoto S. Role of airway epithelial barrier dysfunction in "
           "pathogenesis of asthma. Allergol Int 67:12-17 (2018).",
           "10.1016/j.alit.2017.08.011"),
    "R5": ("Gan H et al. Protein kinase D promotes airway epithelial barrier dysfunction "
           "and permeability through down-regulation of claudin-1. J Biol Chem "
           "288:37343-37354 (2013).", "10.1074/jbc.M113.511527"),
    "R6": ("Ganz T. Antimicrobial polypeptides. J Leukoc Biol 75:34-38 (2004).",
           "10.1189/jlb.0403150"),
    "R7": ("Almuntashiri S et al. Club cell secreted protein CC16: potential applications "
           "in prognosis and therapy for pulmonary diseases. J Clin Med 9:4039 (2020).",
           "10.3390/jcm9124039"),
    "R8": ("Khovidhunkit W et al. Parotid secretory protein is an HDL-associated protein "
           "with anticandidal activity. Am J Physiol Regul Integr Comp Physiol "
           "288:R1306-R1315 (2005).", "10.1152/ajpregu.00007.2004"),
    "R9": ("Kato K, Hansen L, Clausen H. Polypeptide N-acetylgalactosaminyltransferase-"
           "associated phenotypes in mammals. Molecules 26:5504 (2021).",
           "10.3390/molecules26185504"),
    "R10": ("Pishesha N, Harmand TJ, Ploegh HL. A guide to antigen processing and "
            "presentation. Nat Rev Immunol 22:751-764 (2022).",
            "10.1038/s41577-022-00707-2"),
    "R11": ("Mar KB et al. LY6E mediates an evolutionarily conserved enhancement of virus "
            "infection by targeting a late entry step. Nat Commun 9:3603 (2018).",
            "10.1038/s41467-018-06000-y"),
    "R12": ("Hvidberg V et al. Identification of the receptor scavenging hemopexin-heme "
            "complexes [LRP1/CD91]. Blood 106:2572-2579 (2005).",
            "10.1182/blood-2005-03-1185"),
    "R13": ("Grewal IS, Flavell RA. The role of CD40 ligand in costimulation and T-cell "
            "activation. Immunol Rev 153:85-106 (1996).",
            "10.1111/j.1600-065x.1996.tb00921.x"),
    "R14": ("Reuter U et al. Correction of complete interferon-gamma receptor 1 deficiency "
            "by bone marrow transplantation. Blood 100:4234-4235 (2002).",
            "10.1182/blood-2002-02-0433"),
    "R15": ("Blackburn JB et al. Secretory cells are the primary source of pIgR in small "
            "airways. Am J Respir Cell Mol Biol 67:334-345 (2022).",
            "10.1165/rcmb.2021-0548OC"),
    "R16": ("Wyatt TA et al. Malondialdehyde-acetaldehyde adduct formation decreases "
            "immunoglobulin A transport across airway epithelium in smokers who abuse "
            "alcohol. Am J Pathol 191:1732-1742 (2021).",
            "10.1016/j.ajpath.2021.06.007"),
    "R17": ("Tezuka H et al. Regulation of IgA production by naturally occurring "
            "TNF/iNOS-producing dendritic cells. Nature 448:929-933 (2007).",
            "10.1038/nature06033"),
}

ANN = {
    # gene: (function relevant to its programme, [refs])
    "MUC1": ("Membrane-tethered mucin of the periciliary brush; also limits and resolves "
             "airway inflammation during infection", ["R1", "R2"]),
    "MUC4": ("Membrane-tethered mucin of the periciliary layer", ["R1", "R2"]),
    "MUC16": ("Largest membrane-tethered mucin; inner, cell-associated barrier layer",
              ["R1", "R2"]),
    "MUC5AC": ("Gel-forming secreted mucin; concentrated in proximal airways and made by "
               "CCSP-positive secretory and goblet cells", ["R1", "R3"]),
    "MUC5B": ("Gel-forming secreted mucin; dominant mucin of healthy superficial "
              "epithelium and submucosal glands, required for mucociliary clearance",
              ["R1", "R3"]),
    "GALNT4": ("GalNAc-transferase that initiates mucin-type O-glycosylation, the "
               "modification that gives mucins their gel properties", ["R9"]),
    "BPIFA2": ("Secreted BPI-fold protein of the salivary/airway secretome with "
               "anticandidal activity", ["R8"]),
    "LYZ": ("Lysozyme; principal muramidase of airway surface liquid", ["R6"]),
    "SCGB1A1": ("Club cell secretory protein CC16; anti-inflammatory and anti-oxidant "
                "component of airway secretions", ["R7"]),
    "CLDN1": ("Tight-junction claudin; its loss raises paracellular permeability in airway "
              "epithelium", ["R4", "R5"]),
    "CLDN4": ("Tight-junction claudin of the apical junctional complex", ["R4", "R5"]),
    "OCLN": ("Occludin; integral tight-junction protein", ["R4", "R5"]),
    "TJP1": ("ZO-1; scaffold linking tight-junction strands to the actin cytoskeleton",
             ["R4", "R5"]),
    "HLA-A": ("MHC class I heavy chain; presents cytosolic peptides to CD8+ T cells",
              ["R10"]),
    "HLA-B": ("MHC class I heavy chain; presents cytosolic peptides to CD8+ T cells",
              ["R10"]),
    "HLA-C": ("MHC class I heavy chain; presents cytosolic peptides to CD8+ T cells",
              ["R10"]),
    "HLA-DRA": ("MHC class II alpha chain; presents endosomal peptides to CD4+ T cells",
                ["R10"]),
    "HLA-DRB1": ("MHC class II beta chain; presents endosomal peptides to CD4+ T cells",
                 ["R10"]),
    "CD74": ("Invariant chain; chaperones MHC class II and blocks the groove until "
             "endosomal peptide loading", ["R10"]),
    "CD40": ("Costimulatory receptor; CD40-CD40L engagement licenses antigen-presenting "
             "cells and supports B-cell class switching", ["R13"]),
    "IFNGR1": ("IFN-gamma receptor ligand-binding chain; complete deficiency causes "
               "susceptibility to mycobacterial and other intracellular infection",
               ["R14"]),
    "IFNGR2": ("IFN-gamma receptor signalling chain, partner of IFNGR1", ["R14"]),
    "LY6E": ("GPI-anchored interferon-stimulated protein that modulates a late step of "
             "viral entry", ["R11"]),
    "LRP1": ("CD91; endocytic scavenger receptor that internalises ligand complexes for "
             "intracellular processing", ["R12"]),
    "PIGR": ("Polymeric immunoglobulin receptor; transcytoses dimeric IgA to the mucosal "
             "surface and is made predominantly by secretory cells", ["R15"]),
    "TGFB1": ("Switch factor for IgA class-switch recombination; also down-regulates "
              "epithelial pIgR and so reduces secretory IgA transport", ["R16", "R17"]),
    "TNFSF13": ("APRIL; drives T-cell-independent IgA class-switch recombination at "
                "mucosal sites", ["R17"]),
}

R = pd.read_csv(AUD / "axis_audit.csv").set_index("axis")
GC = pd.read_csv(F5 / "gene_celltype.csv")
DET = pd.read_csv(F5 / "per_gene_detection.csv")
ORDER = [a for a in ["Mucus permeability", "Antigen presentation", "IgA support"]
         if a in set(GC.axis)]

rows = []
for a in ORDER:
    sub = GC[GC.axis == a].sort_values("detection_overall", ascending=False)
    for r in sub.itertuples():
        func, refs = ANN.get(r.gene, ("", []))
        rows.append(dict(
            programme=a, gene=r.gene,
            detection=f"{r.detection_overall:.2f}",
            top_cell_type=f"{r.top1} ({r.top1_det:.2f})",
            second_cell_type=f"{r.top2} ({r.top2_det:.2f})" if r.top2 else "",
            function=func, refs=", ".join(refs)))
T = pd.DataFrame(rows)
T.to_csv(OUT / "supp_table1_gene_evidence.csv", index=False)

used = sorted({x for r in T.refs for x in r.split(", ") if x},
              key=lambda k: int(k[1:]))
num = {k: i + 1 for i, k in enumerate(used)}


def esc(s):
    return (s.replace("&", r"\&").replace("%", r"\%").replace("_", r"\_"))


lines = []
for a in ORDER:
    rr = R.loc[a]
    lines.append(r"\multicolumn{5}{l}{\textbf{%s} \;---\; %d genes scored, "
                 r"empirical $P$ = %.3f, in %s} \\[1pt]"
                 % (esc(a), int(rr.n_genes_scored), float(rr.empirical_p),
                    esc(f"{int(rr.n_cells):,} cells")))
    lines.append(r"\midrule")
    for r in T[T.programme == a].itertuples():
        cites = ",".join(str(num[x]) for x in r.refs.split(", ") if x)
        lines.append(r"\textit{%s} & %s & %s & %s & [%s] \\" % (
            esc(r.gene), esc(r.detection), esc(r.top_cell_type),
            esc(r.function), cites))
    lines.append(r"\addlinespace")

reflist = "\n".join(
    r"\item[%d.] %s \textsc{doi}: \texttt{%s}" % (num[k], esc(REFS[k][0]), REFS[k][1])
    for k in used)

tex = r"""\documentclass[10pt]{article}
\usepackage[a4paper,margin=15mm,landscape]{geometry}
\usepackage[T1]{fontenc}
\usepackage{lmodern}
\usepackage{booktabs}
\usepackage{longtable}
\usepackage{array}
\usepackage{enumitem}
\usepackage{microtype}
\setlength{\parindent}{0pt}
\pagestyle{empty}

\begin{document}

{\large\textbf{Supplementary Table 1}\par}\vspace{2pt}
{\textbf{Gene-level evidence for the programmes in Figure 5a.}\par}\vspace{6pt}

{\footnotesize Every gene drawn in Figure 5a, that is, every gene detected in at least 5\%
of the cells its own programme is scored in. \textbf{Detection} is the fraction of those
cells with a non-zero count. \textbf{Principal cell type} is the cell type with the highest
detection of that gene \emph{in this atlas}, computed over the programme's own cell subset
and restricted to types with at least 200 cells; the value in parentheses is detection
within that type. The function column gives one indicative reference per claim and is not an
exhaustive citation. Genes below the detection floor, and the innate-sensing and
cytokine-tone programmes that fail the 60\% measurability gate, are excluded from Figure 5
and from this table.\par}

\vspace{8pt}
{\footnotesize
\begin{longtable}{@{}l r l p{85mm} l@{}}
\toprule
\textbf{Gene} & \textbf{Detection} & \textbf{Principal cell type (detection)} &
\textbf{Function relevant to the programme} & \textbf{Ref.} \\
\midrule
\endhead
PLACEHOLDER_ROWS
\bottomrule
\end{longtable}
}

\vspace{4pt}
{\footnotesize\textbf{References}\par}
{\scriptsize
\begin{itemize}[leftmargin=2.2em,itemsep=1pt,topsep=2pt]
PLACEHOLDER_REFS
\end{itemize}
}

\end{document}
"""
tex = tex.replace("PLACEHOLDER_ROWS", "\n".join(lines))
tex = tex.replace("PLACEHOLDER_REFS", reflist)
(OUT / "supp_table1.tex").write_text(tex, encoding="utf-8")
print(T.to_string(index=False))
print("\nwrote", OUT / "supp_table1_gene_evidence.csv")
print("wrote", OUT / "supp_table1.tex")
