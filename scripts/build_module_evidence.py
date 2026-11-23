"""Regenerate docs/MODULE_EVIDENCE.md from gene_modules.json + module_evidence.csv.

Run from anywhere:  python scripts/build_module_evidence.py
"""
import json
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
R = ROOT / "src" / "nasage" / "resources"
OUT = ROOT / "docs" / "MODULE_EVIDENCE.md"
M = json.loads((R / "gene_modules.json").read_text(encoding="utf-8"))
A = json.loads((R / "module_celltype_applicability.json").read_text(encoding="utf-8"))
ev = pd.read_csv(R / "module_evidence.csv", dtype=str).fillna("")

ROLE = {"definition": "Definition", "function": "Function", "airway_relevance": "Airway / nasal",
        "age_relevance": "Age", "caveat": "Caveat", "method": "Method"}
ROLE_ORDER = list(ROLE)
VER = {"abstract_checked": "abstract", "full_text_checked": "full text",
       "metadata_only": "metadata only"}

GAPS = {
    "M01": ["No verified human reference links airway barrier gene expression to age.",
            "The desmosomal genes (DSP, PKP1, PKP3, JUP) are not specifically covered by the "
            "references below."],
    "M02": ["The FOXJ1 master-regulator evidence is from zebrafish.",
            "Individual axonemal genes were not checked against a published gene list."],
    "M02b": ["SCGB1A1 is a weak secretory marker in nasal (vs bronchial) epithelium."],
    "M03b": ["No verified age reference."],
    "M04": ["TLR4, TLR5 and TLR9 have no airway-epithelium-specific reference here.",
            "TLR7 is barely expressed by airway epithelial cells (cell-line evidence)."],
    "M05": ["That type III interferon dominates at the airway epithelium rests largely on mouse "
            "data."],
    "M06": ["The HLA mapping-bias reference is from DNA sequencing; applying it to scRNA-seq "
            "expression is an extrapolation.",
            "The age row is an inference from the ISG age difference, not a direct MHC-I "
            "measurement."],
    "M07": ["LGMN has no verified supporting reference."],
    "M08": ["CD83, LAMP3, ICOSLG, TNFSF9, CD70 and TNFSF4 are not named in the abstracts checked.",
            "The age reference concerns innate immunity generally, not DCs specifically."],
    "M09": ["No gene-specific reference for MZB1, DERL3, SEC11C, SSR4, FKBP11 or IRF4.",
            "The XBP1, PRDM1 and BCMA references are mouse studies.",
            "No verified age reference."],
    "M10": ["No verified age reference."],
    "M11": ["SenMayo deliberately excludes CDKN1A and CDKN2A; their inclusion here follows the "
            "consensus review, not SenMayo.",
            "SERPINE1 and GDF15 membership of SenMayo could not be checked (supplementary list "
            "not accessible); GDF15 and SERPINs are supported by the ageing-plasma SASP atlas."],
    "M12": ["These cytokines are not specific to ageing; infection is a strong confounder."],
    "M13": ["No verified age reference.",
            "That tissue dissociation induces CD69 could not be verified (no PubMed abstract)."],
}


def cite(r):
    link = f"[PMID {r.pmid}](https://pubmed.ncbi.nlm.nih.gov/{r.pmid}/)" if r.pmid else "no PMID"
    return f"{r.first_author} et al. {r.year}, *{r.journal}* ({link})"


def table(rows):
    rows = rows.assign(_o=rows.evidence_role.map({r: i for i, r in enumerate(ROLE_ORDER)}))
    out = ["| Role | What it supports | Reference | Checked |", "|---|---|---|---|"]
    for r in rows.sort_values("_o", kind="stable").itertuples():
        # drop internal "Claim 10:" / "Claim 1 (partial):" labels from the literature check
        claim = re.sub(r"^Claim \d+(\s*\([^)]*\))?:\s*", "", r.claim_supported).replace("|", "/")
        out.append(f"| {ROLE.get(r.evidence_role, r.evidence_role)} | {claim} | {cite(r)} | "
                   f"{VER.get(r.verification, r.verification)} |")
    return "\n".join(out)


def tested(mod):
    d = A[mod]
    prim = [c for c, v in d.items() if v == 1.0]
    expl = [c for c, v in d.items() if v == 0.5]
    s = f"primary in {', '.join(prim) if prim else 'none'}"
    return s + (f"; exploratory in {', '.join(expl)}" if expl else "")


L = ["# Module evidence",
     "",
     "What each module is meant to read out, and the literature behind it. Generated from "
     "`src/nasage/resources/gene_modules.json` and "
     "`src/nasage/resources/module_evidence.csv`; edit those, not this file.",
     "",
     "**How the references were checked.** Every PMID, DOI, year, journal and title below was "
     "confirmed against PubMed metadata on 2026-11-23. The *Checked* column says how far the "
     "claim itself was checked: against the abstract, the full text, or only the metadata "
     "(the title). Each module also lists the **evidence gaps** found: genes or claims that no "
     "verified reference supports. Treat those as open, not settled.",
     "",
     "Method references cited elsewhere in the docs (Scanpy, muscat, scCODA, propeller, "
     "Splatter, powsimR, DESeq2, Benjamini-Hochberg and similar) are in `references.csv` with "
     "their original provenance; they were not re-checked in this pass.",
     ""]
L += ["## Contents", ""]
for mod, meta in M.items():
    L.append(f"- [{mod.split('_')[0]} {meta['name']}](#{mod.lower()})")
L.append("- [General: cell types, QC and methods](#general)")
L.append("")

for mod, meta in M.items():
    short = mod.split("_")[0]
    L += [f'<a id="{mod.lower()}"></a>', f"## {short} {meta['name']}", "",
          f"**Genes ({len(meta['genes'])}):** {', '.join(meta['genes'])}  ",
          f"**Tested:** {tested(mod)}  ",
          f"**Reading for mucosal vaccination:** {meta.get('interpretation', '')}", "",
          f"**What it measures.** {meta['definition']}", "",
          f"**How to read it.** {meta.get('interpretation_note', '')}", "",
          f"**Caveats.** {meta.get('limitations', '')}", ""]
    rows = ev[ev.module_id == short]
    if len(rows):
        L += ["### Evidence", "", table(rows), ""]
    if short in GAPS:
        L += ["### Evidence gaps", ""] + [f"- {g}" for g in GAPS[short]] + [""]

L += ['<a id="general"></a>', "## General: cell types, QC and methods", "",
      table(ev[ev.module_id == "GENERAL"]), "",
      "### Evidence gaps", "",
      "- That neutrophils have low gene counts in droplet scRNA-seq, and are frequent in healthy "
      "nasal brushings, could not be verified from abstracts; the references show only that "
      "neutrophils are technically hard to capture and that nasal granulocytes rise in COVID-19.",
      "- The MAD-based QC reference supports the general method; its abstract does not mention "
      "MAD thresholds explicitly.", ""]
OUT.write_text("\n".join(L), encoding="utf-8")
print("wrote", OUT, len(L), "lines")
