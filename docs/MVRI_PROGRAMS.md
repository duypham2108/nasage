# Epithelial cell-type programs for the mucosal vaccine readiness screen

Status: **LOCKED** by the commit that adds this file. Gene lists are fixed here and are not
edited after any outcome is seen. Supersedes the unsupervised derivation, which failed
(`results/mvri/DERIVATION_v1_FAILED.md`: the only programs passing the size filter were the
mitochondrial / rRNA artifact block).

**Why curated rather than derived.** The unsupervised attempt recovered a technical axis, and the
one result in this project that has replicated externally — M02 mucociliary axoneme in Figure 5 —
is literature-curated and carries no selection bias. Curation also makes every gene defensible to
a reviewer, which a Leiden cluster is not.

**No outcome data has been read at the time of this commit.** Both validation cohorts
(GSE155237, GSE261277) remain untouched, so choosing curation over derivation cannot have been
influenced by them.

---

## Selection rules applied to every program

1. **Each gene has a stated reason and at least one reference below.** A gene that is only
   "co-expressed" does not qualify.
2. **Cell-type assignment** is the epithelial population where the programme is principally
   expressed in human airway, used to pick the composition proxy it is adjusted for.
3. **Composition-proxy genes are forbidden** and appear in no program: FOXJ1, TPPP3, PIFO, SNTN,
   CAPS (ciliated); KRT5, KRT14, TP63, DLK2, MIR205HG (basal); SCGB1A1, BPIFA1, MSMB, WFDC2,
   SLPI (secretory). This is what allows the primary model to separate programme activity from
   cell-type abundance.
4. **Redundancy with the existing curated modules** M01–M13 is computed, and a program exceeding
   Jaccard 0.5 is dropped before the family cap, per `PREREG_mvri.md` §2.2.
5. Programs are **not** tuned for size; each is as large as its biology and the proxy exclusions
   allow.

---

## P1 — Basal progenitor and repair capacity (Basal)

**Rationale.** Mucosal vaccine take requires an epithelium that can be transiently infected or
perturbed and then repaired. Basal cells are the airway progenitor pool, and their adhesion to
the basement membrane and hemidesmosome machinery is what licenses regeneration.

| Gene | Reason |
|---|---|
| ITGA6, ITGB4 | α6β4 integrin, the basal-cell basement-membrane receptor defining the progenitor niche |
| NGFR | p75NTR, marker of the regenerative basal subset |
| COL17A1 | hemidesmosome collagen; progenitor maintenance |
| LAMB3, LAMA3, LAMC2 | laminin-332, the ligand basal cells adhere to |
| DST, PLEC | hemidesmosome plakins linking integrin to keratin |
| CDH3, BCAM | basal adhesion |
| SOX2 | airway epithelial progenitor transcription factor |
| MKI67, TOP2A | proliferative fraction of the basal pool |

*References:* Gilpin et al. 2016 *Biomaterials* [doi:10.1016/j.biomaterials.2016.08.055](https://doi.org/10.1016/j.biomaterials.2016.08.055);
Wang et al. 2023 *J Allergy Clin Immunol* [doi:10.1016/j.jaci.2023.01.030](https://doi.org/10.1016/j.jaci.2023.01.030).

## P2 — Epithelial junction integrity (Basal)

**Rationale.** Barrier integrity governs both pathogen access and antigen sampling, and
respiratory viruses actively disrupt it. A loose barrier at baseline is a plausible readiness
axis in either direction, which is why the test is two-sided.

| Gene | Reason |
|---|---|
| CLDN1, CLDN3, CLDN4, CLDN7 | airway tight-junction claudins |
| OCLN | occludin |
| TJP1, TJP2, TJP3 | ZO-1/2/3 scaffolds |
| CGN | cingulin |
| F11R | JAM-A |
| CDH1 | E-cadherin, adherens junction |
| DSP, PKP3, JUP | desmosomal adhesion |

*References:* Kast et al. 2017 *Clin Exp Immunol* [doi:10.1111/cei.13042](https://doi.org/10.1111/cei.13042);
Looi et al. 2021 *Int J Mol Sci* [doi:10.3390/ijms22168790](https://doi.org/10.3390/ijms22168790).

## P3 — Secretory antimicrobial capacity (Secretory)

**Rationale.** Secreted antimicrobial proteins are the first non-specific barrier a live
intranasal vaccine or a challenge organism meets, and their levels track airway infection and
exacerbation in independent disease cohorts.

| Gene | Reason |
|---|---|
| LTF | lactoferrin; iron sequestration, directly measured in airway infection studies |
| LYZ | lysozyme; peptidoglycan hydrolysis |
| BPIFB1 | LPLUNC1, airway BPI-fold antimicrobial (BPIFA1 excluded as proxy) |
| LCN2 | lipocalin-2, siderophore sequestration |
| DEFB1, DEFB4A | β-defensins |
| CST3 | cystatin C, protease inhibition |
| PI3 | elafin, antiprotease/antimicrobial |
| CP | ceruloplasmin |
| AZGP1 | secretory zinc-α2-glycoprotein |
| S100A8, S100A9 | calprotectin subunits, metal sequestration |

*References:* Sagel et al. 2009 *Pediatr Pulmonol* [doi:10.1002/ppul.21028](https://doi.org/10.1002/ppul.21028);
Sibila et al. 2019 *Thorax* [doi:10.1136/thoraxjnl-2018-212895](https://doi.org/10.1136/thoraxjnl-2018-212895).

## P4 — IgA transport axis (Secretory) — the most directly vaccine-relevant program

**Rationale.** Secretory IgA is the mucosal effector a nasal vaccine is trying to induce, and
nasal IgA is an established correlate of protection against influenza challenge in humans,
including in volunteers with low serum titres. The epithelium does not make IgA, but it
**transports** it: pIgR is the rate-limiting step, and it is inducible by Th17 signals. An
epithelium with low baseline transport capacity cannot convert a plasma-cell response into
secreted antibody.

| Gene | Reason |
|---|---|
| PIGR | polymeric Ig receptor; the transcytosis step itself, rate-limiting |
| FCGRT | neonatal Fc receptor; IgG transcytosis arm |
| CCL28 | MEC; recruits IgA-secreting plasma cells to mucosa |
| TNFSF13, TNFSF13B | APRIL and BAFF; epithelial support of local plasma cells and class switching |
| IL6 | plasma-cell differentiation support |
| CXCL17 | mucosal chemokine of airway epithelium |
| JCHAIN | J chain; present in tissue, marks the polymeric Ig being transported |
| IL17RA, IL17RC | receptors for the Th17 signal shown to induce pIgR |

*References:* Jaffar et al. 2009 *J Immunol* [doi:10.4049/jimmunol.0900237](https://doi.org/10.4049/jimmunol.0900237);
Gould et al. 2017 *Front Microbiol* [doi:10.3389/fmicb.2017.00900](https://doi.org/10.3389/fmicb.2017.00900);
de Fays et al. 2022 *Cells* [doi:10.3390/cells11081324](https://doi.org/10.3390/cells11081324);
Bean et al. 2024 *mBio* [doi:10.1128/mbio.02372-23](https://doi.org/10.1128/mbio.02372-23);
Gerard et al. 2026 *Int J Mol Sci* [doi:10.3390/ijms27219525](https://doi.org/10.3390/ijms27219525).

## P5 — Gel-forming mucin barrier (Secretory)

**Rationale.** The mucus gel is the physical trap a live vaccine must penetrate to reach its
target cell. Nasal epithelial structure and function differ between children and elderly adults
in culture, which makes this axis plausibly age-dependent as well as readiness-relevant.

| Gene | Reason |
|---|---|
| MUC5AC, MUC5B | the two airway gel-forming mucins |
| MUC4, MUC1, MUC16 | membrane-tethered mucins of the periciliary barrier |
| TFF1, TFF3 | trefoil factors; mucus rheology and repair |
| AGR2 | required for mucin maturation |
| CLCA1 | mucin granule maturation |
| FCGBP | IgG Fc-binding protein, mucus crosslinking |
| ZG16B | secretory lectin |

*References:* Balázs et al. 2022 *Front Immunol* [doi:10.3389/fimmu.2022.822437](https://doi.org/10.3389/fimmu.2022.822437).

## P6 — Ciliated-cell entry-factor availability (Ciliated) — the novel axis

**Rationale.** This is the program with the clearest mechanistic link to a **live** intranasal
vaccine. LAIV must actually infect ciliated cells to take; influenza tropism in human airway
epithelium is for ciliated cells, mediated by α2,6-linked sialic acid, and entry requires
host-protease cleavage of haemagglutinin. Baseline availability of the sialyltransferases and
proteases that create those receptors and perform that cleavage is therefore a candidate
determinant of vaccine take that is **not** an interferon axis and is not represented in M01–M13.

| Gene | Reason |
|---|---|
| ST6GAL1 | makes α2,6-sialic acid, the human influenza receptor |
| ST3GAL4, ST3GAL6 | α2,3-sialylation, the avian-type receptor present in airway |
| ST6GALNAC1 | O-linked sialylation of mucins |
| NEU1, NEU3 | sialidases; remove the receptor, opposing the above |
| TMPRSS2, TMPRSS4, TMPRSS11D, TMPRSS11E | host serine proteases that cleave influenza HA |
| CTSL | cathepsin L, alternative cleavage route |
| FURIN | polybasic-site cleavage |
| ANPEP, BSG | additional respiratory virus entry factors expressed in airway epithelium |

*References:* Thompson et al. 2006 *J Virol* [doi:10.1128/JVI.00384-06](https://doi.org/10.1128/JVI.00384-06);
Laloli et al. 2022 *Front Immunol* [doi:10.3389/fimmu.2022.970325](https://doi.org/10.3389/fimmu.2022.970325).

## P7 — Ciliary beat machinery (Ciliated)

**Rationale.** Mucociliary clearance sets how long an inoculum dwells on the mucosa. Included
despite an expected overlap with the curated M02 module; the redundancy rule decides its fate,
not a judgement made after seeing results. FOXJ1, TPPP3, PIFO, SNTN and CAPS are excluded as
composition proxies, so this program is axonemal structure rather than ciliated-cell identity.

| Gene | Reason |
|---|---|
| DNAH5, DNAH9, DNAI1, DNAI2, DNAL1 | outer dynein arm; primary ciliary dyskinesia genes |
| DNAAF1, DNAAF2 | dynein assembly factors |
| CCDC39, CCDC40 | nexin–dynein regulatory complex ruler |
| RSPH1, RSPH4A, RSPH9 | radial spoke head |
| SPAG6, SPAG1 | central apparatus |
| HYDIN | central pair |
| ARMC4 | outer dynein arm docking |
| CFAP43, CFAP44 | cilia and flagella associated |

*References:* Thompson et al. 2006 *J Virol* (ciliated-cell tropism and clearance context)
[doi:10.1128/JVI.00384-06](https://doi.org/10.1128/JVI.00384-06);
Balázs et al. 2022 *Front Immunol* [doi:10.3389/fimmu.2022.822437](https://doi.org/10.3389/fimmu.2022.822437).

---

## Context references for the endpoints themselves

Not gene sources; they define what "readiness" means in each validation cohort and why these
endpoints were chosen.

- **RSV challenge endpoint** — Habibi et al. 2020 *Science*: neutrophilic inflammation in the
  respiratory mucosa predisposes to RSV infection.
  [doi:10.1126/science.aba9301](https://doi.org/10.1126/science.aba9301)
- **Pneumococcal carriage endpoint** — Morton et al. 2021 *EBioMedicine*, controlled human
  pneumococcal infection. [doi:10.1016/j.ebiom.2021.103579](https://doi.org/10.1016/j.ebiom.2021.103579)
- **Why nasal mucosal state should matter for vaccination** — de Steenhuijsen Piters et al. 2019
  *Nat Commun*, nasal microbiota and pneumococcus in the context of LAIV.
  [doi:10.1038/s41467-019-10814-9](https://doi.org/10.1038/s41467-019-10814-9)
- **Mucosal correlates of protection after influenza challenge** — Bean et al. 2024 *mBio*,
  from the same challenge programme as GSE253686.
  [doi:10.1128/mbio.02372-23](https://doi.org/10.1128/mbio.02372-23)

All references above were retrieved from PubMed.

---

## Pre-declared expectations, recorded before any test

These are **not** hypotheses with locked directions — `PREREG_mvri.md` §4 fixes every test as
two-sided — but recording what would be unsurprising keeps a post-hoc story from being told later.

| Program | Would be unsurprising |
|---|---|
| P4 IgA transport | higher transport capacity → less infection / less carriage |
| P6 entry factors | more entry-factor availability → more infection; and, in the exploratory LAIV endpoint, a **larger** response, since the vaccine must infect to take |
| P3 antimicrobial | higher → less carriage, pneumococcus especially |
| P1, P2, P5, P7 | no strong prior in either direction |

If P6 behaves oppositely between the challenge endpoints (where infection is bad) and the LAIV
endpoint (where infection is the point), that is the predicted pattern and not a contradiction —
stated here so it cannot be claimed as a discovery afterwards.
