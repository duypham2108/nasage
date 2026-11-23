# How to read these outputs (and how not to)

## The one rule
A module score is **relative expression of a gene set within a cell, compared
with that cell's other genes**. It is not a rate, not a protein level, not a
function. Two modules' scores are not comparable to each other, and the same
module is only comparable across donors *within the same cell population*.

## Why no composite "vaccine readiness" score
The programmes measured here pull in different directions, and the direction of
benefit depends on the vaccine platform:

| Programme | Plausible supportive reading | Plausible opposing reading |
|---|---|---|
| Barrier / junctions | keeps pathogens out | keeps antigen out; many mucosal adjuvants work by loosening junctions |
| Mucociliary machinery | clears pathogens | clears the vaccine too, shortening mucosal residence |
| Mucins | traps antigen for sampling | physically impedes delivery and uptake |
| PIGR transport | required for secretory IgA to reach the lumen | reports capacity only, not antigen-specific antibody |
| Innate sensing | responsiveness to TLR/STING adjuvants | reactogenicity |
| Interferon-stimulated genes (type I and III) | antiviral state | restricts live-attenuated replication; marks recent infection; immunopathology |
| MHC-I / MHC-II | antigen presentation capacity | largely interferon-inducible, so partly the same signal as IFN tone |
| Costimulation | limiting factor for naive T-cell priming | may mark ongoing inflammation |
| Plasma-cell secretory machinery | local antibody secretion capacity | depends entirely on which cells the brush happened to recover |
| IgA fraction | secretory IgA is the mucosal antibody that pIgR transports | a fall may reflect a shift to IgG rather than less antibody; says nothing about specificity |
| Senescence arrest | marks cells that may impair repair and responses with age | p21 also rises with DNA damage and interferon; no single marker defines senescence |
| Inflammatory SASP / inflammaging | chronic low-grade inflammation is an ageing hallmark | the same cytokines rise acutely with infection, allergy, smoking and the brushing itself |
| Resident memory T cells | rapid local protection at the mucosa | CD69 also marks activation; CD103 favours CD8 cells |

Collapsing these into one number hides the trade-offs and invents a ranking the
data cannot support. Ranking age groups for vaccination requires vaccination:
immunogenicity or protection endpoints, by platform. This package deliberately
ships no composite score and no ranking function.

## Things that will bite you on real data
1. **Composition masquerading as cell state.** Pooling epithelial cells makes a
   shift in the ciliated:secretory ratio look like a change in cilia or mucin
   programmes. Always compare within annotated populations; use
   `nas.tl.pooled_scores` only to show the contrast.
2. **Absent immune cells are a sampling result.** A brush that returns no
   plasma cells tells you about the brush. Record it as missing, never as zero,
   and never as "the biology is absent".
3. **Interferon tone is exposure.** Baseline ISG scores track recent or ongoing
   respiratory viral infection. Without pathogen testing, an age difference in
   ISG tone is not an age *trait*.
4. **Ambient RNA.** Nasal brushings carry abundant secreted transcripts
   (mucins, secretoglobins, immunoglobulins). Ambient contamination inflates
   exactly the modules anchored on secreted products. Run an ambient-RNA
   correction and check whether M02b, M03, M09 and M10 survive it.
5. **Epithelial MHC class II** is inducible and is not professional
   antigen-presenting function. It is scored as an exploratory combination.
6. **Donors, not cells.** With cells as the unit, every test is anticonservative
   by orders of magnitude. The donor count in each output table is the real
   sample size.
7. **Age is confounded** with season, exposure history, smoking, allergic
   disease and sampling site. Adjusting for batch and sex does not address any
   of these.
8. **Interferon couples modules.** Interferon induces ISGs (M05), the MHC-I
   machinery (M06) and part of innate sensing (M04). Check
   `nas.tl.module_correlation` and report modules that move together across
   donors as one signal.
9. **Subtype mix within a population.** Goblet vs club, mature vs immature DCs,
   CD4 vs CD8 T cells: keep these as your labels so the model can adjust for
   each donor's mix (`subfrac_1/2`). Otherwise a shift in mix reads as a change
   in cell state.
10. **HLA genes.** M06 and M07 contain polymorphic HLA genes that align to the
    reference with allele-dependent efficiency; rerun with `nas.drop_hla()`.
11. **Ceiling and floor.** Check `score_distribution` in `adata.uns`: a module near the
    UCell ceiling cannot show increases, and one that is zero in most cells
    carries little information per cell.
