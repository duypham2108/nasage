# Search for a third validation cohort for the P5 mucin-barrier result

Date: 2026-11-24. Recorded so the same search is not repeated, and so the choice of cohorts is
visible rather than assumed.

**What was wanted.** A human cohort with (a) a **baseline / pre-exposure** nasal transcriptome,
(b) an **infection, carriage or vaccine outcome**, (c) **no participant overlap** with
GSE155237 (Imperial RSV) or the Liverpool EHPC programme, and (d) n meaningfully above 39.

**Searched.** GEO DataSets, human, expression series, on: nasal + challenge; nasal + rhinovirus;
nasal + SARS-CoV-2 + challenge; "controlled human infection"; "human challenge" + nasal +
baseline; influenza + challenge + nasal + shedding; nasal + vaccine + IgA. Plus follow-up on
every plausible hit.

## Result of the search

| Dataset | Nasal baseline | Outcome | Verdict |
|---|---|---|---|
| GSE155237 RSV challenge, Imperial | 39 | infection | **in use** |
| GSE261277 pneumococcal, Liverpool | 39 | carriage | **in use** |
| **GSE253686 NOVA influenza challenge** | **75** | viral shedding | **best available — outcome behind a PMC download gate** |
| GSE166161 RSV D0/D7, Imperial | 13 | infection | **excluded: same donors as GSE155237** |
| GSE117580 Liverpool, LAIV vs TIV | 35 | carriage + vaccine arm | overlap risk with GSE261277 |
| GSE164649 Liverpool, LAIV vs TIV | 44 | carriage + **LAIV nasal IgA** | overlap risk; but see below |
| GSE124949 Liverpool EHPC | 20 | carriage | overlap; too small |
| GSE310241 pertussis CHIM, Dalhousie / La Jolla | **19** | symptomatic / asymptomatic / non-infected | too small; outcome not in GEO |
| GSE43523 allergic rhinitis | 12 | none (atopy only) | no outcome, no age |
| GSE205172 | — | — | mouse |

### The overlap finding, verified not assumed

GSE166161 is described as a separate RSV challenge series (26 samples, D0 and D7, 18 Infected /
8 Uninfected) and looked like a clean third cohort. Its `donor id` field uses the same scheme as
GSE155237 and **shares donors with it** — `AG_MH003` appears in both. It is the D0/D7 arm of the
same Imperial challenge programme whose D0/D3 arm is already in use. Excluded.

### GSE310241, the only genuinely new challenge cohort found

A *Bordetella pertussis* controlled human infection model (NCT05136599), North American, a third
pathogen and a different continent from both cohorts in use. Appealing in principle: participants
were primed in infancy with whole-cell or acellular pertussis vaccine, so it carries a vaccine
variable none of the others has.

It does not work as a primary test:

- of 414 samples, only **143 are nasopharyngeal swabs**, from **33 subjects**;
- only **36 pre-challenge nasal samples, from 19 subjects**;
- the outcome (symptomatic / asymptomatic / non-infected) and the vaccine-priming variable are
  **not in the GEO metadata** — they are in the preprint, Willemsen et al. 2026, medRxiv
  PMC13370606, [doi:10.64898/2026.06.26.26356067](https://doi.org/10.64898/2026.06.26.26356067),
  which is open access (retrieved via PubMed).

At 19 subjects split three ways it is below any useful power. It would be legitimate only as a
**pre-declared direction-only check**, never as a replication, and the sensible time to do that
is after a properly powered cohort exists, not instead of one.

## Second sweep: beyond GEO expression series

The first sweep was GEO-only and transcriptomic, which was too narrow. Extended to:

**PRIDE (proteomics).** Searched `nasal`, `nasopharyngeal`, `nasosorption`, `nasal lining
fluid`. **Nothing usable.** Human nasal proteomics in PRIDE is air–liquid-interface culture
(PXD066146, PXD038893), diesel/particle exposure (PXD065940), chronic rhinosinusitis
(PXD079333), nasopharyngeal carcinoma, or acute COVID swab diagnostics (PXD025277, PXD022889,
PXD029143). `nasosorption` and `nasal lining fluid` return **zero** studies. There is no
deposited nasal proteomic dataset with a baseline-then-outcome design.

This is worth noting because P5 is unusually suited to a proteomic test: MUC5AC, MUC5B, TFF1,
TFF3, FCGBP, AGR2 and ZG16B are all **secreted**, so the programme is in principle measurable in
nasal lining fluid. The data simply do not exist publicly.

**ArrayExpress / BioStudies.** Searched `nasal challenge`, `nasal carriage`, `nasosorption`,
`nasal vaccine human`. No human nasal baseline-plus-outcome study. The only challenge dataset,
E-MTAB-13041 (influenza H3N2 human challenge), is **whole blood**. `nasosorption` returns zero.

**Further GEO terms** (`colonization`, `carriage`, `prospective respiratory infection`,
`vaccination`, `household transmission`) surfaced two more candidates, both rejected:

| Dataset | Why not |
|---|---|
| GSE183071, n = 156 | nasal epithelium in **acute** COVID with severity outcome, not baseline→acquisition; targeted immune panel (GPL30569), so the mucin genes are probably not even measured |
| GSE150133, n = 71 | prospective RTI rates in infancy, but the transcriptome is **whole blood**, not nasal |

## A proteomic note on GSE253686

The NOVA study has **nasal Olink proteomics linked to shedding outcome** (its data files S9 and
S10 report proteins differing between shedders and non-shedders at D3/D5). So nasal proteomics
with a real outcome does exist — it is inside the same gated supplement as the transcriptomic
outcome table.

**But it would not test P5.** Olink inflammation panels measure cytokines and chemokines
(CXCL10, IL-6, CXCL11, TRAIL, PD-L1 and so on), not gel-forming mucins. The NOVA proteomics is
worth having for other reasons; it is not an orthogonal test of the mucin-barrier programme.

## Conclusion and recommendation

**There is no third cohort in GEO that meets the requirement.** Baseline-nasal-plus-outcome human
challenge designs are rare, and the ones that exist are small or already in use. Padding the
analysis with the Liverpool datasets would break the non-independence rule declared in
`PREREG_mvri_engagement.md` §0.1.

Two things are worth doing instead, in order:

1. **GSE253686 (NOVA), 75 participants with baseline SD −1 nasal arrays.** Already downloaded.
   Only the per-participant age and shedding outcome are missing, in supplementary data files
   S1 and S16 behind a PMC download gate — two browser clicks. This roughly doubles the evidence
   base and is the single highest-value unlock available. The plan must be locked before those
   files are opened.
2. ~~**Test P5 against the LAIV nasal IgA response in GSE164649.**~~ **CHECKED 2026-11-24 — not
   possible.** The GEO metadata carries only timepoint, inoculation date, virus species, sex,
   carriage, age and vaccine arm; the string "IgA" appears nowhere in it. The series has one
   supplementary file, the raw counts matrix. The paper (Carniel et al. 2021, *JCI Insight*,
   [doi:10.1172/jci.insight.141088](https://doi.org/10.1172/jci.insight.141088), PMC7934923,
   open access, retrieved via PubMed) has a single supplement — a 1.3 MB **PDF** — and its data
   availability statement says only "all other underlying data are provided in the manuscript",
   i.e. in figures. **Per-participant IgA titres are not public**, the same severed-link
   situation as the Gambia LAIV cohort. The original text of this recommendation follows, struck:

   ~~This is a *vaccine* endpoint,
   and the project currently has none — every outcome tested so far is pathogen challenge, which
   is exactly why the phrase "vaccine readiness" is still unlicensed. The Liverpool overlap
   weakens it as a replication of the *carriage* result but matters much less for a different
   endpoint. It would need its own locked plan, and whether per-participant IgA titres are public
   has not yet been checked.~~

   **Consequence: the project still has no vaccine endpoint of any kind, and no public route to
   one.** Every outcome available is pathogen challenge. The phrase "mucosal vaccine readiness"
   cannot be supported from public data as things stand.

Neither is a substitute for the third challenge cohort that does not exist. The honest position
for the manuscript remains the one in `results/mvri_outcome/PREREG_RESULTS.md` §5: two cohorts,
a likely inflated effect, and a direction that was already visible in prior work.
