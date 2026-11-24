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
2. **Test P5 against the LAIV nasal IgA response in GSE164649.** This is a *vaccine* endpoint,
   and the project currently has none — every outcome tested so far is pathogen challenge, which
   is exactly why the phrase "vaccine readiness" is still unlicensed. The Liverpool overlap
   weakens it as a replication of the *carriage* result but matters much less for a different
   endpoint. It would need its own locked plan, and whether per-participant IgA titres are public
   has not yet been checked.

Neither is a substitute for the third challenge cohort that does not exist. The honest position
for the manuscript remains the one in `results/mvri_outcome/PREREG_RESULTS.md` §5: two cohorts,
a likely inflated effect, and a direction that was already visible in prior work.
