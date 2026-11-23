# Module evidence

What each module is meant to read out, and the literature behind it. Generated from `src/nasage/resources/gene_modules.json` and `src/nasage/resources/module_evidence.csv`; edit those, not this file.

**How the references were checked.** Every PMID, DOI, year, journal and title below was confirmed against PubMed metadata on 2026-11-23. The *Checked* column says how far the claim itself was checked: against the abstract, the full text, or only the metadata (the title). Each module also lists the **evidence gaps** found: genes or claims that no verified reference supports. Treat those as open, not settled.

Method references cited elsewhere in the docs (Scanpy, muscat, scCODA, propeller, Splatter, powsimR, DESeq2, Benjamini-Hochberg and similar) are in `references.csv` with their original provenance; they were not re-checked in this pass.

## Contents

- [M01 Epithelial barrier and junction maintenance](#m01_epithelial_barrier)
- [M02 Motile ciliary (mucociliary) machinery](#m02_mucociliary_axoneme)
- [M02b Airway mucin / secretory product panel (small panel)](#m02b_airway_mucin_panel)
- [M03 Polymeric antibody transport (PIGR gene-level feature)](#m03_pigr_antibody_transport)
- [M03b FcRn IgG transport (gene-level feature)](#m03b_fcgrt_igg_transport)
- [M04 Innate pathogen sensing (constitutive receptors and proximal adaptors)](#m04_innate_sensing)
- [M05 Interferon-stimulated genes (type I and III interferon response)](#m05_interferon_stimulated_genes)
- [M06 MHC class I antigen processing and peptide loading](#m06_mhc_i_processing)
- [M07 MHC class II antigen processing and presentation](#m07_mhc_ii_presentation)
- [M08 Professional APC costimulation and maturation](#m08_apc_costimulation)
- [M09 Plasma-cell secretory machinery](#m09_plasma_secretory_machinery)
- [M10 IgA share of immunoglobulin heavy-chain transcripts (ratio feature)](#m10_iga_fraction)
- [M11 Senescence-associated cell-cycle arrest](#m11_senescence_arrest)
- [M12 Inflammatory cytokines (SASP / inflammaging core)](#m12_inflammatory_sasp)
- [M13 Tissue-resident memory T-cell core programme](#m13_trm_core)
- [General: cell types, QC and methods](#general)

<a id="m01_epithelial_barrier"></a>
## M01 Epithelial barrier and junction maintenance

**Genes (17):** CLDN1, CLDN3, CLDN4, CLDN7, OCLN, TJP1, TJP2, TJP3, CGN, F11R, CDH1, DSP, PKP1, PKP3, JUP, MARVELD2, CLDN23  
**Tested:** primary in basal, ciliated, secretory; exploratory in squamous  
**Reading for mucosal vaccination:** context-dependent

**What it measures.** Apical junctional complex: tight-junction strand proteins (claudins, occludin), cytoplasmic scaffolds (ZO-1/2/3, cingulin), and adherens/desmosomal adhesion molecules that set paracellular permeability of the nasal epithelium.

**How to read it.** A tighter barrier restricts pathogen and allergen penetration, but also restricts access of intranasally delivered antigen to subepithelial APCs; several mucosal adjuvants act by transiently loosening junctions. Higher score must not be read as 'better vaccine response'.

**Caveats.** Transcript abundance is a weak proxy for assembled junction function (no TEER equivalent); junction proteins are heavily post-translationally regulated. Claudin family members are functionally opposed (sealing vs pore-forming), so a single summed score is directionally ambiguous. TJP3 and MARVELD2 returned no Reactome mapping.

### Evidence

| Role | What it supports | Reference | Checked |
|---|---|---|---|
| Definition | Tight junctions form selective paracellular barriers in epithelia, mark the apical/basolateral border and connect to signalling networks and other adhesion complexes. | Zihni et al. 2016, *Nat Rev Mol Cell Biol* ([PMID 27353478](https://pubmed.ncbi.nlm.nih.gov/27353478/)) | abstract |
| Airway / nasal | In human sinonasal epithelium, CRS with nasal polyps shows lower trans-tissue resistance and lower occludin and ZO-1 (TJP1), and IFN-gamma and IL-4 reduce barrier integrity in primary nasal ALI cultures. | Soyka et al. 2012, *J Allergy Clin Immunol* ([PMID 22840853](https://pubmed.ncbi.nlm.nih.gov/22840853/)) | abstract |

### Evidence gaps

- No verified human reference links airway barrier gene expression to age.
- The desmosomal genes (DSP, PKP1, PKP3, JUP) are not specifically covered by the references below.

<a id="m02_mucociliary_axoneme"></a>
## M02 Motile ciliary (mucociliary) machinery

**Genes (20):** FOXJ1, TPPP3, CIMAP3, SPAG6, DNAH5, DNAI1, DNAI2, DNAH9, RSPH1, RSPH4A, RSPH9, CFAP43, CFAP53, HYDIN, SPEF2, TEKT1, CCDC39, CCDC40, DNALI1, IFT81  
**Tested:** primary in ciliated; exploratory in basal  
**Reading for mucosal vaccination:** context-dependent

**What it measures.** Axonemal structural and motility apparatus of multiciliated cells: outer/inner dynein arms, radial spokes, central apparatus, dynein assembly factors, plus the FOXJ1 multiciliogenesis transcriptional program.

**How to read it.** Efficient mucociliary clearance removes pathogens but equally removes intranasally delivered vaccine antigen, shortening mucosal residence time. Direction of benefit depends on formulation and residence requirements.

**Caveats.** Reactome maps only 2/20 of these genes; evidence rests on GO annotations and primary-ciliary-dyskinesia genetics, much of it IEA or model-organism derived. mRNA level does not measure ciliary beat frequency or coordinated transport.

### Evidence

| Role | What it supports | Reference | Checked |
|---|---|---|---|
| Definition | Foxj1 activates genes essential for motile cilia formation and is sufficient for ectopic motile-like cilia (zebrafish). | Yu et al. 2008, *Nat Genet* ([PMID 19011630](https://pubmed.ncbi.nlm.nih.gov/19011630/)) | abstract |
| Function | Mutations that impair motile-cilia beating cause motile ciliopathies including primary ciliary dyskinesia with impaired respiratory mucociliary clearance (the individual gene list is not in the abstract). | Wallmeier et al. 2020, *Nat Rev Dis Primers* ([PMID 32943623](https://pubmed.ncbi.nlm.nih.gov/32943623/)) | abstract |
| Airway / nasal | scRNA-seq of human nasal epithelial cultures traces basal-to-multiciliated differentiation via deuterosomal cells and finds hybrid goblet/multiciliated cells (so ciliary genes can appear in secretory-like cells). | Ruiz García et al. 2019, *Development* ([PMID 31558434](https://pubmed.ncbi.nlm.nih.gov/31558434/)) | abstract |
| Age | In healthy volunteers aged 11-90, nasal ciliary beat frequency fell and mucociliary clearance time and microtubular abnormalities rose with age (functional and ultrastructural data, not transcriptomic). | Ho et al. 2001, *Am J Respir Crit Care Med* ([PMID 11282777](https://pubmed.ncbi.nlm.nih.gov/11282777/)) | abstract |
| Caveat | During COVID-19 nasopharyngeal samples show loss of ciliated and expansion of secretory cells, so ciliary module scores at sample level partly reflect cell composition and infection. | Ziegler et al. 2021, *Cell* ([PMID 34352228](https://pubmed.ncbi.nlm.nih.gov/34352228/)) | abstract |

### Evidence gaps

- The FOXJ1 master-regulator evidence is from zebrafish.
- Individual axonemal genes were not checked against a published gene list.

<a id="m02b_airway_mucin_panel"></a>
## M02b Airway mucin / secretory product panel (small panel)

**Genes (5):** MUC5AC, MUC5B, SCGB1A1, BPIFA1, BPIFB1  
**Tested:** primary in secretory; exploratory in ciliated  
**Reading for mucosal vaccination:** context-dependent

**What it measures.** Gel-forming mucins and abundant secretory-cell products of the airway surface liquid.

**How to read it.** Mucus traps and clears antigen; it can also physically impede antigen uptake and vaccine delivery.

**Caveats.** Only five genes, with extreme dynamic range and strong secretory-cell restriction; reported as a small panel, not as an enrichment signature. SCGB1A1 has no Reactome mapping. Counts for MUC5AC/MUC5B are poorly captured by 3'-biased scRNA-seq in real data. In the nasal epithelium SCGB1A1 is a weaker secretory marker than in the bronchi (nasal secretory cells express MUC4, PI3 and SIX3 instead), so it contributes little in nasal brushings.

### Evidence

| Role | What it supports | Reference | Checked |
|---|---|---|---|
| Function | MUC5AC and MUC5B encode the principal airway mucins; Muc5b (not Muc5ac) is required for mucociliary clearance and airway defence in mice, and MUC5B at human airway surfaces can fall by up to 90% in disease. | Roy et al. 2014, *Nature* ([PMID 24317696](https://pubmed.ncbi.nlm.nih.gov/24317696/)) | abstract |
| Function | Review of airway mucus and mucin function and dysfunction in airway disease (no abstract in PubMed; content not checked). | Fahy et al. 2010, *N Engl J Med* ([PMID 21121836](https://pubmed.ncbi.nlm.nih.gov/21121836/)) | metadata only |
| Airway / nasal | upper airway epithelium contains distinct multiciliated, goblet, club and basal cells; goblet cells can be precursors of multiciliated cells, explaining hybrid goblet/ciliated cells. | Ruiz García et al. 2019, *Development* ([PMID 31558434](https://pubmed.ncbi.nlm.nih.gov/31558434/)) | abstract |
| Airway / nasal | scRNA-seq of human sinonasal tissue and nasal scrapings defined core, healthy, inflamed and polyp secretory-cell signatures, with secretory antimicrobial programmes shifting in type 2 inflammation. | Ordovas-Montanes et al. 2018, *Nature* ([PMID 30135581](https://pubmed.ncbi.nlm.nih.gov/30135581/)) | abstract |
| Airway / nasal | human lung and airway scRNA-seq atlas defining the expression profiles of 58 cell populations; the abstract does not itself state goblet-versus-club differences. | Travaglini et al. 2020, *Nature* ([PMID 33208946](https://pubmed.ncbi.nlm.nih.gov/33208946/)) | abstract |
| Age | In public bulk and single-cell lung data CFTR expression in goblet, club and basal cells fell with age (not ciliated cells); evidence that airway secretory-cell programmes change with age (not mucin genes themselves). | Corcoran et al. 2024, *Sci Rep* ([PMID 39572772](https://pubmed.ncbi.nlm.nih.gov/39572772/)) | abstract |
| Caveat | the same cell types (suprabasal, secretory, multiciliated) differ between nose (MUC4, PI3, SIX3) and tracheobronchial airways (SCGB1A1, TFF3), so SCGB1A1 is a weaker secretory marker in nasal brushings. | Deprez et al. 2020, *Am J Respir Crit Care Med* ([PMID 32726565](https://pubmed.ncbi.nlm.nih.gov/32726565/)) | abstract |

### Evidence gaps

- SCGB1A1 is a weak secretory marker in nasal (vs bronchial) epithelium.

<a id="m03_pigr_antibody_transport"></a>
## M03 Polymeric antibody transport (PIGR gene-level feature)

**Genes (1):** PIGR  
**Tested:** primary in secretory; exploratory in ciliated  
**Reading for mucosal vaccination:** supportive

**What it measures.** Epithelial receptor-mediated transcytosis of dimeric IgA/pentameric IgM to the mucosal surface and release of secretory component.

**How to read it.** Transport capacity is necessary for secretory IgA to reach the lumen, so higher PIGR plausibly supports mucosal antibody delivery. It reports transport capacity only, not antigen-specific antibody.

**Caveats.** Single gene: the 'score' is simply normalised PIGR expression and is reported as a gene-level feature, not an enrichment score. PIGR is induced by IFN-gamma, IL-17 and TNF, so its level also tracks inflammation. Airway-specific functional evidence cited here is from a pIgR-deficient mouse model (PMID 27046438); the human nasal evidence cited (PMID 41822483) is an association between nasal secretory IgA and viral shedding, not a measurement of PIGR transcript.

### Evidence

| Role | What it supports | Reference | Checked |
|---|---|---|---|
| Definition | secretory IgA is mainly dimeric IgA made by local plasma cells and transported by pIgR, the only receptor for epithelial IgA/IgM transport in pIgR-deficient mice; airways receive B cells primed in nasopharynx-associated lymphoid tissue. | Brandtzaeg et al. 2009, *Scand J Immunol* ([PMID 19906191](https://pubmed.ncbi.nlm.nih.gov/19906191/)) | abstract |
| Function | pIgR mediates epithelial transcytosis of polymeric IgA; its cleaved secretory component is released free or as part of SIgA; pIgR expression is regulated by TLR signalling and cytokines. | Kaetzel et al. 2005, *Immunol Rev* ([PMID 16048543](https://pubmed.ncbi.nlm.nih.gov/16048543/)) | abstract |
| Function | only J-chain-containing IgA/IgM polymers bind pIgR with high affinity, and SIgA and free secretory component come from cleavage of the pIgR extracellular domain. | Johansen et al. 2000, *Scand J Immunol* ([PMID 10972899](https://pubmed.ncbi.nlm.nih.gov/10972899/)) | abstract |
| Function | pIgR transports locally produced dimeric IgA across mucosal epithelia and is the precursor of secretory component, with expression regulated by host and microbial factors. | Johansen et al. 2011, *Mucosal Immunol* ([PMID 21956244](https://pubmed.ncbi.nlm.nih.gov/21956244/)) | abstract |
| Airway / nasal | In 143 people with acute SARS-CoV-2, an early rise in nasal-lavage secretory IgA (especially anti-spike NTD) was associated with shorter viral shedding (association only). | Brackett et al. 2026, *Front Immunol* ([PMID 41822483](https://pubmed.ncbi.nlm.nih.gov/41822483/)) | abstract |
| Airway / nasal | In human bronchi normal ciliated epithelium expresses pIgR with surface SIgA, remodelled epithelium has reduced pIgR and local SIgA deficiency, and full epithelial differentiation is needed for pIgR expression and IgA transcytosis. | Polosukhin et al. 2011, *Am J Respir Crit Care Med* ([PMID 21512171](https://pubmed.ncbi.nlm.nih.gov/21512171/)) | abstract |

<a id="m03b_fcgrt_igg_transport"></a>
## M03b FcRn IgG transport (gene-level feature)

**Genes (1):** FCGRT  
**Tested:** primary in none; exploratory in basal, ciliated, secretory, dendritic, macrophage  
**Reading for mucosal vaccination:** context-dependent

**What it measures.** Neonatal Fc receptor heavy chain mediating bidirectional IgG transcytosis and IgG/albumin recycling at epithelial surfaces.

**How to read it.** Relevant to IgG-mediated mucosal protection, which is a different effector arm from secretory IgA; grouping it with PIGR would merge two distinct transport systems.

**Caveats.** Single gene, reported separately from PIGR; beta-2 microglobulin (B2M) is the obligate partner chain but is also an MHC-I component and is scored in M06.

### Evidence

| Role | What it supports | Reference | Checked |
|---|---|---|---|
| Function | FcRn transports IgG (maternal-fetal transfer), protects IgG from degradation and is expressed at various sites in adults. | Roopenian et al. 2007, *Nat Rev Immunol* ([PMID 17703228](https://pubmed.ncbi.nlm.nih.gov/17703228/)) | abstract |
| Airway / nasal | Adult human, non-human primate and mouse bronchial epithelial cells express FcRn, and FcRn transports IgG-Fc across mouse respiratory epithelium in vivo. | Spiekermann et al. 2002, *J Exp Med* ([PMID 12163559](https://pubmed.ncbi.nlm.nih.gov/12163559/)) | abstract |

### Evidence gaps

- No verified age reference.

<a id="m04_innate_sensing"></a>
## M04 Innate pathogen sensing (constitutive receptors and proximal adaptors)

**Genes (15):** TLR2, TLR4, TLR5, TLR7, STING1, NOD1, NOD2, NLRP3, MAVS, TICAM1, MYD88, TBK1, IRF3, TRAF3, TLR9  
**Tested:** primary in basal, ciliated, secretory, dendritic, macrophage, Bcell; exploratory in Tcell  
**Reading for mucosal vaccination:** context-dependent

**What it measures.** Pattern-recognition receptors (TLRs, STING, NLRs) and proximal adaptors/kinases up to IRF3. Sensors that are themselves interferon-stimulated genes (RIGI, IFIH1, DHX58, TLR3, AIM2, CGAS) are excluded so that the module is less of a proxy for interferon tone.

**How to read it.** Sensing capacity may support responsiveness to TLR/STING-agonist adjuvants, but also drives reactogenicity; it is not intrinsically beneficial.

**Caveats.** Removing the canonical ISG sensors reduces but does not remove interferon dependence; check the donor-level correlation with M05 (nas.tl.module_correlation) before reading M04 as independent of ISG tone. TLR7 is expressed mainly by immune cells and is barely detected in airway epithelial cells, so in epithelium the module is carried by TLR2/4/5, NOD1/2, NLRP3 and the adaptors. MYD88, TBK1, TRAF3 and MAVS are broadly expressed. Expression of a receptor does not imply signalling competence.

### Evidence

| Role | What it supports | Reference | Checked |
|---|---|---|---|
| Definition | Innate sensing is mediated by PRR families (TLRs, RIG-I-like receptors, NOD-like receptors, C-type lectins) that signal to inflammatory gene expression. | Takeuchi et al. 2010, *Cell* ([PMID 20303872](https://pubmed.ncbi.nlm.nih.gov/20303872/)) | abstract |
| Function | RIG-I and MDA5 sense viral RNA and signal via their CARDs to induce type I IFN, whereas LGP2 lacks CARDs and acts as a negative regulator (relevant to the excluded RLR genes). | Yoneyama et al. 2005, *J Immunol* ([PMID 16116171](https://pubmed.ncbi.nlm.nih.gov/16116171/)) | abstract |
| Function | The cGAS-STING pathway detects intracellular DNA from pathogens and from self (DAMP) and mediates inflammatory responses, including sterile inflammation. | Ablasser et al. 2019, *Science* ([PMID 30846571](https://pubmed.ncbi.nlm.nih.gov/30846571/)) | abstract |
| Airway / nasal | Airway epithelium expresses many innate immune receptors including TLRs and NOD-like receptors, which drive cytokine/chemokine production and recruit neutrophils and T cells. | Parker et al. 2011, *Am J Respir Cell Mol Biol* ([PMID 21330463](https://pubmed.ncbi.nlm.nih.gov/21330463/)) | abstract |
| Age | Claims 1 and 3: upper-airway scRNA-seq (4 weeks to 77 years) showed higher basal expression of the PRRs MDA5 (IFIH1) and RIG-I (DDX58) in children's epithelial cells, macrophages and DCs than in adults. | Loske et al. 2022, *Nat Biotechnol* ([PMID 34408314](https://pubmed.ncbi.nlm.nih.gov/34408314/)) | abstract |
| Age | With age innate immunity becomes dysregulated: basal inflammation rises (partly from endogenous DAMPs acting on PRRs) while PRR signalling and cell migration needed for pathogen and vaccine responses are impaired. | Shaw et al. 2013, *Nat Rev Immunol* ([PMID 24157572](https://pubmed.ncbi.nlm.nih.gov/24157572/)) | abstract |
| Caveat | TLR7 and TLR8 were not expressed in human A549 epithelial cells, so endosomal TLR7 in M04 will mostly reflect immune cells rather than epithelium (cell-line data). | Tissari et al. 2005, *J Immunol* ([PMID 15778392](https://pubmed.ncbi.nlm.nih.gov/15778392/)) | abstract |
| Caveat | C6orf150/MB21D1 (CGAS) was in the ISG library and was a broadly acting antiviral ISG, so CGAS is interferon-inducible and not purely a constitutive sensor. | Schoggins et al. 2011, *Nature* ([PMID 21478870](https://pubmed.ncbi.nlm.nih.gov/21478870/)) | full text |

### Evidence gaps

- TLR4, TLR5 and TLR9 have no airway-epithelium-specific reference here.
- TLR7 is barely expressed by airway epithelial cells (cell-line evidence).

<a id="m05_interferon_stimulated_genes"></a>
## M05 Interferon-stimulated genes (type I and III interferon response)

**Genes (22):** ISG15, IFI6, IFI44L, IFIT1, IFIT2, IFIT3, MX1, MX2, OAS1, OAS2, OAS3, OASL, RSAD2, STAT1, STAT2, IRF7, IRF9, BST2, XAF1, HERC5, USP18, IFITM3  
**Tested:** primary in basal, ciliated, secretory, dendritic, macrophage, Tcell, Bcell, plasma; exploratory in neutrophil, NK  
**Reading for mucosal vaccination:** context-dependent

**What it measures.** Canonical interferon-stimulated effector and amplifier genes downstream of IFNAR and IFNLR1 JAK-STAT signalling. At the airway epithelium type III interferon (IFN-lambda) induces largely the same genes as type I, so the module cannot tell the two apart.

**How to read it.** Higher baseline ISG tone has been reported in children's upper airways and is antiviral, but type I IFN is also immunopathological in other settings and can restrict replication of live-attenuated vaccines. A high score must not be labelled automatically beneficial.

**Caveats.** Not specific to type I interferon: type III (dominant at the airway epithelium) and partly type II interferon induce the same genes. Baseline ISG tone largely reflects recent or ongoing respiratory viral exposure rather than a stable trait; without viral testing it is unattributable.

### Evidence

| Role | What it supports | Reference | Checked |
|---|---|---|---|
| Definition | ISG products enhance pathogen sensing, negatively regulate JAK-STAT signalling and directly inhibit viruses (review of ISG biology). | Schneider et al. 2014, *Annu Rev Immunol* ([PMID 24555472](https://pubmed.ncbi.nlm.nih.gov/24555472/)) | abstract |
| Definition | AIM2 is an IFN-inducible cytosolic DNA sensor; IFN-beta treatment induced AIM2 protein in human fibroblasts, and AIM2 rose in senescent fibroblasts with IL-1beta. | Duan et al. 2011, *Mol Cancer Res* ([PMID 21471287](https://pubmed.ncbi.nlm.nih.gov/21471287/)) | abstract |
| Definition | LGP2 (DHX58) mRNA is induced by interferon, dsRNA or Sendai virus and acts as a feedback inhibitor of RIG-I/MAVS signalling. | Komuro et al. 2006, *J Virol* ([PMID 17020950](https://pubmed.ncbi.nlm.nih.gov/17020950/)) | abstract |
| Definition | IFN-alpha upregulates TLR3 expression in human A549 epithelial cells and HUVECs, i.e. TLR3 is interferon-inducible. | Tissari et al. 2005, *J Immunol* ([PMID 15778392](https://pubmed.ncbi.nlm.nih.gov/15778392/)) | abstract |
| Definition | overexpression screen of >380 human ISGs; the abstract names RIG-I (DDX58), MDA5 (IFIH1), C6orf150/MB21D1 (CGAS), IFI6, IFI44L, OASL and IFITM3 among ISG effectors. | Schoggins et al. 2011, *Nature* ([PMID 21478870](https://pubmed.ncbi.nlm.nih.gov/21478870/)) | abstract |
| Definition | Interferome v2.0 is a curated database of genes regulated by type I, II and III IFN in human and mouse cells and tissues, usable to check the ISG status of individual genes; the abstract does not list specific genes. | Rusinova et al. 2013, *Nucleic Acids Res* ([PMID 23203888](https://pubmed.ncbi.nlm.nih.gov/23203888/)) | abstract |
| Function | type I and type III IFNs activate shared signalling pathways and transcriptional programs, with type III IFN acting as a front-line defence at epithelial barriers while limiting inflammation. | Lazear et al. 2019, *Immunity* ([PMID 30995506](https://pubmed.ncbi.nlm.nih.gov/30995506/)) | abstract |
| Airway / nasal | in mouse influenza infection IFN-lambdas are the first IFNs produced and act at the respiratory epithelial barrier to suppress initial spread without inflammation, with type I IFN acting later (mouse data). | Galani et al. 2017, *Immunity* ([PMID 28514692](https://pubmed.ncbi.nlm.nih.gov/28514692/)) | abstract |
| Airway / nasal | rhinovirus infection of human bronchial epithelial cells induced IFN-beta, IFN-lambda1/2/3, IRF7, RIG-I and MDA5 mRNA in a TLR3/TRIF-, MDA5- and IRF3-dependent way. | Wang et al. 2009, *J Immunol* ([PMID 19890046](https://pubmed.ncbi.nlm.nih.gov/19890046/)) | abstract |
| Airway / nasal | after intranasal influenza A infection in mice type III IFN levels greatly exceed type I in the lung, and the IFN-lambda receptor is expressed mainly by epithelial cells (mouse data). | Jewell et al. 2010, *J Virol* ([PMID 20739515](https://pubmed.ncbi.nlm.nih.gov/20739515/)) | abstract |
| Airway / nasal | scRNA-seq of nasopharyngeal swabs from 58 people: mild/moderate COVID-19 epithelial cells express antiviral/IFN-responsive genes while severe cases have muted responses; infection also changes epithelial composition (secretory up, ciliated down). | Ziegler et al. 2021, *Cell* ([PMID 34352228](https://pubmed.ncbi.nlm.nih.gov/34352228/)) | abstract |
| Airway / nasal | In serial human nasopharyngeal samples SARS-CoV-2 induced ISGs in the upper respiratory tract, and prior rhinovirus infection accelerated ISG responses and blocked SARS-CoV-2 replication, so baseline nasal ISG state depends on recent infections. | Cheemarla et al. 2021, *J Exp Med* ([PMID 34128960](https://pubmed.ncbi.nlm.nih.gov/34128960/)) | abstract |
| Age | children's upper airway cells are primed for virus sensing, with higher basal IFIH1/DDX58 and stronger early innate antiviral responses than adults. | Loske et al. 2022, *Nat Biotechnol* ([PMID 34408314](https://pubmed.ncbi.nlm.nih.gov/34408314/)) | abstract |
| Age | Older people with severe COVID-19 showed type I and II IFN deficiencies and reduced antiviral defence genes (LY6E, IFNAR1) in nasal secretory cells (integrative re-analysis of public data). | Hou et al. 2022, *Aging Cell* ([PMID 35023286](https://pubmed.ncbi.nlm.nih.gov/35023286/)) | abstract |
| Age | single-cell profiling of nasal, tracheal and bronchial samples showed that airway cells of healthy children are already in an interferon-activated state compared with adults. | Yoshida et al. 2022, *Nature* ([PMID 34937051](https://pubmed.ncbi.nlm.nih.gov/34937051/)) | abstract |
| Caveat | Type I IFNs are protective against viruses but can cause immunopathology in acute infections and immunosuppression in chronic infections, so a high ISG score is not simply 'better immunity'. | McNab et al. 2015, *Nat Rev Immunol* ([PMID 25614319](https://pubmed.ncbi.nlm.nih.gov/25614319/)) | abstract |
| Caveat | Full text: the compiled library of 389 ISGs included PSMB8, PSMB9, MX1, IFIT3, IFI6, IFITM3, BST2, IRF7, RSAD2 and also CDKN1A (toxic on overexpression), so CDKN1A is IFN-inducible and overlaps the ISG and senescence modules. | Schoggins et al. 2011, *Nature* ([PMID 21478870](https://pubmed.ncbi.nlm.nih.gov/21478870/)) | full text |

### Evidence gaps

- That type III interferon dominates at the airway epithelium rests largely on mouse data.

<a id="m06_mhc_i_processing"></a>
## M06 MHC class I antigen processing and peptide loading

**Genes (15):** PSMB8, PSMB9, PSMB10, TAP1, TAP2, TAPBP, B2M, HLA-A, HLA-B, HLA-C, ERAP1, ERAP2, CALR, CANX, NLRC5  
**Tested:** primary in basal, ciliated, secretory, dendritic, macrophage, Tcell, Bcell, plasma; exploratory in NK  
**Reading for mucosal vaccination:** context-dependent

**What it measures.** Immunoproteasome subunits, peptide transport into the ER, the peptide-loading complex, and MHC-I heavy/light chains.

**How to read it.** Relevant to CD8 T-cell priming/recall, but PSMB8, PSMB9, TAP1 and B2M are themselves interferon-inducible, so a high score may simply reflect interferon tone (M05) rather than independent presentation capacity.

**Caveats.** Strong overlap with the interferon module by biology, not by gene sharing; HLA-A/B/C mapping in real data is confounded by allelic reference bias (not reproduced in this simulation). Immunoproteasome, TAP and NLRC5 are interferon-inducible, so M06 tracks M05; check nas.tl.module_correlation. HLA-A/B/C are polymorphic and align to the reference genome with allele-dependent efficiency, so scores partly reflect HLA type and ancestry; rerun with nas.drop_hla() as a sensitivity analysis. HLA-A/B/C and B2M are among the most highly expressed genes in most cells, so UCell scores can sit near their ceiling (see the score_distribution table).

### Evidence

| Role | What it supports | Reference | Checked |
|---|---|---|---|
| Definition | Review of the molecular machinery of MHC class I antigen processing and presentation (proteasome, peptide transport and loading). | Neefjes et al. 2011, *Nat Rev Immunol* ([PMID 22076556](https://pubmed.ncbi.nlm.nih.gov/22076556/)) | abstract |
| Definition | review establishing NLRC5 as the MHC class I transactivator (CITA), the class I counterpart of CIITA. | Kobayashi et al. 2012, *Nat Rev Immunol* ([PMID 23175229](https://pubmed.ncbi.nlm.nih.gov/23175229/)) | abstract |
| Function | full text lists PSMB8 and PSMB9 among the ISGs tested in the type I IFN effector screen. | Schoggins et al. 2011, *Nature* ([PMID 21478870](https://pubmed.ncbi.nlm.nih.gov/21478870/)) | full text |
| Function | NLRC5 is IFN-gamma-inducible, binds and transactivates MHC-I promoters, is required for IFN-gamma-induced MHC-I upregulation, and also induces B2M, TAP and LMP (PSMB) genes. | Meissner et al. 2010, *Proc Natl Acad Sci U S A* ([PMID 20639463](https://pubmed.ncbi.nlm.nih.gov/20639463/)) | abstract |
| Age | Healthy children's airway cells are in an interferon-activated state; because MHC-I processing genes are IFN-inducible (claim 4), age differences in M06 may partly track the IFN state (inference). | Yoshida et al. 2022, *Nature* ([PMID 34937051](https://pubmed.ncbi.nlm.nih.gov/34937051/)) | abstract |
| Caveat | reference mapping of short reads at HLA-A/-B/-C/-DRB1/-DQB1 gave 18.6% incorrect genotype calls with bias towards reference alleles (DNA-seq; the implication for HLA expression quantification in scRNA-seq is an extrapolation). | Brandt et al. 2015, *G3 (Bethesda)* ([PMID 25787242](https://pubmed.ncbi.nlm.nih.gov/25787242/)) | abstract |

### Evidence gaps

- The HLA mapping-bias reference is from DNA sequencing; applying it to scRNA-seq expression is an extrapolation.
- The age row is an inference from the ISG age difference, not a direct MHC-I measurement.

<a id="m07_mhc_ii_presentation"></a>
## M07 MHC class II antigen processing and presentation

**Genes (12):** CD74, HLA-DRA, HLA-DRB1, HLA-DPA1, HLA-DPB1, HLA-DQA1, HLA-DQB1, HLA-DMA, HLA-DMB, CIITA, CTSS, LGMN  
**Tested:** primary in dendritic, macrophage, Bcell; exploratory in basal, ciliated, secretory, Tcell  
**Reading for mucosal vaccination:** supportive in professional APC; context-dependent in epithelium

**What it measures.** Class II alpha/beta chains, invariant chain, HLA-DM editing, CIITA transactivation and the endolysosomal proteases that generate class II peptides.

**How to read it.** In dendritic cells, macrophages and B cells this reports classical CD4 priming machinery. Airway and nasal epithelial cells can express HLA-DR and CD86 under inflammatory/IFN-gamma conditions, but epithelial class II expression is not evidence of professional APC function and is scored and interpreted separately.

**Caveats.** CIITA and CTSS are expressed well below the class II structural genes, so the score is dominated by HLA-D* and CD74 abundance. Epithelial class II is inducible and cannot be distinguished from professional APC function by transcript level alone. HLA class II genes are highly polymorphic; reference-alignment bias makes scores partly reflect HLA type and ancestry. Rerun with nas.drop_hla() as a sensitivity analysis. LGMN has no verified supporting reference in references.csv.

### Evidence

| Role | What it supports | Reference | Checked |
|---|---|---|---|
| Definition | MHC-II expression requires RFX (RFXANK, RFX5, RFXAP) and the master coactivator CIITA, whose defects cause bare lymphocyte syndrome. | Reith et al. 2001, *Annu Rev Immunol* ([PMID 11244040](https://pubmed.ncbi.nlm.nih.gov/11244040/)) | abstract |
| Definition | Review of the molecular machinery of MHC class II antigen presentation. | Neefjes et al. 2011, *Nat Rev Immunol* ([PMID 22076556](https://pubmed.ncbi.nlm.nih.gov/22076556/)) | abstract |
| Function | Cathepsin S (CTSS) and cathepsin L degrade invariant chain (CD74) and process antigen in the MHC-II pathway; the abstract does not discuss legumain (LGMN). | Hsing et al. 2005, *Immunol Rev* ([PMID 16181340](https://pubmed.ncbi.nlm.nih.gov/16181340/)) | abstract |
| Function | CIITA is induced by IFN-gamma, is required for MHC-II induction and is sufficient to activate MHC class II genes, i.e. it regulates both inducible and constitutive MHC-II. | Steimle et al. 1994, *Science* ([PMID 8016643](https://pubmed.ncbi.nlm.nih.gov/8016643/)) | abstract |
| Airway / nasal | Nasal epithelial cells constitutively express HLA-DR, IFN-gamma upregulates it, and cultured nasal epithelial cells can drive antigen-specific autologous T-cell proliferation. | Takizawa et al. 2007, *Clin Exp Allergy* ([PMID 17359392](https://pubmed.ncbi.nlm.nih.gov/17359392/)) | abstract |
| Airway / nasal | MHC class II is highly expressed on lung and intestinal epithelial cells, which may act as non-professional antigen-presenting cells. | Wosen et al. 2018, *Front Immunol* ([PMID 30319613](https://pubmed.ncbi.nlm.nih.gov/30319613/)) | abstract |
| Airway / nasal | Human nasal polyp epithelium stains for HLA-DR in situ, and IFN-gamma increases HLA-DR on cultured nasal epithelial cells. | Papon et al. 2002, *Laryngoscope* ([PMID 12439183](https://pubmed.ncbi.nlm.nih.gov/12439183/)) | abstract |
| Caveat | reference mapping bias at HLA-DRB1 and HLA-DQB1 (DNA-seq), so per-gene HLA class II counts from reference-aligned scRNA-seq may be biased between donors with different alleles. | Brandt et al. 2015, *G3 (Bethesda)* ([PMID 25787242](https://pubmed.ncbi.nlm.nih.gov/25787242/)) | abstract |

### Evidence gaps

- LGMN has no verified supporting reference.

<a id="m08_apc_costimulation"></a>
## M08 Professional APC costimulation and maturation

**Genes (12):** CD80, CD86, CD40, CD83, ICOSLG, TNFSF9, CD70, TNFSF4, RELB, NFKB2, CCR7, LAMP3  
**Tested:** primary in dendritic, macrophage, Bcell  
**Reading for mucosal vaccination:** supportive

**What it measures.** B7-family and TNF-superfamily costimulatory ligands plus the non-canonical NF-kB maturation axis and DC maturation markers.

**How to read it.** Costimulatory ligand availability is a limiting factor for naive T-cell priming, so higher levels are plausibly supportive; however, high baseline costimulation may equally mark ongoing local inflammation.

**Caveats.** Mixed content: CD83, CCR7, LAMP3, RELB and NFKB2 are maturation/migration markers rather than costimulatory ligands, so the module conflates ligand availability with maturation state. CD83, ICOSLG, CCR7 and LAMP3 returned no mapping to the two anchor pathways. Restricted to professional APC populations. CCR7, LAMP3, CD80, CD86 and CD40 together mark the mature/regulatory DC state, so in a dendritic population that pools cDC1, cDC2, pDC and mature DCs the score largely measures the mature-DC fraction. Annotate DC subtypes and keep them as subtypes so the analysis can adjust for subtype fractions.

### Evidence

| Role | What it supports | Reference | Checked |
|---|---|---|---|
| Definition | scRNA-seq in human and mouse lung cancer defined mregDCs co-expressing maturation genes (CD40, CCR7, IL12B) and immunoregulatory genes (CD274, PDCD1LG2, CD200), acquired by DC1 and DC2 after antigen uptake. | Maier et al. 2020, *Nature* ([PMID 32269339](https://pubmed.ncbi.nlm.nih.gov/32269339/)) | abstract |
| Age | Innate immune cell functions such as cell migration and PRR signalling are impaired in aged individuals (review, not DC- or airway-specific). | Shaw et al. 2013, *Nat Rev Immunol* ([PMID 24157572](https://pubmed.ncbi.nlm.nih.gov/24157572/)) | abstract |
| Caveat | Human nasal epithelial cells constitutively express HLA-DR and CD86 but not CD80, and CD86 is upregulated by IL-1beta/TNF-alpha, so M08 scores in epithelium do not mean DC maturation. | Takizawa et al. 2007, *Clin Exp Allergy* ([PMID 17359392](https://pubmed.ncbi.nlm.nih.gov/17359392/)) | abstract |

### Evidence gaps

- CD83, LAMP3, ICOSLG, TNFSF9, CD70 and TNFSF4 are not named in the abstracts checked.
- The age reference concerns innate immunity generally, not DCs specifically.

<a id="m09_plasma_secretory_machinery"></a>
## M09 Plasma-cell secretory machinery

**Genes (9):** PRDM1, XBP1, IRF4, DERL3, MZB1, SEC11C, SSR4, FKBP11, TNFRSF17  
**Tested:** primary in plasma; exploratory in Bcell  
**Reading for mucosal vaccination:** supportive

**What it measures.** Plasma-cell differentiation regulators (PRDM1, XBP1, IRF4) and the ER/secretory apparatus that supports high-rate antibody secretion (MZB1, DERL3, SEC11C, SSR4, FKBP11), plus the plasma-cell survival receptor BCMA (TNFRSF17). Immunoglobulin constant-region genes and JCHAIN are excluded: isotype is scored separately (M10), and Ig transcripts are the commonest ambient contaminant in nasal samples.

**How to read it.** JCHAIN is required for polymeric IgA assembly and therefore couples this module mechanistically to PIGR-mediated transport; both are needed for secretory IgA.

**Caveats.** Measures the secretory programme of the plasma cells recovered, not antibody output or specificity. Plasma cells are rare in brushings, so few donors may reach the minimum cell count. Within plasma cells these genes are very highly expressed, so check the score_distribution table for ceiling effects.

### Evidence

| Role | What it supports | Reference | Checked |
|---|---|---|---|
| Definition | Review: plasmablast and plasma-cell differentiation from B cells is governed by a central gene-regulatory network (the abstract does not name individual genes). | Nutt et al. 2015, *Nat Rev Immunol* ([PMID 25698678](https://pubmed.ncbi.nlm.nih.gov/25698678/)) | abstract |
| Function | BCMA (TNFRSF17) is expressed more highly in plasma cells than in resting B cells and is required for survival of long-lived bone-marrow plasma cells (mouse). | O'Connor et al. 2004, *J Exp Med* ([PMID 14707116](https://pubmed.ncbi.nlm.nih.gov/14707116/)) | abstract |
| Function | Blimp-1 (PRDM1) is required to induce most plasma-cell genes including XBP1, and XBP1 induces secretory-pathway genes and expands the ER (mouse). | Shaffer et al. 2004, *Immunity* ([PMID 15345222](https://pubmed.ncbi.nlm.nih.gov/15345222/)) | abstract |
| Function | B-cell-specific loss of Blimp-1 (PRDM1) abolishes Ig-secreting plasma cells and XBP-1 induction, and XBP-1 cannot rescue plasmacytic differentiation without Blimp-1 (mouse). | Shapiro-Shelef et al. 2003, *Immunity* ([PMID 14563324](https://pubmed.ncbi.nlm.nih.gov/14563324/)) | abstract |
| Airway / nasal | in upper airway mucosa memory B cells differentiate locally into Ig-producing plasma cells whose product is mainly J-chain-containing polymeric IgA, transported by pIgR/secretory component. | Brandtzaeg et al. 1996, *Acta Otolaryngol* ([PMID 8725503](https://pubmed.ncbi.nlm.nih.gov/8725503/)) | abstract |

### Evidence gaps

- No gene-specific reference for MZB1, DERL3, SEC11C, SSR4, FKBP11 or IRF4.
- The XBP1, PRDM1 and BCMA references are mouse studies.
- No verified age reference.

<a id="m10_iga_fraction"></a>
## M10 IgA share of immunoglobulin heavy-chain transcripts (ratio feature)

**Genes (7):** IGHA1, IGHA2, IGHG1, IGHG2, IGHG3, IGHG4, IGHM  
**Tested:** primary in plasma; exploratory in Bcell  
**Reading for mucosal vaccination:** supportive

**What it measures.** Donor-level fraction of IgA heavy-chain counts (IGHA1+IGHA2) among all IgA, IgG and IgM heavy-chain counts in the cell population, reported on the logit scale. Not a UCell score.

**How to read it.** Secretory IgA is the main antibody at the nasal mucosal surface; a shift from IgA towards IgG among local plasma cells changes the type of antibody available there. It does not measure antigen-specific antibody.

**Caveats.** Computed from pooled counts, so it is dominated by the cells with the most Ig transcript. Ambient Ig contamination comes mostly from plasma cells themselves, which limits but does not remove bias. IGHD and IGHE are not included.

### Evidence

| Role | What it supports | Reference | Checked |
|---|---|---|---|
| Definition | Human IgA1 and IgA2 switch recombination was profiled in nasal mucosa and other compartments; local IgA1/IgA2 proportions parallel IgA1+/IgA2+ plasma-cell distributions, supporting summing IGHA1 and IGHA2. | Lin et al. 2014, *Mucosal Immunol* ([PMID 24064668](https://pubmed.ncbi.nlm.nih.gov/24064668/)) | abstract |
| Function | J chain made by mucosal and glandular plasma cells is needed to form polymeric IgA/IgM with high affinity for pIgR, so JCHAIN plus pIgR are key to secretory IgA. | Johansen et al. 2000, *Scand J Immunol* ([PMID 10972899](https://pubmed.ncbi.nlm.nih.gov/10972899/)) | abstract |
| Airway / nasal | locally produced Ig in upper-airway mucosa consists mainly of polymeric IgA, while IgG reaches secretions largely by passive diffusion. | Brandtzaeg et al. 1996, *Acta Otolaryngol* ([PMID 8725503](https://pubmed.ncbi.nlm.nih.gov/8725503/)) | abstract |

### Evidence gaps

- No verified age reference.

<a id="m11_senescence_arrest"></a>
## M11 Senescence-associated cell-cycle arrest

**Genes (5):** CDKN1A, CDKN2A, SERPINE1, IGFBP7, GDF15  
**Tested:** primary in basal, ciliated, secretory; exploratory in macrophage, Tcell  
**Reading for mucosal vaccination:** context-dependent

**What it measures.** Cyclin-dependent kinase inhibitors that enforce senescence-associated arrest (p21/CDKN1A, p16/CDKN2A) plus secreted factors reported in senescent cells and ageing plasma (SERPINE1/PAI-1, IGFBP7, GDF15).

**How to read it.** Senescent epithelial and immune cells accumulate with age and may alter mucosal repair and responses. A transcript score is not a diagnosis of senescence.

**Caveats.** No single transcript defines senescence. CDKN1A is also induced by p53 activation, DNA damage and interferon (it is in the interferon-stimulated gene library), and CDKN2A mRNA is low and unreliably detected in scRNA-seq; for these reasons the SenMayo gene set deliberately excludes both. p16- and p21-positive cells were not found to increase with age in human lung. Read the module as an arrest-associated programme, not as a count of senescent cells.

### Evidence

| Role | What it supports | Reference | Checked |
|---|---|---|---|
| Definition | consensus definition of cellular senescence and recommended biomarkers (the abstract does not list individual genes). | Gorgoulis et al. 2019, *Cell* ([PMID 31675495](https://pubmed.ncbi.nlm.nih.gov/31675495/)) | abstract |
| Definition | The SASP Atlas found GDF15 and SERPINs among proteins secreted by senescent cells that also correlate with age in human plasma (the abstract says 'SERPINs', not SERPINE1 by name). | Basisty et al. 2020, *PLoS Biol* ([PMID 31945054](https://pubmed.ncbi.nlm.nih.gov/31945054/)) | abstract |
| Function | Secreted IGFBP7 is needed for BRAF-V600E oncogene-induced senescence in human fibroblasts and melanocytes. | Wajapeyee et al. 2008, *Cell* ([PMID 18267069](https://pubmed.ncbi.nlm.nih.gov/18267069/)) | abstract |
| Age | There is no universal senescence marker; p16 and p21 positive cells rose with age in several human organs but not in lung, so the age trend is tissue-specific. | Idda et al. 2020, *Aging (Albany NY)* ([PMID 32160592](https://pubmed.ncbi.nlm.nih.gov/32160592/)) | abstract |
| Caveat | International Cell Senescence Association consensus that accurate in vivo detection of senescent cells needs several markers combined (it proposes a multi-marker algorithm), so no single marker is sufficient. | Gorgoulis et al. 2019, *Cell* ([PMID 31675495](https://pubmed.ncbi.nlm.nih.gov/31675495/)) | abstract |
| Caveat | Senescence markers (CDKN2A, CDKN1A, SASP, LMNB1 loss) are non-exclusive and the senescent transcriptome varies by inducer, cell type and time since induction. | Hernandez-Segura et al. 2017, *Curr Biol* ([PMID 28844647](https://pubmed.ncbi.nlm.nih.gov/28844647/)) | abstract |
| Caveat | full text states CDKN2A and CDKN1A were intentionally left out of the 125-gene SenMayo set (used to validate it), and CDKN2A mRNA is low and may not be reliably detected in scRNA-seq; IGFBP7 is a SenMayo SASP gene; SERPINE1 and GDF15 membership could not be confirmed (Supplementary Data 1 not accessible). | Saul et al. 2022, *Nat Commun* ([PMID 35974106](https://pubmed.ncbi.nlm.nih.gov/35974106/)) | full text |

### Evidence gaps

- SenMayo deliberately excludes CDKN1A and CDKN2A; their inclusion here follows the consensus review, not SenMayo.
- SERPINE1 and GDF15 membership of SenMayo could not be checked (supplementary list not accessible); GDF15 and SERPINs are supported by the ageing-plasma SASP atlas.

<a id="m12_inflammatory_sasp"></a>
## M12 Inflammatory cytokines (SASP / inflammaging core)

**Genes (7):** IL6, IL1A, IL1B, CXCL8, CXCL1, CCL2, TNF  
**Tested:** primary in basal, ciliated, secretory, macrophage; exploratory in dendritic, neutrophil  
**Reading for mucosal vaccination:** context-dependent

**What it measures.** Pro-inflammatory cytokines and chemokines that form the core of the senescence-associated secretory phenotype and of age-associated chronic low-grade inflammation (inflammaging).

**How to read it.** Chronic low-grade inflammation is a hallmark of ageing, but the same cytokines rise acutely with infection, allergy, smoking and sampling trauma.

**Caveats.** Cytokine transcripts are sparse and bursty in scRNA-seq, so the per-donor mean is used (aggregation 'auto'). Acute infection and sampling stress are strong confounders; brushing itself can induce immediate-early inflammatory genes.

### Evidence

| Role | What it supports | Reference | Checked |
|---|---|---|---|
| Definition | defined the SASP: senescent fibroblasts and epithelial cells induced by genotoxic stress secrete many inflammatory factors, with IL-6 and IL-8 (CXCL8) as key paracrine mediators. | Coppé et al. 2008, *PLoS Biol* ([PMID 19053174](https://pubmed.ncbi.nlm.nih.gov/19053174/)) | abstract |
| Definition | Full text: IL1A, CXCL8 and CCL2 are highly connected nodes of the SenMayo gene network, and CDKN1A, CCL2 and IL6 rose with age in two human bone cohorts. | Saul et al. 2022, *Nat Commun* ([PMID 35974106](https://pubmed.ncbi.nlm.nih.gov/35974106/)) | full text |
| Age | inflammaging is the chronic, sterile, low-grade inflammation that develops with ageing and contributes to age-related disease. | Franceschi et al. 2018, *Nat Rev Endocrinol* ([PMID 30046148](https://pubmed.ncbi.nlm.nih.gov/30046148/)) | abstract |
| Caveat | Nasopharyngeal epithelial states change strongly with acute viral infection, so current infection status is a confounder for inflammatory and IFN module scores. | Ziegler et al. 2021, *Cell* ([PMID 34352228](https://pubmed.ncbi.nlm.nih.gov/34352228/)) | abstract |
| Caveat | Systemic chronic inflammation is promoted by infections, inactivity, diet, toxicants and stress and underlies many diseases, so inflammatory cytokine expression is not specific to ageing. | Furman et al. 2019, *Nat Med* ([PMID 31806905](https://pubmed.ncbi.nlm.nih.gov/31806905/)) | abstract |

### Evidence gaps

- These cytokines are not specific to ageing; infection is a strong confounder.

<a id="m13_trm_core"></a>
## M13 Tissue-resident memory T-cell core programme

**Genes (7):** CD69, ITGAE, ITGA1, CXCR6, PDCD1, CRTAM, DUSP6  
**Tested:** primary in Tcell; exploratory in NK  
**Reading for mucosal vaccination:** supportive

**What it measures.** Genes of the human tissue-resident memory T-cell (TRM) core signature (Kumar et al. 2017): residency/retention markers (CD69, CD103/ITGAE, CD49a/ITGA1), the CXCR6 chemokine receptor, and PDCD1, CRTAM and DUSP6.

**How to read it.** Resident memory T cells in the airway mucosa provide rapid local protection, and their maintenance changes with age.

**Caveats.** CD69 is also an early activation marker, so residency should not be inferred from CD69 alone. ITGAE (CD103) marks CD8+ TRM much more than CD4+ TRM, so the score leans towards CD8 cells; annotate CD4 and CD8 T cells as subtypes so the model adjusts for their mix. Scored only in T cells.

### Evidence

| Role | What it supports | Reference | Checked |
|---|---|---|---|
| Definition | Full text: the 31-gene core TRM signature (CD69+ vs CD69- memory T cells, lung and spleen, CD4 and CD8) includes upregulated ITGAE, ITGA1, CXCR6, PDCD1, DUSP6, CRTAM, IL2 and IL10 and downregulated S1PR1, KLF2, KLF3, SELL, CX3CR1, RAP1GAP1 and RGS1; CD69 defines the subset. | Kumar et al. 2017, *Cell Rep* ([PMID 28930685](https://pubmed.ncbi.nlm.nih.gov/28930685/)) | full text |
| Function | Nasal-epithelium CD8+ TRM persist long term and block influenza spread from the upper to the lower respiratory tract (mouse). | Pizzolla et al. 2017, *Sci Immunol* ([PMID 28783656](https://pubmed.ncbi.nlm.nih.gov/28783656/)) | abstract |
| Airway / nasal | Human nasal mucosa holds CD38+CD8+ TRM and SARS-CoV-2-specific CD8+ T cells that persisted at least 2 months after viral clearance. | Roukens et al. 2022, *Nat Immunol* ([PMID 34937933](https://pubmed.ncbi.nlm.nih.gov/34937933/)) | abstract |
| Caveat | Tissue dissociation induces gene-expression artefacts in scRNA-seq (from title only; no PubMed abstract, so whether CD69 or DUSP6 are affected was not verified here). | van den Brink et al. 2017, *Nat Methods* ([PMID 28960196](https://pubmed.ncbi.nlm.nih.gov/28960196/)) | metadata only |
| Caveat | Full text: CD101 is not in the 31-gene core; it was a lung-specific transcript, raised on CD8+CD69+ cells by flow but minimally on CD4+ TRM, and CD103 (ITGAE) is mainly CD8 TRM, so these genes bias the score to CD8 TRM. | Kumar et al. 2017, *Cell Rep* ([PMID 28930685](https://pubmed.ncbi.nlm.nih.gov/28930685/)) | full text |

### Evidence gaps

- No verified age reference.
- That tissue dissociation induces CD69 could not be verified (no PubMed abstract).

<a id="general"></a>
## General: cell types, QC and methods

| Role | What it supports | Reference | Checked |
|---|---|---|---|
| Airway / nasal | scRNA-seq of 77,969 cells from nose to 12th airway division in 10 healthy volunteers, including nasal brushings: 89% epithelial, 6% immune and 5% stromal cells. | Deprez et al. 2020, *Am J Respir Crit Care Med* ([PMID 32726565](https://pubmed.ncbi.nlm.nih.gov/32726565/)) | abstract |
| Caveat | scRNA-seq of human neutrophils is technically difficult, and changing the data-analysis pipeline (not the chemistry) greatly increases neutrophil detection; the abstract does not give gene counts. | Wigerblad et al. 2022, *J Immunol* ([PMID 35858733](https://pubmed.ncbi.nlm.nih.gov/35858733/)) | abstract |
| Caveat | nasal granulocytes/neutrophils are present in human nasal mucosal samples and rise during acute COVID-19 (flow cytometry, not scRNA-seq). | Roukens et al. 2022, *Nat Immunol* ([PMID 34937933](https://pubmed.ncbi.nlm.nih.gov/34937933/)) | abstract |
| Method | UCell scores gene signatures in single cells with the Mann-Whitney U statistic and is robust to dataset size and heterogeneity. | Andreatta et al. 2021, *Comput Struct Biotechnol J* ([PMID 34285779](https://pubmed.ncbi.nlm.nih.gov/34285779/)) | abstract |
| Method | overview of the Bioconductor single-cell workflow and the OSCA online book (the usual citation for adaptive MAD-based outlier QC); the abstract itself does not mention MAD thresholds. | Amezquita et al. 2020, *Nat Methods* ([PMID 31792435](https://pubmed.ncbi.nlm.nih.gov/31792435/)) | abstract |
| Method | UCell (R) and pyUCell (Python) implement rank-based single-cell signature scoring compatible with scanpy/scverse. | Andreatta et al. 2026, *Bioinformatics* ([PMID 41665985](https://pubmed.ncbi.nlm.nih.gov/41665985/)) | abstract |
| Method | across 5.5 million cells mtDNA% is higher in human than mouse tissues and a uniform 5% threshold fails in 13 of 44 human tissues, so tissue-specific mitochondrial thresholds are needed. | Osorio et al. 2021, *Bioinformatics* ([PMID 32840568](https://pubmed.ncbi.nlm.nih.gov/32840568/)) | abstract |

### Evidence gaps

- That neutrophils have low gene counts in droplet scRNA-seq, and are frequent in healthy nasal brushings, could not be verified from abstracts; the references show only that neutrophils are technically hard to capture and that nasal granulocytes rise in COVID-19.
- The MAD-based QC reference supports the general method; its abstract does not mention MAD thresholds explicitly.
