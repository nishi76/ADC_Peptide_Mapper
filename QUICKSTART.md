# ADC Peptide Mapper v1.0 — Quick Start Guide

**Goal:** get from zero to a ready-to-import MRM transition list in under 10 minutes.

---

## Step 0 — Install

```r
install.packages("remotes")
remotes::install_github("nishi76/ADC_Peptide_Mapper")
ADCPeptideMapper::run_app()
```

The app opens in your default browser. Background proteome databases (Human, Cyno, Rat) are bundled -- no build step required.

---

## Step 1 — Tab 1: Input & Setup

1. **ADC Name** -- type a short identifier (e.g. `TDM1_batch1`). This label appears in all exports.
2. **FASTA input** -- either:
   - Upload a `.fasta` file with your heavy and light chain sequences, or
   - Paste sequences directly into the text area, or
   - Click **Load Demo (Trastuzumab)** to try the pre-loaded example.
3. **Enzyme** -- select your primary enzyme. Trypsin is the default for most LC-MS/MS workflows. Lys-C or a Trypsin + Lys-C combination is useful for incomplete digestion studies.
4. **Second enzyme (optional)** -- enable for sequential two-enzyme digestion.
5. **Missed cleavages** -- 0 is standard; increase to 1 or 2 to capture missed cleavage peptides.
6. **Peptide length** -- default 6-40 AA. Adjust if you expect very short or long peptides.
7. **Background species** -- tick Human, Cynomolgus Monkey, and/or Rat for uniqueness filtering.
8. Click **Run Digest & Uniqueness Check**.

The digest runs in seconds. A status banner confirms the number of peptides generated and unique peptides found.

---

## Step 2 — Tab 2: Modifications

Configure the modification state before inspecting results. Changes here trigger an automatic re-scoring of the peptide table.

**Fixed modification**
- Carbamidomethylation (C, +57.021 Da) is on by default. Disable only if your sample was not alkylated.

**Variable modifications**
- Tick **Oxidation (M)** for routine samples.
- Tick any **ADCDB drug-linker payload** matching your conjugate (MMAE, DM1, DXd, SN-38, etc.). Each payload adds a selectable mass shift to conjugation-site residues.
- Tick **Linker biotransformations** relevant to your chemistry (e.g. maleimide ring hydrolysis for thiol-maleimide ADCs).

**Special PTMs**
- Add Deamidation (N/Q), Pyroglutamate, Acetylation, or Phosphorylation as needed.

**Custom modifications**
- Click **Add Custom Mod** to define any residue + mass shift combination.

**DAR settings**
- Set conjugation chemistry (Cys thiol-maleimide, Lys NHS-ester, or site-specific) and the maximum DAR level. This drives Tab 4 DAR-level transition generation.

---

## Step 3 — Tab 3: Peptide Results

The table shows every theoretical peptide with:

| Column | Meaning |
|--------|---------|
| Sequence | Amino acid sequence |
| Modified sequence | Sequence with modification annotations |
| UniqueToADC | TRUE if not present in any selected background proteome |
| PeptideLength | Number of residues |
| MissedCleavages | MC count |
| Charge states | Predicted charge states |
| Precursor m/z | Monoisotopic m/z per charge state |

**Filtering tips**
- Tick **Unique peptides only** to focus on ADC-specific peptides.
- Use the **Co-uniqueness** toggle to require uniqueness across all selected species simultaneously (stricter than per-species).
- The **Sequence Coverage Map** card below the table visualises peptide positions on each chain. Colour by Uniqueness, Missed Cleavages, or Length. Download as PNG at 300 dpi.

**Export**
- **Download CSV** exports the filtered table.
- **Download Excel** produces a multi-sheet workbook (Transition List, Peptide Summary, Unique Peptides, Instrument Reference).

---

## Step 4 — Tab 4: Transition List

1. **Instrument** -- select the platform you will run on (Skyline, Thermo, SCIEX, Bruker, Agilent, Waters).
2. **DAR level** -- "All" generates a combined list; individual DAR levels (DAR0, DAR2, etc.) generate level-specific lists for DAR-resolved quantitation.
3. **Unique only** -- tick to restrict the export to ADC-specific peptides.
4. **Top N ions** -- default 5 fragment ions per precursor, ranked by product m/z.
5. Click **Generate Transition List**.
6. Click **Download CSV** to save the instrument-ready file.

**Column layouts per instrument**

| Platform | Key columns |
|----------|------------|
| Skyline | Protein, Peptide, Precursor m/z, Product m/z, CE, Charge, ... |
| Thermo Xcalibur | Compound, Precursor (m/z), Product (m/z), CE, ... |
| SCIEX Analyst | Q1, Q3, CE, Declustering Potential, ... |
| Bruker timsControl | CompoundName, Precursor m/z, Product m/z, CE, ... |
| Agilent MassHunter | Compound Name, Precursor Ion, Product Ion, Collision Energy, ... |
| Waters MassLynx | Parent, Daughter, Cone, Collision Energy, Component Name, ... |

For Skyline users, import the CSV directly via **File > Import > Transition List**.

---

## Step 5 — Tab 5: Heavy Labelling (optional)

Use this tab to generate stable-isotope labelled (SIL) internal standard pairs for quantitative LC-MS/MS.

1. Select a **label preset** (e.g. 13C6 15N2 Lys for K-labelled SILAC peptides).
2. The table shows light and heavy peptide pairs with their mass shifts.
3. Download as CSV or add to your Skyline document alongside the unlabelled transition list.

---

## Step 6 — Tab 6: MS/MS Search (optional)

Use this tab to confirm theoretical peptides against real LC-MS/MS data.

**If you have a search engine installed:**
1. The status badge at the top of the tab shows which engine was detected (MS Amanda 3.0 preferred).
2. Upload your mzML/mzXML/MGF spectral files.
3. Click **Run Search**.
4. Adjust the score threshold slider; results are filtered in real time.
5. The **Sequence Coverage** sub-tab shows which theoretical peptides were confirmed.

**If you have pre-computed results:**
- Upload directly: `.mzid` (MS Amanda), `.pepxml` (Tide/Crux), or `psm.tsv` (FragPipe/MSFragger). Format is auto-detected.

**FDR estimation**
- Target-decoy FDR is calculated automatically. A 1% FDR threshold is highlighted in the table.

For engine installation instructions see [the README](README.md#ms/ms-search-engine-setup).

---

## Step 7 — Tab 7: MRM Assessment (optional)

After generating a transition list, use this tab to rank transitions by predicted signal quality.

- **All Transitions** -- full list with quality scores.
- **DAR-Level Summary** -- per-DAR signal comparison.
- **Peak Plots** -- simulated peak shapes for visual QC.
- **Ranked Export** -- download the top-ranked transitions ready for final method setup.

---

## ADC Design tabs (8-11)

These tabs are independent of the digest workflow and can be used at any time.

**Tab 8 -- ADC Design**
- Browse the linker and payload database with approved ADC reference data.
- Select a linker chemistry and payload to see predicted deconjugation kinetics plotted over time in plasma vs. tumour microenvironment conditions.

**Tab 9 -- Target Biology**
- Enter a gene/protein target name.
- **Internalization** scoring estimates lysosomal delivery efficiency based on receptor biology.
- **Surface accessibility** scores how accessible the epitope is to a large IgG molecule, incorporating shedding and secretion flags.

**Tab 10 -- Antibody**
- **Binding affinity** scoring interprets your KD against optimal ADC ranges for solid tumour vs. haematologic targets.
- **Epitope characterization** flags binding-site barrier risk and cross-reactivity considerations.
- **Half-life prediction** compares wild-type IgG1 Fc against engineering variants (LS, YTE, GASDALIE, XTEND, etc.).
- **Bystander effect** scoring estimates payload permeability and diffusion potential.

**Tab 11 -- PK & Efficacy**
- **PK simulation** models total antibody, conjugated ADC, and free drug plasma concentrations over time.
- **Bystander kill** estimates therapeutic window in antigen-heterogeneous tumour models.

---

## Tab 12 -- AI Assistant

1. Go to the **AI Assistant** tab.
2. Paste your Anthropic API key and click **Save to .Renviron** (persists across sessions).
3. Select a Claude model from the dropdown.
4. Enable **Include digest context** to automatically attach your current ADC name, chains, enzyme, and peptide count to each message.
5. Type your question and press **Send** or Enter.

Example questions:
- "Which of my unique peptides are best suited for MRM quantitation?"
- "Explain why maleimide hydrolysis matters for ADC LC-MS/MS analysis."
- "What DAR level should I target for optimal therapeutic index?"

---

---

## Sequence Analysis tools in depth

The seven Sequence Analysis tabs form a linear pipeline. Each step feeds into the next, but you can jump back at any point to change settings and re-run.

### Tab 1 -- Input & Setup

This is where the analysis starts. The app accepts any multi-chain FASTA: heavy chain (HC) and light chain (LC) are auto-detected by length and sequence composition. You can also paste raw sequence blocks without FASTA headers -- the app assigns chain labels automatically. The demo Trastuzumab sequence is useful for learning the workflow before using your own molecule.

Enzyme selection drives everything downstream. Trypsin is the gold standard for most LC-MS/MS ADC bioanalysis because it produces peptides in the 800-3000 Da range with predictable chromatography. Lys-C produces longer, more hydrophobic peptides useful for complementary coverage. Glu-C (E/D) is common for middle-down workflows. The optional second enzyme enables sequential digestion -- useful for studying incomplete digestion or confirming borderline unique peptides by orthogonal cleavage.

Missed cleavages = 0 keeps the peptide list clean and realistic for a fully digested sample. Increasing to 1 captures the 15-30% of peptides that are typically only partially cleaved in practice, which matters when building MRM methods that need to monitor real biological samples.

The background species selection determines which proteomes are used for uniqueness filtering in Tab 3. If your ADC will be dosed in cynomolgus monkey toxicology studies, tick Monkey -- peptides flagged as unique against human but not monkey can confound cross-species PK/bioanalytical bridging.

### Tab 2 -- Modifications

Modifications determine the exact mass of every peptide and whether a peptide can be detected as the native vs. modified form. Getting this right before reviewing results is important.

**Fixed mods** are applied to every occurrence of the target residue, no exception. Carbamidomethylation on Cys (+57.021 Da) is the standard because iodoacetamide alkylation is nearly quantitative. If your sample prep used a different alkylating agent (NEM, propionamide), switch here.

**Variable mods** are applied in all combinations to generate every possible modified form. Each additional variable mod multiplies the peptide table. Keep the list to only what is relevant for your sample -- over-modifying causes noise. For an ADC, the most important variable mod is the drug-linker payload mass on the conjugation site residue. The ADCDB payload library covers all major approved and clinical-stage payloads with their exact linker-conjugated mass shifts.

**Linker biotransformations** are variable mods representing chemical degradation of the linker after conjugation. These appear in plasma samples and stressed stability samples. Maleimide ring hydrolysis and succinimide ring-opening are the most common for thiol-maleimide linked ADCs and add +18.011 Da to the conjugation-site Cys peptide. Including these in your transition list is important for complete DAR-level monitoring in PK samples.

**DAR settings** tell the app how many conjugation sites exist and what the conjugation chemistry is. The app then generates a separate modification state and precursor mass for each DAR level from DAR0 (naked antibody) through DARn (fully loaded). This is essential for DAR-resolved MRM -- each DAR species is a distinct mass and must be monitored with its own transitions.

### Tab 3 -- Peptide Results

The peptide table is the master output of the digest engine. Every row is one theoretical peptide after all selected modifications are applied. The `UniqueToADC` flag is the most important filter for bioanalytical method development -- only unique peptides can be used as signature peptides for ADC quantitation in complex biological matrices.

The **Sequence Coverage Map** deserves particular attention. It lays out all peptides as colour-coded bars across the full HC and LC sequences. Gaps in coverage indicate regions where no peptide was generated -- often because of very short or very long sequences, or because they fall entirely within missed cleavage regions. Coloring by uniqueness immediately shows which parts of the antibody are covered by unique peptides vs. which regions are shared with background proteins. This is especially useful for Fc-region peptides, which are often not unique against the human background because the Fc sequence is conserved across IgG subclasses.

The co-uniqueness toggle applies a stricter filter: a peptide must be unique against all selected species simultaneously. This is the conservative choice for cross-species bioanalytical methods.

### Tab 4 -- Transition List

The transition list tab converts theoretical peptides into ready-to-import instrument methods. Select your instrument platform first -- the app applies the correct collision energy formula and output column layout for each vendor. CE formulas are empirically derived from published ADC and peptide MRM literature for each instrument family.

The DAR level selector generates DAR-resolved methods. For a DAR2/DAR4 mixture (typical for cysteine-conjugated ADCs), you generate three separate transition lists (DAR0, DAR2, DAR4) and use the precursor masses to set up separate MRM channels per DAR species. This enables DAR-resolved pharmacokinetics -- the most complete picture of ADC catabolism in vivo.

For Skyline users: import the CSV via File > Import > Transition List. Skyline recognises the column headers and populates all precursor and product m/z fields automatically.

Top N filtering (default 5 ions per precursor) is a practical limit for MRM cycle time. In a typical 15-30 min gradient with 20-50 precursors to monitor, 5 transitions each keeps the dwell time reasonable at 20-50 ms per channel.

### Tab 5 -- Heavy Labelling

Stable-isotope labelled (SIL) peptides are the internal standard of choice for absolute quantitation of ADC signature peptides. The SIL peptide and its light counterpart co-elute exactly, share the same ionisation efficiency, and differ only in mass -- making the light/heavy ratio a robust quantitative readout independent of matrix effects and instrument response drift.

Select the isotope label that matches your SIL peptide vendor specification (most commercial SIL peptides for IgG use 13C6 15N2 Lys or 13C6 15N4 Arg). The table shows both the light peptide (from Tab 3) and its heavy counterpart with the mass shift applied. Download this as a separate CSV and add it to your Skyline document alongside the unlabelled transition list.

### Tab 6 -- MS/MS Search

This tab closes the loop between theory and experiment. After running LC-MS/MS, use this tab to confirm which theoretical peptides were actually detected.

Two workflows are supported. If you have a supported search engine installed, upload your raw spectral files (mzML is recommended; the app also accepts mzXML and MGF) and click Run Search. The engine runs locally on your machine -- no data leaves your environment. The result is a PSM table cross-referenced against your theoretical peptide list from Tab 3.

If you already have search results from another tool (MS Amanda, MSFragger, FragPipe, Tide), upload the output file directly. The parser auto-detects mzIdentML, pepXML, and tab-delimited PSM formats.

The score threshold slider and FDR display update in real time. Target-decoy FDR estimation follows the standard Kall 2008 approach -- use the 1% FDR threshold as your starting point and tighten it for high-confidence signature peptide confirmation.

The Sequence Coverage sub-tab overlays confirmed PSMs onto the coverage map from Tab 3, so you can immediately see which theoretical peptides were experimentally confirmed and which remained undetected.

### Tab 7 -- MRM Assessment

After generating a transition list, this tab helps you prioritise which transitions to carry into the final method. Transitions are scored on signal quality metrics including precursor charge state, product ion type (y-ions generally outperform b-ions for ESI), product m/z (avoid low-mass region < 200 Da), and co-elution with known interferences.

The Peak Plots sub-tab simulates the expected peak shape for each transition based on its isotope envelope and charge state -- useful for visual QC before committing to method setup on the instrument. The Ranked Export downloads a filtered list of your highest-confidence transitions sorted by score, ready for instrument programming.

---

## ADC Design tools in depth

Tabs 8-11 are a standalone decision support suite for ADC candidate evaluation. They are independent of the digest workflow and can be used at any point, even before a molecule exists in the lab.

### Tab 8 -- ADC Design

The linker and payload database covers all major approved ADCs (Kadcyla, Enhertu, Padcev, Trodelvy, Zynlonta, Besylomab, etc.) and clinical-stage molecules, with linker chemistry, payload class, DAR, and PK annotations drawn from published literature.

The **Deconjugation Prediction** sub-tab models the chemical stability of your chosen linker in plasma vs. the tumour microenvironment over time. For thiol-maleimide linked ADCs, maleimide hydrolysis and thiosuccinimide ring-opening are the dominant degradation pathways. The model projects the fraction of intact ADC vs. deconjugated species over a 14-day period and flags stability windows that are consistent with the dosing interval.

This is a qualitative guide, not a quantitative pharmacokinetic model -- it assumes literature-derived rate constants and is most useful for comparing relative linker stability rather than predicting absolute plasma half-life.

### Tab 9 -- Target Biology

Target selection is one of the most consequential decisions in ADC development. A target with low internalization, poor surface accessibility, or rapid shedding will undermine even an optimally designed ADC.

The **Internalization** scorer covers 50+ validated ADC targets (HER2, HER3, TROP2, Nectin-4, FRalpha, CD22, CD30, CD33, B7-H3, EGFR, c-Met, and more). Scores reflect the depth and speed of receptor internalization based on receptor biology -- RTKs and CD antigens with constitutive endocytosis score highest. The score is on a 0-100 scale with an accompanying tier label (High / Moderate / Low) and a summary of the receptor biology rationale.

The **Surface Accessibility** scorer assesses how readily a large IgG-drug conjugate can access the epitope. Factors include receptor density, extracellular domain count, transmembrane topology, and whether the target is shed or secreted into circulation. Shed targets create a sink effect that reduces effective ADC delivery to tumour cells.

For novel targets not in the database, both scorers return an "Unknown" tier with suggested characterisation experiments.

### Tab 10 -- Antibody

**Binding affinity** scoring contextualises your measured KD in terms of ADC-specific considerations. Unlike naked antibody therapeutics, an ADC that binds too tightly can cause a binding-site barrier effect -- the ADC is retained in the tumour periphery around blood vessels and fails to penetrate to antigen-low cells deeper in the mass. The optimal KD for solid tumour ADCs is generally 1-10 nM. Haematologic targets tolerate tighter binding (0.1-1 nM) because diffusion barriers are less of a limiting factor.

**Epitope characterization** flags cross-reactivity risk. Conformational epitopes (discontinuous binding surface) are generally more target-specific than linear epitopes but are harder to reproduce with synthetic peptides for competitive binding assays. Linear epitopes warrant a BLAST search of the epitope sequence against the human proteome to assess off-target binding risk.

**Half-life prediction** covers nine Fc engineering variants including LS (M428L/N434S, ~4x extension), YTE (M252Y/S254T/T256E, ~3.5x), GASDALIE (G236A/S239D/A330L/I332E, enhanced FcgammaRIII), and XTEND. The model uses published FcRn binding affinities and extrapolates from wild-type IgG1 half-life adjusted for the variant's FcRn affinity ratio. Half-life extension matters for ADC dosing frequency and the exposure window needed for target engagement.

**Bystander effect scoring** estimates the likelihood that free payload released after lysosomal degradation will diffuse into neighbouring antigen-negative tumour cells. High-permeability payloads (MMAE, Duocarmycin, PBD) drive strong bystander kill, which is advantageous for antigen-heterogeneous tumours but raises safety concerns in normal tissue adjacent to tumour.

### Tab 11 -- PK & Efficacy

The PK simulation uses a two-compartment model with DAR-dependent clearance. As DAR increases, hydrophobicity increases and clearance accelerates -- DAR2 and DAR4 species have markedly different exposure profiles from DAR0 (naked antibody). The simulation plots total antibody, conjugated ADC, and free drug concentrations over time for the selected dosing regimen.

The bystander kill model estimates the effective concentration of free payload in antigen-negative cells as a function of tumour antigen density, payload permeability, and ADC dose. This is particularly useful for evaluating ADCs against heterogeneous solid tumours where a subset of cells do not express the target antigen.

These simulations use literature-derived PK parameters scaled to your payload/linker selection and are intended for qualitative comparison and hypothesis generation, not for regulatory submission or clinical dosing decisions.

---

## Typical workflows

### Workflow A -- MRM method development from scratch

Tab 1 (digest) > Tab 2 (mods + DAR) > Tab 3 (review unique peptides) > Tab 4 (export for your instrument) > Tab 5 (SIL pairs for quantitation)

### Workflow B -- Confirm identity with existing LC-MS/MS data

Tab 1 > Tab 2 > Tab 6 (search or upload PSMs) > Tab 3 (coverage map with confirmed peptides highlighted) > Tab 4 (export confirmed transitions only)

### Workflow C -- ADC candidate evaluation

Tab 8 (select linker + payload) > Tab 9 (target internalization + accessibility) > Tab 10 (antibody affinity + half-life) > Tab 11 (PK/efficacy simulation)

---

## Troubleshooting

| Problem | Fix |
|---------|-----|
| "Failed to load R/digest.R" | Reinstall the package: `remotes::install_github("nishi76/ADC_Peptide_Mapper")` |
| Background databases not found | Check the package was installed from GitHub, not run from source |
| No engine detected in Tab 6 | Set `MSAMANDA_EXE` or `CRUX_EXE` env var, or enter the path in the Tab 6 UI |
| AI Assistant returns no models | Check your Anthropic API key is valid and has access |
| Very slow digest | Reduce missed cleavages to 0, or narrow the peptide length range |

---

## Getting help

- GitHub issues: https://github.com/nishi76/ADC_Peptide_Mapper/issues
- Email: nishikant.wase@gmail.com
- R help pages: `?ADCPeptideMapper::parse_fasta`, `?ADCPeptideMapper::generate_transition_list`, etc.

---

*ADC Peptide Mapper v1.0 -- Nishikant Wase, PhD -- For research use only.*
