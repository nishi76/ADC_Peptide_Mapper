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
