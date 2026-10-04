# ADC Peptide Mapper v1.0

**12-tab R Shiny application for in-silico ADC peptide mapping and LC-MS/MS method development.**

Covers the full workflow from FASTA upload through proteolytic digest, uniqueness checking, DAR modeling, transition list generation, MRM quality assessment, MS/MS database search, heavy labelling, and ADC design scoring -- all in one tool.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20681412.svg)](https://doi.org/10.5281/zenodo.20681412)

---

## Installation

```r
# install.packages("remotes")
remotes::install_github("nishi76/ADC_Peptide_Mapper")
ADCPeptideMapper::run_app()
```

Dependencies are declared in `DESCRIPTION` and installed automatically. For MS/MS search (Tab 6) see [Search Engine Setup](#ms/ms-search-engine-setup).

---

## Tabs

The application is split into two groups accessible from the sidebar.

### Sequence Analysis

| Tab | Name | What it does |
|-----|------|-------------|
| 1 | **Input & Setup** | Upload ADC FASTA (file or paste); name your ADC; select enzyme(s), missed cleavages, and background species; run digest |
| 2 | **Modifications** | Fixed CAM, variable mods (oxidation, propionamide, NEM), ADCDB drug-linker payloads, linker biotransformation mods, special PTMs, custom mod builder, DAR settings |
| 3 | **Peptide Results** | Browse and filter the full modified peptide table; sequence coverage map; co-uniqueness filter; CSV/Excel export |
| 4 | **Transition List** | Select instrument platform and DAR level; generate and download MRM/DIA transition lists for 6 instrument families |
| 5 | **Heavy Labelling** | Generate SIL-IS light/heavy peptide pairs; 6 isotope label presets plus custom |
| 6 | **MS/MS Search** | Run MS Amanda 3.0 or Tide/Crux locally; upload pre-computed results (mzIdentML, pepXML, psm.tsv); cross-reference PSMs with theoretical peptides; target-decoy FDR estimation |
| 7 | **MRM Assessment** | Score transition quality; peak simulation plots; ranked export by signal confidence |

### ADC Design

| Tab | Name | What it does |
|-----|------|-------------|
| 8 | **ADC Design** | Linker chemistry selector; payload database; deconjugation prediction with stability timeline |
| 9 | **Target Biology** | Internalization scoring; surface accessibility assessment for 50+ known ADC targets |
| 10 | **Antibody** | Binding affinity scoring; epitope characterization; half-life prediction; Fc engineering variants; bystander effect scoring |
| 11 | **PK & Efficacy** | ADC PK simulation; bystander kill modeling; DAR-PK interaction curves |
| 12 | **AI Assistant** | Anthropic Claude chat pre-loaded with ADC/proteomics context; digest context injection; model selector |

---

## Features

### Digest engine
- 11 enzymes: Trypsin, Trypsin/P, Lys-C, Lys-C/P, Lys-N, Asp-N, Glu-C (E/D), Arg-C, Chymotrypsin, Papain, Elastase
- Optional second enzyme for sequential dual-enzyme digestion
- 0, 1, or 2 missed cleavages
- Peptide length filter (configurable)

### Modifications
- **Fixed:** Carbamidomethylation (C, +57.021 Da)
- **Variable:** Oxidation (M), Propionamide (C), NEM (C), ADCDB drug-linker payloads (MMAE, DM1, DXd, SN-38, Duocarmycin, PBD, and more)
- **Linker biotransformations:** maleimide ring hydrolysis (+18.011 Da), succinimide ring-opening (+18.011 Da), thioether/sulfoxide oxidation (+15.995 Da), disulfide loss (-31.990 Da), conjugation-site deamidation (+0.984 Da)
- **Special PTMs:** Deamidation (N/Q), Pyroglutamate (N-term Q/E), Acetylation (K), Phosphorylation (S/T/Y)
- **Custom mod builder:** any residue, any mass shift, N-term / C-term / any-position

### Uniqueness checking
- Pre-built background proteomes: Human, Cynomolgus Monkey, Rat (UniProt Swiss-Prot reviewed + TrEMBL)
- Cross-species co-uniqueness: filter for peptides unique across all selected species simultaneously
- cRAP contaminant database integrated

### Transition list generation
- Full b/y ion series (b2..b(n-1), y2..y(n-1)); singly and doubly charged products
- Instrument-specific collision energy formulas for 6 platforms
- DAR0-DARn level-specific transition lists
- Averagine isotope envelope (Senko 1995) with recommended precursor isotope selection

### Instrument export formats (Tab 4)

| Platform | CE formula | Notes |
|----------|-----------|-------|
| Skyline | Generic | Direct Skyline import |
| Thermo Xcalibur / TSQ Altis | Linear, charge-dependent | HCD/CID optimised |
| SCIEX Analyst / QTRAP / TripleTOF | Empirical, charge-dependent | MRM and SWATH |
| Bruker timsControl / timsTOF | TIMS-adjusted | PASEF compatible |
| Agilent MassHunter / QQQ | Agilent empirical | MRM optimised |
| Waters MassLynx / Xevo TQ | Waters empirical | MRM optimised |

### Heavy label presets (Tab 5)

| Label | Residue | Mass shift |
|-------|---------|-----------|
| 13C6 15N2 Lys | K | +8.014199 Da |
| 13C6 15N4 Arg | R | +10.008269 Da |
| D4 Lys | K | +4.025107 Da |
| D6 Leu | L | +6.031817 Da |
| 13C6 Leu | L | +6.020129 Da |
| 13C9 15N1 Tyr | Y | +10.009369 Da |
| Custom | User-defined | User-defined |

### DAR modeling
- Cysteine thiol-maleimide, lysine NHS-ester/hydrazone, and site-specific conjugation
- Full MRM transition list per DAR level (DAR0 through DARn)
- Linker biotransformation variable mods applied per DAR species

### MRM Assessment (Tab 7)
- Transition quality scoring
- Peak shape simulation plots
- Ranked export by signal confidence
- DAR-level summary view

### ADC Design tabs (8-11)
- Linker and payload database with approved ADC reference data
- Deconjugation prediction with plasma stability timeline
- Internalization and surface accessibility scoring for 50+ targets
- Antibody half-life prediction across Fc engineering variants (LS, GASDALIE, YTE, etc.)
- Bystander effect scoring
- ADC PK simulation with DAR-dependent clearance modeling

### AI Assistant (Tab 12)
- Anthropic Claude API with ADC/proteomics system context
- Digest context injection (ADC name, chains, enzyme, peptide count)
- Dynamic model selector
- API key saved to `.Renviron` for persistence

---

## Setup

### 1. Install from GitHub

```r
install.packages("remotes")
remotes::install_github("nishi76/ADC_Peptide_Mapper")
ADCPeptideMapper::run_app()
```

Background databases (Human, Cyno, Rat) are bundled in the package -- no manual build step needed.

### 2. (Optional) MS/MS Search Engine

Tab 6 requires MS Amanda 3.0 or Tide/Crux installed on your machine. See [Search Engine Setup](#ms/ms-search-engine-setup) below.

### 3. (Optional) Anthropic API Key for AI Assistant

```r
Sys.setenv(ANTHROPIC_API_KEY = "sk-ant-...")
# Or paste it directly into Tab 12 and click "Save to .Renviron"
```

Get a key at: https://console.anthropic.com

---

## MS/MS Search Engine Setup

### MS Amanda 3.0 (recommended)

Free, standalone, no Java required. Developed at the Institute of Molecular Pathology (IMP), Vienna.

**Download:** https://github.com/hgb-bin-proteomics/MSAmanda/releases

| Platform | Binary | Notes |
|----------|--------|-------|
| Windows 10/11 | `MSAmanda.exe` | No installer needed |
| Linux | `MSAmanda` | Requires .NET 6 runtime |
| macOS | `MSAmanda` | Requires .NET 6 runtime |

**.NET 6 runtime (Linux/macOS):**
```bash
# Ubuntu/Debian
sudo apt-get install -y dotnet-runtime-6.0
# macOS
brew install --cask dotnet-runtime
```

### Tide / Crux 4.x (fallback)

**Download:** https://crux.ms/download.html

### Converting raw files (ProteoWizard MSConvert)

```bash
msconvert input.raw --mzML --filter "peakPicking true 1-"
```

**Download:** https://proteowizard.sourceforge.io

Detection order: explicit path in UI > `MSAMANDA_EXE` / `CRUX_EXE` environment variable > auto-scan of common directories > system PATH. MS Amanda is always preferred over Tide.

---

## Package structure

```
ADCPeptideMapper/
+-- R/
|   +-- digest.R                <- 11-enzyme digest engine
|   +-- modifications.R         <- PTM definitions, ADCDB payloads, linker biotransformations
|   +-- transitions.R           <- b/y ion series, CE calculation, DAR transitions
|   +-- isotopes.R              <- Averagine isotope envelope (Senko 1995)
|   +-- export.R                <- 6-platform instrument formatters
|   +-- uniqueness.R            <- Background proteome loading and uniqueness checking
|   +-- dar.R                   <- DAR distribution modeling
|   +-- msearch.R               <- MS Amanda + Tide engine detection, search, FDR, parsing
|   +-- adc_design.R            <- Linker/payload database, deconjugation prediction
|   +-- antibody_characterization.R  <- Binding affinity, epitope, half-life, bystander
|   +-- target_biology.R        <- Internalization scoring, surface accessibility
|   +-- run_app.R               <- run_app() entry point
|   +-- zzz.R                   <- Package-level imports and globalVariables
+-- inst/
|   +-- shiny/                  <- Bundled Shiny app (12 tabs)
|   +-- extdata/                <- Background databases (human, monkey, rat, cRAP)
+-- man/                        <- Auto-generated help pages (58 topics)
+-- DESCRIPTION
+-- NAMESPACE
```

---

## Running tests

```r
setwd("path/to/ADC_Peptide_Mapper_v1.0")
source("tests/test_masses.R")            # 79 unit tests
source("tests/benchmark_mass_accuracy.R") # 10 reference peptides, <= 0.05 mDa
```

---

## Citation

If you use ADC Peptide Mapper in your research, please cite:

```
Wase, N. (2026). ADC Peptide Mapper (Version 1.0) [Software].
https://doi.org/10.5281/zenodo.20681412
```

If using the MS/MS Search tab with MS Amanda, also cite:

```
Dorfer V, et al. MS Amanda, a Universal Identification Algorithm Optimized
for High Accuracy Tandem Mass Spectra. J Proteome Res. 2014;13(8):3679-3684.
doi:10.1021/pr500202e
```

If using Tide/Crux:

```
McIlwain S, et al. Crux: Rapid Open Source Protein Tandem Mass Spectrometry
Analysis. J Proteome Res. 2014;13(10):4488-4491. doi:10.1021/pr500741y
```

Isotope envelope: Senko MW, et al. Determination of monoisotopic masses and ion
populations for large biomolecules from resolved isotopic distributions.
J Am Soc Mass Spectrom. 1995;6(4):229-233.

See `CITATION.cff` for full metadata.

---

## License

MIT -- see `LICENSE` for details.

---

## Author

**Nishikant Wase, PhD** -- [nishikant.wase@gmail.com](mailto:nishikant.wase@gmail.com)  
Portfolio: [nishi76.github.io](https://nishi76.github.io)  
DOI: [10.5281/zenodo.20681412](https://doi.org/10.5281/zenodo.20681412)

*For research use only. Monoisotopic masses throughout. Background databases sourced from UniProt Swiss-Prot reviewed proteomes.*
