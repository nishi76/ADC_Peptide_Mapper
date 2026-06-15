# ADC Peptide Mapper v0.8

**R Shiny application** — In-silico proteolytic digest, uniqueness checking, DAR distribution modeling, linker biotransformation variable modifications, instrument-specific transition list export, heavy labelling, MS/MS search confirmation, and an AI assistant for Antibody-Drug Conjugates.

**Live app:** https://nishiw.shinyapps.io/ADC_Peptide_Mapper/

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20681412.svg)](https://doi.org/10.5281/zenodo.20681412)

---

## What's New in v0.8

| Feature | v0.7 | v0.8 |
|---|---|---|
| Tabs | 6 | 7 (+ AI Assistant) |
| DAR modeling | — | DAR0–DARn distribution with per-level MRM transition lists |
| Conjugation chemistry | — | Cysteine thiol-maleimide, lysine NHS-ester/hydrazone, site-specific classification |
| Linker biotransformations | — | 5 variable mods: maleimide hydrolysis, succinimide ring-opening, thioether oxidation, disulfide loss, deamidation near conjugation site |
| Isotope envelope | — | Averagine-based isotope distribution (Senko 1995) |
| FDR estimation | — | Target-decoy FDR (Käll 2008) in MS/MS Search tab |
| Cross-species filtering | — | Co-species uniqueness check (Human / Cyno / Rat simultaneously) |
| cRAP contaminants | — | Common Repository of Adventitious Proteins integrated into background build |
| Search engine support | MSFragger (commercial) | MS Amanda 3.0 (primary) + Tide/Crux 4.x (fallback) — both free |
| Mass accuracy benchmark | — | Sub-mDa accuracy verified against IgG1 tryptic reference peptides |
| Unit tests | — | 79 tests, 100% pass rate (`tests/test_masses.R`) |
| AI Assistant | — | Tab 7: Anthropic Claude API chat with ADC/proteomics context |
| Sidebar citation | — | Author, DOI, and contact always visible in sidebar |
| Run Digest placement | Separate row below inputs | Embedded in FASTA Input card, eliminating empty space |
| Citation DOI | Placeholder | Real Zenodo DOI: `10.5281/zenodo.20681412` |

---


## Features

- **FASTA upload** — multi-chain ADC (HC + LC auto-detected); demo Trastuzumab sequence pre-loaded
- **11-enzyme digest engine** — Trypsin, Trypsin/P (incl. proline), Lys-C, Lys-C/P (incl. proline), Lys-N, Asp-N, Glu-C (E/D), Arg-C, Chymotrypsin (F/Y/W), Papain (K/R/Q), Elastase (A/V/S/G/T)
- **Optional second enzyme** — sequential dual-enzyme digestion
- **Missed cleavages** — 0, 1, or 2 (selectable per run)
- **Fixed mod** — Carbamidomethylation (CAM, +57.021 Da on C)
- **Variable mods** — Oxidation (M), Propionamide (C), NEM (C), ADCDB drug-linker payloads (MMAE, DM1, DXd, SN-38, and more)
- **Special mods** — Deamidation (N/Q), Pyroglutamate (Q/E, N-term), Acetylation (K), Phosphorylation (S/T/Y)
- **Linker biotransformation mods** — maleimide ring hydrolysis (+18.011 Da), succinimide ring-opening (+18.011 Da), thioether→sulfoxide (+15.995 Da), disulfide loss (−31.990 Da), conjugation-site deamidation (+0.984 Da)
- **DAR distribution modeling** — DAR0–DARn with full MRM transition list per DAR level
- **Custom mod builder** — any residue, any mass shift, N-term / C-term / any-position location
- **Uniqueness check** — vs pre-built Human, Cynomolgus Monkey, and Rat backgrounds (UniProt Swiss-Prot reviewed + TrEMBL)
- **Cross-species co-uniqueness** — filter peptides unique across all selected species simultaneously
- **cRAP contaminant integration** — common laboratory contaminants flagged in background
- **Sequence coverage map (Tab 3)** — visualisation; colour by uniqueness, missed cleavages, or length; PNG download
- **Full b/y ion series** — b2..b(n-1) and y2..y(n-1); singly charged products
- **Averagine isotope envelope** — monoisotopic + isotope distribution per peptide
- **Instrument-specific export (Tab 4)** — per-platform collision energy formulas and CSV column layouts for 6 instrument families
- **Heavy Labelling (Tab 5)** — 6 isotope label presets + custom; light/heavy peptide pairs with mass shifts
- **MS/MS Search (Tab 6)** — MS Amanda 3.0 (primary) or Tide/Crux (fallback); target-decoy FDR estimation; dynamic score filter UI
- **AI Assistant (Tab 7)** — Anthropic Claude API chat with ADC/proteomics system context; dynamic model selection; API key persisted to `.Renviron`
- **Mass accuracy benchmark** — sub-mDa accuracy validated against IgG1 Fc tryptic reference peptides (GPSVFPLAPSSR, ELASGLSFPVGFK, CASIQKFGR, DTLMISR)
- **Unit test suite** — 79 tests covering mass functions, enzyme cleavage, DAR, isotopes, and constants (`tests/test_masses.R`)

---

## Setup

### 1. Install R Packages (one-time)

```r
install.packages(c(
  "shiny", "bs4Dash", "DT", "data.table", "openxlsx",
  "httr2", "stringr", "dplyr", "shinycssloaders", "shinyjs", "htmltools"
))
```

Bioconductor packages (required for mzIdentML / pepXML parsing in Tab 6):

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("Biostrings", "mzR"))
install.packages("XML")   # for pepXML / mzIdentML parsing
```

### 2. Build Background Databases (one-time, ~10 min)

```r
setwd("path/to/ADC_Peptide_Mapper_v0.8")
source("build_background_db.R")
```

Creates three files in `data/`:
- `bg_human.rds` — UniProt Swiss-Prot human reviewed (~20,400 proteins)
- `bg_monkey.rds` — Cynomolgus monkey (~1,200 proteins)
- `bg_rat.rds` — Rat (~8,200 proteins)

> **Note:** If you already have these files from v0.7, copy them into the v0.8 `data/` folder — no rebuild needed.

> **Cloning from GitHub?** The `.rds` files are excluded from the repository (they are 50–150 MB each and exceed GitHub's file size limit). After cloning, run `source("build_background_db.R")` once to generate them locally.

### 3. (Optional) Set Anthropic API Key for AI Assistant

To use the AI Assistant tab, you need a free or paid Anthropic API key:

```r
# In the app — paste key into the AI Assistant tab and click "Save to .Renviron"
# The key is then loaded automatically on every future session.

# Or set it manually before launching:
Sys.setenv(ANTHROPIC_API_KEY = "sk-ant-...")
```

Get a key at: https://console.anthropic.com

### 4. Launch the App

```r
shiny::runApp("app.R")
```

*Tip: In RStudio, open `app.R` and click the **Run App** button.*

---

## Tabs

| Tab | Name | Purpose |
|---|---|---|
| 1 | **Input & Setup** | Upload ADC FASTA, name your ADC, select enzyme(s), missed cleavages, background species, run digest |
| 2 | **Modifications** | Fixed, variable, special, linker biotransformation, and custom PTMs; ADCDB drug-linker payloads; DAR settings |
| 3 | **Peptide Results** | Browse, filter, and export the full modified peptide table; standalone sequence coverage map; co-uniqueness check |
| 4 | **Transition List** | Select instrument platform, DAR level, generate and download MRM/DIA transition lists |
| 5 | **Heavy Labelling** | Generate SIL-IS light/heavy peptide pairs for quantitative LC-MS/MS |
| 6 | **MS/MS Search** | Run MS Amanda 3.0 or Tide/Crux locally; cross-reference PSMs with theoretical peptides; FDR estimation |
| 7 | **AI Assistant** | Claude-powered chat for ADC peptide mapping questions; uses your Anthropic API key |

---

## Tab 7 — AI Assistant

The AI Assistant tab embeds an Anthropic Claude chat interface pre-configured with ADC/proteomics context. It can help interpret results, explain mass spectrometry concepts, suggest experimental designs, and troubleshoot workflows.

**Setup:**
1. Go to the **AI Assistant** tab
2. Paste your Anthropic API key into the key field
3. Click **Save to .Renviron** — the key persists across sessions
4. Select a Claude model from the dropdown (populated from your account's available models)
5. Type your question and click **Send** (or press Enter)

**Context toggle:** enable "Include digest context" to automatically attach the current ADC name, chains detected, enzyme, and peptide count to each message.

---

## Tab 3 — Sequence Coverage Map

After running the digest, a **Sequence Coverage Map** card appears below the peptide results table. It shows all theoretical peptides mapped onto each chain using a greedy lane-assignment algorithm (no overlapping bars).

**Controls:**
- **Show chain** — display all chains or one at a time
- **Colour by** — Uniqueness (navy = unique, grey = non-unique), Missed cleavages (dark → light blue), or Peptide length (viridis scale)
- **Show peptide labels** — overlay sequence text on peptides ≥ 8 AA
- **Download PNG** — 300 dpi export

Coverage percentage is annotated per chain (e.g. "HC: 74.0% covered (333 / 450 AA)").

---

## Tab 6 — MS/MS Search Engine Setup

Tab 6 runs MS Amanda 3.0 (primary) or Tide/Crux (fallback) on your local machine and cross-references PSM results against the theoretical peptide list from Tab 3.

### Engine detection order

1. **MS Amanda 3.0** — checked first
   - Explicit path in the "Engine executable path" field
   - `MSAMANDA_EXE` environment variable
   - Auto-scan: `MSAmanda` / `MSAmanda.exe` in `getwd()`, `~/tools/`, `~/bin/`, `~/MSAmanda/`
   - System PATH

2. **Tide/Crux 4.x** — checked only if MS Amanda not found
   - `CRUX_EXE` environment variable
   - Auto-scan: `crux` / `crux.exe` in `getwd()`, `~/tools/`, `~/bin/`, `~/crux/bin/`
   - System PATH

The status badge in Tab 6 shows which engine is active and links to download pages if neither is found.

### Score filter UI

The score slider label and range change automatically based on the detected engine:
- **MS Amanda:** "Minimum Amanda Score" (0–2000, default 100)
- **Tide/Crux:** "Minimum XCorr" (0–10, default 1.5)

### Path A — upload pre-computed results

Upload any of the following directly (skips the search step):
- `.mzid` / `.mzidentml` — MS Amanda output
- `.pepxml` / `.pep.xml` — Tide/Crux output
- `psm.tsv` — FragPipe / MSFragger legacy
- Amanda summary `.csv`

Format is auto-detected from file extension and content.

---

## MS/MS Search Engine Installation

### MS Amanda 3.0 (Primary — recommended)

MS Amanda is a free, standalone peptide identification engine developed at the Institute of Molecular Pathology (IMP), Vienna. Designed for high-resolution Orbitrap data. No Java, no licence fee.

**Download:** https://github.com/hgb-bin-proteomics/MSAmanda/releases

| Platform | Binary name | Notes |
|---|---|---|
| Windows 10/11 (x64) | `MSAmanda.exe` | Standalone `.exe`; no installer needed |
| Linux (x86_64) | `MSAmanda` | Requires .NET 6 runtime |
| macOS (Intel / Apple Silicon) | `MSAmanda` | Requires .NET 6 runtime |

**.NET 6 runtime (Linux/macOS):**
```bash
# Ubuntu/Debian
sudo apt-get install -y dotnet-runtime-6.0
# macOS (Homebrew)
brew install --cask dotnet-runtime
```

### Tide / Crux 4.x (Fallback)

**Download:** https://crux.ms/download.html — statically linked, no dependencies.

### Converting Raw Files (ProteoWizard MSConvert)

```bash
msconvert input.raw --mzML --filter "peakPicking true 1-"
```

**Download:** https://proteowizard.sourceforge.io

### Requirements Summary

| Component | Required for | Source | Free? |
|---|---|---|---|
| R ≥ 4.2 | App runtime | https://cran.r-project.org | Yes |
| R packages (see §1) | App runtime | CRAN / Bioconductor | Yes |
| MS Amanda 3.0 | Tab 6 (primary) | https://github.com/hgb-bin-proteomics/MSAmanda/releases | Yes |
| .NET 6 runtime | MS Amanda on Linux/Mac | https://dotnet.microsoft.com/download/dotnet/6.0 | Yes |
| Crux 4.x (Tide) | Tab 6 (fallback) | https://crux.ms/download.html | Yes |
| ProteoWizard MSConvert | Raw file conversion | https://proteowizard.sourceforge.io | Yes |
| Anthropic API key | Tab 7 AI Assistant | https://console.anthropic.com | Free tier available |
| Internet access | `build_background_db.R` only | — | — |

---

## Supported Instruments (Tab 4)

| Platform | CE Formula | Notes |
|---|---|---|
| Thermo (Orbitrap/TSQ) | Linear, charge-dependent | HCD/CID optimised |
| SCIEX (QTRAP/TripleTOF) | Empirical, charge-dependent | MRM & SWATH |
| Bruker (timsTOF) | TIMS-adjusted | PASEF compatible |
| Agilent (QQQ/QTOF) | Agilent empirical | MRM optimised |
| Waters (Xevo/Synapt) | Waters empirical | MRM optimised |
| Skyline (generic) | Sciex-style default | Direct Skyline import |

---

## Heavy Label Presets (Tab 5)

| Label | Residue | Mass Shift (Da) |
|---|---|---|
| 13C6 15N2 Lys | K | +8.014199 |
| 13C6 15N4 Arg | R | +10.008269 |
| D4 Lys | K | +4.025107 |
| D6 Leu | L | +6.031817 |
| 13C6 Leu | L | +6.020129 |
| 13C9 15N1 Tyr | Y | +10.009369 |
| Custom | User-defined | User-defined |

---

## DAR Modeling (Tab 2 / Tab 4)

Drug-to-Antibody Ratio (DAR) distribution modeling generates a complete MRM transition list for each DAR species (DAR0 through DARn). Each DAR level adds the appropriate number of drug-linker payload mass units to the conjugated peptide(s), reflecting the statistical distribution of conjugation sites in the ADC drug product.

**Conjugation chemistry supported:**
- Cysteine thiol-maleimide (interchain disulfide reduction)
- Lysine NHS-ester / hydrazone
- Site-specific (engineered cysteines, unnatural amino acids)

**Linker biotransformations modeled as variable modifications:**

| Biotransformation | Mass shift | Residue |
|---|---|---|
| Maleimide ring hydrolysis | +18.011 Da | C |
| Succinimide ring-opening | +18.011 Da | C |
| Thioether → sulfoxide oxidation | +15.995 Da | C |
| Disulfide loss | −31.990 Da | C |
| Deamidation at conjugation site | +0.984 Da | N |

---

## File Structure

```
ADC_Peptide_Mapper_v0.8/
├── app.R                    ← Main Shiny application (7 tabs)
├── build_background_db.R    ← One-time database builder (UniProt + cRAP)
├── deploy.R                 ← shinyapps.io deployment script
├── DESCRIPTION              ← Package metadata (for rsconnect)
├── CITATION.cff             ← Citation metadata (CFF v1.2.0)
├── README.md                ← This file
├── R/
│   ├── digest.R             ← 11-enzyme digest engine
│   ├── modifications.R      ← PTM definitions, ADCDB payloads, linker biotransformations
│   ├── transitions.R        ← b/y ion series + CE calculation + DAR transitions
│   ├── isotopes.R           ← Averagine isotope envelope (Senko 1995)
│   ├── export.R             ← 6-platform instrument formatters
│   ├── uniqueness.R         ← Background proteome loading & uniqueness checking
│   ├── dar.R                ← DAR distribution modeling
│   └── msearch.R            ← MS Amanda + Tide engine detection, search, FDR, result parsing
├── tests/
│   ├── test_masses.R        ← 79 unit tests (100% pass rate)
│   └── benchmark_mass_accuracy.R  ← Sub-mDa accuracy benchmark (10/10 pass)
├── data/
│   ├── bg_human.rds         ← (generated by build_background_db.R)
│   ├── bg_monkey.rds        ← (generated by build_background_db.R)
│   ├── bg_rat.rds           ← (generated by build_background_db.R)
│   └── README.txt
└── www/
    └── custom.css           ← App styling + AI chat UI + sidebar citation styles
```

---

## Deployment (shinyapps.io)

```r
# One-time account setup
install.packages("rsconnect")
rsconnect::setAccountInfo(
  name   = "your-account-name",   # from shinyapps.io → Account → Profile
  token  = "YOUR_TOKEN",          # from shinyapps.io → Account → Tokens
  secret = "YOUR_SECRET"
)

# Deploy
source("deploy.R")
```

After deploying, set your Anthropic API key as a **shinyapps.io environment variable** (app Settings → Environment Variables → `ANTHROPIC_API_KEY`) — never hard-code it or commit it to git.

---

## Running Unit Tests

```r
setwd("path/to/ADC_Peptide_Mapper_v0.8")
source("R/digest.R")
source("R/modifications.R")
source("R/transitions.R")
source("R/isotopes.R")
source("R/dar.R")
source("tests/test_masses.R")         # 79 tests
source("tests/benchmark_mass_accuracy.R")  # 10 reference peptides
```

All 79 tests and all 10 mass accuracy benchmarks (≤ 0.05 mDa) should pass before deploying.

---

## Citation

If you use ADC Peptide Mapper in your research, please cite:

```
Wase, N. (2026). ADC Peptide Mapper (Version 0.8) [Software].
https://doi.org/10.5281/zenodo.20681412
```

If you use the MS/MS Search tab (Tab 6) with MS Amanda, also cite:

```
Dorfer V, et al. MS Amanda, a Universal Identification Algorithm Optimized
for High Accuracy Tandem Mass Spectra. J Proteome Res. 2014;13(8):3679-3684.
doi:10.1021/pr500202e
```

If using Tide/Crux, also cite:

```
McIlwain S, et al. Crux: Rapid Open Source Protein Tandem Mass Spectrometry
Analysis. J Proteome Res. 2014;13(10):4488-4491. doi:10.1021/pr500741y
```

See `CITATION.cff` for full metadata including all 16 scientific references.

---

## License

MIT License — see `LICENSE` for details.

---

## AI Development Disclosure

Portions of this application were developed with the assistance of [Claude](https://www.anthropic.com) (Anthropic), an AI assistant. Specifically, AI assistance was used for code structure, UI layout, documentation, and iterative debugging during development of version 0.8.

All scientific logic — including mass accuracy models, enzymatic digestion rules, DAR distribution modeling, linker biotransformation definitions, uniqueness filtering against reference proteomes, and instrument-specific collision energy formulas — was authored, designed, and independently validated by Nishikant Wase. The 79-unit test suite and 10 mass accuracy benchmarks (≤0.05 mDa vs. Unimod reference values) confirm the correctness of the core mass engine.

This disclosure follows best practices for transparent AI-assisted software development. For research use only — users should independently validate all outputs against their own experimental data.

---

## Author

**Nishikant Wase, PhD** — [nishikant.wase@gmail.com](mailto:nishikant.wase@gmail.com)  
Portfolio: [nishi76.github.io](https://nishi76.github.io)  
DOI: [10.5281/zenodo.20681412](https://doi.org/10.5281/zenodo.20681412)

*For research use only. Monoisotopic masses throughout. Background databases sourced from UniProt Swiss-Prot reviewed proteomes (Human, Cynomolgus Monkey, Rat).*
