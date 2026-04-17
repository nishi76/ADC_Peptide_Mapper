# ADC Peptide Mapper v0.5

**4-tab version** — In-silico tryptic digest, uniqueness checking, and Skyline transition list export.

## Features
- FASTA upload (multi-chain ADC: HC + LC auto-detected)
- Trypsin digest (KP/RP rule, 0 missed cleavages, 6-30 AA)
- Fixed mod: Carbamidomethylation (CAM, +57.021 Da on C)
- Variable mods: Oxidation (M), Propionamide (C), NEM (C), Drug-linker payloads (MMAE/DM1/DXd/SN-38/custom)
- Special mods: Deamidation, Pyroglutamate, Acetylation, Phosphorylation, custom builder
- Uniqueness check vs pre-built Human / Cynomolgus Monkey / Rat backgrounds (MC=0/1/2)
- Skyline CSV export (MRM and HR modes) + Excel summary

## Setup

```r
install.packages(c("shiny","bs4Dash","DT","data.table","openxlsx",
                   "httr2","stringr","dplyr","shinycssloaders","shinyjs"))
source("build_background_db.R")   # one-time, ~10 min
shiny::runApp("app.R")
```

## Tabs
1. Input & Setup
2. Modifications
3. Peptide Results
4. Transition List

## Auto-Deploy (GitHub Actions)
See `.github/workflows/deploy.yml`. Add these three secrets to your GitHub repo:
- `SHINYAPPS_ACCOUNT` — your shinyapps.io username
- `SHINYAPPS_TOKEN`   — from shinyapps.io > Account > Tokens
- `SHINYAPPS_SECRET`  — from shinyapps.io > Account > Tokens
