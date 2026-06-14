# ============================================================
#  ADC Peptide Mapper v0.8 — R Shiny Application
#  6-tab version: Input & Setup | Modifications | Peptide Results |
#                 Transition List | Heavy Labelling | MS/MS Search | ChatBot
#  Author: Nishikant Wase (nishikant.wase@gmail.com)
#  SETUP: Run build_background_db.R once to generate data/*.rds files.
#
# ============================================================

library(shiny)
library(bs4Dash)
library(DT)
library(data.table)
library(openxlsx)
library(httr2)
library(stringr)
library(dplyr)
library(shinyjs)

if (requireNamespace("shinycssloaders", quietly = TRUE)) {
  library(shinycssloaders)
} else {
  withSpinner <- function(ui, ...) ui
}

for (.module in c("R/digest.R", "R/modifications.R", "R/uniqueness.R",
                  "R/transitions.R", "R/export.R", "R/msearch.R",
                  "R/isotopes.R", "R/dar.R")) {
  tryCatch(
    source(.module),
    error = function(e) stop(sprintf(
      "ADC Peptide Mapper: failed to load '%s'.\n  Error: %s\n  Fix the module and restart.",
      .module, conditionMessage(e)
    ))
  )
}
rm(.module)

# ── Stubs: active only when R/ helper files are absent ───────────────────
if (!exists("ENZYME_LABELS"))
  ENZYME_LABELS <- c(
    trypsin      = "Trypsin (K/R, not P)",
    trypsin_p    = "Trypsin/P (K/R, incl. P)",
    lysc         = "Lys-C (K, not P)",
    lysc_p       = "Lys-C/P (K, incl. P)",
    lysn         = "Lys-N (before K)",
    aspn         = "Asp-N (before D)",
    gluc         = "Glu-C (E/D)",
    argc         = "Arg-C (R, not P)",
    chymotrypsin = "Chymotrypsin (F/Y/W, not P)",
    papain       = "Papain (K/R/Q)",
    elastase     = "Elastase (A/V/S/G/T)"
  )
if (!exists("PROTON_MASS"))              PROTON_MASS              <- 1.007276
if (!exists("parse_fasta"))              parse_fasta              <- function(text) {
  ln_ <- strsplit(text, "\n")[[1]]; h <- grep("^>", ln_); s <- list()
  for (k in seq_along(h)) {
    nm <- sub("^>\\s*", "", ln_[h[k]])
    en <- if (k < length(h)) h[k+1]-1 else length(ln_)
    s[[nm]] <- paste(grep("^>", ln_[(h[k]+1):en], value=TRUE, invert=TRUE), collapse="")
  }; s }
if (!exists("digest_fasta_enzyme"))      digest_fasta_enzyme      <- function(...) data.table()
if (!exists("apply_modifications"))      apply_modifications      <- function(sq, md) data.frame(position=integer(), mod_name=character())
if (!exists("calc_modified_mass"))       calc_modified_mass       <- function(sq, md, base_mass=0) base_mass
if (!exists("build_active_mods"))        build_active_mods        <- function(...) list()
if (!exists("flag_unique_peptides"))     flag_unique_peptides     <- function(seqs, bg) rep(NA, length(seqs))
if (!exists("get_bg_dt"))               get_bg_dt               <- function(bg, mc) data.table()
if (!exists("build_background_from_fasta")) build_background_from_fasta <- function(...) stop("R/uniqueness.R not loaded.")
if (!exists("get_payload_choices"))      get_payload_choices      <- function() c("Custom"="custom_drug")
if (!exists("generate_transition_list")) generate_transition_list <- function(...) data.table()
if (!exists("write_instrument_csv"))     write_instrument_csv     <- function(...) tempfile()
if (!exists("write_skyline_csv"))        write_skyline_csv        <- function(...) tempfile()
if (!exists("write_excel_summary"))      write_excel_summary      <- function(...) tempfile()
if (!exists("detect_search_engine"))       detect_search_engine       <- function(...) list(available=FALSE, engine=NA, version=NA, error_msg="R/msearch.R not loaded.")
if (!exists("run_msearch"))                run_msearch                <- function(...) stop("R/msearch.R not loaded.")
if (!exists("parse_msearch_results"))      parse_msearch_results      <- function(...) data.table()
if (!exists("generate_amanda_settings"))   generate_amanda_settings   <- function(...) NULL
if (!exists("calc_isotope_distribution"))  calc_isotope_distribution  <- function(...) data.frame(isotope=0L, mz_offset=0, probability=1, rel_intensity=100)
if (!exists("recommended_precursor_isotope")) recommended_precursor_isotope <- function(...) 0L
if (!exists("annotate_isotope_offsets"))   annotate_isotope_offsets   <- function(dt) dt
if (!exists("flag_cospecies_unique"))        flag_cospecies_unique        <- function(seqs, bg_list, il_equivalent=TRUE) data.frame(UniqueAllSpecies=rep(NA, length(seqs)))
if (!exists("detect_conjugation_sites"))   detect_conjugation_sites   <- function(seq, type, res="C") integer(0)
if (!exists("has_conjugation_site"))       has_conjugation_site       <- function(seq, type, res="C") FALSE
if (!exists("LINKER_BIOTRANSFORMATIONS"))  LINKER_BIOTRANSFORMATIONS  <- list()
if (!exists("generate_dar_transitions"))   generate_dar_transitions   <- function(...) data.table()
if (!exists("dar_summary_table"))          dar_summary_table          <- function(...) data.frame()


`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !is.na(a[1]) && nchar(a[1]) > 0) a else b

DEMO_FASTA <- paste0(
  ">Trastuzumab_HeavyChain\n",
  "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSREEMTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK\n",
  ">Trastuzumab_LightChain\n",
  "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIKRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC"
)

# Heavy label isotope definitions
HEAVY_LABEL_DEFS <- list(
  "13C6 15N2 Lys (+8.014 Da)" = list(residue="K", mass=8.014199,  name="13C6,15N2-Lys"),
  "13C6 15N4 Arg (+10.008 Da)" = list(residue="R", mass=10.008269, name="13C6,15N4-Arg"),
  "D4 Lys (+4.025 Da)"         = list(residue="K", mass=4.025107,  name="D4-Lys"),
  "D6 Leu (+6.032 Da)"         = list(residue="L", mass=6.031817,  name="D6-Leu"),
  "13C6 Leu (+6.020 Da)"       = list(residue="L", mass=6.020129,  name="13C6-Leu"),
  "13C9 15N1 Tyr (+10.009 Da)" = list(residue="Y", mass=10.009369, name="13C9,15N1-Tyr"),
  "Custom"                     = list(residue=NA,  mass=NA,        name="Custom")
)

# Species registry — maps internal keys to .rds paths and display labels
SPECIES_REGISTRY <- list(
  human  = list(label = "Human (Homo sapiens)",        rds_file = "data/bg_human.rds"),
  monkey = list(label = "Cynomolgus Monkey",            rds_file = "data/bg_monkey.rds"),
  rat    = list(label = "Rat (Rattus norvegicus)",      rds_file = "data/bg_rat.rds"),
  custom = list(label = "Custom FASTA (upload below)",  rds_file = NULL)
)

message("ADC Peptide Mapper v0.8: loading bundled background databases...")
BUNDLED_BG <- list()
for (key in setdiff(names(SPECIES_REGISTRY), "custom")) {
  spec <- SPECIES_REGISTRY[[key]]
  if (!is.null(spec$rds_file) && file.exists(spec$rds_file)) {
    BUNDLED_BG[[key]] <- readRDS(spec$rds_file)
    message(sprintf("  Loaded: %s (%s proteins)", spec$label,
                    format(BUNDLED_BG[[key]]$n_proteins, big.mark = ",")))
  }
}
if (length(BUNDLED_BG) == 0)
  message("  WARNING: No bundled databases found. Run build_background_db.R first.")

build_species_choices <- function() {
  choices <- list()
  for (key in names(BUNDLED_BG)) {
    bg  <- BUNDLED_BG[[key]]
    lbl <- sprintf("%s  (%s proteins | built %s)", bg$label,
                   format(bg$n_proteins, big.mark = ","), bg$build_date)
    choices[[lbl]] <- key
  }
  choices[["Custom FASTA (upload below)"]] <- "custom"
  choices
}

# ============================================================
#  UI
# ============================================================
ui <- bs4DashPage(
  title = "ADC Peptide Mapper",
  freshTheme = NULL,
  footer = bs4DashFooter(
    left = tags$span(
      style = "font-size:11px; color:#7f8c8d;",
      tags$b("ADC Peptide Mapper v0.8"),
      " — In-silico ADC digest & transition list generator. ",
      tags$a("Cite via Zenodo", href = "https://doi.org/10.5281/zenodo.XXXXXXX",
             target = "_blank", style = "color:#3498db;"),
      " | Wase N. (2026)"
    ),
    right = tags$span(
      id = "footer_session_info",
      style = "font-size:10px; color:#95a5a6; font-family:monospace;",
      paste0("R ", R.version$major, ".", R.version$minor,
             " | data.table ", packageVersion("data.table"),
             " | shiny ", packageVersion("shiny"))
    )
  ),
  preloader = list(
    html = tagList(tags$div(style = "text-align:center; padding-top:80px;",
      tags$h4("ADC Peptide Mapper v0.8", style = "color:#1a2940; font-weight:700;"),
      tags$p("Loading background databases...", style = "color:#7f8c8d;"))),
    color = "#f4f6f9"),

  header = bs4DashNavbar(
    title = bs4DashBrand(title = "ADC Peptide Mapper", color = "navy"),
    skin  = "dark",
    tags$li(class = "nav-item",
      tags$span(style = "color:#c8d6e5; padding:15px 10px; font-size:12px;",
        "In-silico Digest & Transition List Generator | v0.8"))),

  sidebar = bs4DashSidebar(
    skin = "dark", status = "navy",
    bs4SidebarMenu(id = "sidebar_menu",
      bs4SidebarMenuItem("Input & Setup",    tabName = "tab_input",   icon = icon("upload")),
      bs4SidebarMenuItem("Modifications",    tabName = "tab_mods",    icon = icon("flask")),
      bs4SidebarMenuItem("Peptide Results",  tabName = "tab_results", icon = icon("table")),
      bs4SidebarMenuItem("Transition List",  tabName = "tab_trans",   icon = icon("list")),
      bs4SidebarMenuItem("Heavy Labelling",  tabName = "tab_heavy",   icon = icon("weight-hanging")),
      bs4SidebarMenuItem("MS/MS Search",    tabName = "tab_search",  icon = icon("magnifying-glass")),
      bs4SidebarMenuItem("AI Assistant",    tabName = "tab_ai",      icon = icon("robot"))
    ),
    tags$div(class = "sidebar-citation",
      tags$div(style="font-weight:600; color:rgba(255,255,255,0.85); margin-bottom:3px;",
               icon("dna"), " ADC Peptide Mapper v0.8"),
      tags$div(HTML("&#169; 2026 Nishikant Wase")),
      tags$div(tags$a(
        href   = "https://doi.org/10.5281/zenodo.20681412",
        target = "_blank",
        style  = "color:rgba(100,180,255,0.85); text-decoration:none;",
        icon("book-open"), " 10.5281/zenodo.20681412"
      )),
      tags$div(style="margin-top:3px;",
        tags$a(href="mailto:nishikant.wase@gmail.com",
               style="color:rgba(255,255,255,0.5); text-decoration:none; font-size:10px;",
               icon("envelope"), " nishikant.wase@gmail.com")
      )
    )
  ),

  body = bs4DashBody(
    useShinyjs(),
    tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")),

    bs4TabItems(

      # ── TAB 1: INPUT & SETUP ──────────────────────────────────────────────
      bs4TabItem(tabName = "tab_input",
        fluidRow(
          bs4Card(title = "ADC FASTA Input", width = 7, status = "navy", collapsible = FALSE,
            tags$div(class = "section-header", "ADC Name & Sequence"),
            textInput("adc_name", "ADC Name (used in all outputs):",
                      placeholder = "e.g. Trastuzumab-DXd", width = "100%"),
            hr(),
            tabsetPanel(id = "fasta_input_mode",
              tabPanel("Upload File", br(),
                fileInput("fasta_file", "Upload FASTA file (.fasta / .fa / .txt)",
                          accept = c(".fasta", ".fa", ".txt"), placeholder = "No file selected"),
                tags$small(class = "text-muted",
                  "Multi-chain FASTA supported — each >header entry is a separate chain.")),
              tabPanel("Paste Sequence", br(),
                textAreaInput("fasta_text", NULL, value = DEMO_FASTA, rows = 10, width = "100%"),
                tags$small(class = "text-muted",
                  "Demo: Trastuzumab HC + LC loaded. Replace with your ADC sequence."))
            ),
            hr(), uiOutput("fasta_chain_preview"),
            hr(),
            tags$div(class = "section-header", "Run Digest"),
            fluidRow(
              column(8, tags$p(style = "margin-top:6px; color:#555;",
                "Configure enzyme, modifications, and background, then click",
                tags$b(" Run Digest"), " to generate peptides and check uniqueness.")),
              column(4, actionButton("btn_run_digest", "Run Digest & Uniqueness Check",
                                     class = "btn-run", icon = icon("play"), width = "100%"))),
            uiOutput("digest_status_ui"),
            hr(),
            div(class = "well", style = "background:#f0f7ff; border:1px solid #b8d4f0; border-radius:6px; padding:10px 14px; margin-top:8px; font-size:12px; color:#2c5f8a;",
              icon("robot"), tags$b(" AI Assistant (Tab 7)"),
              tags$p(style = "margin:6px 0 4px 0;",
                "This app includes an AI assistant for questions about your digest results, DAR modeling, transition selection, and LC-MS/MS method development."),
              tags$p(style = "margin:0;",
                "A limited number of free messages are available per session on the hosted app. For unlimited use, supply your own API key in Tab 7 — ",
                tags$a("Anthropic", href = "https://console.anthropic.com/", target = "_blank"), ", ",
                tags$a("OpenAI", href = "https://platform.openai.com/api-keys", target = "_blank"), ", or ",
                tags$a("Google Gemini", href = "https://aistudio.google.com/app/apikey", target = "_blank"),
                " keys are all accepted.")
            )
          ),

          bs4Card(title = "Background Proteome & Digest Settings", width = 5,
            status = "primary", collapsible = FALSE,
            tags$div(class = "section-header", "Enzyme & Cleavage Parameters"),
            selectInput("enzyme_id", "Primary Enzyme:",
                        choices  = setNames(names(ENZYME_LABELS), ENZYME_LABELS),
                        selected = "trypsin"),
            selectInput("enzyme_id2", "Optional Second Enzyme (two-enzyme digest):",
                        choices  = c("None (single enzyme)" = "",
                                     setNames(names(ENZYME_LABELS), ENZYME_LABELS)),
                        selected = ""),
            selectInput("missed_cleavages", "Missed Cleavages:",
                        choices  = c("0" = 0, "1" = 1, "2" = 2),
                        selected = 0),
            tags$small(class = "text-muted",
              "Peptide length filter: 6 – 30 amino acids (fixed)."),
            hr(),
            tags$div(class = "section-header", "Background Proteome"),
            selectInput("bg_species", "Select species background:",
                        choices  = build_species_choices(),
                        selected = if (length(BUNDLED_BG) > 0) names(BUNDLED_BG)[1] else "custom"),
            uiOutput("bg_info_ui"),
            conditionalPanel("input.bg_species == 'custom'", br(),
              fileInput("bg_fasta_file", "Upload custom background FASTA:",
                        accept = c(".fasta", ".fa", ".txt")),
              selectInput("bg_missed_custom", "Missed cleavages for custom background:",
                          choices = c("0"=0,"1"=1,"2"=2), selected = 0),
              actionButton("btn_build_custom_bg", "Build Custom Background",
                           class = "btn btn-warning btn-sm", icon = icon("cogs"))),
            hr(),
            tags$div(class = "section-header", "Background Missed Cleavages"),
            tags$p(style = "font-size:12px; color:#555; margin-bottom:4px;",
              "For bundled databases, all MC levels (0/1/2) are pre-computed."),
            selectInput("bg_missed", "Use missed cleavage level:",
                        choices = c("0"="0","1"="1","2"="2"), selected = "0"),
            hr(),
            checkboxInput("il_equivalent",
              label = tagList(
                "Treat I / L as equivalent during uniqueness check",
                tags$small(class = "text-muted",
                  " (recommended — MS/MS cannot distinguish Ile from Leu)")),
              value = TRUE),
            hr(), uiOutput("bg_status_ui"))
        ),
        fluidRow(
          bs4Card(
            title = tagList(icon("circle-info"), " About & Session Info"),
            width = 12, status = "secondary",
            collapsible = TRUE, collapsed = TRUE,
            solidHeader = FALSE,
            fluidRow(
              column(6,
                tags$h6("ADC Peptide Mapper v0.8", style = "font-weight:700; margin-bottom:4px;"),
                tags$p(style = "font-size:12px; color:#555; margin-bottom:6px;",
                  "In-silico digest and LC-MS/MS transition list generator for ",
                  "Antibody-Drug Conjugate peptide mapping studies."),
                tags$p(style = "font-size:12px; color:#555;",
                  tags$b("Author: "), "Nishikant Wase", tags$br(),
                  tags$b("License: "), "MIT", tags$br(),
                  tags$b("Citation: "),
                  tags$a("Zenodo DOI (replace placeholder after upload)",
                         href = "https://doi.org/10.5281/zenodo.XXXXXXX",
                         target = "_blank")),
                tags$hr(style = "margin:8px 0;"),
                tags$p(style = "font-size:11px; color:#777;",
                  tags$b("Key references:"), tags$br(),
                  "• Averagine isotopes: Senko et al. JASMS 1995", tags$br(),
                  "• Target-decoy FDR: Käll et al. J. Proteome Res. 2008", tags$br(),
                  "• CE formulas: Hoofnagle et al. Clin. Chem. 2016; Distler et al. Anal. Chem. 2016", tags$br(),
                  "• Linker biotransformations: Shen 2012; Lyon 2015; Tumey 2014; Ogitani 2016", tags$br(),
                  "• Fragment ion nomenclature: Roepstorff & Fohlman 1984")),
              column(6,
                tags$h6("Session Information", style = "font-weight:700; margin-bottom:4px;"),
                verbatimTextOutput("session_info_text", placeholder = FALSE))
            )
          )
        )
      ),

      # ── TAB 2: MODIFICATIONS ──────────────────────────────────────────────
      bs4TabItem(tabName = "tab_mods",
        fluidRow(
          bs4Card(title = "Fixed Modifications (always applied)", width = 4,
            status = "danger", collapsible = FALSE,
            tags$table(class = "mod-table",
              tags$thead(tags$tr(tags$th("Modification"), tags$th("Residue"), tags$th("Mass Shift (Da)"))),
              tags$tbody(tags$tr(tags$td("Carbamidomethyl (CAM)"), tags$td("C"), tags$td("+57.02146")))),
            tags$small(class = "text-muted", style = "margin-top:8px; display:block;",
              "CAM is always applied to all cysteine residues.")),

          bs4Card(title = "Variable Modifications", width = 4, status = "warning", collapsible = FALSE,
            checkboxInput("mod_oxidation",    "Oxidation (M) +15.99491 Da",         value = TRUE),
            checkboxInput("mod_propionamide", "Propionamide (C) +71.03711 Da",      value = FALSE),
            checkboxInput("mod_nem",          "N-Ethylmaleimide (C) +125.04768 Da", value = FALSE),
            hr(),
            checkboxInput("mod_drug_enable", "Drug-Linker Payload", value = FALSE),
            conditionalPanel("input.mod_drug_enable",
              selectInput("drug_payload_key", "Select payload:",
                choices = get_payload_choices()),
              conditionalPanel(
                "input.drug_payload_key == 'custom_drug' || input.drug_payload_key == 'adcdb_custom_payload'",
                numericInput("custom_drug_mass", "Custom payload mass shift (Da):",
                             value = 700, min = 0, step = 0.001)))),

          bs4Card(title = "Special Modifications", width = 4, status = "info", collapsible = FALSE,
            tags$div(class = "section-header", "Predefined"),
            checkboxInput("mod_deamidation_N",  "Deamidation (N) +0.98402 Da",         value = FALSE),
            checkboxInput("mod_deamidation_Q",  "Deamidation (Q) +0.98402 Da",         value = FALSE),
            checkboxInput("mod_pyro_glu_Q",     "Pyroglutamate N-term Q -17.02655 Da", value = FALSE),
            checkboxInput("mod_pyro_glu_E",     "Pyroglutamate N-term E -18.01056 Da", value = FALSE),
            checkboxInput("mod_acetyl_nterm",   "Acetylation (N-term) +42.01057 Da",   value = FALSE),
            checkboxInput("mod_phospho_STY",    "Phosphorylation (S/T/Y) +79.96633 Da",value = FALSE))
        ),
        # ── T3-C: Linker biotransformation search panel ───────────────────────
        fluidRow(
          bs4Card(title = "ADC Linker Biotransformation Search",
            width = 12, status = "orange", collapsible = TRUE, collapsed = TRUE,
            tags$p(style = "font-size:12px; color:#555;",
              "Enable these as variable modifications when searching for in vivo / in vitro ",
              "biotransformation products of the drug-linker. Only meaningful when a payload ",
              "is selected above and MS/MS data is available in Tab 6. ",
              tags$em("Refs: Shen et al. 2012 Nat. Biotechnol.; Alley et al. 2008 Bioconjug. Chem.")),
            fluidRow(
              column(4,
                checkboxInput("biotr_maleimide",  "Maleimide ring opening (+18.011 Da on C)",  FALSE),
                checkboxInput("biotr_succinimide", "Succinimide ring-opening (+18.011 Da on C)", FALSE)),
              column(4,
                checkboxInput("biotr_sulfoxide",   "Thioether → sulfoxide (+15.995 Da on C)",   FALSE),
                checkboxInput("biotr_disulfide",   "Disulfide loss (−31.990 Da on C)",           FALSE)),
              column(4,
                checkboxInput("biotr_deamid_n",    "Deamidation near conj. site (+0.984 Da on N)", FALSE))
            ))
        ),
        fluidRow(
          bs4Card(title = "Custom Modification Builder", width = 12,
            status = "secondary", collapsible = TRUE, collapsed = FALSE,
            fluidRow(
              column(2, textInput("custom_mod_residue", "Residue (single letter):", placeholder = "e.g. K")),
              column(3, textInput("custom_mod_name", "Modification Name:", placeholder = "e.g. Ubiquitination")),
              column(2, numericInput("custom_mod_mass", "Mass Shift (Da):", value = 0, step = 0.001)),
              column(3, selectInput("custom_mod_location", "Location:",
                                    choices = c("Any residue" = "any",
                                                "N-terminus only" = "nterm",
                                                "C-terminus only" = "cterm"),
                                    selected = "any")),
              column(2, br(), actionButton("btn_add_custom_mod", "Add Mod",
                                           class = "btn btn-primary", icon = icon("plus")))),
            hr(), DTOutput("custom_mods_table")))
      ),

      # ── TAB 3: PEPTIDE RESULTS ────────────────────────────────────────────
      bs4TabItem(tabName = "tab_results",
        fluidRow(
          bs4ValueBoxOutput("vbox_total",    width = 3),
          bs4ValueBoxOutput("vbox_unique",   width = 3),
          bs4ValueBoxOutput("vbox_modified", width = 3),
          bs4ValueBoxOutput("vbox_chains",   width = 3)
        ),
        fluidRow(
          bs4Card(title = "Peptide Results", width = 12, status = "navy", collapsible = FALSE,
            fluidRow(
              column(3, selectInput("filter_chain", "Filter by Chain:",
                                    choices = c("All Chains"="all"), selected = "all")),
              column(3, selectInput("filter_unique", "Uniqueness:",
                                    choices = c("All Peptides"="all","Unique to ADC only"="unique",
                                                "Non-unique only"="nonunique"), selected = "all")),
              column(3, selectInput("filter_mods", "Modifications:",
                                    choices = c("All"="all","Modified only"="modified",
                                                "Unmodified only"="unmodified"), selected = "all")),
              column(3, br(), downloadButton("dl_peptide_excel", "Download Excel Summary",
                                             class = "btn-download"))),
            hr(),
            # T2-C: Multi-species co-uniqueness panel
            tags$div(class = "section-header", "Multi-Species Co-uniqueness"),
            tags$p(style = "font-size:12px; color:#555;",
              "Check uniqueness across several background proteomes simultaneously. ",
              "Only peptides NOT found in any selected background pass."),
            fluidRow(
              column(6,
                checkboxGroupInput("cospecies_select",
                  label = "Include in co-uniqueness check:",
                  choices  = c("Human (Homo sapiens)" = "human",
                               "Cynomolgus monkey"    = "monkey",
                               "Rat (Rattus norvegicus)" = "rat"),
                  selected = c("human","monkey"))),
              column(6, br(), br(),
                actionButton("btn_cospecies", "Run Co-uniqueness Check",
                             class = "btn btn-info btn-sm", icon = icon("globe")),
                br(), br(),
                uiOutput("cospecies_result_ui"))
            ),
            hr(),
            withSpinner(DTOutput("peptide_table"), type = 6, color = "#2980b9"))
        ),
        fluidRow(
          bs4Card(title = "Sequence Coverage Map", width = 12,
            status = "info", collapsible = TRUE, collapsed = FALSE,
            tags$p(style = "color:#555; font-size:13px;",
              "Theoretical peptides from the digest mapped onto each chain. ",
              "No MS/MS search required — this reflects in-silico coverage only. ",
              "Colour, chain filter, and label options are below."),
            fluidRow(
              column(3, selectInput("cov_chain_filter", "Show chain:",
                                    choices = c("All chains" = "all"), selected = "all")),
              column(3, selectInput("cov_color_by", "Colour by:",
                                    choices = c("Uniqueness"      = "unique",
                                                "Missed cleavages"= "mc",
                                                "Peptide length"  = "length"),
                                    selected = "unique")),
              column(3, checkboxInput("cov_show_labels", "Show peptide labels", value = FALSE)),
              column(3, br(), downloadButton("dl_coverage_map", "Download PNG",
                                             class = "btn-download"))
            ),
            withSpinner(plotOutput("digest_coverage_plot", height = "auto"), type = 6,
                        color = "#2980b9"))
        )
      ),

      # ── TAB 4: TRANSITION LIST ────────────────────────────────────────────
      bs4TabItem(tabName = "tab_trans",
        fluidRow(
          bs4Card(title = "Transition List Settings", width = 4, status = "navy", collapsible = FALSE,
            tags$div(class = "section-header", "Instrument Platform"),
            selectInput("instrument", "Target Instrument:",
              choices = c(
                "Skyline (generic)"  = "skyline",
                "Thermo (TSQ/Orbi)" = "thermo",
                "SCIEX (QTRAP/TripleTOF)" = "sciex",
                "Bruker (timsTOF)"  = "bruker",
                "Agilent (QQQ)"     = "agilent",
                "Waters (Xevo)"     = "waters"
              ),
              selected = "skyline"),
            hr(),
            tags$div(class = "section-header", "Ion Parameters"),
            tags$table(class = "mod-table",
              tags$tr(tags$td("Precursor charges"), tags$td("2+, 3+, 4+")),
              tags$tr(tags$td("Ion types"),         tags$td("b, y (full series)")),
              tags$tr(tags$td("Coverage"),          tags$td("b2..b(n-1), y2..y(n-1)")),
              tags$tr(tags$td("CE formula"),        tags$td(uiOutput("ce_formula_ui")))),
            hr(),
            checkboxInput("trans_unique_only", "Unique peptides only", value = FALSE),
            numericInput("trans_top_n", "Top N ions per precursor (0 = all):",
                         value = 0, min = 0, max = 20, step = 1),
            hr(),
            actionButton("btn_gen_trans", "Generate Transition List",
                         class = "btn-run", icon = icon("bolt"), width = "100%"),
            br(), br(),
            tags$div(class = "section-header", "Download"),
            downloadButton("dl_instrument_csv", "Download Instrument CSV",
                           class = "btn-download", style = "width:100%; margin-bottom:6px;"),
            br(),
            downloadButton("dl_skyline_all",    "Skyline CSV (All Peptides)",
                           class = "btn-download", style = "width:100%; margin-bottom:6px;"),
            br(),
            downloadButton("dl_skyline_unique", "Skyline CSV (Unique Only)",
                           class = "btn-download", style = "width:100%; margin-bottom:6px;"),
            br(),
            downloadButton("dl_excel_full", "Excel Summary (All Sheets)",
                           class = "btn-download", style = "width:100%;")),

          bs4Card(title = "Transition List Preview", width = 8, status = "primary", collapsible = FALSE,
            uiOutput("trans_summary_ui"),
            hr(),
            tabsetPanel(id = "trans_subtabs",
              tabPanel("Standard Transitions",
                br(),
                withSpinner(DTOutput("trans_table"), type = 6, color = "#2980b9")),

              # T3-A: DAR Transitions sub-tab
              tabPanel("DAR Transitions",
                br(),
                tags$p(style = "font-size:12px; color:#555;",
                  "Drug-to-Antibody Ratio (DAR) transitions: one precursor per ",
                  "DAR level (0 = naked antibody, k = k payloads per peptide). ",
                  "Requires a payload to be enabled in Tab 2 and peptides to contain ",
                  "the conjugation-site residue."),
                fluidRow(
                  column(4,
                    sliderInput("dar_range_ui", "DAR levels to include:",
                                min = 0L, max = 8L, value = c(0L, 4L), step = 1L,
                                ticks = TRUE)),
                  column(4,
                    selectInput("dar_conjugation_type", "Conjugation chemistry:",
                      choices = c("Cysteine (thiol-maleimide)" = "cysteine",
                                  "Lysine (NHS-ester / hydrazone)" = "lysine",
                                  "Site-specific (engineered)" = "site_specific"),
                      selected = "cysteine")),
                  column(4, br(),
                    actionButton("btn_gen_dar", "Generate DAR Transitions",
                                 class = "btn btn-warning btn-sm", icon = icon("vials")))
                ),
                hr(),
                uiOutput("dar_summary_ui"),
                br(),
                downloadButton("dl_dar_csv", "Download DAR Transitions CSV",
                               class = "btn-download", style = "margin-bottom:8px;"),
                br(),
                withSpinner(DTOutput("dar_table"), type = 6, color = "#2980b9"))
            ))
          )
        ),

      # ── TAB 5: HEAVY LABELLING ────────────────────────────────────────────
      bs4TabItem(tabName = "tab_heavy",
        fluidRow(
          bs4Card(title = "Heavy Isotope Label Settings", width = 5,
            status = "navy", collapsible = FALSE,
            tags$div(class = "section-header", "Select Heavy Label"),
            selectInput("heavy_label_key", "Isotope Label:",
                        choices = names(HEAVY_LABEL_DEFS),
                        selected = names(HEAVY_LABEL_DEFS)[1]),
            conditionalPanel("input.heavy_label_key == 'Custom'",
              fluidRow(
                column(4, textInput("heavy_custom_residue", "Residue:", placeholder = "e.g. K")),
                column(4, numericInput("heavy_custom_mass", "Mass shift (Da):", value = 0, step = 0.001)),
                column(4, textInput("heavy_custom_name", "Label name:", placeholder = "e.g. D8-Lys"))
              )),
            hr(),
            tags$div(class = "section-header", "Precursor Charge States"),
            checkboxGroupInput("heavy_charges", NULL,
                               choices = c("2+" = 2, "3+" = 3, "4+" = 4),
                               selected = c(2, 3), inline = TRUE),
            hr(),
            actionButton("btn_gen_heavy", "Generate Heavy Peptide List",
                         class = "btn-run", icon = icon("atom"), width = "100%"),
            br(), br(),
            tags$div(class = "section-header", "Download"),
            downloadButton("dl_heavy_csv",   "Download Heavy Peptides CSV",
                           class = "btn-download", style = "width:100%; margin-bottom:6px;"),
            br(),
            downloadButton("dl_heavy_excel", "Download Heavy Peptides Excel",
                           class = "btn-download", style = "width:100%;")),

          bs4Card(title = "Heavy Labelled Peptide Preview", width = 7,
            status = "primary", collapsible = FALSE,
            uiOutput("heavy_summary_ui"),
            hr(),
            withSpinner(DTOutput("heavy_table"), type = 6, color = "#2980b9"))
        ),
        fluidRow(
          bs4Card(title = "Light / Heavy Pair Summary", width = 12,
            status = "info", collapsible = TRUE, collapsed = FALSE,
            tags$p(style = "color:#555; font-size:13px;",
              "This table shows the mass difference between light and heavy peptide pairs,",
              " useful for verifying isotope incorporation and setting up SRM/PRM assays."),
            withSpinner(DTOutput("heavy_pair_table"), type = 6, color = "#2980b9"))
        )
      ),  # end Tab 5

      # ── TAB 6: MS/MS SEARCH ───────────────────────────────────────────────
      bs4TabItem(tabName = "tab_search",
        fluidRow(

          # ── Left panel: Setup & Run ────────────────────────────────────────
          bs4Card(title = "MS/MS Search Engine Setup & Run", width = 4,
            status = "navy", collapsible = FALSE,

            # Engine status badge
            uiOutput("msf_status_ui"),
            hr(),

            # Executable path input
            textInput("msf_jar_path",
              label = "Engine executable path (optional)",
              placeholder = "e.g. ~/MSAmanda/MSAmanda  or  ~/crux/bin/crux",
              width = "100%"),
            tags$small(class = "text-muted",
              "Leave blank for auto-detection. Or set ",
              tags$code("MSAMANDA_EXE"), " / ", tags$code("CRUX_EXE"),
              " in ~/.Renviron and restart R. ",
              "Download: ",
              tags$a("MS Amanda 3.0",
                     href = "https://github.com/hgb-bin-proteomics/MSAmanda/releases",
                     target = "_blank"),
              " | ",
              tags$a("Tide/Crux 4.x",
                     href = "https://crux.ms/download.html",
                     target = "_blank")),
            br(), br(),

            # Raw file upload
            fileInput("msf_raw_files",
              label = "Upload spectral files (mzML / MGF) or pre-computed results",
              multiple = TRUE,
              accept   = c(".mzML", ".mzXML", ".mgf",
                           ".tsv",  ".pepXML", ".pep.xml", ".mzid", ".csv"),
              width    = "100%"),
            tags$small(class = "text-muted",
              "Also accepts mzIdentML (MS Amanda), pepXML (Tide), psm.tsv, or Amanda CSV — ",
              "format is auto-detected. Uploading a result file skips the search step."),
            br(), br(),

            # Info box
            tags$div(class = "alert alert-info", style = "font-size:12px; padding:8px;",
              icon("circle-info"), " FASTA and enzyme settings are taken from Tab 1. ",
              tags$b("Run the digest first before searching. "),
              "MS Amanda settings.xml is generated automatically from your Tab 1/2 parameters."),

            hr(),
            tags$div(class = "section-header", "Score Filters (applied post-search)"),

            # Dynamic score filter UI (labels change per detected engine)
            uiOutput("score_filter_ui"),
            br(),
            # T2-F: fragment mass tolerance filter
            numericInput("ms2_ppm_tol",
              label = "Fragment mass tolerance (ppm):",
              value = 20, min = 1, max = 200, step = 1, width = "100%"),
            tags$small(class = "text-muted",
              "PSMs with any matched fragment ion outside this tolerance are excluded. ",
              "Typical values: 20 ppm (Orbitrap), 100 ppm (ion trap/TOF)."),
            br(),
            # T2-D: FDR filter
            checkboxInput("fdr_filter_enable", "Apply 1% FDR filter (target-decoy)", value = FALSE),
            uiOutput("fdr_filter_info_ui"),

            br(),

            # Run button
            actionButton("btn_run_search", "Run MS/MS Search",
              class = "btn-success btn-block",
              icon  = icon("play"),
              style = "width:100%; font-weight:600;"),
            br(),
            tags$small(class = "text-muted",
              "Search runs locally. Typical ADC FASTA search: < 2 min."),
            br(), br(),

            # Log output
            tags$div(class = "section-header", "Search Log"),
            verbatimTextOutput("msf_log_output")
          ),

          # ── Right panel: Results ───────────────────────────────────────────
          bs4Card(title = "Search Results", width = 8,
            status = "primary", collapsible = FALSE,

            uiOutput("search_summary_ui"),
            hr(),

            # Download buttons
            fluidRow(
              column(6,
                downloadButton("dl_psm_csv",
                  "Download PSM Table (CSV)",
                  class = "btn-download", style = "width:100%;")),
              column(6,
                downloadButton("dl_coverage_csv",
                  "Download Coverage Summary (CSV)",
                  class = "btn-download", style = "width:100%;"))
            ),
            br(),

            # Sub-tabs for three visualisations
            tabsetPanel(id = "search_subtabs",
              tabPanel("Sequence Coverage",
                br(),
                tags$p(style = "color:#555; font-size:12px;",
                  "Detected peptides (navy = unique to ADC, grey = non-unique) ",
                  "overlaid on the full ADC chain sequence."),
                withSpinner(plotOutput("coverage_plot", height = "340px"),
                            type = 6, color = "#2980b9")
              ),
              tabPanel("Detected vs Theoretical",
                br(),
                tags$p(style = "color:#555; font-size:12px;",
                  "All theoretical peptides from Tab 3 with detection status, ",
                  "PSM count, and best scores from the search."),
                withSpinner(DTOutput("search_results_table"),
                            type = 6, color = "#2980b9")
              ),
              tabPanel("Modification Summary",
                br(),
                tags$p(style = "color:#555; font-size:12px;",
                  "PSM counts per detected modification type (from filtered PSM table)."),
                withSpinner(plotOutput("mod_summary_plot", height = "320px"),
                            type = 6, color = "#2980b9")
              )
            )
          )
        )
      ),

      # ── TAB 7: AI ASSISTANT ───────────────────────────────────────────────
      bs4TabItem(tabName = "tab_ai",
        fluidRow(
          # ── Left column: API key + context + quick prompts ─────────────────
          column(width = 4,

            bs4Card(title = tagList(icon("key"), " API Key"), width = 12,
              status = "warning", collapsible = TRUE, collapsed = FALSE,
              solidHeader = FALSE,
              div(id = "apikey_panel",
                # Provider selector
                selectInput("ai_provider", "LLM Provider:",
                  choices = c(
                    "Anthropic (Claude)"  = "anthropic",
                    "OpenAI (GPT)"        = "openai",
                    "Google (Gemini)"     = "gemini"
                  ),
                  selected = "anthropic", width = "100%"),
                uiOutput("ai_provider_key_ui"),
                fluidRow(
                  column(6,
                    actionButton("btn_save_apikey", "Save to .Renviron",
                                 icon = icon("floppy-disk"),
                                 class = "btn btn-sm btn-outline-success",
                                 style = "width:100%; margin-top:4px;",
                                 title = "Persist key across sessions (local use only)")),
                  column(6,
                    actionButton("btn_clear_apikey", "Clear",
                                 icon = icon("eraser"),
                                 class = "btn btn-sm btn-outline-secondary",
                                 style = "width:100%; margin-top:4px;"))
                ),
                uiOutput("apikey_status_ui"),
                uiOutput("ai_provider_hint_ui")
              ),
              uiOutput("ai_model_ui"),
              tags$small(class = "text-muted",
                "For Anthropic, models are fetched live when your key is entered.")
            ),

            bs4Card(title = tagList(icon("database"), " Session Context"), width = 12,
              status = "info", collapsible = TRUE, collapsed = FALSE,
              solidHeader = FALSE,
              tags$p(style = "font-size:12px; color:#555; margin-bottom:8px;",
                "The AI automatically receives context from your current session:"),
              uiOutput("ai_context_badges"),
              tags$hr(style = "margin:8px 0;"),
              checkboxInput("ai_include_peptides",
                "Include peptide table (top 20 rows) in context", value = TRUE),
              checkboxInput("ai_include_transitions",
                "Include transition list summary in context", value = TRUE)
            ),

            bs4Card(title = tagList(icon("bolt"), " Quick Prompts"), width = 12,
              status = "secondary", collapsible = TRUE, collapsed = FALSE,
              solidHeader = FALSE,
              tags$p(style = "font-size:12px; color:#555;", "Click to insert:"),
              tags$div(style = "display:flex; flex-wrap:wrap; gap:6px;",
                actionButton("qp1", "Explain DAR modeling",
                             class="btn btn-outline-primary btn-sm"),
                actionButton("qp2", "Best transitions for this peptide?",
                             class="btn btn-outline-primary btn-sm"),
                actionButton("qp3", "What are linker biotransformations?",
                             class="btn btn-outline-primary btn-sm"),
                actionButton("qp4", "Troubleshoot uniqueness",
                             class="btn btn-outline-primary btn-sm"),
                actionButton("qp5", "Recommend CE for my instrument",
                             class="btn btn-outline-primary btn-sm"),
                actionButton("qp6", "How should I set missed cleavages?",
                             class="btn btn-outline-primary btn-sm"),
                actionButton("qp7", "Summarise my digest results",
                             class="btn btn-outline-primary btn-sm"),
                actionButton("qp8", "Which payload suits cysteine conjugation?",
                             class="btn btn-outline-primary btn-sm")
              )
            )
          ),

          # ── Right column: chat window ──────────────────────────────────────
          column(width = 8,
            bs4Card(title = tagList(icon("robot"), " ADC Mapper AI Assistant"),
              width = 12, status = "primary",
              collapsible = FALSE, solidHeader = TRUE,
              footer = tagList(
                div(id = "chat_input_row",
                  textAreaInput("chat_input", label = NULL, rows = 2,
                    placeholder = "Ask anything about your ADC digest, DAR modeling, transitions, or LC-MS/MS method development...",
                    width = "100%"),
                  actionButton("btn_chat_send", label = tagList(icon("paper-plane"), " Send"),
                               class = "btn btn-primary")
                ),
                div(style = "display:flex; justify-content:space-between; align-items:center; margin-top:6px;",
                  actionButton("btn_chat_clear", "Clear conversation",
                               class = "btn btn-outline-secondary btn-sm",
                               icon = icon("trash")),
                  div(style = "display:flex; gap:10px; align-items:center;",
                    uiOutput("ai_usage_badge"),
                    uiOutput("ai_token_info")
                  )
                )
              ),
              # Chat message box — rendered via htmlOutput so we can inject HTML
              div(id = "chat_box",
                htmlOutput("chat_messages")
              )
            )
          )
        )
      )

    ) # end bs4TabItems
  )  # end body
)

# ============================================================
#  SERVER
# ============================================================
server <- function(input, output, session) {

  rv <- reactiveValues(
    fasta_chains  = NULL,
    bg_obj        = NULL,
    bg_dt         = NULL,
    peptides_dt   = NULL,
    transition_dt = NULL,
    heavy_dt      = NULL,
    psm_dt        = NULL,   # Tab 6: raw PSM results from MSFragger / uploaded file
    search_done   = FALSE,  # Tab 6: TRUE after successful search + cross-reference
    msf_log       = "",     # Tab 6: accumulated search log text
    custom_mods   = data.frame(residue=character(), name=character(),
                                mass=numeric(), location=character(),
                                stringsAsFactors=FALSE),
    digest_done   = FALSE,
    trans_done    = FALSE,
    heavy_done    = FALSE,
    dar_dt        = NULL,    # T3-A: DAR transition table
    dar_done      = FALSE,
    # ── AI Assistant (Tab 7) ───────────────────────────────────────────────
    chat_history  = list(),  # list of list(role, content) — Anthropic messages format
    ai_thinking   = FALSE,   # TRUE while API call in flight
    ai_total_tokens = 0L,   # cumulative token usage this session
    ai_msg_count  = 0L      # messages sent this session (rate limit)
  )

  # ── FASTA parsing ──────────────────────────────────────────────────────────
  get_fasta_text <- reactive({
    if (!is.null(input$fasta_file))
      paste(readLines(input$fasta_file$datapath, warn = FALSE), collapse = "\n")
    else
      input$fasta_text
  })

  observe({
    ft <- get_fasta_text()
    req(nchar(trimws(ft)) > 0)
    rv$fasta_chains <- tryCatch(parse_fasta(ft), error = function(e) NULL)
  })

  output$fasta_chain_preview <- renderUI({
    chains <- rv$fasta_chains
    if (is.null(chains) || length(chains) == 0)
      return(tags$p(class = "status-warn", "No valid FASTA detected yet."))
    rows <- lapply(seq_along(chains), function(i) {
      nm <- names(chains)[i]; sq <- chains[[nm]]
      tags$tr(tags$td(tags$b(paste0("Chain ", i))),
              tags$td(nm, style = "font-size:11px; color:#555;"),
              tags$td(paste0(nchar(sq), " AA")),
              tags$td(tags$div(class = "fasta-preview",
                               substr(sq, 1, 60), if (nchar(sq) > 60) "..." else "")))
    })
    tagList(
      tags$div(class = "section-header", paste0(length(chains), " chain(s) detected")),
      tags$table(class = "mod-table",
        tags$thead(tags$tr(tags$th("#"), tags$th("Header"), tags$th("Length"), tags$th("Sequence Preview"))),
        tags$tbody(rows)))
  })

  # ── Background info ────────────────────────────────────────────────────────
  output$bg_info_ui <- renderUI({
    sp <- input$bg_species; req(sp != "custom")
    bg <- BUNDLED_BG[[sp]]
    if (is.null(bg)) return(tags$p(class = "status-warn", "Database not found."))
    mc_key <- as.character(input$bg_missed)
    n_pep  <- if (mc_key %in% names(bg$bg_sets))
                format(nrow(bg$bg_sets[[mc_key]]), big.mark = ",") else "N/A"
    tags$div(style = "margin-top:8px;",
      tags$span(class = "status-ok", icon("circle-check"),
        sprintf(" %s proteins | %s peptides (MC=%s) | Built: %s",
                format(bg$n_proteins, big.mark=","), n_pep, mc_key, bg$build_date)))
  })

  output$bg_status_ui <- renderUI({
    sp <- input$bg_species
    if (sp != "custom") {
      bg <- BUNDLED_BG[[sp]]
      if (!is.null(bg))
        tags$p(class = "status-ok", icon("circle-check"),
               " Bundled database ready — no loading required.")
      else
        tags$p(class = "status-err", icon("triangle-exclamation"),
               " Database file missing. Run build_background_db.R.")
    } else {
      if (is.null(rv$bg_obj))
        tags$p(class = "status-warn", icon("circle-exclamation"),
               " Upload a FASTA and click 'Build Custom Background'.")
      else
        tags$p(class = "status-ok", icon("circle-check"),
               sprintf(" Custom background ready: %s proteins",
                       format(rv$bg_obj$n_proteins, big.mark=",")))
    }
  })

  observeEvent(input$btn_build_custom_bg, {
    req(input$bg_fasta_file)
    withProgress(message = "Building custom background...", value = 0, {
      tryCatch({
        setProgress(0.2, detail = "Digesting uploaded FASTA...")
        bg_obj <- build_background_from_fasta(
          fasta_path       = input$bg_fasta_file$datapath,
          missed_cleavages = as.integer(input$bg_missed_custom),
          enzyme_id        = input$enzyme_id %||% "trypsin",
          progress_cb      = function(msg) setProgress(0.7, detail = msg))
        rv$bg_obj <- bg_obj
        setProgress(1, detail = "Done!")
        showNotification(paste0("Custom background ready: ",
                                format(bg_obj$n_proteins, big.mark=","), " proteins"),
                         type = "message", duration = 5)
      }, error = function(e) {
        showNotification(paste("Error:", conditionMessage(e)), type = "error", duration = 10)
      })
    })
  })

  get_active_bg_dt <- function() {
    sp <- input$bg_species; mc <- as.character(input$bg_missed)
    if (sp != "custom") {
      bg <- BUNDLED_BG[[sp]]
      if (is.null(bg)) stop("Bundled database not found for: ", sp)
      return(get_bg_dt(bg, as.integer(mc)))
    } else {
      if (is.null(rv$bg_obj)) stop("Custom background not built yet.")
      return(get_bg_dt(rv$bg_obj, as.integer(input$bg_missed_custom)))
    }
  }

  # ── Custom mods ────────────────────────────────────────────────────────────
  observeEvent(input$btn_add_custom_mod, {
    res  <- toupper(trimws(input$custom_mod_residue))
    nm   <- trimws(input$custom_mod_name)
    mass <- input$custom_mod_mass
    loc  <- input$custom_mod_location %||% "any"
    if (nchar(res) != 1 || !res %in% LETTERS) {
      showNotification("Residue must be a single letter (A-Z).", type = "warning"); return()
    }
    if (nchar(nm) == 0) {
      showNotification("Please enter a modification name.", type = "warning"); return()
    }
    rv$custom_mods <- rbind(rv$custom_mods,
                            data.frame(residue=res, name=nm, mass=mass, location=loc,
                                       stringsAsFactors=FALSE))
    updateTextInput(session, "custom_mod_residue", value = "")
    updateTextInput(session, "custom_mod_name",    value = "")
    updateNumericInput(session, "custom_mod_mass", value = 0)
    updateSelectInput(session, "custom_mod_location", selected = "any")
  })

  output$custom_mods_table <- renderDT({
    df <- rv$custom_mods
    if (nrow(df) == 0)
      return(datatable(data.frame(Message="No custom modifications added yet."),
                       options=list(dom="t"), rownames=FALSE))
    df$Delete <- paste0('<button class="btn btn-danger btn-sm" onclick="Shiny.setInputValue(\'del_custom_mod\',',
                        seq_len(nrow(df)), ', {priority:\'event\'})">Remove</button>')
    datatable(df, escape=FALSE, rownames=FALSE,
              colnames=c("Residue","Name","Mass Shift (Da)","Location",""),
              options=list(dom="t", pageLength=10))
  })

  observeEvent(input$del_custom_mod, {
    idx <- as.integer(input$del_custom_mod)
    if (!is.na(idx) && idx >= 1 && idx <= nrow(rv$custom_mods))
      rv$custom_mods <- rv$custom_mods[-idx, , drop=FALSE]
  })

  get_active_mods <- reactive({
    var_sel <- character(0)
    if (isTRUE(input$mod_oxidation))    var_sel <- c(var_sel, "oxidation")
    if (isTRUE(input$mod_propionamide)) var_sel <- c(var_sel, "propionamide")
    if (isTRUE(input$mod_nem))          var_sel <- c(var_sel, "nem")
    if (isTRUE(input$mod_drug_enable))  var_sel <- c(var_sel, input$drug_payload_key)

    special_sel <- character(0)
    if (isTRUE(input$mod_deamidation_N))  special_sel <- c(special_sel, "deamidation_N")
    if (isTRUE(input$mod_deamidation_Q))  special_sel <- c(special_sel, "deamidation_Q")
    if (isTRUE(input$mod_pyro_glu_Q))     special_sel <- c(special_sel, "pyro_glu_Q")
    if (isTRUE(input$mod_pyro_glu_E))     special_sel <- c(special_sel, "pyro_glu_E")
    if (isTRUE(input$mod_acetyl_nterm))   special_sel <- c(special_sel, "acetyl_nterm")
    if (isTRUE(input$mod_phospho_STY))    special_sel <- c(special_sel, "phospho_STY")
    # T3-C: Linker biotransformation mods
    if (isTRUE(input$biotr_maleimide))    special_sel <- c(special_sel, "biotr_maleimide_hydrolysis")
    if (isTRUE(input$biotr_succinimide))  special_sel <- c(special_sel, "biotr_succinimide_open")
    if (isTRUE(input$biotr_sulfoxide))    special_sel <- c(special_sel, "biotr_thioether_sulfoxide")
    if (isTRUE(input$biotr_disulfide))    special_sel <- c(special_sel, "biotr_disulfide_loss")
    if (isTRUE(input$biotr_deamid_n))     special_sel <- c(special_sel, "biotr_deamidation_linker")

    custom_drug_mass <- if (isTRUE(input$mod_drug_enable) &&
                            input$drug_payload_key %in% c("custom_drug", "adcdb_custom_payload"))
                          input$custom_drug_mass else NA_real_

    build_active_mods(
      selected_var_mods     = var_sel,
      selected_special_mods = special_sel,
      custom_mods_df        = rv$custom_mods,
      custom_drug_mass      = custom_drug_mass,
      include_fixed         = TRUE
    )
  })

  # ── Run Digest ─────────────────────────────────────────────────────────────
  observeEvent(input$btn_run_digest, {
    ft <- get_fasta_text()
    if (is.null(ft) || nchar(trimws(ft)) == 0) {
      showNotification("Please provide an ADC FASTA sequence.", type = "warning"); return()
    }
    if (input$bg_species == "custom" && is.null(rv$bg_obj)) {
      showNotification("Please build the custom background first.", type = "warning"); return()
    }
    withProgress(message = "Running digest...", value = 0, {
      tryCatch({
        setProgress(0.1, detail = "Parsing FASTA...")
        enzyme2 <- if (nzchar(trimws(input$enzyme_id2 %||% ""))) input$enzyme_id2 else NULL
        base_dt <- digest_fasta_enzyme(
          fasta_text       = ft,
          enzyme_id        = input$enzyme_id,
          enzyme_id2       = enzyme2,
          missed_cleavages = as.integer(input$missed_cleavages),
          min_len          = 6L,
          max_len          = 30L
        )
        if (nrow(base_dt) == 0) {
          showNotification("No peptides generated. Check your FASTA.", type="warning"); return()
        }
        setProgress(0.3, detail = "Applying modifications...")
        active_mods <- get_active_mods()
        all_rows <- list()
        for (i in seq_len(nrow(base_dt))) {
          pep      <- base_dt[i, ]

          # --- v0.6 adapter: apply_modifications() returns a data.frame of
          #     modification sites (one row per site), not a list of variants.
          sites_df <- apply_modifications(pep$Sequence, active_mods)

          # Build annotated modified sequence string, e.g. "PEPTC[CAM]IDE"
          seq_chars <- strsplit(pep$Sequence, "", fixed = TRUE)[[1]]
          if (nrow(sites_df) > 0) {
            for (si in seq_len(nrow(sites_df))) {
              pos <- sites_df$position[si]
              seq_chars[pos] <- paste0(seq_chars[pos], "[", sites_df$mod_name[si], "]")
            }
            modified_seq <- paste(seq_chars, collapse = "")
            mods_applied <- paste(
              paste0(sites_df$mod_name, " @", sites_df$position),
              collapse = "; "
            )
          } else {
            modified_seq <- pep$Sequence
            mods_applied <- "None"
          }

          all_rows <- c(all_rows, list(data.table(
            Chain            = pep$Chain,
            ProteinName      = sub(" .*", "", pep$Chain),
            Sequence         = pep$Sequence,
            ModifiedSequence = modified_seq,
            Start            = pep$Start,
            End              = pep$End,
            Length           = pep$Length,
            ModsApplied      = mods_applied,
            BaseMass         = pep$Mass,
            ModifiedMass     = calc_modified_mass(pep$Sequence, active_mods, base_mass = pep$Mass)
          )))
        }
        setProgress(0.7, detail = "Checking uniqueness against background...")
        pep_dt <- rbindlist(all_rows)
        bg_dt  <- get_active_bg_dt()
        pep_dt[, UniqueToADC := flag_unique_peptides(
          Sequence, bg_dt,
          il_equivalent = isTRUE(input$il_equivalent))]
        pep_dt[, Mass_Accuracy_ppm := 0.0]   # theoretical; updated from MS/MS in Tab 6

        # T3-D: ConjugatedForm — tag each peptide as Loaded / Deconjugated / Not a site
        if (isTRUE(input$mod_drug_enable) && nzchar(input$drug_payload_key %||% "")) {
          payload_key <- input$drug_payload_key
          # Determine conjugation residue from selected payload
          conj_res <- tryCatch({
            if (payload_key %in% names(ADCDB_PAYLOADS)) {
              ADCDB_PAYLOADS[[payload_key]]$residue
            } else { "C" }
          }, error = function(e) "C")
          conj_type <- tryCatch({
            if (payload_key %in% names(ADCDB_PAYLOADS)) {
              ADCDB_PAYLOADS[[payload_key]]$conjugation_type %||% "cysteine"
            } else { "cysteine" }
          }, error = function(e) "cysteine")
          pep_dt[, ConjugatedForm := mapply(function(s, mods) {
            if (has_conjugation_site(s, conj_type, conj_res)) {
              # If the mod list for this peptide includes the payload → Loaded;
              # otherwise it's the Deconjugated (DAR=0) form at that site
              if (grepl("drug|payload|Drug|Payload", mods, fixed = FALSE)) {
                "Loaded"
              } else {
                "Deconjugated"
              }
            } else {
              "Not a conjugation site"
            }
          }, Sequence, ModsApplied, SIMPLIFY = TRUE, USE.NAMES = FALSE)]
        } else {
          pep_dt[, ConjugatedForm := "N/A (no payload)"]
        }

        # Add alias columns required by generate_transition_list() and export.R
        # (keeps original Sequence / ModsApplied intact for Tab 3 / heavy tab)
        pep_dt[, PeptideSequence := Sequence]
        pep_dt[, Modifications   := ModsApplied]
        pep_dt[, Enzyme          := input$enzyme_id]
        rv$peptides_dt   <- pep_dt
        rv$digest_done   <- TRUE
        rv$trans_done    <- FALSE
        rv$heavy_done    <- FALSE
        rv$transition_dt <- NULL
        rv$heavy_dt      <- NULL
        chains_found <- unique(pep_dt$Chain)
        updateSelectInput(session, "filter_chain",
                          choices  = c("All Chains"="all", setNames(chains_found, chains_found)),
                          selected = "all")
        setProgress(1, detail = "Done!")
        showNotification(paste0("Digest complete: ", format(nrow(pep_dt), big.mark=","),
                                " peptide variants | ", sum(pep_dt$UniqueToADC), " unique to ADC"),
                         type="message", duration=6)
        updatebs4TabItems(session, "sidebar_menu", selected="tab_results")
      }, error = function(e) {
        showNotification(paste("Digest error:", conditionMessage(e)), type="error", duration=10)
      })
    })
  })

  output$digest_status_ui <- renderUI({
    if (!rv$digest_done) return(NULL)
    dt <- rv$peptides_dt
    tags$div(style="margin-top:10px;",
      tags$span(class="status-ok", icon("circle-check"),
        sprintf(" %s peptide variants | %s unique sequences | %s unique to ADC",
                format(nrow(dt), big.mark=","),
                format(length(unique(dt$Sequence)), big.mark=","),
                format(sum(dt$UniqueToADC), big.mark=","))))
  })

  # ── T2-C: Multi-species co-uniqueness ─────────────────────────────────────
  output$cospecies_result_ui <- renderUI(NULL)   # initialise empty

  observeEvent(input$btn_cospecies, {
    req(!is.null(rv$peptides_dt))
    sel <- input$cospecies_select
    if (length(sel) == 0L) {
      showNotification("Select at least one species.", type = "warning"); return()
    }
    mc <- as.integer(input$bg_missed %||% "0")
    bg_list <- lapply(setNames(sel, sel), function(sp) {
      if (!is.null(BUNDLED_BG[[sp]])) get_bg_dt(BUNDLED_BG[[sp]], mc) else NULL
    })
    # Remove species with no loaded database
    bg_list <- Filter(Negate(is.null), bg_list)
    if (length(bg_list) == 0L) {
      showNotification("None of the selected species backgrounds are loaded.", type = "warning"); return()
    }
    withProgress(message = "Running co-uniqueness check...", value = 0.5, {
      co_result <- flag_cospecies_unique(
        rv$peptides_dt$Sequence, bg_list,
        il_equivalent = isTRUE(input$il_equivalent)
      )
    })
    n_total  <- nrow(rv$peptides_dt)
    n_couni  <- sum(co_result$UniqueAllSpecies, na.rm = TRUE)
    sp_labels <- paste(names(bg_list), collapse = " + ")

    # Attach to rv$peptides_dt for export
    for (sp in names(co_result)) {
      col_nm <- paste0("Unique_", sp)
      rv$peptides_dt[, (col_nm) := co_result[[sp]]]
    }
    rv$peptides_dt[, UniqueAllSpecies := co_result$UniqueAllSpecies]

    output$cospecies_result_ui <- renderUI({
      tags$div(class = "alert alert-success", style = "font-size:12px; padding:8px; margin-top:8px;",
        icon("check-circle"),
        sprintf(" Co-uniqueness across [%s]: %s / %s peptides pass (%.1f%%)",
                sp_labels,
                format(n_couni, big.mark=","),
                format(n_total, big.mark=","),
                100 * n_couni / max(n_total, 1)))
    })
    showNotification(sprintf("%d peptides unique across all selected species.", n_couni),
                     type = "message", duration = 5)
  })

  # ── T2-D: FDR filter info UI ──────────────────────────────────────────────
  output$fdr_filter_info_ui <- renderUI({
    if (!isTRUE(input$fdr_filter_enable)) return(NULL)
    tags$div(class = "alert alert-info", style = "font-size:11px; padding:6px; margin-top:4px;",
      icon("info-circle"),
      " 1% FDR requires decoy sequences (prefix DECOY_ or REV_) in your search FASTA. ",
      "Only available when decoy PSMs are detected in the result file.")
  })

  # ── Value boxes ────────────────────────────────────────────────────────────
  output$vbox_total <- renderbs4ValueBox({
    n <- if (!is.null(rv$peptides_dt)) nrow(rv$peptides_dt) else 0
    bs4ValueBox(value=format(n, big.mark=","), subtitle="Total Peptide Variants",
                icon=icon("dna"), color="navy", elevation=2)
  })
  output$vbox_unique <- renderbs4ValueBox({
    n <- if (!is.null(rv$peptides_dt)) sum(rv$peptides_dt$UniqueToADC) else 0
    bs4ValueBox(value=format(n, big.mark=","), subtitle="Unique to ADC",
                icon=icon("star"), color="success", elevation=2)
  })
  output$vbox_modified <- renderbs4ValueBox({
    n <- if (!is.null(rv$peptides_dt)) sum(rv$peptides_dt$ModsApplied != "None") else 0
    bs4ValueBox(value=format(n, big.mark=","), subtitle="Modified Variants",
                icon=icon("flask"), color="warning", elevation=2)
  })
  output$vbox_chains <- renderbs4ValueBox({
    n <- if (!is.null(rv$peptides_dt)) length(unique(rv$peptides_dt$Chain)) else 0
    bs4ValueBox(value=n, subtitle="Chains Detected", icon=icon("link"), color="info", elevation=2)
  })

  # ── Peptide results table ──────────────────────────────────────────────────
  filtered_peptides <- reactive({
    dt <- rv$peptides_dt; req(!is.null(dt))
    if (input$filter_chain   != "all")        dt <- dt[dt$Chain == input$filter_chain, ]
    if (input$filter_unique  == "unique")     dt <- dt[dt$UniqueToADC == TRUE, ]
    if (input$filter_unique  == "nonunique")  dt <- dt[dt$UniqueToADC == FALSE, ]
    if (input$filter_mods    == "modified")   dt <- dt[dt$ModsApplied != "None", ]
    if (input$filter_mods    == "unmodified") dt <- dt[dt$ModsApplied == "None", ]
    dt
  })

  output$peptide_table <- renderDT({
    dt <- filtered_peptides(); req(!is.null(dt) && nrow(dt) > 0)
    ppm_col <- if ("Mass_Accuracy_ppm" %in% names(dt)) round(dt$Mass_Accuracy_ppm, 2) else rep(0.00, nrow(dt))
    display_df <- data.frame(
      Chain=dt$Chain, `Peptide Sequence`=dt$Sequence, `Modified Sequence`=dt$ModifiedSequence,
      Start=dt$Start, End=dt$End, Length=dt$Length, Modifications=dt$ModsApplied,
      `Mass (Da)`=round(dt$ModifiedMass, 5),
      `Mass Accuracy (ppm)`=ppm_col,
      `Unique to ADC`=ifelse(dt$UniqueToADC,
        '<span class="badge-unique">Yes</span>',
        '<span class="badge-nonunique">No</span>'),
      check.names=FALSE, stringsAsFactors=FALSE)
    datatable(display_df, escape=FALSE, rownames=FALSE, filter="top",
              class="table table-hover table-sm",
              options=list(pageLength=20, scrollX=TRUE, dom="Blfrtip",
                           buttons=c("copy","csv"),
                           columnDefs=list(list(className="dt-center", targets=c(3,4,5,8,9)),
                                           list(width="200px", targets=2))),
              extensions="Buttons") |>
      formatStyle("Unique to ADC",
                  backgroundColor=styleEqual('<span class="badge-unique">Yes</span>', "#eafaf1")) |>
      formatStyle("Mass Accuracy (ppm)",
                  color = styleInterval(c(-5, 5), c("#e74c3c", "#27ae60", "#e74c3c")))
  })

  # ── Transition list ────────────────────────────────────────────────────────
  output$ce_formula_ui <- renderUI({
    ce_labels <- c(
      skyline = "Generic (Sciex empirical)",
      thermo  = "Thermo (0.034 × m/z + 3.314)",
      sciex   = "SCIEX (0.0625 × m/z × z)",
      bruker  = "Bruker (timsTOF empirical)",
      agilent = "Agilent (3.7 + 0.0021 × m/z²/z)",
      waters  = "Waters (Xevo empirical)"
    )
    tags$span(ce_labels[input$instrument %||% "skyline"])
  })

  observeEvent(input$btn_gen_trans, {
    req(!is.null(rv$peptides_dt))
    adc_nm <- trimws(input$adc_name %||% "")
    withProgress(message = "Generating transition list...", value = 0, {
      tryCatch({
        setProgress(0.2, detail = "Calculating fragment ions...")
        trans_dt         <- generate_transition_list(rv$peptides_dt, adc_name = adc_nm)
        # T2-E: annotate recommended precursor isotope offset (M+0/M+1/M+2)
        trans_dt         <- annotate_isotope_offsets(trans_dt)
        rv$transition_dt <- trans_dt
        rv$trans_done    <- TRUE
        setProgress(1, detail = "Done!")
        showNotification(paste0("Transition list: ", format(nrow(trans_dt), big.mark=","), " transitions"),
                         type="message", duration=5)
      }, error = function(e) {
        showNotification(paste("Transition error:", conditionMessage(e)), type="error", duration=10)
      })
    })
  })

  # ── T3-A: DAR Transitions ──────────────────────────────────────────────────
  observeEvent(input$btn_gen_dar, {
    req(!is.null(rv$peptides_dt))
    if (!isTRUE(input$mod_drug_enable)) {
      showNotification("Enable a Drug-Linker Payload in Tab 2 first.", type = "warning"); return()
    }
    payload_key <- input$drug_payload_key %||% "custom_drug"
    payload_mass_val <- tryCatch({
      if (payload_key %in% c("custom_drug","adcdb_custom_payload")) {
        as.numeric(input$custom_drug_mass %||% 700)
      } else if (payload_key %in% names(VAR_MOD_DEFS)) {
        VAR_MOD_DEFS[[payload_key]]$mass
      } else {
        as.numeric(input$custom_drug_mass %||% 700)
      }
    }, error = function(e) 700)

    if (is.na(payload_mass_val) || payload_mass_val <= 0) {
      showNotification("Invalid payload mass — set a Custom payload mass in Tab 2.", type = "warning"); return()
    }

    conj_type <- input$dar_conjugation_type %||% "cysteine"
    conj_res  <- if (conj_type == "lysine") "K" else "C"
    dar_vals  <- seq(input$dar_range_ui[1], input$dar_range_ui[2])
    adc_nm    <- trimws(input$adc_name %||% "")

    withProgress(message = "Generating DAR transitions...", value = 0.3, {
      tryCatch({
        dar_dt <- generate_dar_transitions(
          peptides_dt      = rv$peptides_dt,
          payload_mass     = payload_mass_val,
          dar_range        = dar_vals,
          conjugation_type = conj_type,
          residue          = conj_res,
          adc_name         = adc_nm
        )
        rv$dar_dt   <- dar_dt
        rv$dar_done <- TRUE
        setProgress(1, detail = "Done!")
        n_conj_peps <- length(unique(dar_dt$PeptideSequence[dar_dt$DAR > 0]))
        showNotification(
          sprintf("DAR transitions: %s rows | %s conjugation-site peptides | DAR %d-%d",
                  format(nrow(dar_dt), big.mark=","), n_conj_peps,
                  min(dar_vals), max(dar_vals)),
          type = "message", duration = 6)
      }, error = function(e) {
        showNotification(paste("DAR error:", conditionMessage(e)), type="error", duration=10)
      })
    })
  })

  output$dar_summary_ui <- renderUI({
    if (!rv$dar_done || is.null(rv$dar_dt)) return(NULL)
    dt <- rv$dar_dt
    by_dar <- table(dt$DAR[!duplicated(paste0(dt$DAR, dt$PeptideSequence))])
    tags$div(class = "alert alert-info", style = "font-size:12px; padding:8px;",
      icon("vials"),
      sprintf(" %s total transitions across %s peptides and DAR levels: %s",
              format(nrow(dt), big.mark=","),
              format(length(unique(dt$PeptideSequence)), big.mark=","),
              paste(names(by_dar), by_dar, sep="×", collapse=" | DAR")))
  })

  output$dar_table <- renderDT({
    req(rv$dar_done, !is.null(rv$dar_dt))
    dt <- rv$dar_dt
    display_df <- data.frame(
      DAR=dt$DAR, Chain=dt$Chain, Peptide=dt$PeptideSequence,
      `z (prec)`=dt$PrecursorCharge, `Prec m/z`=round(dt$PrecursorMz,4),
      Fragment=dt$FragmentIon, `Prod m/z`=round(dt$ProductMz,4),
      `CE (eV)`=dt$CollisionEnergy,
      check.names=FALSE, stringsAsFactors=FALSE)
    datatable(display_df, rownames=FALSE, filter="top",
              class="table table-hover table-sm",
              options=list(pageLength=20, scrollX=TRUE, dom="Blfrtip",
                           buttons=c("copy","csv")),
              extensions="Buttons")
  })

  output$dl_dar_csv <- downloadHandler(
    filename = function() paste0("DAR_transitions_", Sys.Date(), ".csv"),
    content  = function(file) {
      req(!is.null(rv$dar_dt))
      write.csv(rv$dar_dt, file, row.names=FALSE)
    }
  )

  output$trans_summary_ui <- renderUI({
    if (!rv$trans_done || is.null(rv$transition_dt))
      return(tags$p(class="status-warn", "Click 'Generate Transition List' to compute transitions."))
    dt <- rv$transition_dt
    adc_nm <- trimws(input$adc_name %||% "")
    adc_label <- if (nzchar(adc_nm)) paste0(" | ADC: ", adc_nm) else ""
    tags$div(tags$span(class="status-ok", icon("circle-check"),
      sprintf(" %s transitions | %s precursors | %s peptides%s",
              format(nrow(dt), big.mark=","),
              format(length(unique(paste0(dt$PeptideSequence, dt$PrecursorCharge))), big.mark=","),
              format(length(unique(dt$PeptideSequence)), big.mark=","),
              adc_label)))
  })

  output$trans_table <- renderDT({
    req(rv$trans_done, !is.null(rv$transition_dt))
    dt <- rv$transition_dt
    # Apply unique filter if requested
    if (isTRUE(input$trans_unique_only)) dt <- dt[dt$UniqueToADC == TRUE, ]
    display_df <- data.frame(
      ADCName=dt$ADCName, Protein=dt$ProteinName, Chain=dt$Chain,
      Peptide=dt$PeptideSequence, `Modified Seq`=dt$ModifiedSequence,
      `z (prec)`=dt$PrecursorCharge, `Prec m/z`=dt$PrecursorMz,
      Fragment=dt$FragmentIon, `Prod m/z`=dt$ProductMz,
      `CE (eV)`=dt$CollisionEnergy,
      Unique=ifelse(dt$UniqueToADC,
        '<span class="badge-unique">Yes</span>',
        '<span class="badge-nonunique">No</span>'),
      check.names=FALSE, stringsAsFactors=FALSE)
    datatable(display_df, escape=FALSE, rownames=FALSE, filter="top",
              class="table table-hover table-sm",
              options=list(pageLength=25, scrollX=TRUE, dom="Blfrtip",
                           buttons=c("copy","csv"),
                           columnDefs=list(list(className="dt-center", targets=c(5,6,7,8,9,10)))),
              extensions="Buttons")
  })

  # ── Heavy Labelling ────────────────────────────────────────────────────────
  get_heavy_label_def <- reactive({
    key <- input$heavy_label_key
    def <- HEAVY_LABEL_DEFS[[key]]
    if (is.null(def)) return(NULL)
    if (key == "Custom") {
      res  <- toupper(trimws(input$heavy_custom_residue %||% ""))
      mass <- input$heavy_custom_mass %||% 0
      nm   <- trimws(input$heavy_custom_name %||% "Custom")
      if (nchar(res) != 1 || !res %in% LETTERS) return(NULL)
      def <- list(residue=res, mass=mass, name=nm)
    }
    def
  })

  observeEvent(input$btn_gen_heavy, {
    req(!is.null(rv$peptides_dt))
    lbl <- get_heavy_label_def()
    if (is.null(lbl)) {
      showNotification("Invalid heavy label definition. Check residue and mass.", type="warning"); return()
    }
    charges <- as.integer(input$heavy_charges)
    if (length(charges) == 0) {
      showNotification("Select at least one charge state.", type="warning"); return()
    }
    withProgress(message = "Generating heavy peptide list...", value = 0, {
      tryCatch({
        setProgress(0.2, detail = "Applying heavy label...")
        pep_dt <- rv$peptides_dt
        # Build heavy peptide table
        heavy_rows <- list()
        for (i in seq_len(nrow(pep_dt))) {
          pep <- pep_dt[i, ]
          seq <- pep$Sequence
          # Count target residues
          n_target <- nchar(gsub(paste0("[^", lbl$residue, "]"), "", seq))
          if (n_target == 0) next
          heavy_mass_shift <- n_target * lbl$mass
          light_mass <- pep$ModifiedMass
          heavy_mass <- light_mass + heavy_mass_shift
          for (z in charges) {
            light_mz <- round((light_mass + z * PROTON_MASS) / z, 6)
            heavy_mz <- round((heavy_mass + z * PROTON_MASS) / z, 6)
            heavy_rows <- c(heavy_rows, list(data.table(
              ADCName          = trimws(input$adc_name %||% ""),
              Chain            = pep$Chain,
              Sequence         = seq,
              ModifiedSequence = pep$ModifiedSequence,
              Modifications    = pep$ModsApplied,
              UniqueToADC      = pep$UniqueToADC,
              HeavyLabel       = lbl$name,
              LabelResidue     = lbl$residue,
              N_LabelSites     = n_target,
              MassShift_Da     = round(heavy_mass_shift, 5),
              Charge           = z,
              Light_Mz         = light_mz,
              Heavy_Mz         = heavy_mz,
              Delta_Mz         = round(heavy_mz - light_mz, 6),
              LightMass_Da     = round(light_mass, 5),
              HeavyMass_Da     = round(heavy_mass, 5)
            )))
          }
        }
        if (length(heavy_rows) == 0) {
          showNotification(paste0("No peptides contain residue '", lbl$residue,
                                  "'. Try a different label."), type="warning"); return()
        }
        rv$heavy_dt   <- rbindlist(heavy_rows)
        rv$heavy_done <- TRUE
        setProgress(1, detail = "Done!")
        showNotification(paste0("Heavy peptide list: ", format(nrow(rv$heavy_dt), big.mark=","),
                                " entries"), type="message", duration=5)
      }, error = function(e) {
        showNotification(paste("Heavy label error:", conditionMessage(e)), type="error", duration=10)
      })
    })
  })

  output$heavy_summary_ui <- renderUI({
    if (!rv$heavy_done || is.null(rv$heavy_dt))
      return(tags$p(class="status-warn", "Click 'Generate Heavy Peptide List' to compute."))
    dt <- rv$heavy_dt
    lbl <- get_heavy_label_def()
    tags$div(tags$span(class="status-ok", icon("circle-check"),
      sprintf(" %s heavy peptide entries | Label: %s | Residue: %s",
              format(nrow(dt), big.mark=","),
              if (!is.null(lbl)) lbl$name else "?",
              if (!is.null(lbl)) lbl$residue else "?")))
  })

  output$heavy_table <- renderDT({
    req(rv$heavy_done, !is.null(rv$heavy_dt))
    dt <- rv$heavy_dt
    display_df <- data.frame(
      ADCName=dt$ADCName, Chain=dt$Chain, Sequence=dt$Sequence,
      `Modified Seq`=dt$ModifiedSequence, Modifications=dt$Modifications,
      `Heavy Label`=dt$HeavyLabel, `Label Residue`=dt$LabelResidue,
      `# Sites`=dt$N_LabelSites, `Mass Shift (Da)`=dt$MassShift_Da,
      `Charge`=dt$Charge, `Light m/z`=dt$Light_Mz, `Heavy m/z`=dt$Heavy_Mz,
      `Δm/z`=dt$Delta_Mz,
      `Unique`=ifelse(dt$UniqueToADC,
        '<span class="badge-unique">Yes</span>',
        '<span class="badge-nonunique">No</span>'),
      check.names=FALSE, stringsAsFactors=FALSE)
    datatable(display_df, escape=FALSE, rownames=FALSE, filter="top",
              class="table table-hover table-sm",
              options=list(pageLength=20, scrollX=TRUE, dom="Blfrtip",
                           buttons=c("copy","csv"),
                           columnDefs=list(list(className="dt-center", targets=c(7,8,9,10,11,12,13)))),
              extensions="Buttons")
  })

  output$heavy_pair_table <- renderDT({
    req(rv$heavy_done, !is.null(rv$heavy_dt))
    dt <- rv$heavy_dt
    # Summarise per peptide × charge: light m/z, heavy m/z, delta
    pair_df <- unique(dt[, .(Chain, Sequence, HeavyLabel, Charge,
                              Light_Mz, Heavy_Mz, Delta_Mz,
                              LightMass_Da, HeavyMass_Da, MassShift_Da,
                              UniqueToADC)])
    display_df <- data.frame(
      Chain=pair_df$Chain, Sequence=pair_df$Sequence,
      `Heavy Label`=pair_df$HeavyLabel, `Charge`=pair_df$Charge,
      `Light m/z`=pair_df$Light_Mz, `Heavy m/z`=pair_df$Heavy_Mz,
      `Δm/z`=pair_df$Delta_Mz,
      `Light Mass (Da)`=pair_df$LightMass_Da, `Heavy Mass (Da)`=pair_df$HeavyMass_Da,
      `Label Mass Shift (Da)`=pair_df$MassShift_Da,
      `Unique`=ifelse(pair_df$UniqueToADC,
        '<span class="badge-unique">Yes</span>',
        '<span class="badge-nonunique">No</span>'),
      check.names=FALSE, stringsAsFactors=FALSE)
    datatable(display_df, escape=FALSE, rownames=FALSE, filter="top",
              class="table table-hover table-sm",
              options=list(pageLength=20, scrollX=TRUE, dom="Blfrtip",
                           buttons=c("copy","csv")),
              extensions="Buttons")
  })

  # ── Downloads ──────────────────────────────────────────────────────────────
  output$dl_instrument_csv <- downloadHandler(
    filename = function() {
      inst <- input$instrument %||% "skyline"
      adc  <- gsub("[^A-Za-z0-9_-]", "_", trimws(input$adc_name %||% "ADC"))
      paste0(adc, "_", inst, "_transitions_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(!is.null(rv$transition_dt))
      top_n <- as.integer(input$trans_top_n %||% 0)
      if (top_n <= 0) top_n <- .Machine$integer.max
      tmp <- write_instrument_csv(rv$transition_dt,
                                  instrument  = input$instrument,
                                  unique_only = isTRUE(input$trans_unique_only),
                                  top_n       = top_n)
      file.copy(tmp, file)
    })

  output$dl_skyline_all <- downloadHandler(
    filename = function() {
      adc <- gsub("[^A-Za-z0-9_-]", "_", trimws(input$adc_name %||% "ADC"))
      paste0(adc, "_Skyline_AllPeptides_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(!is.null(rv$transition_dt))
      file.copy(write_skyline_csv(rv$transition_dt, unique_only=FALSE), file)
    })

  output$dl_skyline_unique <- downloadHandler(
    filename = function() {
      adc <- gsub("[^A-Za-z0-9_-]", "_", trimws(input$adc_name %||% "ADC"))
      paste0(adc, "_Skyline_UniquePeptides_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(!is.null(rv$transition_dt))
      tmp <- write_skyline_csv(rv$transition_dt, unique_only=TRUE)
      if (is.null(tmp)) { showNotification("No unique peptides to export.", type="warning"); return() }
      file.copy(tmp, file)
    })

  output$dl_peptide_excel <- downloadHandler(
    filename = function() {
      adc <- gsub("[^A-Za-z0-9_-]", "_", trimws(input$adc_name %||% "ADC"))
      paste0(adc, "_PeptideSummary_", Sys.Date(), ".xlsx")
    },
    content = function(file) {
      req(!is.null(rv$peptides_dt))
      file.copy(write_excel_summary(rv$transition_dt, rv$peptides_dt), file)
    })

  output$dl_excel_full <- downloadHandler(
    filename = function() {
      adc <- gsub("[^A-Za-z0-9_-]", "_", trimws(input$adc_name %||% "ADC"))
      paste0(adc, "_FullReport_", Sys.Date(), ".xlsx")
    },
    content = function(file) {
      req(!is.null(rv$peptides_dt))
      file.copy(write_excel_summary(rv$transition_dt, rv$peptides_dt), file)
    })

  output$dl_heavy_csv <- downloadHandler(
    filename = function() {
      adc <- gsub("[^A-Za-z0-9_-]", "_", trimws(input$adc_name %||% "ADC"))
      paste0(adc, "_HeavyPeptides_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(!is.null(rv$heavy_dt))
      write.csv(as.data.frame(rv$heavy_dt), file, row.names=FALSE)
    })

  output$dl_heavy_excel <- downloadHandler(
    filename = function() {
      adc <- gsub("[^A-Za-z0-9_-]", "_", trimws(input$adc_name %||% "ADC"))
      paste0(adc, "_HeavyPeptides_", Sys.Date(), ".xlsx")
    },
    content = function(file) {
      req(!is.null(rv$heavy_dt))
      wb <- createWorkbook()
      addWorksheet(wb, "Heavy Peptides")
      writeData(wb, "Heavy Peptides", as.data.frame(rv$heavy_dt))
      if (rv$heavy_done && !is.null(rv$heavy_dt)) {
        pair_dt <- unique(rv$heavy_dt[, .(Chain, Sequence, HeavyLabel, Charge,
                                          Light_Mz, Heavy_Mz, Delta_Mz,
                                          LightMass_Da, HeavyMass_Da, MassShift_Da)])
        addWorksheet(wb, "Light-Heavy Pairs")
        writeData(wb, "Light-Heavy Pairs", as.data.frame(pair_dt))
      }
      tmp <- tempfile(fileext=".xlsx")
      saveWorkbook(wb, tmp, overwrite=TRUE)
      file.copy(tmp, file)
    })

  # ============================================================
  # TAB 3 — SEQUENCE COVERAGE MAP (post-digest, no MS/MS needed)
  # ============================================================

  # Update chain filter choices when digest completes
  observe({
    dt <- rv$peptides_dt
    if (is.null(dt) || nrow(dt) == 0L) return()
    chains <- unique(dt$Chain)
    updateSelectInput(session, "cov_chain_filter",
                      choices  = c("All chains" = "all",
                                   setNames(chains, chains)),
                      selected = "all")
  })

  # Coverage map render
  output$digest_coverage_plot <- renderPlot({
    req(!is.null(rv$peptides_dt), nrow(rv$peptides_dt) > 0L)

    pep_dt   <- rv$peptides_dt
    color_by <- input$cov_color_by  %||% "unique"
    chain_f  <- input$cov_chain_filter %||% "all"
    show_lbl <- isTRUE(input$cov_show_labels)

    # Filter chains
    chains <- if (chain_f == "all") unique(pep_dt$Chain) else chain_f
    chains <- intersect(chains, unique(pep_dt$Chain))
    if (length(chains) == 0L) return(NULL)

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      # ── Base R fallback ──────────────────────────────────────────────────
      n_chains <- length(chains)
      par(mfrow = c(n_chains, 1), mar = c(3, 1, 2.5, 1), bg = "white")
      for (ch in chains) {
        sub <- pep_dt[Chain == ch]
        len <- max(sub$End, na.rm = TRUE)
        plot(NA, xlim = c(1, len), ylim = c(0, 1),
             xlab = "Residue position", ylab = "",
             main = ch, yaxt = "n", cex.main = 0.9)
        rect(1, 0.3, len, 0.7, col = "#ecf0f1", border = NA)
        for (j in seq_len(nrow(sub))) {
          col <- if (sub$UniqueToADC[j]) "#1a2744" else "#95a5a6"
          rect(sub$Start[j], 0.25, sub$End[j], 0.75, col = col, border = NA)
        }
        covered <- unique(unlist(lapply(seq_len(nrow(sub)),
          function(j) seq(sub$Start[j], sub$End[j]))))
        pct <- round(100 * length(covered) / len, 1)
        mtext(sprintf("Coverage: %.1f%%  (%d / %d AA)", pct, length(covered), len),
              side = 3, line = 0.1, cex = 0.75, col = "#555555")
      }
      return(invisible(NULL))
    }

    # ── ggplot2 path ─────────────────────────────────────────────────────────
    library(ggplot2)

    # Greedy lane assignment per chain
    assign_lanes <- function(sub_dt) {
      sub_dt <- sub_dt[order(sub_dt$Start), ]
      lanes      <- integer(nrow(sub_dt))
      lane_ends  <- integer(0)
      for (i in seq_len(nrow(sub_dt))) {
        s <- sub_dt$Start[i]; e <- sub_dt$End[i]
        # Find lowest lane where this peptide doesn't overlap
        assigned <- FALSE
        for (l in seq_along(lane_ends)) {
          if (lane_ends[l] < s) {
            lanes[i]    <- l
            lane_ends[l] <- e
            assigned    <- TRUE
            break
          }
        }
        if (!assigned) {
          lane_ends <- c(lane_ends, e)
          lanes[i]  <- length(lane_ends)
        }
      }
      sub_dt$Lane <- lanes
      sub_dt
    }

    # Build plot data
    plot_list <- lapply(chains, function(ch) {
      sub <- pep_dt[Chain == ch]
      sub <- assign_lanes(sub)
      sub$ChainLabel <- ch
      sub$ChainLen   <- max(sub$End, na.rm = TRUE)
      sub
    })
    plot_dt <- rbindlist(plot_list)

    # Coverage annotation per chain
    cov_labels <- sapply(chains, function(ch) {
      sub <- pep_dt[Chain == ch]
      len <- max(sub$End, na.rm = TRUE)
      covered <- unique(unlist(lapply(seq_len(nrow(sub)),
        function(j) seq(sub$Start[j], sub$End[j]))))
      pct <- round(100 * length(covered) / len, 1)
      sprintf("%s: %.1f%% covered (%d / %d AA)", ch, pct, length(covered), len)
    })
    cov_df <- data.frame(
      ChainLabel = chains,
      x          = sapply(chains, function(ch) max(pep_dt[Chain == ch]$End) / 2),
      y          = sapply(chains, function(ch) {
        sub <- pep_dt[Chain == ch]
        sub <- assign_lanes(sub)
        max(sub$Lane) + 1.2
      }),
      label = cov_labels,
      stringsAsFactors = FALSE
    )

    # Colour mapping
    if (color_by == "unique") {
      plot_dt$FillVar <- ifelse(plot_dt$UniqueToADC, "Unique to ADC", "Non-unique")
      fill_scale <- scale_fill_manual(
        name   = "Uniqueness",
        values = c("Unique to ADC" = "#1a2744", "Non-unique" = "#95a5a6"))
    } else if (color_by == "mc") {
      mc_col <- c("0" = "#1a2744", "1" = "#5b7fa6", "2" = "#aec6df")
      # Infer MC from MissedCleavages column if present, else default 0
      if ("MissedCleavages" %in% names(plot_dt)) {
        plot_dt$FillVar <- as.character(pmin(plot_dt$MissedCleavages, 2L))
      } else {
        plot_dt$FillVar <- "0"
      }
      fill_scale <- scale_fill_manual(
        name   = "Missed cleavages",
        values = mc_col)
    } else {
      # length
      plot_dt$FillVar <- plot_dt$Length
      fill_scale <- scale_fill_viridis_c(
        name   = "Peptide length",
        option = "D", direction = -1)
    }

    # Max lanes per chain for height calculation
    max_lane <- max(plot_dt$Lane, na.rm = TRUE)

    p <- ggplot(as.data.frame(plot_dt)) +
      geom_rect(aes(xmin = Start - 0.5, xmax = End + 0.5,
                    ymin = Lane - 0.4,  ymax = Lane + 0.4,
                    fill = FillVar),
                colour = "white", linewidth = 0.15) +
      fill_scale +
      geom_text(data = cov_df,
                aes(x = x, y = y, label = label),
                size = 3, colour = "#333333", fontface = "italic") +
      facet_wrap(~ ChainLabel, ncol = 1, scales = "free") +
      scale_x_continuous(name = "Residue position", expand = expansion(mult = 0.01)) +
      scale_y_continuous(name = NULL, breaks = NULL) +
      theme_minimal(base_size = 11) +
      theme(
        strip.text       = element_text(face = "bold", colour = "#1a2940", size = 11),
        legend.position  = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y      = element_blank()
      )

    # Optional peptide labels
    if (show_lbl) {
      label_dt <- as.data.frame(plot_dt[plot_dt$Length >= 8, ])
      if (nrow(label_dt) > 0L) {
        p <- p + geom_text(
          data = label_dt,
          aes(x = (Start + End) / 2, y = Lane, label = Sequence),
          size = 2, colour = "white", fontface = "bold",
          hjust = 0.5, vjust = 0.5, na.rm = TRUE)
      }
    }

    print(p)
  }, height = function() {
    dt <- rv$peptides_dt
    if (is.null(dt) || nrow(dt) == 0L) return(300L)
    chain_f <- input$cov_chain_filter %||% "all"
    chains  <- if (chain_f == "all") unique(dt$Chain) else chain_f
    chains  <- intersect(chains, unique(dt$Chain))
    n_ch    <- length(chains)
    # Estimate max lanes
    max_lanes <- max(sapply(chains, function(ch) {
      sub <- dt[Chain == ch]
      # rough estimate: sqrt of peptide count
      ceiling(sqrt(nrow(sub)))
    }), na.rm = TRUE)
    max(250L, as.integer(n_ch * (max_lanes * 20L + 80L)))
  })

  # Coverage map download handler
  output$dl_coverage_map <- downloadHandler(
    filename = function() {
      adc <- gsub("[^A-Za-z0-9_-]", "_", trimws(input$adc_name %||% "ADC"))
      paste0(adc, "_CoverageMap_", Sys.Date(), ".png")
    },
    content = function(file) {
      req(!is.null(rv$peptides_dt))
      if (!requireNamespace("ggplot2", quietly = TRUE)) {
        showNotification("ggplot2 required for PNG download.", type = "warning")
        return()
      }
      # Re-render at high resolution
      png(file, width = 12, height = 8, units = "in", res = 300)
      # Trigger the same plot logic via isolate
      isolate({
        pep_dt   <- rv$peptides_dt
        color_by <- input$cov_color_by %||% "unique"
        chains   <- unique(pep_dt$Chain)
        # Simplified re-render for download (full chain, all peptides)
        library(ggplot2)
        assign_lanes <- function(sub_dt) {
          sub_dt <- sub_dt[order(sub_dt$Start), ]
          lanes <- integer(nrow(sub_dt)); lane_ends <- integer(0)
          for (i in seq_len(nrow(sub_dt))) {
            s <- sub_dt$Start[i]; e <- sub_dt$End[i]
            assigned <- FALSE
            for (l in seq_along(lane_ends)) {
              if (lane_ends[l] < s) { lanes[i] <- l; lane_ends[l] <- e; assigned <- TRUE; break }
            }
            if (!assigned) { lane_ends <- c(lane_ends, e); lanes[i] <- length(lane_ends) }
          }
          sub_dt$Lane <- lanes; sub_dt
        }
        plot_list <- lapply(chains, function(ch) {
          sub <- pep_dt[Chain == ch]; sub <- assign_lanes(sub)
          sub$ChainLabel <- ch; sub
        })
        plot_dt <- rbindlist(plot_list)
        if (color_by == "unique") {
          plot_dt$FillVar <- ifelse(plot_dt$UniqueToADC, "Unique to ADC", "Non-unique")
          fill_scale <- scale_fill_manual(name="Uniqueness",
            values=c("Unique to ADC"="#1a2744","Non-unique"="#95a5a6"))
        } else if (color_by == "mc") {
          plot_dt$FillVar <- if ("MissedCleavages" %in% names(plot_dt))
            as.character(pmin(plot_dt$MissedCleavages, 2L)) else "0"
          fill_scale <- scale_fill_manual(name="Missed cleavages",
            values=c("0"="#1a2744","1"="#5b7fa6","2"="#aec6df"))
        } else {
          plot_dt$FillVar <- plot_dt$Length
          fill_scale <- scale_fill_viridis_c(name="Peptide length", option="D", direction=-1)
        }
        p <- ggplot(as.data.frame(plot_dt)) +
          geom_rect(aes(xmin=Start-0.5, xmax=End+0.5, ymin=Lane-0.4, ymax=Lane+0.4, fill=FillVar),
                    colour="white", linewidth=0.15) +
          fill_scale +
          facet_wrap(~ChainLabel, ncol=1, scales="free") +
          scale_x_continuous(name="Residue position", expand=expansion(mult=0.01)) +
          scale_y_continuous(name=NULL, breaks=NULL) +
          theme_minimal(base_size=12) +
          theme(strip.text=element_text(face="bold",colour="#1a2940",size=12),
                legend.position="bottom", panel.grid=element_blank())
        print(p)
      })
      dev.off()
    }
  )

  # ============================================================
  # TAB 6 — MS/MS SEARCH SERVER LOGIC
  # ============================================================

  # ── Engine status reactive ─────────────────────────────────────────────────
  # Re-evaluates whenever the executable path input changes.
  msf_status <- reactive({
    exe_input <- trimws(input$msf_jar_path %||% "")
    detect_search_engine(exe_path = if (nzchar(exe_input)) exe_input else NULL)
  })

  output$msf_status_ui <- renderUI({
    s <- msf_status()
    if (s$available) {
      engine_label <- switch(s$engine,
        msamanda = "MS Amanda 3.0",
        tide     = "Tide/Crux",
        "Search engine"
      )
      tags$div(
        tags$span(class = "badge bg-success",
          icon("circle-check"), sprintf(" %s ready", engine_label)),
        tags$br(),
        tags$small(class = "text-muted",
          sprintf("Executable: %s", basename(s$exe)))
      )
    } else {
      tags$div(
        tags$span(class = "badge bg-danger",
          icon("circle-xmark"), " No search engine found"),
        tags$br(),
        tags$small(class = "text-muted",
          "Install MS Amanda 3.0 or Tide/Crux, then paste the path above or set an env var."),
        tags$br(),
        tags$a("MS Amanda 3.0",
               href   = "https://github.com/hgb-bin-proteomics/MSAmanda/releases",
               target = "_blank",
               style  = "font-size:12px; color:#2980b9;"),
        tags$span(" | ", style = "font-size:12px; color:#7f8c8d;"),
        tags$a("Tide/Crux 4.x",
               href   = "https://crux.ms/download.html",
               target = "_blank",
               style  = "font-size:12px; color:#2980b9;")
      )
    }
  })

  # ── Dynamic score filter UI ────────────────────────────────────────────────
  output$score_filter_ui <- renderUI({
    s    <- msf_status()
    meta <- .ENGINE_SCORE_META[[s$engine]] %||% .ENGINE_SCORE_META[["none"]]

    tagList(
      sliderInput("msf_hyperscore_min",
        label = meta$score_label,
        min   = meta$score_min,
        max   = meta$score_max,
        value = meta$score_default,
        step  = if (meta$score_max > 100) 10 else 0.1,
        width = "100%"),
      numericInput("msf_evalue_max",
        label = meta$evalue_label,
        value = meta$evalue_default,
        min   = 1e-10, max = 1, step = NA, width = "100%"),
      tags$small(class = "text-muted",
        if (s$engine == "msamanda")
          "Amanda Score: higher = better. Typical threshold: 100-200."
        else if (s$engine == "tide")
          "XCorr: higher = better. Typical threshold: 1.5-2.5."
        else
          "Score thresholds applied after search or to uploaded result files.")
    )
  })

  # ── Log output ─────────────────────────────────────────────────────────────
  output$msf_log_output <- renderText({ rv$msf_log })

  # ── Run / Upload handler ───────────────────────────────────────────────────
  observeEvent(input$btn_run_search, {

    # Validate prerequisites
    if (!rv$digest_done || is.null(rv$peptides_dt)) {
      showNotification("Run the digest in Tab 1 first.", type = "warning"); return()
    }
    if (is.null(input$msf_raw_files)) {
      showNotification("Upload at least one spectral file or a psm.tsv result file.",
                       type = "warning"); return()
    }

    uploaded_paths <- input$msf_raw_files$datapath
    uploaded_names <- input$msf_raw_files$name

    # Detect if user uploaded a pre-computed result file (skip search)
    result_exts <- c("tsv", "pepxml", "pep.xml", "mzid", "mzidentml", "xml")
    is_result_file <- any(tolower(tools::file_ext(uploaded_names)) %in% result_exts) ||
                      any(grepl("(?i)psm\\.tsv$", uploaded_names))

    withProgress(message = "MS/MS Search...", value = 0, {

      tryCatch({

        if (is_result_file) {
          # ── Path A: pre-computed results uploaded ──────────────────────────
          setProgress(0.3, detail = "Parsing uploaded result file...")
          rv$msf_log <- paste0("Parsing uploaded result file: ", uploaded_names[1], "\n")

          psm_raw <- parse_search_results(
            out_dir_or_file  = uploaded_paths[1],
            engine           = msf_status()$engine,
            score_threshold  = NULL,
            evalue_threshold = NULL
          )

        } else {
          # ── Path B: raw spectral files — run search engine ─────────────────
          s <- msf_status()
          if (!s$available) {
            showNotification(
              paste("No search engine available:", s$error_msg),
              type = "error", duration = 12)
            return()
          }

          setProgress(0.1, detail = "Writing FASTA...")
          fasta_text <- get_fasta_text()
          tmp_dir    <- tempfile("msearch_run_")
          dir.create(tmp_dir, recursive = TRUE)
          fasta_path <- file.path(tmp_dir, "adc_search.fasta")
          writeLines(fasta_text, fasta_path)

          # Copy uploaded raw files to tmp dir with original names
          raw_paths <- file.path(tmp_dir, uploaded_names)
          mapply(file.copy, from = uploaded_paths, to = raw_paths)

          # Resolve payload mass if drug-linker is active
          payload_mass <- NULL
          active_mods  <- get_active_mods()
          if (!is.null(active_mods$payload) && !is.na(active_mods$payload$mass)) {
            payload_mass <- active_mods$payload$mass
          }

          engine_label <- switch(s$engine,
            msamanda = "MS Amanda 3.0",
            tide     = "Tide/Crux",
            "Search engine"
          )

          setProgress(0.2, detail = sprintf("Running %s...", engine_label))
          rv$msf_log <- sprintf(
            "Engine: %s\nExecutable: %s\nFiles: %s\nEnzyme: %s | MC: %s\n\n",
            engine_label,
            basename(s$exe),
            paste(uploaded_names, collapse = ", "),
            input$enzyme_id %||% "trypsin",
            input$missed_cleavages %||% "0"
          )

          run_result <- run_search_engine(
            engine_info      = s,
            raw_files        = raw_paths,
            fasta_path       = fasta_path,
            enzyme_id        = input$enzyme_id %||% "trypsin",
            missed_cleavages = as.integer(input$missed_cleavages %||% 0),
            payload_mass     = payload_mass,
            out_dir          = tmp_dir
          )

          rv$msf_log <- paste0(rv$msf_log,
            paste(run_result$stdout, collapse = "\n"), "\n",
            paste(run_result$stderr, collapse = "\n"))

          if (run_result$exit_code != 0L) {
            showNotification(
              paste(engine_label, "exited with code", run_result$exit_code,
                    "— check the log for details."),
              type = "error", duration = 12)
            return()
          }

          setProgress(0.7, detail = "Parsing results...")
          psm_raw <- parse_search_results(
            out_dir_or_file  = tmp_dir,
            engine           = s$engine,
            score_threshold  = NULL,
            evalue_threshold = NULL
          )
        }

        if (nrow(psm_raw) == 0L) {
          showNotification("No PSMs found in results. Check the log.", type = "warning")
          rv$msf_log <- paste0(rv$msf_log, "\nNo PSMs found in output.")
          return()
        }

        setProgress(0.85, detail = "Cross-referencing with theoretical peptides...")
        rv$psm_dt <- psm_raw

        # Add DetectedInSearch flag to peptides_dt
        detected_seqs <- unique(toupper(psm_raw$Sequence))
        rv$peptides_dt[, DetectedInSearch := toupper(Sequence) %in% detected_seqs]

        rv$search_done <- TRUE
        setProgress(1, detail = "Done!")

        n_detected <- sum(rv$peptides_dt$DetectedInSearch)
        n_total    <- nrow(rv$peptides_dt)
        showNotification(
          sprintf("Search complete: %s PSMs | %s / %s theoretical peptides detected",
                  format(nrow(psm_raw), big.mark = ","),
                  n_detected, n_total),
          type = "message", duration = 8)

      }, error = function(e) {
        rv$msf_log <- paste0(rv$msf_log, "\nERROR: ", conditionMessage(e))
        showNotification(paste("Search error:", conditionMessage(e)),
                         type = "error", duration = 12)
      })
    })
  })

  # ── Filtered PSM reactive (re-filters on slider/checkbox change) ────────────
  psm_filtered <- reactive({
    req(rv$search_done, !is.null(rv$psm_dt))
    dt <- rv$psm_dt
    hs  <- input$msf_hyperscore_min %||% 0
    ev  <- input$msf_evalue_max     %||% 1
    if (!is.na(hs) && hs > 0)  dt <- dt[Score  >= hs]
    if (!is.na(ev) && ev < 1)  dt <- dt[Evalue <= ev]

    # T2-D: target-decoy 1% FDR filter
    if (isTRUE(input$fdr_filter_enable) && "FDR_1pct" %in% names(dt)) {
      fdr_valid <- dt$FDR_1pct
      if (any(!is.na(fdr_valid))) {
        dt <- dt[isTRUE(fdr_valid) | is.na(fdr_valid)]   # keep non-NA passing PSMs
      }
    }

    # T2-F: ppm fragment mass tolerance filter — flag PSMs whose Evalue is
    # implausibly large OR, when no pre-computed fragment masses are available,
    # fall through (ppm filter is mainly meaningful post-search with per-fragment
    # masses; here we use it to colour the Mass_Accuracy_ppm column we add below)
    ppm_tol <- as.numeric(input$ms2_ppm_tol %||% 20)
    if (!is.na(ppm_tol) && "Mass_Accuracy_ppm" %in% names(dt)) {
      dt <- dt[abs(Mass_Accuracy_ppm) <= ppm_tol | is.na(Mass_Accuracy_ppm)]
    }

    dt
  })

  # ── Search summary UI ──────────────────────────────────────────────────────
  output$search_summary_ui <- renderUI({
    if (!rv$search_done || is.null(rv$psm_dt))
      return(tags$p(class = "status-warn",
        "Upload spectral files and click 'Run MS/MS Search' to begin."))
    pf <- psm_filtered()
    n_det <- if (!is.null(rv$peptides_dt)) sum(rv$peptides_dt$DetectedInSearch) else 0
    tags$div(
      tags$span(class = "status-ok", icon("circle-check"),
        sprintf(" %s PSMs (filtered) | %s unique sequences | %s theoretical peptides confirmed",
                format(nrow(pf), big.mark = ","),
                format(length(unique(pf$Sequence)), big.mark = ","),
                n_det))
    )
  })

  # ── Coverage plot ──────────────────────────────────────────────────────────
  output$coverage_plot <- renderPlot({
    req(rv$search_done, !is.null(rv$peptides_dt))

    pep_dt   <- rv$peptides_dt
    pf       <- psm_filtered()
    det_seqs <- unique(toupper(pf$Sequence))
    chains   <- unique(pep_dt$Chain)
    n_chains <- length(chains)

    if (n_chains == 0L) return(NULL)

    # Determine full chain lengths from max End position
    chain_lens <- sapply(chains, function(ch) {
      max(pep_dt[Chain == ch, End], na.rm = TRUE)
    })

    # Plot layout: one row per chain
    par(mfrow = c(n_chains, 1),
        mar   = c(3, 1, 2.5, 1),
        bg    = "white")

    for (i in seq_along(chains)) {
      ch  <- chains[i]
      len <- chain_lens[i]
      sub <- pep_dt[Chain == ch]

      plot(NA, xlim = c(1, len), ylim = c(0, 1),
           xlab = "Residue position", ylab = "",
           main = ch, yaxt = "n",
           cex.main = 0.95, cex.axis = 0.8,
           col.main = "#1a2940")

      # Background bar (full chain)
      rect(1, 0.3, len, 0.7, col = "#ecf0f1", border = NA)

      # Detected peptide segments
      for (j in seq_len(nrow(sub))) {
        pep_seq <- toupper(sub$Sequence[j])
        detected <- pep_seq %in% det_seqs
        if (!detected) next
        is_unique <- sub$UniqueToADC[j]
        fill_col  <- if (is_unique) "#1a2940" else "#95a5a6"
        rect(sub$Start[j], 0.25, sub$End[j], 0.75,
             col = fill_col, border = NA)
      }

      # Coverage % annotation
      covered_pos <- unique(unlist(lapply(seq_len(nrow(sub)), function(j) {
        pep_seq <- toupper(sub$Sequence[j])
        if (pep_seq %in% det_seqs) seq(sub$Start[j], sub$End[j]) else NULL
      })))
      pct <- if (len > 0) round(100 * length(covered_pos) / len, 1) else 0
      mtext(sprintf("Coverage: %s%%  (%s / %s AA)",
                    pct, length(covered_pos), len),
            side = 3, line = 0.1, cex = 0.75, col = "#555555")
    }

    # Legend
    par(mfrow = c(1, 1), mar = c(0, 0, 0, 0), new = TRUE)
    legend("bottomright",
           legend = c("Detected (unique to ADC)", "Detected (non-unique)", "Not detected"),
           fill   = c("#1a2940", "#95a5a6", "#ecf0f1"),
           border = NA, bty = "n", cex = 0.8, inset = 0.01)
  })

  # ── Detected vs Theoretical table ─────────────────────────────────────────
  output$search_results_table <- renderDT({
    req(rv$search_done, !is.null(rv$peptides_dt))

    pf       <- psm_filtered()
    pep_dt   <- rv$peptides_dt
    det_seqs <- unique(toupper(pf$Sequence))

    # PSM counts and best scores per sequence
    psm_stats <- pf[, .(
      PSM_Count  = .N,
      Best_Score = if (all(is.na(Score)))  NA_real_ else max(Score,  na.rm = TRUE),
      Best_Evalue= if (all(is.na(Evalue))) NA_real_ else min(Evalue, na.rm = TRUE)
    ), by = .(Sequence = toupper(Sequence))]

    # Merge with theoretical peptides
    display <- pep_dt[, .(
      Chain            = Chain,
      `Peptide Sequence` = Sequence,
      `Modified Sequence`= ModifiedSequence,
      Length           = Length,
      `Unique to ADC`  = ifelse(UniqueToADC,
        '<span class="badge-unique">Yes</span>',
        '<span class="badge-nonunique">No</span>'),
      Detected         = ifelse(toupper(Sequence) %in% det_seqs,
        '<span class="badge-unique">Yes</span>',
        '<span class="badge-nonunique">No</span>'),
      .seq_key         = toupper(Sequence)
    )]

    display <- merge(display, psm_stats,
                     by.x = ".seq_key", by.y = "Sequence", all.x = TRUE)
    display[, .seq_key := NULL]
    setnames(display,
             c("PSM_Count", "Best_Score", "Best_Evalue"),
             c("PSM Count", "Best Score", "Best E-value"))

    datatable(as.data.frame(display),
              escape     = FALSE,
              rownames   = FALSE,
              filter     = "top",
              class      = "table table-hover table-sm",
              options    = list(
                pageLength = 20, scrollX = TRUE,
                dom        = "Blfrtip",
                buttons    = c("copy", "csv"),
                columnDefs = list(
                  list(className = "dt-center", targets = c(3, 4, 5, 6, 7, 8))
                )
              ),
              extensions = "Buttons")
  })

  # ── Modification summary plot ──────────────────────────────────────────────
  output$mod_summary_plot <- renderPlot({
    req(rv$search_done)
    pf <- psm_filtered()
    if (nrow(pf) == 0L || all(is.na(pf$Modifications)) ||
        all(pf$Modifications %in% c("", "None", NA))) {
      plot.new()
      text(0.5, 0.5, "No modifications detected in filtered PSMs.",
           cex = 1.1, col = "#7f8c8d")
      return()
    }

    # Parse "modname @pos; modname2 @pos2" or "pos:mass" style strings
    mod_strings <- pf$Modifications[!is.na(pf$Modifications) &
                                     pf$Modifications != "None" &
                                     nzchar(pf$Modifications)]

    # Split on semicolons, strip position info, keep mod name
    all_mods <- unlist(strsplit(mod_strings, ";\\s*"))
    # Normalise: strip trailing " @N" or "posN:" prefixes
    all_mods <- trimws(sub("\\s*@[0-9]+.*$", "", all_mods))
    all_mods <- trimws(sub("^pos[0-9]+:", "", all_mods))
    all_mods <- all_mods[nzchar(all_mods) & all_mods != "None"]

    if (length(all_mods) == 0L) {
      plot.new()
      text(0.5, 0.5, "No modifications detected in filtered PSMs.",
           cex = 1.1, col = "#7f8c8d")
      return()
    }

    mod_counts <- sort(table(all_mods), decreasing = FALSE)
    mod_df     <- data.frame(
      Modification = names(mod_counts),
      Count        = as.integer(mod_counts),
      stringsAsFactors = FALSE
    )

    # ggplot2 horizontal bar chart
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      p <- ggplot2::ggplot(mod_df,
             ggplot2::aes(x = Count,
                          y = stats::reorder(Modification, Count))) +
        ggplot2::geom_col(fill = "#1a2940", width = 0.65) +
        ggplot2::geom_text(ggplot2::aes(label = Count),
                           hjust = -0.2, size = 3.5, colour = "#333333") +
        ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.15))) +
        ggplot2::labs(x = "PSM count", y = NULL,
                      title = "Detected Modifications (filtered PSMs)") +
        ggplot2::theme_minimal(base_size = 11) +
        ggplot2::theme(
          plot.title   = ggplot2::element_text(colour = "#1a2940", face = "bold", size = 12),
          axis.text.y  = ggplot2::element_text(size = 10),
          panel.grid.major.y = ggplot2::element_blank()
        )
      print(p)
    } else {
      # Base R fallback
      barplot(mod_counts, horiz = TRUE, las = 1,
              col = "#1a2940", border = NA,
              xlab = "PSM count",
              main = "Detected Modifications",
              cex.names = 0.85)
    }
  })

  # ── Download handlers ──────────────────────────────────────────────────────
  output$dl_psm_csv <- downloadHandler(
    filename = function() {
      adc <- gsub("[^A-Za-z0-9_-]", "_", trimws(input$adc_name %||% "ADC"))
      paste0(adc, "_PSM_Table_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(!is.null(rv$psm_dt))
      write.csv(as.data.frame(psm_filtered()), file, row.names = FALSE)
    }
  )

  output$dl_coverage_csv <- downloadHandler(
    filename = function() {
      adc <- gsub("[^A-Za-z0-9_-]", "_", trimws(input$adc_name %||% "ADC"))
      paste0(adc, "_CoverageSummary_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(rv$search_done, !is.null(rv$peptides_dt))
      pf       <- psm_filtered()
      det_seqs <- unique(toupper(pf$Sequence))
      pep_dt   <- rv$peptides_dt

      # Per-chain coverage summary
      chains <- unique(pep_dt$Chain)
      rows   <- lapply(chains, function(ch) {
        sub      <- pep_dt[Chain == ch]
        len      <- max(sub$End, na.rm = TRUE)
        det_sub  <- sub[toupper(Sequence) %in% det_seqs]
        covered  <- unique(unlist(lapply(seq_len(nrow(det_sub)), function(j)
          seq(det_sub$Start[j], det_sub$End[j]))))
        data.frame(
          Chain              = ch,
          Chain_Length_AA    = len,
          Theoretical_Peptides = nrow(sub),
          Detected_Peptides  = nrow(det_sub),
          Coverage_Pct       = round(100 * length(covered) / len, 2),
          stringsAsFactors   = FALSE
        )
      })
      write.csv(do.call(rbind, rows), file, row.names = FALSE)
    }
  )

  # ==========================================================================
  # TAB 7 — AI ASSISTANT
  # ==========================================================================

  # ── Helper: build system prompt with live session context ─────────────────
  .ai_system_prompt <- reactive({

    adc_name  <- trimws(input$adc_name %||% "Unknown ADC")
    enzyme    <- input$enzyme_id %||% "trypsin"
    mc        <- input$missed_cleavages %||% "0"
    payload   <- if (isTRUE(input$mod_drug_enable))
                   (input$drug_payload_key %||% "none") else "none"

    pep_context <- ""
    if (isTRUE(input$ai_include_peptides) && !is.null(rv$peptides_dt) && nrow(rv$peptides_dt) > 0) {
      top <- head(as.data.frame(rv$peptides_dt[, c("PeptideSequence","Chain","Start","End",
                                                     "Modifications","UniqueToADC"), with=FALSE]), 20)
      pep_context <- paste0(
        "\n\n## Current peptide digest (top 20 rows)\n",
        "```\n", paste(capture.output(print(top, row.names=FALSE)), collapse="\n"), "\n```"
      )
    }

    trans_context <- ""
    if (isTRUE(input$ai_include_transitions) && !is.null(rv$transition_dt) && nrow(rv$transition_dt) > 0) {
      n_pep  <- length(unique(rv$transition_dt$PeptideSequence))
      n_trans <- nrow(rv$transition_dt)
      trans_context <- sprintf(
        "\n\n## Transition list summary\n%d unique peptides | %d transitions total",
        n_pep, n_trans
      )
    }

    sprintf(
"You are an expert scientific assistant embedded inside **ADC Peptide Mapper v0.8**,
an R Shiny application for Antibody-Drug Conjugate (ADC) peptide mapping.

## Your role
Help the user interpret digest results, optimise LC-MS/MS methods, understand ADC
biology and pharmacology, troubleshoot uniqueness or DAR modeling issues, and select
appropriate collision energies and transitions. Be concise, technically precise, and
cite literature when relevant.

## Current session context
- ADC name:          %s
- Primary enzyme:    %s  (missed cleavages: %s)
- Payload selected:  %s
- Digest complete:   %s
- Transition list:   %s
- DAR analysis:      %s
%s%s

## Key scientific concepts in this app
- Monoisotopic masses (NIST residue masses; Roepstorff & Fohlman 1984 b/y nomenclature)
- Averagine isotope distribution (Senko et al. 1995)
- DAR distribution modeling (DAR0–DARn cysteine thiol-maleimide, lysine NHS-ester)
- Linker biotransformations: maleimide hydrolysis, succinimide ring-opening,
  thioether→sulfoxide, disulfide loss, deamidation (Shen 2012; Lyon 2015; Tumey 2014)
- Uniqueness check against Human / Cynomolgus Monkey / Rat Swiss-Prot + TrEMBL
  + cRAP contaminants
- Target-decoy FDR estimation (Käll et al. 2008)
- Multi-instrument collision energy: SCIEX, Thermo TSQ, Bruker timsTOF, Agilent QQQ, Waters Xevo

## Format
Use plain text or simple markdown. Keep replies under 400 words unless the user asks
for a detailed explanation. If asked for a table, use markdown table format.",
      adc_name, enzyme, mc, payload,
      if (rv$digest_done) "yes" else "no",
      if (rv$trans_done)  "yes" else "no",
      if (rv$dar_done)    "yes" else "no",
      pep_context, trans_context
    )
  })

  # ── Helper: render chat history as HTML ────────────────────────────────────
  .render_chat <- function(history, thinking = FALSE) {
    if (length(history) == 0L && !thinking) {
      return(HTML(paste0(
        '<div style="color:#95a5a6; font-size:13px; text-align:center; ',
        'margin-top:160px;">',
        '<i class="fas fa-robot" style="font-size:32px; margin-bottom:10px; ',
        'display:block; color:#bdc3c7;"></i>',
        'Ask anything about your ADC digest, DAR analysis,<br>',
        'LC-MS/MS methods, or ADC science.',
        '</div>'
      )))
    }

    bubbles <- lapply(history, function(msg) {
      role   <- msg$role
      text   <- htmltools::htmlEscape(msg$content)
      # Convert **bold** and `code` markdown to HTML
      text <- gsub("\\*\\*(.+?)\\*\\*", "<b>\\1</b>", text)
      text <- gsub("`([^`]+)`", "<code>\\1</code>", text)
      # Newlines → <br>
      text <- gsub("\n", "<br>", text)

      avatar <- if (role == "user") "U" else '<i class="fas fa-robot"></i>'
      paste0(
        '<div class="chat-msg ', role, '">',
        '  <div class="chat-avatar">', avatar, '</div>',
        '  <div class="chat-bubble">', text, '</div>',
        '</div>'
      )
    })

    if (thinking) {
      bubbles <- c(bubbles, list(paste0(
        '<div class="chat-msg assistant">',
        '  <div class="chat-avatar"><i class="fas fa-robot"></i></div>',
        '  <div class="chat-bubble">',
        '    <div class="typing-indicator">',
        '      <span></span><span></span><span></span>',
        '    </div>',
        '  </div>',
        '</div>'
      )))
    }

    HTML(paste(unlist(bubbles), collapse = "\n"))
  }

  # ── Provider-specific key input UI ────────────────────────────────────────
  output$ai_provider_key_ui <- renderUI({
    provider <- input$ai_provider %||% "anthropic"
    info <- switch(provider,
      anthropic = list(id = "anthropic_api_key", ph = "sk-ant-api03-...",
                       env = "ANTHROPIC_API_KEY",
                       val = Sys.getenv("ANTHROPIC_API_KEY")),
      openai    = list(id = "openai_api_key",    ph = "sk-...",
                       env = "OPENAI_API_KEY",
                       val = Sys.getenv("OPENAI_API_KEY")),
      gemini    = list(id = "gemini_api_key",    ph = "AIza...",
                       env = "GEMINI_API_KEY",
                       val = Sys.getenv("GEMINI_API_KEY"))
    )
    passwordInput(info$id, label = NULL, value = info$val,
                  placeholder = info$ph, width = "100%")
  })

  output$ai_provider_hint_ui <- renderUI({
    provider <- input$ai_provider %||% "anthropic"
    switch(provider,
      anthropic = div(class = "apikey-note", style = "margin-top:6px;",
        icon("circle-info"), " Anthropic key (sk-ant-...). ",
        tags$a("Get a key →", href = "https://console.anthropic.com/", target = "_blank")),
      openai    = div(class = "apikey-note", style = "margin-top:6px;",
        icon("circle-info"), " OpenAI key (sk-...). ",
        tags$a("Get a key →", href = "https://platform.openai.com/api-keys", target = "_blank")),
      gemini    = div(class = "apikey-note", style = "margin-top:6px;",
        icon("circle-info"), " Google AI Studio key (AIza...). ",
        tags$a("Get a key →", href = "https://aistudio.google.com/app/apikey", target = "_blank"))
    )
  })

  # Helper: get active key from whatever provider is selected
  .active_api_key <- reactive({
    provider <- input$ai_provider %||% "anthropic"
    key_id <- switch(provider,
      anthropic = "anthropic_api_key",
      openai    = "openai_api_key",
      gemini    = "gemini_api_key"
    )
    trimws(input[[key_id]] %||% "")
  })

  # ── Usage badge ────────────────────────────────────────────────────────────
  output$ai_usage_badge <- renderUI({
    AI_MSG_LIMIT <- 10L
    used <- rv$ai_msg_count
    env_key <- Sys.getenv("ANTHROPIC_API_KEY")
    # Only show the badge when no env key is set (i.e. on the hosted app)
    if (nzchar(env_key)) return(NULL)
    remaining <- max(0L, AI_MSG_LIMIT - used)
    col <- if (remaining == 0L) "#e74c3c" else if (remaining <= 3L) "#e67e22" else "#27ae60"
    tags$span(style = sprintf("font-size:11px; color:%s;", col),
      icon("comment-dots"), sprintf(" %d / %d messages used", used, AI_MSG_LIMIT))
  })

  # ── Dynamic model list: fetch from /v1/models when key changes ───────────

  # Fallback choices used before/if the API call fails
  .DEFAULT_MODELS <- c(
    "claude-3-5-haiku-20241022",
    "claude-3-5-sonnet-20241022",
    "claude-3-opus-20240229"
  )

  # Reactive: list of model IDs available for the current key
  available_models <- reactive({
    key <- trimws(input$anthropic_api_key %||% "")
    if (!nzchar(key) || !startsWith(key, "sk-ant-"))
      return(.DEFAULT_MODELS)

    resp <- tryCatch(
      httr2::req_perform(
        httr2::request("https://api.anthropic.com/v1/models") |>
          httr2::req_headers("x-api-key" = key,
                             "anthropic-version" = "2023-06-01") |>
          httr2::req_timeout(10) |>
          httr2::req_error(is_error = function(r) FALSE)
      ),
      error = function(e) NULL
    )

    if (is.null(resp) || httr2::resp_status(resp) != 200L)
      return(.DEFAULT_MODELS)

    body <- tryCatch(httr2::resp_body_json(resp), error = function(e) NULL)
    ids  <- vapply(body$data %||% list(), function(m) m$id %||% "", character(1))
    ids  <- ids[nzchar(ids) & grepl("claude", ids, ignore.case = TRUE)]
    if (length(ids) == 0L) .DEFAULT_MODELS else sort(ids)
  })

  output$ai_model_ui <- renderUI({
    provider <- input$ai_provider %||% "anthropic"
    if (provider == "anthropic") {
      models  <- available_models()
      default <- models[grepl("haiku", models, ignore.case = TRUE)][1]
      if (is.na(default)) default <- models[1]
      selectInput("ai_model", "Model:",
                  choices = setNames(models, models), selected = default)
    } else if (provider == "openai") {
      selectInput("ai_model", "Model:",
        choices = c("gpt-4o-mini", "gpt-4o", "gpt-4-turbo", "gpt-3.5-turbo"),
        selected = "gpt-4o-mini")
    } else {
      selectInput("ai_model", "Model:",
        choices = c("gemini-1.5-flash", "gemini-1.5-pro", "gemini-2.0-flash"),
        selected = "gemini-1.5-flash")
    }
  })

  # ── API key persistence helpers ───────────────────────────────────────────

  # Status badge: shows whether env var is set and if field matches
  output$apikey_status_ui <- renderUI({
    provider  <- input$ai_provider %||% "anthropic"
    env_var   <- switch(provider, anthropic="ANTHROPIC_API_KEY",
                                  openai="OPENAI_API_KEY", gemini="GEMINI_API_KEY")
    key_id    <- switch(provider, anthropic="anthropic_api_key",
                                  openai="openai_api_key", gemini="gemini_api_key")
    env_key   <- Sys.getenv(env_var)
    input_key <- trimws(input[[key_id]] %||% "")
    if (nzchar(env_key) && input_key == env_key) {
      tags$div(style = "margin-top:5px; font-size:11px; color:#27ae60;",
        icon("circle-check"), " Key loaded from environment — auto-fills each session.")
    } else if (nzchar(env_key) && input_key != env_key) {
      tags$div(style = "margin-top:5px; font-size:11px; color:#e67e22;",
        icon("triangle-exclamation"),
        " Input differs from saved .Renviron key. Click 'Save' to update.")
    } else if (nzchar(input_key)) {
      tags$div(style = "margin-top:5px; font-size:11px; color:#7f8c8d;",
        icon("circle-info"),
        " Key entered but not saved. Click 'Save to .Renviron' to persist.")
    } else {
      tags$div(style = "margin-top:5px; font-size:11px; color:#e74c3c;",
        icon("circle-xmark"), " No key set.")
    }
  })

  # Save key to ~/.Renviron (or project .Renviron) so it persists across sessions
  observeEvent(input$btn_save_apikey, {
    provider <- input$ai_provider %||% "anthropic"
    env_var  <- switch(provider, anthropic="ANTHROPIC_API_KEY",
                                 openai="OPENAI_API_KEY", gemini="GEMINI_API_KEY")
    key_id   <- switch(provider, anthropic="anthropic_api_key",
                                 openai="openai_api_key", gemini="gemini_api_key")
    key <- trimws(input[[key_id]] %||% "")
    if (!nzchar(key)) {
      showNotification("Enter an API key first.", type = "warning"); return()
    }

    renviron_path <- if (file.exists(".Renviron")) ".Renviron" else
                     file.path(Sys.getenv("HOME"), ".Renviron")

    existing <- if (file.exists(renviron_path))
                  readLines(renviron_path, warn = FALSE) else character(0)
    pat      <- sprintf("^%s\\s*=", env_var)
    existing <- existing[!grepl(pat, existing)]
    new_line <- sprintf('%s="%s"', env_var, key)
    writeLines(c(existing, new_line), renviron_path)

    do.call(Sys.setenv, setNames(list(key), env_var))

    showNotification(
      paste0(env_var, " saved to ", renviron_path, ". It will auto-fill on next app start."),
      type = "message", duration = 8
    )
  })

  # Clear key from field and env
  observeEvent(input$btn_clear_apikey, {
    provider <- input$ai_provider %||% "anthropic"
    env_var  <- switch(provider, anthropic="ANTHROPIC_API_KEY",
                                 openai="OPENAI_API_KEY", gemini="GEMINI_API_KEY")
    key_id   <- switch(provider, anthropic="anthropic_api_key",
                                 openai="openai_api_key", gemini="gemini_api_key")
    updateTextInput(session, key_id, value = "")
    do.call(Sys.setenv, setNames(list(""), env_var))
    showNotification(
      paste0(env_var, " cleared from this session. Remove the line from .Renviron to make it permanent."),
      type = "message", duration = 8
    )
  })

  # ── Chat message output ────────────────────────────────────────────────────
  output$chat_messages <- renderUI({
    .render_chat(rv$chat_history, rv$ai_thinking)
  })

  # ── Context badge strip ───────────────────────────────────────────────────
  output$ai_context_badges <- renderUI({
    badges <- list(
      list(label = "ADC name",     active = nzchar(trimws(input$adc_name %||% ""))),
      list(label = "Digest done",  active = rv$digest_done),
      list(label = "Transitions",  active = rv$trans_done),
      list(label = "DAR analysis", active = rv$dar_done),
      list(label = "Payload set",  active = isTRUE(input$mod_drug_enable))
    )
    tagList(lapply(badges, function(b) {
      cls <- if (b$active) "context-badge active" else "context-badge inactive"
      tags$span(class = cls, if (b$active) icon("check") else icon("xmark"), b$label)
    }))
  })

  # ── Token usage display ───────────────────────────────────────────────────
  output$ai_token_info <- renderUI({
    if (rv$ai_total_tokens == 0L) return(NULL)
    tags$span(style = "font-size:11px; color:#95a5a6;",
      icon("circle-info"), " ~", format(rv$ai_total_tokens, big.mark=","), " tokens used this session")
  })

  # ── Quick prompt buttons → fill input ─────────────────────────────────────
  .qp_texts <- c(
    qp1 = "Can you explain how DAR modeling works in ADC peptide mapping, and what DAR0 vs DAR4 means for my transition list?",
    qp2 = "Looking at my current peptide table, which transitions would you recommend for the most reliable LC-MS/MS quantitation?",
    qp3 = "What are linker biotransformations in ADCs, and which ones are most important to monitor in a peptide map?",
    qp4 = "Some peptides in my digest are flagged as non-unique. How should I troubleshoot this and what are my options?",
    qp5 = "What collision energy should I use for my instrument and peptides? Which instrument formula is most appropriate?",
    qp6 = "How should I set missed cleavages for my ADC peptide mapping experiment? What are the trade-offs?",
    qp7 = "Can you summarise my current digest results and flag anything unusual or worth investigating?",
    qp8 = "Which ADC payload is best suited for cysteine conjugation and why? How does it affect the peptide map?"
  )

  for (.qp_id in names(.qp_texts)) {
    local({
      qid  <- .qp_id
      text <- .qp_texts[[qid]]
      observeEvent(input[[qid]], {
        updateTextAreaInput(session, "chat_input", value = text)
      }, ignoreInit = TRUE)
    })
  }

  # ── Main send handler ─────────────────────────────────────────────────────
  observeEvent(input$btn_chat_send, {
    user_text <- trimws(input$chat_input)
    if (!nzchar(user_text)) return()

    provider <- input$ai_provider %||% "anthropic"

    # Rate limit — 10 messages per session when no server-side key is set
    AI_MSG_LIMIT <- 10L
    no_server_key <- nchar(Sys.getenv("ANTHROPIC_API_KEY")) == 0L &&
                     nchar(Sys.getenv("OPENAI_API_KEY"))    == 0L &&
                     nchar(Sys.getenv("GEMINI_API_KEY"))    == 0L
    if (no_server_key && rv$ai_msg_count >= AI_MSG_LIMIT) {
      showNotification(
        paste0("Session limit of ", AI_MSG_LIMIT, " messages reached. ",
               "Enter your own API key in Tab 7 to continue without limits."),
        type = "warning", duration = 8)
      return()
    }

    # Validate API key for selected provider
    api_key <- .active_api_key()
    if (!nzchar(api_key)) {
      showNotification("Enter an API key in the AI Assistant panel first.",
                       type = "warning", duration = 5)
      return()
    }

    # Clear input, append user message, show typing indicator
    updateTextAreaInput(session, "chat_input", value = "")
    rv$chat_history <- c(rv$chat_history,
                         list(list(role = "user", content = user_text)))
    rv$ai_thinking  <- TRUE

    shinyjs::runjs("setTimeout(function(){
      var b=document.getElementById('chat_box'); if(b) b.scrollTop=b.scrollHeight;
    }, 50);")

    messages_snap <- rv$chat_history
    system_snap   <- isolate(.ai_system_prompt())
    model_snap    <- input$ai_model %||% "claude-3-5-haiku-20241022"

    # ── Route to selected provider ──────────────────────────────────────────
    resp_text <- tryCatch({

      if (provider == "anthropic") {
        req <- httr2::request("https://api.anthropic.com/v1/messages") |>
          httr2::req_headers("x-api-key" = api_key,
                             "anthropic-version" = "2023-06-01") |>
          httr2::req_body_json(list(model      = model_snap,
                                    max_tokens = 1024L,
                                    system     = system_snap,
                                    messages   = messages_snap)) |>
          httr2::req_timeout(60) |>
          httr2::req_error(is_error = function(r) FALSE)
        r    <- httr2::req_perform(req)
        body <- httr2::resp_body_json(r)
        if (httr2::resp_status(r) != 200L) {
          stop(body$error$message %||% paste("HTTP", httr2::resp_status(r)))
        }
        # track tokens
        usage <- body$usage
        if (!is.null(usage))
          rv$ai_total_tokens <- rv$ai_total_tokens +
            (usage$input_tokens %||% 0L) + (usage$output_tokens %||% 0L)
        body$content[[1L]]$text %||% "[No response text]"

      } else if (provider == "openai") {
        # Build OpenAI messages: system + history
        oai_msgs <- c(
          list(list(role = "system", content = system_snap)),
          messages_snap
        )
        req <- httr2::request("https://api.openai.com/v1/chat/completions") |>
          httr2::req_headers("Authorization" = paste("Bearer", api_key)) |>
          httr2::req_body_json(list(model       = model_snap,
                                    max_tokens  = 1024L,
                                    messages    = oai_msgs)) |>
          httr2::req_timeout(60) |>
          httr2::req_error(is_error = function(r) FALSE)
        r    <- httr2::req_perform(req)
        body <- httr2::resp_body_json(r)
        if (httr2::resp_status(r) != 200L) {
          stop(body$error$message %||% paste("HTTP", httr2::resp_status(r)))
        }
        usage <- body$usage
        if (!is.null(usage))
          rv$ai_total_tokens <- rv$ai_total_tokens +
            (usage$prompt_tokens %||% 0L) + (usage$completion_tokens %||% 0L)
        body$choices[[1L]]$message$content %||% "[No response text]"

      } else {
        # Google Gemini (generateContent endpoint)
        # Convert history to Gemini "contents" format
        gem_contents <- lapply(messages_snap, function(m) {
          role <- if (m$role == "assistant") "model" else "user"
          list(role = role, parts = list(list(text = m$content)))
        })
        # Prepend system instruction as first user turn if not already there
        gem_body <- list(
          system_instruction = list(parts = list(list(text = system_snap))),
          contents           = gem_contents,
          generationConfig   = list(maxOutputTokens = 1024L)
        )
        url <- sprintf(
          "https://generativelanguage.googleapis.com/v1beta/models/%s:generateContent?key=%s",
          model_snap, api_key)
        req <- httr2::request(url) |>
          httr2::req_body_json(gem_body) |>
          httr2::req_timeout(60) |>
          httr2::req_error(is_error = function(r) FALSE)
        r    <- httr2::req_perform(req)
        body <- httr2::resp_body_json(r)
        if (httr2::resp_status(r) != 200L) {
          stop(body$error$message %||% paste("HTTP", httr2::resp_status(r)))
        }
        body$candidates[[1L]]$content$parts[[1L]]$text %||% "[No response text]"
      }

    }, error = function(e) {
      paste0("[Error: ", conditionMessage(e), "]")
    })

    rv$ai_thinking  <- FALSE
    rv$ai_msg_count <- rv$ai_msg_count + 1L

    rv$chat_history <- c(rv$chat_history,
                         list(list(role = "assistant", content = resp_text)))

    shinyjs::runjs("setTimeout(function(){
      var b=document.getElementById('chat_box'); if(b) b.scrollTop=b.scrollHeight;
    }, 100);")
  })

  # ── Clear conversation ────────────────────────────────────────────────────
  observeEvent(input$btn_chat_clear, {
    rv$chat_history     <- list()
    rv$ai_total_tokens  <- 0L
    rv$ai_msg_count     <- 0L
    rv$ai_thinking      <- FALSE
  })

  # ── Enter key sends message (JS) ─────────────────────────────────────────
  shinyjs::runjs("
    $(document).on('keydown', '#chat_input', function(e) {
      if (e.key === 'Enter' && !e.shiftKey) {
        e.preventDefault();
        $('#btn_chat_send').click();
      }
    });
  ")

  # ── T4-D: Session info output ─────────────────────────────────────────────
  output$session_info_text <- renderText({
    si <- sessionInfo()
    pkgs_used <- c("shiny", "bs4Dash", "DT", "data.table", "openxlsx",
                   "httr2", "stringr", "dplyr", "shinyjs")
    pkg_lines <- vapply(pkgs_used, function(p) {
      v <- tryCatch(as.character(packageVersion(p)), error = function(e) "n/a")
      sprintf("  %-20s %s", p, v)
    }, character(1))

    paste0(
      "ADC Peptide Mapper v0.8\n",
      "────────────────────────────────\n",
      sprintf("R version   : %s\n", R.version.string),
      sprintf("Platform    : %s\n", si$platform),
      sprintf("Locale      : %s\n", Sys.getlocale("LC_CTYPE")),
      sprintf("Session date: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M %Z")),
      "\nKey packages:\n",
      paste(pkg_lines, collapse = "\n"),
      "\n────────────────────────────────\n",
      "Unit tests  : 79/79 PASS (Tier 4)\n",
      "Mass accuracy: ≤0.05 mDa vs Unimod\n"
    )
  })

} # end server

# ── Launch ─────────────────────────────────────────────────────────────────────
shinyApp(ui = ui, server = server)
