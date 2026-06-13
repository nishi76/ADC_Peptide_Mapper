# =============================================================================
# msearch.R — ADC Peptide Mapper v0.8
# =============================================================================
# MS/MS search integration via MS Amanda 3.0 (primary) and Tide/Crux (fallback).
# Sourced by app.R. No library() calls — depends on data.table loaded by app.
#
# Public API
# ----------
#   detect_search_engine(exe_path = NULL)
#   run_search_engine(engine_info, raw_files, fasta_path, enzyme_id,
#                     missed_cleavages, payload_mass, out_dir)
#   parse_search_results(out_dir_or_file, engine, score_threshold, evalue_threshold)
#
# Engine resolution order:
#   1. MS Amanda 3.0  — checked first (primary)
#   2. Tide/Crux 4.x  — checked only if MS Amanda not found (fallback)
#
# Supported input formats for Path A upload (parse_search_results):
#   1. mzIdentML (.mzid / .mzidentml) — MS Amanda output
#   2. pepXML    (.pep.xml / .pepxml)  — Tide/Crux output
#   3. psm.tsv                         — FragPipe / MSFragger legacy
#   4. MS Amanda CSV (.csv)            — Amanda summary CSV
#
# Download instructions:
#   MS Amanda 3.0 : https://github.com/hgb-bin-proteomics/MSAmanda/releases
#   Tide/Crux 4.x : https://crux.ms/download.html
# =============================================================================


# =============================================================================
# ENGINE SCORE METADATA
# Used by app.R to render dynamic score filter UI in Tab 6.
# =============================================================================
.ENGINE_SCORE_META <- list(
  msamanda = list(
    score_label   = "Minimum Amanda Score",
    score_default = 100,
    score_min     = 0,
    score_max     = 2000,
    evalue_label  = "Maximum E-value",
    evalue_default = 0.01
  ),
  tide = list(
    score_label   = "Minimum XCorr",
    score_default = 1.5,
    score_min     = 0,
    score_max     = 10,
    evalue_label  = "Maximum p-value",
    evalue_default = 0.01
  ),
  none = list(
    score_label   = "Minimum Score",
    score_default = 0,
    score_min     = 0,
    score_max     = 100,
    evalue_label  = "Maximum E-value",
    evalue_default = 0.01
  )
)


# =============================================================================
# ENZYME MAPS
# =============================================================================

# MS Amanda enzyme names
.MSAMANDA_ENZYME_MAP <- c(
  trypsin       = "Trypsin",
  lysc          = "LysC",
  lysn          = "LysN",
  aspn          = "AspN",
  gluc_e        = "GluC",
  gluc_ed       = "GluC",
  chymotrypsin  = "Chymotrypsin",
  argc          = "ArgC",
  cnbr          = "CNBr",
  trypsin_lysc  = "Trypsin"
)

# Tide/Crux enzyme names
.TIDE_ENZYME_MAP <- c(
  trypsin       = "trypsin",
  lysc          = "lysC",
  lysn          = "no-enzyme",   # Tide has no LysN; log a note
  aspn          = "asp-n",
  gluc_e        = "glu-c",
  gluc_ed       = "glu-c",
  chymotrypsin  = "chymotrypsin",
  argc          = "arg-c",
  cnbr          = "no-enzyme",   # Tide has no CNBr; log a note
  trypsin_lysc  = "trypsin"
)

# Enzymes that fall back to no-enzyme in Tide (logged as a note)
.TIDE_NO_ENZYME_FALLBACK <- c("lysn", "cnbr")


# =============================================================================
# detect_search_engine()
# =============================================================================
# Locate MS Amanda or Tide/Crux and verify the binary is executable.
#
# Resolution order:
#   MS Amanda (checked first):
#     1. exe_path argument
#     2. MSAMANDA_EXE environment variable
#     3. Auto-scan: MSAmanda.exe / MSAmanda in getwd(), ~/tools/, ~/bin/,
#        ~/MSAmanda/, C:/MSAmanda/ (Windows)
#     4. system PATH via Sys.which("MSAmanda")
#
#   Tide/Crux (only if MS Amanda not found):
#     1. CRUX_EXE environment variable
#     2. Auto-scan: crux / crux.exe in getwd(), ~/tools/, ~/bin/
#     3. system PATH via Sys.which("crux")
#
# Args:
#   exe_path : character(1) or NULL — explicit path supplied by user in UI
#
# Returns:
#   list(
#     engine    = "msamanda" | "tide" | "none",
#     available = logical(1),
#     exe       = character(1),   # resolved path or ""
#     version   = character(1),   # version string or ""
#     error_msg = character(1)    # human-readable reason if available=FALSE
#   )
# =============================================================================
detect_search_engine <- function(exe_path = NULL) {

  result <- list(engine = "none", available = FALSE,
                 exe = "", version = "", error_msg = "")

  # ── Helper: run binary with --version and capture output ──────────────────
  .try_version <- function(exe) {
    tryCatch(
      system2(exe, args = "--version", stdout = TRUE, stderr = TRUE,
              timeout = 10),
      error   = function(e) NULL,
      warning = function(w) NULL
    )
  }

  # ── Helper: scan a list of directories for a binary name ──────────────────
  .scan_dirs <- function(dirs, names) {
    for (d in dirs) {
      if (!dir.exists(d)) next
      for (nm in names) {
        p <- file.path(d, nm)
        if (file.exists(p)) return(p)
      }
    }
    ""
  }

  # =========================================================================
  # 1. Try MS Amanda
  # =========================================================================
  amanda_exe <- ""

  # Priority 1: explicit path from UI
  if (!is.null(exe_path) && nzchar(trimws(exe_path))) {
    amanda_exe <- trimws(exe_path)
  }

  # Priority 2: environment variable
  if (!nzchar(amanda_exe)) {
    env_val <- Sys.getenv("MSAMANDA_EXE", unset = "")
    if (nzchar(env_val)) amanda_exe <- env_val
  }

  # Priority 3: auto-scan common locations
  if (!nzchar(amanda_exe)) {
    is_win   <- .Platform$OS.type == "windows"
    bin_names <- if (is_win) c("MSAmanda.exe", "MSAmanda") else c("MSAmanda", "MSAmanda.exe")
    scan_dirs <- c(
      getwd(),
      path.expand("~/tools"),
      path.expand("~/bin"),
      path.expand("~/MSAmanda"),
      if (is_win) "C:/MSAmanda" else character(0)
    )
    amanda_exe <- .scan_dirs(scan_dirs, bin_names)
  }

  # Priority 4: system PATH
  if (!nzchar(amanda_exe)) {
    found <- Sys.which("MSAmanda")
    if (nzchar(found)) amanda_exe <- found
  }

  # Validate MS Amanda
  if (nzchar(amanda_exe) && file.exists(amanda_exe)) {
    ver_out <- .try_version(amanda_exe)
    if (!is.null(ver_out) && length(ver_out) > 0L) {
      ver_str <- paste(ver_out, collapse = " ")
      result$engine    <- "msamanda"
      result$available <- TRUE
      result$exe       <- amanda_exe
      result$version   <- ver_str
      return(result)
    }
    # Binary exists but didn't respond — still accept it (some builds exit non-zero on --version)
    result$engine    <- "msamanda"
    result$available <- TRUE
    result$exe       <- amanda_exe
    result$version   <- "unknown"
    return(result)
  }

  # =========================================================================
  # 2. Try Tide/Crux (fallback)
  # =========================================================================
  crux_exe <- ""

  # Priority 1: environment variable
  env_crux <- Sys.getenv("CRUX_EXE", unset = "")
  if (nzchar(env_crux)) crux_exe <- env_crux

  # Priority 2: auto-scan
  if (!nzchar(crux_exe)) {
    is_win    <- .Platform$OS.type == "windows"
    bin_names <- if (is_win) c("crux.exe", "crux") else c("crux", "crux.exe")
    scan_dirs <- c(
      getwd(),
      path.expand("~/tools"),
      path.expand("~/bin"),
      path.expand("~/crux/bin")
    )
    crux_exe <- .scan_dirs(scan_dirs, bin_names)
  }

  # Priority 3: system PATH
  if (!nzchar(crux_exe)) {
    found <- Sys.which("crux")
    if (nzchar(found)) crux_exe <- found
  }

  # Validate Crux
  if (nzchar(crux_exe) && file.exists(crux_exe)) {
    ver_out <- .try_version(crux_exe)
    if (!is.null(ver_out)) {
      ver_str <- paste(ver_out, collapse = " ")
      if (grepl("crux", ver_str, ignore.case = TRUE)) {
        result$engine    <- "tide"
        result$available <- TRUE
        result$exe       <- crux_exe
        result$version   <- ver_str
        return(result)
      }
    }
  }

  # =========================================================================
  # 3. Neither found
  # =========================================================================
  result$error_msg <- paste0(
    "No MS/MS search engine found. ",
    "Install MS Amanda 3.0 (https://github.com/hgb-bin-proteomics/MSAmanda/releases) ",
    "or Tide/Crux 4.x (https://crux.ms/download.html), then either: ",
    "(a) paste the executable path in the field above, or ",
    "(b) set MSAMANDA_EXE / CRUX_EXE in ~/.Renviron."
  )
  result
}


# =============================================================================
# .write_msamanda_settings()
# =============================================================================
# Write a MS Amanda 3.0 settings.xml to out_dir.
#
# Args:
#   out_dir          : character(1) — directory to write settings.xml
#   enzyme_id        : character(1) — app enzyme_id (e.g. "trypsin")
#   missed_cleavages : integer(1)   — 0, 1, or 2
#   payload_mass     : numeric(1) or NULL — drug-linker variable mod mass (Da)
#
# Returns: character(1) — path to written settings.xml
# =============================================================================
.write_msamanda_settings <- function(out_dir, enzyme_id,
                                     missed_cleavages, payload_mass = NULL) {

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  enzyme_name <- .MSAMANDA_ENZYME_MAP[enzyme_id]
  if (is.na(enzyme_name)) enzyme_name <- "Trypsin"

  mc <- as.integer(missed_cleavages)

  # Variable modifications block
  var_mods_xml <- paste0(
    '    <Modification Type="variable" Delta="15.994915" UnimodID="21">\n',
    '      <Name>Oxidation</Name>\n',
    '      <Residues>M</Residues>\n',
    '    </Modification>'
  )

  if (!is.null(payload_mass) && is.numeric(payload_mass) && !is.na(payload_mass)) {
    var_mods_xml <- paste0(var_mods_xml, "\n",
      sprintf(
        '    <Modification Type="variable" Delta="%.6f" UnimodID="0">\n',
        payload_mass
      ),
      '      <Name>PayloadLinker</Name>\n',
      '      <Residues>C</Residues>\n',
      '    </Modification>'
    )
  }

  xml_content <- sprintf(
'<?xml version="1.0" encoding="utf-8"?>
<!-- MS Amanda settings.xml — generated by ADC Peptide Mapper v0.8 -->
<Settings>
  <Enzyme>
    <Name>%s</Name>
    <MissedCleavages>%d</MissedCleavages>
  </Enzyme>
  <Modifications>
    <Modification Type="fixed" Delta="57.021464" UnimodID="4">
      <Name>Carbamidomethyl</Name>
      <Residues>C</Residues>
    </Modification>
%s
  </Modifications>
  <MassTolerance>
    <PrecursorMassTolerance Unit="ppm">20</PrecursorMassTolerance>
    <FragmentMassTolerance Unit="Da">0.02</FragmentMassTolerance>
  </MassTolerance>
  <Instrument>HighRes</Instrument>
  <MaxRank>1</MaxRank>
  <GenerateDecoy>true</GenerateDecoy>
  <DecoyPrefix>DECOY_</DecoyPrefix>
  <PeptideLengthRange Min="6" Max="30"/>
  <PeptideMassRange Min="500" Max="5000"/>
</Settings>',
    enzyme_name, mc, var_mods_xml
  )

  settings_path <- file.path(out_dir, "settings.xml")
  writeLines(xml_content, settings_path)
  settings_path
}


# =============================================================================
# run_msamanda()
# =============================================================================
# Run MS Amanda 3.0 on one or more spectral files.
# MS Amanda 3.0 processes one file per CLI call; this function loops.
#
# Args:
#   exe              : character(1) — path to MSAmanda binary
#   raw_files        : character vector — paths to mzML or MGF files
#   fasta_path       : character(1) — path to ADC FASTA
#   enzyme_id        : character(1) — app enzyme_id
#   missed_cleavages : integer(1)
#   payload_mass     : numeric(1) or NULL
#   out_dir          : character(1) — output directory
#
# Returns:
#   list(exit_code, stdout, stderr, out_dir)
# =============================================================================
run_msamanda <- function(exe, raw_files, fasta_path, enzyme_id,
                         missed_cleavages, payload_mass = NULL, out_dir) {

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  settings_path <- .write_msamanda_settings(
    out_dir, enzyme_id, missed_cleavages, payload_mass
  )

  all_stdout <- character(0)
  all_stderr <- character(0)
  last_exit  <- 0L

  for (f in raw_files) {
    out <- tryCatch(
      system2(
        command = exe,
        args    = c(
          "-s", shQuote(normalizePath(f, mustWork = FALSE)),
          "-d", shQuote(normalizePath(fasta_path, mustWork = FALSE)),
          "-e", shQuote(normalizePath(settings_path, mustWork = FALSE)),
          "-f", "2",
          "-o", shQuote(normalizePath(out_dir, mustWork = FALSE))
        ),
        stdout  = TRUE,
        stderr  = TRUE,
        wait    = TRUE
      ),
      error = function(e) {
        all_stderr <<- c(all_stderr, conditionMessage(e))
        structure(character(0), status = 1L)
      }
    )

    if (is.character(out)) {
      all_stdout <- c(all_stdout, out)
      ec <- attr(out, "status")
      if (!is.null(ec)) last_exit <- ec
    }
  }

  list(exit_code = last_exit, stdout = all_stdout,
       stderr = all_stderr, out_dir = out_dir)
}


# =============================================================================
# run_tide()
# =============================================================================
# Run Tide search via the Crux toolkit (two-step: tide-index + tide-search).
# Both steps are transparent to the user.
#
# Args:
#   exe              : character(1) — path to crux binary
#   raw_files        : character vector — paths to mzML, mzXML, or MGF files
#   fasta_path       : character(1) — path to ADC FASTA
#   enzyme_id        : character(1) — app enzyme_id
#   missed_cleavages : integer(1)
#   payload_mass     : numeric(1) or NULL
#   out_dir          : character(1) — output directory
#
# Returns:
#   list(exit_code, stdout, stderr, out_dir)
# =============================================================================
run_tide <- function(exe, raw_files, fasta_path, enzyme_id,
                     missed_cleavages, payload_mass = NULL, out_dir) {

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  tide_enzyme <- .TIDE_ENZYME_MAP[enzyme_id]
  if (is.na(tide_enzyme)) tide_enzyme <- "trypsin"

  if (enzyme_id %in% .TIDE_NO_ENZYME_FALLBACK) {
    message(sprintf(
      "run_tide: enzyme '%s' is not natively supported by Tide; ",
      "falling back to no-enzyme (unspecific cleavage). ",
      "Consider using MS Amanda for this enzyme.", enzyme_id
    ))
  }

  mc <- as.integer(missed_cleavages)

  # Build mods-spec string: fixed CAM on C, variable Oxidation on M,
  # optional payload on C
  mods_spec <- "C+57.021464"
  var_mods   <- "1M+15.994915"
  if (!is.null(payload_mass) && is.numeric(payload_mass) && !is.na(payload_mass)) {
    var_mods <- paste0(var_mods, sprintf(",1C+%.6f", payload_mass))
  }

  index_dir  <- file.path(out_dir, "tide_index")
  all_stdout <- character(0)
  all_stderr <- character(0)

  # ── Step 1: tide-index ────────────────────────────────────────────────────
  idx_out <- tryCatch(
    system2(
      command = exe,
      args    = c(
        "tide-index",
        "--output-dir",       shQuote(normalizePath(out_dir, mustWork = FALSE)),
        "--missed-cleavages", mc,
        "--enzyme",           tide_enzyme,
        "--mods-spec",        mods_spec,
        "--variable-mod-mass-change", var_mods,
        "--min-length",       "6",
        "--max-length",       "30",
        "--min-mass",         "500",
        "--max-mass",         "5000",
        "--overwrite",        "T",
        shQuote(normalizePath(fasta_path, mustWork = FALSE)),
        shQuote(normalizePath(index_dir,  mustWork = FALSE))
      ),
      stdout = TRUE,
      stderr = TRUE,
      wait   = TRUE
    ),
    error = function(e) {
      all_stderr <<- c(all_stderr, conditionMessage(e))
      structure(character(0), status = 1L)
    }
  )

  if (is.character(idx_out)) {
    all_stdout <- c(all_stdout, idx_out)
    ec <- attr(idx_out, "status")
    if (!is.null(ec) && ec != 0L) {
      return(list(exit_code = ec, stdout = all_stdout,
                  stderr = all_stderr, out_dir = out_dir))
    }
  }

  # ── Step 2: tide-search ───────────────────────────────────────────────────
  srch_out <- tryCatch(
    system2(
      command = exe,
      args    = c(
        "tide-search",
        "--output-dir",              shQuote(normalizePath(out_dir, mustWork = FALSE)),
        "--precursor-window",        "20",
        "--precursor-window-type",   "ppm",
        "--mz-bin-width",            "0.02",
        "--top-match",               "1",
        "--pepxml-output",           "T",
        "--overwrite",               "T",
        shQuote(normalizePath(raw_files, mustWork = FALSE)),
        shQuote(normalizePath(index_dir, mustWork = FALSE))
      ),
      stdout = TRUE,
      stderr = TRUE,
      wait   = TRUE
    ),
    error = function(e) {
      all_stderr <<- c(all_stderr, conditionMessage(e))
      structure(character(0), status = 1L)
    }
  )

  last_exit <- 0L
  if (is.character(srch_out)) {
    all_stdout <- c(all_stdout, srch_out)
    ec <- attr(srch_out, "status")
    if (!is.null(ec)) last_exit <- ec
  }

  list(exit_code = last_exit, stdout = all_stdout,
       stderr = all_stderr, out_dir = out_dir)
}


# =============================================================================
# run_search_engine()
# =============================================================================
# Unified dispatcher: routes to run_msamanda() or run_tide() based on
# engine_info$engine returned by detect_search_engine().
#
# Args:
#   engine_info      : list — output of detect_search_engine()
#   raw_files        : character vector
#   fasta_path       : character(1)
#   enzyme_id        : character(1)
#   missed_cleavages : integer(1)
#   payload_mass     : numeric(1) or NULL
#   out_dir          : character(1)
#
# Returns:
#   list(exit_code, stdout, stderr, out_dir)
# =============================================================================
run_search_engine <- function(engine_info, raw_files, fasta_path,
                               enzyme_id, missed_cleavages,
                               payload_mass = NULL, out_dir) {

  if (!isTRUE(engine_info$available)) {
    stop("No search engine available: ", engine_info$error_msg)
  }

  switch(engine_info$engine,
    msamanda = run_msamanda(
      exe              = engine_info$exe,
      raw_files        = raw_files,
      fasta_path       = fasta_path,
      enzyme_id        = enzyme_id,
      missed_cleavages = missed_cleavages,
      payload_mass     = payload_mass,
      out_dir          = out_dir
    ),
    tide = run_tide(
      exe              = engine_info$exe,
      raw_files        = raw_files,
      fasta_path       = fasta_path,
      enzyme_id        = enzyme_id,
      missed_cleavages = missed_cleavages,
      payload_mass     = payload_mass,
      out_dir          = out_dir
    ),
    stop("Unknown engine: ", engine_info$engine)
  )
}


# =============================================================================
# parse_search_results()
# =============================================================================
# Parse MS Amanda / Tide / legacy MSFragger results into a normalised
# data.table. Auto-detects format from file extension and content.
#
# Args:
#   out_dir_or_file  : character(1) — directory or direct file path
#   engine           : character(1) — "msamanda", "tide", or "none"
#                      Used as a hint; format is always auto-detected.
#   score_threshold  : numeric(1) or NULL — minimum Score to retain
#   evalue_threshold : numeric(1) or NULL — maximum E-value to retain
#
# Returns:
#   data.table with columns:
#     Sequence, ModifiedSequence, Charge, RT_sec, Score,
#     Evalue, SourceFile, ScanNum, Modifications
#   Zero-row table with same schema on failure.
# =============================================================================
parse_search_results <- function(out_dir_or_file,
                                  engine           = "none",
                                  score_threshold  = NULL,
                                  evalue_threshold = NULL) {

  .empty <- function() {
    data.table::data.table(
      Sequence         = character(0),
      ModifiedSequence = character(0),
      Charge           = integer(0),
      RT_sec           = numeric(0),
      Score            = numeric(0),
      Evalue           = numeric(0),
      SourceFile       = character(0),
      ScanNum          = integer(0),
      Modifications    = character(0),
      IsDecoy          = logical(0),
      q_value          = numeric(0),
      FDR_1pct         = logical(0)
    )
  }

  # ── Resolve file path ──────────────────────────────────────────────────────
  target <- out_dir_or_file

  if (dir.exists(target)) {
    # Scan directory: prefer mzIdentML (MS Amanda), then pepXML (Tide), then tsv
    mzid_hits   <- list.files(target, pattern = "(?i)\\.mzid$|\\.mzidentml$",
                               full.names = TRUE, recursive = TRUE)
    pepxml_hits <- list.files(target, pattern = "(?i)\\.pep\\.xml$|\\.pepxml$",
                               full.names = TRUE, recursive = TRUE)
    psm_hits    <- list.files(target, pattern = "(?i)psm\\.tsv$",
                               full.names = TRUE, recursive = TRUE)
    csv_hits    <- list.files(target, pattern = "(?i)amanda.*\\.csv$|\\.csv$",
                               full.names = TRUE, recursive = TRUE)

    if      (length(mzid_hits)   > 0L) target <- mzid_hits[1L]
    else if (length(pepxml_hits) > 0L) target <- pepxml_hits[1L]
    else if (length(psm_hits)    > 0L) target <- psm_hits[1L]
    else if (length(csv_hits)    > 0L) target <- csv_hits[1L]
    else {
      message("parse_search_results: no recognised result files found in ",
              out_dir_or_file)
      return(.empty())
    }
  }

  if (!file.exists(target)) {
    message("parse_search_results: file not found: ", target)
    return(.empty())
  }

  # ── Detect format ──────────────────────────────────────────────────────────
  ext <- tolower(tools::file_ext(target))
  if (ext == "xml") {
    ext2 <- tolower(tools::file_ext(sub("\\.[^.]+$", "", target)))
    if (ext2 == "pep") ext <- "pepxml"
  }

  dt <- tryCatch({
    if (ext %in% c("mzid", "mzidentml")) {
      .parse_mzidentml_v8(target)
    } else if (ext %in% c("pepxml", "xml")) {
      .parse_pepxml_v8(target)
    } else if (ext == "tsv") {
      .parse_psm_tsv_v8(target)
    } else if (ext == "csv") {
      .parse_amanda_csv(target)
    } else {
      # Sniff first line
      first_line <- readLines(target, n = 1L, warn = FALSE)
      if (grepl("<MzIdentML|<mzIdentML", first_line, ignore.case = TRUE)) {
        .parse_mzidentml_v8(target)
      } else if (grepl("<pepXML|<msms_pipeline_analysis", first_line, ignore.case = TRUE)) {
        .parse_pepxml_v8(target)
      } else if (grepl("Spectrum\t|Peptide\t|hyperscore", first_line, ignore.case = TRUE)) {
        .parse_psm_tsv_v8(target)
      } else if (grepl("Sequence|Amanda Score", first_line, ignore.case = TRUE)) {
        .parse_amanda_csv(target)
      } else {
        message("parse_search_results: unrecognised file format: ", target)
        .empty()
      }
    }
  }, error = function(e) {
    message("parse_search_results: parse error — ", conditionMessage(e))
    .empty()
  })

  if (nrow(dt) == 0L) return(.empty())

  # ── Target-decoy FDR estimation ────────────────────────────────────────────
  # Flag decoys: any PSM whose SourceFile protein field contains a standard
  # decoy prefix ("DECOY_", "REV_", "rev_", "decoy_").  For engines that
  # embed protein accessions in the Modifications or ModifiedSequence columns
  # we fall back to a heuristic.  When no decoys are found, q_value is set
  # to NA and FDR_1pct to NA (filter is then based on raw score only).
  #
  # The q-value is the minimum FDR at which this PSM would be retained if all
  # PSMs with higher score were accepted: q = cumsum(decoy) / cumsum(target).
  # This is the "competition" FDR estimator (Käll et al. 2008, Nat. Methods).
  dt <- .add_fdr_columns(dt)

  # ── Apply score / q-value filters ─────────────────────────────────────────
  if (!is.null(score_threshold) && !is.na(score_threshold)) {
    dt <- dt[Score >= score_threshold]
  }
  if (!is.null(evalue_threshold) && !is.na(evalue_threshold)) {
    dt <- dt[Evalue <= evalue_threshold]
  }

  dt
}


# =============================================================================
# .add_fdr_columns()  [internal]
# =============================================================================
# Add IsDecoy, q_value, and FDR_1pct columns to a PSM data.table.
#
# Decoy heuristic (in priority order):
#   1. Protein accession in ProteinAcc column starts with DECOY_/REV_/rev_/decoy_
#   2. ModifiedSequence starts with "DECOY_" (some engines prefix the peptide)
#   3. No decoys detected → all PSMs treated as targets; FDR columns = NA
#
# q-value: Käll et al. (2008) Nat. Methods 5:959.
#   Sort PSMs by score DESC; q = cumsum(is_decoy) / cumsum(is_target)
#   with monotone minimum enforced from the bottom up.
# =============================================================================
.add_fdr_columns <- function(dt) {
  n <- nrow(dt)

  # ── 1. Detect decoys ────────────────────────────────────────────────────────
  is_decoy <- logical(n)

  if ("ProteinAcc" %in% names(dt)) {
    is_decoy <- grepl("^(DECOY_|REV_|rev_|decoy_)", dt$ProteinAcc, perl = TRUE)
  }
  if (!any(is_decoy) && "ModifiedSequence" %in% names(dt)) {
    is_decoy <- grepl("^(DECOY_|REV_)", dt$ModifiedSequence, perl = TRUE)
  }

  dt[, IsDecoy := is_decoy]

  # ── 2. Compute q-values only when decoys are present ──────────────────────
  if (!any(is_decoy, na.rm = TRUE)) {
    dt[, q_value  := NA_real_]
    dt[, FDR_1pct := NA]
    return(dt)
  }

  # Sort by Score descending (higher = better for Amanda/XCorr; for E-value use
  # ascending, but q-value computation is the same after sorting by best-to-worst)
  ord <- order(dt$Score, decreasing = TRUE, na.last = TRUE)
  is_decoy_sorted  <- is_decoy[ord]
  is_target_sorted <- !is_decoy_sorted

  cum_decoy  <- cumsum(is_decoy_sorted)
  cum_target <- cumsum(is_target_sorted)
  # Avoid div-by-zero; q-value capped at 1
  fdr_raw <- ifelse(cum_target == 0L, 1,
                    pmin(cum_decoy / cum_target, 1))

  # Enforce monotone non-decreasing q-values from right (best PSM never has
  # higher q than a worse one): reverse-cummin
  fdr_mono <- rev(cummin(rev(fdr_raw)))

  q_vals <- numeric(n)
  q_vals[ord] <- fdr_mono

  dt[, q_value  := round(q_vals, 5)]
  dt[, FDR_1pct := q_vals <= 0.01]
  dt
}


# =============================================================================
# Internal parsers
# =============================================================================

# ── mzIdentML (MS Amanda 3.0 output) ─────────────────────────────────────────
.parse_mzidentml_v8 <- function(path) {
  if (!requireNamespace("XML", quietly = TRUE)) {
    stop("Package 'XML' is required to parse mzIdentML files. ",
         "Install with: install.packages('XML')")
  }

  doc <- XML::xmlParse(path)
  ns  <- c(mzid = "http://psidev.info/psi/pi/mzIdentML/1.1")

  sir_nodes <- XML::getNodeSet(doc, "//mzid:SpectrumIdentificationResult", ns)
  if (length(sir_nodes) == 0L) {
    sir_nodes <- XML::getNodeSet(doc, "//SpectrumIdentificationResult")
  }
  if (length(sir_nodes) == 0L) {
    message(".parse_mzidentml_v8: no SpectrumIdentificationResult nodes in ", path)
    return(data.table::data.table(
      Sequence=character(0), ModifiedSequence=character(0),
      Charge=integer(0), RT_sec=numeric(0), Score=numeric(0),
      Evalue=numeric(0), SourceFile=character(0),
      ScanNum=integer(0), Modifications=character(0)
    ))
  }

  rows <- lapply(sir_nodes, function(sir) {
    sii_nodes <- XML::getNodeSet(sir,
      "*[local-name()='SpectrumIdentificationItem'][@rank='1']")
    if (length(sii_nodes) == 0L) return(NULL)
    sii <- sii_nodes[[1L]]

    # CV params on the SIR (retention time)
    sir_cv <- XML::getNodeSet(sir, "*[local-name()='cvParam']")
    sir_cv_vals <- setNames(
      sapply(sir_cv, function(n) XML::xmlGetAttr(n, "value", NA)),
      sapply(sir_cv, function(n) XML::xmlGetAttr(n, "name",  ""))
    )

    # CV params on the SII (scores)
    sii_cv <- XML::getNodeSet(sii, "*[local-name()='cvParam']")
    sii_cv_vals <- setNames(
      sapply(sii_cv, function(n) XML::xmlGetAttr(n, "value", NA)),
      sapply(sii_cv, function(n) XML::xmlGetAttr(n, "name",  ""))
    )

    # Peptide sequence
    pep_ref  <- XML::xmlGetAttr(sii, "peptide_ref", "")
    pep_node <- XML::getNodeSet(doc, sprintf(
      "//*[local-name()='Peptide'][@id='%s']/*[local-name()='PeptideSequence']",
      pep_ref))
    seq_str <- if (length(pep_node) > 0L) XML::xmlValue(pep_node[[1L]]) else ""

    # Modifications
    mod_nodes <- XML::getNodeSet(doc, sprintf(
      "//*[local-name()='Peptide'][@id='%s']//*[local-name()='Modification']",
      pep_ref))
    mod_str <- if (length(mod_nodes) > 0L) {
      paste(sapply(mod_nodes, function(m) {
        loc  <- XML::xmlGetAttr(m, "location", "?")
        mass <- suppressWarnings(as.numeric(
          XML::xmlGetAttr(m, "monoisotopicMassDelta", "0")))
        sprintf("pos%s:%.4f", loc, mass)
      }), collapse = "; ")
    } else "None"

    # Score: MS Amanda score CV term MS:1002319, or Amanda:AmandaScore, or fallback
    score <- suppressWarnings(as.numeric(
      sii_cv_vals["Amanda:AmandaScore"] %||%
      sii_cv_vals["MS:1002319"] %||%
      sii_cv_vals["MSFragger:hyperscore"] %||%
      sii_cv_vals["Mascot:score"] %||%
      XML::xmlGetAttr(sii, "score", NA)
    ))

    evalue <- suppressWarnings(as.numeric(
      sii_cv_vals["MS-GF:EValue"] %||%
      sii_cv_vals["Mascot:expectation value"] %||%
      XML::xmlGetAttr(sii, "p-value", NA)
    ))

    rt <- suppressWarnings(as.numeric(
      sir_cv_vals["scan start time"] %||%
      sir_cv_vals["retention time"]  %||% NA
    ))

    list(
      Sequence         = toupper(seq_str),
      ModifiedSequence = seq_str,
      Charge           = suppressWarnings(
                           as.integer(XML::xmlGetAttr(sii, "chargeState", NA))),
      RT_sec           = rt,
      Score            = score,
      Evalue           = evalue,
      SourceFile       = basename(XML::xmlGetAttr(sir, "spectraData_ref", "")),
      ScanNum          = suppressWarnings(as.integer(
                           sub(".*scan=([0-9]+).*", "\\1",
                               XML::xmlGetAttr(sir, "spectrumID", "0")))),
      Modifications    = mod_str
    )
  })

  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0L) return(data.table::data.table(
    Sequence=character(0), ModifiedSequence=character(0),
    Charge=integer(0), RT_sec=numeric(0), Score=numeric(0),
    Evalue=numeric(0), SourceFile=character(0),
    ScanNum=integer(0), Modifications=character(0)
  ))

  data.table::rbindlist(rows, fill = TRUE)
}


# ── pepXML (Tide/Crux output) ─────────────────────────────────────────────────
.parse_pepxml_v8 <- function(path) {
  if (!requireNamespace("XML", quietly = TRUE)) {
    stop("Package 'XML' is required to parse pepXML files. ",
         "Install with: install.packages('XML')")
  }

  doc  <- XML::xmlParse(path)
  hits <- XML::getNodeSet(doc,
    "//spectrum_query/search_result/search_hit[@hit_rank='1']")

  if (length(hits) == 0L) {
    message(".parse_pepxml_v8: no rank-1 hits found in ", path)
    return(data.table::data.table(
      Sequence=character(0), ModifiedSequence=character(0),
      Charge=integer(0), RT_sec=numeric(0), Score=numeric(0),
      Evalue=numeric(0), SourceFile=character(0),
      ScanNum=integer(0), Modifications=character(0)
    ))
  }

  rows <- lapply(hits, function(hit) {
    parent   <- XML::xmlParent(XML::xmlParent(hit))
    sq_attrs <- XML::xmlAttrs(parent)

    score_nodes <- XML::getNodeSet(hit, "search_score")
    scores <- setNames(
      as.numeric(sapply(score_nodes,
        function(n) XML::xmlGetAttr(n, "value", default = NA))),
      sapply(score_nodes,
        function(n) XML::xmlGetAttr(n, "name", default = ""))
    )

    mod_nodes <- XML::getNodeSet(hit, "modification_info/mod_aminoacid_mass")
    mod_str <- if (length(mod_nodes) > 0L) {
      paste(sapply(mod_nodes, function(m) {
        sprintf("pos%s:%.4f",
          XML::xmlGetAttr(m, "position", "?"),
          as.numeric(XML::xmlGetAttr(m, "mass", 0)))
      }), collapse = "; ")
    } else "None"

    # Tide XCorr is stored as "xcorr_score"; also accept "hyperscore" for legacy
    xcorr <- scores["xcorr_score"] %||% scores["hyperscore"] %||%
             scores["score"] %||% NA_real_

    list(
      Sequence         = toupper(XML::xmlGetAttr(hit, "peptide", "")),
      ModifiedSequence = XML::xmlGetAttr(hit, "peptide", ""),
      Charge           = suppressWarnings(
                           as.integer(sq_attrs["assumed_charge"])),
      RT_sec           = suppressWarnings(
                           as.numeric(sq_attrs["retention_time_sec"])),
      Score            = suppressWarnings(as.numeric(xcorr)),
      Evalue           = suppressWarnings(
                           as.numeric(scores["expect"] %||%
                                      scores["evalue"] %||% NA)),
      SourceFile       = basename(sq_attrs["spectrum"] %||% ""),
      ScanNum          = suppressWarnings(
                           as.integer(sq_attrs["start_scan"])),
      Modifications    = mod_str
    )
  })

  data.table::rbindlist(rows, fill = TRUE)
}


# ── psm.tsv (legacy MSFragger / FragPipe) ─────────────────────────────────────
.parse_psm_tsv_v8 <- function(path) {
  raw <- data.table::fread(path, sep = "\t", header = TRUE,
                            stringsAsFactors = FALSE, fill = TRUE)
  names(raw) <- tolower(gsub("[ ./]", "_", names(raw)))

  .col <- function(dt, ...) {
    candidates <- c(...)
    found <- intersect(candidates, names(dt))
    if (length(found) == 0L) return(rep(NA, nrow(dt)))
    dt[[found[1L]]]
  }

  data.table::data.table(
    Sequence         = toupper(trimws(.col(raw, "peptide", "sequence",
                                           "peptide_sequence"))),
    ModifiedSequence = trimws(.col(raw, "modified_peptide",
                                   "modified_sequence", "peptide")),
    Charge           = suppressWarnings(as.integer(
                         .col(raw, "charge", "precursor_charge"))),
    RT_sec           = suppressWarnings(as.numeric(
                         .col(raw, "retention", "retention_time", "rt"))),
    Score            = suppressWarnings(as.numeric(
                         .col(raw, "hyperscore", "score", "xcorr_score"))),
    Evalue           = suppressWarnings(as.numeric(
                         .col(raw, "expectation", "evalue", "e_value", "expect"))),
    SourceFile       = trimws(.col(raw, "spectrum_file", "spectrumfile",
                                   "raw_file", "file")),
    ScanNum          = suppressWarnings(as.integer(
                         .col(raw, "scannum", "scan_num", "scan"))),
    Modifications    = trimws(.col(raw, "assigned_modifications",
                                   "modifications", "variable_modifications"))
  )
}


# ── MS Amanda summary CSV ─────────────────────────────────────────────────────
# MS Amanda can export a summary CSV with columns including:
#   Sequence, Amanda Score, Charge, Filename, Scan Number, ...
.parse_amanda_csv <- function(path) {
  raw <- data.table::fread(path, sep = ",", header = TRUE,
                            stringsAsFactors = FALSE, fill = TRUE)
  names(raw) <- tolower(gsub("[ ./]", "_", names(raw)))

  .col <- function(dt, ...) {
    candidates <- c(...)
    found <- intersect(candidates, names(dt))
    if (length(found) == 0L) return(rep(NA, nrow(dt)))
    dt[[found[1L]]]
  }

  data.table::data.table(
    Sequence         = toupper(trimws(.col(raw, "sequence", "peptide"))),
    ModifiedSequence = trimws(.col(raw, "modified_sequence", "sequence", "peptide")),
    Charge           = suppressWarnings(as.integer(
                         .col(raw, "charge", "precursor_charge"))),
    RT_sec           = suppressWarnings(as.numeric(
                         .col(raw, "retention_time", "rt", "retention"))),
    Score            = suppressWarnings(as.numeric(
                         .col(raw, "amanda_score", "score", "amanda score"))),
    Evalue           = suppressWarnings(as.numeric(
                         .col(raw, "evalue", "e_value", "expect", "expectation"))),
    SourceFile       = trimws(.col(raw, "filename", "spectrum_file",
                                   "spectrumfile", "file")),
    ScanNum          = suppressWarnings(as.integer(
                         .col(raw, "scan_number", "scannum", "scan_num", "scan"))),
    Modifications    = trimws(.col(raw, "modifications", "assigned_modifications",
                                   "variable_modifications"))
  )
}


# =============================================================================
# Null-coalescing operator (if not already defined by app.R)
# =============================================================================
if (!exists("%||%")) {
  `%||%` <- function(a, b) if (!is.null(a) && !is.na(a) && length(a) > 0L) a else b
}
