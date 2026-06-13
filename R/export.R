# =============================================================================
# export.R — ADC Peptide Mapper v0.8
# =============================================================================
# Sourced by the Shiny app (server.R). No library() calls here.
# Depends on:
#   data.table      — loaded by app
#   openxlsx        — loaded by app (for write_excel_summary)
#   transitions.R   — sourced before this file; provides select_top_ions(),
#                     calc_ce_instrument()
#
# Public API
# ----------
#   format_skyline_csv(transition_dt, unique_only)
#   write_skyline_csv(transition_dt, unique_only)
#   format_thermo(transition_dt, unique_only, top_n)
#   format_sciex(transition_dt, unique_only, top_n)
#   format_bruker(transition_dt, unique_only, top_n)
#   format_agilent(transition_dt, unique_only, top_n)
#   format_waters(transition_dt, unique_only, top_n)
#   write_instrument_csv(transition_dt, instrument, unique_only, top_n)
#   write_excel_summary(transition_dt, peptides_dt, unique_only)
#
# v0.8 changes vs v0.6
# ---------------------
#   - format_skyline_csv / write_skyline_csv: added ADCName column
#   - Added format_thermo(), format_sciex(), format_bruker(),
#     format_agilent(), format_waters()
#   - Added write_instrument_csv() dispatcher
#   - write_excel_summary(): added 4th sheet "Instrument Exports"
# =============================================================================


# =============================================================================
# INTERNAL HELPERS
# =============================================================================

# -----------------------------------------------------------------------------
# .apply_unique_filter()  [internal]
# -----------------------------------------------------------------------------
# Optionally filter a transition data.table to unique-to-ADC peptides only.
#
# Args:
#   dt          : data.table — transition list (must have UniqueToADC column)
#   unique_only : logical(1) — if TRUE, keep only rows where UniqueToADC == TRUE
#
# Returns:
#   data.table — filtered (or unchanged) copy
# -----------------------------------------------------------------------------
.apply_unique_filter <- function(dt, unique_only) {
  if (isTRUE(unique_only) && "UniqueToADC" %in% names(dt)) {
    dt <- dt[UniqueToADC == TRUE]
  }
  dt
}


# -----------------------------------------------------------------------------
# .apply_top_n_per_precursor()  [internal]
# -----------------------------------------------------------------------------
# Keep the top N fragment ions per (PeptideSequence × PrecursorCharge) group,
# ranked by ProductMz descending (largest product ions first).
#
# This is the standard top-N filtering applied by all instrument formatters.
# Ranking by largest m/z preferentially retains high-mass y-ions, which
# typically give the best signal-to-noise in MRM experiments.
#
# Args:
#   dt    : data.table — must have columns PeptideSequence, PrecursorCharge,
#                        ProductMz
#   top_n : integer(1) — maximum ions per precursor group
#
# Returns:
#   data.table — filtered, preserving original column order
# -----------------------------------------------------------------------------
.apply_top_n_per_precursor <- function(dt, top_n) {

  top_n <- as.integer(top_n)
  if (top_n <= 0L || nrow(dt) == 0L) return(dt)

  # Use data.table grouping for efficiency
  # .SD[order(-ProductMz)][seq_len(min(.N, top_n))] selects top N per group
  dt[
    dt[, .I[order(-ProductMz)[seq_len(min(.N, top_n))]],
       by = .(PeptideSequence, PrecursorCharge)]$V1
  ]
}


# =============================================================================
# SKYLINE FORMAT
# =============================================================================

# -----------------------------------------------------------------------------
# format_skyline_csv()
# -----------------------------------------------------------------------------
# Format the transition list for import into Skyline (Proteomics Research
# Software, MacCoss Lab). Returns a data.frame ready for write.csv().
#
# Skyline expects a specific column set for its "Transition List" import.
# See: https://skyline.ms/wiki/home/software/Skyline/page.view?name=import_transition_list
#
# Args:
#   transition_dt : data.table — full transition list from generate_transition_list()
#   unique_only   : logical(1) — if TRUE, export only UniqueToADC == TRUE rows
#
# Returns:
#   data.frame with Skyline transition list columns (no top-N filtering —
#   Skyline handles its own peak picking)
# -----------------------------------------------------------------------------
format_skyline_csv <- function(transition_dt, unique_only = FALSE) {

  dt <- .apply_unique_filter(transition_dt, unique_only)

  if (nrow(dt) == 0L) {
    warning("format_skyline_csv: no rows after filtering — returning empty data.frame.")
    return(data.frame())
  }

  # Build Skyline-compatible data.frame
  # Column names match Skyline's expected transition list headers exactly
  out <- data.frame(
    # ADC product name (v0.6 addition)
    "ADC Name"           = dt$ADCName,

    # Protein / molecule identifiers
    "Protein Name"       = dt$ProteinName,
    "Peptide Sequence"   = dt$PeptideSequence,
    "Modified Sequence"  = dt$ModifiedSequence,

    # Precursor
    "Precursor Charge"   = dt$PrecursorCharge,
    "Precursor Mz"       = round(dt$PrecursorMz, 6),

    # Product / fragment
    "Product Charge"     = dt$ProductCharge,
    "Product Mz"         = round(dt$ProductMz, 6),
    "Fragment Ion"       = dt$FragmentIon,

    # Acquisition parameters
    "Collision Energy"   = dt$CollisionEnergy,

    # Annotation columns
    "Modifications"      = dt$Modifications,
    "Unique To ADC"      = dt$UniqueToADC,
    "Peptide Length"     = dt$PeptideLength,
    "Start"              = dt$Start,
    "End"                = dt$End,
    "Enzyme"             = dt$Enzyme,

    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  out
}


# -----------------------------------------------------------------------------
# write_skyline_csv()
# -----------------------------------------------------------------------------
# Write the Skyline transition list to a temporary CSV file.
#
# Args:
#   transition_dt : data.table — full transition list
#   unique_only   : logical(1) — filter to unique peptides only
#
# Returns:
#   character(1) — path to the temporary CSV file
# -----------------------------------------------------------------------------
write_skyline_csv <- function(transition_dt, unique_only = FALSE) {

  out_df   <- format_skyline_csv(transition_dt, unique_only = unique_only)
  tmp_path <- tempfile(pattern = "skyline_transitions_", fileext = ".csv")

  write.csv(out_df, file = tmp_path, row.names = FALSE, quote = TRUE)

  message(sprintf("write_skyline_csv: wrote %d rows to %s", nrow(out_df), tmp_path))
  tmp_path
}


# =============================================================================
# THERMO FORMAT
# =============================================================================

# -----------------------------------------------------------------------------
# format_thermo()
# -----------------------------------------------------------------------------
# Format the transition list for Thermo Xcalibur / TSQ Altis / TSQ Quantis.
#
# The output matches the "Compound List" (SRM method) format expected by
# Xcalibur's TSQ Method Editor. Each row represents one MRM transition.
#
# Reaction Monitor: a semicolon-delimited list of all fragment ion labels
# for a given precursor (grouped by PeptideSequence + PrecursorCharge).
# This is written once per unique precursor row.
#
# Args:
#   transition_dt : data.table — full transition list
#   unique_only   : logical(1) — filter to UniqueToADC == TRUE rows
#   top_n         : integer(1) — keep top N ions per precursor (default 5),
#                   ranked by ProductMz descending
#
# Returns:
#   data.frame with Thermo Xcalibur MRM columns
# -----------------------------------------------------------------------------
format_thermo <- function(transition_dt, unique_only = FALSE, top_n = 5L) {

  dt <- .apply_unique_filter(transition_dt, unique_only)
  dt <- .apply_top_n_per_precursor(dt, top_n)

  if (nrow(dt) == 0L) {
    warning("format_thermo: no rows after filtering — returning empty data.frame.")
    return(data.frame())
  }

  # Recalculate CE using Thermo-specific formula
  dt <- data.table::copy(dt)
  dt[, CollisionEnergy := mapply(
    calc_ce_instrument,
    precursor_mz = PrecursorMz,
    charge       = PrecursorCharge,
    MoreArgs     = list(instrument = "thermo")
  )]

  # Compound Name: "SEQUENCE_z2", "SEQUENCE_z3", etc.
  dt[, CompoundName := paste(PeptideSequence, PrecursorCharge, sep = "_z")]

  # Reaction Monitor: collapse all fragment ion labels per precursor into one
  # semicolon-delimited string (one row per transition, same value per group)
  dt[, ReactionMonitor := paste(FragmentIon, collapse = ";"),
     by = .(PeptideSequence, PrecursorCharge)]

  out <- data.frame(
    "Compound Name"      = dt$CompoundName,
    "Precursor (m/z)"    = round(dt$PrecursorMz, 6),
    "Product (m/z)"      = round(dt$ProductMz, 6),
    "Collision Energy"   = dt$CollisionEnergy,
    "Start Time (min)"   = 0,          # placeholder — set in Xcalibur method
    "Stop Time (min)"    = 0,          # placeholder
    "Polarity"           = "Positive",
    "Trigger"            = 0.0,
    "Reaction Monitor"   = dt$ReactionMonitor,
    check.names          = FALSE,
    stringsAsFactors     = FALSE
  )

  out
}


# =============================================================================
# SCIEX FORMAT
# =============================================================================

# -----------------------------------------------------------------------------
# format_sciex()
# -----------------------------------------------------------------------------
# Format the transition list for SCIEX Analyst MRM method import.
#
# The output matches the "MRM Transition" table format used by SCIEX Analyst
# software for TripleTOF, QTRAP, and Triple Quad instruments.
#
# Args:
#   transition_dt : data.table — full transition list
#   unique_only   : logical(1) — filter to UniqueToADC == TRUE rows
#   top_n         : integer(1) — keep top N ions per precursor (default 5),
#                   ranked by ProductMz descending
#
# Returns:
#   data.frame with SCIEX Analyst MRM columns
# -----------------------------------------------------------------------------
format_sciex <- function(transition_dt, unique_only = FALSE, top_n = 5L) {

  dt <- .apply_unique_filter(transition_dt, unique_only)
  dt <- .apply_top_n_per_precursor(dt, top_n)

  if (nrow(dt) == 0L) {
    warning("format_sciex: no rows after filtering — returning empty data.frame.")
    return(data.frame())
  }

  # Recalculate CE using SCIEX-specific formula (charge-dependent)
  dt <- data.table::copy(dt)
  dt[, CollisionEnergy := mapply(
    calc_ce_instrument,
    precursor_mz = PrecursorMz,
    charge       = PrecursorCharge,
    MoreArgs     = list(instrument = "sciex")
  )]

  # Transition name: "SEQUENCE_FragmentIon_z2"
  dt[, TransitionName := paste(PeptideSequence, FragmentIon, PrecursorCharge, sep = "_")]

  out <- data.frame(
    "Q1 Mass (Da)"               = round(dt$PrecursorMz, 6),
    "Q3 Mass (Da)"               = round(dt$ProductMz, 6),
    "Time (msec)"                = 100,          # default dwell time (ms)
    "Name"                       = dt$TransitionName,
    "Collision Energy (V)"       = dt$CollisionEnergy,
    "Declustering Potential (V)" = 80,           # default DP for peptides
    "Cell Exit Potential (V)"    = 15,           # default CXP
    check.names                  = FALSE,
    stringsAsFactors             = FALSE
  )

  out
}


# =============================================================================
# BRUKER FORMAT
# =============================================================================

# -----------------------------------------------------------------------------
# format_bruker()
# -----------------------------------------------------------------------------
# Format the transition list for Bruker timsControl / QTOF / EVOQ instruments.
#
# The output matches the "Compound List" format used by Bruker's timsControl
# and otofControl software for targeted MRM / PRM acquisition.
#
# Args:
#   transition_dt : data.table — full transition list
#   unique_only   : logical(1) — filter to UniqueToADC == TRUE rows
#   top_n         : integer(1) — keep top N ions per precursor (default 5),
#                   ranked by ProductMz descending
#
# Returns:
#   data.frame with Bruker timsControl compound list columns
# -----------------------------------------------------------------------------
format_bruker <- function(transition_dt, unique_only = FALSE, top_n = 5L) {

  dt <- .apply_unique_filter(transition_dt, unique_only)
  dt <- .apply_top_n_per_precursor(dt, top_n)

  if (nrow(dt) == 0L) {
    warning("format_bruker: no rows after filtering — returning empty data.frame.")
    return(data.frame())
  }

  # Recalculate CE using Bruker-specific formula
  dt <- data.table::copy(dt)
  dt[, CollisionEnergy := mapply(
    calc_ce_instrument,
    precursor_mz = PrecursorMz,
    charge       = PrecursorCharge,
    MoreArgs     = list(instrument = "bruker")
  )]

  # Adduct notation: [M+2H]2+, [M+3H]3+, [M+4H]4+
  dt[, Adduct := paste0("[M+", PrecursorCharge, "H]", PrecursorCharge, "+")]

  out <- data.frame(
    "Compound"     = dt$PeptideSequence,
    "Formula"      = "",                        # molecular formula unknown
    "Adduct"       = dt$Adduct,
    "m/z"          = round(dt$PrecursorMz, 6),
    "z"            = dt$PrecursorCharge,
    "CID"          = dt$CollisionEnergy,        # CID = collision-induced dissociation energy
    "Fragment m/z" = round(dt$ProductMz, 6),
    "Fragment Name"= dt$FragmentIon,
    check.names    = FALSE,
    stringsAsFactors = FALSE
  )

  out
}


# =============================================================================
# AGILENT FORMAT
# =============================================================================

# -----------------------------------------------------------------------------
# format_agilent()
# -----------------------------------------------------------------------------
# Format the transition list for Agilent MassHunter / Triple Quad instruments.
#
# The output matches the "MRM Transition" table format used by Agilent
# MassHunter Workstation for 6400-series QQQ instruments.
#
# Args:
#   transition_dt : data.table — full transition list
#   unique_only   : logical(1) — filter to UniqueToADC == TRUE rows
#   top_n         : integer(1) — keep top N ions per precursor (default 5),
#                   ranked by ProductMz descending
#
# Returns:
#   data.frame with Agilent MassHunter MRM columns
# -----------------------------------------------------------------------------
format_agilent <- function(transition_dt, unique_only = FALSE, top_n = 5L) {

  dt <- .apply_unique_filter(transition_dt, unique_only)
  dt <- .apply_top_n_per_precursor(dt, top_n)

  if (nrow(dt) == 0L) {
    warning("format_agilent: no rows after filtering — returning empty data.frame.")
    return(data.frame())
  }

  # Recalculate CE using Agilent-specific formula
  dt <- data.table::copy(dt)
  dt[, CollisionEnergy := mapply(
    calc_ce_instrument,
    precursor_mz = PrecursorMz,
    charge       = PrecursorCharge,
    MoreArgs     = list(instrument = "agilent")
  )]

  # Transition name: "SEQUENCE_FragmentIon"
  dt[, TransitionName := paste(PeptideSequence, FragmentIon, sep = "_")]

  out <- data.frame(
    "Precursor Ion"           = round(dt$PrecursorMz, 6),
    "MS1 Res"                 = "Unit",          # unit resolution on Q1
    "Product Ion"             = round(dt$ProductMz, 6),
    "MS2 Res"                 = "Unit",          # unit resolution on Q3
    "Dwell"                   = 20,              # dwell time in ms
    "Fragmentor"              = 135,             # fragmentor voltage (V), typical for peptides
    "Collision Energy"        = dt$CollisionEnergy,
    "Cell Accelerator Voltage"= 4,               # CAV (V), Agilent default
    "Polarity"                = "Positive",
    "Name"                    = dt$TransitionName,
    check.names               = FALSE,
    stringsAsFactors          = FALSE
  )

  out
}


# =============================================================================
# WATERS FORMAT
# =============================================================================

# -----------------------------------------------------------------------------
# format_waters()
# -----------------------------------------------------------------------------
# Format the transition list for Waters MassLynx / Xevo TQ instruments.
#
# The output matches the "MRM Transition" table format used by Waters
# MassLynx software for Xevo TQ-S and Xevo TQ-XS instruments.
#
# Args:
#   transition_dt : data.table — full transition list
#   unique_only   : logical(1) — filter to UniqueToADC == TRUE rows
#   top_n         : integer(1) — keep top N ions per precursor (default 5),
#                   ranked by ProductMz descending
#
# Returns:
#   data.frame with Waters MassLynx MRM columns
# -----------------------------------------------------------------------------
format_waters <- function(transition_dt, unique_only = FALSE, top_n = 5L) {

  dt <- .apply_unique_filter(transition_dt, unique_only)
  dt <- .apply_top_n_per_precursor(dt, top_n)

  if (nrow(dt) == 0L) {
    warning("format_waters: no rows after filtering — returning empty data.frame.")
    return(data.frame())
  }

  # Recalculate CE using Waters-specific formula
  dt <- data.table::copy(dt)
  dt[, CollisionEnergy := mapply(
    calc_ce_instrument,
    precursor_mz = PrecursorMz,
    charge       = PrecursorCharge,
    MoreArgs     = list(instrument = "waters")
  )]

  # Component name: "SEQUENCE_FragmentIon_z2"
  dt[, ComponentName := paste(PeptideSequence, FragmentIon, PrecursorCharge, sep = "_")]

  out <- data.frame(
    "Parent"          = round(dt$PrecursorMz, 6),
    "Daughter"        = round(dt$ProductMz, 6),
    "Dwell Time"      = 0.05,           # dwell time in seconds (50 ms)
    "Cone"            = 35,             # cone voltage (V), typical for peptides
    "Collision Energy"= dt$CollisionEnergy,
    "Function"        = 1,              # MassLynx function number (1 = MRM)
    "Component Name"  = dt$ComponentName,
    check.names       = FALSE,
    stringsAsFactors  = FALSE
  )

  out
}


# =============================================================================
# UNIFIED INSTRUMENT CSV WRITER
# =============================================================================

# -----------------------------------------------------------------------------
# write_instrument_csv()
# -----------------------------------------------------------------------------
# Dispatcher: format the transition list for a specified instrument and write
# to a temporary CSV file. Returns the path to the file for Shiny download.
#
# Args:
#   transition_dt : data.table — full transition list from generate_transition_list()
#   instrument    : character(1) — target instrument platform. One of:
#                     "skyline"  Skyline (no top-N filtering)
#                     "thermo"   Thermo Xcalibur / TSQ Altis
#                     "sciex"    SCIEX Analyst / TripleTOF / QTRAP
#                     "bruker"   Bruker timsControl / QTOF / EVOQ
#                     "agilent"  Agilent MassHunter / QQQ
#                     "waters"   Waters MassLynx / Xevo TQ
#   unique_only   : logical(1) — if TRUE, export only UniqueToADC == TRUE rows
#   top_n         : integer(1) — top N fragment ions per precursor (default 5).
#                   Ignored for "skyline" (no top-N filtering applied).
#
# Returns:
#   character(1) — absolute path to the temporary CSV file.
#                  The file is created in the system temp directory and will
#                  be cleaned up by the OS or Shiny session end.
#
# Errors:
#   Stops with an informative message if `instrument` is not recognised.
# -----------------------------------------------------------------------------
write_instrument_csv <- function(transition_dt,
                                 instrument,
                                 unique_only = FALSE,
                                 top_n       = 5L) {

  instrument <- tolower(trimws(instrument))

  # Supported instruments and their formatter functions
  supported <- c("skyline", "thermo", "sciex", "bruker", "agilent", "waters")

  if (!instrument %in% supported) {
    stop(sprintf(
      "write_instrument_csv: unknown instrument '%s'. Supported: %s",
      instrument, paste(supported, collapse = ", ")
    ))
  }

  # Format the data using the appropriate formatter
  out_df <- switch(instrument,
    "skyline" = format_skyline_csv(transition_dt, unique_only = unique_only),
    "thermo"  = format_thermo(transition_dt,  unique_only = unique_only, top_n = top_n),
    "sciex"   = format_sciex(transition_dt,   unique_only = unique_only, top_n = top_n),
    "bruker"  = format_bruker(transition_dt,  unique_only = unique_only, top_n = top_n),
    "agilent" = format_agilent(transition_dt, unique_only = unique_only, top_n = top_n),
    "waters"  = format_waters(transition_dt,  unique_only = unique_only, top_n = top_n)
  )

  # Write to a named temp file
  tmp_path <- tempfile(
    pattern = paste0(instrument, "_transitions_"),
    fileext = ".csv"
  )

  write.csv(out_df, file = tmp_path, row.names = FALSE, quote = TRUE)

  message(sprintf(
    "write_instrument_csv [%s]: wrote %d rows to %s",
    instrument, nrow(out_df), tmp_path
  ))

  tmp_path
}


# =============================================================================
# EXCEL SUMMARY EXPORT
# =============================================================================

# -----------------------------------------------------------------------------
# write_excel_summary()
# -----------------------------------------------------------------------------
# Write a multi-sheet Excel workbook summarising the ADC peptide mapping
# results and transition list.
#
# Sheets:
#   1. "Transition List"    — full transition list (all ions, all charges)
#   2. "Peptide Summary"    — one row per unique peptide with coverage stats
#   3. "Unique Peptides"    — subset of Peptide Summary for UniqueToADC == TRUE
#   4. "Instrument Exports" — reference table of supported instruments and
#                             their export column formats (v0.6 addition)
#
# Args:
#   transition_dt : data.table — full transition list from generate_transition_list()
#   peptides_dt   : data.table — peptide table from the digest step (one row
#                   per peptide; used to build the Peptide Summary sheet)
#   unique_only   : logical(1) — if TRUE, the "Transition List" sheet contains
#                   only UniqueToADC == TRUE rows (Peptide Summary sheets are
#                   always complete)
#
# Returns:
#   character(1) — path to the temporary .xlsx file
# -----------------------------------------------------------------------------
write_excel_summary <- function(transition_dt,
                                peptides_dt,
                                unique_only = FALSE) {

  # ---- Sheet 1: Transition List --------------------------------------------

  tl_sheet <- .apply_unique_filter(transition_dt, unique_only)
  tl_df    <- as.data.frame(tl_sheet)

  # ---- Sheet 2: Peptide Summary --------------------------------------------
  # One row per unique peptide; summarise charge states and ion counts

  pep_summary <- transition_dt[, .(
    ProteinName      = ProteinName[1L],
    Chain            = Chain[1L],
    ModifiedSequence = ModifiedSequence[1L],
    Modifications    = Modifications[1L],
    UniqueToADC      = UniqueToADC[1L],
    PeptideLength    = PeptideLength[1L],
    Start            = Start[1L],
    End              = End[1L],
    Enzyme           = Enzyme[1L],
    ChargeStates     = paste(sort(unique(PrecursorCharge)), collapse = "/"),
    TotalTransitions = .N,
    UniqueFragments  = length(unique(FragmentIon))
  ), by = .(ADCName, PeptideSequence)]

  pep_df <- as.data.frame(pep_summary)

  # ---- Sheet 3: Unique Peptides --------------------------------------------

  uniq_df <- pep_df[pep_df$UniqueToADC == TRUE, , drop = FALSE]

  # ---- Sheet 4: Instrument Exports reference table -------------------------

  instrument_ref <- data.frame(
    Instrument = c(
      "Skyline",
      "Thermo Xcalibur / TSQ Altis",
      "SCIEX Analyst / TripleTOF / QTRAP",
      "Bruker timsControl / QTOF / EVOQ",
      "Agilent MassHunter / QQQ",
      "Waters MassLynx / Xevo TQ"
    ),
    Key = c(
      "skyline", "thermo", "sciex", "bruker", "agilent", "waters"
    ),
    TopNFiltering = c(
      "None (all ions exported)",
      "Top N by ProductMz desc",
      "Top N by ProductMz desc",
      "Top N by ProductMz desc",
      "Top N by ProductMz desc",
      "Top N by ProductMz desc"
    ),
    DefaultTopN = c(
      "N/A", "5", "5", "5", "5", "5"
    ),
    CEFormula = c(
      "SCIEX (charge-dependent)",
      "0.0340 * mz + 3.0",
      "0.0448/0.0533/0.0580 * mz - 2.0 (z=2/3/4)",
      "0.0380 * mz + 2.5",
      "0.0360 * mz + 4.0",
      "0.0350 * mz + 3.5"
    ),
    OutputColumns = c(
      paste(c("ADC Name", "Protein Name", "Peptide Sequence", "Modified Sequence",
              "Precursor Charge", "Precursor Mz", "Product Charge", "Product Mz",
              "Fragment Ion", "Collision Energy", "Modifications", "Unique To ADC",
              "Peptide Length", "Start", "End", "Enzyme"), collapse = " | "),
      paste(c("Compound Name", "Precursor (m/z)", "Product (m/z)",
              "Collision Energy", "Start Time (min)", "Stop Time (min)",
              "Polarity", "Trigger", "Reaction Monitor"), collapse = " | "),
      paste(c("Q1 Mass (Da)", "Q3 Mass (Da)", "Time (msec)", "Name",
              "Collision Energy (V)", "Declustering Potential (V)",
              "Cell Exit Potential (V)"), collapse = " | "),
      paste(c("Compound", "Formula", "Adduct", "m/z", "z",
              "CID", "Fragment m/z", "Fragment Name"), collapse = " | "),
      paste(c("Precursor Ion", "MS1 Res", "Product Ion", "MS2 Res",
              "Dwell", "Fragmentor", "Collision Energy",
              "Cell Accelerator Voltage", "Polarity", "Name"), collapse = " | "),
      paste(c("Parent", "Daughter", "Dwell Time", "Cone",
              "Collision Energy", "Function", "Component Name"), collapse = " | ")
    ),
    Notes = c(
      "No top-N filtering; Skyline performs its own peak picking. ADCName column added in v0.6.",
      "Reaction Monitor column contains semicolon-delimited fragment ion labels per precursor. Start/Stop Time = 0 (placeholder).",
      "Declustering Potential = 80 V (default). Cell Exit Potential = 15 V (default). Dwell = 100 ms.",
      "Formula column is empty (molecular formula not computed). Adduct notation: [M+zH]z+.",
      "Fragmentor = 135 V (default). Cell Accelerator Voltage = 4 V. Dwell = 20 ms.",
      "Cone = 35 V (default). Dwell Time = 0.05 s (50 ms). Function = 1 (MRM)."
    ),
    stringsAsFactors = FALSE
  )

  # ---- Assemble workbook ---------------------------------------------------

  wb <- openxlsx::createWorkbook()

  # Helper: add a sheet with auto-column widths and a frozen header row
  .add_sheet <- function(wb, sheet_name, df) {
    openxlsx::addWorksheet(wb, sheet_name)
    openxlsx::writeData(wb, sheet_name, df, headerStyle = openxlsx::createStyle(
      textDecoration = "bold",
      fgFill         = "#D9E1F2",
      border         = "Bottom",
      borderColour   = "#4472C4"
    ))
    openxlsx::freezePane(wb, sheet_name, firstRow = TRUE)
    openxlsx::setColWidths(wb, sheet_name, cols = seq_len(ncol(df)), widths = "auto")
  }

  .add_sheet(wb, "Transition List",    tl_df)
  .add_sheet(wb, "Peptide Summary",    pep_df)
  .add_sheet(wb, "Unique Peptides",    uniq_df)
  .add_sheet(wb, "Instrument Exports", instrument_ref)

  # Write to temp file
  tmp_path <- tempfile(pattern = "adc_peptide_mapper_", fileext = ".xlsx")
  openxlsx::saveWorkbook(wb, file = tmp_path, overwrite = TRUE)

  message(sprintf(
    "write_excel_summary: wrote workbook (%d transitions, %d peptides) to %s",
    nrow(tl_df), nrow(pep_df), tmp_path
  ))

  tmp_path
}
