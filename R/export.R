# export.R — Skyline CSV and Excel export formatters

# ── Skyline Transition List CSV ───────────────────────────────────────────────
# Formats transition_dt into Skyline-importable CSV
# Skyline expects specific column names for transition list import
format_skyline_csv <- function(transition_dt, unique_only = FALSE) {
  if (is.null(transition_dt) || nrow(transition_dt) == 0) return(NULL)

  dt <- transition_dt
  if (unique_only) {
    dt <- dt[dt$UniqueToADC == TRUE, ]
  }
  if (nrow(dt) == 0) return(NULL)

  # Skyline transition list column mapping
  skyline_df <- data.frame(
    `Protein Name`       = dt$ProteinName,
    `Peptide Sequence`   = dt$PeptideSequence,
    `Modified Sequence`  = dt$ModifiedSequence,
    `Precursor Charge`   = dt$PrecursorCharge,
    `Precursor Mz`       = dt$PrecursorMz,
    `Product Charge`     = dt$ProductCharge,
    `Product Mz`         = dt$ProductMz,
    `Fragment Ion`       = dt$FragmentIon,
    `Collision Energy`   = dt$CollisionEnergy,
    `Chain`              = dt$Chain,
    `Modifications`      = dt$Modifications,
    `Unique To ADC`      = ifelse(dt$UniqueToADC, "True", "False"),
    `Peptide Length`     = dt$PeptideLength,
    `Start Position`     = dt$Start,
    `End Position`       = dt$End,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  skyline_df
}

# ── Write Skyline CSV to temp file ────────────────────────────────────────────
write_skyline_csv <- function(transition_dt, unique_only = FALSE) {
  df <- format_skyline_csv(transition_dt, unique_only)
  if (is.null(df)) return(NULL)

  tmp <- tempfile(fileext = ".csv")
  write.csv(df, tmp, row.names = FALSE, quote = TRUE)
  tmp
}

# ── Excel Peptide Summary ─────────────────────────────────────────────────────
# Writes a multi-sheet Excel workbook:
#   Sheet 1: All peptides summary
#   Sheet 2: Unique peptides only
#   Sheet 3: Transition list
write_excel_summary <- function(peptides_dt, transition_dt) {
  if (is.null(peptides_dt) || nrow(peptides_dt) == 0) return(NULL)

  tmp <- tempfile(fileext = ".xlsx")

  wb <- openxlsx::createWorkbook()

  # ── Sheet 1: All Peptides ──
  openxlsx::addWorksheet(wb, "All Peptides")
  pep_df <- as.data.frame(peptides_dt[, c(
    "Chain", "ProteinName", "Sequence", "ModifiedSequence",
    "Start", "End", "Length", "ModsApplied", "ModifiedMass", "UniqueToADC"
  )])
  colnames(pep_df) <- c(
    "Chain", "Protein", "Peptide Sequence", "Modified Sequence",
    "Start", "End", "Length", "Modifications", "Monoisotopic Mass (Da)", "Unique to ADC"
  )
  pep_df$`Unique to ADC` <- ifelse(pep_df$`Unique to ADC`, "Yes", "No")

  openxlsx::writeDataTable(wb, "All Peptides", pep_df, tableStyle = "TableStyleMedium9")

  # Style header
  header_style <- openxlsx::createStyle(
    fontColour = "#FFFFFF", bgFill = "#2C3E50",
    halign = "CENTER", textDecoration = "Bold"
  )
  openxlsx::addStyle(wb, "All Peptides", header_style,
                     rows = 1, cols = seq_len(ncol(pep_df)), gridExpand = TRUE)

  # Highlight unique peptides in green
  unique_rows <- which(pep_df$`Unique to ADC` == "Yes") + 1  # +1 for header
  if (length(unique_rows) > 0) {
    green_style <- openxlsx::createStyle(bgFill = "#D5F5E3")
    openxlsx::addStyle(wb, "All Peptides", green_style,
                       rows = unique_rows, cols = seq_len(ncol(pep_df)),
                       gridExpand = TRUE, stack = TRUE)
  }

  # ── Sheet 2: Unique Peptides ──
  openxlsx::addWorksheet(wb, "Unique Peptides")
  unique_df <- pep_df[pep_df$`Unique to ADC` == "Yes", ]
  if (nrow(unique_df) > 0) {
    openxlsx::writeDataTable(wb, "Unique Peptides", unique_df,
                             tableStyle = "TableStyleMedium4")
    openxlsx::addStyle(wb, "Unique Peptides", header_style,
                       rows = 1, cols = seq_len(ncol(unique_df)), gridExpand = TRUE)
  } else {
    openxlsx::writeData(wb, "Unique Peptides",
                        data.frame(Message = "No unique peptides found."))
  }

  # ── Sheet 3: Transition List ──
  openxlsx::addWorksheet(wb, "Transition List")
  if (!is.null(transition_dt) && nrow(transition_dt) > 0) {
    trans_df <- format_skyline_csv(transition_dt, unique_only = FALSE)
    openxlsx::writeDataTable(wb, "Transition List", trans_df,
                             tableStyle = "TableStyleMedium2")
    openxlsx::addStyle(wb, "Transition List", header_style,
                       rows = 1, cols = seq_len(ncol(trans_df)), gridExpand = TRUE)
  } else {
    openxlsx::writeData(wb, "Transition List",
                        data.frame(Message = "No transitions generated yet."))
  }

  # Auto-fit columns
  for (sheet in c("All Peptides", "Unique Peptides", "Transition List")) {
    openxlsx::setColWidths(wb, sheet, cols = 1:20, widths = "auto")
  }

  openxlsx::saveWorkbook(wb, tmp, overwrite = TRUE)
  tmp
}
