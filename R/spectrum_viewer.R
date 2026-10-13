#' @title Read mzML index
#' @description Opens an mzML file with mzR and returns a header data.frame
#'   with one row per MS2 scan (retentionTime, precursorMZ, seqNum, etc.).
#'
#' @param path Character. Path to an mzML file.
#' @return data.frame with columns seqNum, retentionTime, precursorMZ,
#'   precursorIntensity, msLevel filtered to MS2 scans only.
#' @export
#' @examples
#' \dontrun{
#'   idx <- read_mzml_index("sample.mzML")
#'   head(idx)
#' }
read_mzml_index <- function(path) {
  if (!requireNamespace("mzR", quietly = TRUE))
    stop("mzR is required. Install via: BiocManager::install('mzR')")
  ms <- mzR::openMSfile(path)
  on.exit(mzR::close(ms))
  hdr <- mzR::header(ms)
  hdr[hdr$msLevel == 2, c("seqNum", "retentionTime", "precursorMZ",
                            "precursorIntensity", "peaksCount")]
}

#' @title Extract one MS2 scan
#' @description Returns the peak list for a single scan from an mzML file.
#'
#' @param path Character. Path to the mzML file.
#' @param scan_idx Integer. The seqNum of the scan to extract.
#' @return data.frame with columns mz, intensity, rel_intensity (0-100).
#' @export
#' @examples
#' \dontrun{
#'   pk <- extract_scan("sample.mzML", scan_idx = 5L)
#'   head(pk)
#' }
extract_scan <- function(path, scan_idx) {
  if (!requireNamespace("mzR", quietly = TRUE))
    stop("mzR is required. Install via: BiocManager::install('mzR')")
  ms <- mzR::openMSfile(path)
  on.exit(mzR::close(ms))
  pk <- mzR::peaks(ms, scan_idx)
  if (is.null(pk) || nrow(pk) == 0)
    return(data.frame(mz = numeric(0), intensity = numeric(0),
                      rel_intensity = numeric(0)))
  df <- as.data.frame(pk)
  names(df) <- c("mz", "intensity")
  max_int <- max(df$intensity, na.rm = TRUE)
  df$rel_intensity <- if (max_int > 0) 100 * df$intensity / max_int else 0
  df
}

#' @title Match theoretical ions to observed spectrum
#' @description Annotates observed peaks by matching them against a theoretical
#'   ion table (from \code{\link{calc_fragment_ions}}) within a ppm tolerance.
#'   Also includes neutral losses: b-H2O, b-NH3, y-H2O, y-NH3.
#'
#' @param ions_dt data.table. Output of \code{calc_fragment_ions()} with
#'   columns ion_type, ion_number, charge, mz, label.
#' @param spectrum_df data.frame. Output of \code{extract_scan()} with columns
#'   mz, intensity, rel_intensity.
#' @param ppm_tol Numeric. Match tolerance in ppm (default 10).
#' @return data.frame with all spectrum peaks plus columns: matched (logical),
#'   ion_label (character), ion_series (b/y/b_loss/y_loss).
#' @export
#' @examples
#' \dontrun{
#'   ions <- calc_fragment_ions("GPSVFPLAPSSR")
#'   spec  <- extract_scan("sample.mzML", 12L)
#'   ann   <- match_ions(ions, spec, ppm_tol = 10)
#' }
match_ions <- function(ions_dt, spectrum_df, ppm_tol = 10) {
  # Add neutral losses to theoretical ion table
  ions_aug <- as.data.frame(ions_dt)

  h2o_mass <- 18.01056
  nh3_mass <- 17.02655

  loss_rows <- list()
  for (i in seq_len(nrow(ions_aug))) {
    row <- ions_aug[i, ]
    z   <- row$charge
    for (loss_info in list(
      list(mass = h2o_mass, suffix = "-H2O"),
      list(mass = nh3_mass, suffix = "-NH3")
    )) {
      new_mz <- row$mz - loss_info$mass / z
      if (new_mz > 0)
        loss_rows <- c(loss_rows, list(data.frame(
          ion_type   = row$ion_type,
          ion_number = row$ion_number,
          charge     = z,
          mz         = new_mz,
          label      = paste0(row$label, loss_info$suffix),
          stringsAsFactors = FALSE
        )))
    }
  }
  if (length(loss_rows) > 0)
    ions_aug <- rbind(ions_aug, do.call(rbind, loss_rows))

  # Match each observed peak to the closest theoretical ion within ppm_tol
  result <- spectrum_df
  result$matched   <- FALSE
  result$ion_label <- NA_character_
  result$ion_series <- NA_character_

  for (i in seq_len(nrow(result))) {
    obs_mz <- result$mz[i]
    ppm_err <- abs(ions_aug$mz - obs_mz) / obs_mz * 1e6
    best    <- which.min(ppm_err)
    if (ppm_err[best] <= ppm_tol) {
      result$matched[i]    <- TRUE
      result$ion_label[i]  <- ions_aug$label[best]
      # classify series
      lbl <- ions_aug$ion_type[best]
      result$ion_series[i] <- if (grepl("^b", lbl)) "b" else "y"
    }
  }
  result
}

#' @title Plot annotated MS2 spectrum
#' @description Builds an interactive plotly bar chart of an annotated MS2
#'   spectrum. b-ions (blue), y-ions (red), unmatched (grey). Neutral losses
#'   are shown in lighter shades.
#'
#' @param matched_df data.frame. Output of \code{match_ions()}.
#' @param peptide_seq Character. Amino acid sequence (for title).
#' @param precursor_mz Numeric. Precursor m/z (for subtitle).
#' @param precursor_charge Integer. Precursor charge (for subtitle).
#' @return A plotly figure object.
#' @export
#' @examples
#' \dontrun{
#'   fig <- plot_annotated_spectrum(matched_df, "GPSVFPLAPSSR",
#'                                  precursor_mz = 617.35, precursor_charge = 2L)
#'   fig
#' }
plot_annotated_spectrum <- function(matched_df, peptide_seq = "",
                                    precursor_mz = NA_real_,
                                    precursor_charge = NA_integer_) {
  if (!requireNamespace("plotly", quietly = TRUE))
    stop("plotly is required. Install with: install.packages('plotly')")

  df <- matched_df
  df$color <- "#cccccc"   # unmatched grey
  df$color[df$matched & df$ion_series == "b"] <- "#2196F3"  # blue b-ions
  df$color[df$matched & df$ion_series == "y"] <- "#e53935"  # red y-ions

  # Lighter shades for neutral losses
  is_loss <- df$matched & grepl("(-H2O|-NH3)$", df$ion_label)
  df$color[is_loss & df$ion_series == "b"] <- "#90CAF9"
  df$color[is_loss & df$ion_series == "y"] <- "#EF9A9A"

  # Annotation text: only show matched labels with rel_intensity > 5
  df$text_label <- ifelse(
    df$matched & df$rel_intensity > 5,
    df$ion_label, ""
  )

  title_str <- if (nchar(peptide_seq) > 0) {
    mz_str <- if (!is.na(precursor_mz))
      sprintf(" | Precursor: %.4f (z=%s)", precursor_mz,
              if (!is.na(precursor_charge)) as.character(precursor_charge) else "?")
    else ""
    paste0(peptide_seq, mz_str)
  } else "MS2 Spectrum"

  plotly::plot_ly(
    data        = df,
    x           = ~mz,
    y           = ~rel_intensity,
    type        = "bar",
    marker      = list(color = ~color, line = list(width = 0)),
    text        = ~text_label,
    textposition = "outside",
    hovertemplate = paste0(
      "<b>m/z:</b> %{x:.4f}<br>",
      "<b>Rel. Int.:</b> %{y:.1f}%<br>",
      "<b>Ion:</b> %{text}<extra></extra>"
    )
  ) |>
    plotly::layout(
      title      = list(text = title_str, font = list(size = 13)),
      xaxis      = list(title = "m/z", zeroline = FALSE),
      yaxis      = list(title = "Relative Intensity (%)", range = c(0, 120)),
      bargap     = 0,
      showlegend = FALSE,
      plot_bgcolor  = "#fafafa",
      paper_bgcolor = "#ffffff",
      margin     = list(t = 60, b = 50, l = 60, r = 20)
    )
}

#' @title Find closest precursor scan
#' @description Searches an mzML index for the MS2 scan whose precursor m/z
#'   is closest to a target value within a given ppm window.
#'
#' @param mzml_index data.frame. Output of \code{read_mzml_index()}.
#' @param target_mz Numeric. Target precursor m/z.
#' @param charge Integer. Expected precursor charge state.
#' @param ppm_tol Numeric. Tolerance in ppm (default 20).
#' @return Integer seqNum of the best matching scan, or NA if none found.
#' @export
#' @examples
#' \dontrun{
#'   idx  <- read_mzml_index("sample.mzML")
#'   scan <- find_precursor_scan(idx, target_mz = 617.35, charge = 2L)
#' }
find_precursor_scan <- function(mzml_index, target_mz, charge = NA_integer_,
                                ppm_tol = 20) {
  if (nrow(mzml_index) == 0 || is.na(target_mz)) return(NA_integer_)
  ppm_err <- abs(mzml_index$precursorMZ - target_mz) / target_mz * 1e6
  best    <- which.min(ppm_err)
  if (ppm_err[best] <= ppm_tol) mzml_index$seqNum[best] else NA_integer_
}
