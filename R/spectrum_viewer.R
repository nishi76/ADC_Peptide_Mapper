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

# ── Full-run index & TIC ─────────────────────────────────────────────────────

#' @title Read full mzML scan index (MS1 + MS2)
#' @description Returns a data.frame with one row per scan for all MS levels,
#'   suitable for building a TIC/BPC chromatogram and navigating the raw file.
#'
#' @param path Character. Path to an mzML file.
#' @return data.frame with columns seqNum, msLevel, rt_min, totIonCurrent,
#'   basePeakMZ, basePeakIntensity, precursorMZ, precursorCharge, peaksCount.
#' @export
read_mzml_full_index <- function(path) {
  if (!requireNamespace("mzR", quietly = TRUE))
    stop("mzR is required. Install via: BiocManager::install('mzR')")
  ms <- mzR::openMSfile(path)
  on.exit(mzR::close(ms))
  h <- mzR::header(ms)
  data.frame(
    seqNum            = h$seqNum,
    msLevel           = h$msLevel,
    rt_min            = h$retentionTime / 60,
    totIonCurrent     = h$totIonCurrent,
    basePeakMZ        = h$basePeakMZ,
    basePeakIntensity = h$basePeakIntensity,
    precursorMZ       = ifelse(h$msLevel > 1L, h$precursorMZ,      NA_real_),
    precursorCharge   = ifelse(h$msLevel > 1L, h$precursorCharge,  NA_integer_),
    precursorScanNum  = ifelse(h$msLevel > 1L, h$precursorScanNum, NA_integer_),
    peaksCount        = h$peaksCount,
    stringsAsFactors  = FALSE
  )
}

#' @title Plot base-peak chromatogram with MS2 scan markers
#' @description Builds an interactive plotly BPC/TIC trace. MS2 events are
#'   drawn as triangle markers at y = 0. Clicking any point selects that
#'   region. A vertical dashed line marks the selected scan.
#'
#' @param tic_df data.frame. Output of \code{read_mzml_full_index()}.
#' @param selected_scan Integer or NULL. seqNum of the currently selected scan.
#' @param matched_scans Integer vector or NULL. seqNums to highlight in red.
#' @param source_id Character. Plotly source id for click events.
#' @return A plotly figure object.
#' @export
plot_tic <- function(tic_df, selected_scan = NULL, matched_scans = NULL,
                     source_id = "sv2_tic") {
  if (!requireNamespace("plotly", quietly = TRUE))
    stop("plotly is required.")

  ms1 <- tic_df[tic_df$msLevel == 1L, ]
  ms2 <- tic_df[tic_df$msLevel == 2L, ]

  bpi_max <- max(ms1$basePeakIntensity, na.rm = TRUE)
  if (!is.finite(bpi_max) || bpi_max == 0) bpi_max <- 1
  ms1$bpi_norm <- ms1$basePeakIntensity / bpi_max * 100

  # MS2 triangle colours
  ms2_col <- rep("rgba(150,150,150,0.30)", nrow(ms2))
  if (!is.null(matched_scans) && length(matched_scans) > 0)
    ms2_col[ms2$seqNum %in% matched_scans] <- "rgba(231,76,60,0.75)"
  if (!is.null(selected_scan) && !is.na(selected_scan) && length(selected_scan) > 0)
    ms2_col[ms2$seqNum == selected_scan] <- "#c0392b"

  fig <- plotly::plot_ly(source = source_id) |>
    plotly::add_trace(
      data = ms1, x = ~rt_min, y = ~bpi_norm,
      type = "scatter", mode = "lines",
      line = list(color = "#2980b9", width = 1.5),
      name = "BPC", showlegend = FALSE,
      hovertemplate = "RT: %{x:.2f} min | BPI: %{text}<extra></extra>",
      text = format(ms1$basePeakIntensity, scientific = TRUE, digits = 3)
    ) |>
    plotly::add_trace(
      data = ms2, x = ~rt_min, y = rep(0, nrow(ms2)),
      type = "scatter", mode = "markers",
      marker = list(symbol = "triangle-up", size = 7, color = ms2_col,
                    line = list(width = 0)),
      customdata = ms2$seqNum,
      name = "MS2", showlegend = FALSE,
      hovertemplate = paste0(
        "Scan #%{customdata}<br>",
        "RT: %{x:.2f} min<br>",
        "Prec: %{text} m/z<extra></extra>"
      ),
      text = round(ms2$precursorMZ, 4)
    )

  shapes <- list()
  if (!is.null(selected_scan) && !is.na(selected_scan) && length(selected_scan) > 0) {
    sel_row <- tic_df[tic_df$seqNum == selected_scan, ]
    if (nrow(sel_row) > 0)
      shapes <- list(list(
        type = "line",
        x0 = sel_row$rt_min[1], x1 = sel_row$rt_min[1],
        y0 = 0, y1 = 1, yref = "paper",
        line = list(color = "#c0392b", dash = "dash", width = 1.5)
      ))
  }

  fig |> plotly::layout(
    xaxis  = list(title = "Retention Time (min)", showgrid = FALSE, zeroline = FALSE),
    yaxis  = list(title = "BPC (% max)", range = c(-8, 110),
                  showgrid = TRUE, gridcolor = "#f0f0f0", zeroline = FALSE),
    margin = list(l = 60, r = 15, t = 12, b = 45),
    plot_bgcolor  = "#ffffff",
    paper_bgcolor = "#ffffff",
    shapes = shapes
  )
}

#' @title Plot annotated MS2 spectrum (publication style)
#' @description Vertical-stick plot with colour-coded b/y ion series and
#'   rotated text labels above major annotated peaks. Suitable for the
#'   dedicated Spectrum Viewer tab.
#'
#' @param matched_df data.frame. Output of \code{match_ions()}.
#' @param peptide_seq Character. Bare amino acid sequence for title.
#' @param precursor_mz Numeric. Observed precursor m/z.
#' @param precursor_charge Integer. Precursor charge state.
#' @param scan_rt Numeric. Retention time in minutes.
#' @param label_threshold Numeric. Minimum relative intensity (0-100) for a
#'   matched peak to receive a visible label (default 3).
#' @return A plotly figure object.
#' @export
plot_ms2_spectrum <- function(matched_df, peptide_seq = "",
                               precursor_mz     = NA_real_,
                               precursor_charge = NA_integer_,
                               scan_rt          = NA_real_,
                               label_threshold  = 3) {
  if (!requireNamespace("plotly", quietly = TRUE))
    stop("plotly is required.")

  df <- matched_df
  if (nrow(df) == 0)
    return(plotly::plot_ly() |> plotly::layout(title = "No peaks in spectrum"))

  # Classify series
  df$series <- "unmatched"
  df$series[df$matched & df$ion_series == "b" & !grepl("(-H2O|-NH3)$", df$ion_label)] <- "b"
  df$series[df$matched & df$ion_series == "y" & !grepl("(-H2O|-NH3)$", df$ion_label)] <- "y"
  df$series[df$matched & df$ion_series == "b" &  grepl("(-H2O|-NH3)$", df$ion_label)] <- "b_loss"
  df$series[df$matched & df$ion_series == "y" &  grepl("(-H2O|-NH3)$", df$ion_label)] <- "y_loss"

  col_map <- c(
    "b"         = "#1a73e8",
    "y"         = "#d93025",
    "b_loss"    = "#93c4f5",
    "y_loss"    = "#f5a3a3",
    "unmatched" = "#c8c8c8"
  )
  lbl_col_map <- c(
    "b" = "#1565c0", "y" = "#b71c1c",
    "b_loss" = "#1565c0", "y_loss" = "#b71c1c",
    "unmatched" = "#aaaaaa"
  )
  legend_names <- c(
    "b" = "b-ions", "y" = "y-ions",
    "b_loss" = "b-H₂O/NH₃",
    "y_loss" = "y-H₂O/NH₃",
    "unmatched" = "Unmatched"
  )

  # Title
  prec_str <- if (!is.na(precursor_mz)) {
    z_str <- if (!is.na(precursor_charge)) sprintf("+%d", as.integer(precursor_charge)) else "?"
    sprintf("%.4f (%s)", precursor_mz, z_str)
  } else NULL
  rt_str <- if (!is.na(scan_rt)) sprintf("RT %.2f min", scan_rt) else NULL
  pep_str <- if (nchar(trimws(peptide_seq)) > 0) peptide_seq else NULL
  title_str <- paste(c(pep_str, prec_str, rt_str), collapse = "  |  ")
  if (nchar(title_str) == 0) title_str <- "MS2 Spectrum"

  # Build traces (vertical sticks per series group)
  series_order <- c("unmatched", "b_loss", "y_loss", "b", "y")
  fig <- plotly::plot_ly()

  for (ser in series_order) {
    sdf <- df[df$series == ser, ]
    if (nrow(sdf) == 0) next

    # Vertical stick: (mz,0) → (mz, intensity) → NA per peak
    mz_v  <- as.vector(rbind(sdf$mz, sdf$mz, NA))
    int_v <- as.vector(rbind(rep(0, nrow(sdf)), sdf$rel_intensity, NA))

    lw <- if (ser == "unmatched") 1 else 1.6
    fig <- plotly::add_trace(fig,
      x = mz_v, y = int_v,
      type = "scatter", mode = "lines",
      line = list(color = col_map[ser], width = lw),
      name = legend_names[ser],
      showlegend = (ser != "unmatched"),
      legendgroup = ser,
      hoverinfo = "skip"
    )

    # Invisible hover points at peak tips
    hover_text <- sprintf(
      "<b>Ion:</b> %s<br><b>m/z:</b> %.4f<br><b>Rel. Int.:</b> %.1f%%",
      ifelse(is.na(sdf$ion_label), "—", sdf$ion_label),
      sdf$mz, sdf$rel_intensity
    )
    fig <- plotly::add_trace(fig,
      x = sdf$mz, y = sdf$rel_intensity,
      type = "scatter", mode = "markers",
      marker = list(size = 8, opacity = 0, color = col_map[ser]),
      name = legend_names[ser],
      showlegend = FALSE,
      legendgroup = ser,
      text = hover_text,
      hovertemplate = "%{text}<extra></extra>"
    )
  }

  # On-chart labels for major matched peaks
  ann_df <- df[df$matched & !is.na(df$ion_label) &
               df$rel_intensity >= label_threshold, ]
  ann_list <- if (nrow(ann_df) > 0) {
    lapply(seq_len(nrow(ann_df)), function(i) {
      list(
        x         = ann_df$mz[i],
        y         = ann_df$rel_intensity[i],
        text      = ann_df$ion_label[i],
        showarrow = FALSE,
        xanchor   = "center",
        yanchor   = "bottom",
        textangle = -90,
        font      = list(size = 9, color = lbl_col_map[ann_df$series[i]]),
        yshift    = 4
      )
    })
  } else list()

  fig |> plotly::layout(
    title  = list(text = title_str, font = list(size = 12), x = 0.02, xanchor = "left"),
    xaxis  = list(title = "m/z", showgrid = FALSE, zeroline = FALSE),
    yaxis  = list(
      title     = "Relative Intensity (%)",
      range     = c(0, 125),
      showgrid  = TRUE, gridcolor = "#f0f0f0", zeroline = FALSE
    ),
    legend = list(orientation = "h", x = 0.5, xanchor = "center",
                  y = -0.13, font = list(size = 10)),
    margin        = list(l = 65, r = 15, t = 50, b = 65),
    plot_bgcolor  = "#ffffff",
    paper_bgcolor = "#ffffff",
    annotations   = ann_list
  )
}
