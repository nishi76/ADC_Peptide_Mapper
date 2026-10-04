# =============================================================================
# digest.R  -  ADC Peptide Mapper v0.8
# =============================================================================
# Sourced by app.R.  No library() calls here - all packages (data.table, etc.)
# are loaded in app.R before this file is sourced.
#
# Contents
# --------
#   1. Mass constants
#   2. ENZYME_LABELS  - human-readable names for the 9 supported enzymes
#   3. parse_fasta()              - FASTA text -> named character vector
#   4. calc_peptide_mass()        - monoisotopic neutral mass of a peptide
#   5. filter_peptides_by_length()- keep peptides within [min_len, max_len]
#   6. enzyme_cleave()            - generic in-silico digest (9 enzymes)
#   7. trypsin_cleave()           - backward-compat alias -> enzyme_cleave()
#   8. digest_fasta_enzyme()      - digest a full FASTA with one or two enzymes
#   9. digest_fasta()             - backward-compat wrapper -> digest_fasta_enzyme()
# =============================================================================


# -----------------------------------------------------------------------------
# 1.  Mass constants
# -----------------------------------------------------------------------------

#' Monoisotopic residue masses for the 20 standard amino acids
#'
#' @description Named numeric vector of monoisotopic residue masses (Da) for
#'   all 20 standard amino acids. Keys are single-letter codes (uppercase).
#'   Used throughout the mass calculation functions in this package.
#'
#' @format A named numeric vector of length 20.
#' @export
#' @examples
#' AA_MONO_MASS["K"]  # lysine residue mass
#' sum(AA_MONO_MASS[c("A","K","L")]) + WATER_MASS  # AKL peptide neutral mass
AA_MONO_MASS <- c(
  A =  71.03711,  R = 156.10111,  N = 114.04293,  D = 115.02694,
  C = 103.00919,  E = 129.04259,  Q = 128.05858,  G =  57.02146,
  H = 137.05891,  I = 113.08406,  L = 113.08406,  K = 128.09496,
  M = 131.04049,  F = 147.06841,  P =  97.05276,  S =  87.03203,
  T = 101.04768,  W = 186.07931,  Y = 163.06333,  V =  99.06841
)

#' Mass of water (Da)
#' @description Monoisotopic mass of water (H2O), used to compute peptide
#'   neutral masses and y-ion series.
#' @format A numeric scalar.
#' @export
WATER_MASS  <- 18.01056

#' Mass of a proton (Da)
#' @description Monoisotopic mass of a proton (H+), used to convert neutral
#'   masses to m/z values for protonated ions.
#' @format A numeric scalar.
#' @export
PROTON_MASS <-  1.007276


# -----------------------------------------------------------------------------
# 2.  ENZYME_LABELS
# -----------------------------------------------------------------------------

#' Human-readable labels for all supported proteolytic enzymes
#'
#' @description Named character vector mapping internal enzyme identifiers to
#'   human-readable display names. Used to populate the enzyme selector UI
#'   and to label the \code{Enzyme} column in digest output tables.
#'
#' @format A named character vector with 11 elements.
#' @export
#' @examples
#' ENZYME_LABELS["trypsin"]
#' names(ENZYME_LABELS)
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


# -----------------------------------------------------------------------------
# 3.  parse_fasta()
# -----------------------------------------------------------------------------

#' Parse a FASTA text string into a named character vector
#'
#' @description Splits the contents of a FASTA file (provided as a single
#'   string) into a named character vector where names are the sequence
#'   headers and values are the uppercase amino acid sequences.
#'
#' @param fasta_text character(1). Full contents of a FASTA file as a single
#'   string. Accepts both Unix (\code{\\n}) and Windows (\code{\\r\\n}) line endings.
#'
#' @return A named character vector. Names are the header lines (without the
#'   leading \code{>}), values are uppercase sequences with whitespace removed.
#'   Returns a zero-length named character vector if the input is empty or
#'   contains no valid headers.
#'
#' @export
#' @examples
#' fasta <- ">Heavy chain\nEVQLVESGGGLVQPGGSLR\n>Light chain\nDIQMTQSPS"
#' seqs <- parse_fasta(fasta)
#' names(seqs)
#' nchar(seqs)
parse_fasta <- function(fasta_text) {
  if (is.null(fasta_text) || length(fasta_text) == 0L ||
      !nzchar(trimws(fasta_text))) {
    return(setNames(character(0L), character(0L)))
  }
  lines    <- strsplit(gsub("\r\n", "\n", fasta_text), "\n", fixed = TRUE)[[1L]]
  is_header <- startsWith(lines, ">")
  if (!any(is_header)) return(setNames(character(0L), character(0L)))
  group_id    <- cumsum(is_header)
  header_lines <- lines[is_header]
  chain_names  <- trimws(sub("^>", "", header_lines))
  sequences <- vapply(seq_along(chain_names), function(i) {
    seq_lines <- lines[!is_header & group_id == i]
    paste0(toupper(gsub("[[:space:]]", "", seq_lines)), collapse = "")
  }, character(1L))
  setNames(sequences, chain_names)
}


# -----------------------------------------------------------------------------
# 4.  calc_peptide_mass()
# -----------------------------------------------------------------------------

#' Compute the monoisotopic neutral mass of a peptide
#'
#' @description Calculates the neutral monoisotopic mass of a peptide from its
#'   single-letter amino acid sequence using the residue masses in
#'   \code{\link{AA_MONO_MASS}} and adding one water molecule for the
#'   N- and C-terminal groups.
#'
#' @param sequence character(1). Single-letter amino acid sequence (uppercase).
#'
#' @return numeric(1). Monoisotopic neutral mass in Da. Returns \code{NA_real_}
#'   for \code{NULL} or zero-length input. Unknown amino acid codes contribute
#'   a mass of 0.
#'
#' @export
#' @examples
#' calc_peptide_mass("PEPTIDE")   # ~799.36 Da
#' calc_peptide_mass("EVQLVESGGG")
calc_peptide_mass <- function(sequence) {
  if (is.null(sequence) || !nzchar(sequence)) return(NA_real_)
  residues <- strsplit(sequence, "", fixed = TRUE)[[1L]]
  masses   <- AA_MONO_MASS[residues]
  masses[is.na(masses)] <- 0
  sum(masses) + WATER_MASS
}


# -----------------------------------------------------------------------------
# 5.  filter_peptides_by_length()
# -----------------------------------------------------------------------------

#' Filter a peptide data.table by sequence length
#'
#' @description Keeps only rows whose \code{Length} column falls within
#'   \code{[min_len, max_len]}. Passes through \code{NULL} or empty tables
#'   unchanged.
#'
#' @param dt data.table. Must contain an integer \code{Length} column.
#' @param min_len integer(1). Minimum peptide length (inclusive). Default 6.
#' @param max_len integer(1). Maximum peptide length (inclusive). Default 30.
#'
#' @return data.table. Subset of \code{dt} containing only rows where
#'   \code{min_len <= Length <= max_len}.
#'
#' @export
#' @examples
#' \dontrun{
#'   pep_dt <- digest_fasta_enzyme(fasta, min_len = 1, max_len = 50)
#'   filtered <- filter_peptides_by_length(pep_dt, min_len = 7, max_len = 25)
#' }
filter_peptides_by_length <- function(dt, min_len = 6L, max_len = 30L) {
  if (is.null(dt) || nrow(dt) == 0L) return(dt)
  dt[Length >= min_len & Length <= max_len]
}


# -----------------------------------------------------------------------------
# 6.  enzyme_cleave()
# -----------------------------------------------------------------------------

#' In-silico proteolytic digest of a single protein sequence
#'
#' @description Performs a theoretical (in-silico) proteolytic digest of a
#'   single amino acid sequence using one of 11 supported enzymes. Returns a
#'   data.table of peptide fragments, optionally including missed cleavages.
#'
#' @param sequence character(1). Single-letter amino acid sequence (uppercase).
#' @param enzyme_id character(1). Enzyme identifier. One of:
#'   \code{"trypsin"}, \code{"trypsin_p"}, \code{"lysc"}, \code{"lysc_p"},
#'   \code{"lysn"}, \code{"aspn"}, \code{"gluc"}, \code{"argc"},
#'   \code{"chymotrypsin"}, \code{"papain"}, \code{"elastase"}.
#'   See \code{\link{ENZYME_LABELS}} for human-readable names.
#' @param missed_cleavages integer(1). Number of missed cleavages to include
#'   (0 = fully digested, 1 = one internal cleavage site missed, etc.).
#'   Default 0.
#'
#' @return data.table with columns:
#'   \describe{
#'     \item{Sequence}{character. Peptide amino acid sequence.}
#'     \item{Start}{integer. 1-based start position in the parent protein.}
#'     \item{End}{integer. 1-based end position.}
#'     \item{Length}{integer. Number of residues.}
#'     \item{MC}{integer. Number of missed cleavages for this peptide.}
#'   }
#'
#' @export
#' @examples
#' enzyme_cleave("MEPEPTIDEKRPEPTIDE", enzyme_id = "trypsin", missed_cleavages = 0)
#' enzyme_cleave("MEPEPTIDEKRPEPTIDE", enzyme_id = "trypsin", missed_cleavages = 1)
enzyme_cleave <- function(sequence, enzyme_id = "trypsin", missed_cleavages = 0L) {

  empty_dt <- data.table(
    Sequence = character(0L), Start = integer(0L), End = integer(0L),
    Length   = integer(0L),   MC    = integer(0L)
  )

  if (is.null(sequence) || length(sequence) == 0L || !nzchar(sequence)) return(empty_dt)

  enzyme_id        <- tolower(trimws(enzyme_id))
  missed_cleavages <- max(0L, as.integer(missed_cleavages))
  n                <- nchar(sequence)

  if (n == 1L) {
    return(data.table(Sequence=sequence, Start=1L, End=1L, Length=1L, MC=0L))
  }

  residues   <- strsplit(sequence, "", fixed = TRUE)[[1L]]
  cut_after  <- logical(n)
  cut_before <- logical(n)

  if (enzyme_id == "trypsin") {
    for (i in seq_len(n - 1L)) {
      if (residues[i] %in% c("K","R") && residues[i+1L] != "P") cut_after[i] <- TRUE
    }
    if (residues[n] %in% c("K","R")) cut_after[n] <- TRUE

  } else if (enzyme_id == "trypsin_p") {
    # Trypsin/P: cleaves after K/R regardless of following P
    for (i in seq_len(n)) {
      if (residues[i] %in% c("K","R")) cut_after[i] <- TRUE
    }

  } else if (enzyme_id == "lysc_p") {
    # Lys-C/P: cleaves after K regardless of following P
    for (i in seq_len(n)) {
      if (residues[i] == "K") cut_after[i] <- TRUE
    }

  } else if (enzyme_id %in% c("lysn", "lys_n")) {
    # Lys-N: cleaves before K (lys_n kept as internal alias)
    for (i in seq(2L, n)) {
      if (residues[i] == "K") cut_before[i] <- TRUE
    }

  } else if (enzyme_id == "chymotrypsin") {
    for (i in seq_len(n - 1L)) {
      if (residues[i] %in% c("F","Y","W") && residues[i+1L] != "P") cut_after[i] <- TRUE
    }
    if (residues[n] %in% c("F","Y","W")) cut_after[n] <- TRUE

  } else if (enzyme_id == "argc") {
    for (i in seq_len(n - 1L)) {
      if (residues[i] == "R" && residues[i+1L] != "P") cut_after[i] <- TRUE
    }
    if (residues[n] == "R") cut_after[n] <- TRUE

  } else if (enzyme_id == "lysc") {
    for (i in seq_len(n - 1L)) {
      if (residues[i] == "K" && residues[i+1L] != "P") cut_after[i] <- TRUE
    }
    if (residues[n] == "K") cut_after[n] <- TRUE

  } else if (enzyme_id == "papain") {
    for (i in seq_len(n)) {
      if (residues[i] %in% c("K","R","Q")) cut_after[i] <- TRUE
    }

  } else if (enzyme_id == "elastase") {
    for (i in seq_len(n)) {
      if (residues[i] %in% c("A","V","S","G","T")) cut_after[i] <- TRUE
    }

  } else if (enzyme_id == "gluc") {
    for (i in seq_len(n)) {
      if (residues[i] %in% c("E","D")) cut_after[i] <- TRUE
    }

  } else if (enzyme_id == "aspn") {
    for (i in seq(2L, n)) {
      if (residues[i] == "D") cut_before[i] <- TRUE
    }

  } else {
    stop(sprintf("enzyme_cleave(): unknown enzyme_id '%s'. Valid: %s",
                 enzyme_id, paste(names(ENZYME_LABELS), collapse=", ")))
  }

  new_pep_starts <- 1L
  after_positions <- which(cut_after)
  after_positions <- after_positions[after_positions < n]
  if (length(after_positions) > 0L) new_pep_starts <- c(new_pep_starts, after_positions + 1L)
  before_positions <- which(cut_before)
  if (length(before_positions) > 0L) new_pep_starts <- c(new_pep_starts, before_positions)
  new_pep_starts <- sort(unique(new_pep_starts))

  num_peps   <- length(new_pep_starts)
  pep_ends   <- c(new_pep_starts[-1L] - 1L, n)
  pep_starts <- new_pep_starts
  pep_seqs   <- mapply(function(s,e) substr(sequence,s,e), pep_starts, pep_ends, USE.NAMES=FALSE)
  pep_lens   <- pep_ends - pep_starts + 1L

  valid      <- pep_lens > 0L
  pep_seqs   <- pep_seqs[valid]; pep_starts <- pep_starts[valid]
  pep_ends   <- pep_ends[valid]; pep_lens   <- pep_lens[valid]
  num_peps   <- sum(valid)

  base_dt <- data.table(Sequence=pep_seqs, Start=pep_starts, End=pep_ends,
                        Length=pep_lens, MC=0L)

  if (missed_cleavages > 0L && num_peps > 1L) {
    mc_rows <- vector("list", missed_cleavages * (num_peps - 1L))
    idx <- 0L
    for (mc in seq_len(missed_cleavages)) {
      for (i in seq_len(num_peps - mc)) {
        j      <- i + mc
        s      <- pep_starts[i]; e <- pep_ends[j]
        idx    <- idx + 1L
        mc_rows[[idx]] <- data.table(
          Sequence = substr(sequence, s, e), Start = s, End = e,
          Length   = e - s + 1L,            MC    = mc
        )
      }
    }
    mc_rows <- mc_rows[seq_len(idx)]
    if (length(mc_rows) > 0L)
      base_dt <- rbindlist(c(list(base_dt), mc_rows), use.names=TRUE)
  }

  setorder(base_dt, Start, MC)
  base_dt
}


# -----------------------------------------------------------------------------
# 7.  trypsin_cleave()  - backward-compatibility alias
# -----------------------------------------------------------------------------

#' Trypsin digest (backward-compatibility alias)
#'
#' @description Convenience wrapper that calls \code{\link{enzyme_cleave}} with
#'   \code{enzyme_id = "trypsin"}. Kept for backward compatibility.
#'
#' @param sequence character(1). Single-letter amino acid sequence.
#' @param missed_cleavages integer(1). Number of missed cleavages. Default 0.
#'
#' @return data.table. Same schema as \code{\link{enzyme_cleave}}.
#'
#' @export
#' @seealso \code{\link{enzyme_cleave}}, \code{\link{digest_fasta_enzyme}}
#' @examples
#' trypsin_cleave("MEPEPTIDEKRPEPTIDE")
trypsin_cleave <- function(sequence, missed_cleavages = 0L) {
  enzyme_cleave(sequence, enzyme_id="trypsin", missed_cleavages=missed_cleavages)
}


# -----------------------------------------------------------------------------
# 8.  digest_fasta_enzyme()
# -----------------------------------------------------------------------------

#' Digest all sequences in a FASTA string with one or two enzymes
#'
#' @description Parses a multi-sequence FASTA text and performs in-silico
#'   proteolytic digestion of every chain. Supports single-enzyme digestion and
#'   sequential two-enzyme digestion (pass 1 with \code{enzyme_id}, pass 2 with
#'   \code{enzyme_id2}).
#'
#' @param fasta_text character(1). Full FASTA file contents as a single string.
#' @param enzyme_id character(1). Primary enzyme identifier (see
#'   \code{\link{ENZYME_LABELS}}). Default \code{"trypsin"}.
#' @param enzyme_id2 character(1) or NULL. Optional second enzyme for
#'   sequential two-enzyme digestion. \code{NULL} (default) uses only
#'   \code{enzyme_id}.
#' @param missed_cleavages integer(1). Missed cleavages allowed. Default 0.
#' @param min_len integer(1). Minimum peptide length to retain. Default 6.
#' @param max_len integer(1). Maximum peptide length to retain. Default 30.
#'
#' @return data.table with columns: \code{Chain}, \code{Sequence},
#'   \code{Start}, \code{End}, \code{Length}, \code{Mass}, \code{Enzyme},
#'   and optionally \code{MC}. Returns an empty table if no sequences or no
#'   peptides pass length filters.
#'
#' @export
#' @examples
#' fasta <- paste0(">Heavy\nEVQLVESGGGLVQPGGSLRLSCAAS",
#'                 "\n>Light\nDIQMTQSPSSLSASVGDRVTITC")
#' digest_fasta_enzyme(fasta, enzyme_id = "trypsin", missed_cleavages = 1)
digest_fasta_enzyme <- function(fasta_text,
                                enzyme_id        = "trypsin",
                                enzyme_id2       = NULL,
                                missed_cleavages = 0L,
                                min_len          = 6L,
                                max_len          = 30L) {

  empty_result <- data.table(
    Chain=character(0L), Sequence=character(0L), Start=integer(0L),
    End=integer(0L), Length=integer(0L), Mass=numeric(0L), Enzyme=character(0L)
  )

  seqs <- parse_fasta(fasta_text)
  if (length(seqs) == 0L) return(empty_result)

  enzyme_id        <- tolower(trimws(enzyme_id))
  missed_cleavages <- max(0L, as.integer(missed_cleavages))
  min_len          <- as.integer(min_len)
  max_len          <- as.integer(max_len)

  use_two_enzyme <- !is.null(enzyme_id2) && nzchar(trimws(enzyme_id2))
  if (use_two_enzyme) {
    enzyme_id2   <- tolower(trimws(enzyme_id2))
    enzyme_label <- paste0(enzyme_id, "+", enzyme_id2)
  } else {
    enzyme_label <- enzyme_id
  }

  chain_results <- vector("list", length(seqs))

  for (ci in seq_along(seqs)) {
    chain_name <- names(seqs)[ci]
    chain_seq  <- seqs[[ci]]
    if (!nzchar(chain_seq)) next

    if (!use_two_enzyme) {
      pep_dt <- enzyme_cleave(chain_seq, enzyme_id=enzyme_id,
                              missed_cleavages=missed_cleavages)
    } else {
      pass1_dt <- enzyme_cleave(chain_seq, enzyme_id=enzyme_id, missed_cleavages=0L)
      if (nrow(pass1_dt) == 0L) next
      pass2_list <- vector("list", nrow(pass1_dt))
      for (pi in seq_len(nrow(pass1_dt))) {
        sub_seq   <- pass1_dt$Sequence[pi]
        sub_start <- pass1_dt$Start[pi]
        sub_dt    <- enzyme_cleave(sub_seq, enzyme_id=enzyme_id2,
                                   missed_cleavages=missed_cleavages)
        if (nrow(sub_dt) == 0L) next
        sub_dt[, Start := Start + sub_start - 1L]
        sub_dt[, End   := End   + sub_start - 1L]
        pass2_list[[pi]] <- sub_dt
      }
      pass2_list <- Filter(Negate(is.null), pass2_list)
      if (length(pass2_list) == 0L) next
      pep_dt <- rbindlist(pass2_list, use.names=TRUE)
      # Deduplicate on (Sequence, Start, End) to preserve positional information.
      # Two identical sequences at different positions are kept as separate rows.
      pep_dt <- unique(pep_dt, by=c("Sequence","Start","End"))
    }

    pep_dt <- filter_peptides_by_length(pep_dt, min_len=min_len, max_len=max_len)
    if (nrow(pep_dt) == 0L) next

    pep_dt[, Mass   := vapply(Sequence, calc_peptide_mass, numeric(1L))]
    pep_dt[, Chain  := chain_name]
    pep_dt[, Enzyme := enzyme_label]
    chain_results[[ci]] <- pep_dt
  }

  chain_results <- Filter(Negate(is.null), chain_results)
  if (length(chain_results) == 0L) return(empty_result)

  result    <- rbindlist(chain_results, use.names=TRUE, fill=TRUE)
  core_cols <- c("Chain","Sequence","Start","End","Length","Mass","Enzyme")
  extra_cols <- setdiff(names(result), core_cols)
  setcolorder(result, c(core_cols, extra_cols))
  if ("MC" %in% names(result)) setorder(result, Chain, Start, MC) else setorder(result, Chain, Start)
  result
}


# -----------------------------------------------------------------------------
# 9.  digest_fasta()  - backward-compatibility wrapper
# -----------------------------------------------------------------------------

#' Tryptic digest of a FASTA string (backward-compatibility wrapper)
#'
#' @description Calls \code{\link{digest_fasta_enzyme}} with
#'   \code{enzyme_id = "trypsin"}. Preserved for backward compatibility.
#'
#' @param fasta_text character(1). FASTA file contents as a single string.
#' @param missed_cleavages integer(1). Default 0.
#' @param min_len integer(1). Minimum peptide length. Default 6.
#' @param max_len integer(1). Maximum peptide length. Default 30.
#'
#' @return data.table. Same schema as \code{\link{digest_fasta_enzyme}}.
#'
#' @export
#' @seealso \code{\link{digest_fasta_enzyme}}
digest_fasta <- function(fasta_text, missed_cleavages=0L, min_len=6L, max_len=30L) {
  digest_fasta_enzyme(fasta_text, enzyme_id="trypsin", enzyme_id2=NULL,
                      missed_cleavages=missed_cleavages, min_len=min_len, max_len=max_len)
}

# =============================================================================
# End of digest.R
# =============================================================================
