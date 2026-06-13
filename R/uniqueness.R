# =============================================================================
# uniqueness.R — ADC Peptide Mapper v0.8
# =============================================================================
# Functions for checking peptide uniqueness against background proteomes.
# Sourced by app.R. No library() calls — depends on data.table loaded by app.
#
# Public API
# ----------
#   flag_unique_peptides(sequences, bg_dt, il_equivalent)
#   flag_cospecies_unique(sequences, bg_list, il_equivalent)
#   get_bg_dt(bg_obj, missed_cleavages)
#   build_background_from_fasta(fasta_path, missed_cleavages, enzyme_id, progress_cb)
#
# Background object structure (.rds files):
#   list(
#     label        = character(1),
#     n_proteins   = integer(1),
#     build_date   = character(1),
#     bg_sets      = list(
#       "0" = data.table(Sequence = character),   # MC=0 peptides
#       "1" = data.table(Sequence = character),   # MC=0+1 peptides
#       "2" = data.table(Sequence = character)    # MC=0+1+2 peptides
#     )
#   )
#
# Note: build_background_db.R stores peptides_mc0/1/2 as character vectors.
#       get_bg_dt() handles both formats transparently.
# =============================================================================


# =============================================================================
# flag_unique_peptides()
# =============================================================================
# Flag peptides that are NOT present in the background proteome.
#
# Args:
#   sequences    : character vector — peptide sequences to check (uppercase)
#   bg_dt        : data.table with column "Sequence" — background peptide set
#                  (may be NULL or zero-row, in which case all are flagged unique)
#   il_equivalent: logical(1) — if TRUE (default), treat isoleucine (I) and
#                  leucine (L) as identical during comparison, matching the
#                  behaviour of database search engines that cannot distinguish
#                  I/L by MS/MS alone. Both query sequences and background
#                  sequences are normalised to "L" before the set-membership
#                  test. Flag defaults to TRUE; set FALSE only when an
#                  orthogonal method (e.g. ECD/ETD, Peptide Retention Time
#                  Prediction) is available to differentiate I from L.
#
# Returns:
#   logical vector of same length as sequences:
#     TRUE  = peptide is unique to the ADC (not in background)
#     FALSE = peptide found in background
# =============================================================================
flag_unique_peptides <- function(sequences, bg_dt, il_equivalent = TRUE) {
  if (is.null(bg_dt) || nrow(bg_dt) == 0L) {
    return(rep(TRUE, length(sequences)))
  }

  # Ensure uppercase for comparison
  seqs_upper <- toupper(sequences)
  bg_seqs    <- toupper(bg_dt$Sequence)

  # I/L equivalence: normalise both query and background to "L"
  if (isTRUE(il_equivalent)) {
    seqs_upper <- chartr("I", "L", seqs_upper)
    bg_seqs    <- chartr("I", "L", bg_seqs)
  }

  !seqs_upper %in% bg_seqs
}


# =============================================================================
# flag_cospecies_unique()
# =============================================================================
# Cross-species co-uniqueness check: test each peptide against multiple
# background proteomes simultaneously.  This is the correct approach for
# cross-species PK studies (e.g. human / cyno / rat in the same experiment),
# where a peptide is only analytically useful if it is unique in ALL tested
# species at once.
#
# Args:
#   sequences    : character vector — peptide sequences to check (uppercase)
#   bg_list      : named list of data.tables, each with column "Sequence".
#                  Names become column names in the result (e.g. "human",
#                  "monkey", "rat"). May include NULL elements (species not
#                  loaded); those columns are filled with NA.
#   il_equivalent: logical(1) — treat I/L as identical (default TRUE).
#                  Forwarded to flag_unique_peptides() for each species.
#
# Returns:
#   data.frame with one column per species (logical: TRUE = unique) PLUS:
#     UniqueAllSpecies : logical — TRUE iff unique in every non-NULL species
#
# Example:
#   bg_list <- list(human = get_bg_dt(BUNDLED_BG$human, 0),
#                   monkey = get_bg_dt(BUNDLED_BG$monkey, 0))
#   result  <- flag_cospecies_unique(pep_dt$Sequence, bg_list)
#   pep_dt  <- cbind(pep_dt, result)
# =============================================================================
flag_cospecies_unique <- function(sequences, bg_list, il_equivalent = TRUE) {
  n <- length(sequences)
  sp_names <- names(bg_list)

  result <- as.data.frame(
    setNames(
      lapply(sp_names, function(sp) {
        bg <- bg_list[[sp]]
        if (is.null(bg) || (is.data.frame(bg) && nrow(bg) == 0L)) {
          return(rep(NA, n))
        }
        flag_unique_peptides(sequences, bg, il_equivalent = il_equivalent)
      }),
      sp_names
    ),
    stringsAsFactors = FALSE
  )

  # A peptide is co-unique iff it is unique in every species with a non-NA result
  result$UniqueAllSpecies <- apply(result, 1, function(row) {
    non_na <- row[!is.na(row)]
    if (length(non_na) == 0L) return(NA)
    all(as.logical(non_na))
  })

  result
}


# =============================================================================
# get_bg_dt()
# =============================================================================
# Extract a data.table of background peptide sequences from a background
# object at the requested missed-cleavage level.
#
# Handles two .rds formats:
#   Format A (build_background_db.R): bg_obj$peptides_mc0, $peptides_mc1, $peptides_mc2
#   Format B (bg_sets):               bg_obj$bg_sets[["0"]], [["1"]], [["2"]]
#
# Args:
#   bg_obj          : list — background object loaded from .rds
#   missed_cleavages: integer(1) — 0, 1, or 2
#
# Returns:
#   data.table with column "Sequence" (character)
# =============================================================================
get_bg_dt <- function(bg_obj, missed_cleavages) {
  mc <- as.integer(missed_cleavages)
  mc <- max(0L, min(2L, mc))   # clamp to 0-2

  # ── Format A: peptides_mc0 / peptides_mc1 / peptides_mc2 ─────────────────
  mc_key_a <- paste0("peptides_mc", mc)
  if (!is.null(bg_obj[[mc_key_a]])) {
    peps <- bg_obj[[mc_key_a]]
    if (is.character(peps)) {
      return(data.table::data.table(Sequence = toupper(peps)))
    }
    if (is.data.frame(peps) || data.table::is.data.table(peps)) {
      if ("Sequence" %in% names(peps)) return(peps)
    }
  }

  # ── Format B: bg_sets[["0"]] / [["1"]] / [["2"]] ─────────────────────────
  if (!is.null(bg_obj$bg_sets)) {
    mc_key_b <- as.character(mc)
    if (!is.null(bg_obj$bg_sets[[mc_key_b]])) {
      dt <- bg_obj$bg_sets[[mc_key_b]]
      if (is.character(dt)) {
        return(data.table::data.table(Sequence = toupper(dt)))
      }
      if (data.table::is.data.table(dt) || is.data.frame(dt)) {
        if ("Sequence" %in% names(dt)) return(data.table::as.data.table(dt))
      }
    }
  }

  # ── Fallback: return empty table ──────────────────────────────────────────
  message(sprintf(
    "get_bg_dt: could not find peptides for MC=%d in background object. ",
    mc, "Returning empty set (all peptides will be flagged unique)."
  ))
  data.table::data.table(Sequence = character(0))
}


# =============================================================================
# build_background_from_fasta()
# =============================================================================
# Build a background proteome object from a user-uploaded FASTA file.
# Used for the "Custom FASTA" background option in Tab 1.
#
# The background is digested with the SAME enzyme the user chose for their
# ADC sequence, so uniqueness comparisons are enzyme-consistent.
#
# Args:
#   fasta_path       : character(1) — path to FASTA file
#   missed_cleavages : integer(1)   — 0, 1, or 2 (digest MC level)
#   enzyme_id        : character(1) — enzyme identifier passed to
#                      enzyme_cleave() (default "trypsin"). Must match
#                      one of the keys in ENZYME_LABELS.
#   progress_cb      : function(msg) or NULL — optional progress callback
#
# Returns:
#   list(label, n_proteins, build_date, enzyme_id, bg_sets)
#   bg_sets is a named list: "0", "1", "2" → data.table(Sequence)
# =============================================================================
build_background_from_fasta <- function(fasta_path,
                                         missed_cleavages = 0L,
                                         enzyme_id        = "trypsin",
                                         progress_cb      = NULL) {

  .progress <- function(msg) {
    if (!is.null(progress_cb)) progress_cb(msg)
    message(msg)
  }

  mc         <- as.integer(missed_cleavages)
  enzyme_id  <- tolower(trimws(enzyme_id))

  # ── Parse FASTA ────────────────────────────────────────────────────────────
  .progress("Parsing FASTA...")
  fasta_text <- paste(readLines(fasta_path, warn = FALSE), collapse = "\n")
  chains     <- parse_fasta(fasta_text)

  if (length(chains) == 0L) {
    stop("No valid sequences found in the uploaded FASTA file.")
  }

  n_proteins <- length(chains)
  .progress(sprintf("Parsed %d sequences.", n_proteins))

  # ── Digest at requested enzyme + MC level ─────────────────────────────────
  .progress(sprintf("Digesting with %s at MC=%d...", enzyme_id, mc))

  .digest_one_mc <- function(mc_level) {
    all_peps <- character(0)
    for (seq in chains) {
      pep_dt <- tryCatch(
        enzyme_cleave(seq, enzyme_id = enzyme_id, missed_cleavages = mc_level),
        error = function(e) data.table::data.table(Sequence = character(0))
      )
      # enzyme_cleave() returns a data.table — extract the Sequence column
      peps <- pep_dt$Sequence
      peps <- peps[nzchar(peps) & nchar(peps) >= 6L & nchar(peps) <= 30L]
      all_peps <- c(all_peps, peps)
    }
    data.table::data.table(Sequence = unique(toupper(all_peps)))
  }

  # Build bg_sets for MC 0, 1, 2 (always build all three for consistency)
  bg_sets <- list()
  for (mc_level in 0L:2L) {
    .progress(sprintf("  Digesting MC=%d...", mc_level))
    bg_sets[[as.character(mc_level)]] <- .digest_one_mc(mc_level)
  }

  list(
    label      = sprintf("Custom FASTA (%d proteins, %s)", n_proteins, enzyme_id),
    n_proteins = n_proteins,
    build_date = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    enzyme_id  = enzyme_id,
    bg_sets    = bg_sets
  )
}
