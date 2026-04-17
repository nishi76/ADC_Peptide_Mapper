# uniqueness.R — Background proteome loading and uniqueness flagging
# Background databases are pre-built by build_background_db.R and
# stored as .rds files in data/. No internet access needed at runtime.

# ── Species registry ──────────────────────────────────────────────────────────
# ── Resolve data directory ────────────────────────────────────────────────────
# Priority: (1) data/ relative to app working dir, (2) /workspace/adc_peptide_mapper/data/
# This allows the app to work both from /workspace and when deployed elsewhere.
.resolve_data_dir <- function() {
  candidates <- c(
    file.path(getwd(), "data"),
    "/workspace/adc_peptide_mapper/data"
  )
  for (d in candidates) {
    if (dir.exists(d) &&
        any(file.exists(file.path(d, c("bg_human.rds","bg_monkey.rds","bg_rat.rds"))))) {
      return(d)
    }
  }
  # Fallback: relative data/ (will show "not found" gracefully in UI)
  file.path(getwd(), "data")
}
.DATA_DIR <- .resolve_data_dir()

SPECIES_REGISTRY <- list(
  human = list(
    label    = "Human (Homo sapiens)",
    rds_file = file.path(.DATA_DIR, "bg_human.rds")
  ),
  monkey = list(
    label    = "Cynomolgus Monkey (Macaca fascicularis)",
    rds_file = file.path(.DATA_DIR, "bg_monkey.rds")
  ),
  rat = list(
    label    = "Rat (Rattus norvegicus)",
    rds_file = file.path(.DATA_DIR, "bg_rat.rds")
  ),
  custom = list(
    label    = "Custom FASTA (user upload)",
    rds_file = NULL
  )
)

# ── Check which bundled databases are available ───────────────────────────────
get_available_species <- function() {
  avail <- list()
  for (key in names(SPECIES_REGISTRY)) {
    spec <- SPECIES_REGISTRY[[key]]
    if (key == "custom") {
      avail[[key]] <- spec$label
    } else if (!is.null(spec$rds_file) && file.exists(spec$rds_file)) {
      avail[[key]] <- spec$label
    }
  }
  avail
}

# ── Load a pre-built background set from RDS ─────────────────────────────────
# Returns the bg_object (list with $bg_sets, $n_proteins, $label, etc.)
load_background_rds <- function(species_key) {
  spec <- SPECIES_REGISTRY[[species_key]]
  if (is.null(spec) || is.null(spec$rds_file)) {
    stop("No RDS file defined for species: ", species_key)
  }
  if (!file.exists(spec$rds_file)) {
    stop("Background database not found: ", spec$rds_file,
         "\nRun build_background_db.R to generate it.")
  }
  message("Loading background: ", spec$rds_file)
  readRDS(spec$rds_file)
}

# ── Build background set from a custom uploaded FASTA ────────────────────────
# Used when user uploads their own background FASTA
# Returns same structure as load_background_rds()
build_background_from_fasta <- function(fasta_path,
                                         missed_cleavages = 0,
                                         min_len = 6, max_len = 30,
                                         progress_cb = NULL) {
  if (!file.exists(fasta_path)) stop("FASTA file not found: ", fasta_path)

  if (!is.null(progress_cb)) progress_cb("Reading custom FASTA...")
  fasta_text <- paste(readLines(fasta_path, warn = FALSE), collapse = "\n")

  if (!is.null(progress_cb)) progress_cb("Parsing sequences...")
  chains     <- parse_fasta(fasta_text)
  n_proteins <- length(chains)

  if (!is.null(progress_cb))
    progress_cb(paste0("Digesting ", n_proteins, " proteins (MC=", missed_cleavages, ")..."))

  all_seqs <- character(0)
  for (i in seq_along(chains)) {
    peptides <- trypsin_cleave(chains[[i]], missed_cleavages)
    peptides <- filter_peptides_by_length(peptides, min_len, max_len)
    if (length(peptides) > 0)
      all_seqs <- c(all_seqs, vapply(peptides, `[[`, character(1), "sequence"))
  }

  dt <- data.table::data.table(Sequence = unique(all_seqs))
  data.table::setkey(dt, Sequence)

  if (!is.null(progress_cb)) progress_cb("Done.")

  # Return same structure as pre-built RDS
  list(
    species     = "custom",
    label       = paste0("Custom FASTA (", n_proteins, " proteins)"),
    n_proteins  = n_proteins,
    build_date  = Sys.Date(),
    bg_sets     = setNames(list(dt), as.character(missed_cleavages))
  )
}

# ── Get the correct peptide set for a given missed cleavage level ─────────────
# bg_obj: result of load_background_rds() or build_background_from_fasta()
# mc: integer missed cleavages (0, 1, or 2)
# For custom uploads only one MC level is built; falls back to closest available
get_bg_dt <- function(bg_obj, mc) {
  mc_key <- as.character(mc)
  if (mc_key %in% names(bg_obj$bg_sets)) {
    return(bg_obj$bg_sets[[mc_key]])
  }
  # Fallback: use the first available set
  warning("MC=", mc, " not found in background set; using MC=",
          names(bg_obj$bg_sets)[1])
  bg_obj$bg_sets[[1]]
}

# ── Flag ADC peptide uniqueness ───────────────────────────────────────────────
# peptide_sequences: character vector of bare AA sequences
# bg_dt: data.table with key column "Sequence" (from get_bg_dt())
# Returns: logical vector, TRUE = unique to ADC (not in background)
flag_unique_peptides <- function(peptide_sequences, bg_dt) {
  !peptide_sequences %in% bg_dt$Sequence
}
