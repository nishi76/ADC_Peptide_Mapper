# ============================================================
#  build_background_db.R
#  Run ONCE to download and pre-digest background proteomes.
#  Saves .rds files to data/ — bundled with the app.
#  Usage: Rscript build_background_db.R
#         or source("build_background_db.R") from RStudio
# ============================================================

library(data.table)
library(httr2)

source("R/digest.R")

# ── Species definitions ───────────────────────────────────────────────────────
# Output directory: write to /workspace/adc_peptide_mapper/data/ (local filesystem)
# The app reads from this location via the auto-resolver in R/uniqueness.R
DATA_DIR <- "/workspace/adc_peptide_mapper/data"
if (!dir.exists(DATA_DIR)) dir.create(DATA_DIR, recursive = TRUE)

SPECIES <- list(
  human = list(
    label       = "Human (Homo sapiens)",
    organism_id = "9606",
    rds_file    = file.path(DATA_DIR, "bg_human.rds")
  ),
  monkey = list(
    label       = "Cynomolgus Monkey (Macaca fascicularis)",
    organism_id = "9541",
    rds_file    = file.path(DATA_DIR, "bg_monkey.rds")
  ),
  rat = list(
    label       = "Rat (Rattus norvegicus)",
    organism_id = "10116",
    rds_file    = file.path(DATA_DIR, "bg_rat.rds")
  )
)

# ── UniProt download ──────────────────────────────────────────────────────────
download_uniprot_fasta <- function(organism_id, label) {
  url <- paste0(
    "https://rest.uniprot.org/uniprotkb/stream",
    "?query=reviewed%3Atrue%20AND%20organism_id%3A", organism_id,
    "&format=fasta&compressed=false"
  )
  cat(sprintf("[%s] Downloading UniProt Swiss-Prot FASTA (organism %s)...\n",
              label, organism_id))
  resp <- httr2::request(url) |>
    httr2::req_timeout(600) |>
    httr2::req_retry(max_tries = 3, backoff = ~ 10) |>
    httr2::req_perform()
  httr2::resp_body_string(resp)
}

# ── Build and save background set ────────────────────────────────────────────
build_and_save <- function(key, spec, missed_cleavages_list = c(0, 1, 2),
                            min_len = 6, max_len = 30) {
  rds_path <- spec$rds_file

  if (file.exists(rds_path) && file.size(rds_path) > 1000) {
    cat(sprintf("[%s] Already exists: %s — skipping.\n", spec$label, rds_path))
    return(invisible(NULL))
  }

  # Download FASTA
  fasta_text <- download_uniprot_fasta(spec$organism_id, spec$label)

  # Count proteins
  n_proteins <- length(gregexpr(">", fasta_text, fixed = TRUE)[[1]])
  cat(sprintf("[%s] %d proteins downloaded.\n", spec$label, n_proteins))

  # Parse chains
  cat(sprintf("[%s] Parsing sequences...\n", spec$label))
  chains <- parse_fasta(fasta_text)

  # Build peptide sets for each missed cleavage level
  bg_sets <- list()

  for (mc in missed_cleavages_list) {
    cat(sprintf("[%s] Digesting (missed cleavages = %d)...\n", spec$label, mc))

    all_seqs <- character(0)
    n_chains  <- length(chains)

    for (i in seq_along(chains)) {
      seq_str  <- chains[[i]]
      peptides <- trypsin_cleave(seq_str, mc)
      peptides <- filter_peptides_by_length(peptides, min_len, max_len)
      if (length(peptides) > 0) {
        all_seqs <- c(all_seqs, vapply(peptides, `[[`, character(1), "sequence"))
      }
      if (i %% 1000 == 0) {
        cat(sprintf("  ... %d / %d proteins processed\n", i, n_chains))
      }
    }

    # Build keyed data.table for O(1) lookup
    dt <- data.table(Sequence = unique(all_seqs))
    setkey(dt, Sequence)

    bg_sets[[as.character(mc)]] <- dt
    cat(sprintf("[%s] MC=%d: %s unique peptides\n",
                spec$label, mc, format(nrow(dt), big.mark = ",")))
  }

  # Save as RDS
  result <- list(
    species      = key,
    label        = spec$label,
    organism_id  = spec$organism_id,
    n_proteins   = n_proteins,
    build_date   = Sys.Date(),
    bg_sets      = bg_sets   # named list: "0", "1", "2"
  )

  saveRDS(result, rds_path, compress = "xz")
  cat(sprintf("[%s] Saved to %s (%.1f MB)\n",
              spec$label, rds_path,
              file.size(rds_path) / 1e6))
}

# ── Run for all species ───────────────────────────────────────────────────────
cat("=== ADC Peptide Mapper — Background Database Builder ===\n")
cat(sprintf("Started: %s\n\n", Sys.time()))

for (key in names(SPECIES)) {
  tryCatch(
    build_and_save(key, SPECIES[[key]]),
    error = function(e) {
      cat(sprintf("ERROR building %s: %s\n", key, conditionMessage(e)))
    }
  )
}

cat(sprintf("\nDone: %s\n", Sys.time()))
cat("Background databases saved to data/\n")
cat("Files:\n")
for (f in list.files("data", pattern = "\\.rds$", full.names = TRUE)) {
  cat(sprintf("  %s  (%.1f MB)\n", f, file.size(f) / 1e6))
}
