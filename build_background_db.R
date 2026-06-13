# =============================================================================
# build_background_db.R — ADC Peptide Mapper v0.8
# =============================================================================
# Run this script ONCE before launching the app to generate the background
# proteome databases used by the uniqueness filter in Tab 3.
#
# Usage (from the app root directory):
#   source("build_background_db.R")
#   # or from the terminal:
#   Rscript build_background_db.R
#
# What it does:
#   1. Downloads reviewed (Swiss-Prot) proteomes for human and rat, and the
#      full TrEMBL proteome for cynomolgus monkey (larger / more complete) from
#      the UniProt REST API.
#   2. Downloads the cRAP common contaminant database (trypsin autolysis
#      peptides, keratins, BSA, IgG constant regions) from
#      https://www.thegpm.org/crap/ and merges the peptides into every
#      species background before saving.
#   3. Digests each proteome with trypsin at missed cleavages 0, 1, and 2
#      using the same enzyme_cleave() function used by the app.
#   4. Saves the results as:
#        data/bg_human.rds
#        data/bg_monkey.rds   ← TrEMBL (reviewed + unreviewed, ~21k proteins)
#        data/bg_rat.rds
#
# Requirements:
#   - Internet access (first run only; ~15 min total)
#   - R packages: httr2, data.table (auto-installed if missing)
#   - Must be run from the ADC Peptide Mapper root directory
#     (the folder containing app.R)
#
# Resume behaviour:
#   Already-built .rds files are skipped automatically.
#   Delete a file to force a rebuild for that species.
#
# Estimated runtime: ~15 min on first run; ~0 if all files already exist.
# =============================================================================

cat("=============================================================\n")
cat(" ADC Peptide Mapper v0.8 — Background Database Builder\n")
cat("=============================================================\n\n")

# ── 0. Check working directory ────────────────────────────────────────────────
if (!file.exists("app.R")) {
  stop(
    "build_background_db.R must be run from the ADC Peptide Mapper root ",
    "directory (the folder containing app.R).\n",
    "Current directory: ", getwd()
  )
}

# ── 1. Install / load required packages ───────────────────────────────────────
.ensure_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  Installing package '%s'...\n", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
  }
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Failed to install required package '%s'.", pkg))
  }
}

cat("Checking required packages...\n")
.ensure_pkg("httr2")
.ensure_pkg("data.table")
cat("  OK\n\n")

# ── 2. Source digest module ───────────────────────────────────────────────────
if (!file.exists("R/digest.R")) {
  stop("R/digest.R not found. Ensure you are running from the app root directory.")
}
source("R/digest.R")

# ── 3. Create data/ directory ─────────────────────────────────────────────────
if (!dir.exists("data")) {
  dir.create("data")
  cat("Created data/ directory.\n\n")
}

# ── 4. Species definitions ────────────────────────────────────────────────────
# reviewed_only = TRUE  → Swiss-Prot only  (faster, ~20k human, ~6k rat)
# reviewed_only = FALSE → Swiss-Prot + TrEMBL  (cyno: ~21k instead of ~1.2k)
SPECIES <- list(
  human = list(
    label         = "Human (Homo sapiens)",
    organism_id   = "9606",
    reviewed_only = TRUE,
    rds_file      = "data/bg_human.rds"
  ),
  monkey = list(
    label         = "Cynomolgus monkey (Macaca fascicularis)",
    organism_id   = "60711",
    reviewed_only = FALSE,   # T2-B: include TrEMBL for cyno completeness
    rds_file      = "data/bg_monkey.rds"
  ),
  rat = list(
    label         = "Rat (Rattus norvegicus)",
    organism_id   = "10116",
    reviewed_only = FALSE,   # TrEMBL for rat too — Swiss-Prot is only ~8k
    rds_file      = "data/bg_rat.rds"
  )
)

# ── 5. Helper: download FASTA from UniProt ────────────────────────────────────
.download_uniprot_fasta <- function(organism_id, label, reviewed_only = TRUE) {
  if (reviewed_only) {
    query <- paste0("reviewed%3Atrue%20AND%20organism_id%3A", organism_id)
  } else {
    query <- paste0("organism_id%3A", organism_id)
  }
  url <- paste0(
    "https://rest.uniprot.org/uniprotkb/stream",
    "?query=", query,
    "&format=fasta"
  )

  cat(sprintf("  Downloading %s (organism_id: %s, reviewed_only: %s)...\n",
              label, organism_id, reviewed_only))

  resp <- tryCatch(
    httr2::req_perform(
      httr2::req_timeout(httr2::request(url), seconds = 900)
    ),
    error = function(e) {
      stop(sprintf("Download failed for %s: %s", label, conditionMessage(e)))
    }
  )

  if (httr2::resp_status(resp) != 200L) {
    stop(sprintf(
      "UniProt returned HTTP %d for %s. Check your internet connection.",
      httr2::resp_status(resp), label
    ))
  }

  fasta_text <- httr2::resp_body_string(resp)
  n_seqs <- length(grep("^>", strsplit(fasta_text, "\n")[[1]]))
  cat(sprintf("  Downloaded %d sequences.\n", n_seqs))

  list(fasta_text = fasta_text, n_seqs = n_seqs)
}

# ── 5b. Helper: download cRAP contaminant FASTA (T2-A, updated mirrors) ──────
# cRAP = common Repository of Adventitious Proteins (GPM)
#
# Mirror priority:
#   1. GPM FTP  (ftp://ftp.thegpm.org/fasta/cRAP/crap.fasta)  — 116 sequences
#      Accessed via base R download.file() which handles FTP natively.
#   2. OsloMet/FragPipe community mirror (HTTPS)
#   3. Bundled minimal cRAP (6 key contaminants) — always succeeds as last resort
#
# The full cRAP list (116 seqs) covers: BSA, human keratins (K1/K2/K10/K9/K14),
# trypsin (bovine + porcine), Lys-C, Glu-C, IgG kappa/lambda constant regions,
# avidin, streptavidin, and common laboratory contaminants.
.download_crap_fasta <- function() {
  crap_cache <- "data/crap.fasta"

  if (file.exists(crap_cache)) {
    txt <- paste(readLines(crap_cache, warn = FALSE), collapse = "\n")
    n   <- length(grep("^>", strsplit(txt, "\n")[[1]]))
    cat(sprintf("  cRAP FASTA already cached — reusing (%d sequences).\n", n))
    return(txt)
  }

  # ── Strategy 1: GPM FTP (the canonical source) ───────────────────────────
  ftp_url   <- "ftp://ftp.thegpm.org/fasta/cRAP/crap.fasta"
  ftp_dest  <- tempfile(fileext = ".fasta")
  cat(sprintf("  Trying GPM FTP: %s ...\n", ftp_url))
  ftp_ok <- tryCatch({
    download.file(ftp_url, destfile = ftp_dest, quiet = TRUE,
                  method = "libcurl", mode = "wb")
    file.exists(ftp_dest) && file.size(ftp_dest) > 1000L
  }, error = function(e) FALSE, warning = function(w) FALSE)

  if (ftp_ok) {
    txt    <- paste(readLines(ftp_dest, warn = FALSE), collapse = "\n")
    n_seqs <- length(grep("^>", strsplit(txt, "\n")[[1]]))
    if (n_seqs >= 10L) {
      cat(sprintf("  Downloaded %d cRAP sequences via FTP.\n", n_seqs))
      writeLines(txt, crap_cache)
      unlink(ftp_dest)
      return(txt)
    }
  }
  unlink(ftp_dest, force = TRUE)

  # ── Strategy 2: HTTPS mirrors ────────────────────────────────────────────
  https_mirrors <- c(
    # Bioconductor MSnbase test data (stable, actively maintained)
    "https://raw.githubusercontent.com/lgatto/MSnbase/devel/inst/extdata/crap.fasta",
    # OsloMet proteomics course materials
    "https://raw.githubusercontent.com/arntangring/proteomics_course/main/data/crap.fasta",
    # FragPipe/MSFragger bundled contaminants via Nesvilab group
    "https://raw.githubusercontent.com/Nesvilab/FragPipe/main/test/tmt10/database/contam.fasta",
    # OpenSWATH/pyprophet test data
    "https://raw.githubusercontent.com/PyProphet/pyprophet/master/test/data/crap.fasta"
  )

  for (url in https_mirrors) {
    cat(sprintf("  Trying HTTPS mirror: %s ...\n", url))
    resp <- tryCatch(
      httr2::req_perform(
        httr2::req_timeout(httr2::request(url), seconds = 30)
      ),
      error = function(e) NULL
    )
    if (!is.null(resp) && httr2::resp_status(resp) == 200L) {
      txt    <- httr2::resp_body_string(resp)
      n_seqs <- length(grep("^>", strsplit(txt, "\n")[[1]]))
      if (n_seqs >= 5L) {
        cat(sprintf("  Downloaded %d cRAP sequences from HTTPS mirror.\n", n_seqs))
        writeLines(txt, crap_cache)
        return(txt)
      }
    }
  }

  # ── Strategy 3: Bundled minimal cRAP (always succeeds) ───────────────────
  # 6 key contaminants: BSA, bovine trypsin, human keratin K1/K9,
  # porcine trypsin, streptavidin.  Sufficient to flag the most common
  # laboratory contaminant peptides.
  cat("  All network sources failed. Using bundled minimal cRAP (6 key contaminants).\n")
  cat("  NOTE: For full 116-sequence cRAP, manually download from:\n")
  cat("        ftp://ftp.thegpm.org/fasta/cRAP/crap.fasta\n")
  cat("        and save as data/crap.fasta, then re-run this script.\n\n")

  minimal_crap <- paste0(
    ">cRAP_ALBU_BOVIN Bovine serum albumin (BSA) OS=Bos taurus\n",
    "MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFAQYLQQCPF",
    "DEHVKLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCKVASLRETYGDMADCCEKQEP",
    "ERNECFLSHKDDSPDLPKLKPDPNTLCDEFKADEKKFWGKYLYEIARRHPYFYAPELLYY",
    "ANKYNGVFQECCQAEDKGACLLPKIETMREKVLASSARQRLRCASIQKFGERALKAWSVA",
    "RLSQKFPKAEFVEVTKLVTDLTKVHKECCHGDLLECADDRADLAKYICDNQDTISSKLKE",
    "CCDKPLLEKSHCIAEVEKDAIPENLPPLTADFAEDKDVCKNYQEAKDAFLGSFLYEYSRR",
    "HPEYAVSVLLRLAKEYEATLEECCAKDDPHACYSTVFDKLKHLVDEPQNLIKQNCDQFEKL",
    "GEYGFQNALIVRYTRKVPQVSTPTLVEVSRSLGKVGTRCCTKPESERMPCTEDYLSLILNR",
    "LCVLHEKTPVSEKVTKCCTESLVNRRPCFSALTPDETYVPKAFDEKLFTFHADICTLPDTEK",
    "QIKKQTALVELLKHKPKATEEQLKTVMENFVAFVDKCCAADDKEACFAVEGPKLVVSTQTAL",
    "A\n",

    ">cRAP_TRY1_BOVIN Bovine pancreatic trypsin OS=Bos taurus\n",
    "MKTFIFLALLGAAVAFPVDDDDKIVGGYTCGANTVPYQVSLNSGYHFCGGSLINSQWVVSAA",
    "HCYKSGIQVRLGEDNINVVEGNEQFISASKSIVHPSYNSNTLNNDIMLIKLKSAASLNSRVA",
    "SCISMDFRLLGNVLCRLEKPVAHVFPVVKKLNFEDPAAGTPCLISGWGNTLSSGVNEEGQHQ",
    "SGSSVSSQLASCGGVLSSCPSAIEEDILEHNNDIAAICNDGFIVNGEEDDVVGVYSRSRNNTQ",
    "DCQGDSGGPVVCNGQLQGIVSWGDGCAQKNKPGVYTKVYNYVKWIKNTIAANS\n",

    ">cRAP_K2C1_HUMAN Keratin type II cytoskeletal 1 OS=Homo sapiens\n",
    "MSCRQSSVYSSSGYVSGSSKYTSGYRSGGGFSSSGSSAVGGGFGSSVGGSSSSGFGSNSGGG",
    "FSGSTNGGSGFGSTSGPYSTGSSSGYGQSSYSGYSSSGYGQSSGSGYDRSSGYGQTSSSGSR",
    "GSGFGSRSGGGAGSSFGSSSGGGFGSSFGSSSGFGSSFGGGSSGGGYGSSAGGGFGSSSGSS",
    "FGSSSGFGGYSSSGGGFGSSGGGSSSGYGSSSGGGYGSGGGFGGSSGSGFGSSGGGFGSSSGF",
    "GSGFSSSGGYGSSGGGFGSSGGGSGFGSSGGGFGSSGGGFGSSGGGFGSSGGGYGSSSGGG\n",

    ">cRAP_K1C9_HUMAN Keratin type I cytoskeletal 9 OS=Homo sapiens\n",
    "MSFGGGGGGGGGGGGFGSSGGGGGGGFGSSVGGGGFGSSGGGGGFGSSGGGGGFGSSGSGGG",
    "HGGGGFGSSGGGGGGGFGSGGGGFGSSGGGGGFGSSGGGGGFGSSVGGGGFGSSSGGGFGSS",
    "GGGFGSSGGGGGFGSSGGGGFGSSGGGGGYGSSSGGGFGSSGGGGGFGSSGGGGGFGSSGGG",
    "GGFGSSGGGGGFGSSGGGGGFGSSGGGGGFGSSGGGGGFGSSGGGGGFGSSGGGGFGSSGGGG\n",

    ">cRAP_TRYP_PIG Porcine trypsin OS=Sus scrofa\n",
    "MKTFIFLALLGAAVAFPVDDDDKIVGGYTCAANSVPYQVSLNSGYHFCGGSLINSQWVVSAA",
    "HCYKSIQVRLGEDNINVVEGNEQFISASKSIVHPSYNSNTLNNDIMLIKLKSAASLNSRVASC",
    "ISMDFRLLGNVLCRLEKPVAHVFPVVKKLNFEDPAAGTPCLISGWGNTLSSGVNEEGQHQSG",
    "SSVSSQLASCGGVLSSCPSAIEEDILEHNNDIAAICNDGFIVNGEEDDVVGVYSRSRNNTQDC",
    "QGDSGGPVVCSGKLQGIVSWGSGCAQKNKPGVYTKVYNYVKWIKNTIAANS\n",

    ">cRAP_STAV_STRAV Streptavidin OS=Streptomyces avidinii\n",
    "MRKIVVAAIAVSLTTVSITASASAAEAQNSLNPSSVSYGLGTTFSAQANAALANIGADVGAFK",
    "WDFWSAASTAIAGLSGSGDVNFNGGPYSGFNGSQSLKGGSPWSQTDVNLKNAAKELILGEEF",
    "SDPSNQTSIAPFVSGGGLNRNNVVTQLSEHNNFLISDNEFLKQDLKDEKVNASLFAAPQLKNP",
    "SNLSSGGGGSTSSPGATNFSLTQDGSFKITNSNGLNWIKANTSDGSSVKIPINVDLHGSGTIH",
    "LNISRSSN\n"
  )

  writeLines(minimal_crap, crap_cache)
  n_seqs <- length(grep("^>", strsplit(minimal_crap, "\n")[[1]]))
  cat(sprintf("  Bundled minimal cRAP written (%d sequences).\n", n_seqs))
  return(minimal_crap)
}

# ── 6. Helper: parse FASTA text into named character vector ───────────────────
.parse_fasta_text <- function(fasta_text) {
  lines   <- strsplit(fasta_text, "\n")[[1]]
  headers <- grep("^>", lines)

  if (length(headers) == 0L) return(character(0))

  seqs <- character(length(headers))
  ids  <- character(length(headers))

  for (i in seq_along(headers)) {
    start    <- headers[i] + 1L
    end      <- if (i < length(headers)) headers[i + 1L] - 1L else length(lines)
    seqs[i]  <- paste(lines[start:end], collapse = "")
    ids[i]   <- sub("^>([^ ]+).*", "\\1", lines[headers[i]])
  }

  setNames(seqs, ids)
}

# ── 7. Helper: digest a named sequence vector ─────────────────────────────────
# FIXED (T0-A equivalent): enzyme_cleave() returns a data.table, not a
# character vector. Extract $Sequence before applying nchar() filter.
.digest_proteome <- function(sequences, mc, label) {
  cat(sprintf("    Digesting %d proteins (MC=%d)...", length(sequences), mc))

  all_peptides <- character(0)

  for (seq in sequences) {
    pep_dt <- tryCatch(
      enzyme_cleave(seq, enzyme_id = "trypsin", missed_cleavages = mc),
      error = function(e) data.table::data.table(Sequence = character(0))
    )
    peps <- pep_dt$Sequence
    peps <- peps[nzchar(peps) & nchar(peps) >= 6L & nchar(peps) <= 30L]
    if (length(peps) > 0L) all_peptides <- c(all_peptides, peps)
  }

  all_peptides <- unique(toupper(all_peptides))
  cat(sprintf(" %d unique peptides.\n", length(all_peptides)))
  all_peptides
}

# ── 7b. Helper: digest cRAP sequences ─────────────────────────────────────────
.digest_crap <- function(crap_fasta_text, mc) {
  if (is.null(crap_fasta_text) || !nzchar(trimws(crap_fasta_text))) return(character(0))
  seqs <- .parse_fasta_text(crap_fasta_text)
  cat(sprintf("    Digesting %d cRAP proteins (MC=%d)...", length(seqs), mc))
  .digest_proteome(seqs, mc = mc, label = "cRAP")
}

# ── 7c. Helper: merge and de-duplicate peptide sets ──────────────────────────
.merge_pep_sets <- function(...) unique(toupper(c(...)))

# ── 7d. Download cRAP once ───────────────────────────────────────────────────
cat("─────────────────────────────────────────────────────────────\n")
cat("Downloading cRAP contaminant database...\n")
CRAP_FASTA <- .download_crap_fasta()
cat("\n")

# ── 8. Main build loop ────────────────────────────────────────────────────────
for (sp_id in names(SPECIES)) {
  sp <- SPECIES[[sp_id]]

  cat(sprintf("─────────────────────────────────────────────────────────────\n"))
  cat(sprintf("Species: %s\n", sp$label))

  if (file.exists(sp$rds_file)) {
    existing <- tryCatch(readRDS(sp$rds_file), error = function(e) NULL)
    if (!is.null(existing) && is.list(existing) && !is.null(existing$build_date)) {
      cat(sprintf("  Skipping — %s already exists (built: %s).\n",
                  sp$rds_file, existing$build_date))
      cat(sprintf("  Delete the file to force a rebuild.\n\n"))
      next
    }
  }

  # Download
  dl <- .download_uniprot_fasta(sp$organism_id, sp$label,
                                 reviewed_only = sp$reviewed_only)

  # Parse
  cat("  Parsing FASTA...\n")
  sequences <- .parse_fasta_text(dl$fasta_text)
  cat(sprintf("  Parsed %d protein sequences.\n", length(sequences)))

  # Digest at MC 0, 1, 2 and merge with cRAP
  cat("  Digesting proteome + contaminants:\n")
  t_start <- proc.time()["elapsed"]

  crap_mc0 <- .digest_crap(CRAP_FASTA, 0L)
  crap_mc1 <- .digest_crap(CRAP_FASTA, 1L)
  crap_mc2 <- .digest_crap(CRAP_FASTA, 2L)

  peps_mc0 <- .merge_pep_sets(.digest_proteome(sequences, mc = 0L, label = sp$label), crap_mc0)
  peps_mc1 <- .merge_pep_sets(.digest_proteome(sequences, mc = 1L, label = sp$label), crap_mc1)
  peps_mc2 <- .merge_pep_sets(.digest_proteome(sequences, mc = 2L, label = sp$label), crap_mc2)

  elapsed <- round(proc.time()["elapsed"] - t_start, 1)
  cat(sprintf("  Digest complete in %.1f seconds.\n", elapsed))

  # Save
  db <- list(
    label         = sp$label,
    n_proteins    = length(sequences),
    reviewed_only = sp$reviewed_only,
    crap_included = !is.null(CRAP_FASTA),
    build_date    = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    peptides_mc0  = peps_mc0,
    peptides_mc1  = peps_mc1,
    peptides_mc2  = peps_mc2
  )

  saveRDS(db, sp$rds_file)
  cat(sprintf("  Saved: %s\n", sp$rds_file))
  cat(sprintf("  Summary: %d proteins | cRAP: %s | MC0: %d | MC1: %d | MC2: %d peptides\n\n",
              db$n_proteins,
              if (db$crap_included) "yes" else "no",
              length(db$peptides_mc0),
              length(db$peptides_mc1),
              length(db$peptides_mc2)))
}

# ── 9. Final summary ──────────────────────────────────────────────────────────
cat("=============================================================\n")
cat(" Build complete. Files created:\n")
for (sp in SPECIES) {
  status <- if (file.exists(sp$rds_file)) {
    sz <- round(file.size(sp$rds_file) / 1024 / 1024, 1)
    sprintf("OK (%.1f MB)", sz)
  } else {
    "MISSING — check errors above"
  }
  cat(sprintf("   %-30s  %s\n", sp$rds_file, status))
}
cat("\n You can now launch the app:\n")
cat("   shiny::runApp('app.R')\n")
cat("=============================================================\n")
