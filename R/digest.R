# digest.R — Trypsin digestion engine for ADC Peptide Mapper
# Monoisotopic residue masses (Da) — standard values
AA_MONO_MASS <- c(
  A = 71.03711,  R = 156.10111, N = 114.04293, D = 115.02694,
  C = 103.00919, E = 129.04259, Q = 128.05858, G = 57.02146,
  H = 137.05891, I = 113.08406, L = 113.08406, K = 128.09496,
  M = 131.04049, F = 147.06841, P = 97.05276,  S = 87.03203,
  T = 101.04768, W = 186.07931, Y = 163.06333, V = 99.06841
)

WATER_MASS  <- 18.01056   # H2O added to complete peptide
PROTON_MASS <- 1.007276   # proton mass for m/z calculation

# ── FASTA Parser ─────────────────────────────────────────────────────────────
# Returns a named list: list(chain_name = "SEQUENCE", ...)
parse_fasta <- function(fasta_text) {
  lines <- strsplit(fasta_text, "\n")[[1]]
  lines <- trimws(lines)
  lines <- lines[nchar(lines) > 0]

  chains    <- list()
  cur_name  <- NULL
  cur_seq   <- character(0)

  for (line in lines) {
    if (startsWith(line, ">")) {
      if (!is.null(cur_name)) {
        chains[[cur_name]] <- paste(cur_seq, collapse = "")
      }
      # Use full header (strip ">") as chain name
      cur_name <- sub("^>\\s*", "", line)
      cur_seq  <- character(0)
    } else {
      cur_seq <- c(cur_seq, toupper(gsub("[^A-Za-z]", "", line)))
    }
  }
  # Save last chain
  if (!is.null(cur_name) && length(cur_seq) > 0) {
    chains[[cur_name]] <- paste(cur_seq, collapse = "")
  }
  chains
}

# ── Trypsin Cleavage ─────────────────────────────────────────────────────────
# Cleave after K or R, NOT if next residue is P (KP/RP rule)
# Returns character vector of peptide sequences
trypsin_cleave <- function(sequence, missed_cleavages = 0) {
  seq_vec <- strsplit(sequence, "")[[1]]
  n       <- length(seq_vec)

  # Find cleavage sites: after K or R, not followed by P
  cut_after <- integer(0)
  for (i in seq_len(n - 1)) {
    if (seq_vec[i] %in% c("K", "R") && seq_vec[i + 1] != "P") {
      cut_after <- c(cut_after, i)
    }
  }
  # Always include end of sequence
  cut_after <- c(cut_after, n)

  # Build base peptides (0 missed cleavages)
  starts <- c(1, cut_after[-length(cut_after)] + 1)
  ends   <- cut_after

  base_peptides <- mapply(function(s, e) {
    list(
      sequence = paste(seq_vec[s:e], collapse = ""),
      start    = s,
      end      = e
    )
  }, starts, ends, SIMPLIFY = FALSE)

  # Add missed cleavage peptides
  all_peptides <- base_peptides
  if (missed_cleavages > 0 && length(base_peptides) > 1) {
    for (mc in seq_len(missed_cleavages)) {
      for (i in seq_len(length(base_peptides) - mc)) {
        j <- i + mc
        merged_seq   <- paste(sapply(base_peptides[i:j], `[[`, "sequence"), collapse = "")
        merged_start <- base_peptides[[i]]$start
        merged_end   <- base_peptides[[j]]$end
        all_peptides <- c(all_peptides, list(list(
          sequence = merged_seq,
          start    = merged_start,
          end      = merged_end
        )))
      }
    }
  }

  all_peptides
}

# ── Peptide Length Filter ────────────────────────────────────────────────────
filter_peptides_by_length <- function(peptide_list, min_len = 6, max_len = 30) {
  keep <- sapply(peptide_list, function(p) {
    len <- nchar(p$sequence)
    len >= min_len & len <= max_len
  })
  peptide_list[keep]
}

# ── Monoisotopic Neutral Mass ────────────────────────────────────────────────
# Calculates neutral monoisotopic mass of a bare peptide sequence
calc_peptide_mass <- function(sequence) {
  residues <- strsplit(sequence, "")[[1]]
  valid    <- residues[residues %in% names(AA_MONO_MASS)]
  if (length(valid) == 0) return(NA_real_)
  sum(AA_MONO_MASS[valid]) + WATER_MASS
}

# ── Full Digest Wrapper ──────────────────────────────────────────────────────
# Returns a data.table with columns:
#   Chain, Sequence, Start, End, Length, Mass
digest_fasta <- function(fasta_text, missed_cleavages = 0,
                          min_len = 6, max_len = 30) {
  chains <- parse_fasta(fasta_text)
  if (length(chains) == 0) return(data.table::data.table())

  result_list <- lapply(names(chains), function(chain_name) {
    seq_str  <- chains[[chain_name]]
    peptides <- trypsin_cleave(seq_str, missed_cleavages)
    peptides <- filter_peptides_by_length(peptides, min_len, max_len)

    if (length(peptides) == 0) return(NULL)

    data.table::data.table(
      Chain    = chain_name,
      Sequence = sapply(peptides, `[[`, "sequence"),
      Start    = sapply(peptides, `[[`, "start"),
      End      = sapply(peptides, `[[`, "end"),
      Length   = nchar(sapply(peptides, `[[`, "sequence")),
      Mass     = sapply(sapply(peptides, `[[`, "sequence"), calc_peptide_mass)
    )
  })

  result_list <- Filter(Negate(is.null), result_list)
  if (length(result_list) == 0) return(data.table::data.table())
  data.table::rbindlist(result_list)
}
