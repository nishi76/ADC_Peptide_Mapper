# transitions.R — Fragment ion calculation and transition list generation

# ── Ion Series Calculation ────────────────────────────────────────────────────
# Calculates b, y (and optionally a) ion m/z values for a bare peptide sequence
# All product ions are singly charged (z=1)
# Returns a data.frame with columns: ion_type, ion_number, product_mz

calc_fragment_ions <- function(sequence, include_a_ions = FALSE) {
  residues <- strsplit(sequence, "")[[1]]
  n        <- length(residues)
  if (n < 2) return(data.frame())

  # Residue masses (bare, no mods — mods handled via mass_delta on precursor)
  res_masses <- AA_MONO_MASS[residues]
  res_masses[is.na(res_masses)] <- 0

  # b-ions: sum of N-terminal residues + proton
  # b_i = sum(res[1..i]) + PROTON
  b_mz <- numeric(n - 1)
  for (i in seq_len(n - 1)) {
    b_mz[i] <- sum(res_masses[1:i]) + PROTON_MASS
  }

  # y-ions: sum of C-terminal residues + water + proton
  # y_i = sum(res[(n-i+1)..n]) + WATER + PROTON
  y_mz <- numeric(n - 1)
  for (i in seq_len(n - 1)) {
    y_mz[i] <- sum(res_masses[(n - i + 1):n]) + WATER_MASS + PROTON_MASS
  }

  ions <- data.frame(
    ion_type   = c(rep("b", n - 1), rep("y", n - 1)),
    ion_number = c(seq_len(n - 1), seq_len(n - 1)),
    product_mz = c(b_mz, y_mz),
    stringsAsFactors = FALSE
  )

  # a-ions: b-ion - CO (27.99491 Da)
  if (include_a_ions) {
    a_ions <- data.frame(
      ion_type   = rep("a", n - 1),
      ion_number = seq_len(n - 1),
      product_mz = b_mz - 27.99491,
      stringsAsFactors = FALSE
    )
    ions <- rbind(ions, a_ions)
  }

  # Remove ions with m/z <= 0 or single-residue terminal ions (ion_number == 1 for b/a, == 1 for y)
  ions <- ions[ions$product_mz > 0, ]
  ions <- ions[!(ions$ion_number == 1 & ions$ion_type %in% c("b", "a")), ]
  ions <- ions[!(ions$ion_number == 1 & ions$ion_type == "y"), ]

  ions
}

# ── Top-N Ion Selection ───────────────────────────────────────────────────────
# Selects top N ions per ion type, preferring larger fragment numbers
# (larger fragments are generally more informative and less noisy)
select_top_ions <- function(ions_df, top_n = 5) {
  if (nrow(ions_df) == 0) return(ions_df)

  result <- do.call(rbind, lapply(unique(ions_df$ion_type), function(itype) {
    sub_df <- ions_df[ions_df$ion_type == itype, ]
    # Sort by ion_number descending (prefer larger fragments)
    sub_df <- sub_df[order(sub_df$ion_number, decreasing = TRUE), ]
    head(sub_df, top_n)
  }))

  result[order(result$ion_type, result$ion_number), ]
}

# ── Precursor m/z Calculation ─────────────────────────────────────────────────
calc_precursor_mz <- function(neutral_mass, charge) {
  (neutral_mass + charge * PROTON_MASS) / charge
}

# ── Collision Energy (Sciex empirical formula) ────────────────────────────────
calc_ce <- function(precursor_mz, charge) {
  if (charge == 2) {
    ce <- 0.0448 * precursor_mz - 2.0
  } else if (charge == 3) {
    ce <- 0.0533 * precursor_mz - 2.0
  } else if (charge == 4) {
    ce <- 0.0580 * precursor_mz - 2.0
  } else {
    ce <- 0.0448 * precursor_mz - 2.0
  }
  round(max(ce, 10), 1)  # floor at 10 eV
}

# ── Generate Full Transition List ─────────────────────────────────────────────
# peptides_dt: data.table from digest + modifications pipeline
#   Required cols: Chain, Sequence, ModifiedSequence, Start, End, Length,
#                  ModifiedMass, ModsApplied, UniqueToADC, ProteinName
# mode: "mrm" (2+/3+, b/y, top5) or "hr" (2+/3+/4+, b/y/a, top6)
# Returns a data.table of transitions

generate_transition_list <- function(peptides_dt, mode = "mrm") {
  if (nrow(peptides_dt) == 0) return(data.table::data.table())

  if (mode == "mrm") {
    charges    <- c(2L, 3L)
    top_n      <- 5L
    incl_a     <- FALSE
  } else {
    charges    <- c(2L, 3L, 4L)
    top_n      <- 6L
    incl_a     <- TRUE
  }

  rows <- list()

  for (i in seq_len(nrow(peptides_dt))) {
    pep       <- peptides_dt[i, ]
    bare_seq  <- pep$Sequence
    mod_seq   <- pep$ModifiedSequence
    neutral_m <- pep$ModifiedMass
    chain     <- pep$Chain
    protein   <- pep$ProteinName
    unique_f  <- pep$UniqueToADC
    mods_app  <- pep$ModsApplied

    # Fragment ions (based on bare sequence — mod mass captured in precursor)
    ions <- calc_fragment_ions(bare_seq, include_a_ions = incl_a)
    ions <- select_top_ions(ions, top_n = top_n)

    if (nrow(ions) == 0) next

    for (z in charges) {
      prec_mz <- calc_precursor_mz(neutral_m, z)
      ce_val  <- calc_ce(prec_mz, z)

      for (j in seq_len(nrow(ions))) {
        rows <- c(rows, list(data.table::data.table(
          ProteinName      = protein,
          Chain            = chain,
          PeptideSequence  = bare_seq,
          ModifiedSequence = mod_seq,
          PrecursorCharge  = z,
          PrecursorMz      = round(prec_mz, 5),
          ProductCharge    = 1L,
          ProductMz        = round(ions$product_mz[j], 5),
          FragmentIon      = paste0(toupper(ions$ion_type[j]), ions$ion_number[j]),
          CollisionEnergy  = ce_val,
          Modifications    = mods_app,
          UniqueToADC      = unique_f,
          PeptideLength    = pep$Length,
          Start            = pep$Start,
          End              = pep$End
        )))
      }
    }
  }

  if (length(rows) == 0) return(data.table::data.table())
  data.table::rbindlist(rows)
}
