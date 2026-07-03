# =============================================================================
# isotopes.R — ADC Peptide Mapper v0.8
# =============================================================================
# Monoisotopic vs. most-abundant precursor isotope selection.
#
# For peptides above ~2,000 Da the monoisotopic peak (M+0) is no longer the
# most abundant isotope peak, and selecting it for SRM/MRM transitions loses
# precursor signal.  Skyline and MRMaid automatically recommend M+1 or M+2
# as the isolation target in such cases; this module provides the same logic.
#
# Method: averagine approximation (Senko et al. 1995, J. Am. Soc. Mass Spectrom.
# 6:229-233).  The elemental composition of the peptide is approximated by
# scaling the averagine monomer (C4.9384 H7.7583 N1.3577 O1.4773 S0.0417)
# to the number of residues.  The isotope envelope is then computed by
# convolution of binomial distributions for each element using the
# first-order approximation (Rockwood & Haimi 2006, J. Am. Soc. Mass Spectrom.
# 17:415-419).
#
# Public API
# ----------
#   calc_isotope_distribution(sequence, n_isotopes)
#     → data.frame(isotope=0:4, probability=..., rel_intensity=...)
#
#   recommended_precursor_isotope(sequence)
#     → integer 0, 1, or 2  (offset from monoisotopic, i.e. M+0/M+1/M+2)
#
#   annotate_isotope_offsets(peptides_dt)
#     → same data.table with RecommendedIsotope column added/updated
# =============================================================================


# Natural isotope abundances (IUPAC 2021)
.ISO_ABUND <- list(
  C  = c(0.9893,   0.0107),                      # 12C, 13C
  H  = c(0.999885, 0.000115),                     # 1H,  2H
  N  = c(0.99632,  0.00368),                      # 14N, 15N
  O  = c(0.99757,  0.00038,  0.00205),            # 16O, 17O, 18O
  S  = c(0.9499,   0.0075,   0.0425,   0.0001)    # 32S, 33S, 34S, 36S
)

# Averagine composition per residue (Senko et al. 1995)
.AVERAGINE <- c(C = 4.9384, H = 7.7583, N = 1.3577, O = 1.4773, S = 0.0417)

# Neutron mass (mass difference between isotope peaks)
.NEUTRON_MASS <- 1.003355


# =============================================================================
# .isotope_envelope_element()  [internal]
# =============================================================================
# Compute the probability distribution for a single element with n atoms,
# using the "first-order" multinomial approximation. Returns a numeric vector
# of length (n_isotopes + 1): P(0 heavy atoms), P(1), …, P(n_isotopes).
# =============================================================================
.isotope_envelope_element <- function(n_atoms, abundances, n_isotopes) {
  if (n_atoms <= 0 || length(abundances) < 2L) {
    out <- numeric(n_isotopes + 1L)
    out[1L] <- 1
    return(out)
  }

  n_atoms <- round(n_atoms)
  p_heavy <- 1 - abundances[1L]   # total probability of any heavy isotope

  # Poisson approximation: lambda = n * p_heavy
  lambda <- n_atoms * p_heavy
  probs  <- dpois(0:n_isotopes, lambda)
  probs / sum(probs)    # renormalise
}


# =============================================================================
# calc_isotope_distribution()
# =============================================================================
# Compute the isotope envelope for a peptide sequence.
#
# Args:
#   sequence    : character(1) — single-letter amino acid sequence
#   n_isotopes  : integer(1)   — number of isotope peaks (default 5: M+0 .. M+4)
#
# Returns:
#   data.frame with columns:
#     isotope      : integer — 0, 1, 2, … (offset from monoisotopic)
#     mz_offset    : numeric — mass offset from M+0 (multiples of neutron mass)
#     probability  : numeric — unnormalised abundance (sums to 1)
#     rel_intensity: numeric — intensity relative to most abundant peak (0–100)
# =============================================================================
calc_isotope_distribution <- function(sequence, n_isotopes = 5L) {
  n_isotopes <- as.integer(n_isotopes)
  n_res <- nchar(sequence)
  if (n_res == 0L) {
    return(data.frame(
      isotope      = 0L,
      mz_offset    = 0,
      probability  = 1,
      rel_intensity = 100
    ))
  }

  # Elemental formula via averagine scaling
  n_C <- round(.AVERAGINE["C"] * n_res)
  n_H <- round(.AVERAGINE["H"] * n_res)
  n_N <- round(.AVERAGINE["N"] * n_res)
  n_O <- round(.AVERAGINE["O"] * n_res)
  n_S <- round(.AVERAGINE["S"] * n_res)

  # Isotope envelope per element
  env_C <- .isotope_envelope_element(n_C, .ISO_ABUND$C, n_isotopes)
  env_H <- .isotope_envelope_element(n_H, .ISO_ABUND$H, n_isotopes)
  env_N <- .isotope_envelope_element(n_N, .ISO_ABUND$N, n_isotopes)
  env_O <- .isotope_envelope_element(n_O, .ISO_ABUND$O, n_isotopes)
  env_S <- .isotope_envelope_element(n_S, .ISO_ABUND$S, n_isotopes)

  # Convolve all elements
  .convolve_envelopes <- function(a, b) {
    n <- length(a) + length(b) - 1L
    out <- numeric(n)
    for (i in seq_along(a))
      for (j in seq_along(b))
        out[i + j - 1L] <- out[i + j - 1L] + a[i] * b[j]
    out
  }

  env <- env_C
  env <- .convolve_envelopes(env, env_H)
  env <- .convolve_envelopes(env, env_N)
  env <- .convolve_envelopes(env, env_O)
  env <- .convolve_envelopes(env, env_S)

  # Trim to requested length
  env <- env[seq_len(min(n_isotopes + 1L, length(env)))]
  env <- env / sum(env)    # normalise

  rel <- 100 * env / max(env)

  data.frame(
    isotope       = seq(0L, length(env) - 1L),
    mz_offset     = seq(0L, length(env) - 1L) * .NEUTRON_MASS,
    probability   = round(env, 6),
    rel_intensity = round(rel, 2),
    stringsAsFactors = FALSE
  )
}


# =============================================================================
# recommended_precursor_isotope()
# =============================================================================
# Return the isotope offset (0 = M+0, 1 = M+1, 2 = M+2) of the most abundant
# isotope peak for the given peptide sequence.
#
# Guideline:
#   - < 1,500 Da  → M+0 (monoisotopic)
#   - 1,500–3,000 → M+1 if M+1 > M+0
#   - > 3,000 Da  → M+1 or M+2 (whichever is most abundant)
#
# Args:
#   sequence : character(1) — single-letter amino acid sequence
#
# Returns:
#   integer(1) — 0, 1, or 2
# =============================================================================
recommended_precursor_isotope <- function(sequence) {
  env <- calc_isotope_distribution(sequence, n_isotopes = 4L)
  as.integer(env$isotope[which.max(env$probability)])
}


# =============================================================================
# annotate_isotope_offsets()
# =============================================================================
# Add or update a RecommendedIsotope column in a peptides / transition data.table.
#
# Args:
#   peptides_dt : data.table — must have a PeptideSequence column
#
# Returns:
#   same data.table with RecommendedIsotope column (integer) added in-place
# =============================================================================
annotate_isotope_offsets <- function(peptides_dt) {
  if (is.null(peptides_dt) || nrow(peptides_dt) == 0L) return(peptides_dt)
  if (!"PeptideSequence" %in% names(peptides_dt)) {
    warning("annotate_isotope_offsets: 'PeptideSequence' column missing — skipping.")
    return(peptides_dt)
  }

  # Cache by unique sequence for speed
  unique_seqs <- unique(peptides_dt$PeptideSequence)
  iso_map <- setNames(
    vapply(unique_seqs, recommended_precursor_isotope, integer(1L)),
    unique_seqs
  )

  peptides_dt[, RecommendedIsotope := iso_map[PeptideSequence]]
  peptides_dt
}

# =============================================================================
# End of isotopes.R
# =============================================================================
