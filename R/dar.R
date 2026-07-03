# =============================================================================
# dar.R — ADC Peptide Mapper v0.8
# =============================================================================
# Drug-to-Antibody Ratio (DAR) distribution modeling for ADC peptide mapping.
#
# Real ADC batches are heterogeneous mixtures of DAR0 (naked antibody) through
# DARn (fully loaded) species.  For a conventional cysteine-conjugated IgG1,
# the theoretical maximum is DAR8 (8 interchain disulfide cysteines).
# For lysine-conjugated ADCs the distribution is even broader.
#
# This module generates transition lists for every DAR level, enabling
# targeted LC-MS/MS detection of each DAR species in a peptide map.
#
# Concept
# -------
# A tryptic peptide that contains a conjugation-site residue appears in two
# forms per site:
#   DAR=0 (naked):   no payload → mass = peptide_base_mass
#   DAR=1 (loaded):  one payload → mass = peptide_base_mass + payload_mass
#
# Peptides with multiple conjugation residues can carry 0, 1, 2, … payloads.
# The overall antibody DAR is the *sum* of loaded sites across all such
# peptides — this module models the per-peptide DAR contribution.
#
# Public API
# ----------
#   calc_dar_peptide_masses(sequence, payload_mass, conjugation_sites, dar_range)
#     → data.table: one row per DAR level showing mass and which sites are loaded
#
#   generate_dar_transitions(peptides_dt, payload_mass, dar_range,
#                            conjugation_type, residue, adc_name)
#     → data.table: full transition list with a DAR column
#
#   dar_summary_table(peptides_dt, payload_mass, dar_range,
#                     conjugation_type, residue)
#     → data.frame: summary of how many peptides contribute at each DAR level
# =============================================================================


# =============================================================================
# calc_dar_peptide_masses()
# =============================================================================
# Enumerate the mass of a peptide at each DAR level (0 to max conjugation sites).
#
# Args:
#   sequence          : character(1) — single-letter amino acid sequence
#   payload_mass      : numeric(1)   — mass shift (Da) per conjugated payload
#   conjugation_sites : integer vector — 1-based positions of conjugation residues
#                       within the peptide.  Use detect_conjugation_sites() to
#                       obtain this from ADCDB_PAYLOADS fields.
#   dar_range         : integer vector — DAR levels to include (default 0:4).
#                       Levels above length(conjugation_sites) are silently dropped.
#   base_mass         : numeric(1) or NULL — pre-computed neutral peptide mass.
#                       If NULL, computed internally via calc_peptide_mass().
#
# Returns:
#   data.table with columns:
#     Sequence         : character — peptide sequence
#     DAR              : integer   — drug-to-antibody ratio contribution
#                                    from this peptide (0 = naked, k = k drugs)
#     LoadedSites      : character — comma-separated positions carrying payload
#                                    (empty string for DAR=0)
#     NConjSites       : integer   — total conjugation sites in this peptide
#     NakedMass        : numeric   — neutral mass of the unloaded peptide
#     LoadedMass       : numeric   — neutral mass with k * payload_mass added
#     MassDelta        : numeric   — LoadedMass − NakedMass
# =============================================================================
calc_dar_peptide_masses <- function(sequence,
                                    payload_mass,
                                    conjugation_sites,
                                    dar_range  = 0L:4L,
                                    base_mass  = NULL) {
  if (is.null(base_mass)) {
    base_mass <- calc_peptide_mass(sequence)
  }

  n_sites <- length(conjugation_sites)
  # Clamp DAR range to what is achievable for this peptide
  valid_dar <- as.integer(dar_range[dar_range >= 0L & dar_range <= n_sites])

  if (length(valid_dar) == 0L) {
    return(data.table::data.table(
      Sequence = character(0), DAR = integer(0),
      LoadedSites = character(0), NConjSites = integer(0),
      NakedMass = numeric(0), LoadedMass = numeric(0), MassDelta = numeric(0)
    ))
  }

  rows <- lapply(valid_dar, function(k) {
    # For simplicity use the first k conjugation sites
    # (combinatorics for k-of-n with specific sites is deferred to future work)
    loaded <- if (k > 0L) conjugation_sites[seq_len(k)] else integer(0L)
    loaded_str <- if (k > 0L) paste(sort(loaded), collapse = ",") else ""

    data.table::data.table(
      Sequence     = sequence,
      DAR          = k,
      LoadedSites  = loaded_str,
      NConjSites   = n_sites,
      NakedMass    = round(base_mass, 6),
      LoadedMass   = round(base_mass + k * payload_mass, 6),
      MassDelta    = round(k * payload_mass, 6)
    )
  })

  data.table::rbindlist(rows)
}


# =============================================================================
# generate_dar_transitions()
# =============================================================================
# Build the full MRM transition list across all requested DAR levels for
# a set of peptides that contain conjugation sites.
#
# Strategy:
#   1. Identify which peptides carry ≥1 conjugation site.
#   2. For each such peptide × DAR level:
#      a. Compute the loaded peptide mass (base + k * payload_mass).
#      b. Build a synthetic mod_list that adds payload_mass to each loaded site.
#      c. Call generate_transition_list() to produce the full b/y transition rows.
#      d. Append a DAR column.
#   3. Peptides with 0 conjugation sites are included once at DAR=0 only.
#
# Args:
#   peptides_dt      : data.table — standard peptide table (from digest + mods)
#   payload_mass     : numeric(1) — mass shift per payload (Da)
#   dar_range        : integer vector — DAR levels to generate (default 0:4)
#   conjugation_type : character(1) — "cysteine" | "lysine" | "site_specific"
#   residue          : character(1) — conjugation residue ("C" or "K")
#   adc_name         : character(1) — written to ADCName column (default "")
#
# Returns:
#   data.table — same schema as generate_transition_list() plus a `DAR` column.
#   Returns an empty table if peptides_dt has no conjugation-site peptides.
# =============================================================================
generate_dar_transitions <- function(peptides_dt,
                                     payload_mass,
                                     dar_range        = 0L:4L,
                                     conjugation_type = "cysteine",
                                     residue          = "C",
                                     adc_name         = "") {

  if (is.null(peptides_dt) || nrow(peptides_dt) == 0L)
    return(.empty_dar_table())

  required <- c("PeptideSequence", "ModifiedSequence", "ProteinName",
                "Chain", "Modifications", "UniqueToADC", "Start", "End", "Enzyme")
  missing_cols <- setdiff(required, names(peptides_dt))
  if (length(missing_cols) > 0L)
    stop(sprintf("generate_dar_transitions: missing columns: %s",
                 paste(missing_cols, collapse = ", ")))

  dar_range     <- as.integer(dar_range)
  payload_mass  <- as.numeric(payload_mass)
  residue       <- toupper(trimws(residue))

  all_blocks <- vector("list", nrow(peptides_dt) * length(dar_range))
  block_idx  <- 0L

  for (pi in seq_len(nrow(peptides_dt))) {
    pep_row <- peptides_dt[pi]
    seq_str <- pep_row$PeptideSequence

    # Identify conjugation sites within this peptide
    conj_sites <- tryCatch(
      detect_conjugation_sites(seq_str, conjugation_type, residue),
      error = function(e) integer(0L)
    )

    # DAR levels achievable for this peptide
    max_dar  <- length(conj_sites)
    pep_dars <- if (max_dar == 0L) 0L else dar_range[dar_range <= max_dar]
    if (length(pep_dars) == 0L) pep_dars <- 0L

    for (k in pep_dars) {
      # Build modification list for k loaded sites
      base_mod_list <- if ("mod_list" %in% names(pep_row)) pep_row$mod_list[[1L]] else NULL

      if (k > 0L && length(conj_sites) > 0L) {
        loaded_pos <- conj_sites[seq_len(k)]
        payload_shifts <- setNames(rep(payload_mass, k), as.character(loaded_pos))
        if (!is.null(base_mod_list) && length(base_mod_list) > 0L) {
          mod_list_k <- c(unlist(base_mod_list), payload_shifts)
        } else {
          mod_list_k <- payload_shifts
        }
      } else {
        mod_list_k <- base_mod_list
      }

      # Temporarily build a one-row peptides_dt with the DAR-specific mod_list
      pep_row_k <- data.table::copy(pep_row)
      pep_row_k$mod_list <- list(mod_list_k)

      # Generate transitions for this peptide at DAR=k
      trans_k <- tryCatch(
        generate_transition_list(pep_row_k, adc_name = adc_name),
        error = function(e) NULL
      )

      if (is.null(trans_k) || nrow(trans_k) == 0L) next

      trans_k[, DAR := k]
      block_idx              <- block_idx + 1L
      all_blocks[[block_idx]] <- trans_k
    }
  }

  if (block_idx == 0L) return(.empty_dar_table())

  result <- data.table::rbindlist(all_blocks[seq_len(block_idx)], fill = TRUE)
  result[, DAR := as.integer(DAR)]
  data.table::setcolorder(result, c("DAR", setdiff(names(result), "DAR")))
  result
}


# =============================================================================
# dar_summary_table()
# =============================================================================
# Produce a concise summary: how many peptides at each DAR level, and the
# expected precursor mass range per DAR.
#
# Args:
#   peptides_dt      : data.table — standard peptide table
#   payload_mass     : numeric(1)
#   dar_range        : integer vector
#   conjugation_type : character(1)
#   residue          : character(1)
#
# Returns:
#   data.frame with columns: DAR, N_Peptides, MassMin_Da, MassMax_Da,
#                             TotalPayloadMass_Da
# =============================================================================
dar_summary_table <- function(peptides_dt,
                               payload_mass,
                               dar_range        = 0L:4L,
                               conjugation_type = "cysteine",
                               residue          = "C") {
  if (is.null(peptides_dt) || nrow(peptides_dt) == 0L)
    return(data.frame(DAR=integer(0), N_Peptides=integer(0),
                      MassMin_Da=numeric(0), MassMax_Da=numeric(0),
                      TotalPayloadMass_Da=numeric(0)))

  residue  <- toupper(trimws(residue))
  dar_range <- as.integer(dar_range)

  rows <- lapply(dar_range, function(k) {
    peps_at_k <- Filter(function(i) {
      sites <- tryCatch(
        detect_conjugation_sites(peptides_dt$PeptideSequence[i], conjugation_type, residue),
        error = function(e) integer(0L)
      )
      k <= length(sites)
    }, seq_len(nrow(peptides_dt)))

    n <- length(peps_at_k)
    if (n == 0L) {
      data.frame(DAR = k, N_Peptides = 0L,
                 MassMin_Da = NA_real_, MassMax_Da = NA_real_,
                 TotalPayloadMass_Da = round(k * payload_mass, 4))
    } else {
      base_masses <- if ("ModifiedMass" %in% names(peptides_dt)) {
        peptides_dt$ModifiedMass[peps_at_k]
      } else {
        vapply(peptides_dt$PeptideSequence[peps_at_k], calc_peptide_mass, numeric(1L))
      }
      loaded_masses <- base_masses + k * payload_mass
      data.frame(DAR = k, N_Peptides = n,
                 MassMin_Da = round(min(loaded_masses), 4),
                 MassMax_Da = round(max(loaded_masses), 4),
                 TotalPayloadMass_Da = round(k * payload_mass, 4))
    }
  })

  do.call(rbind, rows)
}


# =============================================================================
# .empty_dar_table()  [internal]
# =============================================================================
.empty_dar_table <- function() {
  base <- data.table::data.table(
    DAR              = integer(0),
    ADCName          = character(0),
    ProteinName      = character(0),
    Chain            = character(0),
    PeptideSequence  = character(0),
    ModifiedSequence = character(0),
    PrecursorCharge  = integer(0),
    PrecursorMz      = numeric(0),
    ProductCharge    = integer(0),
    ProductMz        = numeric(0),
    FragmentIon      = character(0),
    CollisionEnergy  = numeric(0),
    Modifications    = character(0),
    UniqueToADC      = logical(0),
    PeptideLength    = integer(0),
    Start            = integer(0),
    End              = integer(0),
    Enzyme           = character(0)
  )
  base
}

# =============================================================================
# End of dar.R
# =============================================================================
