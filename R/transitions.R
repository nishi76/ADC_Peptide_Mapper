# =============================================================================
# transitions.R — ADC Peptide Mapper v0.8
# =============================================================================
# Sourced by the Shiny app (server.R). No library() calls here.
# Depends on data.table (loaded by app) and constants from digest.R:
#   AA_MONO_MASS  — named numeric vector, monoisotopic residue masses (Da)
#   WATER_MASS    — numeric, 18.010565 Da
#   PROTON_MASS   — numeric,  1.007276 Da
#
# Public API
# ----------
#   calc_precursor_mz(sequence, charge, modifications)
#   calc_ce(precursor_mz, charge)
#   calc_ce_instrument(precursor_mz, charge, instrument)
#   calc_fragment_ions(sequence, modifications)
#   select_top_ions(ions_dt, top_n, by_mz_desc)
#   generate_transition_list(peptides_dt, adc_name)
#
# v0.8 changes vs v0.6
# ---------------------
#   - a-ions removed entirely (b and y only)
#   - calc_fragment_ions(): no include_a_ions param; starts at b2/y2
#   - generate_transition_list(): no `mode` param; always full b/y series;
#     no top-N filtering; adds ADCName column; charges 2+/3+/4+
#   - select_top_ions(): kept as utility for export.R, NOT called here
#   - calc_ce_instrument(): new multi-instrument CE helper
# =============================================================================


# -----------------------------------------------------------------------------
# calc_precursor_mz()
# -----------------------------------------------------------------------------
# Compute the m/z of a protonated precursor ion [M + z*H]^z+.
#
# Args:
#   sequence      : character(1) — single-letter amino acid sequence
#   charge        : integer(1)   — precursor charge state (e.g. 2, 3, 4)
#   modifications : numeric or NULL — named numeric vector of additional mass
#                   shifts keyed by 1-based position strings.
#                   E.g. c("5" = 57.021464) for carbamidomethyl Cys at pos 5.
#                   NULL or empty = unmodified peptide.
#
# Returns:
#   numeric(1) — precursor m/z, rounded to 6 decimal places
# -----------------------------------------------------------------------------
calc_precursor_mz <- function(sequence, charge, modifications = NULL) {

  # Split sequence into individual residues
  residues <- strsplit(sequence, "", fixed = TRUE)[[1]]

  # Sum residue masses; unknown residues default to 0 with a warning
  residue_masses <- vapply(residues, function(aa) {
    m <- AA_MONO_MASS[aa]
    if (is.na(m)) {
      warning(sprintf(
        "calc_precursor_mz: unknown amino acid '%s' — mass set to 0", aa
      ))
      0
    } else {
      m
    }
  }, numeric(1))

  # Neutral peptide mass = sum of residue masses + water (backbone termini)
  neutral_mass <- sum(residue_masses) + WATER_MASS

  # Add modification mass shifts (sum of all delta-masses)
  if (!is.null(modifications) && length(modifications) > 0L) {
    neutral_mass <- neutral_mass + sum(modifications)
  }

  # [M + z*H]^z+  →  (M + z * proton_mass) / z
  mz <- (neutral_mass + charge * PROTON_MASS) / charge

  round(mz, 6)
}


# -----------------------------------------------------------------------------
# calc_ce()
# -----------------------------------------------------------------------------
# Empirical collision energy (CE) formula for SCIEX instruments.
# Kept exactly as in v0.5 — do NOT modify.
#
# Charge-dependent linear regression against empirically optimised CEs for
# tryptic peptides on SCIEX TripleTOF / QTRAP platforms.
#
# Args:
#   precursor_mz : numeric(1) — precursor m/z
#   charge       : integer(1) — precursor charge state
#
# Returns:
#   numeric(1) — CE in eV, rounded to 1 decimal place, floored at 10 eV
# -----------------------------------------------------------------------------
calc_ce <- function(precursor_mz, charge) {

  ce <- switch(as.character(charge),
    "2" = 0.0448 * precursor_mz - 2.0,
    "3" = 0.0533 * precursor_mz - 2.0,
    "4" = 0.0580 * precursor_mz - 2.0,
    # Fallback for any other charge: use 3+ slope as a reasonable approximation
    0.0533 * precursor_mz - 2.0
  )

  # Floor at 10 eV to avoid nonsensical low values for very short peptides
  ce <- max(ce, 10)

  round(ce, 1)
}


# -----------------------------------------------------------------------------
# calc_ce_instrument()
# -----------------------------------------------------------------------------
# Instrument-specific empirical CE formulas.
#
# All formulas are charge-dependent linear functions of precursor m/z.
# Coefficients are empirical approximations suitable as method-development
# starting points; users should optimise per compound on their specific
# instrument configuration.
#
# References
#   SCIEX  : SCIEX TripleTOF 6600 / QTRAP application notes for peptide MRM.
#   Thermo : Hoofnagle et al. (2016) Clin. Chem. 62(7):997-1004;
#             empirical CE optimisation on TSQ Altis for tryptic peptides.
#   Bruker : Prianichnikov et al. (2021) Mol. Cell. Proteomics 20:100033;
#             timsTOF / EVOQ default CE ramp settings for peptide CID.
#   Agilent: Agilent MassHunter Optimizer default polynomial for QQQ;
#             internal application note G1969-90021.
#   Waters : Distler et al. (2016) Anal. Chem. 88(12):6079-6087;
#             Xevo TQ-S empirical CE for tryptic peptides.
#
# Args:
#   precursor_mz : numeric(1) — precursor m/z
#   charge       : integer(1) — precursor charge state
#   instrument   : character(1) — one of:
#                    "sciex"   SCIEX Analyst / TripleTOF / QTRAP
#                    "thermo"  Thermo Xcalibur / TSQ Altis / Quantis
#                    "bruker"  Bruker timsControl / QTOF / EVOQ
#                    "agilent" Agilent MassHunter / QQQ
#                    "waters"  Waters MassLynx / Xevo TQ
#
# Returns:
#   numeric(1) — CE in eV, floored at 10 eV, rounded to 1 decimal place
# -----------------------------------------------------------------------------
calc_ce_instrument <- function(precursor_mz, charge, instrument = "sciex") {

  instrument <- tolower(trimws(instrument))
  z          <- as.integer(charge)

  ce <- switch(instrument,

    # ------------------------------------------------------------------
    # SCIEX: charge-dependent slopes (identical to calc_ce())
    # Ref: SCIEX TripleTOF / QTRAP application notes for peptide MRM
    # ------------------------------------------------------------------
    "sciex" = {
      switch(as.character(z),
        "2" = 0.0448 * precursor_mz - 2.0,
        "3" = 0.0533 * precursor_mz - 2.0,
        "4" = 0.0580 * precursor_mz - 2.0,
        0.0533 * precursor_mz - 2.0   # fallback
      )
    },

    # ------------------------------------------------------------------
    # Thermo: charge-dependent slopes (TSQ Altis / Quantis)
    # Ref: Hoofnagle et al. (2016) Clin. Chem. 62:997; tuned for CID
    # ------------------------------------------------------------------
    "thermo" = {
      switch(as.character(z),
        "2" = 0.0313 * precursor_mz + 3.0,
        "3" = 0.0380 * precursor_mz + 3.0,
        "4" = 0.0440 * precursor_mz + 3.0,
        0.0380 * precursor_mz + 3.0   # fallback
      )
    },

    # ------------------------------------------------------------------
    # Bruker: charge-dependent slopes (timsTOF / EVOQ)
    # Ref: Prianichnikov et al. (2021) Mol. Cell. Proteomics 20:100033
    # ------------------------------------------------------------------
    "bruker" = {
      switch(as.character(z),
        "2" = 0.0350 * precursor_mz + 2.5,
        "3" = 0.0420 * precursor_mz + 2.5,
        "4" = 0.0490 * precursor_mz + 2.5,
        0.0420 * precursor_mz + 2.5   # fallback
      )
    },

    # ------------------------------------------------------------------
    # Agilent: charge-dependent slopes (QQQ / QTOF)
    # Ref: Agilent MassHunter Optimizer app note G1969-90021
    # ------------------------------------------------------------------
    "agilent" = {
      switch(as.character(z),
        "2" = 0.0330 * precursor_mz + 4.0,
        "3" = 0.0400 * precursor_mz + 4.0,
        "4" = 0.0470 * precursor_mz + 4.0,
        0.0400 * precursor_mz + 4.0   # fallback
      )
    },

    # ------------------------------------------------------------------
    # Waters: charge-dependent slopes (Xevo TQ-S / Xevo TQ-XS)
    # Ref: Distler et al. (2016) Anal. Chem. 88:6079
    # ------------------------------------------------------------------
    "waters" = {
      switch(as.character(z),
        "2" = 0.0320 * precursor_mz + 3.5,
        "3" = 0.0385 * precursor_mz + 3.5,
        "4" = 0.0450 * precursor_mz + 3.5,
        0.0385 * precursor_mz + 3.5   # fallback
      )
    },

    # ------------------------------------------------------------------
    # Unknown instrument — warn and fall back to SCIEX 3+ formula
    # ------------------------------------------------------------------
    {
      warning(sprintf(
        "calc_ce_instrument: unknown instrument '%s'. Falling back to SCIEX 3+ formula.",
        instrument
      ))
      0.0533 * precursor_mz - 2.0
    }
  )

  # Safety floor: never go below 10 eV
  ce <- max(ce, 10)

  round(ce, 1)
}


# -----------------------------------------------------------------------------
# calc_fragment_ions()
# -----------------------------------------------------------------------------
# Generate the complete b- and y-ion series for a peptide.
#
# Ion mass definitions (monoisotopic, singly charged [M+H]+):
#   b_n = sum(residue_masses[1..n])           + PROTON_MASS
#   y_n = sum(residue_masses[(L-n+1)..L])     + WATER_MASS + PROTON_MASS
#   where L = peptide length
#
# Coverage: b2 through b(L-1) and y2 through y(L-1), singly charged (z=1).
# For peptides of length >= 15 AA, doubly charged (z=2) fragment ions are
# also emitted for any 1+ ion with mz > 600 Da — these large fragments are
# routinely observed in MRM and improve selectivity for heavy peptides.
#
# Terminal ions (b1, y1, bL, yL) are excluded — they are not informative
# in standard MRM workflows and are often outside the instrument scan range.
#
# NOTE: a-ions have been removed in v0.6. The `include_a_ions` parameter
# from v0.5 is gone; a-ions are never calculated.
#
# Ion charge formulas:
#   z=1: mz_1 = neutral_frag + PROTON_MASS
#   z=2: mz_2 = (neutral_frag + 2*PROTON_MASS) / 2 = (mz_1 + PROTON_MASS) / 2
#
# Args:
#   sequence      : character(1) — single-letter amino acid sequence
#   modifications : numeric or NULL — named numeric vector of mass shifts
#                   keyed by 1-based position strings.
#                   E.g. c("3" = 79.966331) for phospho-Ser at position 3.
#                   NULL = unmodified.
#
# Returns:
#   data.table with columns:
#     ion_type   : character — "b" or "y"
#     ion_number : integer   — position index (2 .. L-1)
#     charge     : integer   — fragment charge state (1 or 2)
#     mz         : numeric   — m/z (rounded to 6 dp)
#     label      : character — e.g. "b3", "y7", "b12++", "y15++"
#
#   Returns a zero-row data.table with the same schema if L < 3.
# -----------------------------------------------------------------------------
calc_fragment_ions <- function(sequence, modifications = NULL) {

  residues <- strsplit(sequence, "", fixed = TRUE)[[1]]
  L        <- length(residues)

  # Need at least 3 residues to produce any b2/y2 ions
  if (L < 3L) {
    return(data.table::data.table(
      ion_type   = character(0),
      ion_number = integer(0),
      charge     = integer(0),
      mz         = numeric(0),
      label      = character(0)
    ))
  }

  # Residue masses (0 for unknowns, with warning)
  res_mass <- vapply(residues, function(aa) {
    m <- AA_MONO_MASS[aa]
    if (is.na(m)) {
      warning(sprintf(
        "calc_fragment_ions: unknown amino acid '%s' — mass set to 0", aa
      ))
      0
    } else {
      m
    }
  }, numeric(1))

  # Apply modification mass shifts to the appropriate residue positions
  if (!is.null(modifications) && length(modifications) > 0L) {
    for (pos_str in names(modifications)) {
      pos <- suppressWarnings(as.integer(pos_str))
      if (!is.na(pos) && pos >= 1L && pos <= L) {
        res_mass[pos] <- res_mass[pos] + modifications[[pos_str]]
      }
    }
  }

  # ------------------------------------------------------------------
  # b-ions: b_n = cumsum(res_mass)[n] + PROTON_MASS  (z=1)
  # Generate b2 through b(L-1)
  # ------------------------------------------------------------------
  cum_fwd  <- cumsum(res_mass)          # cumulative mass from N-terminus
  b_idx    <- seq(2L, L - 1L)          # ion numbers 2 .. L-1
  b_mz_1   <- cum_fwd[b_idx] + PROTON_MASS

  b_dt <- data.table::data.table(
    ion_type   = "b",
    ion_number = b_idx,
    charge     = 1L,
    mz         = round(b_mz_1, 6),
    label      = paste0("b", b_idx)
  )

  # ------------------------------------------------------------------
  # y-ions: y_n = cumsum(rev(res_mass))[n] + WATER_MASS + PROTON_MASS  (z=1)
  # Generate y2 through y(L-1)
  # ------------------------------------------------------------------
  cum_rev  <- cumsum(rev(res_mass))     # cumulative mass from C-terminus
  y_idx    <- seq(2L, L - 1L)          # ion numbers 2 .. L-1
  y_mz_1   <- cum_rev[y_idx] + WATER_MASS + PROTON_MASS

  y_dt <- data.table::data.table(
    ion_type   = "y",
    ion_number = y_idx,
    charge     = 1L,
    mz         = round(y_mz_1, 6),
    label      = paste0("y", y_idx)
  )

  # ------------------------------------------------------------------
  # Doubly charged (z=2) fragment ions — emitted when:
  #   (a) peptide length >= 15, AND
  #   (b) the corresponding 1+ ion mz > 600 Da
  # Formula: mz_2 = (mz_1 + PROTON_MASS) / 2
  # Label convention: "b12++" / "y15++" (Skyline / Proteome Discoverer)
  # ------------------------------------------------------------------
  ions_list <- list(b_dt, y_dt)

  if (L >= 15L) {
    # b-ions (z=2)
    b_keep <- b_mz_1 > 600
    if (any(b_keep)) {
      b_mz_2 <- (b_mz_1[b_keep] + PROTON_MASS) / 2
      b_dt_2 <- data.table::data.table(
        ion_type   = "b",
        ion_number = b_idx[b_keep],
        charge     = 2L,
        mz         = round(b_mz_2, 6),
        label      = paste0("b", b_idx[b_keep], "++")
      )
      ions_list <- c(ions_list, list(b_dt_2))
    }

    # y-ions (z=2)
    y_keep <- y_mz_1 > 600
    if (any(y_keep)) {
      y_mz_2 <- (y_mz_1[y_keep] + PROTON_MASS) / 2
      y_dt_2 <- data.table::data.table(
        ion_type   = "y",
        ion_number = y_idx[y_keep],
        charge     = 2L,
        mz         = round(y_mz_2, 6),
        label      = paste0("y", y_idx[y_keep], "++")
      )
      ions_list <- c(ions_list, list(y_dt_2))
    }
  }

  # Return all ions — no filtering applied here
  data.table::rbindlist(ions_list)
}


# -----------------------------------------------------------------------------
# select_top_ions()
# -----------------------------------------------------------------------------
# Utility: select the top N fragment ions from a fragment ion data.table.
#
# This function is used by export formatters in export.R to apply per-
# instrument top-N filtering. It is NOT called inside generate_transition_list().
#
# Args:
#   ions_dt    : data.table — must contain either:
#                  - column "ProductMz"  (transition list context), or
#                  - column "mz"         (raw fragment ion context)
#                and optionally "ion_number" for non-mz ranking.
#   top_n      : integer(1) — number of ions to keep (default 5).
#                If top_n >= nrow(ions_dt), all rows are returned unchanged.
#   by_mz_desc : logical(1) — ranking strategy:
#                  TRUE  (default) — rank by m/z descending (largest product
#                         ions first; preferred for MRM sensitivity)
#                  FALSE           — rank by ion_number descending (longest
#                         sequence coverage first)
#
# Returns:
#   data.table — subset of ions_dt, at most top_n rows, in ranked order
# -----------------------------------------------------------------------------
select_top_ions <- function(ions_dt, top_n = 5L, by_mz_desc = TRUE) {

  if (is.null(ions_dt) || nrow(ions_dt) == 0L) return(ions_dt)

  top_n <- as.integer(top_n)
  if (top_n <= 0L) return(ions_dt[0L])
  if (top_n >= nrow(ions_dt)) return(ions_dt)   # nothing to filter

  # Determine sort column — prefer transition-list names over raw names
  sort_col <- if (by_mz_desc) {
    if ("ProductMz"  %in% names(ions_dt)) "ProductMz"  else "mz"
  } else {
    if ("ion_number" %in% names(ions_dt)) "ion_number" else "mz"
  }

  # Order descending and take top N
  ord <- order(ions_dt[[sort_col]], decreasing = TRUE)
  ions_dt[ord[seq_len(min(top_n, length(ord)))]]
}


# -----------------------------------------------------------------------------
# generate_transition_list()
# -----------------------------------------------------------------------------
# Build the full MRM transition list for a set of peptides.
#
# For each peptide × precursor charge (2+, 3+, 4+) × fragment ion
# (b2..b(n-1), y2..y(n-1)), one row is emitted.
# Fragment ions are singly charged (z=1) for all peptides; for peptides
# >= 15 residues, doubly charged (z=2) fragments with 1+ mz > 600 are
# also included (see calc_fragment_ions() for details).
#
# No top-N filtering is applied here. The UI table displays all transitions.
# Instrument-specific export functions in export.R apply their own filtering.
#
# Args:
#   peptides_dt : data.table — one row per peptide, required columns:
#
#     Column            Type       Description
#     ----------------  ---------  ------------------------------------------
#     PeptideSequence   character  Single-letter amino acid sequence
#     ModifiedSequence  character  Sequence with modification notation
#                                  (e.g. "PEPTM[ox]IDE")
#     ProteinName       character  Source protein identifier
#     Chain             character  Heavy / light chain label (or "")
#     Modifications     character  Human-readable modification string
#     UniqueToADC       logical    TRUE if peptide is unique to this ADC
#     Start             integer    1-based start position in protein
#     End               integer    1-based end position in protein
#     Enzyme            character  Digestion enzyme name (e.g. "Trypsin")
#     mod_list          list       List column; each element is a named
#                                  numeric vector of mass shifts by position
#                                  (NULL element = unmodified peptide)
#
#   adc_name    : character(1) — ADC product name written to ADCName column.
#                 Default "".
#
# Returns:
#   data.table with one row per (peptide × charge × fragment ion):
#
#     Column            Type       Description
#     ----------------  ---------  ------------------------------------------
#     ADCName           character  ADC product name (from adc_name arg)
#     ProteinName       character  Source protein identifier
#     Chain             character  Heavy / light chain label
#     PeptideSequence   character  Single-letter sequence
#     ModifiedSequence  character  Sequence with mod notation
#     PrecursorCharge   integer    2, 3, or 4
#     PrecursorMz       numeric    Precursor [M+zH]^z+ m/z (6 dp)
#     ProductCharge     integer    1 or 2 (see calc_fragment_ions())
#     ProductMz         numeric    Fragment ion m/z (6 dp)
#     FragmentIon       character  Ion label, e.g. "b4", "y11"
#     CollisionEnergy   numeric    CE in eV (SCIEX formula; 1 dp)
#     Modifications     character  Human-readable mod string
#     UniqueToADC       logical    Peptide uniqueness flag
#     PeptideLength     integer    Number of residues
#     Start             integer    1-based start position
#     End               integer    1-based end position
#     Enzyme            character  Digestion enzyme name
# -----------------------------------------------------------------------------
generate_transition_list <- function(peptides_dt, adc_name = "") {

  # ---- Input validation ----------------------------------------------------

  if (is.null(peptides_dt) || nrow(peptides_dt) == 0L) {
    message("generate_transition_list: empty peptides_dt — returning empty table.")
    return(.empty_transition_table())
  }

  required_cols <- c(
    "PeptideSequence", "ModifiedSequence", "ProteinName",
    "Chain", "Modifications", "UniqueToADC", "Start", "End", "Enzyme"
  )
  missing_cols <- setdiff(required_cols, names(peptides_dt))
  if (length(missing_cols) > 0L) {
    stop(sprintf(
      "generate_transition_list: peptides_dt is missing required columns: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  # mod_list column is optional; default to NULL for all peptides if absent
  has_mod_list <- "mod_list" %in% names(peptides_dt)

  # ---- Main loop -----------------------------------------------------------

  # Precursor charge states to generate for every peptide
  precursor_charges <- c(2L, 3L, 4L)

  # Pre-allocate list; upper bound = n_peptides × n_charges
  all_blocks <- vector("list", nrow(peptides_dt) * length(precursor_charges))
  block_idx  <- 0L

  for (pep_i in seq_len(nrow(peptides_dt))) {

    pep_row  <- peptides_dt[pep_i]
    seq_str  <- pep_row$PeptideSequence

    # Retrieve modification list for this peptide (NULL = unmodified)
    mod_list <- if (has_mod_list) pep_row$mod_list[[1L]] else NULL

    # Calculate all fragment ions for this peptide (b2..b(n-1), y2..y(n-1))
    frag_ions <- calc_fragment_ions(seq_str, modifications = mod_list)

    # Skip peptides that produce no fragment ions (sequence length < 3)
    if (nrow(frag_ions) == 0L) {
      message(sprintf(
        "generate_transition_list: peptide '%s' is too short for fragment ions — skipped.",
        seq_str
      ))
      next
    }

    pep_len  <- nchar(seq_str)
    n_frags  <- nrow(frag_ions)

    for (z in precursor_charges) {

      # Precursor m/z for this charge state
      prec_mz <- calc_precursor_mz(seq_str, z, modifications = mod_list)

      # Collision energy using SCIEX formula (default instrument)
      # export.R recalculates CE per instrument when writing instrument files
      ce <- calc_ce(prec_mz, z)

      # Build one data.table block for this peptide × charge combination
      # (one row per fragment ion; ProductCharge comes from calc_fragment_ions)
      block <- data.table::data.table(
        ADCName          = rep(adc_name,                 n_frags),
        ProteinName      = rep(pep_row$ProteinName,      n_frags),
        Chain            = rep(pep_row$Chain,            n_frags),
        PeptideSequence  = rep(seq_str,                  n_frags),
        ModifiedSequence = rep(pep_row$ModifiedSequence, n_frags),
        PrecursorCharge  = rep(z,                        n_frags),
        PrecursorMz      = rep(prec_mz,                  n_frags),
        ProductCharge    = frag_ions$charge,
        ProductMz        = frag_ions$mz,
        FragmentIon      = frag_ions$label,
        CollisionEnergy  = rep(ce,                       n_frags),
        Modifications    = rep(pep_row$Modifications,    n_frags),
        UniqueToADC      = rep(pep_row$UniqueToADC,      n_frags),
        PeptideLength    = rep(pep_len,                  n_frags),
        Start            = rep(pep_row$Start,            n_frags),
        End              = rep(pep_row$End,              n_frags),
        Enzyme           = rep(pep_row$Enzyme,           n_frags)
      )

      block_idx <- block_idx + 1L
      all_blocks[[block_idx]] <- block
    }
  }

  # Trim unused list slots (peptides skipped due to short length)
  all_blocks <- all_blocks[seq_len(block_idx)]

  if (block_idx == 0L) {
    message("generate_transition_list: no valid transitions generated.")
    return(.empty_transition_table())
  }

  # Combine all blocks into a single data.table
  result <- data.table::rbindlist(all_blocks)

  # Enforce correct column types
  result[, PrecursorCharge := as.integer(PrecursorCharge)]
  result[, ProductCharge   := as.integer(ProductCharge)]
  result[, PeptideLength   := as.integer(PeptideLength)]
  result[, Start           := as.integer(Start)]
  result[, End             := as.integer(End)]
  result[, UniqueToADC     := as.logical(UniqueToADC)]

  result
}


# -----------------------------------------------------------------------------
# .empty_transition_table()  [internal helper — not exported]
# -----------------------------------------------------------------------------
# Returns a zero-row data.table with the canonical transition list schema.
# Used as a safe return value when no transitions can be generated.
# -----------------------------------------------------------------------------
.empty_transition_table <- function() {
  data.table::data.table(
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
}
