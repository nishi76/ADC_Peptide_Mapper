# modifications.R — Modification mass tables and application logic

# ── Fixed Modifications (always applied) ─────────────────────────────────────
FIXED_MODS <- list(
  list(name = "Carbamidomethyl", residue = "C", mass = 57.02146, abbrev = "CAM")
)

# ── Variable Modification Definitions ────────────────────────────────────────
VAR_MOD_DEFS <- list(
  oxidation    = list(name = "Oxidation",     residue = "M",  mass = 15.99491,  abbrev = "Ox"),
  propionamide = list(name = "Propionamide",  residue = "C",  mass = 71.03711,  abbrev = "Prop"),
  nem          = list(name = "NEM",           residue = "C",  mass = 125.04768, abbrev = "NEM"),
  # Drug-linker payloads (on C, replaces CAM when enabled)
  mmae         = list(name = "MMAE-linker",   residue = "C",  mass = 715.3,     abbrev = "MMAE"),
  dm1          = list(name = "DM1-linker",    residue = "C",  mass = 738.0,     abbrev = "DM1"),
  dxd          = list(name = "DXd-linker",    residue = "C",  mass = 519.2,     abbrev = "DXd"),
  sn38         = list(name = "SN-38-linker",  residue = "C",  mass = 392.1,     abbrev = "SN38"),
  custom_drug  = list(name = "Custom-linker", residue = "C",  mass = NA_real_,  abbrev = "Drug")
)

# ── Predefined Special Modifications ─────────────────────────────────────────
SPECIAL_MOD_DEFS <- list(
  deamid_N     = list(name = "Deamidation",      residue = "N",      mass = 0.98402,   abbrev = "Deam",  nterm = FALSE),
  deamid_Q     = list(name = "Deamidation",      residue = "Q",      mass = 0.98402,   abbrev = "Deam",  nterm = FALSE),
  pyroglu_Q    = list(name = "Pyroglutamate",    residue = "Q",      mass = -17.02655, abbrev = "PyroQ", nterm = TRUE),
  pyroglu_E    = list(name = "Pyroglutamate",    residue = "E",      mass = -18.01056, abbrev = "PyroE", nterm = TRUE),
  acetyl_K     = list(name = "Acetylation",      residue = "K",      mass = 42.01057,  abbrev = "Ac",    nterm = FALSE),
  phospho_S    = list(name = "Phosphorylation",  residue = "S",      mass = 79.96633,  abbrev = "Phos",  nterm = FALSE),
  phospho_T    = list(name = "Phosphorylation",  residue = "T",      mass = 79.96633,  abbrev = "Phos",  nterm = FALSE),
  phospho_Y    = list(name = "Phosphorylation",  residue = "Y",      mass = 79.96633,  abbrev = "Phos",  nterm = FALSE)
)

# ── Build Active Modification List ───────────────────────────────────────────
# Combines fixed + selected variable + selected special + custom mods
# Returns a list of mod definitions to apply
build_active_mods <- function(
    var_mods_selected,       # character vector: e.g. c("oxidation","nem")
    drug_payload_key = NULL, # one of: "mmae","dm1","dxd","sn38","custom_drug"
    custom_drug_mass = NA,   # numeric, used when drug_payload_key == "custom_drug"
    special_mods_selected,   # character vector: e.g. c("deamid_N","pyroglu_Q")
    custom_mods_df = NULL    # data.frame with cols: residue, name, mass
) {
  mods <- FIXED_MODS  # always include CAM

  # Variable mods
  for (key in var_mods_selected) {
    if (key %in% names(VAR_MOD_DEFS)) {
      m <- VAR_MOD_DEFS[[key]]
      if (key == "custom_drug" && !is.na(custom_drug_mass)) {
        m$mass <- as.numeric(custom_drug_mass)
      }
      if (!is.na(m$mass)) mods <- c(mods, list(m))
    }
  }

  # Drug payload (if enabled and not already in var_mods_selected)
  if (!is.null(drug_payload_key) && !drug_payload_key %in% var_mods_selected) {
    m <- VAR_MOD_DEFS[[drug_payload_key]]
    if (!is.null(m)) {
      if (drug_payload_key == "custom_drug" && !is.na(custom_drug_mass)) {
        m$mass <- as.numeric(custom_drug_mass)
      }
      if (!is.na(m$mass)) mods <- c(mods, list(m))
    }
  }

  # Special mods
  for (key in special_mods_selected) {
    if (key %in% names(SPECIAL_MOD_DEFS)) {
      mods <- c(mods, list(SPECIAL_MOD_DEFS[[key]]))
    }
  }

  # Custom mods from user-built table
  if (!is.null(custom_mods_df) && nrow(custom_mods_df) > 0) {
    for (i in seq_len(nrow(custom_mods_df))) {
      mods <- c(mods, list(list(
        name    = as.character(custom_mods_df$name[i]),
        residue = toupper(as.character(custom_mods_df$residue[i])),
        mass    = as.numeric(custom_mods_df$mass[i]),
        abbrev  = as.character(custom_mods_df$name[i]),
        nterm   = FALSE
      )))
    }
  }

  mods
}

# ── Apply Modifications to a Peptide ─────────────────────────────────────────
# Given a bare sequence and active mod list, returns a list of modified variants.
# Each variant: list(modified_seq, mod_mass_delta, mods_applied)
# Strategy: fixed mods always applied; variable mods generate combinations
apply_modifications <- function(sequence, active_mods) {
  residues <- strsplit(sequence, "")[[1]]
  n        <- length(residues)

  # Separate fixed vs variable mods
  fixed_mods <- Filter(function(m) {
    # CAM is fixed; identify by abbrev
    m$abbrev == "CAM"
  }, active_mods)

  variable_mods <- Filter(function(m) m$abbrev != "CAM", active_mods)

  # Build base modified sequence with fixed mods applied
  base_seq_annotated <- residues
  base_mass_delta    <- 0
  base_mods_applied  <- character(0)

  for (mod in fixed_mods) {
    positions <- which(residues == mod$residue)
    if (length(positions) > 0) {
      for (pos in positions) {
        base_seq_annotated[pos] <- paste0(residues[pos], "[", mod$abbrev, "]")
        base_mass_delta         <- base_mass_delta + mod$mass
        base_mods_applied       <- c(base_mods_applied,
                                     paste0(mod$name, "@", residues[pos], pos))
      }
    }
  }

  # For variable mods: find all applicable positions
  # Build list of (position, mod) pairs
  var_sites <- list()
  for (mod in variable_mods) {
    is_nterm <- isTRUE(mod$nterm)
    if (is_nterm) {
      # Only apply to first residue if it matches
      if (residues[1] == mod$residue) {
        var_sites <- c(var_sites, list(list(pos = 1, mod = mod)))
      }
    } else {
      positions <- which(residues == mod$residue)
      for (pos in positions) {
        var_sites <- c(var_sites, list(list(pos = pos, mod = mod)))
      }
    }
  }

  # Generate all combinations: each site is either modified or not
  # To avoid combinatorial explosion, cap at 4 variable sites
  if (length(var_sites) > 4) var_sites <- var_sites[seq_len(4)]

  n_sites <- length(var_sites)
  combos  <- expand.grid(rep(list(c(FALSE, TRUE)), n_sites))

  variants <- list()
  for (ci in seq_len(nrow(combos))) {
    combo          <- as.logical(combos[ci, ])
    seq_ann        <- base_seq_annotated
    mass_delta     <- base_mass_delta
    mods_applied   <- base_mods_applied

    for (si in seq_along(var_sites)) {
      if (combo[si]) {
        site <- var_sites[[si]]
        pos  <- site$pos
        mod  <- site$mod
        # Replace residue annotation (may already have fixed mod annotation)
        current <- seq_ann[pos]
        # Insert variable mod tag
        seq_ann[pos]  <- sub(paste0("^", residues[pos]),
                             paste0(residues[pos], "[", mod$abbrev, "]"),
                             current)
        mass_delta    <- mass_delta + mod$mass
        mods_applied  <- c(mods_applied,
                           paste0(mod$name, "@", residues[pos], pos))
      }
    }

    variants <- c(variants, list(list(
      modified_seq  = paste(seq_ann, collapse = ""),
      mass_delta    = mass_delta,
      mods_applied  = if (length(mods_applied) > 0)
                        paste(mods_applied, collapse = "; ")
                      else "None"
    )))
  }

  # If no variable sites exist, combos has 1 row (all FALSE) — still returns base variant
  # Deduplicate variants by modified_seq
  seen <- character(0)
  unique_variants <- list()
  for (v in variants) {
    if (!v$modified_seq %in% seen) {
      seen            <- c(seen, v$modified_seq)
      unique_variants <- c(unique_variants, list(v))
    }
  }

  # Safety: always return at least the base (fixed-mods-only) variant
  if (length(unique_variants) == 0) {
    unique_variants <- list(list(
      modified_seq = paste(base_seq_annotated, collapse = ""),
      mass_delta   = base_mass_delta,
      mods_applied = if (length(base_mods_applied) > 0)
                       paste(base_mods_applied, collapse = "; ")
                     else "None"
    ))
  }

  unique_variants
}

# ── Calculate Modified Peptide Mass ──────────────────────────────────────────
calc_modified_mass <- function(base_mass, mass_delta) {
  base_mass + mass_delta
}
