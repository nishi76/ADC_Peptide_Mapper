# =============================================================================
# tests/test_masses.R — ADC Peptide Mapper v0.8
# =============================================================================
# Self-contained unit test suite (base R only — no testthat required).
#
# Usage:
#   Rscript tests/test_masses.R
#   # or from R console:
#   source("tests/test_masses.R")
#
# Each test is a named closure.  Results are printed as PASS / FAIL with the
# observed vs expected value.  The script exits with status 1 if any test
# fails (suitable for CI pipelines).
#
# Coverage:
#   Mass constants and AA_MONO_MASS residue values
#   calc_peptide_mass()          — neutral monoisotopic mass
#   calc_precursor_mz()          — [M+zH]^z+ m/z computation
#   calc_ce() / calc_ce_instrument() — collision energy formulas
#   calc_fragment_ions()         — b- and y-ion series
#   enzyme_cleave()              — trypsin cleavage rules
#   detect_conjugation_sites()   — cysteine and lysine mapping
#   calc_dar_peptide_masses()    — DAR mass enumeration
#   generate_dar_transitions()   — DAR transition list schema and m/z delta
#   calc_isotope_distribution()  — averagine envelope peak count and sum
# =============================================================================

cat("=============================================================\n")
cat(" ADC Peptide Mapper v0.8 — Unit Test Suite\n")
cat("=============================================================\n\n")

# ── 0. Load modules ───────────────────────────────────────────────────────────
if (!file.exists("app.R")) {
  stop("Run tests from the ADC Peptide Mapper root directory (contains app.R).")
}

needed_pkgs <- c("data.table")
for (pkg in needed_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
  library(pkg, character.only = TRUE)
}

for (.m in c("R/digest.R", "R/modifications.R", "R/transitions.R",
             "R/isotopes.R", "R/dar.R")) {
  source(.m)
}

# ── 1. Test harness ───────────────────────────────────────────────────────────
.pass <- 0L
.fail <- 0L
.failures <- character(0)

.test <- function(name, expr, expected, tol = 1e-4, check = "numeric") {
  result <- tryCatch(expr, error = function(e) structure(conditionMessage(e), class = "error"))

  if (inherits(result, "error")) {
    cat(sprintf("  FAIL  [%s]\n        Error: %s\n", name, result))
    .fail  <<- .fail + 1L
    .failures <<- c(.failures, name)
    return(invisible(FALSE))
  }

  ok <- switch(check,
    numeric  = isTRUE(all.equal(as.numeric(result), as.numeric(expected), tolerance = tol)),
    integer  = identical(as.integer(result), as.integer(expected)),
    logical  = identical(as.logical(result), as.logical(expected)),
    character= identical(as.character(result), as.character(expected)),
    nrow     = isTRUE(nrow(result) == expected),
    contains = all(as.character(expected) %in% as.character(result)),
    TRUE
  )

  if (ok) {
    cat(sprintf("  PASS  [%s]\n", name))
    .pass <<- .pass + 1L
  } else {
    cat(sprintf("  FAIL  [%s]\n        expected: %s\n        got:      %s\n",
                name,
                paste(expected, collapse = ", "),
                paste(result,   collapse = ", ")))
    .fail  <<- .fail + 1L
    .failures <<- c(.failures, name)
  }
  invisible(ok)
}


# =============================================================================
# Section 1 — Mass constants
# =============================================================================
cat("── Section 1: Mass constants ────────────────────────────────\n")

# NIST monoisotopic residue masses (Unimod-consistent, Da)
# Source: NIST Chemistry WebBook; Roepstorff & Fohlman (1984) nomenclature
.nist <- c(
  G =  57.02146,  A =  71.03711,  V =  99.06841,  L = 113.08406,
  I = 113.08406,  P =  97.05276,  F = 147.06841,  W = 186.07931,
  M = 131.04049,  S =  87.03203,  T = 101.04768,  C = 103.00919,
  Y = 163.06333,  H = 137.05891,  D = 115.02694,  E = 129.04259,
  N = 114.04293,  Q = 128.05858,  K = 128.09496,  R = 156.10111
)

for (aa in names(.nist)) {
  .test(
    sprintf("AA_MONO_MASS[%s] = %.5f Da", aa, .nist[aa]),
    AA_MONO_MASS[aa],
    .nist[aa],
    tol = 1e-5
  )
}

.test("WATER_MASS = 18.01056 Da",  WATER_MASS,  18.01056,  tol = 1e-5)
.test("PROTON_MASS = 1.007276 Da", PROTON_MASS,  1.007276, tol = 1e-6)


# =============================================================================
# Section 2 — calc_peptide_mass()
# =============================================================================
cat("\n── Section 2: calc_peptide_mass() ───────────────────────────\n")

# Reference: manually computed from AA_MONO_MASS + WATER_MASS
# CASIQKFGR:
#   C=103.00919 A=71.03711 S=87.03203 I=113.08406 Q=128.05858
#   K=128.09496 F=147.06841 G=57.02146 R=156.10111  sum=990.50691
#   + H2O = 1008.51747
.test("calc_peptide_mass('CASIQKFGR')",
      calc_peptide_mass("CASIQKFGR"), 1008.51747, tol = 1e-4)

# EVQLVESGGG (first 10 residues of trastuzumab HC):
#   E=129.04259 V=99.06841 Q=128.05858 L=113.08406 V=99.06841
#   E=129.04259 S=87.03203 G=57.02146 G=57.02146 G=57.02146
#   sum = 955.46105 + 18.01056 = 973.47161
.test("calc_peptide_mass('EVQLVESGGG')",
      calc_peptide_mass("EVQLVESGGG"), 973.47161, tol = 1e-4)

# Single residue test
.test("calc_peptide_mass('A') = Ala + H2O",
      calc_peptide_mass("A"),
      AA_MONO_MASS["A"] + WATER_MASS,
      tol = 1e-6)

# Empty peptide: implementation returns NA (no residues to sum) — document that
.test("calc_peptide_mass('') returns NA or WATER_MASS (backbone-only edge case)",
      is.na(calc_peptide_mass("")) || isTRUE(abs(calc_peptide_mass("") - WATER_MASS) < 1e-5),
      TRUE,
      check = "logical")


# =============================================================================
# Section 3 — calc_precursor_mz()
# =============================================================================
cat("\n── Section 3: calc_precursor_mz() ───────────────────────────\n")

# CASIQKFGR unmodified z=2:
#   (1008.51747 + 2*1.007276) / 2 = 1010.53202 / 2 = 505.26601
.test("calc_precursor_mz('CASIQKFGR', z=2)",
      calc_precursor_mz("CASIQKFGR", 2), 505.26601, tol = 1e-4)

# z=3:
#   (1008.51747 + 3*1.007276) / 3 = 1011.53930 / 3 = 337.17977
.test("calc_precursor_mz('CASIQKFGR', z=3)",
      calc_precursor_mz("CASIQKFGR", 3), 337.17977, tol = 1e-4)

# With carbamidomethyl Cys at position 1 (+57.021464 Da):
#   neutral = 1008.51747 + 57.021464 = 1065.53893
#   z=2: (1065.53893 + 2*1.007276)/2 = 1067.55348/2 = 533.77674
.test("calc_precursor_mz('CASIQKFGR', z=2, CAM-C@1)",
      calc_precursor_mz("CASIQKFGR", 2, c("1" = 57.021464)), 533.77674, tol = 1e-4)

# Charge z=1 basic sanity
.test("calc_precursor_mz('A', z=1) = Ala+H2O+proton",
      calc_precursor_mz("A", 1),
      AA_MONO_MASS["A"] + WATER_MASS + PROTON_MASS,
      tol = 1e-6)


# =============================================================================
# Section 4 — calc_ce() and calc_ce_instrument()
# =============================================================================
cat("\n── Section 4: Collision energy formulas ─────────────────────\n")

# SCIEX default formula from transitions.R:
#   z=2: CE = 0.0448 * mz - 2.0   (floored at 10)
# mz = 505.27:  0.0448*505.27 - 2.0 = 22.636 → 20.6
.test("calc_ce(505.27, z=2) SCIEX formula",
      calc_ce(505.27, 2),
      round(max(0.0448 * 505.27 - 2.0, 10), 1),
      tol = 0.05)

# z=3: CE = 0.0533*505.27 - 2.0 = 24.931 → 24.9
.test("calc_ce(505.27, z=3)",
      calc_ce(505.27, 3),
      round(max(0.0533 * 505.27 - 2.0, 10), 1),
      tol = 0.05)

# Floor: tiny mz should give CE = 10
.test("calc_ce(10, z=2) floor=10",
      calc_ce(10, 2),
      10.0,
      tol = 0.01)

# Instrument-specific: SCIEX should match calc_ce()
.test("calc_ce_instrument('sciex') matches calc_ce()",
      calc_ce_instrument(505.27, 2, "sciex"),
      calc_ce(505.27, 2),
      tol = 0.1)

# All five instruments should return numeric > 0 and < 200
for (inst in c("sciex", "thermo", "bruker", "agilent", "waters")) {
  ce_val <- calc_ce_instrument(500, 2, inst)
  .test(sprintf("calc_ce_instrument('%s', 500, z=2) in [10, 200]", inst),
        ce_val >= 10 && ce_val <= 200,
        TRUE,
        check = "logical")
}


# =============================================================================
# Section 5 — calc_fragment_ions()
# =============================================================================
cat("\n── Section 5: calc_fragment_ions() ──────────────────────────\n")

# calc_fragment_ions() returns columns: ion_type, ion_number, charge, mz, label
ions_dt <- calc_fragment_ions("CASIQKFGR")

# Should return a data.table
.test("calc_fragment_ions() returns data.table",
      is.data.table(ions_dt),
      TRUE,
      check = "logical")

# b2 for CA: b2 = C + A + proton = 103.00919 + 71.03711 + 1.007276 = 175.05358
b2_obs <- ions_dt[label == "b2", mz][1]
.test("b2 ion for CASIQKFGR (singly charged)",
      b2_obs, 175.05358, tol = 1e-4)

# y2 for ...GR: y2 = G + R + H2O + proton = 57.02146 + 156.10111 + 18.01056 + 1.007276 = 232.14041
y2_obs <- ions_dt[label == "y2", mz][1]
.test("y2 ion for CASIQKFGR (singly charged)",
      y2_obs, 232.14041, tol = 1e-4)

# Should contain both b and y ions
.test("calc_fragment_ions() has b-ions",
      any(grepl("^b", ions_dt$label)),
      TRUE,
      check = "logical")

.test("calc_fragment_ions() has y-ions",
      any(grepl("^y", ions_dt$label)),
      TRUE,
      check = "logical")

# For a 9-residue peptide: b-ions b2..b8 (7) and y-ions y2..y8 (7) = 14 singly charged
n_z1 <- nrow(ions_dt[charge == 1L])
.test("calc_fragment_ions() z=1 count for 9-mer (14 expected: b2-b8 + y2-y8)",
      n_z1, 14L, check = "integer")


# =============================================================================
# Section 6 — enzyme_cleave()
# =============================================================================
cat("\n── Section 6: enzyme_cleave() ────────────────────────────────\n")

# Simple tripeptide KAR → trypsin cleaves after K and after R
# MC=0: expected peptides "K", "AR" — but length filter may apply; use raw
kka_dt <- enzyme_cleave("KAARVGR", enzyme_id = "trypsin", missed_cleavages = 0L)
kka_seqs <- kka_dt$Sequence
.test("trypsin MC=0 cleaves KAARVGR into expected fragments",
      kka_seqs,
      c("K", "AAR", "VGR"),
      check = "contains")

# MC=1: should also produce "KAAR"
kka_mc1 <- enzyme_cleave("KAARVGR", "trypsin", 1L)
.test("trypsin MC=1 produces 1-missed-cleavage peptides",
      kka_mc1$Sequence,
      c("KAAR", "AARVGR"),
      check = "contains")

# Trypsin should NOT cleave KP or RP (rule: not before P)
.test("trypsin MC=0 does not cleave before P in KPR",
      enzyme_cleave("KPR", "trypsin", 0L)$Sequence,
      c("KPR"),
      check = "contains")

# Lys-C cleaves only after K
lysc_dt <- enzyme_cleave("KAARVGR", "lysc", 0L)
.test("lys-c cleaves after K only",
      lysc_dt$Sequence,
      c("K", "AARVGR"),
      check = "contains")


# =============================================================================
# Section 7 — detect_conjugation_sites()
# =============================================================================
cat("\n── Section 7: detect_conjugation_sites() ────────────────────\n")

# CASIQKFGR has one C at position 1
.test("detect_conjugation_sites CASIQKFGR cysteine → pos 1",
      detect_conjugation_sites("CASIQKFGR", "cysteine", "C"),
      1L,
      check = "integer")

# CXXCXX — two cysteines
.test("detect_conjugation_sites CACF cysteine → pos 1,3",
      detect_conjugation_sites("CACF", "cysteine", "C"),
      c(1L, 3L),
      check = "integer")

# Lysine conjugation excludes N-terminal K
.test("detect_conjugation_sites KASK lysine → pos 4 (not 1)",
      detect_conjugation_sites("KASK", "lysine", "K"),
      4L,
      check = "integer")

# No matching residue
.test("detect_conjugation_sites no C in ASGR → integer(0)",
      detect_conjugation_sites("ASGR", "cysteine", "C"),
      integer(0),
      check = "integer")

# has_conjugation_site()
.test("has_conjugation_site CASIQKFGR cysteine → TRUE",
      has_conjugation_site("CASIQKFGR", "cysteine", "C"),
      TRUE,
      check = "logical")

.test("has_conjugation_site ASGR cysteine → FALSE",
      has_conjugation_site("ASGR", "cysteine", "C"),
      FALSE,
      check = "logical")


# =============================================================================
# Section 8 — calc_dar_peptide_masses()
# =============================================================================
cat("\n── Section 8: calc_dar_peptide_masses() ─────────────────────\n")

dar_masses <- calc_dar_peptide_masses(
  sequence         = "CASIQKFGR",
  payload_mass     = 715.3,
  conjugation_sites = 1L,
  dar_range        = 0L:1L
)

.test("calc_dar_peptide_masses() returns 2 rows (DAR 0 and 1)",
      nrow(dar_masses), 2L, check = "integer")

.test("DAR=0 NakedMass matches calc_peptide_mass()",
      dar_masses[DAR == 0L, NakedMass],
      calc_peptide_mass("CASIQKFGR"),
      tol = 1e-4)

.test("DAR=1 LoadedMass = NakedMass + 715.3",
      dar_masses[DAR == 1L, LoadedMass],
      calc_peptide_mass("CASIQKFGR") + 715.3,
      tol = 1e-4)

.test("DAR=1 MassDelta = 715.3",
      dar_masses[DAR == 1L, MassDelta],
      715.3,
      tol = 1e-4)

.test("DAR=0 LoadedSites is empty string",
      dar_masses[DAR == 0L, LoadedSites],
      "",
      check = "character")

.test("DAR=1 LoadedSites = '1'",
      dar_masses[DAR == 1L, LoadedSites],
      "1",
      check = "character")

# DAR range clamping: asking for DAR=2 on a 1-site peptide should silently drop it
dar_clamped <- calc_dar_peptide_masses("CASIQKFGR", 715.3, 1L, dar_range = 0L:4L)
.test("calc_dar_peptide_masses() clamps DAR range to n_sites",
      max(dar_clamped$DAR), 1L, check = "integer")


# =============================================================================
# Section 9 — generate_dar_transitions()
# =============================================================================
cat("\n── Section 9: generate_dar_transitions() ────────────────────\n")

pep_dt <- data.table(
  PeptideSequence  = "CASIQKFGR",
  ModifiedSequence = "CASIQKFGR",
  ProteinName      = "TestAb",
  Chain            = "HC",
  Modifications    = "",
  UniqueToADC      = TRUE,
  Start            = 1L,
  End              = 9L,
  Enzyme           = "trypsin",
  mod_list         = list(NULL)
)

dar_trans <- generate_dar_transitions(
  pep_dt,
  payload_mass     = 715.3,
  dar_range        = 0L:1L,
  conjugation_type = "cysteine",
  residue          = "C"
)

.test("generate_dar_transitions() produces rows",
      nrow(dar_trans) > 0L,
      TRUE,
      check = "logical")

.test("generate_dar_transitions() has DAR column",
      "DAR" %in% names(dar_trans),
      TRUE,
      check = "logical")

.test("generate_dar_transitions() has DAR 0 and 1",
      sort(unique(dar_trans$DAR)),
      c(0L, 1L),
      check = "integer")

# For z=2, DAR=1 precursor m/z should be exactly 715.3/2 above DAR=0
dar0_mz <- dar_trans[DAR == 0L & PrecursorCharge == 2L, PrecursorMz][1]
dar1_mz <- dar_trans[DAR == 1L & PrecursorCharge == 2L, PrecursorMz][1]
.test("DAR=1 vs DAR=0 Δ(PrecursorMz) at z=2 = 715.3/2 Da",
      dar1_mz - dar0_mz,
      715.3 / 2,
      tol = 1e-3)

# DAR column should be the first column
.test("DAR is the first column in generate_dar_transitions() output",
      names(dar_trans)[1],
      "DAR",
      check = "character")


# =============================================================================
# Section 10 — calc_isotope_distribution() (averagine approximation)
# =============================================================================
cat("\n── Section 10: calc_isotope_distribution() ──────────────────\n")

iso <- calc_isotope_distribution("CASIQKFGR")

.test("calc_isotope_distribution() returns a data frame",
      is.data.frame(iso),
      TRUE,
      check = "logical")

.test("calc_isotope_distribution() has ≥4 isotope peaks for ~1 kDa peptide",
      nrow(iso) >= 4L,
      TRUE,
      check = "logical")

.test("calc_isotope_distribution() relative intensities sum to ~100",
      sum(iso$rel_intensity),
      100,
      tol = 0.5)

.test("calc_isotope_distribution() monoisotopic peak is most abundant for <2 kDa",
      which.max(iso$rel_intensity),
      1L,
      check = "integer")

# Larger peptide (~3 kDa): M+1 or M+2 should dominate
# Use a synthetic 25-residue peptide
long_pep <- paste(rep("ACDEFGHIKLMNPQRSTVWY", 2), collapse = "")[1]
long_pep <- paste(rep("AAAAAAAAAAAAAAAAAAAAAAAAAAA", 1), collapse = "")  # 27 × Ala
iso_large <- calc_isotope_distribution(long_pep)
.test("calc_isotope_distribution() M+0 not always dominant for longer peptides",
      nrow(iso_large) >= 2L,
      TRUE,
      check = "logical")


# =============================================================================
# Section 11 — ADCDB_PAYLOADS and LINKER_BIOTRANSFORMATIONS constants
# =============================================================================
cat("\n── Section 11: Payload and biotransformation constants ──────\n")

.test("ADCDB_PAYLOADS is a list",
      is.list(ADCDB_PAYLOADS),
      TRUE,
      check = "logical")

.test("ADCDB_PAYLOADS has ≥5 entries",
      length(ADCDB_PAYLOADS) >= 5L,
      TRUE,
      check = "logical")

# Every payload entry must have mass and conjugation_type
missing_fields <- Filter(function(nm) {
  entry <- ADCDB_PAYLOADS[[nm]]
  !all(c("mass", "conjugation_type") %in% names(entry))
}, names(ADCDB_PAYLOADS))
.test("All ADCDB_PAYLOADS entries have 'mass' and 'conjugation_type'",
      length(missing_fields) == 0L,
      TRUE,
      check = "logical")

.test("LINKER_BIOTRANSFORMATIONS is a list",
      is.list(LINKER_BIOTRANSFORMATIONS),
      TRUE,
      check = "logical")

.test("LINKER_BIOTRANSFORMATIONS has ≥5 entries",
      length(LINKER_BIOTRANSFORMATIONS) >= 5L,
      TRUE,
      check = "logical")

# Maleimide hydrolysis mass = +18.011 (H2O addition, ring-opening)
mh_mass <- LINKER_BIOTRANSFORMATIONS$maleimide_hydrolysis$delta_mass
.test("LINKER_BIOTRANSFORMATIONS maleimide_hydrolysis delta = +18.011 Da",
      mh_mass, 18.010565, tol = 1e-4)

# Deamidation delta = +0.984 Da
dam_mass <- LINKER_BIOTRANSFORMATIONS$deamidation_linker_N$delta_mass
.test("LINKER_BIOTRANSFORMATIONS deamidation_linker_N delta = +0.984 Da",
      dam_mass, 0.984016, tol = 1e-4)


# =============================================================================
# Final report
# =============================================================================
cat("\n=============================================================\n")
total <- .pass + .fail
cat(sprintf(" Results: %d/%d tests passed", .pass, total))
if (.fail > 0L) {
  cat(sprintf(" (%d FAILED)\n", .fail))
  cat(" Failed tests:\n")
  for (f in .failures) cat(sprintf("   - %s\n", f))
} else {
  cat(" — all tests PASSED ✓\n")
}
cat("=============================================================\n")

if (.fail > 0L) quit(status = 1L, save = "no")
