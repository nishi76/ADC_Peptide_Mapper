# =============================================================================
# modifications.R — ADC Peptide Mapper v0.8
# =============================================================================
# Sourced by the Shiny app (app.R / server.R). No library() calls; base R only.
#
# Constants AA_MONO_MASS, WATER_MASS, and PROTON_MASS are defined in digest.R
# and are available in the shared environment when this file is sourced.
# Do NOT redefine them here.
#
# Contents:
#   1.  FIXED_MODS                  — always-on modifications (e.g. CAM on Cys)
#   2.  ADCDB_PAYLOADS              — curated ADC payload database (with conjugation_type)
#   3.  ADCDB_LINKERS               — curated ADC linker database (ADCDB-sourced)
#   3b. LINKER_BIOTRANSFORMATIONS   — in vivo/vitro linker biotransformation products
#   4.  VAR_MOD_DEFS                — variable modification definitions (incl. ADCDB payloads
#                                     and biotr_ entries from LINKER_BIOTRANSFORMATIONS)
#   5.  SPECIAL_MOD_DEFS            — special / non-residue modifications
#   6.  calc_modified_mass()        — compute peptide mass with active mods applied
#   7.  build_active_mods()         — assemble the active mod list from UI selections
#   8.  apply_modifications()       — apply mods to a peptide sequence, return mod table
#   9.  get_payload_choices()       — helper for Shiny selectInput() grouped choices
#  10.  detect_conjugation_sites()  — T3-B: find conjugation-capable residues in a peptide
#  11.  has_conjugation_site()      — T3-B: convenience wrapper
# =============================================================================


# -----------------------------------------------------------------------------
# 1. FIXED MODIFICATIONS
# -----------------------------------------------------------------------------
# Fixed mods are applied unconditionally to every matching residue.
# Each entry:
#   residue  — single-letter AA code the mod targets
#   name     — human-readable label
#   mass     — monoisotopic mass shift (Da) added to the residue
#   nterm    — logical; TRUE means only apply at peptide N-terminus
#   location — positional constraint: "any" | "nterm" | "cterm" | "<int>"
#              (built-in fixed mods default to "any")

FIXED_MODS <- list(
  CAM = list(
    residue  = "C",
    name     = "Carbamidomethyl (CAM)",
    mass     = 57.02146,   # C2H3NO — iodoacetamide alkylation of Cys thiol
    nterm    = FALSE,
    location = "any"
  )
)


# -----------------------------------------------------------------------------
# 2. ADCDB PAYLOAD DATABASE
# -----------------------------------------------------------------------------
# Curated from https://adcdb.idrblab.net — the ADC Database.
# Each entry represents the cytotoxic payload component of a clinical or
# approved ADC.  The `mass` field is the monoisotopic (or best-available)
# mass shift (Da) that the payload contributes at the conjugation site on
# the antibody/peptide.
#
# Fields:
#   name             — full descriptive name
#   abbrev           — short abbreviation used in labels / column headers
#   mass             — mass shift (Da) at conjugation site (monoisotopic where known)
#   residue          — amino acid targeted for conjugation ("C" = Cys, "K" = Lys)
#   conjugation_type — T3-B: one of "cysteine" | "lysine" | "site_specific"
#                        "cysteine"     → maleimide thiol-chemistry at interchain
#                                         disulfide Cys residues (random, DAR 0–8)
#                        "lysine"       → NHS-ester or hydrazone at solvent-exposed Lys
#                                         (statistical, DAR 0–8)
#                        "site_specific"→ engineered/unnatural amino acid or
#                                         enzymatic conjugation at a unique defined site
#   mechanism        — pharmacological mechanism of action
#   example_adc      — representative approved ADC using this payload

ADCDB_PAYLOADS <- list(

  # --- Topoisomerase I inhibitors -------------------------------------------
  dxd = list(
    name             = "DXd (Deruxtecan payload)",
    abbrev           = "DXd",
    mass             = 519.2,          # exatecan derivative; MW ~519 Da at conjugation
    residue          = "C",
    conjugation_type = "cysteine",
    mechanism        = "Topoisomerase I inhibitor",
    example_adc      = "Trastuzumab deruxtecan (Enhertu)"
  ),

  sn38 = list(
    name             = "SN-38 (Irinotecan metabolite)",
    abbrev           = "SN38",
    mass             = 392.1,          # C20H20N2O5; monoisotopic ~392.1 Da
    residue          = "C",
    conjugation_type = "cysteine",
    mechanism        = "Topoisomerase I inhibitor",
    example_adc      = "Sacituzumab govitecan (Trodelvy)"
  ),

  dxd_trop2 = list(
    name             = "DXd (TROP2-targeted)",
    abbrev           = "DXd-T2",
    mass             = 519.2,          # same DXd payload, different antibody target
    residue          = "C",
    conjugation_type = "cysteine",
    mechanism        = "Topoisomerase I inhibitor",
    example_adc      = "Datopotamab deruxtecan (Datroway)"
  ),

  # --- Auristatins (microtubule disruption) ----------------------------------
  mmae = list(
    name             = "MMAE (Monomethyl auristatin E)",
    abbrev           = "MMAE",
    mass             = 715.3,          # C39H67N5O7; approx monoisotopic 715.3 Da
    residue          = "C",
    conjugation_type = "cysteine",
    mechanism        = "Microtubule disruption (auristatin)",
    example_adc      = "Brentuximab vedotin (Adcetris)"
  ),

  mmaf = list(
    name             = "MMAF (Monomethyl auristatin F)",
    abbrev           = "MMAF",
    mass             = 731.3,          # phenylalanine C-terminus variant of MMAE
    residue          = "C",
    conjugation_type = "cysteine",
    mechanism        = "Microtubule disruption (auristatin)",
    example_adc      = "Belantamab mafodotin (Blenrep)"
  ),

  # --- Maytansinoids (microtubule disruption) --------------------------------
  dm1 = list(
    name             = "DM1 (Emtansine)",
    abbrev           = "DM1",
    mass             = 738.0,          # C35H48ClNO10S; average MW ~738 Da
    residue          = "C",
    conjugation_type = "cysteine",
    mechanism        = "Microtubule disruption (maytansinoid)",
    example_adc      = "Ado-trastuzumab emtansine (Kadcyla)"
  ),

  dm4 = list(
    name             = "DM4 (Ravtansine)",
    abbrev           = "DM4",
    mass             = 752.0,          # DM1 analogue with additional methyl group; ~752 Da
    residue          = "C",
    conjugation_type = "cysteine",
    mechanism        = "Microtubule disruption (maytansinoid)",
    example_adc      = "Mirvetuximab soravtansine (Elahere)"
  ),

  # --- DNA-targeting agents --------------------------------------------------
  pbd = list(
    name             = "PBD dimer (Pyrrolobenzodiazepine)",
    abbrev           = "PBD",
    mass             = 1000.4,         # SG2000-class PBD dimer; ~1000 Da
    residue          = "C",
    conjugation_type = "cysteine",
    mechanism        = "DNA crosslinker",
    example_adc      = "Loncastuximab tesirine (Zynlonta)"
  ),

  calicheamicin = list(
    name             = "Calicheamicin",
    abbrev           = "Calich",
    mass             = 1368.3,         # C55H74IN3O21S4; average MW ~1368 Da
    residue          = "K",            # lysine-conjugated via hydrazone linker
    conjugation_type = "lysine",
    mechanism        = "DNA double-strand break",
    example_adc      = "Gemtuzumab ozogamicin (Mylotarg)"
  ),

  # --- User-defined placeholder ---------------------------------------------
  custom_payload = list(
    name             = "Custom Payload",
    abbrev           = "CustomPL",
    mass             = NA_real_,       # user must supply mass at runtime
    residue          = "C",
    conjugation_type = "cysteine",
    mechanism        = "User-defined",
    example_adc      = "N/A"
  )
)


# -----------------------------------------------------------------------------
# 3. ADCDB LINKER DATABASE
# -----------------------------------------------------------------------------
# Curated from https://adcdb.idrblab.net.
# The `mass` field is the mass contribution (Da) of the LINKER PORTION ONLY,
# i.e. the mass added to the peptide/antibody from the linker chemistry,
# separate from the payload mass.  Total ADC conjugate mass shift =
#   linker$mass + payload$mass
#
# Fields:
#   name               — full descriptive name
#   abbrev             — short abbreviation
#   mass               — linker-only mass contribution (Da)
#   linker_type        — "cleavable" | "non-cleavable" | "custom"
#   cleavage_chemistry — enzyme/condition responsible for linker cleavage
#   example_adc        — representative ADC using this linker

ADCDB_LINKERS <- list(

  # --- Cathepsin-cleavable linkers ------------------------------------------
  mc_ggfg = list(
    name               = "Mc-Gly-Gly-Phe-Gly (cathepsin-cleavable)",
    abbrev             = "McGGFG",
    mass               = 529.2,   # maleimidocaproyl-GGFG tetrapeptide spacer
    linker_type        = "cleavable",
    cleavage_chemistry = "Cathepsin B/L",
    example_adc        = "Trastuzumab deruxtecan"
  ),

  mc_val_cit = list(
    name               = "Mc-Val-Cit-PABC (cathepsin-cleavable)",
    abbrev             = "McVC",
    mass               = 571.3,   # maleimidocaproyl-Val-Cit-PABC
    linker_type        = "cleavable",
    cleavage_chemistry = "Cathepsin B",
    example_adc        = "Brentuximab vedotin"
  ),

  az_pbd_linker = list(
    name               = "SG3249 (PBD linker, cleavable)",
    abbrev             = "SG3249",
    mass               = 612.3,
    linker_type        = "cleavable",
    cleavage_chemistry = "Cathepsin B",
    example_adc        = "Loncastuximab tesirine"
  ),

  # --- Disulfide / reductive cleavage ---------------------------------------
  spdb = list(
    name               = "SPDB (N-succinimidyl-4-(2-pyridyldithio)butanoate)",
    abbrev             = "SPDB",
    mass               = 278.1,   # disulfide-containing crosslinker
    linker_type        = "cleavable",
    cleavage_chemistry = "Disulfide (reductive)",
    example_adc        = "Mirvetuximab soravtansine"
  ),

  # --- pH-sensitive / acid-labile linkers -----------------------------------
  cl2a = list(
    name               = "CL2A (pH-sensitive carbonate linker)",
    abbrev             = "CL2A",
    mass               = 429.2,
    linker_type        = "cleavable",
    cleavage_chemistry = "pH-sensitive hydrolysis",
    example_adc        = "Sacituzumab govitecan"
  ),

  aca_dhz = list(
    name               = "AcBut-hydrazone (pH-sensitive)",
    abbrev             = "AcBut",
    mass               = 245.1,   # 4-(4-acetylphenoxy)butanoic acid hydrazone
    linker_type        = "cleavable",
    cleavage_chemistry = "Acid-labile hydrazone",
    example_adc        = "Gemtuzumab ozogamicin"
  ),

  # --- Non-cleavable linkers ------------------------------------------------
  smcc = list(
    name               = "SMCC (succinimidyl-4-(N-maleimidomethyl)cyclohexane-1-carboxylate)",
    abbrev             = "SMCC",
    mass               = 334.1,   # thioether bond; stable in circulation
    linker_type        = "non-cleavable",
    cleavage_chemistry = "None",
    example_adc        = "Ado-trastuzumab emtansine (Kadcyla)"
  ),

  sgg = list(
    name               = "SGG (maleimide-based non-cleavable)",
    abbrev             = "SGG",
    mass               = 356.1,
    linker_type        = "non-cleavable",
    cleavage_chemistry = "None",
    example_adc        = "Belantamab mafodotin"
  ),

  # --- User-defined placeholder ---------------------------------------------
  custom_linker = list(
    name               = "Custom Linker",
    abbrev             = "CustomLK",
    mass               = NA_real_, # user must supply mass at runtime
    linker_type        = "custom",
    cleavage_chemistry = "User-defined",
    example_adc        = "N/A"
  )
)


# -----------------------------------------------------------------------------
# 3b.  LINKER BIOTRANSFORMATION PRODUCTS  (T3-C)
# -----------------------------------------------------------------------------
# In vivo and in vitro biotransformations of ADC linkers produce predictable
# mass shifts that need to be in the search space for comprehensive ADME
# characterisation.  Ref: Shen et al. (2012) Nat. Biotechnol. 30:184-189;
# Alley et al. (2008) Bioconjug. Chem. 19:759-765.
#
# Each entry:
#   delta_mass  — monoisotopic mass change (Da) applied to the payload-bearing residue
#   description — human-readable label for UI display
#   applies_to  — character vector: residues or "any" (residue where payload sits)
#   ref         — literature reference

LINKER_BIOTRANSFORMATIONS <- list(

  maleimide_hydrolysis = list(
    delta_mass  = +18.010565,
    description = "Maleimide ring opening / hydrolysis (+H2O)",
    applies_to  = "C",
    ref         = "Shen et al. 2012 Nat. Biotechnol."
  ),

  succinimide_ring_open = list(
    delta_mass  = +18.010565,
    description = "Succinimide ring-opening isomerisation (+H2O)",
    applies_to  = "C",
    ref         = "Alley et al. 2008 Bioconjug. Chem."
  ),

  thioether_sulfoxide = list(
    delta_mass  = +15.994915,
    description = "Thioether oxidation → sulfoxide (+O)",
    applies_to  = "C",
    ref         = "Lyon et al. 2015 Nat. Biotechnol."
  ),

  disulfide_loss = list(
    delta_mass  = -31.990415,
    description = "Sulfone reduction / disulfide scrambling (−2S+O; approx.)",
    applies_to  = "C",
    ref         = "Tumey et al. 2014 J. Med. Chem."
  ),

  payload_release_topo = list(
    delta_mass  = -519.2,      # DXd-class; overridden at runtime by payload mass
    description = "Full linker-payload release (topoisomerase payloads; delta varies)",
    applies_to  = "any",
    ref         = "Ogitani et al. 2016 Clin. Cancer Res."
  ),

  deamidation_linker_N = list(
    delta_mass  = +0.984016,
    description = "Deamidation of Asn residue adjacent to conjugation site (+0.984 Da)",
    applies_to  = "N",
    ref         = "Harris et al. 2001 J. Chromatogr. B"
  )
)


# =============================================================================
# detect_conjugation_sites()  (T3-B)
# =============================================================================
# Identify positions within a peptide sequence that are candidate conjugation
# sites based on the payload's conjugation chemistry.
#
# For cysteine-conjugated payloads: returns positions of all Cys residues.
# For lysine-conjugated payloads:   returns positions of all Lys residues
#                                   (excluding N-terminal Lys, which is less
#                                   reactive in NHS-ester conjugation).
# For site_specific:                returns all matching residue positions
#                                   (user-supplied explicit position takes priority
#                                   at the app level).
#
# Args:
#   sequence         : character(1) — single-letter amino acid sequence (peptide)
#   conjugation_type : character(1) — "cysteine" | "lysine" | "site_specific"
#   residue          : character(1) — target amino acid (from ADCDB_PAYLOADS$residue)
#
# Returns:
#   integer vector — 1-based positions within the peptide where conjugation
#   could occur.  Returns integer(0) if no candidate sites found.
# =============================================================================
detect_conjugation_sites <- function(sequence, conjugation_type, residue = "C") {
  residues <- strsplit(toupper(trimws(sequence)), "", fixed = TRUE)[[1L]]
  n        <- length(residues)

  target_aa <- toupper(trimws(residue))
  positions  <- which(residues == target_aa)

  if (length(positions) == 0L) return(integer(0L))

  if (conjugation_type == "lysine") {
    # Exclude N-terminal Lys (position 1) — lower NHS-ester reactivity
    positions <- positions[positions > 1L]
  }

  as.integer(positions)
}


# =============================================================================
# has_conjugation_site()
# =============================================================================
# Convenience wrapper: TRUE iff the peptide contains at least one candidate
# conjugation site for the given payload.
# =============================================================================
has_conjugation_site <- function(sequence, conjugation_type, residue = "C") {
  length(detect_conjugation_sites(sequence, conjugation_type, residue)) > 0L
}


# -----------------------------------------------------------------------------
# 4. VARIABLE MODIFICATION DEFINITIONS  (VAR_MOD_DEFS)
# -----------------------------------------------------------------------------
# Variable mods are user-selectable; zero or more may be active at once.
# Each entry:
#   residue  — single-letter AA code targeted (NA for whole-peptide mods)
#   name     — display label
#   mass     — monoisotopic mass shift (Da)
#   nterm    — logical; TRUE = only apply at peptide N-terminus
#   location — positional constraint: "any" | "nterm" | "cterm" | "<int>"
#              (all built-in var mods default to "any" unless nterm=TRUE)
#
# ADCDB payload entries are generated programmatically from ADCDB_PAYLOADS
# so that the two databases stay in sync automatically.

VAR_MOD_DEFS <- list(

  # --- Common chemical modifications ----------------------------------------

  oxidation = list(
    residue  = "M",
    name     = "Oxidation (Met)",
    mass     = 15.99491,   # +O; methionine sulfoxide
    nterm    = FALSE,
    location = "any"
  ),

  propionamide = list(
    residue  = "C",
    name     = "Propionamide (Cys)",
    mass     = 71.03711,   # acrylamide adduct on Cys; C3H5NO
    nterm    = FALSE,
    location = "any"
  ),

  nem = list(
    residue  = "C",
    name     = "N-Ethylmaleimide (NEM, Cys)",
    mass     = 125.04768,  # C6H7NO2; maleimide adduct on free Cys
    nterm    = FALSE,
    location = "any"
  ),

  # --- Legacy / pre-ADCDB payload entries (kept for backward compatibility) --
  # These mirror entries in ADCDB_PAYLOADS but are retained here so that
  # existing saved sessions / UI state strings continue to resolve correctly.

  mmae = list(
    residue  = "C",
    name     = "MMAE conjugation (Cys)",
    mass     = ADCDB_PAYLOADS$mmae$mass,   # 715.3 Da — references ADCDB_PAYLOADS
    nterm    = FALSE,
    location = "any"
  ),

  dm1 = list(
    residue  = "C",
    name     = "DM1 conjugation (Cys)",
    mass     = ADCDB_PAYLOADS$dm1$mass,    # 738.0 Da
    nterm    = FALSE,
    location = "any"
  ),

  dxd = list(
    residue  = "C",
    name     = "DXd conjugation (Cys)",
    mass     = ADCDB_PAYLOADS$dxd$mass,    # 519.2 Da
    nterm    = FALSE,
    location = "any"
  ),

  sn38 = list(
    residue  = "C",
    name     = "SN-38 conjugation (Cys)",
    mass     = ADCDB_PAYLOADS$sn38$mass,   # 392.1 Da
    nterm    = FALSE,
    location = "any"
  ),

  custom_drug = list(
    residue  = "C",
    name     = "Custom Drug (Cys, user mass)",
    mass     = NA_real_,   # resolved at runtime from UI input
    nterm    = FALSE,
    location = "any"
  ),

  # --- ADCDB payload entries (auto-generated block) -------------------------
  # Each ADCDB_PAYLOADS entry is mirrored here so it appears in the
  # variable-mod selector.  Keys are prefixed "adcdb_" to avoid collisions
  # with the legacy entries above.

  adcdb_dxd = list(
    residue  = ADCDB_PAYLOADS$dxd$residue,
    name     = paste0("[ADCDB] ", ADCDB_PAYLOADS$dxd$name),
    mass     = ADCDB_PAYLOADS$dxd$mass,
    nterm    = FALSE,
    location = "any"
  ),

  adcdb_mmae = list(
    residue  = ADCDB_PAYLOADS$mmae$residue,
    name     = paste0("[ADCDB] ", ADCDB_PAYLOADS$mmae$name),
    mass     = ADCDB_PAYLOADS$mmae$mass,
    nterm    = FALSE,
    location = "any"
  ),

  adcdb_mmaf = list(
    residue  = ADCDB_PAYLOADS$mmaf$residue,
    name     = paste0("[ADCDB] ", ADCDB_PAYLOADS$mmaf$name),
    mass     = ADCDB_PAYLOADS$mmaf$mass,
    nterm    = FALSE,
    location = "any"
  ),

  adcdb_dm1 = list(
    residue  = ADCDB_PAYLOADS$dm1$residue,
    name     = paste0("[ADCDB] ", ADCDB_PAYLOADS$dm1$name),
    mass     = ADCDB_PAYLOADS$dm1$mass,
    nterm    = FALSE,
    location = "any"
  ),

  adcdb_dm4 = list(
    residue  = ADCDB_PAYLOADS$dm4$residue,
    name     = paste0("[ADCDB] ", ADCDB_PAYLOADS$dm4$name),
    mass     = ADCDB_PAYLOADS$dm4$mass,
    nterm    = FALSE,
    location = "any"
  ),

  adcdb_sn38 = list(
    residue  = ADCDB_PAYLOADS$sn38$residue,
    name     = paste0("[ADCDB] ", ADCDB_PAYLOADS$sn38$name),
    mass     = ADCDB_PAYLOADS$sn38$mass,
    nterm    = FALSE,
    location = "any"
  ),

  adcdb_pbd = list(
    residue  = ADCDB_PAYLOADS$pbd$residue,
    name     = paste0("[ADCDB] ", ADCDB_PAYLOADS$pbd$name),
    mass     = ADCDB_PAYLOADS$pbd$mass,
    nterm    = FALSE,
    location = "any"
  ),

  adcdb_calicheamicin = list(
    residue  = ADCDB_PAYLOADS$calicheamicin$residue,   # "K"
    name     = paste0("[ADCDB] ", ADCDB_PAYLOADS$calicheamicin$name),
    mass     = ADCDB_PAYLOADS$calicheamicin$mass,
    nterm    = FALSE,
    location = "any"
  ),

  adcdb_dxd_trop2 = list(
    residue  = ADCDB_PAYLOADS$dxd_trop2$residue,
    name     = paste0("[ADCDB] ", ADCDB_PAYLOADS$dxd_trop2$name),
    mass     = ADCDB_PAYLOADS$dxd_trop2$mass,
    nterm    = FALSE,
    location = "any"
  ),

  adcdb_custom_payload = list(
    residue  = ADCDB_PAYLOADS$custom_payload$residue,
    name     = paste0("[ADCDB] ", ADCDB_PAYLOADS$custom_payload$name),
    mass     = NA_real_,   # user must supply at runtime
    nterm    = FALSE,
    location = "any"
  )
)


# -----------------------------------------------------------------------------
# 5. SPECIAL MODIFICATION DEFINITIONS  (SPECIAL_MOD_DEFS)
# -----------------------------------------------------------------------------
# Special mods cover whole-peptide or terminus-specific chemistry that does
# not map cleanly to a single residue.  Kept unchanged from v0.5.
#
# Each entry:
#   residue  — NA (whole-peptide) or single-letter AA
#   name     — display label
#   mass     — monoisotopic mass shift (Da)
#   nterm    — logical; TRUE = N-terminal modification
#   location — positional constraint (default "any"; "nterm" for N-term mods)

SPECIAL_MOD_DEFS <- list(

  # Pyro-glutamate from N-terminal glutamine (loss of NH3)
  pyro_glu_Q = list(
    residue  = "Q",
    name     = "Pyro-Glu from Q (N-term)",
    mass     = -17.02655,  # -NH3
    nterm    = TRUE,
    location = "nterm"
  ),

  # Pyro-glutamate from N-terminal glutamic acid (loss of H2O)
  pyro_glu_E = list(
    residue  = "E",
    name     = "Pyro-Glu from E (N-term)",
    mass     = -18.01056,  # -H2O
    nterm    = TRUE,
    location = "nterm"
  ),

  # N-terminal acetylation
  acetyl_nterm = list(
    residue  = NA_character_,  # applies to any N-terminal residue
    name     = "Acetylation (N-term)",
    mass     = 42.01057,   # C2H2O
    nterm    = TRUE,
    location = "nterm"
  ),

  # Deamidation of asparagine
  deamidation_N = list(
    residue  = "N",
    name     = "Deamidation (Asn)",
    mass     = 0.98402,    # N->D conversion; +0.984 Da
    nterm    = FALSE,
    location = "any"
  ),

  # Deamidation of glutamine
  deamidation_Q = list(
    residue  = "Q",
    name     = "Deamidation (Gln)",
    mass     = 0.98402,    # Q->E conversion
    nterm    = FALSE,
    location = "any"
  ),

  # Phosphorylation (Ser/Thr/Tyr)
  # Note: residue is NA here; apply_modifications() handles multi-residue
  # targeting by checking aas against c("S","T","Y") when this mod is active.
  # Callers that need per-residue specificity should use separate entries.
  phospho_STY = list(
    residue  = NA_character_,  # applies to S, T, or Y — see apply_modifications()
    name     = "Phosphorylation (Ser/Thr/Tyr)",
    mass     = 79.96633,   # HPO3
    nterm    = FALSE,
    location = "any"
  ),

  # Methyl ester (C-terminal, e.g. from methanol in sample prep)
  methyl_ester_cterm = list(
    residue  = NA_character_,
    name     = "Methyl ester (C-term)",
    mass     = 14.01565,   # +CH2
    nterm    = FALSE,
    location = "cterm"
  ),

  # Amidation of C-terminus
  amidation_cterm = list(
    residue  = NA_character_,
    name     = "Amidation (C-term)",
    mass     = -0.98402,   # -OH + NH2 = -0.984 Da
    nterm    = FALSE,
    location = "cterm"
  ),

  # ── T3-C: Linker biotransformation products ──────────────────────────────
  # Generated from LINKER_BIOTRANSFORMATIONS; each appears in the mod UI
  # when the Biotransformation Search panel is activated.
  biotr_maleimide_hydrolysis = list(
    residue  = "C",
    name     = "[Biotr] Maleimide ring opening (+18.011 Da)",
    mass     = +18.010565,
    nterm    = FALSE,
    location = "any"
  ),

  biotr_succinimide_open = list(
    residue  = "C",
    name     = "[Biotr] Succinimide ring-opening (+18.011 Da)",
    mass     = +18.010565,
    nterm    = FALSE,
    location = "any"
  ),

  biotr_thioether_sulfoxide = list(
    residue  = "C",
    name     = "[Biotr] Thioether → sulfoxide (+15.995 Da)",
    mass     = +15.994915,
    nterm    = FALSE,
    location = "any"
  ),

  biotr_disulfide_loss = list(
    residue  = "C",
    name     = "[Biotr] Disulfide loss / sulfone (−31.990 Da)",
    mass     = -31.990415,
    nterm    = FALSE,
    location = "any"
  ),

  biotr_deamidation_linker = list(
    residue  = "N",
    name     = "[Biotr] Deamidation near conjugation site (+0.984 Da)",
    mass     = +0.984016,
    nterm    = FALSE,
    location = "any"
  )
)


# =============================================================================
# FUNCTIONS
# =============================================================================


# -----------------------------------------------------------------------------
# .filter_positions_by_location()  [internal helper]
# -----------------------------------------------------------------------------
# Filter a vector of candidate residue positions according to the `location`
# constraint and the legacy `nterm` flag.
#
# Arguments:
#   positions  — integer vector of candidate positions (1-indexed)
#   location   — character(1): "any" | "nterm" | "cterm" | "<positive int>"
#   nterm_flag — logical(1): legacy N-terminus flag (TRUE overrides location)
#   n_res      — integer(1): total number of residues in the peptide
#
# Returns:
#   integer vector of positions that pass the filter.

.filter_positions_by_location <- function(positions, location, nterm_flag, n_res) {

  # Legacy nterm flag takes precedence (backward compatibility with v0.5)
  if (isTRUE(nterm_flag)) {
    return(positions[positions == 1L])
  }

  # Normalise: NULL or missing location treated as "any"
  if (is.null(location) || length(location) == 0L || is.na(location)) {
    return(positions)
  }

  if (location == "any") {
    return(positions)
  }

  if (location == "nterm") {
    return(positions[positions == 1L])
  }

  if (location == "cterm") {
    return(positions[positions == n_res])
  }

  # Numeric position string (e.g. "5" means only position 5)
  pos_int <- suppressWarnings(as.integer(location))
  if (!is.na(pos_int) && pos_int >= 1L) {
    return(positions[positions == pos_int])
  }

  # Unrecognised location value — warn and fall back to "any" (safe default)
  warning(sprintf(
    ".filter_positions_by_location: unrecognised location '%s'; defaulting to 'any'",
    location
  ))
  positions
}


# -----------------------------------------------------------------------------
# calc_modified_mass()
# -----------------------------------------------------------------------------
# Compute the monoisotopic neutral mass of a peptide given its unmodified
# residue masses and a list of active modifications.
#
# Arguments:
#   peptide_seq   — character(1); single-letter AA sequence string
#   active_mods   — list of mod entries (each with $residue, $mass, $nterm,
#                   $location); typically the output of build_active_mods()
#   base_mass     — numeric(1); pre-computed unmodified peptide mass (Da).
#                   If NA (default), mass is computed from peptide_seq using
#                   AA_MONO_MASS + WATER_MASS (defined in digest.R).
#
# Returns:
#   numeric(1) — monoisotopic neutral mass (Da) of the modified peptide.
#                Returns NA if any required mass is NA (e.g. custom_drug with
#                no user-supplied value) — the caller should handle this case.

calc_modified_mass <- function(peptide_seq, active_mods, base_mass = NA_real_) {

  # --- Compute base (unmodified) mass if not supplied -----------------------
  if (is.na(base_mass)) {
    aas <- strsplit(peptide_seq, "")[[1]]
    # Sum residue masses; fall back to 0 for unknown AAs (with a warning)
    residue_masses <- vapply(aas, function(aa) {
      m <- AA_MONO_MASS[[aa]]
      if (is.null(m)) {
        warning(sprintf(
          "calc_modified_mass: unknown amino acid '%s' in sequence '%s' — mass set to 0",
          aa, peptide_seq
        ))
        0
      } else {
        m
      }
    }, numeric(1))
    base_mass <- sum(residue_masses) + WATER_MASS
  }

  if (length(active_mods) == 0L) return(base_mass)

  aas      <- strsplit(peptide_seq, "")[[1]]
  n_res    <- length(aas)
  mod_mass <- 0

  for (mod in active_mods) {
    mod_residue  <- mod$residue
    mod_mass_val <- mod$mass
    mod_nterm    <- if (!is.null(mod$nterm))    mod$nterm    else FALSE
    mod_location <- if (!is.null(mod$location)) mod$location else "any"

    # Skip mods with NA mass (unresolved custom entries)
    if (is.na(mod_mass_val)) next

    # Determine which positions this mod applies to
    if (!is.na(mod_residue)) {
      # Residue-specific mod: find all matching positions
      hit_positions <- which(aas == mod_residue)
    } else {
      # Whole-peptide / terminus mod: all positions are candidates
      # (the location filter below will narrow this down as needed)
      hit_positions <- seq_len(n_res)
    }

    # Apply positional filter
    hit_positions <- .filter_positions_by_location(
      positions  = hit_positions,
      location   = mod_location,
      nterm_flag = mod_nterm,
      n_res      = n_res
    )

    # Accumulate mass shift for each qualifying position
    mod_mass <- mod_mass + length(hit_positions) * mod_mass_val
  }

  base_mass + mod_mass
}


# -----------------------------------------------------------------------------
# build_active_mods()
# -----------------------------------------------------------------------------
# Assemble the list of active modifications from:
#   (a) always-on FIXED_MODS
#   (b) user-selected variable mods (keys into VAR_MOD_DEFS)
#   (c) user-selected special mods  (keys into SPECIAL_MOD_DEFS)
#   (d) user-supplied custom mods   (data.frame with columns residue, name,
#                                    mass, and optionally location)
#
# Arguments:
#   selected_var_mods     — character vector of VAR_MOD_DEFS keys to activate.
#                           Pass character(0) or NULL for none.
#   selected_special_mods — character vector of SPECIAL_MOD_DEFS keys to activate.
#                           Pass character(0) or NULL for none.
#   custom_mods_df        — data.frame with columns:
#                             residue  (character) — single-letter AA code
#                             name     (character) — display label
#                             mass     (numeric)   — mass shift in Da
#                             location (character, OPTIONAL) — positional
#                               constraint; defaults to "any" if column absent
#                           Pass NULL or a zero-row data.frame for none.
#   custom_drug_mass      — numeric(1); mass to substitute for NA entries in
#                           VAR_MOD_DEFS (e.g. custom_drug, adcdb_custom_payload).
#                           Ignored if NA (default).
#   include_fixed         — logical(1); if FALSE, FIXED_MODS are excluded.
#                           Useful for computing theoretical unmodified masses.
#
# Returns:
#   Named list of mod entries ready for use by calc_modified_mass() and
#   apply_modifications().  List names are prefixed by source:
#     "fixed_*"   — from FIXED_MODS
#     "var_*"     — from VAR_MOD_DEFS
#     "special_*" — from SPECIAL_MOD_DEFS
#     "custom_*"  — from custom_mods_df (indexed by row number)

build_active_mods <- function(
    selected_var_mods     = character(0),
    selected_special_mods = character(0),
    custom_mods_df        = NULL,
    custom_drug_mass      = NA_real_,
    include_fixed         = TRUE
) {

  active <- list()

  # (a) Fixed modifications --------------------------------------------------
  if (isTRUE(include_fixed)) {
    for (key in names(FIXED_MODS)) {
      active[[paste0("fixed_", key)]] <- FIXED_MODS[[key]]
    }
  }

  # (b) Variable modifications -----------------------------------------------
  if (!is.null(selected_var_mods) && length(selected_var_mods) > 0L) {
    for (key in selected_var_mods) {
      mod <- VAR_MOD_DEFS[[key]]
      if (is.null(mod)) {
        warning(sprintf(
          "build_active_mods: unknown variable mod key '%s' — skipped", key
        ))
        next
      }
      # Resolve NA mass for custom/user-defined entries
      if (is.na(mod$mass) && !is.na(custom_drug_mass)) {
        mod$mass <- as.numeric(custom_drug_mass)
      }
      active[[paste0("var_", key)]] <- mod
    }
  }

  # (c) Special modifications ------------------------------------------------
  if (!is.null(selected_special_mods) && length(selected_special_mods) > 0L) {
    for (key in selected_special_mods) {
      mod <- SPECIAL_MOD_DEFS[[key]]
      if (is.null(mod)) {
        warning(sprintf(
          "build_active_mods: unknown special mod key '%s' — skipped", key
        ))
        next
      }
      active[[paste0("special_", key)]] <- mod
    }
  }

  # (d) Custom modifications from data.frame ---------------------------------
  if (!is.null(custom_mods_df) && is.data.frame(custom_mods_df) &&
      nrow(custom_mods_df) > 0L) {

    # Backward compatibility: add location column if absent (v0.5 sessions)
    if (!"location" %in% colnames(custom_mods_df)) {
      custom_mods_df$location <- "any"
    }

    # Coerce NA location values to "any"
    custom_mods_df$location[is.na(custom_mods_df$location)] <- "any"

    for (i in seq_len(nrow(custom_mods_df))) {
      row     <- custom_mods_df[i, , drop = FALSE]
      mod_key <- paste0("custom_", i)
      active[[mod_key]] <- list(
        residue  = as.character(row$residue),
        name     = as.character(row$name),
        mass     = as.numeric(row$mass),
        nterm    = FALSE,          # custom mods use location field exclusively
        location = as.character(row$location)
      )
    }
  }

  active
}


# -----------------------------------------------------------------------------
# apply_modifications()
# -----------------------------------------------------------------------------
# Apply a set of active modifications to a peptide sequence and return a
# data.frame describing every modification event (one row per modified site).
#
# This function is used to:
#   * annotate peptide fragments with their modification sites
#   * compute per-site mass shifts for display in the results table
#   * support the location-aware filtering introduced in v0.6
#
# Special handling:
#   * phospho_STY (residue = NA): applies to S, T, and Y residues.
#     When a mod has residue = NA and location = "any", ALL positions are
#     candidates (terminus mods with residue = NA are handled via location).
#   * Mods with NA mass are silently skipped (unresolved custom entries).
#
# Arguments:
#   peptide_seq  — character(1); single-letter AA sequence
#   active_mods  — named list from build_active_mods()
#
# Returns:
#   data.frame with columns:
#     position    — integer; 1-indexed position in peptide_seq
#     residue     — character; amino acid at that position
#     mod_key     — character; key from active_mods list
#     mod_name    — character; human-readable modification name
#     mass_shift  — numeric; mass shift (Da) applied at this position
#     location    — character; location constraint that was applied
#
#   Returns a zero-row data.frame (with the same columns) if no mods apply.

apply_modifications <- function(peptide_seq, active_mods) {

  # Pre-split sequence once for efficiency
  aas   <- strsplit(peptide_seq, "")[[1]]
  n_res <- length(aas)

  # Output accumulator (list of single-row lists, assembled at end)
  result_rows <- vector("list", length = 0L)

  for (mod_key in names(active_mods)) {
    mod <- active_mods[[mod_key]]

    mod_residue  <- mod$residue
    mod_mass_val <- mod$mass
    mod_nterm    <- if (!is.null(mod$nterm))    mod$nterm    else FALSE
    mod_location <- if (!is.null(mod$location)) mod$location else "any"
    mod_name     <- if (!is.null(mod$name))     mod$name     else mod_key

    # Skip mods with NA mass — cannot annotate without a mass value
    if (is.na(mod_mass_val)) next

    # Determine candidate positions ----------------------------------------
    if (!is.na(mod_residue)) {
      # Standard residue-specific mod
      candidate_positions <- which(aas == mod_residue)
    } else {
      # Whole-peptide / terminus mod (residue = NA).
      # All positions are candidates; location filter narrows them down.
      # Exception: phospho_STY targets S, T, Y specifically — detected by name.
      if (grepl("phospho", mod_key, ignore.case = TRUE) ||
          grepl("Phosphorylation", mod_name, ignore.case = TRUE)) {
        candidate_positions <- which(aas %in% c("S", "T", "Y"))
      } else {
        candidate_positions <- seq_len(n_res)
      }
    }

    # Apply location / nterm filter ----------------------------------------
    hit_positions <- .filter_positions_by_location(
      positions  = candidate_positions,
      location   = mod_location,
      nterm_flag = mod_nterm,
      n_res      = n_res
    )

    if (length(hit_positions) == 0L) next

    # Build one row per hit position ----------------------------------------
    for (pos in hit_positions) {
      result_rows[[length(result_rows) + 1L]] <- list(
        position   = pos,
        residue    = aas[pos],
        mod_key    = mod_key,
        mod_name   = mod_name,
        mass_shift = mod_mass_val,
        location   = mod_location
      )
    }
  }

  # Assemble output data.frame -----------------------------------------------
  if (length(result_rows) == 0L) {
    return(data.frame(
      position   = integer(0),
      residue    = character(0),
      mod_key    = character(0),
      mod_name   = character(0),
      mass_shift = numeric(0),
      location   = character(0),
      stringsAsFactors = FALSE
    ))
  }

  out <- data.frame(
    position   = vapply(result_rows, `[[`, integer(1),   "position"),
    residue    = vapply(result_rows, `[[`, character(1), "residue"),
    mod_key    = vapply(result_rows, `[[`, character(1), "mod_key"),
    mod_name   = vapply(result_rows, `[[`, character(1), "mod_name"),
    mass_shift = vapply(result_rows, `[[`, numeric(1),   "mass_shift"),
    location   = vapply(result_rows, `[[`, character(1), "location"),
    stringsAsFactors = FALSE
  )

  # Sort by position for readability in the UI table
  out[order(out$position), ]
}


# -----------------------------------------------------------------------------
# get_payload_choices()
# -----------------------------------------------------------------------------
# Returns a named list of named character vectors suitable for use as the
# `choices` argument of shiny::selectInput() / shiny::selectizeInput().
#
# Shiny renders named lists as optgroup-grouped dropdowns:
#   list("Group A" = c("Label 1" = "val1", "Label 2" = "val2"),
#        "Group B" = c("Label 3" = "val3"))
#
# Groups returned:
#   "Standard Payloads" — legacy VAR_MOD_DEFS entries present before ADCDB
#                         integration (mmae, dm1, dxd, sn38, custom_drug)
#   "ADCDB Payloads"    — all adcdb_* entries in VAR_MOD_DEFS, with enriched
#                         labels drawn from ADCDB_PAYLOADS metadata
#
# Usage in UI:
#   selectInput("payload_sel", "Select Payload",
#               choices = get_payload_choices())
#
# Returns:
#   Named list of named character vectors.
#   Inner vector: names = display labels, values = VAR_MOD_DEFS keys.

get_payload_choices <- function() {

  # --- Helper: build a display label for a VAR_MOD_DEFS entry ---------------
  .var_label <- function(key) {
    mod      <- VAR_MOD_DEFS[[key]]
    if (is.null(mod)) return(NA_character_)
    mass_str <- if (is.na(mod$mass)) "user-defined" else sprintf("%.1f Da", mod$mass)
    sprintf("%s [%s]", mod$name, mass_str)
  }

  # --- Helper: build an enriched label for an ADCDB payload entry -----------
  .adcdb_label <- function(key) {
    mod         <- VAR_MOD_DEFS[[key]]
    if (is.null(mod)) return(NA_character_)
    payload_key <- sub("^adcdb_", "", key)
    pl          <- ADCDB_PAYLOADS[[payload_key]]
    mass_str    <- if (is.na(mod$mass)) "user-defined" else sprintf("%.1f Da", mod$mass)
    if (!is.null(pl)) {
      # Include mechanism for richer context in the dropdown
      sprintf("%s [%s | %s]", pl$name, mass_str, pl$mechanism)
    } else {
      sprintf("%s [%s]", mod$name, mass_str)
    }
  }

  # --- Standard (legacy) payload keys ---------------------------------------
  standard_keys <- c("mmae", "dm1", "dxd", "sn38", "custom_drug")
  standard_labels <- vapply(standard_keys, .var_label, character(1))

  # Drop any keys that failed to resolve (defensive)
  valid_std <- !is.na(standard_labels)
  standard_choices <- setNames(standard_keys[valid_std], standard_labels[valid_std])

  # --- ADCDB payload keys (prefixed "adcdb_" in VAR_MOD_DEFS) ---------------
  adcdb_keys   <- grep("^adcdb_", names(VAR_MOD_DEFS), value = TRUE)
  adcdb_labels <- vapply(adcdb_keys, .adcdb_label, character(1))

  valid_adc <- !is.na(adcdb_labels)
  adcdb_choices <- setNames(adcdb_keys[valid_adc], adcdb_labels[valid_adc])

  # --- Return grouped list --------------------------------------------------
  list(
    "Standard Payloads" = standard_choices,
    "ADCDB Payloads"    = adcdb_choices
  )
}


# =============================================================================
# END OF modifications.R  (ADC Peptide Mapper v0.8)
# =============================================================================
