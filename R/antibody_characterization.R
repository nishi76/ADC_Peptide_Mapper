# antibody_characterization.R  —  ADC Peptide Mapper v1.0
# Binding affinity scoring, epitope mapping, FcRn / half-life prediction

# ── Binding affinity scoring ──────────────────────────────────────────────────

score_binding_affinity <- function(kd_nm, kon = NA, koff = NA) {
  tier <- if (kd_nm < 0.01) "Exceptional (picomolar)"
          else if (kd_nm < 1)   "Excellent"
          else if (kd_nm < 10)  "Good"
          else if (kd_nm < 100) "Borderline"
          else                  "Poor"

  badge_class <- if (kd_nm < 1)   "excellent"
                 else if (kd_nm < 10)  "good"
                 else if (kd_nm < 100) "borderline"
                 else                  "poor"

  adc_note <- if (kd_nm < 0.01)
    "Picomolar affinity may cause a 'binding site barrier': ADC saturates perivascular antigen before penetrating the tumor mass. Consider affinity maturation downward for solid-tumor ADCs."
  else if (kd_nm < 1)
    "Excellent affinity for ADC use. Strong tumor retention and slow off-rate support payload delivery."
  else if (kd_nm < 10)
    "Good affinity range for ADC applications. Balance between tumor penetration and retention is favorable."
  else if (kd_nm < 100)
    "Borderline affinity. Ensure antigen density is high enough to compensate. Consider affinity optimization."
  else
    "KD > 100 nM is typically insufficient for ADC applications. Affinity maturation strongly recommended."

  optimal_range <- "1–10 nM KD is generally optimal for solid tumor ADCs. Hematologic malignancy ADCs can tolerate 0.1–1 nM with lower binding-site barrier risk."

  derived_koff <- if (!is.na(kon) && !is.na(kd_nm)) kd_nm * 1e-9 * kon else koff

  list(
    kd_nm        = kd_nm,
    tier         = tier,
    badge_class  = badge_class,
    adc_note     = adc_note,
    optimal_range = optimal_range,
    kon          = kon,
    koff         = derived_koff,
    residence_time_min = if (!is.na(derived_koff) && derived_koff > 0)
      round(1 / derived_koff / 60, 1) else NA
  )
}

# ── Epitope characterization ──────────────────────────────────────────────────

characterize_epitope <- function(epitope_region, epitope_type = "Unknown",
                                  epitope_length = NA, conformation_sensitive = FALSE) {
  type_clean <- tolower(trimws(epitope_type))

  is_linear        <- grepl("linear|sequential", type_clean)
  is_conformational <- grepl("conform|discontinu|structural", type_clean)
  is_glycan         <- grepl("glycan|carbohydrate|sugar", type_clean)

  if (!is_linear && !is_conformational && !is_glycan) {
    is_linear <- !conformation_sensitive
    is_conformational <- conformation_sensitive
  }

  accessibility_note <- if (is_linear)
    "Linear epitopes are generally accessible when the antigen is surface-exposed. Suitable for peptide-blocking assays and ELISA competition."
  else if (is_conformational)
    "Conformational epitopes require the native folded protein. Denaturation (SDS-PAGE, denaturing ELISA) will abolish binding. Ensure assay conditions preserve antigen structure."
  else if (is_glycan)
    "Glycan-dependent epitopes vary with cell-line glycosylation patterns. Binding may differ between recombinant antigen (CHO/HEK) and endogenous tumor antigen."
  else
    "Epitope type undetermined."

  cross_reactivity_risk <- if (is_conformational)
    "Moderate — conformational epitopes are often more target-specific but may cross-react with homologous folded domains."
  else if (is_linear)
    "Low-to-moderate — linear epitopes may cross-react with sequence homologs. BLAST the epitope sequence against the human proteome."
  else
    "Unknown — characterize with mutagenesis or peptide scanning."

  adc_relevance <- if (is_conformational)
    "Conformational epitopes can be retained after receptor internalization and lysosomal trafficking. Generally favorable for ADC targeting."
  else if (is_glycan)
    "Glycan epitopes may be less stable under lysosomal acidification. Consider antigen density and glycan heterogeneity between patients."
  else
    "Linear epitopes on surface domains are generally compatible with ADC delivery."

  list(
    epitope_region        = epitope_region,
    epitope_type_detected = if (is_conformational) "Conformational"
                            else if (is_glycan) "Glycan-dependent"
                            else "Linear",
    accessibility_note    = accessibility_note,
    cross_reactivity_risk = cross_reactivity_risk,
    adc_relevance         = adc_relevance,
    conformation_sensitive = conformation_sensitive
  )
}

# ── FcRn binding / antibody half-life prediction ──────────────────────────────

.igg_base_halflife <- c(
  "IgG1"    = 21,
  "IgG2"    = 20,
  "IgG3"    = 7,
  "IgG4"    = 21,
  "IgG1 (YTE)"        = 31.5,
  "IgG1 (LS)"         = 27.3,
  "IgG1 (M428L/N434S)"= 31.5,
  "IgG1 (GASDALIE)"   = 14
)

.fc_mutation_effects <- data.table::data.table(
  Mutation     = c("YTE (M252Y/S254T/T256E)", "LS (M428L/N434S)", "M428L/N434S",
                   "GASDALIE (G236A/S239D/A330L/I332E)", "ALAYT", "N434H",
                   "Xtend (M428L/N434S — same as LS)", "Wild-type (none)"),
  Fold_change  = c(1.5, 1.3, 1.5, 0.67, 1.2, 1.3, 1.5, 1.0),
  Mechanism    = c(
    "Increased FcRn binding at pH 6.0; faster recycling",
    "FcRn affinity enhancement; similar to YTE but different residues",
    "Increased FcRn pH 6.0 affinity; prolonged recycling",
    "Enhanced FcγRIII binding for ADCC; slightly shorter half-life due to faster clearance",
    "Moderate FcRn improvement",
    "Enhanced FcRn binding",
    "Same as M428L/N434S",
    "Baseline IgG1 FcRn binding"
  )
)

predict_halflife <- function(igg_subclass = "IgG1", fc_mutations = character(0),
                              dar = 4, conjugation_site = "Cys-engineered") {
  base_hl <- .igg_base_halflife[igg_subclass]
  if (is.na(base_hl)) base_hl <- 21

  mut_multiplier <- 1.0
  mut_effects <- character(0)
  for (m in fc_mutations) {
    row <- .fc_mutation_effects[grepl(m, Mutation, ignore.case = TRUE)]
    if (nrow(row) > 0) {
      mut_multiplier <- mut_multiplier * row$Fold_change[1]
      mut_effects <- c(mut_effects, paste0(row$Mutation[1], " (×", row$Fold_change[1], ")"))
    }
  }

  dar_penalty <- if (dar <= 2) 0.05
                 else if (dar <= 4) 0.15
                 else if (dar <= 6) 0.25
                 else 0.35

  site_modifier <- switch(conjugation_site,
    "Cys-engineered"       = 0.0,
    "Lys (NHS ester / non-specific)" = 0.05,
    "Fab-Cys"              = 0.10,
    "Fc-glycan"            = 0.05,
    "N-term"               = 0.08,
    0.10
  )

  total_penalty <- dar_penalty + site_modifier
  predicted_hl  <- base_hl * mut_multiplier * (1 - total_penalty)

  tier <- if (predicted_hl >= 25) "Excellent (> 25 days)"
          else if (predicted_hl >= 18) "Good (18–25 days)"
          else if (predicted_hl >= 10) "Moderate (10–18 days)"
          else "Short (< 10 days)"

  list(
    igg_subclass          = igg_subclass,
    fc_mutations          = fc_mutations,
    mut_effects_applied   = mut_effects,
    base_halflife_days    = round(base_hl, 1),
    mutation_multiplier   = round(mut_multiplier, 2),
    dar_penalty_pct       = round(dar_penalty * 100, 0),
    site_penalty_pct      = round(site_modifier * 100, 0),
    predicted_halflife    = round(predicted_hl, 1),
    halflife_range        = c(round(predicted_hl * 0.85, 1), round(predicted_hl * 1.15, 1)),
    tier                  = tier,
    dar                   = dar,
    conjugation_site      = conjugation_site,
    note = paste0(
      "Predicted t½ based on ", igg_subclass, " baseline (", round(base_hl, 0), " d), ",
      if (length(mut_effects) > 0) paste0("Fc mutations: ", paste(mut_effects, collapse = ", "), "; ") else "no Fc mutations; ",
      "DAR ", dar, " payload penalty: -", round(dar_penalty * 100), "%; ",
      "site penalty (", conjugation_site, "): -", round(site_modifier * 100), "%."
    )
  )
}

# ── PK & Bystander effect (simple models) ────────────────────────────────────

simulate_adc_pk <- function(dose_mg_kg = 3, bw_kg = 70, halflife_days,
                             dar = 4, deconjugation_halflife_days = 5, days = 21) {
  t    <- seq(0, days, by = 0.25)
  Vd   <- 3.5  # L/kg approximate distribution volume for IgG
  dose_ug_ml <- (dose_mg_kg * 1e3) / Vd  # rough initial plasma conc (µg/mL)

  ke_adc  <- log(2) / halflife_days
  ke_deconj <- log(2) / deconjugation_halflife_days

  intact_adc <- dose_ug_ml * exp(-ke_adc * t)
  naked_ab   <- dose_ug_ml * (exp(-ke_adc * t) - exp(-ke_deconj * t)) *
                (ke_deconj / (ke_deconj - ke_adc + 1e-9))

  data.frame(
    Day        = t,
    Intact_ADC = pmax(intact_adc, 0),
    Naked_Ab   = pmax(naked_ab, 0)
  )
}

score_bystander_effect <- function(payload_name, payload_table = NULL) {
  if (is.null(payload_table)) payload_table <- get_payload_table()
  row <- payload_table[toupper(Name) == toupper(payload_name)]
  if (nrow(row) == 0) return(list(score = NA, tier = "Unknown", notes = "Payload not found."))

  permeable <- row$Cell_permeable_bystander[1]
  ic50      <- row$IC50_nM[1]

  score <- 0L
  if (isTRUE(permeable)) score <- score + 5L
  if (ic50 < 1)   score <- score + 3L
  else if (ic50 < 10) score <- score + 1L

  tier  <- if (score >= 7) "High" else if (score >= 4) "Moderate" else "Low"
  notes <- if (isTRUE(permeable))
    paste0(payload_name, " is cell-permeable → strong bystander killing in antigen-heterogeneous tumors. ",
           "IC50 = ", ic50, " nM. Bystander benefit is highest when tumor antigen expression is heterogeneous.")
  else
    paste0(payload_name, " is NOT cell-permeable → limited bystander effect. ",
           "Best suited for targets with uniformly high antigen density. Consider switching to a permeable payload (MMAE, DXd) if antigen heterogeneity is a concern.")

  list(score = score, tier = tier, permeable = permeable, ic50_nm = ic50, notes = notes)
}
