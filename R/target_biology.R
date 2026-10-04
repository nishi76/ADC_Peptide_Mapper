# target_biology.R  -  ADC Peptide Mapper v1.0
# Internalization scoring and surface accessibility assessment

# -- Internalization database --------------------------------------------------

.internalization_db <- data.table::data.table(
  Target = c(
    "HER2","EGFR","CD22","CD33","CD19","CD79b","Trop2","FRalpha","BCMA","HER3",
    "Nectin-4","c-MET","CD30","CD138","ROR1","ROR2","AXL","MSLN","CEACAM5",
    "FOLR1","LY6E","PTK7","CD70","ENPP3","LRRC15","CD123","FLT3","CD25",
    "EGFR vIII","EphA2","GD2","TAG-72","MUC16","STEAP1","B7-H3","B7-H4",
    "DLL3","SLITRK6","KAAG1","PD-L1"
  ),
  Score = c(
    9, 8, 9, 8, 8, 8, 8, 8, 7, 6,
    8, 8, 8, 6, 7, 7, 8, 5, 6,
    8, 7, 8, 8, 7, 5, 8, 8, 7,
    7, 8, 4, 4, 5, 6, 6, 5,
    7, 7, 6, 4
  ),
  Tier = c(
    "High","High","High","High","High","High","High","High","High","Moderate",
    "High","High","High","Moderate","High","High","High","Moderate","Moderate",
    "High","High","High","High","High","Moderate","High","High","High",
    "High","High","Low","Low","Moderate","Moderate","Moderate","Moderate",
    "High","High","Moderate","Low"
  ),
  Receptor_class = c(
    "RTK","RTK","B-cell co-receptor","Myeloid surface antigen","B-cell antigen",
    "B-cell co-receptor","EpCAM family","Folate receptor","TNFR superfamily","RTK",
    "Nectin family","RTK","TNFR superfamily","Proteoglycan","ROR RTK","ROR RTK",
    "TAM RTK","GPI-anchored","CEA family","GPI folate receptor","LY6 family",
    "Orphan RTK","TNFR superfamily","Ectonucleotidase","Fibronectin domain",
    "IL-3R alpha","RTK (FLT3)","IL-2R alpha","EGFRvIII","Eph RTK",
    "Ganglioside (GPI)","Mucin-associated","Mucin","STEAP metalreductase",
    "Ig superfamily (B7)","Ig superfamily (B7)","Notch ligand","SLIT receptor",
    "Unknown","Checkpoint ligand"
  ),
  Mechanism = c(
    "Receptor-mediated endocytosis (clathrin); ligand-independent constitutive",
    "Ligand-induced clathrin endocytosis; constitutive recycling",
    "Constitutive internalization via clathrin and dynamin",
    "Receptor-mediated endocytosis; lysosomal trafficking",
    "Clathrin-mediated; B-cell receptor complex co-internalization",
    "BCR complex co-internalization",
    "Constitutive clathrin-mediated endocytosis",
    "GPI-anchored; caveolae-dependent internalization",
    "APRIL/BAFF-induced endocytosis; constitutive low-level",
    "EGF-induced; lower rate vs HER2; recycling predominates",
    "PVRL4-induced clathrin-mediated",
    "HGF-induced; rapid lysosomal trafficking",
    "CD30L-induced; rapid internalization to lysosome",
    "Heparan sulfate PG; moderate clathrin internalization",
    "Wnt5a-induced; clathrin-dependent",
    "Wnt5b-induced; similar to ROR1",
    "GAS6/Protein S-induced; rapid lysosomal trafficking",
    "GPI-anchored; limited cytoplasmic domain; slow endocytosis",
    "CEA family; moderate clathrin-mediated",
    "GPI-anchored; caveolae; slow relative to transmembrane receptors",
    "GPI-anchored; moderate",
    "Wnt7b-induced; clathrin-dependent",
    "CD70L-induced; rapid TNF receptor internalization",
    "Constitutive; pyrophosphate substrate-induced",
    "Limited internalization; fibronectin domain stroma marker",
    "Cytokine-induced; rapid clathrin endocytosis",
    "FLT3L-induced rapid endocytosis; AML blast surface",
    "IL-2-induced; moderate rate",
    "Constitutive; truncated EGFR; some endosomal escape",
    "EphrinA1-induced; bidirectional; rapid",
    "GPI-ganglioside; minimal cytoplasmic machinery; low internalization",
    "Shed mucin; low surface retention",
    "Shed mucin; co-internalization limited",
    "Iron reductase; constitutive; moderate",
    "Checkpoint ligand; PDL1 recycling; slow constitutive",
    "Checkpoint ligand; limited cytoplasmic domain recycling",
    "Notch ligand; DLL3-high neuroendocrine; constitutive",
    "LRR; neuronal; constitutive moderate",
    "Unknown; early discovery",
    "Checkpoint; surface shedding reduces effective concentration"
  ),
  Literature_note = c(
    "Clathrin pathway dominant; t1/2 internalization ~30 min (Neve 2004)",
    "Rapid recycling after internalization (Wiley 2003)",
    "t1/2 ~10 min (CD22 ligand-independent) (Shan 2000)",
    "Rapid lysosomal trafficking; gemtuzumab basis (Sutherland 2006)",
    "BCR-dependent co-internalization (Press 1989)",
    "Polyclonal CD79b ADC validated in NHL (Dornan 2009)",
    "Trop2 constitutive internalization supports sacituzumab (Cubas 2019)",
    "GPI-caveolae pathway; slower than TM receptors (Elnakat 2004)",
    "BCMA ADC data from belantamab mafodotin (Tai 2014)",
    "HER3 recycling limits payload delivery (Schoeberl 2009)",
    "Nectin-4 ADC basis for enfortumab vedotin (Challita-Eid 2016)",
    "Rapid MET endocytosis to lysosome (Birchmeier 2003)",
    "One of the fastest-internalizing TNF receptors (Falini 1995)",
    "Heparan sulfate retention reduces internalization speed (Bhatt 2016)",
    "ROR1 ADC studies (Zhang 2019)",
    "ROR2 internalization demonstrated (Daneshmanesh 2013)",
    "AXL rapid internalization post-GAS6 (Sasaki 2006)",
    "GPI anchor limits cytoplasmic machinery (Hassan 2016)",
    "CEA internalization rate lower than RTKs (Beauchemin 1999)",
    "GPI-anchored; folate-induced internalization (Parker 2005)",
    "GPI; LY6E overexpression in triple-negative BC (Ly 2020)",
    "PTK7 ADC preclinical internalization (Damelin 2017)",
    "Rapid TNF receptor internalization (Marschner 2020)",
    "ENPP3 ADC studies (Challita-Eid 2016)",
    "LRRC15 stroma target; limited internalization in fibroblasts (Purcell 2018)",
    "CD123 AML internalization (Testa 2014)",
    "FLT3 ADC data (Kubasch 2021)",
    "IL-2R internalization dynamics (Smith 1990)",
    "EGFRvIII constitutive; limited recycling (Sok 2006)",
    "EphA2 ADC rapid trafficking (Thundimadathil 2012)",
    "Ganglioside GD2; low TM internalization (Cheung 2012)",
    "TAG-72 mucin shedding limits ADC payload delivery (Colcher 1987)",
    "MUC16 shedding; co-internalization limited (Scholler 1999)",
    "STEAP1 iron reductase; ADC in prostate (Zhao 2021)",
    "B7-H3 constitutive cycling (Picarda 2016)",
    "B7-H4 limited cytoplasmic signaling domain (Podojil 2020)",
    "DLL3 neuroendocrine; rovalpituzumab basis (Saunders 2015)",
    "SLITRK6 bladder cancer; ADC target (Matsui 2015)",
    "KAAG1 early-stage research",
    "PD-L1 recycling; shedding; limited lysosomal trafficking (Burr 2017)"
  )
)

# -- Internalization scoring ---------------------------------------------------

#' Score target receptor internalization for ADC payload delivery
#'
#' @description Looks up an ADC target by name in a curated internalization
#'   database (40 validated targets) and returns a numerical score (1-10),
#'   tier classification, receptor class, mechanism, and key literature
#'   references. Partial name matching is used when exact lookup fails.
#'
#' @param target_name character(1). Target antigen name (e.g. \code{"HER2"},
#'   \code{"CD22"}, \code{"Trop2"}). Case-insensitive.
#'
#' @return Named list with: \code{found} (logical), \code{target}
#'   (character), \code{score} (integer, 1-10), \code{tier} ("High" /
#'   "Moderate" / "Low" / "Unknown"), \code{receptor_class}, \code{mechanism},
#'   \code{literature}, \code{source} ("curated_db" or "not_found").
#'
#' @export
#' @examples
#' score_internalization("HER2")
#' score_internalization("Trop2")
#' score_internalization("GD2")    # low score
score_internalization <- function(target_name) {
  target_clean <- toupper(trimws(target_name))
  db <- .internalization_db
  db[, Target_upper := toupper(Target)]

  match_row <- db[Target_upper == target_clean]

  if (nrow(match_row) == 0) {
    partial <- db[grepl(target_clean, Target_upper, fixed = TRUE)]
    if (nrow(partial) > 0) match_row <- partial[1]
  }

  if (nrow(match_row) > 0) {
    list(
      found          = TRUE,
      target         = match_row$Target[1],
      score          = match_row$Score[1],
      tier           = match_row$Tier[1],
      receptor_class = match_row$Receptor_class[1],
      mechanism      = match_row$Mechanism[1],
      literature     = match_row$Literature_note[1],
      source         = "curated_db"
    )
  } else {
    list(
      found          = FALSE,
      target         = target_name,
      score          = NA_integer_,
      tier           = "Unknown",
      receptor_class = "Unknown",
      mechanism      = paste0("Target '", target_name, "' not found in curated database. ",
                              "Please consult internalization assay literature for this target."),
      literature     = "No curated reference available.",
      source         = "not_found"
    )
  }
}

# -- Surface accessibility assessment -----------------------------------------

.surface_db <- data.table::data.table(
  Target = c(
    "HER2","EGFR","CD22","CD33","Trop2","BCMA","Nectin-4","HER3","c-MET",
    "MSLN","FRalpha","B7-H3","DLL3","CEACAM5","MUC16","STEAP1","CD19","CD79b",
    "CD30","CD138","ROR1","AXL","FOLR1","PTK7","CD123","FLT3","PD-L1"
  ),
  Extracellular_domains = c(
    4L, 4L, 2L, 2L, 1L, 1L, 4L, 4L, 4L,
    1L, 1L, 3L, 3L, 7L, 1L, 2L, 1L, 2L,
    1L, 1L, 3L, 3L, 1L, 3L, 2L, 5L, 1L
  ),
  TM_count = c(
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
    0L, 0L, 1L, 1L, 1L, 1L, 6L, 1L, 1L,
    1L, 1L, 1L, 1L, 0L, 1L, 1L, 1L, 1L
  ),
  Surface_expression_level = c(
    "High","High","High","Moderate","High","Moderate","High","Moderate","Moderate",
    "High","Moderate","High","Moderate","Moderate","High","Moderate","High","Moderate",
    "Moderate","High","Moderate","Moderate","Moderate","Moderate","Low-Moderate","Moderate","Moderate"
  ),
  Accessibility_score = c(
    92L, 90L, 85L, 78L, 88L, 72L, 86L, 80L, 82L,
    85L, 70L, 84L, 75L, 68L, 60L, 58L, 87L, 76L,
    80L, 74L, 78L, 80L, 65L, 77L, 70L, 78L, 76L
  ),
  Accessibility_tier = c(
    "High","High","High","High","High","Moderate","High","High","High",
    "High","Moderate","High","High","Moderate","Moderate","Moderate","High","High",
    "High","High","High","High","Moderate","High","Moderate","High","High"
  ),
  Shed_or_secreted = c(
    FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE,
    FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE,
    FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE
  ),
  Notes = c(
    "Large 4-domain ECD; highly accessible at tumor surface; HER2+ amplified cells show very high density",
    "Large ECD; accessible; EGFR amplification increases surface density",
    "Accessible B-cell surface receptor; high on malignant B-cells",
    "Myeloid surface antigen; moderate density on AML blasts",
    "Single ECD domain; strong surface expression on epithelial tumors",
    "Shed by ADAM10/17; surface levels variable; consider serum BCMA as biomarker",
    "Nectin-4 overexpressed in bladder/breast; accessible TM protein",
    "HER3 expression lower than HER2; accessible but expression heterogeneous",
    "High MET expression in NSCLC/GBM; large accessible ECD",
    "GPI-anchored; surface exposed; shed form in plasma can act as decoy",
    "GPI-anchored folate receptor; exposed apical surface; shed form exists",
    "Broad tumor expression; 3 Ig-like domains accessible",
    "Notch ligand; high in SCLC neuroendocrine; accessible ECD",
    "CEA family 7-domain structure; partially shed; lower effective surface density",
    "Mucin; shed MUC16 (CA-125) in serum; reduces effective surface binding",
    "6-TM iron reductase; less accessible than single-TM receptors",
    "Type I TM; single ECD; high on B-ALL/CLL",
    "Ig-like ECD; moderate; part of BCR complex",
    "TNFR superfamily; single ECD; expressed on Hodgkin lymphoma",
    "Heparan sulfate PG; single ECD but glycosaminoglycan chains project outward",
    "CRD + kringle domains; accessible; expressed in hematologic malignancies",
    "TAM RTK 3-domain ECD; accessible; widely expressed in solid tumors",
    "GPI-anchored; accessible on apical surface; shed form (CA-125 not to confuse)",
    "3 CRD domains; accessible on colorectal/breast/cervical cancer",
    "IL-3R alpha; expressed on AML/BPDCN blasts; single ECD",
    "Type III TM; 5-domain ECD; good accessibility on AML",
    "Single Ig-like domain; surface shed in TME reduces effective binding"
  )
)

#' Assess target antigen surface accessibility
#'
#' @description Looks up a target antigen in a curated surface accessibility
#'   database (27 validated targets) and returns an accessibility score (0-100),
#'   tier, number of extracellular domains, transmembrane count, surface
#'   expression level, and shedding/secretion status. Partial name matching
#'   is applied when exact lookup fails.
#'
#' @param target_name character(1). Target antigen name. Case-insensitive.
#'
#' @return Named list with: \code{found} (logical), \code{target},
#'   \code{accessibility_score} (integer, 0-100),
#'   \code{accessibility_tier} ("High" / "Moderate" / "Unknown"),
#'   \code{ec_domains} (integer), \code{tm_count} (integer),
#'   \code{expression_level} (character), \code{shed_secreted} (logical),
#'   \code{notes} (character).
#'
#' @export
#' @examples
#' \dontrun{
#' assess_surface_accessibility("HER2")
#' assess_surface_accessibility("MUC16")   # shed antigen, lower score
#' }
assess_surface_accessibility <- function(target_name) {
  target_clean <- toupper(trimws(target_name))
  db <- .surface_db
  db[, Target_upper := toupper(Target)]

  match_row <- db[Target_upper == target_clean]
  if (nrow(match_row) == 0) {
    partial <- db[grepl(target_clean, Target_upper, fixed = TRUE)]
    if (nrow(partial) > 0) match_row <- partial[1]
  }

  if (nrow(match_row) > 0) {
    list(
      found               = TRUE,
      target              = match_row$Target[1],
      accessibility_score = match_row$Accessibility_score[1],
      accessibility_tier  = match_row$Accessibility_tier[1],
      ec_domains          = match_row$Extracellular_domains[1],
      tm_count            = match_row$TM_count[1],
      expression_level    = match_row$Surface_expression_level[1],
      shed_secreted       = match_row$Shed_or_secreted[1],
      notes               = match_row$Notes[1]
    )
  } else {
    list(
      found               = FALSE,
      target              = target_name,
      accessibility_score = NA_integer_,
      accessibility_tier  = "Unknown",
      ec_domains          = NA_integer_,
      tm_count            = NA_integer_,
      expression_level    = "Unknown",
      shed_secreted       = NA,
      notes               = paste0("Target '", target_name, "' not found in curated surface database.")
    )
  }
}
