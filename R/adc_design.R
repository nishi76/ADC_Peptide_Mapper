# adc_design.R  -  ADC Peptide Mapper v1.0
# Linker chemistry table, payload database, deconjugation prediction

# -- Linker database -----------------------------------------------------------

#' Retrieve the curated ADC linker chemistry reference table
#'
#' @description Returns a data.table containing 14 clinically relevant or
#'   investigational ADC linker chemistries with stability rankings, cleavage
#'   triggers, plasma half-lives, DAR ranges, and approved ADC examples.
#'
#' @return data.table with columns: \code{Name}, \code{Type} ("Cleavable" or
#'   "Non-cleavable"), \code{Chemistry}, \code{Cleavage_trigger},
#'   \code{t_half_plasma_days}, \code{DAR_range}, \code{Stability_rank}
#'   (A+/A/B/C/D), \code{Approved_examples}, \code{Notes}.
#'
#' @export
#' @examples
#' lt <- get_linker_table()
#' lt[Type == "Cleavable", .(Name, Stability_rank, t_half_plasma_days)]
get_linker_table <- function() {
  data.table::data.table(
    Name = c(
      "MC-VC-PABC", "MC-VC", "CL2A", "Sulfo-SPDB", "SMCC",
      "SPDB", "CX1-1", "PEG4-VC-PABC", "Disulfide (SPDP)", "Hydrazone (BMPH)",
      "Oxime (aminooxy)", "NHS-PEG4", "Sortase-tag (LPXTG)", "Click chemistry (DBCO-Az)"
    ),
    Type = c(
      "Cleavable", "Cleavable", "Cleavable", "Cleavable", "Non-cleavable",
      "Cleavable", "Cleavable", "Cleavable", "Cleavable", "Cleavable",
      "Cleavable", "Non-cleavable", "Non-cleavable", "Non-cleavable"
    ),
    Chemistry = c(
      "Maleimide-thiol / cathepsin B dipeptide", "Maleimide-thiol / cathepsin B",
      "Maleimide-thiol / acid-labile hydrazone", "Disulfide / cathepsin B",
      "Maleimide-thiol (thioether)", "Disulfide / cathepsin B",
      "Maleimide-thiol / oxime", "PEGylated maleimide / cathepsin B dipeptide",
      "Disulfide", "Hydrazone",
      "Oxime", "NHS ester / lysine amide",
      "Enzymatic transpeptidation", "Bioorthogonal cycloaddition"
    ),
    Cleavage_trigger = c(
      "Lysosomal cathepsin B (pH 5)", "Lysosomal cathepsin B",
      "Acidic pH (< 5.5)", "Reducing env + cathepsin B",
      "None (hydrolysis-resistant)", "Reducing env + cathepsin B",
      "pH + cathepsin", "Lysosomal cathepsin B",
      "Reducing environment (GSH)", "Acidic pH (< 5.5)",
      "Stable at neutral pH; slow hydrolysis", "None",
      "None (site-specific)", "None (bioorthogonal)"
    ),
    t_half_plasma_days = c(7, 6, 0.5, 1.5, ">30", 2, 5, 9, 1.5, 0.8, 21, ">30", ">30", ">30"),
    DAR_range = c("2-8", "2-8", "2-4", "2-4", "3-4", "2-4", "2-8", "4-8", "2-4", "2-4", "2-4", "4-8", "1-2 (site-specific)", "2 (site-specific)"),
    Stability_rank = c("A", "A", "C", "B", "A+", "B", "A", "A+", "C", "D", "A", "A+", "A+", "A+"),
    Approved_examples = c(
      "Adcetris (brentuximab vedotin), Polivy, Padcev", "Trodelvy (SN-38 variant)",
      "Besylomab (preclinical)", "Kadcyla (ado-trastuzumab emtansine DM4 variant)",
      "Kadcyla (T-DM1)", "Preclinical",
      "Experimental", "Preclinical PEG series",
      "Classic ADC research linker", "Early-generation Mylotarg",
      "Experimental", "Research conjugates",
      "Site-specific research (GlaxoSmithKline)", "Research stage"
    ),
    Notes = c(
      "Gold standard for protease-cleavable ADCs; PABC self-immolative spacer releases free payload",
      "Shorter version without PABC; slower release",
      "Rapid release at endosomal pH; high bystander potential; DAR heterogeneity",
      "Disulfide bond reduced in tumor; more selective than SPDB alone",
      "Stable thioether; payload released only after antibody catabolism in lysosome",
      "Slower than Sulfo-SPDB; moderate stability",
      "Hybrid cleavage mechanism; experimental",
      "Reduced aggregation vs. MC-VC; preferred for hydrophobic payloads",
      "Cleaved by intracellular GSH; loses selectivity in high-GSH plasma",
      "First ADC linker generation; rapid hydrolysis limits half-life",
      "Highly stable at pH 7.4; useful for oxime-based site-specific conjugation",
      "Most stable; DAR heterogeneity on lysines; used for less-toxic payloads",
      "Exact site-specific DAR; requires engineered LPXTG tag on antibody",
      "Emerging bioorthogonal platform; precise DAR control"
    )
  )
}

# -- Payload database ---------------------------------------------------------

#' Retrieve the curated ADC cytotoxic payload reference table
#'
#' @description Returns a data.table with 15 clinically relevant or
#'   investigational ADC payloads covering auristatins, maytansinoids,
#'   camptothecin derivatives, PBD dimers, enediynes, anthracyclines, and
#'   novel mechanisms, with IC50, bystander potential, resistance gene, and
#'   approved ADC information.
#'
#' @return data.table with columns: \code{Name}, \code{Class}, \code{MOA},
#'   \code{Target_molecule}, \code{IC50_nM}, \code{DAR_optimal},
#'   \code{Cell_permeable_bystander} (logical), \code{Resistance_genes},
#'   \code{Approved_ADCs}, \code{Potency_rank} (1-5), \code{Notes}.
#'
#' @export
#' @examples
#' pt <- get_payload_table()
#' pt[Cell_permeable_bystander == TRUE, .(Name, IC50_nM)]
get_payload_table <- function() {
  data.table::data.table(
    Name = c(
      "MMAE", "MMAF", "DM1", "DM4", "DXd (Exatecan derivative)",
      "SN-38", "Calicheamicin", "PBD dimer (SG3249)", "Duocarmycin SA",
      "Pyrrolobenzodiazepine (PBD monomer)", "Doxorubicin", "Topotecan",
      "Camptothecin", "Amatoxin (alpha-amanitin)", "Pseudomonas exotoxin (PE38)"
    ),
    Class = c(
      "Auristatin", "Auristatin", "Maytansinoid", "Maytansinoid", "Camptothecin derivative",
      "Camptothecin", "Enediyne", "PBD dimer", "CBI-seco",
      "PBD monomer", "Anthracycline", "Camptothecin", "Camptothecin",
      "Amatoxin", "Bacterial toxin"
    ),
    MOA = c(
      "Tubulin polymerization inhibitor", "Tubulin polymerization inhibitor",
      "Tubulin polymerization inhibitor", "Tubulin polymerization inhibitor",
      "Topoisomerase I inhibitor", "Topoisomerase I inhibitor",
      "DNA double-strand break (DSB)", "DNA interstrand crosslink",
      "DNA alkylation (minor groove)", "DNA alkylation",
      "Topoisomerase II / intercalation", "Topoisomerase I inhibitor",
      "Topoisomerase I inhibitor", "RNA polymerase II inhibitor",
      "ADP-ribosylation of EF-2 (protein synthesis)"
    ),
    Target_molecule = c(
      "beta-tubulin", "beta-tubulin", "beta-tubulin", "beta-tubulin", "Topoisomerase I",
      "Topoisomerase I", "DNA", "DNA", "DNA", "DNA",
      "Topoisomerase II / DNA", "Topoisomerase I", "Topoisomerase I",
      "RNA Pol II", "Elongation factor 2"
    ),
    IC50_nM = c(0.1, 0.1, 0.01, 0.01, 0.7, 3.0, 0.001, 0.001, 0.01, 0.1,
                10.0, 2.0, 2.0, 0.002, 0.001),
    DAR_optimal = c("4", "4", "3-4", "3-4", "4-8", "4-8",
                    "2-3", "2", "4", "2-4",
                    "4-6", "4-8", "4-8", "2-3", "1 (protein)"),
    Cell_permeable_bystander = c(
      TRUE, FALSE, FALSE, FALSE, TRUE, TRUE,
      FALSE, FALSE, TRUE, FALSE,
      TRUE, TRUE, TRUE, FALSE, FALSE
    ),
    Resistance_genes = c(
      "MDR1 (ABCB1), beta-III tubulin", "MDR1 (lower vs MMAE)", "MDR1, ABCC1",
      "MDR1, ABCC1", "ABCG2, ABCB1", "UGT1A1, ABCG2",
      "MDR1 rare", "Low (DNA crosslink)", "Low", "Low",
      "MDR1, MRP1", "ABCG2", "ABCG2", "Rare (novel MOA)", "Rare"
    ),
    Approved_ADCs = c(
      "Adcetris, Polivy, Padcev, Tivdak, Zynlonta", "None approved",
      "Kadcyla (T-DM1)", "Mirvetuximab soravtansine (ELAHERE)",
      "Enhertu (T-DXd), Dato-DXd, Patritumab deruxtecan",
      "Trodelvy (sacituzumab govitecan)", "Mylotarg, Besylomab",
      "Zynlonta (loncastuximab tesirine, PBD dimer)", "Experimental",
      "Experimental", "None (standalone chemo)", "None", "None",
      "Experimental (HDP-101 BCMA)", "Experimental (RG7787)"
    ),
    Potency_rank = c(3L, 3L, 4L, 4L, 2L, 1L, 5L, 5L, 4L, 4L, 1L, 2L, 2L, 5L, 5L),
    Notes = c(
      "Most-used ADC payload. Cell-permeable -> strong bystander killing. High MDR1 susceptibility.",
      "Non-permeable auristatin. Lower bystander; useful where off-target bystander is risky.",
      "Maytansinoid used in Kadcyla. Non-permeable. Requires high antigen density.",
      "Used in ELAHERE. Disulfide linker preferred. Slightly more lipophilic than DM1.",
      "Potent Topo-I inhibitor. Highly cell-permeable. Wide therapeutic window at DAR 8.",
      "Active metabolite of irinotecan. Used in Trodelvy. Strong bystander effect.",
      "Extreme potency; narrow therapeutic window. Best for high-antigen targets.",
      "DNA interstrand crosslinks. Very potent. Low bystander; used in NHL/AML.",
      "CBI alkylator. Strong bystander potential. Experimental stage.",
      "Single-strand DNA alkylation. Moderate potency. Less studied than PBD dimer.",
      "Classical chemo. Lower potency as ADC payload; used for large TAb conjugates.",
      "Moderate Topo-I inhibitor. Used in preclinical ADCs.",
      "Parent of SN-38 / DXd class. Less stable than derivatives.",
      "Novel MOA (RNA Pol II). Not a substrate for MDR1. Experimental.",
      "Protein synthesis inhibitor. Extremely potent but immunogenic."
    )
  )
}

# -- Deconjugation prediction -------------------------------------------------

#' Predict ADC deconjugation kinetics in plasma
#'
#' @description Models DAR decay over time using a first-order rate equation
#'   parameterised by linker-chemistry-specific plasma half-lives. Returns a
#'   kinetics curve, stability classification, days above therapeutic
#'   threshold, and an interpretive warning message.
#'
#' @param chemistry character(1). Linker chemistry string (must match one of
#'   the \code{switch} cases in the function body; see the source for the full
#'   list, e.g. \code{"Maleimide-thiol (MC-VC-PABC / SMCC)"}).
#' @param site character(1). Conjugation site description (used for
#'   site-specific notes, e.g. \code{"Lys (NHS ester / non-specific)"}).
#' @param initial_dar numeric(1). Starting DAR value (default 4).
#' @param days numeric(1). Simulation time in days (default 21).
#'
#' @return Named list with elements: \code{curve} (data.frame with
#'   \code{Day} and \code{DAR}), \code{stability_class} (character, A+/A/B/C/D),
#'   \code{stability_label} (character), \code{t_half} (numeric, days),
#'   \code{therapeutic_threshold} (numeric), \code{days_above_threshold}
#'   (numeric), \code{warning_msg} (character), \code{chemistry},
#'   \code{site}, \code{initial_dar}.
#'
#' @export
#' @examples
#' \dontrun{
#'   res <- predict_deconjugation("Maleimide-thiol (MC-VC-PABC / SMCC)",
#'                                site = "Cys-engineered")
#'   res$stability_class
#'   head(res$curve)
#' }
predict_deconjugation <- function(chemistry, site, initial_dar = 4, days = 21) {
  t <- seq(0, days, by = 0.5)

  params <- switch(chemistry,
    "Maleimide-thiol (MC-VC-PABC / SMCC)" = list(
      t_half      = if (site == "Lys (NHS ester / non-specific)") 30 else 5.5,
      stability   = "A",
      warning_msg = if (site == "Lys (NHS ester / non-specific)")
        "NHS-ester on lysines is highly stable but produces a heterogeneous DAR distribution."
      else
        "Succinimide ring hydrolysis (retro-Michael) can occur in plasma at physiological pH. Engineered Cys or Fab-Cys sites are more stable than non-specific Cys."
    ),
    "Disulfide (SPDB / Sulfo-SPDB)" = list(
      t_half      = 1.5,
      stability   = "C",
      warning_msg = "Plasma albumin and glutathione can reduce disulfide bonds extracellularly. Short t1/2 in reducing tumor microenvironment can be advantageous for payload release, but reduces systemic stability."
    ),
    "Acid-labile hydrazone (CL2A)" = list(
      t_half      = 0.6,
      stability   = "D",
      warning_msg = "Hydrazone linkers hydrolyze rapidly at acidic pH. Significant payload loss expected within 12-24 h in plasma. Only suitable for highly cytotoxic payloads with minimal off-target effects."
    ),
    "Oxime (aminooxy)" = list(
      t_half      = 21,
      stability   = "A",
      warning_msg = "Oxime linkers are among the most stable. Minimal deconjugation at physiological pH. Cleavage occurs only in acidic lysosomes."
    ),
    "NHS ester / Lysine (SMCC)" = list(
      t_half      = 30,
      stability   = "A+",
      warning_msg = "Thioether bond is extremely stable. DAR heterogeneity from lysine conjugation may affect PK and efficacy. Payload released only after full antibody catabolism in lysosome."
    ),
    "Enzymatic / site-specific (sortase, click)" = list(
      t_half      = 30,
      stability   = "A+",
      warning_msg = "Site-specific conjugation produces homogeneous DAR and superior plasma stability. Requires engineered antibody construct."
    ),
    # default
    list(t_half = 7, stability = "B", warning_msg = "")
  )

  ke <- log(2) / params$t_half
  dar_remaining <- initial_dar * exp(-ke * t)
  therapeutic_threshold <- initial_dar * 0.35

  stability_labels <- c(
    "A+" = "Exceptional (> 30 days)",
    "A"  = "Excellent (7-14 days)",
    "B"  = "Good (3-7 days)",
    "C"  = "Moderate (1-3 days)",
    "D"  = "Poor (< 1 day)"
  )

  list(
    curve        = data.frame(Day = t, DAR = pmax(dar_remaining, 0)),
    stability_class  = params$stability,
    stability_label  = stability_labels[params$stability],
    t_half           = params$t_half,
    therapeutic_threshold = therapeutic_threshold,
    days_above_threshold  = max(t[dar_remaining >= therapeutic_threshold], 0),
    warning_msg      = params$warning_msg,
    chemistry        = chemistry,
    site             = site,
    initial_dar      = initial_dar
  )
}

# -- Deconjugation plot --------------------------------------------------------

#' Plot ADC deconjugation kinetics curve
#'
#' @description Creates a ggplot2 line plot of DAR over time from the output
#'   of \code{\link{predict_deconjugation}}, with a shaded sub-therapeutic
#'   zone and dashed threshold line.
#'
#' @param result list. Output of \code{\link{predict_deconjugation}}.
#'
#' @return A \code{ggplot2} object. Print or save with \code{ggsave()}.
#'
#' @export
#' @examples
#' \dontrun{
#'   res <- predict_deconjugation("Maleimide-thiol (MC-VC-PABC / SMCC)",
#'                                site = "Cys-engineered")
#'   plot_deconjugation(res)
#' }
plot_deconjugation <- function(result) {
  df <- result$curve
  thresh <- result$therapeutic_threshold

  p <- ggplot2::ggplot(df, ggplot2::aes(x = Day, y = DAR)) +
    ggplot2::annotate("rect",
      xmin = 0, xmax = max(df$Day),
      ymin = 0, ymax = thresh,
      fill = "#fef3c7", alpha = 0.6) +
    ggplot2::annotate("text",
      x = max(df$Day) * 0.92, y = thresh * 0.5,
      label = "Sub-therapeutic\nzone", size = 3,
      color = "#92400e", hjust = 1) +
    ggplot2::geom_hline(yintercept = thresh,
      linetype = "dashed", color = "#f59e0b", linewidth = 0.7) +
    ggplot2::geom_line(color = "#0d9488", linewidth = 1.8) +
    ggplot2::geom_point(data = df[df$Day %% 3 == 0, ],
      color = "#0d9488", size = 2.5) +
    ggplot2::labs(
      title = paste0("DAR Decay - ", result$chemistry),
      subtitle = paste0("Site: ", result$site,
        "  |  t1/2 = ", result$t_half, " d",
        "  |  Stability class: ", result$stability_class,
        " (", result$stability_label, ")"),
      x = "Days post-administration",
      y = "Mean DAR (drug-to-antibody ratio)"
    ) +
    ggplot2::scale_y_continuous(limits = c(0, NA), expand = ggplot2::expansion(mult = c(0, 0.08))) +
    ggplot2::scale_x_continuous(breaks = seq(0, max(df$Day), by = 3)) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", color = "#1e293b"),
      plot.subtitle = ggplot2::element_text(color = "#64748b", size = 11),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "#e2e8f0"),
      axis.line = ggplot2::element_line(color = "#cbd5e1")
    )
  p
}
