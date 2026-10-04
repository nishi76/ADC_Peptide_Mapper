#' @import data.table
#' @importFrom stats dpois setNames
#' @importFrom utils write.csv
#' @importFrom shiny runApp
#' @importFrom bs4Dash bs4DashPage
#' @importFrom DT datatable
#' @importFrom dplyr filter
#' @importFrom htmltools tags
#' @importFrom httr2 request
#' @importFrom shinyjs useShinyjs
#' @importFrom stringr str_detect
#' @importFrom ggplot2 ggplot
"_PACKAGE"

utils::globalVariables(c(
  ".", ".I", ".N", ":=",
  "ADCName", "Adduct", "Chain", "CollisionEnergy", "ComponentName",
  "CompoundName", "DAR", "Day", "End", "Enzyme", "Evalue", "FDR_1pct",
  "FragmentIon", "IsDecoy", "Length", "MC", "Mass", "Modifications",
  "ModifiedSequence", "Mutation", "Name", "PeptideLength", "PeptideSequence",
  "PrecursorCharge", "PrecursorMz", "ProductCharge", "ProductMz",
  "ProteinName", "ReactionMonitor", "RecommendedIsotope", "Score",
  "Sequence", "Start", "Target", "Target_upper", "TransitionName",
  "UniqueToADC", "dpois", "q_value", "rbindlist", "setcolorder", "setorder"
))
