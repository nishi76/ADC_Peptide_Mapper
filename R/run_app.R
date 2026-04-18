#' Launch the ADC Peptide Mapper Shiny App
#' @export
run_app <- function(...) {
  app_dir <- system.file(package = "ADCPeptideMapper")
  if (app_dir == "") stop("ADCPeptideMapper package not found.")
  shiny::runApp(app_dir, ...)
}
