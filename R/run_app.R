# run_app.R - ADC Peptide Mapper v1.0

#' Launch the ADC Peptide Mapper Shiny application
#'
#' @description Opens the interactive ADC Peptide Mapper interface. The app
#'   provides a 12-tab workflow for in-silico ADC peptide mapping, MRM
#'   transition list generation, DAR modeling, MRM quality assessment, ADC
#'   design tools, and an AI-powered assistant.
#'
#' @param ... Additional arguments passed to \code{\link[shiny]{runApp}},
#'   such as \code{host}, \code{port}, or \code{launch.browser}.
#'
#' @return Does not return; starts a blocking Shiny server process.
#'
#' @export
#' @examples
#' \dontrun{
#'   ADCPeptideMapper::run_app()
#'
#'   # Run on a specific port
#'   ADCPeptideMapper::run_app(port = 7654)
#' }
run_app <- function(...) {
  app_dir <- system.file("shiny", package = "ADCPeptideMapper")
  if (!nzchar(app_dir)) {
    stop(
      "Shiny app directory not found inside the ADCPeptideMapper package. ",
      "Try reinstalling with: remotes::install_github('nishi76/ADC_Peptide_Mapper_v1.0')"
    )
  }
  shiny::runApp(app_dir, ...)
}
