#' Launch the SCP Shiny app
#'
#' Opens the bundled simulation-based power-analysis app in the user's
#' default browser. The app lets the user pick a bundled pilot, set
#' target power and FDR, and read off the recommended sample size for
#' single-cohort rhythmicity detection.
#'
#' @param ... Additional arguments passed to \code{shiny::runApp()}.
#'
#' @return Invisibly returns \code{NULL}. Called for its side effect of
#'   launching the Shiny app.
#'
#' @examples
#' \dontrun{
#'   launchShiny()
#' }
#' @export
launchShiny <- function(...) {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Install the `shiny` package first: install.packages(\"shiny\").")
  }
  appdir <- system.file("shiny", package = "SCP")
  if (!nzchar(appdir)) {
    stop("Shiny app directory not found. Reinstall SCP from a source that ",
         "includes `inst/shiny/app.R`.")
  }
  shiny::runApp(appdir, ...)
}
