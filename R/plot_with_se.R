#' =======================================================================
#' Plotting Helper: Standard Error Bars
#' =======================================================================
#'
#' add_se_bars() draws +/-1 SE error bars on a base-R plot. It is called
#' internally by plotSingleCohortPower() and plotDiffPower() when rendering
#' their stratified-power panels.


# =====================================================================
# HELPER: add error bars to a plot
# =====================================================================
#' Add +/-1 SE Error Bars to a Base-R Plot
#'
#' @param x Numeric vector. x-coordinates of the bar centres.
#' @param y Numeric vector. y-coordinates (point estimates).
#' @param se Numeric vector. Standard errors (bar half-width).
#' @param col Colour passed to \code{arrows()}.
#' @param bar_width Numeric. Relative width of the bar caps (default 0.3).
#'
#' @return Invisibly returns \code{NULL}. Called for its side-effect of
#'   drawing arrows on the active plot device.
add_se_bars <- function(x, y, se, col, bar_width = 0.3) {
  valid <- !is.na(y) & !is.na(se) & se > 0
  if (!any(valid)) return(invisible())
  arrows(x[valid], (y - se)[valid], x[valid], (y + se)[valid],
         angle = 90, code = 3, length = bar_width * 0.15,
         col = col, lwd = 1.2)
}

