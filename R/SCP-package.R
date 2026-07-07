#' SCP: Simulation-Based Circadian Power Analysis with K-harmonic LRT
#'
#' @description
#' Package-level imports and the dynamic-library registration for the
#' compiled cosinor routines in \code{src/cosinor_fast.cpp}. Kept in its own
#' file (rather than folded into an arbitrary R/*.R file) so it is obvious
#' where NAMESPACE's \code{useDynLib}/\code{import} directives come from.
#'
#' @keywords internal
#' @importFrom Rcpp sourceCpp
#' @import stats
#' @import graphics
#' @import grDevices
#' @import utils
#' @import ggplot2
#' @useDynLib SCP, .registration = TRUE
"_PACKAGE"

## Non-standard-evaluation column names used inside ggplot2::aes()/base
## plotting calls throughout R/plot_*.R and R/utils.R. These are data.frame
## column names resolved at plot time, not missing objects; R CMD check
## cannot tell the difference statically, hence the whitelist below.
## Deliberately NOT included: dp_power_raw, dp_phase_results,
## density_results, diff_rhythmicity_permutation, diff_rhythmicity_bootstrap
## -- those are genuine bugs (see audit notes), not NSE bindings, and must
## keep surfacing in R CMD check until fixed.
utils::globalVariables(c(
  "B", "B_fac", "DM", "DP", "DR", "N", "N_jit",
  "dataset", "label", "n", "n_A", "n_B",
  "omega_fac", "power", "power_sd", "stratum", "sweep_fac", "y"
))
