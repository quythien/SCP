# Developer note: this file holds the package-level NAMESPACE directives
# (useDynLib + imports) in one obvious place; the compiled cosinor routines
# live in src/cosinor_fast.cpp. The ?SCP help page shows the description
# below; its title is inherited from the DESCRIPTION file.

#' @description
#' SCP (Simulation-based Circadian Power) is a framework for sample-size
#' planning in circadian transcriptomic studies. It calibrates gene-level
#' distributions of amplitude, noise, phase, and time-of-day sampling from a
#' user's pilot dataset to capture transcriptome-wide signal heterogeneity,
#' then estimates FDR-controlled power for single-cohort rhythmic-biomarker
#' detection and for two-group differential analyses of rhythmicity (DR),
#' phase (DP), and mesor (DM). A bootstrap layer propagates pilot-estimation
#' uncertainty into confidence intervals on every power estimate. A curated
#' database of public circadian pilots and a Shiny app are bundled for
#' point-and-click planning.
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
