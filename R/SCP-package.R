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
#' @section Main functions:
#' \strong{Pilot database:} \code{\link{scp_pilots}}, \code{\link{scp_load_pilot}},
#'   \code{\link{scp_pilot_search}}.
#'
#' \strong{Calibrate your own pilot:} \code{\link{estCircadianParam}},
#'   \code{\link{estCircadianParam2H}}, \code{\link{estCircadianParamTwoGroup}},
#'   \code{\link{prepCircadianData}}.
#'
#' \strong{Design and options:} \code{\link{CircadianBioOptions}},
#'   \code{\link{CircadianDesignOptions}}, \code{\link{CircadianAnalysisOptions}},
#'   \code{\link{CircadianBootstrapOptions}}.
#'
#' \strong{Run power analysis:} \code{\link{runSimsSingleCohort}},
#'   \code{\link{runDifferentialPower}}, \code{\link{runBootstrapDesignGrid}},
#'   \code{\link{npower}}, \code{\link{circaPowerApproxN80}},
#'   \code{\link{recommendDesign}}.
#'
#' \strong{Plot results:} \code{\link{plotSingleCohortPower}},
#'   \code{\link{plotDiffPower}}, \code{\link{plotBootstrapDesignGrid}}.
#'
#' \strong{Detect rhythmicity:} \code{\link{detect_cosinor}}.
#'
#' \strong{Simulate and helpers:} \code{\link{simCircadianSingleCohort2H}},
#'   \code{\link{makeAdaptiveRStrata}}.
#'
#' \strong{App:} \code{\link{launchShiny}}.
#' @importFrom Rcpp sourceCpp
#' @import stats
#' @import graphics
#' @import grDevices
#' @import utils
#' @import ggplot2
#' @useDynLib SCP, .registration = TRUE
"_PACKAGE"

## These names are data-frame columns referenced inside ggplot2 aes() and base
## plotting calls, resolved when a plot is drawn. R CMD check reads them as
## undefined variables, so listing them here silences that false alarm.
utils::globalVariables(c(
  "B", "B_fac", "DM", "DP", "DR", "N", "N_jit",
  "dataset", "label", "n", "n_A", "n_B",
  "omega_fac", "power", "power_sd", "stratum", "sweep_fac", "y"
))
