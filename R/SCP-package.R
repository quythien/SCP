# Developer note: this file holds the package-level NAMESPACE directives
# (useDynLib + imports) in one obvious place; the compiled cosinor routines
# live in src/cosinor_fast.cpp. The ?SCP help page shows the description
# below; its title is inherited from the DESCRIPTION file.

#' @description
#' SCP (Simulation-based Circadian Power) is a framework for sample-size
#' planning in circadian omics studies, where limited sample availability,
#' uncontrolled collection times, and the multiple-testing burden of
#' genome-wide analysis routinely constrain the design. It calibrates the
#' analysis to a user's pilot dataset, drawing on the pilot's gene-level
#' distributions of amplitude, noise, and phase together with its
#' time-of-day sampling to capture transcriptome-wide signal heterogeneity.
#' From this calibration it estimates FDR-controlled power for two goals:
#' single-cohort detection of rhythmic biomarkers, and two-group
#' differential analysis of rhythmicity, phase, or mesor. A bootstrap layer
#' quantifies the uncertainty that the finite pilot contributes to each
#' power estimate. Detection supports single-harmonic and two-harmonic
#' cosinor models. It bundles curated
#' circadian transcriptomic pilots spanning many tissues and several
#' species, plus a Shiny application for point-and-click power calculation
#' and biomarker exploration.
#'
#' @section GTEx input time:
#' Active GTEx pilots use release v10. The default input time is donor-internal
#' phase (DIP), estimated offline with CHIRAL from clock-gene expression across
#' tissues. Season-adjusted time of death (TOD) is available in the package by
#' setting \code{time_basis = "TOD"} in \code{\link{scp_load_pilot}}. The app
#' uses DIP for bundled GTEx pilots. CHIRAL is not rerun during power simulation.
#' Gene acrophase is estimated relative to the selected input time and is a
#' separate quantity. Use \code{citation("SCP")} for software and source citations.
#'
#' @references
#' Talamanca L, Gobet C, Naef F (2023). Sex-dimorphic and age-dependent
#' organization of 24-hour gene expression rhythms in humans.
#' Science 379:478-483. \doi{10.1126/science.add0846}.
#'
#' GTEx Consortium (2020). The GTEx Consortium atlas of genetic regulatory
#' effects across human tissues. Science 369:1318-1330.
#' \doi{10.1126/science.aaz1776}.
#'
#' @section Main functions:
#' \describe{
#'   \item{Pilot database}{
#'     \itemize{
#'       \item \code{\link{scp_pilots}}
#'       \item \code{\link{scp_load_pilot}}
#'       \item \code{\link{scp_pilot_search}}
#'     }
#'   }
#'   \item{Calibrate your own pilot}{
#'     \itemize{
#'       \item \code{\link{estCircadianParam}}
#'       \item \code{\link{estCircadianParam2H}}
#'       \item \code{\link{estCircadianParamTwoGroup}}
#'       \item \code{\link{prepCircadianData}}
#'     }
#'   }
#'   \item{Design and options}{
#'     \itemize{
#'       \item \code{\link{CircadianBioOptions}}
#'       \item \code{\link{CircadianDesignOptions}}
#'       \item \code{\link{CircadianAnalysisOptions}}
#'       \item \code{\link{CircadianBootstrapOptions}}
#'     }
#'   }
#'   \item{Run power analysis}{
#'     \itemize{
#'       \item \code{\link{runSimsSingleCohort}}
#'       \item \code{\link{runDifferentialPower}}
#'       \item \code{\link{runBootstrapDesignGrid}}
#'       \item \code{\link{npower}}
#'       \item \code{\link{circaPowerApproxN80}}
#'       \item \code{\link{recommendDesign}}
#'     }
#'   }
#'   \item{Plot results}{
#'     \itemize{
#'       \item \code{\link{plotSingleCohortPower}}
#'       \item \code{\link{plotDiffPower}}
#'       \item \code{\link{plotBootstrapDesignGrid}}
#'     }
#'   }
#'   \item{Detect rhythmicity}{
#'     \itemize{
#'       \item \code{\link{detect_cosinor}}
#'     }
#'   }
#'   \item{Simulate and helpers}{
#'     \itemize{
#'       \item \code{\link{simCircadianSingleCohort2H}}
#'       \item \code{\link{makeAdaptiveRStrata}}
#'     }
#'   }
#'   \item{App}{
#'     \itemize{
#'       \item \code{\link{launchShiny}}
#'     }
#'   }
#' }
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
