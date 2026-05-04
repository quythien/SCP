#' CircadianPower Package Setup
#' 
#' @description Load all CircadianPower functions in correct order.
#' This framework provides semiparametric power analysis for circadian
#' rhythm studie
#' SEMIPARAMETRIC FEATURES:
#' - Parametric cosinor model (y = M + A·cos(ωt - φ) + ε)
#' - Empirical parameter distributions from pilot data
#' - Empirical Bayes variance estimation
#' - No assumptions on amplitude, phase, or noise distributions
#'
#' @author Thien Pham
#' @date 2025-02-12


required_packages = c(
  "ggplot2",      # Plotting
  "reshape2",     # Data reshaping for heatmaps
  "limma",        # Linear models and empirical Bayes
  "minpack.lm",   # Levenberg-Marquardt optimization
  "nloptr",       # Non-linear optimization for LR tests
  "parallel"      # Parallel processing
)

optional_packages = c(
  "Matrix",        # Matrix operations (for some DCP functions)
  "MetaCycle",     # For JTK_CYCLE rhythmicity detection (detect_JTK)
  "rain",          # For RAIN rhythmicity detection (detect_RAIN)
  "circacompare",  # For CircaCompare differential analysis (detect_CircaCompare)
  "limorhyde",     # For LimoRhyde differential rhythmicity (detect_LimoRhyde)
  "FMM",           # For FMM waveform simulation (simCircadianFMM)
  "tidyr",         # For plotFMMDifferential() (pivot_longer)
  "patchwork",     # For plotFMMViolation() multi-panel layout
  "cowplot",       # For plotFMMViolation() draw wrapper
  "gtable"         # For plotFMMViolation() grob manipulation
)

cat("Checking required packages...\n")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf(
      "Required package '%s' is not installed.\nPlease install it with: install.packages('%s')",
      pkg, pkg
    ))
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}
cat("✓ All required packages loaded\n\n")



# Order matters - later files depend on earlier ones
source_files = c(
  "options.R",                   # Options and parameter setup
  "simulation.R",                # Core simulation (simCircadianDiff, simCircadianSingleCohort)
  "detection.R",                 # Detection methods and bridge functions (DCP pipeline)
  "estimation.R",                # Empirical parameter estimation + one_cosinor_OLS
  "runner.R",                    # Main simulation loops
  "utils.R",                     # Plotting and reporting
  "plot_dr_power.R",             # Stratified power plotting functions
  "plot_with_se.R",              # SE-bar helpers + DR/DP 6-panel figures
  "plot_single_cohort.R",        # Single-cohort Figure 1 (3-panel)
  "plot_diff.R",                 # Differential power plot (18-panel)
  "plot_bvsm.R",                 # B vs m trade-off plot (Figure 3)
  "plot_fmm.R",                  # FMM waveform robustness plots (Figures 4, 4b)
  "npower.R",                    # npower(): find n for target power at given FDR
  "LRTest_diff_phase.R",         # LR test for differential phase
  "LRTest_diff_amp.R",           # LR test for differential amplitude
  "bootstrap_sim.R",             # Bootstrap design grid
  "fourier_sim.R",               # Fourier deviation simulation
  "design_comparison.R",         # Two-stage vs bootstrap comparison
  "summarize_dcp_pairs.R"        # Cross-pair DCP summary table (prelim supp table)
)

for (file in source_files) {
  if (file.exists(file)) {
    cat(sprintf("  ✓ %s\n", file))
    source(file)
  } else {
    warning(sprintf("  ✗ %s NOT FOUND\n", file))
  }
}

# C++ acceleration (Rcpp/RcppArmadillo) — optional, graceful fallback to R
.CPP_LOADED <- FALSE
if (file.exists("src/cosinor_fast.cpp") &&
    requireNamespace("Rcpp", quietly = TRUE) &&
    requireNamespace("RcppArmadillo", quietly = TRUE)) {
  tryCatch({
    suppressMessages(Rcpp::sourceCpp("src/cosinor_fast.cpp"))
    source("src/cosinor_fast.R")
    .CPP_LOADED <- TRUE
    cat("  ✓ C++ acceleration loaded (cosinor_fast)\n")
  }, error = function(e) {
    cat(sprintf("  ✗ C++ acceleration unavailable: %s\n", conditionMessage(e)))
  })
} else {
  cat("  ✗ C++ acceleration skipped (src/cosinor_fast.cpp not found or Rcpp unavailable)\n")
}
