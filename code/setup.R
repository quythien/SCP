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
  "MetaCycle",     # For JTK_CYCLE rhythmicity detection
  "rain",          # For RAIN rhythmicity detection
  "circacompare"   # For CircaCompare differential analysis
)

cat("Checking required packages...\n")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  Installing %s...\n", pkg))
    install.packages(pkg, quiet = TRUE)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}
cat("✓ All required packages loaded\n\n")



# Order matters - later files depend on earlier ones
source_files = c(
  "options.R",                   # Options and parameter setup
  "simulation.R",                # Core simulation (simCircadian, simCircadianDiff)
  "detection.R",                 # Detection methods and bridge functions (DCP pipeline)
  "estimation.R",                # Empirical parameter estimation from pilot data
  "runner.R",                    # Main simulation loops
  "power.R",                     # Power computation and comparison
  "utils.R",                     # Plotting and reporting
  "plot_dr_power.R",             # Stratified power plotting functions
  "plot_with_se.R",              # SE-bar helpers + DR/DP 6-panel figures
  "plot_single_cohort.R",        # Single-cohort Figure 1 (3-panel)
  "plot_diff.R",                 # Differential power plot (18-panel)
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
