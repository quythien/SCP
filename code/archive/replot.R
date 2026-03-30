#' Replot saved simulation results
#'
#' Load .rds files from a previous run and regenerate all figures.
#' Edit the plot function calls below to change aesthetics without
#' re-running any simulations.
#'
#' Usage:
#'   source("setup.R")        # loads all functions
#'   RDS_DIR <- "../output/smoke_test"   # or your production output folder
#'   source("replot.R")

source("setup.R")

# ── Point this to the folder containing the .rds files ──────────────────────
RDS_DIR <- "../output/smoke_test"   # change to production run folder as needed
PDF_DIR <- RDS_DIR                  # save new PDFs alongside the rds files

cat(sprintf("Loading results from: %s\n", RDS_DIR))

# ── Load ─────────────────────────────────────────────────────────────────────

s5 <- readRDS(file.path(RDS_DIR, "s5_boot_grid.rds"))
s6 <- readRDS(file.path(RDS_DIR, "s6_fourier_result.rds"))
s7 <- readRDS(file.path(RDS_DIR, "s7_comparison.rds"))
s8 <- readRDS(file.path(RDS_DIR, "s8_gt_result.rds"))

cat("All results loaded.\n\n")

# ── Scenario 5: Bootstrap Design Grid ────────────────────────────────────────

cat("Plotting S5: Bootstrap design grid...\n")
summaryBootstrapDesignGrid(s5, test_type = "DR")
plotBootstrapDesignGrid(s5, test_type = "DR",
                        output_file = file.path(PDF_DIR, "s5_bootstrap_design_grid.pdf"))

# ── Scenario 6: Fourier Deviation ────────────────────────────────────────────

cat("Plotting S6: Fourier deviation...\n")
plotFourierDeviation(s6, test_type = "DR",
                     output_file = file.path(PDF_DIR, "s6_fourier_deviation.pdf"))

# ── Scenario 7: Design Comparison ────────────────────────────────────────────

cat("Plotting S7: Design comparison...\n")
plotDesignComparison(s7$comparison, target_power = 0.80,
                     output_file = file.path(PDF_DIR, "s7_design_comparison.pdf"))
# s7$boot_s7 contains the dedicated fixed-B bootstrap used in the comparison

# ── Scenario 8: Ground Truth Calibration ─────────────────────────────────────

cat("Plotting S8: Ground truth calibration...\n")
plotGroundTruthComparison(s8, target_power = 0.80,
                          output_file = file.path(PDF_DIR, "s8_ground_truth_calibration.pdf"))

cat("\nAll figures saved to:", PDF_DIR, "\n")
