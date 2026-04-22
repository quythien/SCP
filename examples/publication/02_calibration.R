#' =======================================================================
#' run_calibration.R — Bootstrap vs Two-Stage Method Comparison
#' =======================================================================
#'
#' PURPOSE
#'   Validates the two power-estimation approaches (two-stage and bootstrap)
#'   using fully synthetic data with known true parameters.
#'
#'   Scientific question:
#'     "Given known true biology, do two-stage and bootstrap agree on the
#'      central power estimate, and does bootstrap additionally quantify
#'      the uncertainty from finite pilot data?"
#'
#'   Both methods use the same synthetic pilot drawn from known opts_bio.
#'   The ONLY difference between them is whether parameter estimation
#'   uncertainty is propagated:
#'     Two-stage:  estimate params once → single power curve (no CI)
#'     Bootstrap:  resample params nboot times → power ± CI
#'
#'   Expected result:
#'     - Central power estimates are similar between methods
#'     - Bootstrap CI is strictly wider than the two-stage point estimate
#'     - CI width decreases as n_pilot increases (less estimation uncertainty)
#'
#' WHAT THIS FILE IS NOT
#'   - Not a check of Type I error / FDR (→ run_validation.R)
#'   - Not a real-data power recommendation (→ run_pipeline.R, run_pipeline2.R)
#'   - Not a design optimisation (→ run_pipeline2.R Section 3)
#'
#' OUTPUTS
#'   output/calibration_<timestamp>/
#'     figures/design_comparison.pdf    (two-stage vs bootstrap power curves)
#'     comparison.rds
#'     calibration_summary.txt
#'
#' USAGE
#'   Rscript examples/run_calibration.R
#'
#' RUNTIME
#'   ~15-30 min at production settings (ngenes=3000, nsims=30, nboot=30).
#'   Reduce NSIMS/NBOOT for a quick check; results are noisier but structure
#'   is preserved.
#'
#' @author Thien Pham

# =====================================================================
# SETTINGS — edit here for quick checks vs production runs
# Smoke test: POWERSIM_SMOKE=1 Rscript examples/publication/02_calibration.R
# =====================================================================

SMOKE <- nzchar(Sys.getenv("POWERSIM_SMOKE"))
NGENES      <- if (SMOKE) 500L  else 3000L
NSIMS       <- if (SMOKE) 5L    else 30L
NBOOT       <- if (SMOKE) 5L    else 30L
NSIMS_INNER <- if (SMOKE) 5L    else 20L
N_PILOT     <- if (SMOKE) 10L   else 30L

# N grid: all divisible by LCM(4) = 4; covers a biologically relevant range
N_GRID      <- c(48L, 72L, 144L)

# Fixed single-B design: B=4 (every 6h). Single B ensures CI reflects
# parameter estimation uncertainty only, not design variation.
S_B         <- 4L

# Pilot TOD distribution (post-mortem spacing; 12 evenly-spaced bins)
PILOT_TIMES <- seq(0, 22, by = 2)
DESIGN_VEC  <- PILOT_TIMES


# =====================================================================
# SETUP
# =====================================================================

POWERSIM_ROOT <- {
  env <- Sys.getenv("POWERSIM_ROOT", unset = "")
  if (nzchar(env)) env else
    "/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim"
}
setwd(POWERSIM_ROOT)
old_wd <- setwd(file.path(getwd(), "code"))
source("setup.R")
setwd(old_wd)

run_tag  <- format(Sys.time(), "%Y%m%d_%H%M")
base_out <- file.path("output", paste0("calibration_", run_tag))
fig_dir  <- file.path(base_out, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

cat("=======================================================================\n")
cat("Bootstrap vs Two-Stage Method Comparison\n")
cat("=======================================================================\n\n")
cat(sprintf("Settings: ngenes=%d, nsims=%d, nboot=%d, nsims_inner=%d\n",
            NGENES, NSIMS, NBOOT, NSIMS_INNER))
cat(sprintf("Pilot:    n=%d synthetic subjects (known true parameters)\n", N_PILOT))
cat(sprintf("Design:   fixed B=%d (every 6h), N grid = %s\n",
            S_B, paste(N_GRID, collapse = ", ")))
cat(sprintf("Output:   %s/\n\n", base_out))

t_start <- proc.time()


# =====================================================================
# SHARED OPTIONS
# =====================================================================
# Realistic circadian scenario: 25% rhythmic, 15% with differential
# rhythmicity (DR). Amplitude distribution from BA11/BA47 empirical
# lookup. No DP or DA — isolates the DR test for a clean comparison.

opts_bio <- CircadianBioOptions(
  ngenes        = NGENES,
  prop_rhythmic = 0.25,
  prop_DR       = 0.15,
  prop_DP       = 0.00,
  phase_diff    = c(0, 0),
  amp_diff      = c(1, 1)
)

opts_analysis <- CircadianAnalysisOptions(reference_n = 72L)

cat("True (known) bio options:\n")
print(opts_bio)
cat("\n")


# =====================================================================
# GENERATE SYNTHETIC PILOT FROM KNOWN TRUTH
# =====================================================================

cat(sprintf("Generating synthetic pilot (n=%d) from known opts_bio...\n\n", N_PILOT))
pilot       <- generatePilotData(opts_bio, N_PILOT, PILOT_TIMES, seed = 42L)
pilot_data  <- pilot$data
pilot_times <- pilot$times


# =====================================================================
# TWO-STAGE: estimate once → single power curve
# =====================================================================

cat("Running two-stage (point estimate)...\n")

opts_design <- CircadianDesignOptions(
  sample_sizes = N_GRID,
  nsims        = NSIMS,
  design       = "active"
)

two_stage <- runTwoStagePower(
  pilot_data    = pilot_data,
  pilot_times   = pilot_times,
  design.opts   = opts_design,
  analysis.opts = opts_analysis,
  bio_diff.opts = opts_bio,
  test_type     = "DR",
  verbose       = TRUE
)


# =====================================================================
# BOOTSTRAP: resample params → power ± CI
# =====================================================================

cat(sprintf("\nRunning bootstrap (B=%d, %d draws × %d inner sims)...\n",
            S_B, NBOOT, NSIMS_INNER))

boot.opts <- CircadianBootstrapOptions(
  design_vector = DESIGN_VEC,
  B_values      = S_B,
  N_values      = N_GRID,
  nboot         = NBOOT,
  nsims_inner   = NSIMS_INNER,
  design        = "active",
  seed          = 42L
)

boot_result <- runBootstrapDesignGrid(
  pilot_data    = pilot_data,
  pilot_times   = pilot_times,
  boot.opts     = boot.opts,
  analysis.opts = opts_analysis,
  bio_diff.opts = opts_bio,
  verbose       = TRUE
)


# =====================================================================
# COMPARISON
# =====================================================================

comparison <- compareDesignApproaches(
  two_stage_result = two_stage,
  bootstrap_result = boot_result,
  test_type        = "DR",
  target_power     = 0.80
)

saveRDS(list(two_stage  = two_stage,
             boot       = boot_result,
             comparison = comparison),
        file.path(base_out, "comparison.rds"))

comp_df   <- comparison$comparison
ci_widths <- comp_df$boot_ci_hi - comp_df$boot_ci_lo

cat("\n--- Results ---\n")
cat(sprintf("  %-6s  %-14s  %-22s\n", "N", "Two-stage", "Bootstrap [95% CI]"))
cat(paste0(rep("-", 50), collapse = ""), "\n")
for (j in seq_len(nrow(comp_df))) {
  cat(sprintf("  N=%-4d  %5.1f%%         %5.1f%% [%5.1f%%, %5.1f%%]\n",
              comp_df$n[j],
              100 * comp_df$two_stage_power[j],
              100 * comp_df$boot_power_mean[j],
              100 * comp_df$boot_ci_lo[j],
              100 * comp_df$boot_ci_hi[j]))
}
cat(paste0(rep("-", 50), collapse = ""), "\n")
cat(sprintf("  Two-stage n80:         %s\n",
            ifelse(is.na(comparison$n80_two_stage), ">max(N)", comparison$n80_two_stage)))
cat(sprintf("  Bootstrap n80 median:  %s\n",
            ifelse(is.na(comparison$n80_boot_median), ">max(N)", comparison$n80_boot_median)))
cat(sprintf("  Bootstrap n80 95%% CI:  [%s, %s]\n",
            ifelse(is.na(comparison$n80_boot_lo), "NA", comparison$n80_boot_lo),
            ifelse(is.na(comparison$n80_boot_hi), "NA", comparison$n80_boot_hi)))
cat(sprintf("  Bootstrap CI widths:   %s\n",
            paste(sprintf("%.3f", ci_widths), collapse = ", ")))
cat(sprintf("  CI > 0 (uncertainty propagated): %s\n",
            ifelse(any(ci_widths > 0, na.rm = TRUE), "YES", "NO")))

plotDesignComparison(
  comparison,
  target_power = 0.80,
  panels       = "A",      # Panel B (n80 bar) not needed for calibration
  output_file  = file.path(fig_dir, "design_comparison.pdf")
)
cat(sprintf("\nFigure: %s\n", file.path(fig_dir, "design_comparison.pdf")))


# =====================================================================
# WRAP-UP
# =====================================================================

t_elapsed <- proc.time() - t_start

summary_lines <- c(
  "Bootstrap vs Two-Stage Method Comparison",
  sprintf("Run: %s", run_tag),
  sprintf("Settings: ngenes=%d, nsims=%d, nboot=%d, nsims_inner=%d, n_pilot=%d",
          NGENES, NSIMS, NBOOT, NSIMS_INNER, N_PILOT),
  sprintf("N grid: %s", paste(N_GRID, collapse = ", ")),
  "",
  sprintf("Two-stage n80:         %s",
          ifelse(is.na(comparison$n80_two_stage), ">max(N)", comparison$n80_two_stage)),
  sprintf("Bootstrap n80 median:  %s",
          ifelse(is.na(comparison$n80_boot_median), ">max(N)", comparison$n80_boot_median)),
  sprintf("Bootstrap n80 95%% CI:  [%s, %s]",
          ifelse(is.na(comparison$n80_boot_lo), "NA", comparison$n80_boot_lo),
          ifelse(is.na(comparison$n80_boot_hi), "NA", comparison$n80_boot_hi)),
  sprintf("CI width > 0:          %s",
          ifelse(any(ci_widths > 0, na.rm = TRUE), "YES", "NO")),
  "",
  sprintf("Runtime: %.1f min", t_elapsed[3] / 60),
  sprintf("Output:  %s/", base_out)
)

writeLines(summary_lines, file.path(base_out, "calibration_summary.txt"))

cat("\n=======================================================================\n")
cat("CALIBRATION COMPLETE\n")
cat("=======================================================================\n\n")
writeLines(summary_lines)
cat("\n")
