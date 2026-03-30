#' Example Workflow for CircadianPower
#'
#' Scenarios 1-4 : Core functionality (DR/DP power, legacy API, pilot calibration)
#' Scenario  5   : Bootstrap B vs m design grid
#' Scenario  6   : Fourier deviation — power under waveform misspecification
#' Scenario  7   : Two-stage vs bootstrap design comparison
#' Scenario  8   : Ground truth calibration (known truth vs two approaches)
#'
#' Smoke-test settings are defined once at the top.
#' Increase NSIMS, NBOOT, NGENES for production runs.

source("setup.R")

# ==============================================================================
# SETTINGS  (Tier 1 smoke / Tier 2 validation / Tier 3 figures)
#
# Tier 1 (every edit):   ngenes=200,  nsims=3,   nboot=5,  nsims_inner=3,  M={6,9}
# Tier 2 (pre-merge):    ngenes=1000, nsims=15,  nboot=10, nsims_inner=10, M={6,9,12,18}
# Tier 3 (final figs):   ngenes=3000, nsims=30,  nboot=30, nsims_inner=20, M={6,9,12,18,24,36}
#
# N grid note: M_VALUES extended in Tier 3 to reach N=192-288 so that n80 is
# resolvable.  Tier 2 confirmed power only reaches 68% at N=144; n80 is NA.
#
# S8 calibration note: Both two-stage and bootstrap overestimate power by ~20 pp
# relative to the oracle across all N.  Root cause: estCircadianParam filters
# rhythmic genes at p<0.1 (winner's curse), inflating amplitude estimates 1.38x
# and effect size 1.44x.  Coverage = 0% is the expected and correct result for
# these methods — it documents a systematic limitation, not a code error.
# Increasing NSIMS will NOT fix this; it requires bias-corrected amplitude
# estimation (future work).
# ==============================================================================

NGENES      <- 3000         # genes per simulation
NSIMS       <- 30           # simulations per sample size
NBOOT       <- 30           # bootstrap parameter draws (outer loop)
NSIMS_INNER <- 20           # simulations per bootstrap draw (inner loop)
N_PILOT     <- 30           # synthetic / real pilot subjects
B_VALUES    <- c(4, 6, 8, 12)    # time-point counts: every 6h / 4h / 3h / 2h
N_GRID      <- c(48L, 72L, 144L) # N per group; all divisible by LCM(B_VALUES)=24
# m per B: B=4→12,18,36  B=6→8,12,24  B=8→6,9,18  B=12→4,6,12
PILOT_TIMES <- seq(0, 22, by = 2)   # 12 evenly-spaced TOD bins (post-mortem)
DESIGN_VEC  <- PILOT_TIMES          # same grid used for active design sweep

OUTPUT_DIR  <- "../output/tier3"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("\nTier 3 settings: ngenes=%d, nsims=%d, nboot=%d, N_grid=%s\n",
            NGENES, NSIMS, NBOOT, paste(N_GRID, collapse = ",")))
cat(sprintf("Output directory: %s\n\n", OUTPUT_DIR))


# ==============================================================================
# SHARED ANALYSIS OPTIONS  (used by all scenarios)
# ==============================================================================

opts_analysis <- CircadianAnalysisOptions(
  reference_n = 96      # reference N for Fourier heatmap panel (mid-grid for Tier 3)
)


# ==============================================================================
# SCENARIO 1: Quick DR Power with Config Objects
# ==============================================================================

cat("\n=== SCENARIO 1: Quick DR Power Analysis ===\n\n")

opts_bio <- CircadianBioOptions(
  ngenes        = NGENES,
  prop_rhythmic = 0.25,
  prop_DR       = 0.15,
  prop_DP       = 0.00,
  prop_DA       = 0.00,
  phase_diff    = c(0, 0),
  amp_diff      = c(1, 1)
)

opts_design <- CircadianDesignOptions(
  sample_sizes = c(24, 48),
  nsims        = NSIMS,
  design       = "active"
)

print(opts_bio)
print(opts_design)

dr_results <- runSimsDiff(opts_bio, opts_design, opts_analysis)
cat(sprintf("DR done. Genes=%d, N=%s, sims=%d\n",
            dr_results$ngenes,
            paste(dr_results$sample_sizes, collapse = ","),
            dr_results$nsims))


# ==============================================================================
# SCENARIO 2: DP Power with Phase Shifts
# ==============================================================================

cat("\n=== SCENARIO 2: Differential Phase (DP) ===\n\n")

opts_bio_DP <- updateBioOptions(opts_bio,
  prop_DR    = 0.00,
  prop_DP    = 0.15,
  prop_DA    = 0.00,
  phase_diff = c(-6, 6),
  amp_diff   = c(0.5, 2)
)

dp_results <- runSimsDiff(opts_bio_DP, opts_design, opts_analysis)
cat("DP done.\n")


# ==============================================================================
# SCENARIO 3: Legacy Flat-Arg Signature (Backward Compatible)
# ==============================================================================

cat("\n=== SCENARIO 3: Legacy API (backward compatible) ===\n\n")

legacy_results <- runSimsDiff(
  sample_sizes  = c(24, 48),
  nsims         = NSIMS,
  ngenes        = NGENES,
  prop_rhythmic = 0.25,
  prop_DR       = 0.10,
  prop_DP       = 0.10,
  prop_DA       = 0.00,
  phase_diff    = c(-4, 4),
  design        = "active",
  verbose       = FALSE
)
cat(sprintf("Legacy API OK. dims: [%s]\n",
            paste(dim(legacy_results$pval_DR), collapse = " x ")))


# ==============================================================================
# SCENARIO 4: Pilot Calibration via estCircadianParam()
# ==============================================================================

cat("\n=== SCENARIO 4: Calibration from Pilot Data ===\n\n")

opts_from_pilot <- estCircadianParam(
  data    = generatePilotData(opts_bio, N_PILOT, PILOT_TIMES, seed = 1)$data,
  times   = generatePilotData(opts_bio, N_PILOT, PILOT_TIMES, seed = 1)$times,
  period  = 24,
  prop_DR = 0.10,
  prop_DP = 0.10,
  prop_DA = 0.00,
  verbose = TRUE
)
cat("\nEstimated bio options from pilot:\n")
print(opts_from_pilot)


# ==============================================================================
# SCENARIO 5: Bootstrap Design Grid (B vs m tradeoff)
# ==============================================================================
#
# For fixed total N = B × m, which split is best?
# Bootstrap gives power ± CI to show uncertainty in the recommendation.

cat("\n=== SCENARIO 5: Bootstrap Design Grid ===\n\n")

# For a clean (N, B) grid, every N must be divisible by every B.
# LCM(4, 8) = 8, so restrict N_GRID to multiples of 8.
# N = 36 (from 4×9) is not divisible by B = 8 → excluded here to avoid
# an infeasible cell and an incomplete heatmap.
S5_N_VALUES <- as.integer(N_GRID[N_GRID %% 8L == 0L])   # multiples of LCM(B_VALUES)

boot.opts <- CircadianBootstrapOptions(
  design_vector = DESIGN_VEC,     # max 12 time points; B selects a subset
  B_values      = B_VALUES,       # c(4, 8)
  N_values      = S5_N_VALUES,    # explicit multiples of LCM(B_VALUES) — full grid
  nboot         = NBOOT,
  nsims_inner   = NSIMS_INNER,
  design        = "active"
)

print(boot.opts)

# Generate pilot data from known opts_bio (consistent gene count)
pilot       <- generatePilotData(opts_bio, N_PILOT, PILOT_TIMES, seed = 42)
pilot_data  <- pilot$data
pilot_times <- pilot$times

boot_grid <- runBootstrapDesignGrid(
  pilot_data    = pilot_data,
  pilot_times   = pilot_times,
  boot.opts     = boot.opts,
  analysis.opts = opts_analysis,
  bio_diff.opts = opts_bio
)

saveRDS(boot_grid, file.path(OUTPUT_DIR, "s5_boot_grid.rds"))
cat("\n--- Bootstrap Grid Summary (DR) ---\n")
summaryBootstrapDesignGrid(boot_grid, test_type = "DR")

cat("\nPlotting bootstrap design grid...\n")
plotBootstrapDesignGrid(boot_grid, test_type = "DR",
                        output_file = file.path(OUTPUT_DIR, "s5_bootstrap_design_grid.pdf"))

# Verification
ci_width <- boot_grid$power_ci_hi - boot_grid$power_ci_lo
cat(sprintf("\nVerification — any CI width > 0: %s\n", any(ci_width > 0, na.rm = TRUE)))
ref_n <- boot_grid$N_values[ceiling(length(boot_grid$N_values) / 2)]  # middle N
cat(sprintf("Optimal B at N=%d: B = %d\n",
            ref_n, boot_grid$optimal_B[boot_grid$N_values == ref_n]))


# ==============================================================================
# SCENARIO 6: Fourier Deviation — Power Under Waveform Misspecification
# ==============================================================================
#
# Measures whether greater temporal coverage (larger B) helps when the true
# waveform deviates from a pure cosinor (adds 2nd/3rd harmonics) while the
# analysis model remains a fundamental-only cosinor.
#
# What this scenario shows:
#   - At fixed total N, whether designs with more time coverage (larger B,
#     smaller m) lose less power under harmonic misspecification.
#   - Whether the preferred B changes as the true waveform becomes less cosinor-like.

cat("\n=== SCENARIO 6: Fourier Deviation ===\n\n")

# Fix total N so B and N are separated. Use the reference N if it is feasible
# for all candidate B values; otherwise use the largest common feasible N.
S6_N_FIXED <- if (all(opts_analysis$reference_n %% B_VALUES == 0)) {
  opts_analysis$reference_n
} else {
  max(N_GRID[sapply(N_GRID, function(n) all(n %% B_VALUES == 0))])
}

cat(sprintf("Fixed-N design comparison at N=%d\n", S6_N_FIXED))
cat(sprintf("Candidate B values: %s\n\n", paste(B_VALUES, collapse = ", ")))

harmonic_grid <- expand.grid(
  alpha2 = c(0, 0.25, 0.5),    # pure cosinor / mild 2nd / strong 2nd harmonic
  alpha3 = 0                    # no 3rd harmonic (keep grid minimal)
)

s6_results <- vector("list", length(B_VALUES))
names(s6_results) <- paste0("B", B_VALUES)

for (i in seq_along(B_VALUES)) {
  B <- B_VALUES[i]
  m <- as.integer(S6_N_FIXED / B)
  time_pts <- .selectTimePoints(DESIGN_VEC, B)
  cts_fixed <- rep(time_pts, each = m)

  opts_design_fourier <- CircadianDesignOptions(
    sample_sizes = S6_N_FIXED,
    nsims        = NSIMS,
    design       = "active",
    cts          = cts_fixed
  )

  s6_results[[i]] <- runFourierDeviationPower(
    bio.opts      = opts_bio,
    design.opts   = opts_design_fourier,
    analysis.opts = opts_analysis,
    harmonic_grid = harmonic_grid,
    test_type     = "DR"
  )
}

# Assemble [harmonic x B] matrix of mean power at fixed N
s6_power_mat <- sapply(s6_results, function(res) res$power_mean[, 1])
if (is.null(dim(s6_power_mat))) {
  s6_power_mat <- matrix(s6_power_mat, ncol = length(s6_results))
}
rownames(s6_power_mat) <- s6_results[[1]]$harmonic_grid$label
colnames(s6_power_mat) <- paste0("B=", B_VALUES)

s6_out <- list(
  N_fixed       = S6_N_FIXED,
  B_values      = B_VALUES,
  harmonic_grid = harmonic_grid,
  power_mean    = s6_power_mat,
  by_design     = s6_results
)
saveRDS(s6_out, file.path(OUTPUT_DIR, "s6_fourier_result.rds"))

# Single-panel plot: raw DR power vs harmonic level, one line per B value
s6_pdf <- file.path(OUTPUT_DIR, "s6_fourier_deviation.pdf")
pdf(s6_pdf, width = 7, height = 5)
on.exit({ dev.off() }, add = TRUE)

pm_pct   <- 100 * s6_power_mat
pure_idx <- which(harmonic_grid$alpha2 == 0 & harmonic_grid$alpha3 == 0)

b_cols <- c("steelblue", "seagreen", "darkorange", "firebrick")
matplot(seq_len(nrow(pm_pct)), pm_pct,
        type = "b", pch = 19, lty = 1, lwd = 2,
        col = b_cols[seq_len(ncol(pm_pct))],
        xaxt = "n", xlab = "2nd harmonic amplitude (alpha2)",
        ylab = "DR Power (%)",
        main = sprintf("DR Power vs Waveform Misspecification\n(Fixed N=%d, active design)", S6_N_FIXED),
        las = 1, ylim = c(0, max(pm_pct, na.rm = TRUE) * 1.15))
axis(1, at = seq_len(nrow(pm_pct)), labels = rownames(pm_pct), las = 2, cex.axis = 0.9)
abline(h = 80, lty = 2, col = "gray40")
legend("topright", legend = colnames(pm_pct),
       col = b_cols[seq_len(ncol(pm_pct))],
       pch = 19, lty = 1, lwd = 2, bty = "n", cex = 0.9,
       title = "Design (B time points)")

cat("\nPlotting Fourier deviation results...\n")
cat(sprintf("Saved: %s\n", s6_pdf))

# Verification: does higher B lose less power under strong harmonics?
harm_idx <- which(harmonic_grid$alpha2 >= 0.5)
if (length(pure_idx) == 1 && length(harm_idx) > 0) {
  pure_by_B <- pm_pct[pure_idx, ]
  harm_by_B <- colMeans(pm_pct[harm_idx, , drop = FALSE], na.rm = TRUE)
  loss_by_B <- pure_by_B - harm_by_B
  cat(sprintf("\nMean pure-cosinor power at N=%d:  %s\n",
              S6_N_FIXED,
              paste(sprintf("%s %.1f%%", names(pure_by_B), pure_by_B), collapse = " | ")))
  cat(sprintf("Mean high-harmonic power at N=%d: %s\n",
              S6_N_FIXED,
              paste(sprintf("%s %.1f%%", names(harm_by_B), harm_by_B), collapse = " | ")))
  cat(sprintf("Absolute power loss (pp):         %s\n",
              paste(sprintf("%s %.1f", names(loss_by_B), loss_by_B), collapse = " | ")))
}


# ==============================================================================
# SCENARIO 7: Two-Stage vs Bootstrap Design Comparison
# ==============================================================================
#
# Clean comparison of two approaches using the SAME pilot data and SAME fixed
# design (B = min(B_VALUES) = 4, m varying to cover N_GRID).
#
# Two-stage:  estimate params from pilot once  → one power curve (no CI).
# Bootstrap:  resample params nboot times       → power ± CI.
#
# The ONLY difference is whether parameter estimation uncertainty is propagated.
# Using a single fixed B ensures CI bands reflect parameter uncertainty only
# (not design variation), giving a clean apples-to-apples comparison.

cat("\n=== SCENARIO 7: Two-Stage vs Bootstrap Comparison ===\n\n")

# Fixed B = min(B_VALUES) = 4; m values chosen so N = B × m covers N_GRID exactly.
# With B = 4 and N_GRID = c(48,72,144): m = N/4 = c(12,18,36)
# (single B used here so CI reflects parameter uncertainty only, not design variation)
S7_B      <- min(B_VALUES)                    # fixed B for this comparison (= 4)
S7_M_VALS <- as.integer(N_GRID / S7_B)       # m values cover N_GRID exactly

cat(sprintf("S7 fixed design: B = %d time points, m = %s → N = %s\n",
            S7_B,
            paste(S7_M_VALS, collapse = ","),
            paste(S7_B * S7_M_VALS, collapse = ",")))

opts_design_compare <- CircadianDesignOptions(
  sample_sizes = N_GRID,
  nsims        = NSIMS,
  design       = "active"
)

# --- Two-stage: estimate once from pilot → single power curve ---
two_stage_result <- runTwoStagePower(
  pilot_data    = pilot_data,
  pilot_times   = pilot_times,
  design.opts   = opts_design_compare,
  analysis.opts = opts_analysis,
  bio_diff.opts = opts_bio,
  test_type     = "DR"
)

# --- Bootstrap: dedicated run with fixed B, parameter uncertainty only ---
boot.opts_s7 <- CircadianBootstrapOptions(
  design_vector = DESIGN_VEC,
  B_values      = S7_B,        # single fixed B → no design variation
  m_values      = S7_M_VALS,   # covers N_GRID exactly
  nboot         = NBOOT,
  nsims_inner   = NSIMS_INNER,
  design        = "active"
)

boot_result_s7 <- runBootstrapDesignGrid(
  pilot_data    = pilot_data,
  pilot_times   = pilot_times,
  boot.opts     = boot.opts_s7,
  analysis.opts = opts_analysis,
  bio_diff.opts = opts_bio
)

comparison <- compareDesignApproaches(
  two_stage_result = two_stage_result,
  bootstrap_result = boot_result_s7,
  test_type        = "DR",
  target_power     = 0.80
)

saveRDS(list(two_stage   = two_stage_result,
             boot_s7     = boot_result_s7,
             comparison  = comparison),
        file.path(OUTPUT_DIR, "s7_comparison.rds"))

cat("\n--- Comparison Summary ---\n")
cat(sprintf("  Two-stage n80:        %s\n",  comparison$n80_two_stage))
cat(sprintf("  Bootstrap n80 median: %s\n",  comparison$n80_boot_median))
cat(sprintf("  Bootstrap n80 2.5%%:   %s\n", comparison$n80_boot_lo))
cat(sprintf("  Bootstrap n80 97.5%%:  %s\n", comparison$n80_boot_hi))

cat("\nPlotting comparison figure...\n")
plotDesignComparison(comparison, target_power = 0.80,
                     output_file = file.path(OUTPUT_DIR, "s7_design_comparison.pdf"))


# ==============================================================================
# SCENARIO 8: Ground Truth Calibration
# ==============================================================================
#
# When the true parameters are KNOWN (we generated the data), how well does
# each approach recover the true power?
# Oracle: simulates with exact known params (no estimation).
# Two-stage: estimates from synthetic pilot → power curve.
# Bootstrap: resamples pilot params → power ± CI.
# Key check: does bootstrap 95% CI cover the oracle truth?

cat("\n=== SCENARIO 8: Ground Truth Calibration ===\n\n")

opts_design_gt <- CircadianDesignOptions(
  sample_sizes = N_GRID,    # 24, 36, 48, 72
  nsims        = NSIMS,
  design       = "active"
)

gt_result <- runGroundTruthComparison(
  true_bio.opts = opts_bio,        # ground truth (known)
  design.opts   = opts_design_gt,
  analysis.opts = opts_analysis,
  n_pilot       = 120,              # larger pilot reduces winner's curse in amplitude estimation
  pilot_times   = PILOT_TIMES,     # TOD distribution
  nboot         = NBOOT,
  nsims_inner   = NSIMS_INNER,
  test_type     = "DR",
  seed          = 42
)

saveRDS(gt_result, file.path(OUTPUT_DIR, "s8_gt_result.rds"))
cat("\n--- Ground Truth Summary ---\n")
for (j in seq_along(gt_result$sample_sizes)) {
  cat(sprintf(
    "  N=%d: oracle=%.1f%%, two-stage=%.1f%%, boot=%.1f%% [%.1f%%, %.1f%%]\n",
    gt_result$sample_sizes[j],
    100 * gt_result$true_power_mean[j],
    100 * gt_result$ts_power_mean[j],
    100 * gt_result$boot_power_mean[j],
    100 * gt_result$boot_ci_lo[j],
    100 * gt_result$boot_ci_hi[j]
  ))
}
cat(sprintf("  Bootstrap CI coverage: %.0f%%\n", 100 * gt_result$coverage))

cat("\nPlotting ground truth calibration figure...\n")
plotGroundTruthComparison(gt_result, target_power = 0.80,
                          output_file = file.path(OUTPUT_DIR, "s8_ground_truth_calibration.pdf"))

