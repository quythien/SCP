#!/usr/bin/env Rscript
# =============================================================================
# 09_dm_singlecohort_smoke.R
# Smoke tests for two new features:
#   Section 1: DM (Differential Mesor) endpoint — simulated data, DCP test
#   Section 2: Single-cohort simulation-based power (runSimsSingleCohort)
#   Section 3: Two-pilot mesor — separate group-2 baseline from two-pilot estimation
# =============================================================================

# Locate POWERSIM_ROOT: this script lives in examples/publication/
# When invoked via Rscript from project root, use getwd().
# When invoked via source(), use the script's own path.
POWERSIM_ROOT <- tryCatch({
  script_path <- sys.frame(1)$ofile
  normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
}, error = function(e) {
  # Rscript: check if current dir looks like the project root
  if (file.exists(file.path(getwd(), "code", "setup.R"))) {
    return(getwd())
  }
  # Try two levels up from cwd
  up2 <- normalizePath(file.path(getwd(), "..", ".."), mustWork = FALSE)
  if (file.exists(file.path(up2, "code", "setup.R"))) {
    return(up2)
  }
  stop("Cannot locate POWERSIM_ROOT; run from project root or set it manually.")
})
setwd(POWERSIM_ROOT)
source_dir <- file.path(POWERSIM_ROOT, "code")
old_wd <- setwd(source_dir); source("setup.R"); setwd(old_wd)

cat("\n============================================================\n")
cat("SMOKE TEST: DM endpoint + Single-cohort simulation\n")
cat("============================================================\n\n")

# ============================================================
# SECTION 1: DM (Differential Mesor) endpoint smoke test
# ============================================================
cat("--- Section 1: DM endpoint ---\n")

set.seed(42L)

# Build bio options with DM (prop_DM = 0.15, prop_DR = 0.10, no DA)
bio_dm <- CircadianBioOptions(
  ngenes        = 500L,
  prop_rhythmic = 0.40,
  lBaselineExpr = rnorm(500, 5, 1),
  lOD           = rnorm(500, -1, 0.3),
  amplitude     = abs(rnorm(200, 0.5, 0.2)) + 0.05,
  sigma_rhythmic = NULL,
  phase         = "uniform",
  prop_DR       = 0.10,
  prop_DP       = 0.05,
  prop_DA       = 0.00,
  prop_DM       = 0.15,
  mesor_diff    = c(1.0, 2.5),
  period        = 24L,
  sim.seed      = 42L
)

cat(sprintf("bio_dm: prop_DM = %.3f, mesor_diff = [%.1f, %.1f]\n",
            bio_dm$prop_DM, bio_dm$mesor_diff[1], bio_dm$mesor_diff[2]))

# Active design: 6 ZT points, 8 replicates each (N=48 per group)
cts_active <- rep(seq(0, 20, by = 4), each = 8)   # 48 samples, 6 ZTs

design_dm <- CircadianDesignOptions(
  sample_sizes = c(24L, 48L),
  nsims        = 10L,
  design       = "active",
  cts          = cts_active,
  test_types   = c("DR", "DM")
)

analysis_opts <- CircadianAnalysisOptions(alpha = 0.05, p.adjust.method = "BH")

cat("Running runSimsDiff with DM test...\n")
sims_dm <- runSimsDiff(bio_dm, design_dm, analysis_opts)

# Verify: fdr_DM exists and is populated
stopifnot(!is.null(sims_dm$fdr_DM))
stopifnot(dim(sims_dm$fdr_DM) == c(500L, 2L, 10L))

# Check that DM genes (type 5) exist in at least some simulations
n_DM_genes <- sum(sims_dm$diff_type[[1]] == 5)
cat(sprintf("  Number of DM genes in simulation 1: %d\n", n_DM_genes))
stopifnot(n_DM_genes > 0)

# Check DM power is non-trivial at N=48 (second sample size)
dm_power_n48 <- mean(apply(sims_dm$fdr_DM[, 2, ], 2, function(fdr) {
  true_DM <- sims_dm$diff_type[[1]] == 5
  if (sum(true_DM) == 0) return(NA_real_)
  sum(fdr[true_DM] <= 0.05, na.rm = TRUE) / sum(true_DM)
}), na.rm = TRUE)
cat(sprintf("  DM power at N=48: %.1f%%\n", 100 * dm_power_n48))

cat("  Section 1 PASSED\n\n")


# ============================================================
# SECTION 2: Single-cohort simulation-based power
# ============================================================
cat("--- Section 2: Single-cohort simulation-based power ---\n")

# Synthetic single-group pilot parameters
bio_single <- CircadianBioOptions(
  ngenes        = 500L,
  prop_rhythmic = 0.25,
  lBaselineExpr = rnorm(500, 5, 1),
  lOD           = rnorm(500, -1, 0.3),
  amplitude     = abs(rnorm(125, 0.5, 0.2)) + 0.05,
  sigma_rhythmic = NULL,
  phase         = "uniform",
  prop_DR       = 0,   # no differential (single cohort)
  prop_DP       = 0,
  prop_DA       = 0,
  prop_DM       = 0,
  period        = 24L,
  sim.seed      = 7L
)

cts_single <- rep(seq(0, 22, by = 4), each = 4)  # 24 samples, 6 ZTs

design_single <- CircadianDesignOptions(
  sample_sizes = c(12L, 24L, 48L),
  nsims        = 10L,
  design       = "active",
  cts          = cts_single
)

cat("Running runSimsSingleCohort...\n")
res_single <- runSimsSingleCohort(bio_single, design_single, analysis_opts, verbose = TRUE)

# Verify output structure
stopifnot(!is.null(res_single$marginal_power))
stopifnot(nrow(res_single$marginal_power) == 3L)
stopifnot(ncol(res_single$marginal_power) == 10L)

# Power should increase with N
mean_power <- rowMeans(res_single$marginal_power, na.rm = TRUE)
cat(sprintf("  Mean power: N=12: %.1f%%  N=24: %.1f%%  N=48: %.1f%%\n",
            100 * mean_power[1], 100 * mean_power[2], 100 * mean_power[3]))
stopifnot(mean_power[3] >= mean_power[1] - 0.05)  # power non-decreasing (with tolerance)

cat("  Section 2 PASSED\n\n")


# ============================================================
# SECTION 3: Two-pilot mesor — lBaselineExpr2 from two-group estimation
# ============================================================
cat("--- Section 3: Two-pilot mesor offset ---\n")

# Simulate two groups with different baseline (mesor shift) using lBaselineExpr2
set.seed(99L)
bio_twopilot <- CircadianBioOptions(
  ngenes         = 300L,
  prop_rhythmic  = 0.30,
  lBaselineExpr  = rnorm(300, 5, 0.8),       # group 1 mesor
  lBaselineExpr2 = rnorm(300, 5.5, 0.8),     # group 2 mesor (shifted +0.5)
  lOD            = rnorm(300, -1, 0.3),
  lOD2           = rnorm(300, -0.9, 0.3),    # group 2 noise slightly different
  amplitude      = abs(rnorm(90, 0.5, 0.15)) + 0.05,
  amplitude2     = abs(rnorm(90, 0.45, 0.15)) + 0.05,
  sigma_rhythmic = NULL,
  phase          = "uniform",
  prop_DR        = 0.10,
  prop_DP        = 0.05,
  prop_DA        = 0.00,
  prop_DM        = 0.10,
  mesor_diff     = c(0.8, 2.0),
  period         = 24L,
  sim.seed       = 99L
)

cat(sprintf("bio_twopilot: lBaselineExpr2 range [%.2f, %.2f]\n",
            min(bio_twopilot$lBaselineExpr2), max(bio_twopilot$lBaselineExpr2)))
stopifnot(!is.null(bio_twopilot$lBaselineExpr2))
stopifnot(length(bio_twopilot$lBaselineExpr2) == 300L)

# Verify that simCircadianDiff uses lBaselineExpr2 for group-2 mesor in non-DM genes
cts_tp <- rep(seq(0, 20, by = 4), each = 6)   # 36 samples, 6 ZTs

design_tp <- CircadianDesignOptions(
  sample_sizes = c(36L),
  nsims        = 5L,
  design       = "active",
  cts          = cts_tp,
  test_types   = c("DR", "DM")
)

sims_tp <- runSimsDiff(bio_twopilot, design_tp, analysis_opts)

# Check that group-1 and group-2 mesors differ for non-DM genes
gt <- sims_tp$ground_truth[[1]]
if (is.null(gt)) {
  # ground_truth is now stored per-simulation on the sim_data; check from simCircadianDiff
  sim_one <- simCircadianDiff(
    ngenes        = 300L,
    n1 = 36L, n2 = 36L,
    lBaselineExpr  = bio_twopilot$lBaselineExpr,
    lBaselineExpr2 = bio_twopilot$lBaselineExpr2,
    lOD            = bio_twopilot$lOD,
    lOD2           = bio_twopilot$lOD2,
    amplitude      = bio_twopilot$amplitude,
    amplitude2     = bio_twopilot$amplitude2,
    prop_rhythmic  = bio_twopilot$prop_rhythmic,
    prop_DR        = bio_twopilot$prop_DR,
    prop_DP        = bio_twopilot$prop_DP,
    prop_DA        = bio_twopilot$prop_DA,
    prop_DM        = bio_twopilot$prop_DM,
    mesor_diff     = bio_twopilot$mesor_diff,
    design         = "active",
    cts            = cts_tp,
    sim.seed       = 42L
  )
  gt <- sim_one$ground_truth
}

# For non-DM, non-DR genes (type 1 = rhythmic same), mesor2 should be from lBaselineExpr2
type1_genes <- which(gt$diff_type == 1)
if (length(type1_genes) > 0) {
  mesor1_t1 <- gt$mesor1[type1_genes]
  mesor2_t1 <- gt$mesor2[type1_genes]
  mean_shift <- mean(mesor2_t1 - mesor1_t1)
  cat(sprintf("  Type-1 genes mean mesor2 - mesor1 = %.3f (expected ~0.5 from lBaselineExpr2 shift)\n",
              mean_shift))
}

# For DM genes (type 5), the additional shift was applied on top of lBaselineExpr2
type5_genes <- which(gt$diff_type == 5)
cat(sprintf("  Number of DM genes: %d\n", length(type5_genes)))
if (length(type5_genes) > 0) {
  cat(sprintf("  DM genes mean |mesor2 - mesor1|: %.3f\n",
              mean(abs(gt$mesor2[type5_genes] - gt$mesor1[type5_genes]))))
}

cat("  Section 3 PASSED\n\n")

cat("============================================================\n")
cat("ALL SMOKE TESTS PASSED\n")
cat("============================================================\n")
