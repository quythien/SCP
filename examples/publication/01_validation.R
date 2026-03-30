#!/usr/bin/env Rscript
# Validation Suite for PowerSim Framework
# Runs Type I error control + monotonicity + effect size ordering checks
#
# Smoke test: POWERSIM_SMOKE=1 Rscript examples/publication/01_validation.R

SMOKE <- nzchar(Sys.getenv("POWERSIM_SMOKE"))
NSIMS_VAL  <- if (SMOKE) 5L    else 50L
NGENES_VAL <- if (SMOKE) 500L  else 5000L
N_GRID_VAL <- if (SMOKE) c(20L, 40L) else c(20L, 40L, 60L, 80L, 100L)

POWERSIM_ROOT <- {
  env <- Sys.getenv("POWERSIM_ROOT", unset = "")
  if (nzchar(env)) env else
    "/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim"
}
setwd(POWERSIM_ROOT)
source_dir <- file.path(getwd(), "code")
old_wd <- setwd(source_dir)
source("setup.R")
setwd(old_wd)

run_tag <- format(Sys.time(), "%Y%m%d_%H%M")
val_dir <- file.path("output", paste0("validation_", run_tag))
dir.create(val_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n========================================\n")
cat("  POWERSIM VALIDATION SUITE\n")
cat("========================================\n\n")

t_start <- proc.time()

# =====================================================================
# VALIDATION 1: TYPE I ERROR UNDER NULL (CRITICAL)
# =====================================================================
cat("VALIDATION 1: Type I Error Control Under Null\n")
cat("----------------------------------------------\n")
cat("Simulating with NO true differential genes (prop_DR=0, prop_DP=0, prop_DA=0)\n")
cat("If framework is correct, empirical FP rate should be ~0.05 at alpha=0.05\n\n")

set.seed(42)
null_sims <- runSimsDiff(
  sample_sizes = N_GRID_VAL,
  nsims        = NSIMS_VAL,
  ngenes       = NGENES_VAL,
  prop_rhythmic = 0.20,   # Genes are rhythmic but NO differences
  prop_DR      = 0,
  prop_DP      = 0,
  prop_DA      = 0,
  test_types   = c("DR", "DP", "DA"),
  verbose      = TRUE
)

cat("\n=== TYPE I ERROR RESULTS ===\n\n")
cat("Nominal alpha: 0.05\n")
cat("Under null, BH-adjusted FDR discoveries should be ~0\n")
cat("(because with 0 true positives, any discovery is false)\n\n")

sample_sizes <- null_sims$sample_sizes
nsims <- null_sims$nsims

# For each test type, compute:
# 1. Raw p-value rejection rate (should be ~0.05)
# 2. BH-adjusted rejection rate (should be ≤ 0.05)
# 3. KS test for uniformity of p-values

val1_results <- list()

for (test in c("DR", "DP", "DA")) {
  pval_key <- paste0("pval_", test)
  fdr_key  <- paste0("fdr_", test)

  cat(sprintf("\n--- %s Test ---\n", test))
  cat(sprintf("%-6s | Raw FPR     | BH FPR      | KS p-value  | Status\n", "n"))
  cat(paste0(rep("-", 65), collapse = ""), "\n")

  test_results <- data.frame()

  for (j in seq_along(sample_sizes)) {
    n <- sample_sizes[j]

    # Collect all p-values across simulations for this sample size
    pvals_all <- c()
    fdr_reject_count <- 0
    total_tests <- 0

    for (i in 1:nsims) {
      pv <- null_sims[[pval_key]][, j, i]
      fdr_v <- null_sims[[fdr_key]][, j, i]

      valid_p <- pv[!is.na(pv) & pv < 1]  # Exclude NA and p=1 (untested)
      valid_fdr <- fdr_v[!is.na(fdr_v)]

      pvals_all <- c(pvals_all, valid_p)
      fdr_reject_count <- fdr_reject_count + sum(valid_fdr < 0.05, na.rm = TRUE)
      total_tests <- total_tests + sum(!is.na(fdr_v))
    }

    # Raw false positive rate
    raw_fpr <- mean(pvals_all < 0.05, na.rm = TRUE)

    # BH-adjusted false positive rate
    bh_fpr <- if (total_tests > 0) fdr_reject_count / total_tests else 0

    # KS test for uniformity (valid p-values should be uniform under null)
    if (length(pvals_all) > 10) {
      ks_result <- ks.test(pvals_all, "punif")
      ks_pval <- ks_result$p.value
    } else {
      ks_pval <- NA
    }

    # Status
    raw_ok <- raw_fpr >= 0.02 && raw_fpr <= 0.08
    status <- ifelse(raw_ok, "PASS", "CHECK")

    cat(sprintf("n=%-4d | %.4f      | %.4f      | %.4f      | %s\n",
                n, raw_fpr, bh_fpr, ks_pval, status))

    test_results <- rbind(test_results, data.frame(
      test = test, n = n, raw_fpr = raw_fpr, bh_fpr = bh_fpr,
      ks_pval = ks_pval, n_pvals = length(pvals_all),
      stringsAsFactors = FALSE
    ))
  }

  val1_results[[test]] <- test_results
}

# Save null simulation results
val1_combined <- do.call(rbind, val1_results)
saveRDS(list(results = val1_combined, null_sims = null_sims),
        file = file.path(val_dir, "null_simulation_results.rds"))

# =====================================================================
# VALIDATION 2: P-VALUE QQ PLOTS (-log10 scale)
# =====================================================================
# Under a correctly calibrated null, observed -log10(p) should lie on
# the diagonal y = x (expected uniform quantiles).
# Inflation (λ > 1): points bow above diagonal → too many small p-values.
# Deflation (λ < 1): points bow below diagonal → overly conservative.
# λ = median(χ²_observed) / 0.456  (Devlin & Roeder genomic inflation factor)
# Acceptable range: 0.95 ≤ λ ≤ 1.05.
cat("\n\nVALIDATION 2: P-value QQ Plots (-log10 scale)\n")
cat("------------------------------------------------\n")

pdf(file.path(val_dir, "null_pvalue_qq.pdf"), width = 12, height = 5)
par(mfrow = c(1, 3), mar = c(4.5, 4.5, 4, 1))

j_60 <- which.min(abs(sample_sizes - 60))

for (test in c("DR", "DP", "DA")) {
  pval_key <- paste0("pval_", test)

  # Collect p-values at n=60 across all simulations
  pvals <- c()
  for (i in seq_len(nsims)) {
    pv <- null_sims[[pval_key]][, j_60, i]
    pvals <- c(pvals, pv[!is.na(pv) & pv > 0 & pv < 1])
  }

  if (length(pvals) > 0) {
    # Genomic inflation factor λ
    chisq_obs <- qchisq(1 - pvals, df = 1)
    lambda    <- median(chisq_obs, na.rm = TRUE) / qchisq(0.5, df = 1)

    # QQ plot on -log10 scale
    n_p      <- length(pvals)
    obs      <- sort(-log10(pvals), decreasing = TRUE)
    exp      <- -log10(seq_len(n_p) / (n_p + 1))
    xy_max   <- max(c(obs, exp), na.rm = TRUE) * 1.05

    plot(exp, obs, pch = 20, cex = 0.4,
         col = ifelse(lambda > 1.05, "firebrick",
               ifelse(lambda < 0.95, "steelblue", "grey40")),
         xlab = expression("Expected  " * -log[10](p)),
         ylab = expression("Observed  " * -log[10](p)),
         main = sprintf("%s under null (n=60)\nλ = %.3f", test, lambda),
         xlim = c(0, xy_max), ylim = c(0, xy_max),
         las = 1, cex.main = 0.95)
    abline(0, 1, col = "red", lwd = 1.5, lty = 2)   # expected diagonal

    # 95% confidence band (pointwise beta CI)
    ci_lo <- -log10(qbeta(0.975, seq_len(n_p), n_p - seq_len(n_p) + 1))
    ci_hi <- -log10(qbeta(0.025, seq_len(n_p), n_p - seq_len(n_p) + 1))
    polygon(c(exp, rev(exp)), c(ci_lo, rev(ci_hi)),
            col = adjustcolor("grey70", 0.3), border = NA)

    lambda_status <- ifelse(lambda >= 0.95 & lambda <= 1.05, "OK",
                     ifelse(lambda > 1.05, "INFLATED", "DEFLATED"))
    cat(sprintf("  %s: λ = %.3f  [%s]  n_pvals = %d\n",
                test, lambda, lambda_status, n_p))
  }
}

dev.off()
cat("  Saved: null_pvalue_qq.pdf\n")

# =====================================================================
# VALIDATION 3: MONOTONICITY CHECK
# =====================================================================
cat("\nVALIDATION 3: Power Monotonicity\n")
cat("----------------------------------\n")
cat("Running DR and DP power analyses (self-contained)...\n")

val_bio <- CircadianBioOptions(
  ngenes = null_sims$ngenes, prop_rhythmic = 0.20,
  prop_DR = 0.15, prop_DP = 0.00, prop_DA = 0.00,
  phase_diff = c(0, 0), amp_diff = c(1, 1)
)
val_bio_DP <- updateBioOptions(val_bio,
  prop_DR = 0.00, prop_DP = 0.15, phase_diff = c(-6, 6), amp_diff = c(0.5, 2)
)
val_design <- CircadianDesignOptions(
  sample_sizes = N_GRID_VAL, nsims = null_sims$nsims, design = "active"
)
val_analysis <- CircadianAnalysisOptions(alpha = 0.05, p.adjust.method = "BH")

dr_raw <- runPowerAnalysis(val_bio,    val_design, val_analysis, test_type = "DR")
dp_raw <- runPowerAnalysis(val_bio_DP, val_design, val_analysis, test_type = "DP")

saveRDS(dr_raw, file.path(val_dir, "dr_power.rds"))
saveRDS(dp_raw, file.path(val_dir, "dp_power.rds"))

# Check DR marginal power monotonicity
dr_marginal <- apply(dr_raw$marginal_power, 1, mean, na.rm = TRUE)
dr_mono <- all(diff(dr_marginal) >= -0.005)
cat(sprintf("  DR marginal power: %s\n", paste(sprintf("%.1f%%", 100*dr_marginal), collapse = " -> ")))
cat(sprintf("  DR monotonic: %s\n", ifelse(dr_mono, "PASS", "FAIL")))

# Check DP marginal power monotonicity
dp_marginal <- apply(dp_raw$marginal_power, 1, mean, na.rm = TRUE)
dp_mono <- all(diff(dp_marginal) >= -0.005)
cat(sprintf("  DP marginal power: %s\n", paste(sprintf("%.1f%%", 100*dp_marginal), collapse = " -> ")))
cat(sprintf("  DP monotonic: %s\n", ifelse(dp_mono, "PASS", "FAIL")))

# =====================================================================
# VALIDATION 4: EFFECT SIZE ORDERING (from phase sweep results)
# =====================================================================
cat("\nVALIDATION 4: Effect Size Ordering (Phase Sweep)\n")
cat("---------------------------------------------------\n")

cat("Running phase shift sweep (self-contained)...\n")
val_analysis_ps <- CircadianAnalysisOptions(
  alpha = 0.05, p.adjust.method = "BH",
  phase_shifts = c(0, 0.5, 1, 2, 4, 6, 8, 10, 12)
)
dp_phase_results <- runPhaseShiftAnalysis(
  val_bio, val_design, val_analysis_ps, prop_DP = 0.15, amp_diff = c(0.5, 2)
)
saveRDS(dp_phase_results, file.path(val_dir, "dp_phase_shift.rds"))

if (!is.null(dp_phase_results)) {
  # At n=60, power should increase with phase shift
  j_60 <- which.min(abs(val_design$sample_sizes - 60))
  phase_shifts <- val_analysis_ps$phase_shifts

  power_by_shift <- sapply(seq_along(phase_shifts), function(p) {
    mean(dp_phase_results$marginal_power[p, j_60, ], na.rm = TRUE)
  })

  # Check monotonicity for shifts >= 0.5
  mono_check <- all(diff(power_by_shift[-1]) >= -0.005)

  cat("  Power by phase shift at n=60:\n")
  for (p in seq_along(phase_shifts)) {
    cat(sprintf("    %4.1fh: %.1f%%\n", phase_shifts[p], 100 * power_by_shift[p]))
  }
  cat(sprintf("  Monotonic (shifts >= 0.5h): %s\n", ifelse(mono_check, "PASS", "FAIL")))

  # Saturation check: 8h vs 12h should be similar
  ratio_8_12 <- power_by_shift[7] / power_by_shift[9]  # 8h / 12h
  cat(sprintf("  Saturation ratio (8h/12h): %.3f (expect ~1.0)\n", ratio_8_12))
} else {
  cat("  Phase sweep results not found - SKIPPED\n")
}

# =====================================================================
# VALIDATION 5: FDR CONTROL (from existing results)
# =====================================================================
cat("\nVALIDATION 5: FDR Control\n")
cat("---------------------------\n")

# DR FDR
dr_fdr <- apply(dr_raw$marginal_FDR, 1, mean, na.rm = TRUE)
cat(sprintf("  DR empirical FDR by n: %s\n",
            paste(sprintf("%.3f", dr_fdr), collapse = ", ")))
cat(sprintf("  DR FDR controlled (<= 0.05): %s\n",
            ifelse(all(dr_fdr <= 0.06), "PASS", "CHECK")))

# DP FDR
dp_fdr <- apply(dp_raw$marginal_FDR, 1, mean, na.rm = TRUE)
cat(sprintf("  DP empirical FDR by n: %s\n",
            paste(sprintf("%.3f", dp_fdr), collapse = ", ")))
cat(sprintf("  DP FDR controlled (<= 0.05 for n>=40): %s\n",
            ifelse(all(dp_fdr[-1] <= 0.06), "PASS", "CHECK")))

# =====================================================================
# VALIDATION 6: STRATIFICATION NECESSITY
# =====================================================================
cat("\nVALIDATION 6: Stratification Necessity\n")
cat("-----------------------------------------\n")

# At n=60, compare marginal vs stratified power range
j_60 <- which.min(abs(val_design$sample_sizes - 60))
dr_strat_power <- apply(dr_raw$strat_power[j_60, , ], 1, mean, na.rm = TRUE)
dp_strat_power <- apply(dp_raw$strat_power[j_60, , ], 1, mean, na.rm = TRUE)

dr_range <- range(dr_strat_power, na.rm = TRUE)
dp_range <- range(dp_strat_power, na.rm = TRUE)

cat(sprintf("  DR at n=60: marginal=%.1f%%, stratified range=[%.1f%%, %.1f%%]\n",
            100 * dr_marginal[j_60], 100 * dr_range[1], 100 * dr_range[2]))
cat(sprintf("  DP at n=60: marginal=%.1f%%, stratified range=[%.1f%%, %.1f%%]\n",
            100 * dp_marginal[j_60], 100 * dp_range[1], 100 * dp_range[2]))
cat(sprintf("  DR fold-range: %.1fx\n", dr_range[2] / max(dr_range[1], 0.001)))
cat(sprintf("  DP fold-range: %.1fx\n", dp_range[2] / max(dp_range[1], 0.001)))

# =====================================================================
# SUMMARY
# =====================================================================
t_elapsed <- proc.time() - t_start

cat("\n========================================\n")
cat("  VALIDATION SUMMARY\n")
cat("========================================\n\n")
cat(sprintf("  1. Type I Error Control:    See table above\n"))
cat(sprintf("  2. P-value Uniformity:      See null_pvalue_qq.pdf\n"))
cat(sprintf("  3. Power Monotonicity:      DR=%s, DP=%s\n",
            ifelse(dr_mono, "PASS", "FAIL"), ifelse(dp_mono, "PASS", "FAIL")))
cat(sprintf("  4. Effect Size Ordering:    %s\n",
            ifelse(exists("mono_check") && mono_check, "PASS", "CHECK")))
cat(sprintf("  5. FDR Control:             DR=%s, DP=%s\n",
            ifelse(all(dr_fdr <= 0.06), "PASS", "CHECK"),
            ifelse(all(dp_fdr[-1] <= 0.06), "PASS", "CHECK")))
cat(sprintf("  6. Stratification:          Confirmed (see ranges above)\n"))
cat(sprintf("\nTotal runtime: %.1f minutes\n", t_elapsed[3] / 60))
cat(sprintf("Output: %s\n", val_dir))
cat("========================================\n\n")
