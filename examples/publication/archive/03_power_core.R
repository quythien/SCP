#' =======================================================================
#' End-to-End Power Analysis Pipeline for Circadian Differential Expression
#' =======================================================================
#'
#' WHAT THIS DOES:
#'   Simulation-based power analysis for detecting circadian differential
#'   expression between two groups (e.g., young vs old), using a
#'   semi-parametric framework with  r = A/sigma stratification
#'   and the DCP pipeline for hypothesis testing.
#'
#' SCIENTIFIC QUESTIONS ANSWERED:
#'   1. Differential Rhythmicity (DR):
#'      "How many subjects per group do we need to detect genes that
#'       gain or lose rhythmicity between conditions?"
#'
#'   2. Differential Phase (DP):
#'      "How many subjects per group do we need to detect genes whose
#'       circadian phase shifts between conditions (e.g., ±6h shift)?"
#'
#'   3. Phase Shift Sensitivity:
#'      "How large must a phase shift be before we can reliably detect it,
#'       and how does this interact with sample size and signal strength?"
#'
#' FRAMEWORK:
#'   - Parametric cosinor model: y = M + A*cos(omega*t - phi) + epsilon
#'   - Empirical parameter distributions from pilot data (semi-parametric)
#'   -  stratification by r = A/sigma (signal-to-noise ratio)
#'   - DCP pipeline for differential testing (LR tests for phase/amplitude)
#'   - Passive design: sampling times drawn from real TOD distribution
#'
#' OUTPUTS (3 PDF figures + console summaries):
#'   1. output/example_pipeline/figures/dr_power.pdf       (6 panels)
#'   2. output/example_pipeline/figures/dp_power.pdf       (6 panels)
#'   3. output/example_pipeline/figures/phase_shift.pdf    (6 panels)
#'
#'
#' USAGE:
#'   Rscript examples/run_pipeline.R
#'
#' @author Thien Quy Pham


# =====================================================================
# SECTION 1: SETUP & CONFIGURATION
# =====================================================================

# Smoke test: POWERSIM_SMOKE=1 Rscript examples/publication/03_power_core.R
SMOKE        <- nzchar(Sys.getenv("POWERSIM_SMOKE"))
NSIMS_CORE   <- if (SMOKE) 5L   else 50L
NGENES_CORE  <- if (SMOKE) 500L else 5000L
N_GRID_CORE  <- if (SMOKE) c(20L, 40L) else c(20L, 40L, 60L, 80L, 100L)

# Set POWERSIM_ROOT as env var for portability:
#   Linux/Mac shell: export POWERSIM_ROOT=/path/to/PowerSim
POWERSIM_ROOT <- {
  env <- Sys.getenv("POWERSIM_ROOT", unset = "")
  if (nzchar(env)) env else
    "/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim"
}

# DATA_HUMAN path (outside PowerSim root; override via env var on server):
#   export DATA_HUMAN=/path/to/combined_data.rds
DATA_HUMAN <- {
  env <- Sys.getenv("DATA_HUMAN", unset = "")
  if (nzchar(env)) env else
    "/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Collaborative/Paper/Congruence/PNAS_aging/data/combined_data.rds"
}

setwd(POWERSIM_ROOT)
source_dir <- file.path(getwd(), "code")
old_wd <- setwd(source_dir)
source("setup.R")
setwd(old_wd)

# Source SE plotting functions
source("code/plot_with_se.R")


# =====================================================================
# SECTION 2: LOAD PILOT DATA
# =====================================================================

cat("Loading pilot expression data (PFC younger: BA11 + BA47)...\n")

COMBINED <- readRDS(DATA_HUMAN)

expr_sample_names <- colnames(COMBINED$expr)
pheno_order <- match(expr_sample_names, COMBINED$pheno$sample_name)
valid_samples <- !is.na(pheno_order)
COMBINED$expr <- COMBINED$expr[, valid_samples]
pheno_order <- pheno_order[valid_samples]
pheno_data <- COMBINED$pheno[pheno_order, ]
pheno_data$tod <- if ("TOD.x" %in% colnames(pheno_data)) pheno_data$TOD.x else pheno_data$TOD.y
pheno_data$age_group_final <- if ("AgeGroup" %in% colnames(pheno_data)) pheno_data$AgeGroup else pheno_data$age_group
complete_samples <- !is.na(pheno_data$age_group_final) & !is.na(pheno_data$tod) &
  pheno_data$age_group_final %in% c("younger", "older")
pheno_clean <- pheno_data[complete_samples, ]
younger_idx <- pheno_clean$age_group_final == "younger"
expr_raw_young <- COMBINED$expr[, complete_samples][, younger_idx]
times_raw      <- pheno_clean$tod[younger_idx]
rm(COMBINED, expr_sample_names, pheno_order, valid_samples, pheno_data,
   complete_samples, pheno_clean, younger_idx)

# Use prepCircadianData for standard validation (log2 already; coerces to matrix)
prep_young   <- prepCircadianData(expr_raw_young, times = times_raw, input_type = "log2")
expr_younger <- prep_young$data
times_young  <- prep_young$times
rm(expr_raw_young, times_raw)

cat(sprintf("  Expression matrix: %d genes x %d samples\n", nrow(expr_younger), ncol(expr_younger)))
cat(sprintf("  Reference TOD: n=%d younger subjects\n\n", length(times_young)))


# =====================================================================
# SECTION 3: ESTIMATE PARAMETERS FROM PILOT DATA
# =====================================================================
# Fit cosinor to each gene in the pilot dataset to learn empirical
# distributions of baseline expression, noise, amplitude, and phase.
# This is the data-driven alternative to specifying arbitrary parametric
# distributions. Users without pilot data can omit this step and rely
# on the built-in "ba11_younger" default.

cat("Estimating circadian parameters from pilot data...\n\n")

opts_bio <- estCircadianParam(
  data          = expr_younger,
  times         = times_young,
  period        = 24,
  prop_DR       = 0.15,
  prop_DP       = 0.10,
  phase_diff    = c(-6, 6),
  amp_diff      = c(0.5, 2)
)

rm(expr_younger)  # free memory

# Override ngenes for simulation (pilot has ~20k, we simulate 5k)
opts_bio <- updateBioOptions(opts_bio, ngenes = NGENES_CORE)


# =====================================================================
# SECTION 4: DESIGN & ANALYSIS CONFIGURATION
# =====================================================================

opts_design <- CircadianDesignOptions(
  sample_sizes = N_GRID_CORE,
  nsims        = NSIMS_CORE,
  design       = "passive",
  cts          = times_young,
  test_types   = c("DR", "DP")
)

opts_analysis <- CircadianAnalysisOptions(
  alpha           = 0.05,
  p.adjust.method = "BH",
  parallel.ncores = 1,
  amp.cutoff      = 0,
  target_effect   = 0.1,
  fdr_thresholds  = c(0.01, 0.05, 0.10, 0.20),
  reference_n     = 60,
  phase_shifts    = c(0, 0.5, 1, 2, 4, 6, 8, 10, 12) # this is for phase shift analysis 
)

# Output directories — timestamped to avoid overwriting previous runs
run_tag <- format(Sys.time(), "%Y%m%d_%H%M")
base_out <- file.path("output", paste0("run_", run_tag))
out_dir  <- file.path(base_out, "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("Run output root: %s/\n", base_out))

t_start <- proc.time()

cat("====================================================================\n")
cat("END-TO-END POWER ANALYSIS PIPELINE\n")
cat("====================================================================\n\n")
cat("Configuration:\n")
print(opts_bio)
cat("\n")
print(opts_design)
cat("\n")
print(opts_analysis)
cat(sprintf("\nOutput: %s/\n\n", out_dir))


# =====================================================================
# SECTION 5: DR ANALYSIS
# =====================================================================
# Question: "How many subjects to detect genes gaining/losing rhythmicity?"
# Setup: 15% of genes have differential rhythmicity (types 2,3).
#        No differential phase or amplitude (isolated DR test).

cat("\n====================================================================\n")
cat("ANALYSIS 1: DIFFERENTIAL RHYTHMICITY (DR)\n")
cat("====================================================================\n\n")

opts_bio_DR <- updateBioOptions(opts_bio,
  prop_DR    = 0.15,
  prop_DP    = 0.00,
  phase_diff = c(0, 0),
  amp_diff   = c(1, 1)
)

dr_power_raw <- runPowerAnalysis(opts_bio_DR, opts_design, opts_analysis,
                                 test_type = "DR")

dr_results_file <- file.path(base_out, "dr_power_raw_pvalues.rds")
save(dr_power_raw, file = dr_results_file)
cat(sprintf("\nDR results saved: %s\n", dr_results_file))

# Generate 6-panel PDF with SE bars
dr_fig <- file.path(out_dir, "dr_power.pdf")
plotWithSE(dr_results_file, dr_fig, test_name = "DR", analysis.opts = opts_analysis)


# =====================================================================
# SECTION 6: DP ANALYSIS
# =====================================================================
# Question: "How many subjects to detect phase-shifted genes (±6h shift)?"
# Setup: 15% of genes have differential phase (type 4).
#        No differential rhythmicity or amplitude (isolated DP test).

cat("====================================================================\n")
cat("ANALYSIS 2: DIFFERENTIAL PHASE (DP)\n")
cat("====================================================================\n\n")

opts_bio_DP <- updateBioOptions(opts_bio,
  prop_DR    = 0.00,
  prop_DP    = 0.15,
  phase_diff = c(-6, 6),
  amp_diff   = c(0.5, 2)
)

dp_power_raw <- runPowerAnalysis(opts_bio_DP, opts_design, opts_analysis,
                                 test_type = "DP")

dp_results_file <- file.path(base_out, "dp_power_raw_pvalues.rds")
save(dp_power_raw, file = dp_results_file)
cat(sprintf("\nDP results saved: %s\n", dp_results_file))

# Generate 6-panel PDF with SE bars
dp_fig <- file.path(out_dir, "dp_power.pdf")
plotWithSE(dp_results_file, dp_fig, test_name = "DP", analysis.opts = opts_analysis)


# =====================================================================
# SECTION 7: PHASE SHIFT SENSITIVITY ANALYSIS
# =====================================================================
# Question: "How large must a phase shift be to detect it?"
# This analysis sweeps over phase shift magnitudes AND sample sizes,
# producing a 4D array: [phase_shift, sample_size, r_stratum, sim].

# cat("ANALYSIS 3: PHASE SHIFT SENSITIVITY\n")
# cat("====================================================================\n\n")
#
dp_phase_results <- runPhaseShiftAnalysis(opts_bio, opts_design, opts_analysis,
                                           prop_DP = 0.15, amp_diff = c(0.5, 2))
ps_results_file <- file.path(base_out, "dp_power_phase_shift_results.rds")
save(dp_phase_results, file = ps_results_file, ascii = TRUE)
ps_fig <- file.path(out_dir, "phase_shift.pdf")
plotPhaseShiftWithSE(ps_results_file, ps_fig, analysis.opts = opts_analysis)

sample_sizes <- opts_design$sample_sizes
phase_shifts <- opts_analysis$phase_shifts
nsims <- opts_design$nsims
idx_n60 <- which.min(abs(sample_sizes - opts_analysis$reference_n))
dp_phase_results <- NULL  # placeholder


# =====================================================================
# SECTION 8: SUMMARY & KEY FINDINGS
# =====================================================================

cat("====================================================================\n")
cat("SUMMARY TABLES\n")
cat("====================================================================\n\n")

# --- DR Summary (PROPER-style: Power, Type I Error, FDR, TD, FD, FDC) ---
cat("--- DR Power Summary (FDR 5%, marginal across r strata) ---\n")
dr_summary <- summaryRunPower(dr_power_raw)
cat("\n")

# --- DP Summary ---
cat("--- DP Power Summary (FDR 5%, ±6h shift, marginal across r strata) ---\n")
dp_summary <- summaryRunPower(dp_power_raw)
cat("\n")

## Phase Shift tables SKIPPED (phase sweep disabled for fast testing)
if (!is.null(dp_phase_results)) {
  cat("--- Phase Shift: Marginal Power (FDR 5%) ---\n")
  cat(sprintf("%-11s |", "Phase (h)"))
  for (n in sample_sizes) cat(sprintf(" n=%-4d |", n))
  cat("\n")
  cat(paste0(rep("-", 14 + length(sample_sizes) * 9), collapse = ""), "\n")
  for (p in seq_along(phase_shifts)) {
    cat(sprintf("%-11s |", sprintf("%.1f", phase_shifts[p])))
    for (j in seq_along(sample_sizes)) {
      mp <- mean(dp_phase_results$marginal_power[p, j, ], na.rm = TRUE)
      cat(sprintf(" %5.1f%% |", 100 * mp))
    }
    cat("\n")
  }
  cat("\n--- Phase Shift: Empirical FDR ---\n")
  cat(sprintf("%-11s |", "Phase (h)"))
  for (n in sample_sizes) cat(sprintf(" n=%-4d |", n))
  cat("\n")
  cat(paste0(rep("-", 14 + length(sample_sizes) * 9), collapse = ""), "\n")
  for (p in seq_along(phase_shifts)) {
    cat(sprintf("%-11s |", sprintf("%.1f", phase_shifts[p])))
    for (j in seq_along(sample_sizes)) {
      mfdr <- mean(dp_phase_results$marginal_FDR[p, j, ], na.rm = TRUE)
      cat(sprintf(" %5.3f  |", mfdr))
    }
    cat("\n")
  }
  cat("\n--- Phase Shift: Average True Discoveries ---\n")
  cat(sprintf("%-11s |", "Phase (h)"))
  for (n in sample_sizes) cat(sprintf(" n=%-4d |", n))
  cat("\n")
  cat(paste0(rep("-", 14 + length(sample_sizes) * 9), collapse = ""), "\n")
  for (p in seq_along(phase_shifts)) {
    cat(sprintf("%-11s |", sprintf("%.1f", phase_shifts[p])))
    for (j in seq_along(sample_sizes)) {
      mtd <- mean(dp_phase_results$marginal_TD[p, j, ], na.rm = TRUE)
      cat(sprintf(" %6.1f |", mtd))
    }
    cat("\n")
  }
} else {
  cat("--- Phase Shift: SKIPPED (fast test mode) ---\n")
}

# --- Key Findings ---
cat("\n====================================================================\n")
cat("KEY FINDINGS\n")
cat("====================================================================\n\n")

# DR: n for 80% power
dr_avg_power <- rowMeans(dr_power_raw$marginal_power, na.rm = TRUE)
for (j in seq_along(sample_sizes)) {
  if (dr_avg_power[j] >= 0.80) {
    cat(sprintf("DR:  n=%d per group needed for >=80%% marginal power at FDR 5%%\n", sample_sizes[j]))
    break
  }
  if (j == length(sample_sizes)) {
    cat(sprintf("DR:  80%% power not reached at n=%d (max tested); best=%.1f%%\n",
                sample_sizes[j], 100 * dr_avg_power[j]))
  }
}

# DR: FDR control check
dr_avg_fdr <- rowMeans(dr_power_raw$marginal_FDR, na.rm = TRUE)
cat(sprintf("DR:  Empirical FDR at n=%d: %.3f (nominal: 0.05)\n",
            sample_sizes[length(sample_sizes)], dr_avg_fdr[length(sample_sizes)]))

# DR: Type I error check
dr_avg_alpha <- rowMeans(dr_power_raw$marginal_alpha, na.rm = TRUE)
cat(sprintf("DR:  Empirical Type I error at n=%d: %.4f\n",
            sample_sizes[length(sample_sizes)], dr_avg_alpha[length(sample_sizes)]))

# DP: n for 80% power (6h shift)
if (!is.null(dp_phase_results)) {
  idx_6h <- which(phase_shifts == 6)
  dp_6h_power <- apply(dp_phase_results$marginal_power[idx_6h, , , drop = FALSE], 2, mean, na.rm = TRUE)
  for (j in seq_along(sample_sizes)) {
    if (!is.na(dp_6h_power[j]) && dp_6h_power[j] >= 0.80) {
      cat(sprintf("DP:  n=%d per group needed for >=80%% marginal power (6h shift) at FDR 5%%\n", sample_sizes[j]))
      break
    }
    if (j == length(sample_sizes)) {
      cat(sprintf("DP:  80%% power not reached at n=%d (6h shift, max tested); best=%.1f%%\n",
                  sample_sizes[j], 100 * dp_6h_power[j]))
    }
  }
} else {
  # Use dp_power_raw marginal power when phase sweep is skipped
  dp_avg_power <- rowMeans(dp_power_raw$marginal_power, na.rm = TRUE)
  for (j in seq_along(sample_sizes)) {
    if (dp_avg_power[j] >= 0.80) {
      cat(sprintf("DP:  n=%d per group needed for >=80%% marginal power at FDR 5%%\n", sample_sizes[j]))
      break
    }
    if (j == length(sample_sizes)) {
      cat(sprintf("DP:  80%% power not reached at n=%d (max tested); best=%.1f%%\n",
                  sample_sizes[j], 100 * dp_avg_power[j]))
    }
  }
}

# DP: FDR control check
dp_avg_fdr <- rowMeans(dp_power_raw$marginal_FDR, na.rm = TRUE)
cat(sprintf("DP:  Empirical FDR at n=%d: %.3f (nominal: 0.05)\n",
            sample_sizes[length(sample_sizes)], dp_avg_fdr[length(sample_sizes)]))

# Phase shift: minimum detectable at reference_n
if (!is.null(dp_phase_results)) {
  n_phase <- length(phase_shifts)
  ps_ref_power <- apply(dp_phase_results$marginal_power[, idx_n60, , drop = FALSE], 1, mean, na.rm = TRUE)
  for (p in seq_along(phase_shifts)) {
    if (phase_shifts[p] > 0 && !is.na(ps_ref_power[p]) && ps_ref_power[p] >= 0.80) {
      cat(sprintf("Phase: %.1fh+ shift needed for >=80%% power at n=%d\n", phase_shifts[p], sample_sizes[idx_n60]))
      break
    }
    if (p == n_phase) {
      cat(sprintf("Phase: 80%% power not reached at n=%d for any tested shift\n", sample_sizes[idx_n60]))
    }
  }
} else {
  cat("Phase: SKIPPED (fast test mode)\n")
}

# FDC summary
cat("\n--- False Discovery Cost (FD/TD) ---\n")
cat(sprintf("DR at n=%d:  FDC = %.2f (%.1f FD per TD)\n",
            sample_sizes[length(sample_sizes)],
            dr_summary$FDC[nrow(dr_summary)],
            dr_summary$FDC[nrow(dr_summary)]))
cat(sprintf("DP at n=%d:  FDC = %.2f (%.1f FD per TD)\n",
            sample_sizes[length(sample_sizes)],
            dp_summary$FDC[nrow(dp_summary)],
            dp_summary$FDC[nrow(dp_summary)]))




# =====================================================================
# SECTION 10: DCP vs CircaCompare METHOD COMPARISON
# =====================================================================
# Compare the two differential testing methods on identical simulated
# data. DCP uses LR tests (asymptotically optimal); CircaCompare uses
# Wald t-tests from per-gene NLS.
#
# Comparable tests:
#   DA:  DCP delta-A LRT  vs  CircaCompare alpha1 Wald
#   DP:  DCP delta-phi LRT vs  CircaCompare phi1 Wald
#   DM:  DCP delta-M LRT  vs  CircaCompare k1 Wald
#   DR:  DCP only (CircaCompare has no DR test)
#
# NOTE: CircaCompare is ~1-3 sec/gene (NLS), so we use reduced
#       nsims and fewer sample sizes for speed.
#
# TO RUN: uncomment this section (remove leading "# " from each line)
# WARNING: very slow — use reduced N_GRID_CORE and NSIMS_CORE

# cat("====================================================================\n")
# cat("ANALYSIS 4: METHOD COMPARISON (DCP vs CircaCompare)\n")
# cat("====================================================================\n\n")
#
# compare_base <- file.path(base_out, "method_comparison")
# compare_dir  <- file.path(compare_base, "figures")
# dir.create(compare_dir, recursive = TRUE, showWarnings = FALSE)
#
# # Reduced design for speed (CircaCompare is slow: ~1-3 sec/gene NLS)
# opts_design_compare <- CircadianDesignOptions(
#   sample_sizes = N_GRID_CORE,
#   nsims        = NSIMS_CORE,
#   design       = "passive",
#   cts          = times_young,
#   test_types   = c("DP", "DM")
# )
#
# # Reduced ngenes for CircaCompare (5000 genes x NLS is infeasible)
# opts_bio_compare <- opts_bio
# opts_bio_compare$ngenes <- 5000
# cat(sprintf("  CircaCompare comparison uses ngenes=%d\n\n", opts_bio_compare$ngenes))
#
# opts_analysis_DCP <- CircadianAnalysisOptions(
#   alpha = 0.05, p.adjust.method = "BH", DCmethod = "DCP"
# )
# opts_analysis_CC <- CircadianAnalysisOptions(
#   alpha = 0.05, p.adjust.method = "BH", DCmethod = "CircaCompare"
# )
#
# # --- 10A: DP comparison (Differential Phase) ---
# cat("\n--- 10A: Differential Phase (DP) ---\n")
# cat("  Setting up: 15% DP genes, phase_diff = [-6, 6]\n\n")
#
# opts_bio_DP_cmp <- updateBioOptions(opts_bio_compare,
#   prop_DR    = 0.00,
#   prop_DP    = 0.15,
#   prop_DA    = 0.00,
#   phase_diff = c(-6, 6),
#   amp_diff   = c(0.5, 2)
# )
#
# cat("  Running DCP...\n")
# dp_dcp <- runSimsDiff(opts_bio_DP_cmp, opts_design_compare, opts_analysis_DCP)
# saveRDS(dp_dcp, file = file.path(compare_base, "dp_dcp.rds"))
# cat(sprintf("  Saved: %s\n", file.path(compare_base, "dp_dcp.rds")))
#
# cat("  Running CircaCompare...\n")
# dp_cc <- runSimsDiff(opts_bio_DP_cmp, opts_design_compare, opts_analysis_CC)
# saveRDS(dp_cc, file = file.path(compare_base, "dp_cc.rds"))
# cat(sprintf("  Saved: %s\n", file.path(compare_base, "dp_cc.rds")))
#
#
# # --- 10C: Compute comparison data frames ---
# cat("\n--- Computing comparison statistics ---\n")
#
# .compute_comparison <- function(dcp_out, cc_out, test_type) {
#   fdr_key <- paste0("fdr_", test_type)
#   sample_sizes_cmp <- dcp_out$sample_sizes
#   nsims_cmp <- dcp_out$nsims
#   ngenes_cmp <- dcp_out$ngenes
#   rows <- list()
#   for (j in seq_along(sample_sizes_cmp)) {
#     dcp_power_vals <- numeric(nsims_cmp)
#     cc_power_vals  <- numeric(nsims_cmp)
#     dcp_fdr_vals   <- numeric(nsims_cmp)
#     cc_fdr_vals    <- numeric(nsims_cmp)
#     dcp_td_vals    <- numeric(nsims_cmp)
#     cc_td_vals     <- numeric(nsims_cmp)
#     dcp_fd_vals    <- numeric(nsims_cmp)
#     cc_fd_vals     <- numeric(nsims_cmp)
#     for (i in 1:nsims_cmp) {
#       diff_type <- dcp_out$diff_type[[i]]
#       if (test_type == "DA") {
#         is_target <- diff_type == 5
#       } else if (test_type == "DP") {
#         is_target <- diff_type == 4
#       } else {
#         is_target <- rep(FALSE, ngenes_cmp)
#       }
#       is_null <- !is_target
#       n_target <- sum(is_target)
#       dcp_disc <- dcp_out[[fdr_key]][, j, i] <= 0.05
#       cc_disc  <- cc_out[[fdr_key]][, j, i] <= 0.05
#       dcp_td <- sum(dcp_disc & is_target, na.rm = TRUE)
#       cc_td  <- sum(cc_disc & is_target, na.rm = TRUE)
#       dcp_fd <- sum(dcp_disc & is_null, na.rm = TRUE)
#       cc_fd  <- sum(cc_disc & is_null, na.rm = TRUE)
#       dcp_td_vals[i] <- dcp_td
#       cc_td_vals[i]  <- cc_td
#       dcp_fd_vals[i] <- dcp_fd
#       cc_fd_vals[i]  <- cc_fd
#       dcp_power_vals[i] <- if (n_target > 0) dcp_td / n_target else NA
#       cc_power_vals[i]  <- if (n_target > 0) cc_td / n_target else NA
#       dcp_total <- dcp_td + dcp_fd
#       cc_total  <- cc_td + cc_fd
#       dcp_fdr_vals[i] <- if (dcp_total > 0) dcp_fd / dcp_total else NA
#       cc_fdr_vals[i]  <- if (cc_total > 0) cc_fd / cc_total else NA
#     }
#     rows[[j]] <- data.frame(
#       test_type    = test_type,
#       n            = sample_sizes_cmp[j],
#       DCP_Power    = mean(dcp_power_vals, na.rm = TRUE),
#       CC_Power     = mean(cc_power_vals, na.rm = TRUE),
#       DCP_Power_SE = sd(dcp_power_vals, na.rm = TRUE) / sqrt(sum(!is.na(dcp_power_vals))),
#       CC_Power_SE  = sd(cc_power_vals, na.rm = TRUE) / sqrt(sum(!is.na(cc_power_vals))),
#       DCP_FDR      = mean(dcp_fdr_vals, na.rm = TRUE),
#       CC_FDR       = mean(cc_fdr_vals, na.rm = TRUE),
#       DCP_FDR_SE   = sd(dcp_fdr_vals, na.rm = TRUE) / sqrt(sum(!is.na(dcp_fdr_vals))),
#       CC_FDR_SE    = sd(cc_fdr_vals, na.rm = TRUE) / sqrt(sum(!is.na(cc_fdr_vals))),
#       DCP_Avg_TD   = mean(dcp_td_vals),
#       CC_Avg_TD    = mean(cc_td_vals),
#       DCP_Avg_FD   = mean(dcp_fd_vals),
#       CC_Avg_FD    = mean(cc_fd_vals),
#       stringsAsFactors = FALSE
#     )
#   }
#   do.call(rbind, rows)
# }
#
# .compute_dm_typeI <- function(dcp_out, cc_out) {
#   sample_sizes_cmp <- dcp_out$sample_sizes
#   nsims_cmp <- dcp_out$nsims
#   ngenes_cmp <- dcp_out$ngenes
#   rows <- list()
#   for (j in seq_along(sample_sizes_cmp)) {
#     dcp_vals <- numeric(nsims_cmp)
#     cc_vals  <- numeric(nsims_cmp)
#     for (i in 1:nsims_cmp) {
#       dcp_vals[i] <- sum(dcp_out$fdr_DM[, j, i] <= 0.05, na.rm = TRUE) / ngenes_cmp
#       cc_vals[i]  <- sum(cc_out$fdr_DM[, j, i] <= 0.05, na.rm = TRUE) / ngenes_cmp
#     }
#     rows[[j]] <- data.frame(
#       test_type    = "DM",
#       n            = sample_sizes_cmp[j],
#       DCP_TypeI    = mean(dcp_vals),
#       CC_TypeI     = mean(cc_vals),
#       DCP_TypeI_SE = sd(dcp_vals) / sqrt(nsims_cmp),
#       CC_TypeI_SE  = sd(cc_vals) / sqrt(nsims_cmp),
#       stringsAsFactors = FALSE
#     )
#   }
#   do.call(rbind, rows)
# }
#
# dp_table <- .compute_comparison(dp_dcp, dp_cc, "DP")
# dm_table <- .compute_dm_typeI(dp_dcp, dp_cc)
#
# # --- 10B: Print tables to console ---
# cat("\n====================================================================\n")
# cat("METHOD COMPARISON RESULTS\n")
# cat("====================================================================\n\n")
# for (tbl_name in c("DP")) {
#   tbl <- dp_table
#   cat(sprintf("--- %s: Power (FDR 5%%) ---\n", tbl_name))
#   cat(sprintf("%-6s | %-16s | %-16s\n", "n", "DCP", "CircaCompare"))
#   cat(paste0(rep("-", 44), collapse = ""), "\n")
#   for (r in seq_len(nrow(tbl))) {
#     cat(sprintf("%-6d | %5.1f%% (+/-%.1f%%) | %5.1f%% (+/-%.1f%%)\n",
#                 tbl$n[r],
#                 100 * tbl$DCP_Power[r], 100 * tbl$DCP_Power_SE[r],
#                 100 * tbl$CC_Power[r], 100 * tbl$CC_Power_SE[r]))
#   }
#   cat("\n")
#   cat(sprintf("--- %s: Empirical FDR (nominal 5%%) ---\n", tbl_name))
#   cat(sprintf("%-6s | %-16s | %-16s\n", "n", "DCP", "CircaCompare"))
#   cat(paste0(rep("-", 44), collapse = ""), "\n")
#   for (r in seq_len(nrow(tbl))) {
#     cat(sprintf("%-6d | %5.3f (+/-%.3f) | %5.3f (+/-%.3f)\n",
#                 tbl$n[r],
#                 tbl$DCP_FDR[r], tbl$DCP_FDR_SE[r],
#                 tbl$CC_FDR[r], tbl$CC_FDR_SE[r]))
#   }
#   cat("\n")
# }
# cat("--- DM: Type I Error (no true DM effects) ---\n")
# cat(sprintf("%-6s | %-16s | %-16s\n", "n", "DCP", "CircaCompare"))
# cat(paste0(rep("-", 44), collapse = ""), "\n")
# for (r in seq_len(nrow(dm_table))) {
#   cat(sprintf("%-6d | %6.4f (+/-%.4f) | %6.4f (+/-%.4f)\n",
#               dm_table$n[r],
#               dm_table$DCP_TypeI[r], dm_table$DCP_TypeI_SE[r],
#               dm_table$CC_TypeI[r], dm_table$CC_TypeI_SE[r]))
# }
# cat("\nNote: DR comparison not shown (CircaCompare has no DR test).\n")
# cat("Note: DM shows type I error only (current simulation has no mesor differences).\n\n")
#
# # --- 10E: Save results as RDS ---
# compare_results <- list(
#   dp_dcp = dp_dcp, dp_cc = dp_cc,
#   dp_table = dp_table, dm_table = dm_table,
#   design = opts_design_compare
# )
# saveRDS(compare_results, file = file.path(compare_base, "comparison_results.rds"))
# cat(sprintf("Results saved: %s\n", file.path(compare_base, "comparison_results.rds")))
#
# # --- 10F: Export tables to xlsx ---
# library(openxlsx)
# wb <- createWorkbook()
# addWorksheet(wb, "DP_Comparison")
# writeData(wb, "DP_Comparison", dp_table)
# addWorksheet(wb, "DM_TypeI_Error")
# writeData(wb, "DM_TypeI_Error", dm_table)
# xlsx_file <- file.path(compare_base, "method_comparison_tables.xlsx")
# saveWorkbook(wb, xlsx_file, overwrite = TRUE)
# cat(sprintf("Tables saved: %s\n", xlsx_file))
#
# # --- 10G: Generate comparison figures ---
# cat("\nGenerating comparison figures...\n")
#
# .to_long <- function(tbl, metric_DCP, metric_CC, metric_name,
#                      se_DCP = NULL, se_CC = NULL) {
#   df_dcp <- data.frame(n = tbl$n, value = tbl[[metric_DCP]],
#     se = if (!is.null(se_DCP)) tbl[[se_DCP]] else 0,
#     Method = "DCP (LRT)", stringsAsFactors = FALSE)
#   df_cc <- data.frame(n = tbl$n, value = tbl[[metric_CC]],
#     se = if (!is.null(se_CC)) tbl[[se_CC]] else 0,
#     Method = "CircaCompare (Wald NLS)", stringsAsFactors = FALSE)
#   out <- rbind(df_dcp, df_cc)
#   out$metric <- metric_name
#   out
# }
#
# # Figure: Power comparison (DP, 1 panel)
# dp_pow_long <- .to_long(dp_table, "DCP_Power", "CC_Power", "DP",
#                         "DCP_Power_SE", "CC_Power_SE")
# pdf(file.path(compare_dir, "power_comparison.pdf"), width = 6, height = 5)
# p_power <- ggplot(dp_pow_long, aes(x = n, y = value, color = Method, shape = Method)) +
#   geom_line(linewidth = 0.8) + geom_point(size = 3) +
#   geom_errorbar(aes(ymin = pmax(value - 1.96*se, 0), ymax = pmin(value + 1.96*se, 1)),
#                 width = 2, linewidth = 0.5) +
#   scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
#   scale_color_manual(values = c("DCP (LRT)" = "#2166AC", "CircaCompare (Wald NLS)" = "#B2182B")) +
#   labs(x = "Sample size per group", y = "Power (FDR 5%)",
#        title = "Power: DCP vs CircaCompare (DP)") +
#   theme_bw(base_size = 12) + theme(legend.position = "bottom")
# print(p_power)
# dev.off()
# cat(sprintf("  Figure: %s\n", file.path(compare_dir, "power_comparison.pdf")))
#
# # Figure: FDR control (DP)
# dp_fdr_long <- .to_long(dp_table, "DCP_FDR", "CC_FDR", "DP",
#                         "DCP_FDR_SE", "CC_FDR_SE")
# pdf(file.path(compare_dir, "fdr_comparison.pdf"), width = 6, height = 5)
# p_fdr <- ggplot(dp_fdr_long, aes(x = n, y = value, color = Method, shape = Method)) +
#   geom_line(linewidth = 0.8) + geom_point(size = 3) +
#   geom_errorbar(aes(ymin = pmax(value - 1.96*se, 0), ymax = value + 1.96*se),
#                 width = 2, linewidth = 0.5) +
#   geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey40") +
#   scale_y_continuous(limits = c(0, NA)) +
#   scale_color_manual(values = c("DCP (LRT)" = "#2166AC", "CircaCompare (Wald NLS)" = "#B2182B")) +
#   labs(x = "Sample size per group", y = "Empirical FDR",
#        title = "FDR Control: DCP vs CircaCompare (nominal 5%)") +
#   theme_bw(base_size = 12) + theme(legend.position = "bottom")
# print(p_fdr)
# dev.off()
# cat(sprintf("  Figure: %s\n", file.path(compare_dir, "fdr_comparison.pdf")))
#
# # Figure: DM Type I error
# dm_long <- .to_long(dm_table, "DCP_TypeI", "CC_TypeI", "DM",
#                     "DCP_TypeI_SE", "CC_TypeI_SE")
# pdf(file.path(compare_dir, "dm_type1_comparison.pdf"), width = 6, height = 5)
# p_dm <- ggplot(dm_long, aes(x = n, y = value, color = Method, shape = Method)) +
#   geom_line(linewidth = 0.8) + geom_point(size = 3) +
#   geom_errorbar(aes(ymin = pmax(value - 1.96*se, 0), ymax = value + 1.96*se),
#                 width = 2, linewidth = 0.5) +
#   geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey40") +
#   scale_y_continuous(limits = c(0, NA)) +
#   scale_color_manual(values = c("DCP (LRT)" = "#2166AC", "CircaCompare (Wald NLS)" = "#B2182B")) +
#   labs(x = "Sample size per group", y = "Type I Error Rate",
#        title = "DM Type I Error (no true mesor differences)") +
#   theme_bw(base_size = 12) + theme(legend.position = "bottom")
# print(p_dm)
# dev.off()
# cat(sprintf("  Figure: %s\n", file.path(compare_dir, "dm_type1_comparison.pdf")))

# =====================================================================
# SECTION 9: WRAP-UP
# =====================================================================

t_elapsed <- proc.time() - t_start
cat(sprintf("\n====================================================================\n"))
cat(sprintf("PIPELINE COMPLETE\n"))
cat(sprintf("====================================================================\n\n"))
cat(sprintf("Total runtime: %.1f minutes (%.1f hours)\n", t_elapsed[3] / 60, t_elapsed[3] / 3600))
cat(sprintf("\nOutput directory: %s/\n", base_out))
cat(sprintf("\nOutput files:\n"))
cat(sprintf("  Results:  %s\n", dr_results_file))
cat(sprintf("  Results:  %s\n", dp_results_file))
if (exists("ps_results_file")) cat(sprintf("  Results:  %s\n", ps_results_file))
cat(sprintf("  Figure:   %s\n", dr_fig))
cat(sprintf("  Figure:   %s\n", dp_fig))
if (exists("ps_fig")) cat(sprintf("  Figure:   %s\n", ps_fig))

# if (exists("compare_base")) {
#   cat(sprintf("  Results:  %s\n", file.path(compare_base, "dp_dcp.rds")))
#   cat(sprintf("  Results:  %s\n", file.path(compare_base, "dp_cc.rds")))
#   cat(sprintf("  Results:  %s\n", file.path(compare_base, "comparison_results.rds")))
#   if (exists("xlsx_file")) cat(sprintf("  Tables:   %s\n", xlsx_file))
#   cat(sprintf("  Figure:   %s\n", file.path(compare_base, "figures/power_comparison.pdf")))
#   cat(sprintf("  Figure:   %s\n", file.path(compare_base, "figures/fdr_comparison.pdf")))
#   cat(sprintf("  Figure:   %s\n", file.path(compare_base, "figures/dm_type1_comparison.pdf")))
# }

cat("\n")
sessionInfo()
