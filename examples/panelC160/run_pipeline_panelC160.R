#' =======================================================================
#' End-to-End Power Analysis (Panel C extended to n=160)
#' =======================================================================
#'
#' Runs DR/DP up to n=160, but plots only Panel C with n up to 160.
#' Other panels are limited to n <= 100.
#'
#' USAGE:
#'   Rscript examples/run_pipeline_panelC160.R
#'
#' @author Thien Pham

# =====================================================================
# SECTION 1: SETUP & CONFIGURATION
# =====================================================================

set.seed(12345)

setwd("/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")
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

COMBINED <- readRDS("/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Collaborative/Paper/Congruence/PNAS_aging/data/combined_data.rds")

expr_sample_names <- colnames(COMBINED$expr)
pheno_order <- match(expr_sample_names, COMBINED$pheno$sample_name)
valid_samples <- !is.na(pheno_order)
expr_sample_names <- expr_sample_names[valid_samples]
COMBINED$expr <- COMBINED$expr[, valid_samples]
pheno_order <- pheno_order[valid_samples]
pheno_data <- COMBINED$pheno[pheno_order, ]
pheno_data$tod <- if ("TOD.x" %in% colnames(pheno_data)) pheno_data$TOD.x else pheno_data$TOD.y
pheno_data$age_group_final <- if ("AgeGroup" %in% colnames(pheno_data)) pheno_data$AgeGroup else pheno_data$age_group
complete_samples <- !is.na(pheno_data$age_group_final) & !is.na(pheno_data$tod) &
  pheno_data$age_group_final %in% c("younger", "older")
pheno_clean <- pheno_data[complete_samples, ]
younger_idx <- pheno_clean$age_group_final == "younger"
expr_younger <- COMBINED$expr[, complete_samples][, younger_idx]
times_young <- pheno_clean$tod[younger_idx]

cat(sprintf("  Expression matrix: %d genes x %d samples\n", nrow(expr_younger), ncol(expr_younger)))
cat(sprintf("  Reference TOD: n=%d younger subjects\n\n", length(times_young)))
rm(COMBINED, expr_sample_names, pheno_order, valid_samples, pheno_data,
   complete_samples, pheno_clean, younger_idx)


# =====================================================================
# SECTION 3: ESTIMATE PARAMETERS FROM PILOT DATA
# =====================================================================

cat("Estimating circadian parameters from pilot data...\n\n")

opts_bio <- estCircadianParam(
  data          = expr_younger,
  times         = times_young,
  period        = 24,
  prop_DR       = 0.15,
  prop_DP       = 0.10,
  phase_diff    = c(-6, 6),
  amp_diff      = c(2, 4)
)

rm(expr_younger)
opts_bio <- updateBioOptions(opts_bio, ngenes = 5000)


# =====================================================================
# SECTION 4: DESIGN & ANALYSIS CONFIGURATION
# =====================================================================

opts_design <- CircadianDesignOptions(
  sample_sizes = c(20, 40, 60, 80, 100, 120, 140, 160),
  nsims        = 50,
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
  reference_n     = 60
)

base_out <- file.path("output", "run_panelC160")
out_dir  <- file.path(base_out, "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


# =====================================================================
# SECTION 5: CUSTOM PLOTTING (Panel C up to n=160 only)
# =====================================================================

plotWithSE_panelC160 <- function(results_obj, output_file, test_name = "DR",
                                 analysis.opts = NULL,
                                 panel_c_max_n = 160,
                                 other_max_n = 100) {
  if (exists("dp_power_raw")) results_obj <- dp_power_raw
  dr_power_raw <- results_obj

  sample_sizes  <- dr_power_raw$sample_sizes
  nsims         <- dr_power_raw$nsims
  strata_labels <- dr_power_raw$strata_labels
  r_strata      <- dr_power_raw$r_strata

  if (!is.null(analysis.opts) && inherits(analysis.opts, "CircadianAnalysisOptions")) {
    fdr_thresholds <- analysis.opts$fdr_thresholds
    reference_n    <- analysis.opts$reference_n
  } else {
    fdr_thresholds <- c(0.01, 0.05, 0.10, 0.20)
    reference_n    <- 60
  }
  n_sizes      <- length(sample_sizes)
  n_thresholds <- length(fdr_thresholds)
  n_strata     <- length(strata_labels)

  nested_gt <- is.list(dr_power_raw$is_target_list[[1]]) && length(dr_power_raw$is_target_list) == n_sizes

  power_arr <- array(NA, dim = c(n_sizes, n_strata, n_thresholds, nsims))
  TD_arr    <- array(NA, dim = c(n_sizes, n_strata, n_thresholds, nsims))
  FD_arr    <- array(NA, dim = c(n_sizes, n_strata, n_thresholds, nsims))

  for (j in 1:n_sizes) {
    for (s in 1:nsims) {
      pvals <- dr_power_raw$pvalues[j, , s]

      if (nested_gt) {
        r_vec     <- dr_power_raw$r_values_list[[j]][[s]]
        is_target <- dr_power_raw$is_target_list[[j]][[s]]
        is_null   <- dr_power_raw$is_null_list[[j]][[s]]
      } else {
        r_vec     <- dr_power_raw$r_values_list[[s]]
        is_target <- dr_power_raw$is_target_list[[s]]
        is_null   <- dr_power_raw$is_null_list[[s]]
      }

      xgr <- cut(r_vec, breaks = r_strata, include.lowest = TRUE, labels = FALSE)

      qvals <- rep(1, length(pvals))
      tested <- pvals < 1
      if (sum(tested) > 0) qvals[tested] <- p.adjust(pvals[tested], method = "BH")

      for (t in 1:n_thresholds) {
        discoveries <- qvals <= fdr_thresholds[t]
        for (k in 1:n_strata) {
          in_stratum <- xgr == k
          td <- sum(discoveries & is_target & in_stratum, na.rm = TRUE)
          fd <- sum(discoveries & is_null & in_stratum, na.rm = TRUE)
          nt <- sum(is_target & in_stratum, na.rm = TRUE)
          TD_arr[j, k, t, s] <- td
          FD_arr[j, k, t, s] <- fd
          power_arr[j, k, t, s] <- if (nt > 0) td / nt else NA
        }
      }
    }
  }

  # Summary across sims
  marginal_mean <- matrix(NA, nrow = n_sizes, ncol = n_thresholds)
  marginal_se   <- matrix(NA, nrow = n_sizes, ncol = n_thresholds)

  for (j in 1:n_sizes) {
    for (t in 1:n_thresholds) {
      td <- apply(TD_arr[j, , t, ], 2, sum, na.rm = TRUE)
      nt <- apply(dr_power_raw$strat_n_targets[j, , ], 2, sum, na.rm = TRUE)
      power <- ifelse(nt > 0, td / nt, NA_real_)
      marginal_mean[j, t] <- mean(power, na.rm = TRUE)
      marginal_se[j, t]   <- sd(power, na.rm = TRUE) / sqrt(sum(!is.na(power)))
    }
  }

  # Indices for plotting
  idx_other <- which(sample_sizes <= other_max_n)
  idx_c     <- which(sample_sizes <= panel_c_max_n)

  sample_sizes_other <- sample_sizes[idx_other]
  sample_sizes_c     <- sample_sizes[idx_c]

  # Colors
  threshold_colors <- c("darkgreen", "steelblue", "orange", "red")
  threshold_labels <- paste0("FDR ", fdr_thresholds * 100, "%")
  size_colors <- rainbow(length(sample_sizes_other), s = 0.6, v = 0.8)

  pdf(output_file, width = 12, height = 8)
  par(mfrow = c(2, 3), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))

  # --- Panel A: Power by r (FDR 5%), lines per n (<=100) ---
  idx_fdr5 <- 2
  mean_power_A <- apply(power_arr[idx_other, , idx_fdr5, ], c(1, 2), mean, na.rm = TRUE)
  se_power_A <- apply(power_arr[idx_other, , idx_fdr5, ], c(1, 2),
                      function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))

  y_max_A <- max(mean_power_A + se_power_A, na.rm = TRUE)
  matplot(1:n_strata, t(mean_power_A),
          type = "l", lwd = 2, col = size_colors, lty = 1,
          xlim = c(0.5, n_strata + 0.5), ylim = c(0, y_max_A),
          xlab = expression(r == A/sigma), ylab = "Power",
          main = sprintf("%s Power by r (FDR 5%%)", test_name), xaxt = "n")
  axis(1, at = 1:n_strata, labels = strata_labels, las = 2, cex.axis = 0.6)
  for (j in 1:length(sample_sizes_other)) {
    points(1:n_strata, mean_power_A[j, ], pch = 19, col = size_colors[j], cex = 0.6)
    add_se_bars(1:n_strata, mean_power_A[j, ], se_power_A[j, ], col = size_colors[j])
  }
  grid()
  legend("bottomright", paste0("n=", sample_sizes_other), col = size_colors, lty = 1, lwd = 2, cex = 0.6)
  mtext("A", side = 3, at = -0.02, font = 2)

  # --- Panel B: Power by r at different FDR thresholds (n=60) ---
  idx_n60 <- which.min(abs(sample_sizes_other - reference_n))
  j0 <- idx_other[idx_n60]
  mean_power_B <- apply(power_arr[j0, , , ], c(1, 2), mean, na.rm = TRUE)
  se_power_B <- apply(power_arr[j0, , , ], c(1, 2),
                      function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
  y_max_B <- max(mean_power_B + se_power_B, na.rm = TRUE)
  matplot(1:n_strata, mean_power_B,
          type = "l", lwd = 2, col = threshold_colors, lty = 1,
          xlim = c(0.5, n_strata + 0.5), ylim = c(0, y_max_B),
          xlab = expression(r == A/sigma), ylab = "Power",
          main = sprintf("%s Power by r (n=%d)", test_name, reference_n), xaxt = "n")
  axis(1, at = 1:n_strata, labels = strata_labels, las = 2, cex.axis = 0.6)
  for (t in 1:n_thresholds) {
    points(1:n_strata, mean_power_B[, t], pch = 19, col = threshold_colors[t], cex = 0.6)
    add_se_bars(1:n_strata, mean_power_B[, t], se_power_B[, t], col = threshold_colors[t])
  }
  grid()
  legend("bottomright", threshold_labels, col = threshold_colors, lty = 1, pch = 19, lwd = 2, cex = 0.7)
  mtext("B", side = 3, at = -0.02, font = 2)

  # --- Panel C: Marginal Power vs Sample Size (FDR 5%, extended to n=160) ---
  plot(sample_sizes_c, 100 * marginal_mean[idx_c, idx_fdr5],
       type = "b", pch = 19, lwd = 2, col = "steelblue",
       xlim = c(0, max(sample_sizes_c) * 1.1), ylim = c(0, 100),
       xlab = "Sample Size (per group)", ylab = "Marginal Power (%)",
       xaxt = "n",
       main = sprintf("%s Marginal Power vs Sample Size (FDR 5%%)", test_name))
  axis(1, at = sample_sizes_c, labels = sample_sizes_c)
  add_se_bars(sample_sizes_c, 100 * marginal_mean[idx_c, idx_fdr5],
              100 * marginal_se[idx_c, idx_fdr5], col = "steelblue")
  abline(h = 80, lty = 2, col = "gray")
  abline(v = 60, lty = 3, col = "darkgreen", lwd = 1.5)
  grid()
  text(sample_sizes_c, 100 * marginal_mean[idx_c, idx_fdr5] + 4,
       sprintf("%.1f%%", 100 * marginal_mean[idx_c, idx_fdr5]), cex = 0.65, col = "steelblue")
  text(63, 5, "n=60", cex = 0.6, col = "darkgreen", pos = 4)
  mtext("C", side = 3, at = -0.02, font = 2)

  # --- Panel D: Marginal Power vs Sample Size (all FDR thresholds, <=100) ---
  matplot(sample_sizes_other, 100 * marginal_mean[idx_other, ],
          type = "b", pch = 19, lwd = 2, col = threshold_colors, lty = 1,
          xlim = c(0, max(sample_sizes_other) * 1.1), ylim = c(0, 100),
          xlab = "Sample Size (per group)", ylab = "Marginal Power (%)",
          main = sprintf("%s Power vs Sample Size by FDR", test_name))
  for (t in 1:n_thresholds) {
    add_se_bars(sample_sizes_other, 100 * marginal_mean[idx_other, t],
                100 * marginal_se[idx_other, t], col = threshold_colors[t])
  }
  abline(h = 80, lty = 2, col = "gray")
  grid()
  legend("bottomright", threshold_labels, col = threshold_colors, lty = 1, pch = 19, lwd = 2, cex = 0.7)
  mtext("D", side = 3, at = -0.02, font = 2)

  # --- Panel E: True Discoveries by r stratum (FDR 5%, lines per n <=100) ---
  mean_TD_E <- apply(dr_power_raw$strat_TD[idx_other, , ], c(1, 2), mean, na.rm = TRUE)
  se_TD_E   <- apply(dr_power_raw$strat_TD[idx_other, , ], c(1, 2),
                     function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))

  y_max_E <- max(mean_TD_E + se_TD_E, na.rm = TRUE) * 1.1
  matplot(1:n_strata, t(mean_TD_E),
          type = "l", lwd = 2, col = size_colors, lty = 1,
          xlim = c(0.5, n_strata + 0.5), ylim = c(0, y_max_E),
          xlab = expression(r == A/sigma), ylab = "Mean True Discoveries (per sim)",
          main = sprintf("%s True Discoveries by r (FDR 5%%)", test_name), xaxt = "n")
  axis(1, at = 1:n_strata, labels = strata_labels, las = 2, cex.axis = 0.6)
  for (j in 1:length(sample_sizes_other)) {
    points(1:n_strata, mean_TD_E[j, ], pch = 19, col = size_colors[j], cex = 0.6)
    add_se_bars(1:n_strata, mean_TD_E[j, ], se_TD_E[j, ], col = size_colors[j])
  }
  grid()
  legend("topright", paste0("n=", sample_sizes_other), col = size_colors, lty = 1, lwd = 2, cex = 0.6)
  mtext("E", side = 3, at = -0.02, font = 2)

  # --- Panel F: Distribution of target genes across r strata ---
  if (nested_gt) {
    r_vec     <- dr_power_raw$r_values_list[[1]][[1]]
    is_target <- dr_power_raw$is_target_list[[1]][[1]]
  } else {
    r_vec     <- dr_power_raw$r_values_list[[1]]
    is_target <- dr_power_raw$is_target_list[[1]]
  }
  xgr <- cut(r_vec[is_target], breaks = r_strata, include.lowest = TRUE, labels = FALSE)
  target_counts <- tabulate(xgr, nbins = n_strata)

  barplot(target_counts, names.arg = strata_labels,
          col = "lightblue", border = "gray40",
          xlab = expression(r == A/sigma), ylab = "Number of Target Genes",
          main = sprintf("Distribution of %s Target Genes by r", test_name),
          las = 2, cex.names = 0.6)
  mtext("F", side = 3, at = -0.02, font = 2)

  mtext(sprintf("%s Power Analysis - Stratified by Signal-to-Noise Ratio", test_name),
        outer = TRUE, font = 2, cex = 1.2)
  dev.off()
  cat(sprintf("Figure saved: %s\n", output_file))
}


# =====================================================================
# SECTION 6: DR ANALYSIS
# =====================================================================

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

dr_fig <- file.path(out_dir, "dr_power.pdf")
plotWithSE_panelC160(dr_power_raw, dr_fig, test_name = "DR", analysis.opts = opts_analysis,
                     panel_c_max_n = 160, other_max_n = 100)


# =====================================================================
# SECTION 7: DP ANALYSIS
# =====================================================================

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

dp_fig <- file.path(out_dir, "dp_power.pdf")
plotWithSE_panelC160(dp_power_raw, dp_fig, test_name = "DP", analysis.opts = opts_analysis,
                     panel_c_max_n = 160, other_max_n = 100)

# =====================================================================
# SECTION 8: DA ANALYSIS (extended to n=160)
# =====================================================================

cat("====================================================================\n")
cat("ANALYSIS 3: DIFFERENTIAL AMPLITUDE (DA)\n")
cat("====================================================================\n\n")

opts_bio_DA <- updateBioOptions(opts_bio,
  prop_DR    = 0.00,
  prop_DP    = 0.00,
  phase_diff = c(0, 0),
  amp_diff   = c(2, 4)
)

da_power_raw <- runPowerAnalysis(opts_bio_DA, opts_design, opts_analysis,
                                 test_type = "DA")

da_results_file <- file.path(base_out, "da_power_raw_pvalues.rds")
save(da_power_raw, file = da_results_file)
cat(sprintf("\nDA results saved: %s\n", da_results_file))

da_fig <- file.path(out_dir, "da_power.pdf")
plotWithSE_panelC160(da_power_raw, da_fig, test_name = "DA", analysis.opts = opts_analysis,
                     panel_c_max_n = 160, other_max_n = 100)

cat("\nDONE. Output directory: ", base_out, "\n")
