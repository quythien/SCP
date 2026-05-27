#!/usr/bin/env Rscript
#' =======================================================================
#' 05_method_comparison.R
#'    DCP vs CircaCompare Power & FDR Comparison
#' =======================================================================
#'
#' PURPOSE
#'   Compares two differential circadian testing methods on identical
#'   simulated data drawn from the human PFC younger pilot:
#'     DCP          — LR tests (asymptotically optimal)
#'     CircaCompare — Wald t-tests from per-gene NLS
#'
#'   Comparable tests:
#'     DP:  DCP delta-phi LRT  vs  CircaCompare phi1 Wald
#'     DM:  DCP delta-M LRT   vs  CircaCompare k1 Wald  (type I only)
#'     DR:  DCP only (CircaCompare has no DR test)
#'
#' NOTE
#'   CircaCompare is ~1-3 sec/gene (NLS). NGENES_CMP and NSIMS_CMP are
#'   deliberately reduced relative to 03_power_core.R.
#'   Smoke test: POWERSIM_SMOKE=1 Rscript examples/publication/05_method_comparison.R
#'
#' OUTPUTS
#'   output/05_method_comparison_<ts>/
#'     dp_dcp.rds
#'     dp_cc.rds
#'     comparison_results.rds
#'     method_comparison_tables.xlsx
#'     figures/
#'       power_comparison.pdf
#'       fdr_comparison.pdf
#'       dm_type1_comparison.pdf
#'
#' @author Thien Quy Pham

# =====================================================================
# SECTION 1: SETUP & CONFIGURATION
# =====================================================================

SMOKE       <- nzchar(Sys.getenv("POWERSIM_SMOKE"))
NSIMS_CMP   <- if (SMOKE) 5L   else 20L
NGENES_CMP  <- if (SMOKE) 100L else 500L
N_GRID_CMP  <- if (SMOKE) c(20L, 40L) else c(20L, 40L, 60L, 80L, 100L)

POWERSIM_ROOT <- {
  env <- Sys.getenv("POWERSIM_ROOT", unset = "")
  if (nzchar(env)) env else
    "/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim"
}

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

cat(sprintf("Mode: %s | NGENES=%d | NSIMS=%d\n",
            if (SMOKE) "SMOKE" else "PRODUCTION",
            NGENES_CMP, NSIMS_CMP))
cat("NOTE: CircaCompare uses NLS (~1-3 sec/gene) — runtime will be long.\n\n")

run_tag     <- format(Sys.time(), "%Y%m%d_%H%M")
base_out    <- file.path("output", paste0("05_method_comparison_", run_tag))
compare_dir <- file.path(base_out, "figures")
dir.create(compare_dir, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("Output: %s/\n\n", base_out))


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

prep_young   <- prepCircadianData(expr_raw_young, times = times_raw, input_type = "log2")
expr_younger <- prep_young$data
times_young  <- prep_young$times
rm(expr_raw_young, times_raw)

cat(sprintf("  Expression matrix: %d genes x %d samples\n", nrow(expr_younger), ncol(expr_younger)))
cat(sprintf("  Reference TOD: n=%d younger subjects\n\n", length(times_young)))


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
  amp_diff      = c(0.5, 2)
)
rm(expr_younger)

opts_bio <- updateBioOptions(opts_bio, ngenes = NGENES_CMP)
cat(sprintf("  Simulation ngenes=%d (reduced for CircaCompare speed)\n\n", NGENES_CMP))


# =====================================================================
# SECTION 4: ANALYSIS OPTIONS
# =====================================================================

opts_analysis_DCP <- CircadianAnalysisOptions(
  alpha = 0.05, p.adjust.method = "BH", DCmethod = "DCP"
)
opts_analysis_CC <- CircadianAnalysisOptions(
  alpha = 0.05, p.adjust.method = "BH", DCmethod = "CircaCompare"
)

opts_design_compare <- CircadianDesignOptions(
  sample_sizes = N_GRID_CMP,
  nsims        = NSIMS_CMP,
  design       = "passive",
  cts          = times_young,
  test_types   = c("DP", "DM")
)


# =====================================================================
# SECTION 5: DP COMPARISON (Differential Phase)
# =====================================================================

cat("====================================================================\n")
cat("SECTION 5: Differential Phase (DP) — DCP vs CircaCompare\n")
cat("====================================================================\n\n")
cat("  Setting up: 15% DP genes, phase_diff = [-6, 6]\n\n")

opts_bio_DP_cmp <- updateBioOptions(opts_bio,
  prop_DR    = 0.00,
  prop_DP    = 0.15,
  phase_diff = c(-6, 6),
  amp_diff   = c(0.5, 2)
)

cat("  Running DCP...\n")
dp_dcp <- runSimsDiff(opts_bio_DP_cmp, opts_design_compare, opts_analysis_DCP)
saveRDS(dp_dcp, file.path(base_out, "dp_dcp.rds"))
cat(sprintf("  Saved: %s\n", file.path(base_out, "dp_dcp.rds")))

cat("  Running CircaCompare...\n")
dp_cc <- runSimsDiff(opts_bio_DP_cmp, opts_design_compare, opts_analysis_CC)
saveRDS(dp_cc, file.path(base_out, "dp_cc.rds"))
cat(sprintf("  Saved: %s\n\n", file.path(base_out, "dp_cc.rds")))


# =====================================================================
# SECTION 6: COMPUTE COMPARISON STATISTICS
# =====================================================================

cat("Computing comparison statistics...\n")

.compute_comparison <- function(dcp_out, cc_out, test_type) {
  fdr_key <- paste0("fdr_", test_type)
  sample_sizes_cmp <- dcp_out$sample_sizes
  nsims_cmp <- dcp_out$nsims
  ngenes_cmp <- dcp_out$ngenes
  rows <- list()
  for (j in seq_along(sample_sizes_cmp)) {
    dcp_power_vals <- numeric(nsims_cmp)
    cc_power_vals  <- numeric(nsims_cmp)
    dcp_fdr_vals   <- numeric(nsims_cmp)
    cc_fdr_vals    <- numeric(nsims_cmp)
    dcp_td_vals    <- numeric(nsims_cmp)
    cc_td_vals     <- numeric(nsims_cmp)
    dcp_fd_vals    <- numeric(nsims_cmp)
    cc_fd_vals     <- numeric(nsims_cmp)
    for (i in seq_len(nsims_cmp)) {
      diff_type <- dcp_out$diff_type[[i]]
      is_target <- if (test_type == "DP") diff_type == 4 else rep(FALSE, ngenes_cmp)
      is_null   <- !is_target
      n_target  <- sum(is_target)
      dcp_disc <- dcp_out[[fdr_key]][, j, i] <= 0.05
      cc_disc  <- cc_out[[fdr_key]][, j, i] <= 0.05
      dcp_td <- sum(dcp_disc & is_target, na.rm = TRUE)
      cc_td  <- sum(cc_disc  & is_target, na.rm = TRUE)
      dcp_fd <- sum(dcp_disc & is_null,   na.rm = TRUE)
      cc_fd  <- sum(cc_disc  & is_null,   na.rm = TRUE)
      dcp_td_vals[i] <- dcp_td; cc_td_vals[i] <- cc_td
      dcp_fd_vals[i] <- dcp_fd; cc_fd_vals[i] <- cc_fd
      dcp_power_vals[i] <- if (n_target > 0) dcp_td / n_target else NA
      cc_power_vals[i]  <- if (n_target > 0) cc_td  / n_target else NA
      dcp_total <- dcp_td + dcp_fd; cc_total <- cc_td + cc_fd
      dcp_fdr_vals[i] <- if (dcp_total > 0) dcp_fd / dcp_total else NA
      cc_fdr_vals[i]  <- if (cc_total  > 0) cc_fd  / cc_total  else NA
    }
    rows[[j]] <- data.frame(
      test_type    = test_type,
      n            = sample_sizes_cmp[j],
      DCP_Power    = mean(dcp_power_vals, na.rm = TRUE),
      CC_Power     = mean(cc_power_vals,  na.rm = TRUE),
      DCP_Power_SE = sd(dcp_power_vals, na.rm = TRUE) / sqrt(sum(!is.na(dcp_power_vals))),
      CC_Power_SE  = sd(cc_power_vals,  na.rm = TRUE) / sqrt(sum(!is.na(cc_power_vals))),
      DCP_FDR      = mean(dcp_fdr_vals, na.rm = TRUE),
      CC_FDR       = mean(cc_fdr_vals,  na.rm = TRUE),
      DCP_FDR_SE   = sd(dcp_fdr_vals, na.rm = TRUE) / sqrt(sum(!is.na(dcp_fdr_vals))),
      CC_FDR_SE    = sd(cc_fdr_vals,  na.rm = TRUE) / sqrt(sum(!is.na(cc_fdr_vals))),
      DCP_Avg_TD   = mean(dcp_td_vals),
      CC_Avg_TD    = mean(cc_td_vals),
      DCP_Avg_FD   = mean(dcp_fd_vals),
      CC_Avg_FD    = mean(cc_fd_vals),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

.compute_dm_typeI <- function(dcp_out, cc_out) {
  sample_sizes_cmp <- dcp_out$sample_sizes
  nsims_cmp <- dcp_out$nsims
  ngenes_cmp <- dcp_out$ngenes
  rows <- list()
  for (j in seq_along(sample_sizes_cmp)) {
    dcp_vals <- numeric(nsims_cmp)
    cc_vals  <- numeric(nsims_cmp)
    for (i in seq_len(nsims_cmp)) {
      dcp_vals[i] <- sum(dcp_out$fdr_DM[, j, i] <= 0.05, na.rm = TRUE) / ngenes_cmp
      cc_vals[i]  <- sum(cc_out$fdr_DM[, j, i]  <= 0.05, na.rm = TRUE) / ngenes_cmp
    }
    rows[[j]] <- data.frame(
      test_type    = "DM",
      n            = sample_sizes_cmp[j],
      DCP_TypeI    = mean(dcp_vals),
      CC_TypeI     = mean(cc_vals),
      DCP_TypeI_SE = sd(dcp_vals) / sqrt(nsims_cmp),
      CC_TypeI_SE  = sd(cc_vals)  / sqrt(nsims_cmp),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

dp_table <- .compute_comparison(dp_dcp, dp_cc, "DP")
dm_table <- .compute_dm_typeI(dp_dcp, dp_cc)


# =====================================================================
# SECTION 7: PRINT TABLES
# =====================================================================

cat("\n====================================================================\n")
cat("METHOD COMPARISON RESULTS\n")
cat("====================================================================\n\n")

cat("--- DP: Power (FDR 5%) ---\n")
cat(sprintf("%-6s | %-20s | %-20s\n", "n", "DCP (LRT)", "CircaCompare (Wald)"))
cat(paste0(rep("-", 52), collapse = ""), "\n")
for (r in seq_len(nrow(dp_table))) {
  cat(sprintf("%-6d | %5.1f%% (+/-%.1f%%)      | %5.1f%% (+/-%.1f%%)\n",
              dp_table$n[r],
              100 * dp_table$DCP_Power[r], 100 * dp_table$DCP_Power_SE[r],
              100 * dp_table$CC_Power[r],  100 * dp_table$CC_Power_SE[r]))
}

cat("\n--- DP: Empirical FDR (nominal 5%) ---\n")
cat(sprintf("%-6s | %-20s | %-20s\n", "n", "DCP (LRT)", "CircaCompare (Wald)"))
cat(paste0(rep("-", 52), collapse = ""), "\n")
for (r in seq_len(nrow(dp_table))) {
  cat(sprintf("%-6d | %5.3f (+/-%.3f)       | %5.3f (+/-%.3f)\n",
              dp_table$n[r],
              dp_table$DCP_FDR[r], dp_table$DCP_FDR_SE[r],
              dp_table$CC_FDR[r],  dp_table$CC_FDR_SE[r]))
}

cat("\n--- DM: Type I Error (no true mesor differences) ---\n")
cat(sprintf("%-6s | %-20s | %-20s\n", "n", "DCP (LRT)", "CircaCompare (Wald)"))
cat(paste0(rep("-", 52), collapse = ""), "\n")
for (r in seq_len(nrow(dm_table))) {
  cat(sprintf("%-6d | %6.4f (+/-%.4f)     | %6.4f (+/-%.4f)\n",
              dm_table$n[r],
              dm_table$DCP_TypeI[r], dm_table$DCP_TypeI_SE[r],
              dm_table$CC_TypeI[r],  dm_table$CC_TypeI_SE[r]))
}
cat("\nNote: DR comparison not shown (CircaCompare has no DR test).\n")
cat("Note: DM shows type I error only (simulation has no true mesor differences).\n\n")


# =====================================================================
# SECTION 8: SAVE RESULTS
# =====================================================================

compare_results <- list(
  dp_dcp = dp_dcp, dp_cc = dp_cc,
  dp_table = dp_table, dm_table = dm_table,
  design = opts_design_compare
)
saveRDS(compare_results, file.path(base_out, "comparison_results.rds"))
cat(sprintf("Results saved: %s\n", file.path(base_out, "comparison_results.rds")))

if (requireNamespace("openxlsx", quietly = TRUE)) {
  library(openxlsx)
  wb <- createWorkbook()
  addWorksheet(wb, "DP_Comparison")
  writeData(wb, "DP_Comparison", dp_table)
  addWorksheet(wb, "DM_TypeI_Error")
  writeData(wb, "DM_TypeI_Error", dm_table)
  xlsx_file <- file.path(base_out, "method_comparison_tables.xlsx")
  saveWorkbook(wb, xlsx_file, overwrite = TRUE)
  cat(sprintf("Tables saved: %s\n", xlsx_file))
} else {
  cat("openxlsx not available — skipping xlsx export.\n")
}


# =====================================================================
# SECTION 9: FIGURES
# =====================================================================

cat("\nGenerating comparison figures...\n")

.to_long <- function(tbl, metric_DCP, metric_CC, metric_name,
                     se_DCP = NULL, se_CC = NULL) {
  df_dcp <- data.frame(n = tbl$n, value = tbl[[metric_DCP]],
    se = if (!is.null(se_DCP)) tbl[[se_DCP]] else 0,
    Method = "DCP (LRT)", stringsAsFactors = FALSE)
  df_cc <- data.frame(n = tbl$n, value = tbl[[metric_CC]],
    se = if (!is.null(se_CC)) tbl[[se_CC]] else 0,
    Method = "CircaCompare (Wald NLS)", stringsAsFactors = FALSE)
  out <- rbind(df_dcp, df_cc)
  out$metric <- metric_name
  out
}

cols_method <- c("DCP (LRT)" = "#2166AC", "CircaCompare (Wald NLS)" = "#B2182B")

# Power comparison (DP)
dp_pow_long <- .to_long(dp_table, "DCP_Power", "CC_Power", "DP",
                        "DCP_Power_SE", "CC_Power_SE")
pdf(file.path(compare_dir, "power_comparison.pdf"), width = 6, height = 5)
p_power <- ggplot2::ggplot(dp_pow_long,
    ggplot2::aes(x = n, y = value, color = Method, shape = Method)) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(size = 3) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = pmax(value - 1.96*se, 0), ymax = pmin(value + 1.96*se, 1)),
    width = 2, linewidth = 0.5) +
  ggplot2::geom_hline(yintercept = 0.80, linetype = "dashed", color = "grey40") +
  ggplot2::scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  ggplot2::scale_color_manual(values = cols_method) +
  ggplot2::labs(x = "Sample size per group", y = "Power (FDR 5%)",
       title = "Power: DCP vs CircaCompare (DP, phase_diff=[-6,6])") +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(legend.position = "bottom")
print(p_power)
dev.off()
cat(sprintf("  Figure: %s\n", file.path(compare_dir, "power_comparison.pdf")))

# FDR control (DP)
dp_fdr_long <- .to_long(dp_table, "DCP_FDR", "CC_FDR", "DP",
                        "DCP_FDR_SE", "CC_FDR_SE")
pdf(file.path(compare_dir, "fdr_comparison.pdf"), width = 6, height = 5)
p_fdr <- ggplot2::ggplot(dp_fdr_long,
    ggplot2::aes(x = n, y = value, color = Method, shape = Method)) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(size = 3) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = pmax(value - 1.96*se, 0), ymax = value + 1.96*se),
    width = 2, linewidth = 0.5) +
  ggplot2::geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey40") +
  ggplot2::scale_y_continuous(limits = c(0, NA)) +
  ggplot2::scale_color_manual(values = cols_method) +
  ggplot2::labs(x = "Sample size per group", y = "Empirical FDR",
       title = "FDR Control: DCP vs CircaCompare (nominal 5%)") +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(legend.position = "bottom")
print(p_fdr)
dev.off()
cat(sprintf("  Figure: %s\n", file.path(compare_dir, "fdr_comparison.pdf")))

# DM Type I error
dm_long <- .to_long(dm_table, "DCP_TypeI", "CC_TypeI", "DM",
                    "DCP_TypeI_SE", "CC_TypeI_SE")
pdf(file.path(compare_dir, "dm_type1_comparison.pdf"), width = 6, height = 5)
p_dm <- ggplot2::ggplot(dm_long,
    ggplot2::aes(x = n, y = value, color = Method, shape = Method)) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(size = 3) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = pmax(value - 1.96*se, 0), ymax = value + 1.96*se),
    width = 2, linewidth = 0.5) +
  ggplot2::geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey40") +
  ggplot2::scale_y_continuous(limits = c(0, NA)) +
  ggplot2::scale_color_manual(values = cols_method) +
  ggplot2::labs(x = "Sample size per group", y = "Type I Error Rate",
       title = "DM Type I Error (no true mesor differences)") +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(legend.position = "bottom")
print(p_dm)
dev.off()
cat(sprintf("  Figure: %s\n", file.path(compare_dir, "dm_type1_comparison.pdf")))


# =====================================================================
# WRAP-UP
# =====================================================================

cat("\n====================================================================\n")
cat("05_method_comparison COMPLETE\n")
cat("====================================================================\n\n")
cat(sprintf("Output: %s/\n", base_out))
cat("Done.\n")
