#' =======================================================================
#' 08b_bootstrap_d1d2.R — Two-Stage vs Bootstrap: Mouse D1D2 (n=45)
#' =======================================================================
#' Split from 08_two_stage_vs_bootstrap_realdata.R. Section 2 only.
#' See 08a_bootstrap_baboon.R header for parallel launch instructions.

POWERSIM_ROOT <- {
  env <- Sys.getenv("POWERSIM_ROOT", unset = "")
  if (nzchar(env)) env else
    "/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim"
}
setwd(POWERSIM_ROOT)

SMOKE       <- nzchar(Sys.getenv("POWERSIM_SMOKE"))
NGENES      <- if (SMOKE) 300L else 3000L
NBOOT       <- if (SMOKE) 5L  else 50L
NSIMS       <- if (SMOKE) 5L  else 30L
NSIMS_INNER <- if (SMOKE) 5L  else 20L
RHYTHM_PVAL <- 0.05

B_D1D2        <- 6L
N_GRID_D1D2   <- if (SMOKE) c(24L, 48L, 72L) else c(24L, 36L, 48L, 60L, 72L, 96L, 120L)
DATA_D1D2_PHENO <- "data/mouse_clinicalinfo_03082021_rmOutliers.csv"
DATA_D1D2_EXPR  <- "data/mouse_D1D2_logCPMfiltered_counts.csv"

RUN_TAG <- Sys.getenv("RUN_TAG", format(Sys.Date(), "%Y%m%d"))
out_dir <- file.path("output", paste0("08_two_stage_vs_bootstrap_", RUN_TAG))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("Mode: %s | NGENES=%d | NBOOT=%d | NSIMS=%d | NSIMS_INNER=%d\n",
            if (SMOKE) "SMOKE" else "PRODUCTION", NGENES, NBOOT, NSIMS, NSIMS_INNER))
cat(sprintf("Output -> %s/\n\n", out_dir))

suppressPackageStartupMessages(library(readxl))
source_dir <- file.path(getwd(), "code")
old_wd <- setwd(source_dir); source("setup.R"); setwd(old_wd)

opts_analysis <- CircadianAnalysisOptions(
  alpha = 0.05, p.adjust.method = "BH",
  fdr_thresholds = c(0.05), reference_n = 60L
)

.run_comparison <- function(pilot_data, pilot_times, bio_diff_opts,
                             N_grid, B_val, design_type, label, color,
                             out_prefix) {
  cat(sprintf("\n--- %s ---\n", label))
  cat(sprintf("    Pilot: n=%d | design=%s | B=%d\n", ncol(pilot_data), design_type, B_val))

  design_vec <- if (design_type == "active") {
    seq(0, 24, length.out = B_val + 1)[seq_len(B_val)]
  } else {
    pilot_times
  }

  cat("  Running two-stage...\n")
  design_opts_ts <- CircadianDesignOptions(
    sample_sizes = N_grid, nsims = NSIMS, design = design_type,
    cts = if (design_type == "active") design_vec else pilot_times,
    test_types = "DR"
  )
  ts_result <- tryCatch(
    runTwoStagePower(pilot_data = pilot_data, pilot_times = pilot_times,
                     design.opts = design_opts_ts, analysis.opts = opts_analysis,
                     bio_diff.opts = bio_diff_opts, min_rhythm_pval = RHYTHM_PVAL,
                     test_type = "DR", verbose = FALSE),
    error = function(e) { warning(sprintf("Two-stage failed: %s", e$message)); NULL }
  )

  cat(sprintf("  Running bootstrap (B=%d, %d draws)...\n", B_val, NBOOT))
  boot_opts <- CircadianBootstrapOptions(
    design_vector = design_vec, B_values = B_val, N_values = N_grid,
    nboot = NBOOT, nsims_inner = NSIMS_INNER, design = design_type, seed = 42L
  )
  boot_result <- tryCatch(
    runBootstrapDesignGrid(pilot_data = pilot_data, pilot_times = pilot_times,
                           boot.opts = boot_opts, analysis.opts = opts_analysis,
                           bio_diff.opts = bio_diff_opts, verbose = FALSE),
    error = function(e) { warning(sprintf("Bootstrap failed: %s", e$message)); NULL }
  )

  if (is.null(ts_result) || is.null(boot_result)) return(NULL)

  comparison <- compareDesignApproaches(
    two_stage_result = ts_result, bootstrap_result = boot_result,
    test_type = "DR", target_power = 0.80
  )

  saveRDS(list(two_stage = ts_result, boot = boot_result, comparison = comparison),
          paste0(out_prefix, ".rds"))

  comp_df   <- comparison$comparison
  ci_widths <- comp_df$boot_ci_hi - comp_df$boot_ci_lo

  cat(sprintf("  Two-stage n80:        %s\n",
              ifelse(is.na(comparison$n80_two_stage), ">max(N)", comparison$n80_two_stage)))
  cat(sprintf("  Bootstrap n80 median: %s  [95%% CI: %s, %s]\n",
              ifelse(is.na(comparison$n80_boot_median), ">max(N)", round(comparison$n80_boot_median)),
              ifelse(is.na(comparison$n80_boot_lo), "NA", comparison$n80_boot_lo),
              ifelse(is.na(comparison$n80_boot_hi), "NA", comparison$n80_boot_hi)))
  cat(sprintf("  Mean CI width: %.0f pp\n", 100 * mean(ci_widths, na.rm = TRUE)))

  plotDesignComparison(comparison, target_power = 0.80, panels = "A",
                       output_file = paste0(out_prefix, ".pdf"))
  cat(sprintf("  Figure: %s.pdf\n", out_prefix))

  list(label = label, color = color, n_pilot = ncol(pilot_data),
       design_type = design_type, B = B_val, N_grid = N_grid,
       comparison = comparison, ci_widths = ci_widths,
       ts_result = ts_result, boot_result = boot_result)
}

# =======================================================================
# SECTION 2: MOUSE D1D2 — D1 pilot (n=45)
# =======================================================================
cat("====================================================================\n")
cat("SECTION 2: Mouse D1 vs D2 (n=45 pilot)\n")
cat("====================================================================\n\n")

pheno_d1d2 <- read.csv(DATA_D1D2_PHENO, row.names = 1)
prep_d1d2  <- prepCircadianData(DATA_D1D2_EXPR, times = "time", input_type = "counts",
                                pheno = pheno_d1d2, sample_col = "sample")
log_d1d2   <- prep_d1d2$data

d1_samp  <- pheno_d1d2$sample[pheno_d1d2$cell == "D1"]
d2_samp  <- pheno_d1d2$sample[pheno_d1d2$cell == "D2"]
mat_d1   <- log_d1d2[, colnames(log_d1d2) %in% d1_samp, drop = FALSE]
mat_d2   <- log_d1d2[, colnames(log_d1d2) %in% d2_samp, drop = FALSE]
tod_d1   <- pheno_d1d2$time[match(colnames(mat_d1), pheno_d1d2$sample)]
tod_d2   <- pheno_d1d2$time[match(colnames(mat_d2), pheno_d1d2$sample)]
keep_d1  <- rowSums(mat_d1 > 1) >= floor(ncol(mat_d1) / 2)
mat_d1   <- mat_d1[keep_d1, , drop = FALSE]
mat_d2   <- mat_d2[keep_d1, , drop = FALSE]
set.seed(4)
g_idx_d  <- sample(nrow(mat_d1), min(NGENES, nrow(mat_d1)))
mat_d1_s <- mat_d1[g_idx_d, , drop = FALSE]
mat_d2_s <- mat_d2[g_idx_d, , drop = FALSE]

fit_d1   <- fitCosinorAll(mat_d1_s, times = tod_d1, period = 24)
fit_d2   <- fitCosinorAll(mat_d2_s, times = tod_d2, period = 24)
rhy_d1   <- !is.na(fit_d1$pvalue) & fit_d1$pvalue < RHYTHM_PVAL
rhy_d2   <- !is.na(fit_d2$pvalue) & fit_d2$pvalue < RHYTHM_PVAL
prop_DR_d <- mean(xor(rhy_d1, rhy_d2))

bio_d1d2 <- estCircadianParam(
  data = mat_d1_s, times = tod_d1, period = 24,
  prop_DR = prop_DR_d, prop_DP = 0, prop_DA = 0,
  min_rhythm_pval = RHYTHM_PVAL, verbose = FALSE
)
cat(sprintf("D1D2: prop_DR=%.1f%%  r_med=%.2f\n",
            100 * prop_DR_d,
            median(as.numeric(fit_d1$A[rhy_d1]) /
                   as.numeric(fit_d1$sigma[rhy_d1]), na.rm = TRUE)))

s2 <- .run_comparison(
  pilot_data = mat_d1_s, pilot_times = tod_d1, bio_diff_opts = bio_d1d2,
  N_grid = N_GRID_D1D2, B_val = B_D1D2, design_type = "active",
  label = sprintf("D1D2 D1 (n=%d)", ncol(mat_d1_s)), color = "forestgreen",
  out_prefix = file.path(out_dir, "s2_d1d2_comparison")
)

cat(sprintf("\nDone. Output: %s/\n", out_dir))
