#' =======================================================================
#' 08a_bootstrap_baboon.R — Two-Stage vs Bootstrap: Baboon LUN (n=12)
#' =======================================================================
#' Split from 08_two_stage_vs_bootstrap_realdata.R for parallel execution.
#' Runs Section 1 (Baboon, smallest pilot → widest CI) only.
#'
#' Run all three in parallel, then run 08d for the summary figure:
#'   ROOT=/path/to/PowerSim
#'   TAG=$(date +%Y%m%d)
#'   screen -S boot_baboon -dm bash -c "POWERSIM_ROOT=$ROOT RUN_TAG=$TAG Rscript examples/publication/08a_bootstrap_baboon.R > output/08a.log 2>&1"
#'   screen -S boot_d1d2   -dm bash -c "POWERSIM_ROOT=$ROOT RUN_TAG=$TAG Rscript examples/publication/08b_bootstrap_d1d2.R   > output/08b.log 2>&1"
#'   screen -S boot_seney  -dm bash -c "POWERSIM_ROOT=$ROOT RUN_TAG=$TAG Rscript examples/publication/08c_bootstrap_seney.R  > output/08c.log 2>&1"
#'   # after all three finish:
#'   POWERSIM_ROOT=$ROOT RUN_TAG=$TAG Rscript examples/publication/08d_bootstrap_summary.R

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
NCORES      <- if (SMOKE) 2L  else 20L
RHYTHM_PVAL <- 0.05

B_BABOON      <- 12L
N_GRID_BABOON <- if (SMOKE) c(12L, 24L, 36L) else c(12L, 24L, 36L, 48L, 60L, 72L, 96L)
DATA_BABOON   <- "data/CAMO_PRC_hmb.RData"

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

# -----------------------------------------------------------------------
# Helper
# -----------------------------------------------------------------------
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
                           bio_diff.opts = bio_diff_opts, mode = "differential",
                           verbose = FALSE, mc.cores = NCORES),
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
# SECTION 1: BABOON — LUN pilot (n=12)
# =======================================================================
cat("====================================================================\n")
cat("SECTION 1: Baboon LUN vs CER (n=12 pilot)\n")
cat("====================================================================\n\n")

load(DATA_BABOON)
bab_expr <- baboon_withTOD$baboon
bab_tod  <- baboon_withTOD$tod

prep_lun  <- prepCircadianData(bab_expr[["LUN"]], times = bab_tod[["LUN"]], input_type = "cpm")
prep_cerb <- prepCircadianData(bab_expr[["CER"]], times = bab_tod[["CER"]], input_type = "cpm")
mat_lun   <- prep_lun$data[rowSums(prep_lun$data > 0) >= 6, , drop = FALSE]
mat_cerb  <- prep_cerb$data[rowSums(prep_cerb$data > 0) >= 6, , drop = FALSE]
tod_lun   <- prep_lun$times
common_b  <- intersect(rownames(mat_lun), rownames(mat_cerb))
set.seed(2)
g_idx_b   <- sample(common_b, min(NGENES, length(common_b)))
mat_lun_s  <- mat_lun[g_idx_b, , drop = FALSE]
mat_cerb_s <- mat_cerb[g_idx_b, , drop = FALSE]

fit_lun  <- fitCosinorAll(mat_lun_s,  times = tod_lun,              period = 24)
fit_cerb <- fitCosinorAll(mat_cerb_s, times = prep_cerb$times,      period = 24)
rhy_lun  <- !is.na(fit_lun$pvalue)  & fit_lun$pvalue  < RHYTHM_PVAL
rhy_cerb <- !is.na(fit_cerb$pvalue) & fit_cerb$pvalue < RHYTHM_PVAL
prop_DR_b <- mean(xor(rhy_lun, rhy_cerb))

bio_bab <- estCircadianParam(
  data = mat_lun_s, times = tod_lun, period = 24,
  prop_DR = prop_DR_b, prop_DP = 0,
  min_rhythm_pval = RHYTHM_PVAL, verbose = FALSE
)
cat(sprintf("Baboon: prop_DR=%.1f%%  r_med=%.2f\n",
            100 * prop_DR_b,
            median(as.numeric(fit_lun$A[rhy_lun]) /
                   as.numeric(fit_lun$sigma[rhy_lun]), na.rm = TRUE)))

s1 <- .run_comparison(
  pilot_data = mat_lun_s, pilot_times = tod_lun, bio_diff_opts = bio_bab,
  N_grid = N_GRID_BABOON, B_val = B_BABOON, design_type = "active",
  label = sprintf("Baboon LUN (n=%d)", ncol(mat_lun_s)), color = "darkorange",
  out_prefix = file.path(out_dir, "s1_baboon_comparison")
)

cat(sprintf("\nDone. Output: %s/\n", out_dir))
