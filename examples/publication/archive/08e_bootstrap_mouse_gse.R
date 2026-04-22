#' =======================================================================
#' 08e_bootstrap_mouse_gse.R — Two-Stage vs Bootstrap: Mouse LIV (n=8, GSE54651)
#' =======================================================================
#' Mirrors 08a/b/c structure. Runs two-stage vs bootstrap comparison for
#' Mouse LIV vs CER (Zhang et al. 2014, GSE54651). Strong signal (r~2.9)
#' provides the high-effect-size anchor for cross-dataset comparison in 08d.
#'
#' Screen launch (parallel with 08a/b/c):
#'   screen -S boot_mouse -dm bash -c "POWERSIM_ROOT=$ROOT RUN_TAG=$TAG \
#'     Rscript examples/publication/08e_bootstrap_mouse_gse.R > output/08e.log 2>&1"

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

B_MOUSE      <- 8L
N_GRID_MOUSE <- if (SMOKE) c(16L, 24L, 48L) else c(16L, 24L, 32L, 48L, 64L, 72L, 96L)
DATA_MOUSE   <- "data/mice_GSE54651_CPM.RData"

RUN_TAG <- Sys.getenv("RUN_TAG", format(Sys.Date(), "%Y%m%d"))
out_dir <- file.path("output", paste0("08_two_stage_vs_bootstrap_", RUN_TAG))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("Mode: %s | NGENES=%d | NBOOT=%d | NSIMS=%d | NSIMS_INNER=%d\n",
            if (SMOKE) "SMOKE" else "PRODUCTION", NGENES, NBOOT, NSIMS, NSIMS_INNER))
cat(sprintf("Output -> %s/\n\n", out_dir))

source_dir <- file.path(getwd(), "code")
old_wd <- setwd(source_dir); source("setup.R"); setwd(old_wd)

opts_analysis <- CircadianAnalysisOptions(
  alpha = 0.05, p.adjust.method = "BH",
  fdr_thresholds = c(0.05), reference_n = 24L
)

# -----------------------------------------------------------------------
# Helper (same as 08a)
# -----------------------------------------------------------------------
.run_comparison <- function(pilot_data, pilot_times, bio_diff_opts,
                             N_grid, B_val, design_type, label, color,
                             out_prefix) {
  cat(sprintf("\n--- %s ---\n", label))
  cat(sprintf("    Pilot: n=%d | design=%s | B=%d\n",
              ncol(pilot_data), design_type, B_val))

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
                           bio_diff.opts = bio_diff_opts, verbose = FALSE,
                           mc.cores = NCORES),
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
# SECTION 1: MOUSE LIV vs CER — GSE54651 (n=8 per tissue, B=8)
# =======================================================================
cat("====================================================================\n")
cat("SECTION 1: Mouse LIV vs CER, GSE54651 (n=8 pilot)\n")
cat("====================================================================\n\n")

dat_mouse <- readRDS(DATA_MOUSE)
ct_liv    <- dat_mouse$tod[["LIV"]]
ct_cer    <- dat_mouse$tod[["CER"]]

prep_liv  <- prepCircadianData(dat_mouse$count_clean[["LIV"]],
                                times = ct_liv, input_type = "log2")
prep_cer  <- prepCircadianData(dat_mouse$count_clean[["CER"]],
                                times = ct_cer, input_type = "log2")
rm(dat_mouse)

keep_liv <- rowSums(prep_liv$data > 0) >= 4
keep_cer <- rowSums(prep_cer$data > 0) >= 4
common_g <- intersect(rownames(prep_liv$data)[keep_liv],
                       rownames(prep_cer$data)[keep_cer])
set.seed(2)
g_idx    <- sample(common_g, min(NGENES, length(common_g)))
mat_liv  <- prep_liv$data[g_idx, , drop = FALSE]
mat_cer  <- prep_cer$data[g_idx, , drop = FALSE]
tod_liv  <- prep_liv$times
tod_cer  <- prep_cer$times
rm(prep_liv, prep_cer)

fit_liv  <- fitCosinorAll(mat_liv, times = tod_liv, period = 24)
fit_cer  <- fitCosinorAll(mat_cer, times = tod_cer, period = 24)
rhy_liv  <- !is.na(fit_liv$pvalue) & fit_liv$pvalue < RHYTHM_PVAL
rhy_cer  <- !is.na(fit_cer$pvalue) & fit_cer$pvalue < RHYTHM_PVAL
prop_DR_m <- mean(xor(rhy_liv, rhy_cer))
r_liv_med <- median(as.numeric(fit_liv$A[rhy_liv]) /
                    as.numeric(fit_liv$sigma[rhy_liv]), na.rm = TRUE)

bio_mouse <- estCircadianParam(
  data = mat_liv, times = tod_liv, period = 24,
  prop_DR = prop_DR_m, prop_DP = 0,
  min_rhythm_pval = RHYTHM_PVAL, verbose = FALSE
)
cat(sprintf("Mouse LIV: prop_DR=%.1f%%  r_med=%.2f\n",
            100 * prop_DR_m, r_liv_med))

s1 <- .run_comparison(
  pilot_data    = mat_liv,
  pilot_times   = tod_liv,
  bio_diff_opts = bio_mouse,
  N_grid        = N_GRID_MOUSE,
  B_val         = B_MOUSE,
  design_type   = "active",
  label         = sprintf("Mouse LIV vs CER (n=%d)", ncol(mat_liv)),
  color         = "steelblue",
  out_prefix    = file.path(out_dir, "s4_mouse_gse_comparison")
)

cat(sprintf("\nDone. Output: %s/\n", out_dir))
