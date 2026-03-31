#' =======================================================================
#' 07a_fourier_mouse.R — Fourier Robustness: Mouse GSE54651
#' =======================================================================
#'
#' Split from 07_fourier_robustness.R for parallel execution.
#' Runs Section 1 (Mouse GSE54651 LIV vs CER) only.
#' Saves RDS to shared output dir so 07d_fourier_summary.R can aggregate.
#'
#' Run all three in parallel, then run 07d for the summary figure:
#'   screen -S fourier_mouse  -dm bash -c "POWERSIM_ROOT=$ROOT Rscript examples/exploratory/07a_fourier_mouse.R  > output/07a.log 2>&1"
#'   screen -S fourier_baboon -dm bash -c "POWERSIM_ROOT=$ROOT Rscript examples/exploratory/07b_fourier_baboon.R > output/07b.log 2>&1"
#'   screen -S fourier_d1d2   -dm bash -c "POWERSIM_ROOT=$ROOT Rscript examples/exploratory/07c_fourier_d1d2.R   > output/07c.log 2>&1"
#'   # after all three finish:
#'   Rscript examples/exploratory/07d_fourier_summary.R
#'
#' SHARED OUTPUT DIR: output/07_fourier_robustness_<RUN_TAG>/
#'   Set RUN_TAG env var to coordinate across a/b/c/d (default: today's date).

# -----------------------------------------------------------------------
# 0. Path + settings
# -----------------------------------------------------------------------
POWERSIM_ROOT <- {
  env <- Sys.getenv("POWERSIM_ROOT", unset = "")
  if (nzchar(env)) env else
    "/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim"
}
setwd(POWERSIM_ROOT)

SMOKE       <- nzchar(Sys.getenv("POWERSIM_SMOKE"))
NGENES      <- if (SMOKE) 300L else 3000L
NSIMS       <- if (SMOKE) 5L  else 20L
RHYTHM_PVAL <- 0.05

HARM_GRID <- if (SMOKE) {
  expand.grid(alpha2 = c(0, 0.5), alpha3 = c(0, 0.5))
} else {
  expand.grid(alpha2 = c(0, 0.25, 0.5), alpha3 = c(0, 0.25, 0.5))
}

N_GRID_MOUSE <- if (SMOKE) c(16L, 24L, 32L) else c(8L, 16L, 24L, 48L, 72L)
B_VALS_MOUSE <- c(4L, 8L)
REF_N_MOUSE  <- 24L
DATA_MOUSE   <- "data/mice_GSE54651_CPM.RData"

RUN_TAG <- Sys.getenv("RUN_TAG", format(Sys.Date(), "%Y%m%d"))
out_dir <- file.path("output", paste0("07_fourier_robustness_", RUN_TAG))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("Mode: %s | NGENES=%d | NSIMS=%d | harm_combos=%d\n",
            if (SMOKE) "SMOKE" else "PRODUCTION", NGENES, NSIMS, nrow(HARM_GRID)))
cat(sprintf("Output -> %s/\n\n", out_dir))

source_dir <- file.path(getwd(), "code")
old_wd <- setwd(source_dir); source("setup.R"); setwd(old_wd)

opts_analysis <- CircadianAnalysisOptions(
  alpha = 0.05, p.adjust.method = "BH",
  fdr_thresholds = c(0.05), reference_n = 48L
)

# -----------------------------------------------------------------------
# Helper: one B value
# -----------------------------------------------------------------------
.run_fourier_one_B <- function(bio_opts, N_grid, B, pilot_tod_full,
                                design_type = "active", label = "") {
  cts_B <- seq(0, 24, length.out = B + 1)[seq_len(B)]
  design_opts <- CircadianDesignOptions(
    sample_sizes = N_grid, nsims = NSIMS, design = design_type,
    cts = if (design_type == "active") cts_B else pilot_tod_full,
    test_types = "DR"
  )
  cat(sprintf("    B=%d: %s\n", B, label))
  tryCatch(
    runFourierDeviationPower(bio_opts, design_opts, opts_analysis,
                             harmonic_grid = HARM_GRID, test_type = "DR",
                             verbose = FALSE),
    error = function(e) { warning(sprintf("B=%d failed: %s", B, e$message)); NULL }
  )
}

# -----------------------------------------------------------------------
# Helper: plot
# -----------------------------------------------------------------------
.plot_dataset_fourier <- function(res_low, res_high, B_low, B_high,
                                   dataset_label, ref_n, output_file) {
  pdf(output_file, width = 6, height = 5)
  par(mfrow = c(1, 1), mar = c(4.5, 4.5, 3.5, 1.5))
  if (!is.null(res_low) && !is.null(res_high)) {
    hg <- res_low$harmonic_grid
    a3_0 <- hg$alpha3 == 0
    a2_vals <- sort(unique(hg$alpha2[a3_0]))
    ref_idx_lo <- which.min(abs(res_low$sample_sizes  - ref_n))
    ref_idx_hi <- which.min(abs(res_high$sample_sizes - ref_n))
    pow_lo <- res_low$power_mean[a3_0,  ref_idx_lo]
    pow_hi <- res_high$power_mean[a3_0, ref_idx_hi]
    se_lo  <- res_low$power_se[a3_0,   ref_idx_lo]
    se_hi  <- res_high$power_se[a3_0,  ref_idx_hi]
    plot(a2_vals, pow_lo, type = "b", pch = 16, lwd = 2, col = "steelblue",
         ylim = c(0, 1), xlab = expression(alpha[2] ~ "(2nd harmonic)"),
         ylab = "Power (FDR 5%)", las = 1,
         main = sprintf("%s\nB vs harmonics at N=%d (a3=0)", dataset_label, ref_n))
    arrows(a2_vals, pow_lo - se_lo, a2_vals, pow_lo + se_lo,
           length = 0.05, angle = 90, code = 3, col = "steelblue", lwd = 1.2)
    lines(a2_vals, pow_hi, type = "b", pch = 17, lwd = 2, col = "darkorange", lty = 2)
    arrows(a2_vals, pow_hi - se_hi, a2_vals, pow_hi + se_hi,
           length = 0.05, angle = 90, code = 3, col = "darkorange", lwd = 1.2)
    abline(h = 0.80, lty = 3, col = "gray50")
    legend("bottomleft", legend = sprintf("B=%d", c(B_low, B_high)),
           col = c("steelblue", "darkorange"), lwd = 2,
           pch = c(16, 17), lty = c(1, 2), bty = "n")
  } else {
    plot.new(); text(0.5, 0.5, "Comparison unavailable", cex = 1.2)
  }
  dev.off()
  cat(sprintf("  Figure: %s\n", output_file))
}

# =======================================================================
# SECTION 1: MOUSE GSE54651 — LIV pilot
# =======================================================================
cat("====================================================================\n")
cat("SECTION 1: Mouse GSE54651 — LIV vs CER Fourier robustness\n")
cat("====================================================================\n\n")

dat_mouse <- readRDS(DATA_MOUSE)
ct_liv    <- dat_mouse$tod[["LIV"]]
prep_liv  <- prepCircadianData(dat_mouse$count_clean[["LIV"]], times = ct_liv, input_type = "log2")
mat_liv   <- prep_liv$data
mat_liv   <- mat_liv[rowSums(mat_liv > 0) >= 4, , drop = FALSE]
set.seed(1)
g_idx_m   <- sample(nrow(mat_liv), min(NGENES, nrow(mat_liv)))
mat_liv_s <- mat_liv[g_idx_m, , drop = FALSE]

prep_cer <- prepCircadianData(dat_mouse$count_clean[["CER"]], times = dat_mouse$tod[["CER"]], input_type = "log2")
mat_cer  <- prep_cer$data
mat_cer  <- mat_cer[rowSums(mat_cer > 0) >= 4, , drop = FALSE]
common_m <- intersect(rownames(mat_liv_s), rownames(mat_cer))

fit_liv  <- fitCosinorAll(mat_liv_s[common_m, ], times = ct_liv, period = 24)
fit_cer  <- fitCosinorAll(mat_cer[common_m, ], times = prep_cer$times, period = 24)
rhy_liv  <- !is.na(fit_liv$pvalue) & fit_liv$pvalue < RHYTHM_PVAL
rhy_cer  <- !is.na(fit_cer$pvalue) & fit_cer$pvalue < RHYTHM_PVAL
prop_DR_m <- mean(xor(rhy_liv, rhy_cer))

bio_mouse <- estCircadianParam(
  data = mat_liv_s, times = ct_liv, period = 24,
  prop_DR = prop_DR_m, prop_DP = 0, prop_DA = 0,
  min_rhythm_pval = RHYTHM_PVAL, verbose = FALSE
)
cat(sprintf("Mouse pilot: n=%d, prop_DR=%.1f%%, r_med=%.2f\n",
            ncol(mat_liv_s), 100 * prop_DR_m,
            median(as.numeric(fit_liv$A[rhy_liv]) /
                   as.numeric(fit_liv$sigma[rhy_liv]), na.rm = TRUE)))

cat("Running Fourier analysis (B=4 and B=8)...\n")
s1_B4 <- .run_fourier_one_B(bio_mouse, N_GRID_MOUSE, B_VALS_MOUSE[1], ct_liv, label = "sparse (every 6h)")
s1_B8 <- .run_fourier_one_B(bio_mouse, N_GRID_MOUSE, B_VALS_MOUSE[2], ct_liv, label = "full 8 CT")

.plot_dataset_fourier(s1_B4, s1_B8, B_VALS_MOUSE[1], B_VALS_MOUSE[2],
                       "GSE54651 LIV vs CER", REF_N_MOUSE,
                       file.path(out_dir, "s1_mouse_fourier.pdf"))

saveRDS(list(B4 = s1_B4, B8 = s1_B8,
             label = "GSE54651 (r~2.9)", col = "steelblue",
             B_low = B_VALS_MOUSE[1], B_high = B_VALS_MOUSE[2], ref_n = REF_N_MOUSE),
        file.path(out_dir, "s1_mouse_fourier.rds"))

cat(sprintf("\nDone. Output: %s/\n", out_dir))
