#' =======================================================================
#' 07e_a_fourier_extreme_mouse.R — Extreme Harmonics: Mouse GSE54651
#' =======================================================================
#'
#' Split from 07e_fourier_extreme.R for parallel execution.
#' Runs Section 1 (Mouse GSE54651 LIV vs CER) only.
#' alpha2 grid: 0, 0.25, 0.50, 0.75, 1.00 (alpha3=0 throughout).
#'
#' Run all three in parallel, then run 07e_d for the summary figure:
#'   screen -S 7ea -dm bash -c "POWERSIM_ROOT=$ROOT RUN_TAG=$TAG Rscript examples/exploratory/07e_a_fourier_extreme_mouse.R  > output/07e_a.log 2>&1"
#'   screen -S 7eb -dm bash -c "POWERSIM_ROOT=$ROOT RUN_TAG=$TAG Rscript examples/exploratory/07e_b_fourier_extreme_baboon.R > output/07e_b.log 2>&1"
#'   screen -S 7ec -dm bash -c "POWERSIM_ROOT=$ROOT RUN_TAG=$TAG Rscript examples/exploratory/07e_c_fourier_extreme_d1d2.R   > output/07e_c.log 2>&1"
#'   # after all three finish:
#'   POWERSIM_ROOT=$ROOT RUN_TAG=$TAG Rscript examples/exploratory/07e_d_fourier_extreme_summary.R
#'
#' SHARED OUTPUT DIR: output/07_fourier_extreme_<RUN_TAG>/

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
  expand.grid(alpha2 = c(0, 0.5, 1.0), alpha3 = 0)
} else {
  expand.grid(alpha2 = c(0, 0.25, 0.5, 0.75, 1.0), alpha3 = 0)
}

N_GRID_MOUSE  <- if (SMOKE) c(16L, 24L, 32L) else c(8L, 16L, 24L, 48L, 72L)
B_VALS_MOUSE  <- c(4L, 8L)
REF_N_MOUSE   <- 24L
DATA_MOUSE    <- "data/mice_GSE54651_CPM.RData"

RUN_TAG <- Sys.getenv("RUN_TAG", format(Sys.Date(), "%Y%m%d"))
out_dir <- file.path("output", paste0("07_fourier_extreme_", RUN_TAG))
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

.run_fourier_one_B <- function(bio_opts, N_grid, B, design_type = "active", label = "") {
  cts_B <- seq(0, 24, length.out = B + 1)[seq_len(B)]
  design_opts <- CircadianDesignOptions(
    sample_sizes = N_grid, nsims = NSIMS, design = design_type,
    cts = cts_B, test_types = "DR"
  )
  cat(sprintf("    B=%d: %s\n", B, label))
  tryCatch(
    runFourierDeviationPower(bio_opts, design_opts, opts_analysis,
                             harmonic_grid = HARM_GRID, test_type = "DR",
                             verbose = FALSE),
    error = function(e) { warning(sprintf("B=%d failed: %s", B, e$message)); NULL }
  )
}

.plot_extreme <- function(res_low, res_high, B_low, B_high,
                           dataset_label, ref_n, col, output_file) {
  pdf(output_file, width = 7, height = 5)
  par(mar = c(4.5, 4.5, 3.5, 1.5))
  a2_vals <- sort(unique(HARM_GRID$alpha2))

  if (!is.null(res_low) && !is.null(res_high)) {
    ref_idx_lo <- which.min(abs(res_low$sample_sizes  - ref_n))
    ref_idx_hi <- which.min(abs(res_high$sample_sizes - ref_n))
    pow_lo <- res_low$power_mean[,  ref_idx_lo]
    pow_hi <- res_high$power_mean[, ref_idx_hi]
    se_lo  <- res_low$power_se[,   ref_idx_lo]
    se_hi  <- res_high$power_se[,  ref_idx_hi]
    ylim   <- c(0, min(1, max(c(pow_lo, pow_hi) + c(se_lo, se_hi), na.rm = TRUE) * 1.1))

    plot(a2_vals, pow_lo, type = "b", pch = 16, lwd = 2.5, col = col,
         ylim = ylim, xlab = expression(alpha[2] ~ "(2nd harmonic relative amplitude)"),
         ylab = "Power (FDR 5%)", las = 1,
         main = sprintf("%s — extreme harmonics\nN=%d; solid=B=%d, dashed=B=%d",
                        dataset_label, ref_n, B_low, B_high))
    polygon(c(a2_vals, rev(a2_vals)), c(pow_lo - se_lo, rev(pow_lo + se_lo)),
            col = adjustcolor(col, 0.12), border = NA)
    lines(a2_vals, pow_hi, type = "b", pch = 17, lwd = 2.5, col = col, lty = 2)
    polygon(c(a2_vals, rev(a2_vals)), c(pow_hi - se_hi, rev(pow_hi + se_hi)),
            col = adjustcolor(col, 0.08), border = NA)
    abline(h = 0.80, lty = 3, col = "gray50")
    abline(v = 0.5,  lty = 2, col = "gray70")
    text(0.51, ylim[2] * 0.95, "standard\nmax", col = "gray50", cex = 0.7, adj = 0)
    legend("topright", legend = sprintf("B=%d", c(B_low, B_high)),
           col = col, lwd = 2.5, pch = c(16, 17), lty = c(1, 2), bty = "n")
  } else {
    plot.new(); text(0.5, 0.5, "Results unavailable", cex = 1.2)
  }
  dev.off()
  cat(sprintf("  Figure: %s\n", output_file))
}

# =======================================================================
# SECTION 1: MOUSE GSE54651
# =======================================================================
cat("====================================================================\n")
cat("SECTION 1: Mouse GSE54651 — extreme harmonics\n")
cat("====================================================================\n\n")

dat_mouse <- readRDS(DATA_MOUSE)
ct_liv    <- dat_mouse$tod[["LIV"]]
prep_liv  <- prepCircadianData(dat_mouse$count_clean[["LIV"]], times = ct_liv, input_type = "log2")
mat_liv   <- prep_liv$data[rowSums(prep_liv$data > 0) >= 4, , drop = FALSE]
set.seed(1)
mat_liv_s <- mat_liv[sample(nrow(mat_liv), min(NGENES, nrow(mat_liv))), , drop = FALSE]

prep_cer <- prepCircadianData(dat_mouse$count_clean[["CER"]], times = dat_mouse$tod[["CER"]], input_type = "log2")
mat_cer  <- prep_cer$data[rowSums(prep_cer$data > 0) >= 4, , drop = FALSE]
common_m <- intersect(rownames(mat_liv_s), rownames(mat_cer))
fit_liv  <- fitCosinorAll(mat_liv_s[common_m, ], times = ct_liv,         period = 24)
fit_cer  <- fitCosinorAll(mat_cer[common_m, ],   times = prep_cer$times, period = 24)
rhy_liv  <- !is.na(fit_liv$pvalue) & fit_liv$pvalue < RHYTHM_PVAL
rhy_cer  <- !is.na(fit_cer$pvalue) & fit_cer$pvalue < RHYTHM_PVAL

bio_mouse <- estCircadianParam(
  data = mat_liv_s, times = ct_liv, period = 24,
  prop_DR = mean(xor(rhy_liv, rhy_cer)), prop_DP = 0,
  min_rhythm_pval = RHYTHM_PVAL, verbose = FALSE
)
cat(sprintf("Mouse: r_med=%.2f, prop_DR=%.1f%%\n",
            median(as.numeric(fit_liv$A[rhy_liv]) / as.numeric(fit_liv$sigma[rhy_liv]), na.rm = TRUE),
            100 * mean(xor(rhy_liv, rhy_cer))))

s1_lo <- .run_fourier_one_B(bio_mouse, N_GRID_MOUSE, B_VALS_MOUSE[1], label = "sparse B=4")
s1_hi <- .run_fourier_one_B(bio_mouse, N_GRID_MOUSE, B_VALS_MOUSE[2], label = "dense B=8")

.plot_extreme(s1_lo, s1_hi, B_VALS_MOUSE[1], B_VALS_MOUSE[2],
              "GSE54651 (r~2.9)", REF_N_MOUSE, "steelblue",
              file.path(out_dir, "s1_mouse_extreme.pdf"))

saveRDS(list(lo = s1_lo, hi = s1_hi, B_low = B_VALS_MOUSE[1], B_high = B_VALS_MOUSE[2],
             ref_n = REF_N_MOUSE, label = "GSE54651 (r~2.9)", col = "steelblue"),
        file.path(out_dir, "s1_mouse_extreme.rds"))

a2_vals <- sort(unique(HARM_GRID$alpha2))
for (nm in c("s1_lo", "s1_hi")) {
  res <- get(nm); B <- if (nm == "s1_lo") B_VALS_MOUSE[1] else B_VALS_MOUSE[2]
  if (!is.null(res)) {
    ref_idx <- which.min(abs(res$sample_sizes - REF_N_MOUSE))
    p0 <- res$power_mean[a2_vals == 0,        ref_idx]
    p1 <- res$power_mean[which.max(a2_vals),  ref_idx]
    cat(sprintf("GSE54651 B=%d (ref_n=%d): power %.0f%% (a2=0) -> %.0f%% (a2=1.0)  delta=%.0f pp\n",
                B, REF_N_MOUSE, 100 * p0, 100 * p1, 100 * (p0 - p1)))
  }
}

cat(sprintf("\nDone. Output: %s/\n", out_dir))
