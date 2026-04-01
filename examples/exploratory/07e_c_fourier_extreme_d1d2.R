#' =======================================================================
#' 07e_c_fourier_extreme_d1d2.R — Extreme Harmonics: Mouse D1 vs D2
#' =======================================================================
#'
#' Split from 07e_fourier_extreme.R for parallel execution.
#' Runs Section 3 (Mouse D1D2) only.
#' See 07e_a_fourier_extreme_mouse.R header for parallel launch instructions.
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

N_GRID_D1D2   <- if (SMOKE) c(24L, 48L, 72L) else c(24L, 48L, 72L, 96L, 120L)
B_VALS_D1D2   <- c(3L, 6L)
REF_N_D1D2    <- 72L
DATA_D1D2_PHENO <- "data/mouse_clinicalinfo_03082021_rmOutliers.csv"
DATA_D1D2_EXPR  <- "data/mouse_D1D2_logCPMfiltered_counts.csv"

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
# SECTION 3: MOUSE D1 vs D2
# =======================================================================
cat("====================================================================\n")
cat("SECTION 3: Mouse D1 vs D2 — extreme harmonics\n")
cat("====================================================================\n\n")

pheno_d1d2 <- read.csv(DATA_D1D2_PHENO, row.names = 1)
prep_d1d2  <- prepCircadianData(DATA_D1D2_EXPR, times = "time", input_type = "counts",
                                pheno = pheno_d1d2, sample_col = "sample")
log_d1d2   <- prep_d1d2$data
d1_samp <- pheno_d1d2$sample[pheno_d1d2$cell == "D1"]
mat_d1  <- log_d1d2[, colnames(log_d1d2) %in% d1_samp, drop = FALSE]
mat_d2  <- log_d1d2[, colnames(log_d1d2) %in% pheno_d1d2$sample[pheno_d1d2$cell == "D2"], drop = FALSE]
tod_d1  <- pheno_d1d2$time[match(colnames(mat_d1), pheno_d1d2$sample)]
keep_d1 <- rowSums(mat_d1 > 1) >= floor(ncol(mat_d1) / 2)
mat_d1  <- mat_d1[keep_d1, , drop = FALSE]; mat_d2 <- mat_d2[keep_d1, , drop = FALSE]
set.seed(4); g_idx_d <- sample(nrow(mat_d1), min(NGENES, nrow(mat_d1)))
mat_d1_s <- mat_d1[g_idx_d, , drop = FALSE]; mat_d2_s <- mat_d2[g_idx_d, , drop = FALSE]
fit_d1   <- fitCosinorAll(mat_d1_s, times = tod_d1, period = 24)
rhy_d1   <- !is.na(fit_d1$pvalue) & fit_d1$pvalue < RHYTHM_PVAL
tod_d2   <- pheno_d1d2$time[match(colnames(mat_d2_s), pheno_d1d2$sample)]
fit_d2   <- fitCosinorAll(mat_d2_s, times = tod_d2, period = 24)
rhy_d2   <- !is.na(fit_d2$pvalue) & fit_d2$pvalue < RHYTHM_PVAL

bio_d1d2 <- estCircadianParam(
  data = mat_d1_s, times = tod_d1, period = 24,
  prop_DR = mean(xor(rhy_d1, rhy_d2)), prop_DP = 0, prop_DA = 0,
  min_rhythm_pval = RHYTHM_PVAL, verbose = FALSE
)
cat(sprintf("D1D2: r_med=%.2f, prop_DR=%.1f%%\n",
            median(as.numeric(fit_d1$A[rhy_d1]) / as.numeric(fit_d1$sigma[rhy_d1]), na.rm = TRUE),
            100 * mean(xor(rhy_d1, rhy_d2))))

s3_lo <- .run_fourier_one_B(bio_d1d2, N_GRID_D1D2, B_VALS_D1D2[1], label = "sparse B=3")
s3_hi <- .run_fourier_one_B(bio_d1d2, N_GRID_D1D2, B_VALS_D1D2[2], label = "dense B=6")

.plot_extreme(s3_lo, s3_hi, B_VALS_D1D2[1], B_VALS_D1D2[2],
              "D1D2 (r~0.66)", REF_N_D1D2, "forestgreen",
              file.path(out_dir, "s3_d1d2_extreme.pdf"))

saveRDS(list(lo = s3_lo, hi = s3_hi, B_low = B_VALS_D1D2[1], B_high = B_VALS_D1D2[2],
             ref_n = REF_N_D1D2, label = "D1D2 (r~0.66)", col = "forestgreen"),
        file.path(out_dir, "s3_d1d2_extreme.rds"))

a2_vals <- sort(unique(HARM_GRID$alpha2))
for (nm in c("s3_lo", "s3_hi")) {
  res <- get(nm); B <- if (nm == "s3_lo") B_VALS_D1D2[1] else B_VALS_D1D2[2]
  if (!is.null(res)) {
    ref_idx <- which.min(abs(res$sample_sizes - REF_N_D1D2))
    p0 <- res$power_mean[a2_vals == 0,        ref_idx]
    p1 <- res$power_mean[which.max(a2_vals),  ref_idx]
    cat(sprintf("D1D2 B=%d (ref_n=%d): power %.0f%% (a2=0) -> %.0f%% (a2=1.0)  delta=%.0f pp\n",
                B, REF_N_D1D2, 100 * p0, 100 * p1, 100 * (p0 - p1)))
  }
}

cat(sprintf("\nDone. Output: %s/\n", out_dir))
