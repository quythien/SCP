#' =======================================================================
#' 07e_b_fourier_extreme_baboon.R — Extreme Harmonics: Baboon LUN vs CER
#' =======================================================================
#'
#' Split from 07e_fourier_extreme.R for parallel execution.
#' Runs Section 2 (Baboon LUN vs CER) only.
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

N_GRID_BABOON <- if (SMOKE) c(24L, 36L, 48L) else c(12L, 24L, 36L, 48L, 72L, 96L)
B_VALS_BABOON <- c(4L, 12L)
REF_N_BABOON  <- 48L
DATA_BABOON   <- "data/CAMO_PRC_hmb.RData"

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
# SECTION 2: BABOON LUN vs CER
# =======================================================================
cat("====================================================================\n")
cat("SECTION 2: Baboon LUN vs CER — extreme harmonics\n")
cat("====================================================================\n\n")

load(DATA_BABOON)
bab_expr  <- baboon_withTOD$baboon; bab_tod <- baboon_withTOD$tod
prep_lun  <- prepCircadianData(bab_expr[["LUN"]], times = bab_tod[["LUN"]], input_type = "cpm")
prep_cerb <- prepCircadianData(bab_expr[["CER"]], times = bab_tod[["CER"]], input_type = "cpm")
mat_lun   <- prep_lun$data[rowSums(prep_lun$data > 0) >= 6, , drop = FALSE]
mat_cerb  <- prep_cerb$data[rowSums(prep_cerb$data > 0) >= 6, , drop = FALSE]
tod_lun   <- prep_lun$times; tod_cerb <- prep_cerb$times
common_b  <- intersect(rownames(mat_lun), rownames(mat_cerb))
set.seed(2); g_idx_b <- sample(common_b, min(NGENES, length(common_b)))
mat_lun_s  <- mat_lun[g_idx_b,  , drop = FALSE]
mat_cerb_s <- mat_cerb[g_idx_b, , drop = FALSE]
fit_lun  <- fitCosinorAll(mat_lun_s,  times = tod_lun,  period = 24)
fit_cerb <- fitCosinorAll(mat_cerb_s, times = tod_cerb, period = 24)
rhy_lun  <- !is.na(fit_lun$pvalue)  & fit_lun$pvalue  < RHYTHM_PVAL
rhy_cerb <- !is.na(fit_cerb$pvalue) & fit_cerb$pvalue < RHYTHM_PVAL

bio_bab <- estCircadianParam(
  data = mat_lun_s, times = tod_lun, period = 24,
  prop_DR = mean(xor(rhy_lun, rhy_cerb)), prop_DP = 0,
  min_rhythm_pval = RHYTHM_PVAL, verbose = FALSE
)
cat(sprintf("Baboon: r_med=%.2f, prop_DR=%.1f%%\n",
            median(as.numeric(fit_lun$A[rhy_lun]) / as.numeric(fit_lun$sigma[rhy_lun]), na.rm = TRUE),
            100 * mean(xor(rhy_lun, rhy_cerb))))

s2_lo <- .run_fourier_one_B(bio_bab, N_GRID_BABOON, B_VALS_BABOON[1], label = "sparse B=4")
s2_hi <- .run_fourier_one_B(bio_bab, N_GRID_BABOON, B_VALS_BABOON[2], label = "dense B=12")

.plot_extreme(s2_lo, s2_hi, B_VALS_BABOON[1], B_VALS_BABOON[2],
              "Baboon (r~1.7)", REF_N_BABOON, "darkorange",
              file.path(out_dir, "s2_baboon_extreme.pdf"))

saveRDS(list(lo = s2_lo, hi = s2_hi, B_low = B_VALS_BABOON[1], B_high = B_VALS_BABOON[2],
             ref_n = REF_N_BABOON, label = "Baboon (r~1.7)", col = "darkorange"),
        file.path(out_dir, "s2_baboon_extreme.rds"))

a2_vals <- sort(unique(HARM_GRID$alpha2))
for (nm in c("s2_lo", "s2_hi")) {
  res <- get(nm); B <- if (nm == "s2_lo") B_VALS_BABOON[1] else B_VALS_BABOON[2]
  if (!is.null(res)) {
    ref_idx <- which.min(abs(res$sample_sizes - REF_N_BABOON))
    p0 <- res$power_mean[a2_vals == 0,        ref_idx]
    p1 <- res$power_mean[which.max(a2_vals),  ref_idx]
    cat(sprintf("Baboon B=%d (ref_n=%d): power %.0f%% (a2=0) -> %.0f%% (a2=1.0)  delta=%.0f pp\n",
                B, REF_N_BABOON, 100 * p0, 100 * p1, 100 * (p0 - p1)))
  }
}

cat(sprintf("\nDone. Output: %s/\n", out_dir))
