#' =======================================================================
#' 07f_c_bvsm_d1d2.R — B vs m Tradeoff Under Waveform Distortion: D1D2
#' =======================================================================
#'
#' Sweeps B in {3,4,6,12} at fixed N in {12,24,48,72} across alpha2 grid.
#' All (B, N, alpha2) combinations run in parallel via mclapply.
#'
#' Pilot: Mouse D1D2 (n=45, B=6, active design).
#' Output: output/07_bvsm_<RUN_TAG>/s3_d1d2_bvsm.rds + .pdf

POWERSIM_ROOT <- {
  env <- Sys.getenv("POWERSIM_ROOT", unset = "")
  if (nzchar(env)) env else
    "/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim"
}
setwd(POWERSIM_ROOT)

SMOKE       <- nzchar(Sys.getenv("POWERSIM_SMOKE"))
NGENES      <- if (SMOKE) 300L  else 3000L
NSIMS       <- if (SMOKE) 5L    else 20L
NCORES      <- if (SMOKE) 2L    else 20L
RHYTHM_PVAL <- 0.05

B_VALS  <- if (SMOKE) c(3L, 6L)   else c(3L, 4L, 6L, 12L)
N_GRID  <- if (SMOKE) c(12L, 24L) else c(12L, 24L, 48L, 72L)
HARM_GRID <- expand.grid(alpha2 = c(0, 0.25, 0.5, 0.75, 1.0), alpha3 = 0)

DATA_D1D2     <- "data/mouse_D1D2_logCPMfiltered_counts.csv"
PILOT_N       <- 45L
DATASET_LABEL <- "D1D2 (r~0.64)"
COL           <- "forestgreen"

RUN_TAG <- Sys.getenv("RUN_TAG", format(Sys.Date(), "%Y%m%d"))
out_dir <- file.path("output", paste0("07_bvsm_", RUN_TAG))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("Mode: %s | NGENES=%d | NSIMS=%d | NCORES=%d\n",
            if (SMOKE) "SMOKE" else "PRODUCTION", NGENES, NSIMS, NCORES))
cat(sprintf("B grid: %s\n", paste(B_VALS, collapse = ", ")))
cat(sprintf("N grid: %s\n", paste(N_GRID, collapse = ", ")))
cat(sprintf("Output -> %s/\n\n", out_dir))

source_dir <- file.path(getwd(), "code")
old_wd <- setwd(source_dir); source("setup.R"); setwd(old_wd)

opts_analysis <- CircadianAnalysisOptions(
  alpha = 0.05, p.adjust.method = "BH",
  fdr_thresholds = c(0.05), reference_n = PILOT_N
)

# =======================================================================
# Load and prep pilot data
# =======================================================================
cat("Loading Mouse D1D2...\n")
pheno_d1d2 <- read.csv("data/mouse_clinicalinfo_03082021_rmOutliers.csv", row.names = 1)
prep_d1d2  <- prepCircadianData(DATA_D1D2, times = "time", input_type = "counts",
                                pheno = pheno_d1d2, sample_col = "sample")
log_d1d2   <- prep_d1d2$data

d1_samp <- pheno_d1d2$sample[pheno_d1d2$cell == "D1"]
d2_samp <- pheno_d1d2$sample[pheno_d1d2$cell == "D2"]
mat_d1  <- log_d1d2[, colnames(log_d1d2) %in% d1_samp, drop = FALSE]
mat_d2  <- log_d1d2[, colnames(log_d1d2) %in% d2_samp, drop = FALSE]
tod_d1  <- pheno_d1d2$time[match(colnames(mat_d1), pheno_d1d2$sample)]
tod_d2  <- pheno_d1d2$time[match(colnames(mat_d2), pheno_d1d2$sample)]

keep_d1 <- rowSums(mat_d1 > 1) >= floor(ncol(mat_d1) / 2)
mat_d1  <- mat_d1[keep_d1, , drop = FALSE]
mat_d2  <- mat_d2[keep_d1, , drop = FALSE]
set.seed(4)
g_idx   <- sample(nrow(mat_d1), min(NGENES, nrow(mat_d1)))
mat_d1_s <- mat_d1[g_idx, , drop = FALSE]
mat_d2_s <- mat_d2[g_idx, , drop = FALSE]

fit_d1  <- fitCosinorAll(mat_d1_s, times = tod_d1, period = 24)
fit_d2  <- fitCosinorAll(mat_d2_s, times = tod_d2, period = 24)
rhy_d1  <- !is.na(fit_d1$pvalue) & fit_d1$pvalue < RHYTHM_PVAL
rhy_d2  <- !is.na(fit_d2$pvalue) & fit_d2$pvalue < RHYTHM_PVAL

bio_opts <- estCircadianParam(
  data = mat_d1_s, times = tod_d1, period = 24,
  prop_DR = mean(xor(rhy_d1, rhy_d2)), prop_DP = 0, prop_DA = 0,
  min_rhythm_pval = RHYTHM_PVAL, verbose = FALSE
)
cat(sprintf("D1D2: r_med=%.2f, prop_DR=%.1f%%\n",
            median(as.numeric(fit_d1$A[rhy_d1]) / as.numeric(fit_d1$sigma[rhy_d1]), na.rm = TRUE),
            100 * mean(xor(rhy_d1, rhy_d2))))

# =======================================================================
# Parallel sweep
# =======================================================================
cat(sprintf("\nBuilding combination grid: %d B x %d N x %d alpha2 = %d total\n",
            length(B_VALS), length(N_GRID), nrow(HARM_GRID),
            length(B_VALS) * length(N_GRID) * nrow(HARM_GRID)))

combo_grid <- expand.grid(
  b_idx = seq_along(B_VALS),
  n_idx = seq_along(N_GRID),
  h_idx = seq_len(nrow(HARM_GRID))
)
cat(sprintf("Launching %d combinations on %d cores...\n\n", nrow(combo_grid), NCORES))

run_one <- function(i) {
  B      <- B_VALS[combo_grid$b_idx[i]]
  N      <- N_GRID[combo_grid$n_idx[i]]
  alpha2 <- HARM_GRID$alpha2[combo_grid$h_idx[i]]
  alpha3 <- HARM_GRID$alpha3[combo_grid$h_idx[i]]

  cts_B <- seq(0, 24, length.out = B + 1)[seq_len(B)]
  cts_N <- rep(cts_B, each = N / B)

  design_opts <- CircadianDesignOptions(
    sample_sizes = N, nsims = NSIMS, design = "active",
    cts = cts_N, test_types = "DR"
  )

  tryCatch({
    sim_out <- runSimsDiff(bio_opts, design_opts, opts_analysis,
                           harmonics = c(alpha2, alpha3))
    fdr_arr <- sim_out$fdr_DR
    power_vec <- vapply(seq_len(NSIMS), function(s) {
      fdr_vec   <- fdr_arr[, 1, s]
      diff_type <- sim_out$diff_type[[s]]
      target    <- diff_type %in% c(2L, 3L)
      if (sum(target) == 0) return(NA_real_)
      sum(fdr_vec[target] <= 0.05, na.rm = TRUE) / sum(target)
    }, numeric(1))
    list(B = B, N = N, alpha2 = alpha2, alpha3 = alpha3,
         power_mean = mean(power_vec, na.rm = TRUE),
         power_se   = sd(power_vec,   na.rm = TRUE),
         power_vec  = power_vec)
  }, error = function(e) {
    message(sprintf("  FAILED B=%d N=%d a2=%.2f: %s", B, N, alpha2, e$message))
    list(B = B, N = N, alpha2 = alpha2, alpha3 = alpha3,
         power_mean = NA_real_, power_se = NA_real_, power_vec = rep(NA_real_, NSIMS))
  })
}

results_list <- parallel::mclapply(
  seq_len(nrow(combo_grid)), run_one,
  mc.cores = NCORES, mc.set.seed = TRUE
)

# =======================================================================
# Assemble
# =======================================================================
power_mean <- array(NA_real_,
  dim      = c(length(B_VALS), length(N_GRID), nrow(HARM_GRID)),
  dimnames = list(paste0("B", B_VALS), paste0("N", N_GRID),
                  sprintf("a2=%.2f", HARM_GRID$alpha2)))
power_se <- power_mean

for (i in seq_along(results_list)) {
  r     <- results_list[[i]]
  b_idx <- combo_grid$b_idx[i]
  n_idx <- combo_grid$n_idx[i]
  h_idx <- combo_grid$h_idx[i]
  power_mean[b_idx, n_idx, h_idx] <- r$power_mean
  power_se[b_idx,   n_idx, h_idx] <- r$power_se
}

cat("\nPower at alpha2=0 across B and N (Proposition 1 check):\n")
h0 <- which(HARM_GRID$alpha2 == 0)
print(round(power_mean[, , h0], 3))

out <- list(
  power_mean  = power_mean,
  power_se    = power_se,
  B_vals      = B_VALS,
  N_grid      = N_GRID,
  harm_grid   = HARM_GRID,
  pilot_n     = PILOT_N,
  label       = DATASET_LABEL,
  col         = COL,
  nsims       = NSIMS,
  ngenes      = NGENES
)
saveRDS(out, file.path(out_dir, "s3_d1d2_bvsm.rds"))
cat(sprintf("Saved: %s\n", file.path(out_dir, "s3_d1d2_bvsm.rds")))

# =======================================================================
# Diagnostic plot
# =======================================================================
pdf(file.path(out_dir, "s3_d1d2_bvsm.pdf"),
    width = 4 * length(N_GRID), height = 4.5)
par(mfrow = c(1, length(N_GRID)), mar = c(4.5, 4.5, 3, 1))
a2_vals <- HARM_GRID$alpha2
b_cols  <- colorRampPalette(c("gray70", COL))(length(B_VALS))
b_pchs  <- c(16, 17, 15, 18)[seq_along(B_VALS)]

for (ni in seq_along(N_GRID)) {
  N <- N_GRID[ni]
  plot(NA, xlim = range(a2_vals), ylim = c(0, 1), las = 1,
       xlab = expression(alpha[2]), ylab = "Power (DR, FDR 5%)",
       main = sprintf("%s\nN = %d", DATASET_LABEL, N))
  abline(h = 0.80, lty = 3, col = "gray70")
  abline(v = 0.50, lty = 3, col = "gray70")
  for (bi in seq_along(B_VALS)) {
    pm <- power_mean[bi, ni, ]
    se <- power_se[bi, ni, ]
    lines(a2_vals, pm, col = b_cols[bi], lwd = 2)
    points(a2_vals, pm, col = b_cols[bi], pch = b_pchs[bi], cex = 1.1)
    polygon(c(a2_vals, rev(a2_vals)),
            c(pm - se, rev(pm + se)),
            col = adjustcolor(b_cols[bi], 0.12), border = NA)
  }
  if (ni == 1)
    legend("topright", legend = sprintf("B=%d (m=%d)", B_VALS, N / B_VALS),
           col = b_cols, lwd = 2, pch = b_pchs, bty = "n", cex = 0.8)
}
dev.off()
cat(sprintf("Figure: %s\n", file.path(out_dir, "s3_d1d2_bvsm.pdf")))
cat(sprintf("\nDone. Output: %s/\n", out_dir))
