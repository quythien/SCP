#' =======================================================================
#' 07f_b_bvsm_baboon.R — B vs m Tradeoff Under Waveform Distortion: Baboon
#' =======================================================================
#'
#' Sweeps B in {3,4,6,12} at fixed N in {12,24,48,72} across alpha2 grid.
#' All (B, N, alpha2) combinations run in parallel via mclapply.
#'
#' Pilot: Baboon LUN (n=12, B=12, m=1 — no replication in original study).
#' Output: output/07_bvsm_<RUN_TAG>/s2_baboon_bvsm.rds + .pdf

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

DATA_BABOON   <- "data/CAMO_PRC_hmb.RData"
PILOT_N       <- 12L
DATASET_LABEL <- "Baboon LUN vs CER (r~1.7)"
COL           <- "darkorange"

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
cat("Loading Baboon LUN vs CER...\n")
load(DATA_BABOON)
bab_expr <- baboon_withTOD$baboon
bab_tod  <- baboon_withTOD$tod

prep_lun  <- prepCircadianData(bab_expr[["LUN"]], times = bab_tod[["LUN"]], input_type = "cpm")
prep_cerb <- prepCircadianData(bab_expr[["CER"]], times = bab_tod[["CER"]], input_type = "cpm")

fit_lun  <- fitCosinorAll(prep_lun$data,  times = prep_lun$times,  period = 24)
fit_cerb <- fitCosinorAll(prep_cerb$data, times = prep_cerb$times, period = 24)
rhy_lun  <- !is.na(fit_lun$pvalue)  & fit_lun$pvalue  < RHYTHM_PVAL
rhy_cerb <- !is.na(fit_cerb$pvalue) & fit_cerb$pvalue < RHYTHM_PVAL

bio_opts <- estCircadianParam(
  data = prep_lun$data, times = prep_lun$times, period = 24,
  prop_DR = mean(xor(rhy_lun, rhy_cerb)), prop_DP = 0, prop_DA = 0,
  min_rhythm_pval = RHYTHM_PVAL, verbose = FALSE
)
cat(sprintf("Baboon: r_med=%.2f, prop_DR=%.1f%%\n",
            median(as.numeric(fit_lun$A[rhy_lun]) / as.numeric(fit_lun$sigma[rhy_lun]), na.rm = TRUE),
            100 * mean(xor(rhy_lun, rhy_cerb))))

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
saveRDS(out, file.path(out_dir, "s2_baboon_bvsm.rds"))
cat(sprintf("Saved: %s\n", file.path(out_dir, "s2_baboon_bvsm.rds")))

# =======================================================================
# Diagnostic plot
# =======================================================================
pdf(file.path(out_dir, "s2_baboon_bvsm.pdf"),
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
cat(sprintf("Figure: %s\n", file.path(out_dir, "s2_baboon_bvsm.pdf")))
cat(sprintf("\nDone. Output: %s/\n", out_dir))
