#!/usr/bin/env Rscript
# =============================================================================
# Fig 6 v8 - Active B vs m trade-off, USING THE PACKAGE SIMULATOR (no custom)
#   Panel A: runSingleCohortPower on K=1 paired Liver pilot, active design per B
#   Panel B: simCircadianSingleCohort2H on K=2 paired Liver pilot, active per B
# Matches Fig 1B simulator for Panel A; matches Fig 5 simulator for Panel B.
# =============================================================================
setwd("/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")
sink(stderr())
source("code/utils.R"); source("code/options.R"); source("code/simulation.R")
source("code/estimation.R"); source("code/bootstrap_sim.R"); source("code/detection.R")
source("code/runner.R")
sink()
suppressPackageStartupMessages(library(parallel))

bio_K1 <- readRDS("data/gtex_Liver_single_pilot.rds")
bio_K1$paired_sigma <- TRUE          # match Fig 1B paired Panel A
psi_K2 <- readRDS("output/two_harmonic/results/pilot_2h_GTExLiver_topK300.rds")
cat(sprintf("K=1 pilot: n_amp=%d, n_sigma=%d, prop_rhy=%.3f, paired_sigma=TRUE\n",
            length(bio_K1$amplitude), length(bio_K1$sigma_rhythmic),
            bio_K1$prop_rhythmic))
cat(sprintf("K=2 pilot: n=%d, A2/A1 med=%.2f, prop_rhy=%.3f\n",
            attr(psi_K2, "n_pilot"),
            psi_K2$diagnostics$A2_over_A1_med, psi_K2$prop_rhythmic))

PERIOD <- 24; OMEGA0 <- 2*pi/PERIOD
NSIMS <- 100; NCORES <- 12; NGENES_SIM <- 2000; SEED_BASE <- 20260527
DISPLAY_N <- c(20L, 40L, 60L, 80L, 100L, 120L, 150L, 200L, 250L)
B_GRID    <- c(3L, 4L, 6L, 8L, 12L, 24L)

# Panel A: runSingleCohortPower per B with active equispaced cts
cat("\n== Panel A: K=1 via runSingleCohortPower, active per B ==\n")
pwr_K1 <- matrix(NA_real_, length(DISPLAY_N), length(B_GRID),
                  dimnames = list(paste0("N=",DISPLAY_N), paste0("B=",B_GRID)))
se_K1 <- pwr_K1
for (k in seq_along(B_GRID)) {
  B <- B_GRID[k]
  cts_active <- seq(0, PERIOD - PERIOD/B, length.out = B)
  design <- CircadianDesignOptions(
    sample_sizes = DISPLAY_N, nsims = NSIMS,
    design = "active", cts = cts_active
  )
  analysis <- CircadianAnalysisOptions(
    alpha = 0.05, p.adjust.method = "BH",
    r_strata = makeAdaptiveRStrata(bio_K1, bin_width = 0.25)
  )
  res <- runSingleCohortPower(bio_K1, design, analysis,
                                methods = "DCP", plot = FALSE,
                                verbose = FALSE, mc.cores = NCORES)
  # Compute marginal power per N at FDR=0.05 across sims
  for (j in seq_along(DISPLAY_N)) {
    rep_pwr <- numeric(NSIMS)
    for (s in seq_len(NSIMS)) {
      rv <- res$r_values_list[[j]][[s]]; rhy <- rv > 0
      pv <- res$pvalues[j, , s]; pv[is.na(pv)] <- 1
      q  <- p.adjust(pv, "BH")
      rep_pwr[s] <- if (sum(rhy) > 0) sum(q <= 0.05 & rhy) / sum(rhy) else NA
    }
    pwr_K1[j, k] <- mean(rep_pwr, na.rm = TRUE)
    se_K1[j, k]  <- sd(rep_pwr, na.rm = TRUE) / sqrt(sum(is.finite(rep_pwr)))
  }
  cat(sprintf("  B=%2d: %s\n", B,
              paste(sprintf("%5.1f%%", 100*pwr_K1[, k]), collapse = " ")))
}

# Panel B: simCircadianSingleCohort2H per B with active equispaced ts (matches Fig 5)
cat("\n== Panel B: K=2 via simCircadianSingleCohort2H, active per B ==\n")
B_K2 <- B_GRID
pwr_K2 <- matrix(NA_real_, length(DISPLAY_N), length(B_K2),
                  dimnames = list(paste0("N=",DISPLAY_N), paste0("B=",B_K2)))
se_K2 <- pwr_K2

active_times <- function(N, B) {
  m <- N %/% B
  ts <- rep(seq(0, PERIOD - PERIOD/B, length.out = B), each = m)
  if (length(ts) < N) ts <- c(ts, sample(seq(0, PERIOD - PERIOD/B, length.out = B), N - length(ts), replace = TRUE))
  if (length(ts) > N) ts <- ts[seq_len(N)]
  ts
}

sim_K2_one <- function(N, B, seed) {
  if (N < B || B < 5) return(NA_real_)
  ts <- active_times(N, B)
  if (length(unique(ts)) < 5) return(NA_real_)
  set.seed(seed)
  psi_local <- psi_K2; psi_local$ngenes <- NGENES_SIM
  sim <- simCircadianSingleCohort2H(psi_local, ts, seed = seed)
  Xs <- cbind(1, cos(OMEGA0*ts), sin(OMEGA0*ts),
              cos(2*OMEGA0*ts), sin(2*OMEGA0*ts))
  df_full <- length(ts) - 5L
  if (df_full <= 0) return(NA_real_)
  G <- nrow(sim$expr); pv <- numeric(G)
  for (g in seq_len(G)) {
    y <- sim$expr[g, ]; R0 <- sum((y - mean(y))^2)
    if (R0 <= 0) { pv[g] <- 1; next }
    f <- tryCatch(lm.fit(Xs, y), error = function(e) NULL)
    if (is.null(f)) { pv[g] <- 1; next }
    R <- sum(f$residuals^2)
    if (R <= 0 || any(is.na(f$coefficients))) { pv[g] <- 1; next }
    pv[g] <- pf(((R0-R)/4) / (R/df_full), 4, df_full, lower.tail = FALSE)
  }
  q <- p.adjust(pv, "BH"); rhy <- sim$is_rhythmic
  sum(q <= 0.05 & rhy) / max(1, sum(rhy))
}

for (k in seq_along(B_K2)) {
  B <- B_K2[k]
  for (j in seq_along(DISPLAY_N)) {
    N <- DISPLAY_N[j]
    v <- mclapply(seq_len(NSIMS), function(i)
      sim_K2_one(N, B, SEED_BASE + (N*1000L + B + i)*7919L),
      mc.cores = NCORES, mc.preschedule = FALSE)
    v <- unlist(v); v <- v[is.finite(v)]
    if (length(v) > 0) {
      pwr_K2[j, k] <- mean(v)
      se_K2[j, k]  <- sd(v) / sqrt(length(v))
    }
  }
  cat(sprintf("  B=%2d: %s\n", B,
              paste(sprintf("%5.1f%%", 100*pwr_K2[, k]), collapse = " ")))
}

saveRDS(list(N = DISPLAY_N, B = B_GRID,
              pwr_K1 = pwr_K1, se_K1 = se_K1,
              pwr_K2 = pwr_K2, se_K2 = se_K2),
        "output/two_harmonic/results/fig6_v8_data.rds")

cols <- c("#440154","#3b528b","#21918c","#5ec962","#fde725","#7d028c")[seq_along(B_GRID)]

draw_panel <- function(N, M, SE, B_grid, cols, letter, main_text, jitter_step = 2.2) {
  matplot(N, 100 * M, type = "n", ylim = c(0, 100), xlim = c(0, max(N) * 1.05),
          xlab = "Sample size (n)", ylab = "Power (%)", main = "")
  title(main = sprintf("%s   %s", letter, main_text),
        adj = 0.5, font.main = 2, cex.main = 1.10, line = 0.3)
  abline(h = 80, lty = 2, col = "grey50", lwd = 1.2); grid()
  for (k in seq_along(B_grid)) {
    xj <- N + (k - (length(B_grid) + 1)/2) * jitter_step
    lo <- 100 * pmax(0, M[, k] - SE[, k])
    hi <- 100 * pmin(1, M[, k] + SE[, k])
    ok <- is.finite(lo) & is.finite(hi) & SE[, k] > 0
    if (any(ok))
      arrows(xj[ok], lo[ok], xj[ok], hi[ok],
             code = 3, angle = 90, length = 0.025,
             col = cols[k], lwd = 1.0)
    M_k <- M[, k]; M_k[!is.finite(M_k)] <- 0
    lines(xj, 100 * M_k, type = "o", pch = 19, lwd = 1.6, col = cols[k], cex = 0.55)
  }
  legend("bottomright", paste0("B=", B_grid), col = cols[seq_along(B_grid)],
         lty = 1, pch = 19, lwd = 1.4, cex = 0.62,
         bty = "o", box.col = "grey70", box.lwd = 0.5, inset = 0.01,
         y.intersp = 0.82, bg = "white", seg.len = 1.4)
}

pdf("submission/figures/Fig6_active_BvsM.pdf", width = 9.0, height = 4.2)
par(mfrow = c(1, 2), mai = c(0.95, 0.85, 0.45, 0.15),
    mgp = c(2.5, 0.55, 0), oma = c(0, 0, 2.0, 0),
    cex.axis = 0.95, cex.lab = 1.10, cex.main = 1.10, font.main = 2)
draw_panel(DISPLAY_N, pwr_K1, se_K1, B_GRID, cols,
           "A", "K=1 Power, FDR 5% (K=1 pilot)")
draw_panel(DISPLAY_N, pwr_K2, se_K2, B_K2, cols,
           "B", "K=2 Power, FDR 5% (K=2 pilot)")
mtext("Active design B vs m trade-off (GTEx Liver pilot)",
      outer = TRUE, side = 3, line = 0.5, font = 2, cex = 1.20)
dev.off()
cat("\nSaved: submission/figures/Fig6_active_BvsM.pdf\n")
