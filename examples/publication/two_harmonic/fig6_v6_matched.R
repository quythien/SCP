#!/usr/bin/env Rscript
# =============================================================================
# Fig 6 v6 - 2-panel matched detector + simulator
#   Panel A: K=1 detector on K=1 simulator (A2 = 0), B-stratified active design
#   Panel B: K=2 detector on K=2 simulator (full),   B-stratified active design
# Jittered N offsets to separate overlapping curves, +/- 1 SE error bars.
# =============================================================================
setwd("/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")
sink(stderr())
source("code/utils.R"); source("code/options.R"); source("code/simulation.R")
source("code/estimation.R"); source("code/bootstrap_sim.R"); source("code/detection.R")
sink()
suppressPackageStartupMessages(library(parallel))

psi <- readRDS("output/two_harmonic/results/pilot_2h_GTExLiver_topK300.rds")
cat(sprintf("Pilot: GTEx Liver K=2 (A2/A1 med=%.2f)\n",
            psi$diagnostics$A2_over_A1_med))

PERIOD <- 24; OMEGA0 <- 2*pi/PERIOD
NSIMS <- 100; NCORES <- 12; NGENES <- 2000; SEED_BASE <- 20260527
DISPLAY_N <- c(20L, 40L, 60L, 80L, 100L, 120L, 150L, 200L, 250L)
B_K1 <- c(2L, 3L, 6L, 8L, 12L, 24L)   # B=2 below identifiability for K=1
B_K2 <- c(3L, 4L, 6L, 8L, 12L, 24L)   # B=3,4 below identifiability for K=2

make_K1_pilot <- function(psi) {
  psi$amplitude2 <- rep(0, length(psi$amplitude))
  psi$phase2     <- rep(0, length(psi$phase))
  psi
}

active_times <- function(N, B) {
  if (N < B) return(sort(sample(seq(0, PERIOD - PERIOD/B, length.out = B), N, replace = FALSE)))
  m <- N %/% B
  ts <- rep(seq(0, PERIOD - PERIOD/B, length.out = B), each = m)
  if (length(ts) < N) ts <- c(ts, sample(seq(0, PERIOD - PERIOD/B, length.out = B), N - length(ts), replace = TRUE))
  if (length(ts) > N) ts <- ts[seq_len(N)]
  ts
}

sim_one <- function(N, B, seed, pilot, K) {
  ts <- active_times(N, B); set.seed(seed)
  psi_local <- pilot; psi_local$ngenes <- NGENES
  sim <- simCircadianSingleCohort2H(psi_local, ts, seed = seed)
  Xs <- if (K == 1L) cbind(1, cos(OMEGA0*ts), sin(OMEGA0*ts))
        else cbind(1, cos(OMEGA0*ts), sin(OMEGA0*ts),
                   cos(2*OMEGA0*ts), sin(2*OMEGA0*ts))
  df_full <- length(ts) - (K*2 + 1L); if (df_full <= 0) return(NA_real_)
  G <- nrow(sim$expr); pv <- numeric(G)
  for (g in seq_len(G)) {
    y <- sim$expr[g, ]; R0 <- sum((y - mean(y))^2)
    if (R0 <= 0) { pv[g] <- 1; next }
    f <- tryCatch(lm.fit(Xs, y), error = function(e) NULL)
    if (is.null(f)) { pv[g] <- 1; next }
    R <- sum(f$residuals^2)
    if (R <= 0 || any(is.na(f$coefficients))) { pv[g] <- 1; next }
    pv[g] <- pf(((R0-R)/(K*2)) / (R/df_full), K*2, df_full, lower.tail = FALSE)
  }
  q <- p.adjust(pv, "BH"); rhy <- sim$is_rhythmic
  sum(q <= 0.05 & rhy) / max(1, sum(rhy))
}

run_grid <- function(N_grid, B_grid, K, pilot) {
  M  <- matrix(NA_real_, length(N_grid), length(B_grid))
  SE <- matrix(NA_real_, length(N_grid), length(B_grid))
  for (j in seq_along(N_grid)) {
    for (k in seq_along(B_grid)) {
      N <- N_grid[j]; B <- B_grid[k]; if (N < B) next
      vals <- mclapply(seq_len(NSIMS), function(i)
        sim_one(N, B, SEED_BASE + (N*1000L + B + i)*7919L, pilot, K),
        mc.cores = NCORES, mc.preschedule = FALSE)
      v <- unlist(vals); v <- v[is.finite(v)]
      if (length(v) > 0) {
        M[j, k]  <- mean(v)
        SE[j, k] <- sd(v) / sqrt(length(v))
      }
    }
    cat(sprintf("  N=%3d (K=%d): %s\n", N_grid[j], K,
                paste(sprintf("%5.1f%%", 100*M[j,]), collapse=" ")))
  }
  list(M = M, SE = SE)
}

psi_K1 <- make_K1_pilot(psi)
cat("\n== Panel A: K=1 detector + K=1 simulator ==\n")
res_K1 <- run_grid(DISPLAY_N, B_K1, K = 1L, psi_K1)
cat("\n== Panel B: K=2 detector + K=2 simulator ==\n")
res_K2 <- run_grid(DISPLAY_N, B_K2, K = 2L, psi)

saveRDS(list(N = DISPLAY_N, B_K1 = B_K1, B_K2 = B_K2,
              pwr_K1 = res_K1$M, se_K1 = res_K1$SE,
              pwr_K2 = res_K2$M, se_K2 = res_K2$SE),
        "output/two_harmonic/results/fig6_v6_data.rds")

cols_K1 <- c("#7d028c","#440154","#3b528b","#21918c","#5ec962","#fde725")
cols_K2 <- c("#7d028c","#440154","#3b528b","#21918c","#5ec962","#fde725")

draw_panel <- function(N, M, SE, B_grid, cols, letter, main_text, jitter_step = 1.8) {
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
    lines(xj, 100 * M[, k], type = "o", pch = 19, lwd = 1.6,
          col = cols[k], cex = 0.55)
  }
  legend("bottomright", paste0("B=", B_grid), col = cols[seq_along(B_grid)],
         lty = 1, pch = 19, lwd = 1.4, cex = 0.62,
         bty = "o", box.col = "grey70", box.lwd = 0.5, inset = 0.01,
         y.intersp = 0.82, bg = "white", seg.len = 1.4)
}

pdf("submission/figures/Fig6_active_BvsM.pdf", width = 8.5, height = 4.0)
par(mfrow = c(1, 2), mai = c(0.95, 0.85, 0.45, 0.15),
    mgp = c(2.5, 0.55, 0), oma = c(0, 0, 2.0, 0),
    cex.axis = 0.95, cex.lab = 1.10, cex.main = 1.10, font.main = 2)

draw_panel(DISPLAY_N, res_K1$M, res_K1$SE, B_K1, cols_K1,
           "A", "K=1 Power (K=1 simulator)", jitter_step = 2.2)
draw_panel(DISPLAY_N, res_K2$M, res_K2$SE, B_K2, cols_K2,
           "B", "K=2 Power (K=2 simulator)", jitter_step = 2.2)

mtext("Active design B vs m trade-off (GTEx Liver pilot)",
      outer = TRUE, side = 3, line = 0.5, font = 2, cex = 1.20)
dev.off()
cat("\nSaved: submission/figures/Fig6_active_BvsM.pdf\n")
