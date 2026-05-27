#!/usr/bin/env Rscript
# =============================================================================
# Fig 6 v4: B-invariance + K=2 phase advantage (GTEx Liver K=2)
#
# Panel A: K=1 power vs N, B-stratified active     — B-invariance theorem
# Panel B: K=2 power vs N, B-stratified active     — K=2 also B-invariant
# Panel C: Phase MSE — 4 curves: K=1 active, K=2 active, K=1 passive, K=2 passive
#          (fixed B=12 for active; passive uses empirical Liver TOD distribution)
# =============================================================================
setwd("/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")
sink(stderr())
source("code/utils.R"); source("code/options.R"); source("code/simulation.R")
source("code/estimation.R"); source("code/bootstrap_sim.R"); source("code/detection.R")
sink()
suppressPackageStartupMessages(library(parallel))

psi <- readRDS("output/two_harmonic/results/pilot_2h_GTExLiver_topK300.rds")
cat(sprintf("GTEx Liver K=2 pilot: A2/A1 med=%.2f, prop_rhythmic=%.3f, n=%d\n",
            psi$diagnostics$A2_over_A1_med, psi$prop_rhythmic,
            attr(psi, "n_pilot")))

PERIOD <- 24; OMEGA0 <- 2*pi/PERIOD
DISPLAY_N_AB <- c(20L, 40L, 60L, 80L, 100L, 120L, 150L, 200L, 250L)
DISPLAY_N_C  <- c(20L, 40L, 60L, 80L, 100L, 150L, 200L, 250L)
B_K1 <- c(3L, 6L, 8L, 12L, 24L)
B_K2 <- c(6L, 8L, 12L, 24L)
B_FIXED_C <- 12L
NSIMS <- 100; NCORES <- 12; SEED_BASE <- 20260527; NGENES <- 2000

active_times <- function(N, B) {
  if (N < B) return(sort(sample(seq(0, PERIOD - PERIOD/B, length.out = B), N, replace = FALSE)))
  m <- N %/% B
  ts <- rep(seq(0, PERIOD - PERIOD/B, length.out = B), each = m)
  if (length(ts) < N) ts <- c(ts, sample(seq(0, PERIOD - PERIOD/B, length.out = B),
                                          N - length(ts), replace = TRUE))
  if (length(ts) > N) ts <- ts[seq_len(N)]
  ts
}
passive_times <- function(N, pilot_cts) sample(pilot_cts, N, replace = TRUE)

circ_mse <- function(phi_hat, phi_true) {
  d <- (phi_hat - phi_true) %% PERIOD; d <- pmin(d, PERIOD - d)
  median(d^2, na.rm = TRUE)
}

# Single sim: returns power and phase MSE for K=1 and K=2 detectors
sim_one <- function(N, ts, seed) {
  set.seed(seed)
  psi_local <- psi; psi_local$ngenes <- NGENES
  sim <- simCircadianSingleCohort2H(psi_local, ts, seed = seed)
  G <- nrow(sim$expr)
  X1 <- cbind(1, cos(OMEGA0*ts), sin(OMEGA0*ts))
  X2 <- cbind(X1, cos(2*OMEGA0*ts), sin(2*OMEGA0*ts))
  df1 <- length(ts) - 3L; df2 <- length(ts) - 5L
  if (df1 <= 0 || df2 <= 0) return(NULL)
  pv1 <- pv2 <- phi1_K1 <- phi1_K2 <- numeric(G)
  for (g in seq_len(G)) {
    y <- sim$expr[g, ]; R0 <- sum((y - mean(y))^2)
    if (R0 <= 0) { pv1[g] <- pv2[g] <- 1; phi1_K1[g] <- phi1_K2[g] <- NA; next }
    f1 <- lm.fit(X1, y); R1 <- sum(f1$residuals^2)
    f2 <- lm.fit(X2, y); R2 <- sum(f2$residuals^2)
    if (R1 <= 0 || R2 <= 0) { pv1[g] <- pv2[g] <- 1; phi1_K1[g] <- phi1_K2[g] <- NA; next }
    pv1[g] <- pf(((R0-R1)/2) / (R1/df1), 2, df1, lower.tail = FALSE)
    pv2[g] <- pf(((R0-R2)/4) / (R2/df2), 4, df2, lower.tail = FALSE)
    phi1_K1[g] <- (atan2(f1$coefficients[3], f1$coefficients[2]) / OMEGA0) %% PERIOD
    phi1_K2[g] <- (atan2(f2$coefficients[3], f2$coefficients[2]) / OMEGA0) %% PERIOD
  }
  q1 <- p.adjust(pv1, "BH"); q2 <- p.adjust(pv2, "BH")
  rhy <- sim$is_rhythmic
  list(power_K1 = sum(q1 <= 0.05 & rhy) / max(1, sum(rhy)),
       power_K2 = sum(q2 <= 0.05 & rhy) / max(1, sum(rhy)),
       mse_K1   = if (sum(rhy) > 0) circ_mse(phi1_K1[rhy], sim$phi1[rhy]) else NA,
       mse_K2   = if (sum(rhy) > 0) circ_mse(phi1_K2[rhy], sim$phi1[rhy]) else NA)
}

# Run grid of (N, B) → power[K1/K2]
run_active_power <- function(N_grid, B_grid, K) {
  M <- matrix(NA_real_, length(N_grid), length(B_grid),
              dimnames = list(paste0("N=", N_grid), paste0("B=", B_grid)))
  for (j in seq_along(N_grid)) {
    for (k in seq_along(B_grid)) {
      N <- N_grid[j]; B <- B_grid[k]
      if (N < B) next
      res <- mclapply(seq_len(NSIMS), function(i)
        sim_one(N, active_times(N, B), SEED_BASE + (N*1000L + B + i)*7919L),
        mc.cores = NCORES, mc.preschedule = FALSE)
      res <- res[vapply(res, is.list, logical(1))]
      if (length(res) > 0) {
        M[j, k] <- mean(sapply(res, function(x)
          if (K == 1L) x$power_K1 else x$power_K2), na.rm = TRUE)
      }
    }
  }
  M
}

cat("\n== Panel A: K=1 power (B-stratified active) ==\n")
pwr_K1 <- run_active_power(DISPLAY_N_AB, B_K1, K = 1L)
for (j in seq_along(DISPLAY_N_AB))
  cat(sprintf("  N=%3d: %s\n", DISPLAY_N_AB[j],
              paste(sprintf("%5.1f%%", 100*pwr_K1[j,]), collapse = " ")))

cat("\n== Panel B: K=2 power (B-stratified active) ==\n")
pwr_K2 <- run_active_power(DISPLAY_N_AB, B_K2, K = 2L)
for (j in seq_along(DISPLAY_N_AB))
  cat(sprintf("  N=%3d: %s\n", DISPLAY_N_AB[j],
              paste(sprintf("%5.1f%%", 100*pwr_K2[j,]), collapse = " ")))

cat("\n== Panel C: Phase MSE — active vs passive, K=1 vs K=2 ==\n")
mse_active_K1  <- mse_active_K2  <- numeric(length(DISPLAY_N_C))
mse_passive_K1 <- mse_passive_K2 <- numeric(length(DISPLAY_N_C))
for (j in seq_along(DISPLAY_N_C)) {
  N <- DISPLAY_N_C[j]
  res_a <- mclapply(seq_len(NSIMS), function(i)
    sim_one(N, active_times(N, B_FIXED_C), SEED_BASE + (N*1000L + 99L + i)*7919L),
    mc.cores = NCORES, mc.preschedule = FALSE)
  res_p <- mclapply(seq_len(NSIMS), function(i)
    sim_one(N, passive_times(N, psi$cts),  SEED_BASE + (N*1000L + 88L + i)*7919L),
    mc.cores = NCORES, mc.preschedule = FALSE)
  res_a <- res_a[vapply(res_a, is.list, logical(1))]
  res_p <- res_p[vapply(res_p, is.list, logical(1))]
  if (length(res_a) > 0) {
    mse_active_K1[j] <- mean(sapply(res_a, function(x) x$mse_K1), na.rm = TRUE)
    mse_active_K2[j] <- mean(sapply(res_a, function(x) x$mse_K2), na.rm = TRUE)
  }
  if (length(res_p) > 0) {
    mse_passive_K1[j] <- mean(sapply(res_p, function(x) x$mse_K1), na.rm = TRUE)
    mse_passive_K2[j] <- mean(sapply(res_p, function(x) x$mse_K2), na.rm = TRUE)
  }
  cat(sprintf("  N=%3d: active K1=%.2f K2=%.2f | passive K1=%.2f K2=%.2f\n",
              N, mse_active_K1[j], mse_active_K2[j],
              mse_passive_K1[j], mse_passive_K2[j]))
}

saveRDS(list(pwr_K1 = pwr_K1, pwr_K2 = pwr_K2,
              N_AB = DISPLAY_N_AB, N_C = DISPLAY_N_C,
              B_K1 = B_K1, B_K2 = B_K2,
              mse_active_K1 = mse_active_K1,
              mse_active_K2 = mse_active_K2,
              mse_passive_K1 = mse_passive_K1,
              mse_passive_K2 = mse_passive_K2),
        "output/two_harmonic/results/fig6_v4_data.rds")

cols <- c("#440154","#3b528b","#21918c","#5ec962","#fde725")
pdf("submission/figures/Fig6_v4_phase_story.pdf", width = 11.0, height = 3.7)
par(mfrow = c(1, 3), mai = c(0.95, 0.85, 0.50, 0.10),
    mgp = c(2.6, 0.55, 0), oma = c(0, 0, 2.0, 0),
    cex.axis = 0.90, cex.lab = 1.05, cex.main = 1.05, font.main = 2)

# Panel A
matplot(DISPLAY_N_AB, 100 * pwr_K1, type = "o", pch = 19, lwd = 1.8,
        col = cols[seq_along(B_K1)], lty = 1, cex = 0.55,
        ylim = c(0, 100), xlim = c(0, max(DISPLAY_N_AB) * 1.05),
        xlab = "Sample size (n)", ylab = "Power (%)", main = "")
title(main = "A   K=1 Power (active)",
      adj = 0.5, font.main = 2, cex.main = 1.05, line = 0.3)
abline(h = 80, lty = 2, col = "grey50"); grid()
legend("bottomright", paste0("B=", B_K1), col = cols[seq_along(B_K1)],
       lty = 1, pch = 19, lwd = 1.4, cex = 0.55,
       bty = "o", box.col = "grey70", inset = 0.01, y.intersp = 0.78, bg = "white")

# Panel B
matplot(DISPLAY_N_AB, 100 * pwr_K2, type = "o", pch = 19, lwd = 1.8,
        col = cols[seq_along(B_K2)], lty = 1, cex = 0.55,
        ylim = c(0, 100), xlim = c(0, max(DISPLAY_N_AB) * 1.05),
        xlab = "Sample size (n)", ylab = "Power (%)", main = "")
title(main = "B   K=2 Power (active)",
      adj = 0.5, font.main = 2, cex.main = 1.05, line = 0.3)
abline(h = 80, lty = 2, col = "grey50"); grid()
legend("bottomright", paste0("B=", B_K2), col = cols[seq_along(B_K2)],
       lty = 1, pch = 19, lwd = 1.4, cex = 0.55,
       bty = "o", box.col = "grey70", inset = 0.01, y.intersp = 0.78, bg = "white")

# Panel C: 4-curve Phase MSE comparison
ymax_C <- max(c(mse_active_K1, mse_active_K2,
                mse_passive_K1, mse_passive_K2), na.rm = TRUE) * 1.05
plot(DISPLAY_N_C, mse_active_K1, type = "o", pch = 19, lwd = 1.8,
     col = "darkorange", cex = 0.55, lty = 1,
     ylim = c(0, ymax_C), xlim = c(0, max(DISPLAY_N_C) * 1.05),
     xlab = "Sample size (n)",
     ylab = expression(median ~ "(" * hat(phi)[1] - phi[1*","*true] * ")"^2 ~ "(h"^2*")"),
     main = "")
lines(DISPLAY_N_C, mse_active_K2,  type = "o", pch = 19, lwd = 1.8,
      col = "steelblue", cex = 0.55, lty = 1)
lines(DISPLAY_N_C, mse_passive_K1, type = "o", pch = 17, lwd = 1.8,
      col = "darkorange", cex = 0.65, lty = 2)
lines(DISPLAY_N_C, mse_passive_K2, type = "o", pch = 17, lwd = 1.8,
      col = "steelblue", cex = 0.65, lty = 2)
title(main = "C   Phase MSE: active vs passive",
      adj = 0.5, font.main = 2, cex.main = 1.05, line = 0.3)
grid()
legend("topright",
       c("K=1 active (B=12)", "K=2 active (B=12)",
         "K=1 passive (empirical TOD)", "K=2 passive (empirical TOD)"),
       col = c("darkorange", "steelblue", "darkorange", "steelblue"),
       pch = c(19, 19, 17, 17), lty = c(1, 1, 2, 2), lwd = 1.6,
       cex = 0.55, bty = "o", box.col = "grey70", box.lwd = 0.5,
       inset = 0.01, y.intersp = 0.85, bg = "white", seg.len = 1.5)

mtext("Active design B-invariance + K=2 phase advantage (GTEx Liver K=2)",
      outer = TRUE, side = 3, line = 0.2, font = 2, cex = 1.10)
dev.off()
cat("\nSaved: submission/figures/Fig6_v4_phase_story.pdf\n")
