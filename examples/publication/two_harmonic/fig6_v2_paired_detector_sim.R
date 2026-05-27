#!/usr/bin/env Rscript
# =============================================================================
# Fig 6 v2: 3 panels, each panel uses matched detector + simulator.
#   Panel A: K=1 detector + K=1 simulator (A2 = 0)
#   Panel B: K=2 detector + K=2 simulator (full two-harmonic data)
#   Panel C: K=1 phase MSE on K=1 simulator data
#
# Two pilot variants rendered for visual comparison:
#   - GTEx Liver (low SNR, passive-like)
#   - Mure 2018 Baboon KIM (active design)
# =============================================================================

setwd("/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")
sink(stderr())
source("code/utils.R"); source("code/options.R"); source("code/simulation.R")
source("code/estimation.R"); source("code/bootstrap_sim.R"); source("code/detection.R")
sink()
suppressPackageStartupMessages(library(parallel))

PERIOD <- 24; OMEGA0 <- 2*pi/PERIOD
NSIMS  <- 100
NCORES <- 12
NGENES_SIM <- 2000
SEED_BASE <- 20260526

# Build a K=1-only variant of a K=2 pilot by zeroing out the second harmonic.
make_K1_pilot <- function(psi) {
  psi$amplitude2 <- rep(0, length(psi$amplitude))
  psi$phase2     <- rep(0, length(psi$phase))
  psi
}

active_times <- function(N, B) {
  if (N < B) return(sort(sample(seq(0, PERIOD - PERIOD/B, length.out = B), N, replace = FALSE)))
  m  <- N %/% B
  ts <- rep(seq(0, PERIOD - PERIOD/B, length.out = B), each = m)
  if (length(ts) < N) ts <- c(ts, sample(seq(0, PERIOD - PERIOD/B, length.out = B),
                                          N - length(ts), replace = TRUE))
  if (length(ts) > N) ts <- ts[seq_len(N)]
  ts
}

# Generic sim: given pilot psi (already K=1 or K=2 form) and detector K, return
# pvalue, phase estimates, and is_rhythmic per gene.
sim_one <- function(N, B, seed, psi, detector_K) {
  ts <- active_times(N, B)
  set.seed(seed)
  psi_local <- psi; psi_local$ngenes <- NGENES_SIM
  sim <- simCircadianSingleCohort2H(psi_local, ts, seed = seed)
  if (detector_K == 1L) {
    Xs <- cbind(1, cos(OMEGA0*ts), sin(OMEGA0*ts))
    df_full <- length(ts) - 3L
  } else {
    Xs <- cbind(1, cos(OMEGA0*ts), sin(OMEGA0*ts),
                cos(2*OMEGA0*ts), sin(2*OMEGA0*ts))
    df_full <- length(ts) - 5L
  }
  if (df_full <= 0) return(NULL)
  G <- nrow(sim$expr); pv <- numeric(G); phi_hat <- rep(NA_real_, G)
  for (g in seq_len(G)) {
    y <- sim$expr[g, ]; R0 <- sum((y-mean(y))^2)
    if (R0 <= 0) { pv[g] <- 1; next }
    f <- lm.fit(Xs, y); R <- sum(f$residuals^2)
    if (R <= 0) { pv[g] <- 1; next }
    pv[g] <- pf(((R0-R)/(detector_K*2)) / (R/df_full),
                detector_K*2, df_full, lower.tail = FALSE)
    if (detector_K == 1L)
      phi_hat[g] <- (atan2(f$coefficients[3], f$coefficients[2]) / OMEGA0) %% PERIOD
  }
  list(pvalue = pv, phi_hat = phi_hat,
       phi_true = sim$phi1, is_rhythmic = sim$is_rhythmic)
}

# Circular phase MSE in hours^2
circ_mse <- function(phi_hat, phi_true) {
  d <- (phi_hat - phi_true) %% PERIOD
  d <- pmin(d, PERIOD - d)
  median(d^2, na.rm = TRUE)
}

run_power <- function(psi, N_grid, B_grid, detector_K, label) {
  cat(sprintf("\n== Power %s ==\n", label))
  M <- matrix(NA_real_, nrow = length(N_grid), ncol = length(B_grid),
              dimnames = list(paste0("N=", N_grid), paste0("B=", B_grid)))
  for (j in seq_along(N_grid)) {
    for (k in seq_along(B_grid)) {
      N <- N_grid[j]; B <- B_grid[k]
      if (N < B) next
      sims <- mclapply(seq_len(NSIMS), function(i)
        sim_one(N, B, SEED_BASE + (N*1000L + B + i)*7919L, psi, detector_K),
        mc.cores = NCORES, mc.preschedule = FALSE)
      sims <- sims[vapply(sims, is.list, logical(1))]
      if (length(sims) > 0)
        M[j,k] <- mean(sapply(sims, function(s) {
          q <- p.adjust(s$pvalue, "BH")
          sum(q <= 0.05 & s$is_rhythmic) / max(1, sum(s$is_rhythmic))
        }), na.rm = TRUE)
    }
    cat(sprintf("  N=%3d: %s\n", N_grid[j],
                paste(sprintf("%5.1f%%", 100*M[j,]), collapse = " ")))
  }
  M
}

run_phase_mse <- function(psi, N_grid, B_grid, label) {
  cat(sprintf("\n== Phase MSE %s ==\n", label))
  M <- matrix(NA_real_, nrow = length(N_grid), ncol = length(B_grid),
              dimnames = list(paste0("N=", N_grid), paste0("B=", B_grid)))
  for (j in seq_along(N_grid)) {
    for (k in seq_along(B_grid)) {
      N <- N_grid[j]; B <- B_grid[k]
      if (N < B) next
      sims <- mclapply(seq_len(NSIMS), function(i)
        sim_one(N, B, SEED_BASE + (N*1000L + B + i)*7919L, psi, detector_K = 1L),
        mc.cores = NCORES, mc.preschedule = FALSE)
      sims <- sims[vapply(sims, is.list, logical(1))]
      if (length(sims) > 0)
        M[j,k] <- mean(sapply(sims, function(s) {
          rhy <- s$is_rhythmic
          if (sum(rhy) == 0) return(NA_real_)
          circ_mse(s$phi_hat[rhy], s$phi_true[rhy])
        }), na.rm = TRUE)
    }
    cat(sprintf("  N=%3d: %s\n", N_grid[j],
                paste(sprintf("%5.2f", M[j,]), collapse = " ")))
  }
  M
}

render_fig6 <- function(pwr_K1, pwr_K2, mse_K1, N_grid, B_grid_K1, B_grid_K2,
                        out_pdf, source_label) {
  pdf(out_pdf, width = 11.5, height = 3.8)
  par(mfrow = c(1, 3), mai = c(0.85, 0.85, 0.55, 0.15),
      mgp = c(2.6, 0.55, 0), oma = c(0, 0, 1.4, 0),
      cex.axis = 0.95, cex.lab = 1.10, cex.main = 1.05, font.main = 2)
  cols_K1 <- c("#440154","#3b528b","#21918c","#5ec962","#fde725")[seq_along(B_grid_K1)]
  cols_K2 <- c("#440154","#3b528b","#21918c","#5ec962","#fde725")[seq_along(B_grid_K2)]

  # Panel A: K=1 power
  matplot(N_grid, 100 * pwr_K1, type = "b", pch = 19, lwd = 1.8,
          col = cols_K1, lty = 1, ylim = c(0, 100),
          xlim = c(0, max(N_grid) * 1.05),
          xlab = "Sample size (n)", ylab = "Power (%)", main = "")
  title(main = "A   K=1 Power (K=1 simulator)",
        adj = 0.5, font.main = 2, cex.main = 1.05, line = 0.3)
  abline(h = 80, lty = 2, col = "grey50", lwd = 1.2); grid()
  legend("bottomright", paste0("B=", B_grid_K1),
         col = cols_K1, lty = 1, pch = 19, lwd = 1.5, cex = 0.58,
         bty = "o", box.col = "grey70", inset = 0.01, y.intersp = 0.78, bg = "white")

  # Panel B: K=2 power
  matplot(N_grid, 100 * pwr_K2, type = "b", pch = 19, lwd = 1.8,
          col = cols_K2, lty = 1, ylim = c(0, 100),
          xlim = c(0, max(N_grid) * 1.05),
          xlab = "Sample size (n)", ylab = "Power (%)", main = "")
  title(main = "B   K=2 Power (K=2 simulator)",
        adj = 0.5, font.main = 2, cex.main = 1.05, line = 0.3)
  abline(h = 80, lty = 2, col = "grey50", lwd = 1.2); grid()
  legend("bottomright", paste0("B=", B_grid_K2),
         col = cols_K2, lty = 1, pch = 19, lwd = 1.5, cex = 0.58,
         bty = "o", box.col = "grey70", inset = 0.01, y.intersp = 0.78, bg = "white")

  # Panel C: K=1 phase MSE
  matplot(N_grid, mse_K1, type = "b", pch = 19, lwd = 1.8,
          col = cols_K1, lty = 1,
          xlim = c(0, max(N_grid) * 1.05),
          ylim = c(0, max(mse_K1, na.rm = TRUE) * 1.05),
          xlab = "Sample size (n)",
          ylab = expression(median ~ "(" * hat(phi) - phi[true] * ")"^2 ~ "(h"^2*")"),
          main = "")
  title(main = "C   Phase MSE (K=1 simulator)",
        adj = 0.5, font.main = 2, cex.main = 1.05, line = 0.3)
  grid()
  legend("topright", paste0("B=", B_grid_K1),
         col = cols_K1, lty = 1, pch = 19, lwd = 1.5, cex = 0.58,
         bty = "o", box.col = "grey70", inset = 0.01, y.intersp = 0.78, bg = "white")

  mtext(sprintf("Active design B-invariance and phase estimation: %s", source_label),
        outer = TRUE, side = 3, line = 0.1, font = 2, cex = 1.00)
  dev.off()
  cat(sprintf("Wrote: %s\n", out_pdf))
}

run_pilot <- function(pilot_path, N_grid, B_grid_K1, B_grid_K2, label, out_pdf) {
  cat(sprintf("\n###### Pilot: %s ######\n", label))
  psi_K2 <- readRDS(pilot_path)
  psi_K1 <- make_K1_pilot(psi_K2)
  cat(sprintf("Pilot: top_k=%d, A2/A1 med=%.2f, prop_rhythmic=%.3f\n",
              psi_K2$diagnostics$top_k_used,
              psi_K2$diagnostics$A2_over_A1_med,
              psi_K2$prop_rhythmic))
  pwr_K1 <- run_power(psi_K1, N_grid, B_grid_K1, detector_K = 1L, paste(label, "K=1"))
  pwr_K2 <- run_power(psi_K2, N_grid, B_grid_K2, detector_K = 2L, paste(label, "K=2"))
  mse_K1 <- run_phase_mse(psi_K1, N_grid, B_grid_K1, paste(label, "phase MSE"))
  render_fig6(pwr_K1, pwr_K2, mse_K1, N_grid, B_grid_K1, B_grid_K2,
              out_pdf, label)
  list(pwr_K1 = pwr_K1, pwr_K2 = pwr_K2, mse_K1 = mse_K1,
       N = N_grid, B_K1 = B_grid_K1, B_K2 = B_grid_K2)
}

# GTEx Liver (passive-like, wider N range)
res_liver <- run_pilot(
  "output/two_harmonic/results/pilot_2h_GTExLiver_topK300.rds",
  N_grid = c(20L, 40L, 60L, 80L, 100L, 120L, 150L, 200L, 250L),
  B_grid_K1 = c(3L, 6L, 8L, 12L, 24L),
  B_grid_K2 = c(6L, 8L, 12L, 24L),
  label = "GTEx Liver",
  out_pdf = "submission/figures/Fig6_v2_GTExLiver.pdf"
)

# Baboon KIM (active design, smaller N range)
res_baboon <- run_pilot(
  "output/two_harmonic/results/pilot_2h_BaboonKIM.rds",
  N_grid = c(12L, 16L, 20L, 24L, 32L, 48L, 72L, 96L),
  B_grid_K1 = c(3L, 6L, 8L, 12L, 24L),
  B_grid_K2 = c(6L, 8L, 12L, 24L),
  label = "Mure 2018 Baboon KIM",
  out_pdf = "submission/figures/Fig6_v2_BaboonKIM.pdf"
)

saveRDS(list(GTExLiver = res_liver, BaboonKIM = res_baboon),
        "output/two_harmonic/results/fig6_v2_data.rds")
cat("\n=== Done ===\n")
