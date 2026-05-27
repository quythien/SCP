#!/usr/bin/env Rscript
# =============================================================================
# Fig 6 v3: B vs m trade-off on Putamen Control (n=59, passive-derived, brain).
#   Panel A: K=1 detector + K=1 simulator (matched)
#   Panel B: K=2 detector + K=2 simulator with A2=0 (matched-detector framing,
#            no real 2nd harmonic in Putamen pilot)
#   Panel C: K=1 phase MSE
# =============================================================================
setwd("/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")
sink(stderr())
source("code/utils.R"); source("code/options.R"); source("code/simulation.R")
source("code/estimation.R"); source("code/bootstrap_sim.R"); source("code/detection.R")
sink()
suppressPackageStartupMessages(library(parallel))

bio1 <- readRDS("data/gse160521_putamen_ctrl_pilot.rds")
cat(sprintf("Putamen Control K=1 pilot: n=%d, prop_rhythmic=%.3f, median r~=%.2f\n",
            length(bio1$cts), bio1$prop_rhythmic,
            median(bio1$amplitude / median(bio1$sigma_rhythmic, na.rm = TRUE),
                    na.rm = TRUE)))

PERIOD <- 24; OMEGA0 <- 2*pi/PERIOD
NGENES <- 2000; NSIMS <- 100; NCORES <- 12; SEED_BASE <- 20260527

active_times <- function(N, B) {
  if (N < B) return(sort(sample(seq(0, PERIOD - PERIOD/B, length.out = B), N, replace = FALSE)))
  m  <- N %/% B
  ts <- rep(seq(0, PERIOD - PERIOD/B, length.out = B), each = m)
  if (length(ts) < N) ts <- c(ts, sample(seq(0, PERIOD - PERIOD/B, length.out = B),
                                          N - length(ts), replace = TRUE))
  if (length(ts) > N) ts <- ts[seq_len(N)]
  ts
}

# K=1 simulator: draw genes from K=1 pilot, generate cosinor data
sim_K1 <- function(N, B, seed) {
  ts <- active_times(N, B); set.seed(seed)
  n_rhy <- round(NGENES * bio1$prop_rhythmic)
  ji <- sample.int(length(bio1$amplitude), n_rhy, replace = TRUE)
  A   <- bio1$amplitude[ji]
  phi <- sample(bio1$phase, n_rhy, replace = TRUE)
  sg  <- sample(bio1$sigma_rhythmic, n_rhy, replace = TRUE)
  mu  <- rnorm(NGENES, mean = median(bio1$lBaselineExpr, na.rm = TRUE),
               sd = sd(bio1$lBaselineExpr, na.rm = TRUE))
  is_rhy <- c(rep(TRUE, n_rhy), rep(FALSE, NGENES - n_rhy))
  s_full <- c(sg, sample(exp(bio1$lOD), NGENES - n_rhy, replace = TRUE))
  A_full <- c(A, rep(0, NGENES - n_rhy)); phi_full <- c(phi, rep(0, NGENES - n_rhy))
  expr <- matrix(rnorm(NGENES * length(ts), sd = 1) * s_full, NGENES, length(ts))
  for (g in seq_len(n_rhy))
    expr[g, ] <- expr[g, ] + mu[g] + A_full[g] * cos(OMEGA0 * (ts - phi_full[g]))
  for (g in (n_rhy+1):NGENES) expr[g, ] <- expr[g, ] + mu[g]
  list(expr = expr, ts = ts, is_rhythmic = is_rhy, phi1 = phi_full)
}

# Detection given expression matrix and K
detect_K <- function(expr, ts, K) {
  Xs <- if (K == 1L) cbind(1, cos(OMEGA0*ts), sin(OMEGA0*ts))
        else cbind(1, cos(OMEGA0*ts), sin(OMEGA0*ts),
                   cos(2*OMEGA0*ts), sin(2*OMEGA0*ts))
  df_full <- length(ts) - (K * 2 + 1L)
  if (df_full <= 0) return(list(pvalue = rep(1, nrow(expr)), phi_hat = rep(NA, nrow(expr))))
  pv <- numeric(nrow(expr)); phi_hat <- rep(NA_real_, nrow(expr))
  for (g in seq_len(nrow(expr))) {
    y <- expr[g, ]; R0 <- sum((y - mean(y))^2)
    if (R0 <= 0) { pv[g] <- 1; next }
    f <- lm.fit(Xs, y); R <- sum(f$residuals^2)
    if (R <= 0) { pv[g] <- 1; next }
    pv[g] <- pf(((R0-R) / (K*2)) / (R/df_full), K*2, df_full, lower.tail = FALSE)
    if (K == 1L)
      phi_hat[g] <- (atan2(f$coefficients[3], f$coefficients[2]) / OMEGA0) %% PERIOD
  }
  list(pvalue = pv, phi_hat = phi_hat)
}

circ_mse <- function(phi_hat, phi_true) {
  d <- (phi_hat - phi_true) %% PERIOD; d <- pmin(d, PERIOD - d)
  median(d^2, na.rm = TRUE)
}

DISPLAY_N <- c(20L, 30L, 40L, 50L, 60L, 80L, 100L, 120L, 150L, 200L)
B_K1 <- c(3L, 6L, 8L, 12L, 24L)
B_K2 <- c(6L, 8L, 12L, 24L)

run_grid <- function(N_grid, B_grid, sim_fn, K, measure = c("power", "phase")) {
  measure <- match.arg(measure)
  M <- matrix(NA_real_, length(N_grid), length(B_grid),
              dimnames = list(paste0("N=", N_grid), paste0("B=", B_grid)))
  for (j in seq_along(N_grid)) {
    for (k in seq_along(B_grid)) {
      N <- N_grid[j]; B <- B_grid[k]
      if (N < B) next
      vals <- mclapply(seq_len(NSIMS), function(i) {
        sim <- sim_fn(N, B, SEED_BASE + (N*1000L + B + i) * 7919L)
        det <- detect_K(sim$expr, sim$ts, K)
        if (measure == "power") {
          q <- p.adjust(det$pvalue, "BH")
          sum(q <= 0.05 & sim$is_rhythmic) / max(1, sum(sim$is_rhythmic))
        } else {
          rhy <- sim$is_rhythmic
          if (sum(rhy) == 0) return(NA_real_)
          circ_mse(det$phi_hat[rhy], sim$phi1[rhy])
        }
      }, mc.cores = NCORES, mc.preschedule = FALSE)
      vals <- unlist(vals[vapply(vals, is.numeric, logical(1))])
      M[j, k] <- if (length(vals) > 0) mean(vals, na.rm = TRUE) else NA
    }
    cat(sprintf("  N=%3d %s K=%d: %s\n", N_grid[j], measure, K,
                paste(sprintf("%6.3f", M[j,]), collapse = " ")))
  }
  M
}

cat("\n== Panel A: K=1 power ==\n")
pwr_K1 <- run_grid(DISPLAY_N, B_K1, sim_K1, K = 1L, measure = "power")
cat("\n== Panel B: K=2 power on K=1 sim ==\n")
pwr_K2 <- run_grid(DISPLAY_N, B_K2, sim_K1, K = 2L, measure = "power")
cat("\n== Panel C: K=1 phase MSE ==\n")
mse_K1 <- run_grid(DISPLAY_N, B_K1, sim_K1, K = 1L, measure = "phase")

saveRDS(list(pwr_K1 = pwr_K1, pwr_K2 = pwr_K2, mse_K1 = mse_K1,
              N = DISPLAY_N, B_K1 = B_K1, B_K2 = B_K2,
              source = "Putamen Control"),
        "output/two_harmonic/results/fig6_v3_PutamenCtrl.rds")

cols <- c("#440154","#3b528b","#21918c","#5ec962","#fde725")
pdf("submission/figures/Fig6_v3_PutamenCtrl.pdf", width = 10.5, height = 3.7)
par(mfrow = c(1, 3), mai = c(0.95, 0.85, 0.50, 0.10),
    mgp = c(2.6, 0.55, 0), oma = c(0, 0, 2.0, 0),
    cex.axis = 0.90, cex.lab = 1.05, cex.main = 1.05, font.main = 2)

# Panel A
matplot(DISPLAY_N, 100 * pwr_K1, type = "o", pch = 19, lwd = 1.8,
        col = cols[seq_along(B_K1)], lty = 1, cex = 0.55,
        ylim = c(0, 100), xlim = c(0, max(DISPLAY_N) * 1.05),
        xlab = "Sample size (n)", ylab = "Power (%)", main = "")
title(main = "A   K=1 Power (K=1 simulator)",
      adj = 0.5, font.main = 2, cex.main = 1.05, line = 0.3)
abline(h = 80, lty = 2, col = "grey50"); grid()
legend("bottomright", paste0("B=", B_K1), col = cols[seq_along(B_K1)],
       lty = 1, pch = 19, lwd = 1.4, cex = 0.55,
       bty = "o", box.col = "grey70", inset = 0.01,
       y.intersp = 0.78, bg = "white", seg.len = 1.2)

# Panel B
matplot(DISPLAY_N, 100 * pwr_K2, type = "o", pch = 19, lwd = 1.8,
        col = cols[seq_along(B_K2)], lty = 1, cex = 0.55,
        ylim = c(0, 100), xlim = c(0, max(DISPLAY_N) * 1.05),
        xlab = "Sample size (n)", ylab = "Power (%)", main = "")
title(main = "B   K=2 Power (K=2 simulator)",
      adj = 0.5, font.main = 2, cex.main = 1.05, line = 0.3)
abline(h = 80, lty = 2, col = "grey50"); grid()
legend("bottomright", paste0("B=", B_K2), col = cols[seq_along(B_K2)],
       lty = 1, pch = 19, lwd = 1.4, cex = 0.55,
       bty = "o", box.col = "grey70", inset = 0.01,
       y.intersp = 0.78, bg = "white", seg.len = 1.2)

# Panel C
matplot(DISPLAY_N, mse_K1, type = "o", pch = 19, lwd = 1.8,
        col = cols[seq_along(B_K1)], lty = 1, cex = 0.55,
        xlim = c(0, max(DISPLAY_N) * 1.05),
        ylim = c(0, max(mse_K1, na.rm = TRUE) * 1.05),
        xlab = "Sample size (n)",
        ylab = expression(median ~ "(" * hat(phi) - phi[true] * ")"^2 ~ "(h"^2*")"),
        main = "")
title(main = "C   Phase MSE (K=1 simulator)",
      adj = 0.5, font.main = 2, cex.main = 1.05, line = 0.3)
grid()
legend("topright", paste0("B=", B_K1), col = cols[seq_along(B_K1)],
       lty = 1, pch = 19, lwd = 1.4, cex = 0.55,
       bty = "o", box.col = "grey70", inset = 0.01,
       y.intersp = 0.78, bg = "white", seg.len = 1.2)

mtext("B vs m Trade-off: Putamen Control (n=59, passive-derived pilot)",
      outer = TRUE, side = 3, line = 0.2, font = 2, cex = 1.10)
dev.off()
cat("\nSaved: submission/figures/Fig6_v3_PutamenCtrl.pdf\n")
