#!/usr/bin/env Rscript
# =============================================================================
# Fig 6 - Active design + phase estimation on Hughes 2009 mouse liver
#
# 1x3 layout (matching prior FMM-era Fig 6 design):
#   A: K=1 (cosinor) power vs N, B in {4, 6, 8, 12, 24}
#   B: Phase estimation - Median phase MSE of phi_1_hat under K=1 cosinor fit
#      vs N, B in {4, 6, 8, 12, 24}, with SE vertical bars
#   C: K=2 power vs N, B in {6, 8, 12, 24} (B>=5 identifiability bound)
#
# Pilot: Hughes 2009 mouse liver, active B=24, N=48.
# Outputs:
#   - output/two_harmonic/results/fig6_sims_Hughes.rds  (per-replicate data)
#   - submission/figures/Fig6_twoharm_active.pdf
# =============================================================================

RUN_SIMS <- TRUE

setwd("/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")
sink(stderr())
source("code/utils.R"); source("code/options.R"); source("code/simulation.R")
source("code/estimation.R"); source("code/bootstrap_sim.R"); source("code/detection.R")
sink()
suppressPackageStartupMessages(library(parallel))

PILOT_PATH <- "output/two_harmonic/results/pilot_2h_Hughes2009.rds"
SIMS_PATH  <- "output/two_harmonic/results/fig6_sims_Hughes.rds"
FIG_PATH   <- "submission/figures/Fig6_twoharm_active.pdf"

PERIOD       <- 24
OMEGA0       <- 2 * pi / PERIOD
DISPLAY_N    <- c(20L, 40L, 60L, 80L, 100L)
B_GRID_K1    <- c(4L, 6L, 8L, 12L, 24L)
B_GRID_K2    <- c(6L, 8L, 12L, 24L)
B_GRID_PHASE <- c(4L, 6L, 8L, 12L, 24L)
NSIMS        <- 100
NGENES_SIM   <- 2000
NCORES       <- min(12, parallel::detectCores() - 4)
SEED_BASE    <- 20260526

psi <- readRDS(PILOT_PATH)
cat(sprintf("Pilot: %s   n=%d   top_k=%d   A2/A1 med=%.2f   prop_rhythmic=%.3f\n",
            attr(psi, "pilot_label"),
            attr(psi, "n_pilot"),
            psi$diagnostics$top_k_used,
            psi$diagnostics$A2_over_A1_med,
            psi$prop_rhythmic))

# Active-design time grid
active_times <- function(N, B) {
  if (N < B)
    return(sort(sample(seq(0, PERIOD - PERIOD/B, length.out = B), N,
                        replace = FALSE)))
  m  <- N %/% B
  ts <- rep(seq(0, PERIOD - PERIOD/B, length.out = B), each = m)
  if (length(ts) < N) ts <- c(ts, sample(seq(0, PERIOD - PERIOD/B,
                                              length.out = B),
                                          N - length(ts), replace = TRUE))
  if (length(ts) > N) ts <- ts[seq_len(N)]
  ts
}

# Circular MSE of phi_hat against phi_true (both in hours, period 24)
circ_mse <- function(phi_hat, phi_true) {
  d <- (phi_hat - phi_true) %% PERIOD
  d <- pmin(d, PERIOD - d)
  median(d^2, na.rm = TRUE)
}

# One simulated replicate -> per-gene p-values + phase estimates
sim_one <- function(N, B, seed, psi, do_K1 = TRUE, do_K2 = TRUE) {
  ts <- active_times(N, B)
  psi_local <- psi; psi_local$ngenes <- NGENES_SIM
  sim <- simCircadianSingleCohort2H(psi_local, ts, seed = seed)
  G <- nrow(sim$expr)

  out <- list(is_rhythmic = sim$is_rhythmic,
              phi1_true   = sim$phi1)

  if (do_K1 && length(unique(ts)) >= 3L) {
    X1  <- cbind(1, cos(OMEGA0 * ts), sin(OMEGA0 * ts))
    df1 <- length(ts) - 3L
    p1  <- phi1_hat_K1 <- numeric(G)
    for (g in seq_len(G)) {
      y <- sim$expr[g, ]
      R0 <- sum((y - mean(y))^2); if (R0 <= 0) { p1[g] <- 1; next }
      f <- lm.fit(X1, y)
      R <- sum(f$residuals^2); if (R <= 0) { p1[g] <- 1; next }
      p1[g] <- pf(((R0 - R)/2) / (R / df1), 2, df1, lower.tail = FALSE)
      phi1_hat_K1[g] <- (atan2(f$coefficients[3], f$coefficients[2]) / OMEGA0) %% PERIOD
    }
    out$p_K1        <- p1
    out$phi1_hat_K1 <- phi1_hat_K1
  }

  if (do_K2 && length(unique(ts)) >= 5L) {
    X2  <- cbind(1, cos(OMEGA0*ts), sin(OMEGA0*ts),
                 cos(2*OMEGA0*ts), sin(2*OMEGA0*ts))
    df2 <- length(ts) - 5L
    p2 <- numeric(G)
    for (g in seq_len(G)) {
      y <- sim$expr[g, ]
      R0 <- sum((y - mean(y))^2); if (R0 <= 0) { p2[g] <- 1; next }
      f <- lm.fit(X2, y)
      R <- sum(f$residuals^2); if (R <= 0) { p2[g] <- 1; next }
      p2[g] <- pf(((R0 - R)/4) / (R / df2), 4, df2, lower.tail = FALSE)
    }
    out$p_K2 <- p2
  }
  out
}

# ----- Run sims over full (N, B) grid -----
if (RUN_SIMS || !file.exists(SIMS_PATH)) {
  cat(sprintf("Running sims: %d N values x %d B values x %d sims x %d cores\n",
              length(DISPLAY_N), length(B_GRID_PHASE), NSIMS, NCORES))
  t0 <- Sys.time()
  B_union <- sort(unique(c(B_GRID_K1, B_GRID_K2, B_GRID_PHASE)))

  # Stats per (N, B): K=1 power, K=2 power, K=1 phase MSE (with per-rep arrays)
  power_K1 <- power_K2 <- mse_K1 <- matrix(NA_real_,
                                            nrow = length(DISPLAY_N),
                                            ncol = length(B_union),
                                            dimnames = list(N = DISPLAY_N, B = B_union))
  mse_K1_se <- mse_K1_sd <- mse_K1
  power_K1_se <- power_K2_se <- mse_K1

  for (b_idx in seq_along(B_union)) {
    B <- B_union[b_idx]
    do_K2 <- B >= 5L
    for (n_idx in seq_along(DISPLAY_N)) {
      N <- DISPLAY_N[n_idx]
      sims <- mclapply(seq_len(NSIMS), function(i) {
        tryCatch(sim_one(N, B,
                          seed = (SEED_BASE + (B*1000L + N + i)*7919L) %% 2147483647L,
                          psi  = psi, do_K1 = TRUE, do_K2 = do_K2),
                 error = function(e) list(err = conditionMessage(e)))
      }, mc.cores = NCORES, mc.preschedule = FALSE)
      sims <- sims[sapply(sims, function(s) is.list(s) && is.null(s$err))]
      if (length(sims) == 0) next

      # K=1 power per replicate
      p1_rep <- sapply(sims, function(s) {
        if (is.null(s$p_K1)) return(NA_real_)
        q <- p.adjust(s$p_K1, "BH")
        sum(q <= 0.05 & s$is_rhythmic) / max(1, sum(s$is_rhythmic))
      })
      power_K1[n_idx, b_idx]    <- mean(p1_rep, na.rm = TRUE)
      power_K1_se[n_idx, b_idx] <- sd(p1_rep, na.rm = TRUE) / sqrt(sum(!is.na(p1_rep)))

      # K=2 power per replicate (only if do_K2)
      if (do_K2) {
        p2_rep <- sapply(sims, function(s) {
          if (is.null(s$p_K2)) return(NA_real_)
          q <- p.adjust(s$p_K2, "BH")
          sum(q <= 0.05 & s$is_rhythmic) / max(1, sum(s$is_rhythmic))
        })
        power_K2[n_idx, b_idx]    <- mean(p2_rep, na.rm = TRUE)
        power_K2_se[n_idx, b_idx] <- sd(p2_rep, na.rm = TRUE) / sqrt(sum(!is.na(p2_rep)))
      }

      # K=1 phase MSE per replicate (across rhythmic genes)
      mse_rep <- sapply(sims, function(s) {
        if (is.null(s$phi1_hat_K1)) return(NA_real_)
        is_R <- s$is_rhythmic & s$phi1_true > 0
        if (!any(is_R)) return(NA_real_)
        circ_mse(s$phi1_hat_K1[is_R], s$phi1_true[is_R])
      })
      mse_K1[n_idx, b_idx]    <- mean(mse_rep, na.rm = TRUE)
      mse_K1_se[n_idx, b_idx] <- sd(mse_rep, na.rm = TRUE) / sqrt(sum(!is.na(mse_rep)))
    }
    cat(sprintf("  B=%d done at %.0f s\n", B,
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
  saveRDS(list(N = DISPLAY_N, B = B_union,
               power_K1 = power_K1, power_K1_se = power_K1_se,
               power_K2 = power_K2, power_K2_se = power_K2_se,
               mse_K1   = mse_K1,   mse_K1_se   = mse_K1_se,
               B_K1 = B_GRID_K1, B_K2 = B_GRID_K2, B_phase = B_GRID_PHASE),
          SIMS_PATH)
}

res <- readRDS(SIMS_PATH)

# ----- Render -----
pdf(FIG_PATH, width = 13.5, height = 4.6)
par(mfrow = c(1, 3), mai = c(0.95, 1.05, 0.80, 0.20),
    mgp = c(3.0, 0.7, 0), oma = c(0, 0, 1.6, 0),
    cex.axis = 1.20, cex.lab = 1.40, font.main = 2)

# Viridis-style palette for B (matches old Fig 6 design)
viridis_pal <- function(n) {
  cols <- c("#440154", "#3b528b", "#21918c", "#5ec962", "#fde725")
  if (n <= length(cols)) cols[seq_len(n)]
  else colorRampPalette(cols)(n)
}

# ----- Panel A: K=1 power vs N at B in B_GRID_K1 -----
b_idx_K1 <- match(B_GRID_K1, res$B)
pwr_K1   <- 100 * res$power_K1[, b_idx_K1]
se_K1    <- 100 * res$power_K1_se[, b_idx_K1]
cols_K1  <- viridis_pal(length(B_GRID_K1))

matplot(res$N, pwr_K1, type = "b", pch = 19, lwd = 2,
        col = cols_K1, lty = 1, cex = 0.9,
        xlim = c(0, max(res$N) * 1.05), ylim = c(0, 100),
        xlab = "N (total samples)", ylab = "Power (K=1, FDR = 0.05)",
        main = "")
title(main = "Active KIM,  cosinor F-test",
      adj = 0.5, font.main = 2, cex.main = 1.30, line = 1.0)
mtext("A", side = 3, line = 1.6, at = par("usr")[1],
      adj = 0, font = 2, cex = 1.20, col = "grey20")
for (k in seq_along(B_GRID_K1)) {
  arrows(res$N, pwr_K1[, k] - 1.96 * se_K1[, k],
         res$N, pwr_K1[, k] + 1.96 * se_K1[, k],
         code = 3, angle = 90, length = 0.04, col = cols_K1[k], lwd = 1.3)
}
abline(h = 80, lty = 3, col = "grey60", lwd = 1)
grid()
legend("bottomright", paste0("B = ", B_GRID_K1),
       col = cols_K1, lty = 1, pch = 19, lwd = 2,
       cex = 0.80, bty = "o", box.col = "grey70", box.lwd = 0.5,
       inset = 0.02, y.intersp = 0.90, bg = "white")

# ----- Panel B: Phase MSE vs N at B in B_GRID_PHASE (K=1 only) -----
b_idx_ph <- match(B_GRID_PHASE, res$B)
mse_med  <- res$mse_K1[, b_idx_ph]
mse_se   <- res$mse_K1_se[, b_idx_ph]
cols_ph  <- viridis_pal(length(B_GRID_PHASE))

ylim_mse <- c(0, max(mse_med + 1.96 * mse_se, na.rm = TRUE) * 1.05)
matplot(res$N, mse_med, type = "b", pch = 19, lwd = 2,
        col = cols_ph, lty = 1, cex = 0.9,
        xlim = c(0, max(res$N) * 1.05), ylim = ylim_mse,
        xlab = "N (total samples)",
        ylab = expression("Median phase MSE (h"^2 * ")"),
        main = "")
title(main = expression("Phase estimation  (" * hat(phi)[g]^cos * ")"),
      adj = 0.5, font.main = 2, cex.main = 1.30, line = 1.0)
mtext("B", side = 3, line = 1.6, at = par("usr")[1],
      adj = 0, font = 2, cex = 1.20, col = "grey20")
for (k in seq_along(B_GRID_PHASE)) {
  arrows(res$N, mse_med[, k] - 1.96 * mse_se[, k],
         res$N, mse_med[, k] + 1.96 * mse_se[, k],
         code = 3, angle = 90, length = 0.04, col = cols_ph[k], lwd = 1.3)
}
grid()
legend("topright", paste0("B = ", B_GRID_PHASE),
       col = cols_ph, lty = 1, pch = 19, lwd = 2,
       cex = 0.80, bty = "o", box.col = "grey70", box.lwd = 0.5,
       inset = 0.02, y.intersp = 0.90, bg = "white")

# ----- Panel C: K=2 power vs N at B in B_GRID_K2 -----
b_idx_K2 <- match(B_GRID_K2, res$B)
pwr_K2   <- 100 * res$power_K2[, b_idx_K2]
se_K2_   <- 100 * res$power_K2_se[, b_idx_K2]
cols_K2  <- viridis_pal(length(B_GRID_K2))

matplot(res$N, pwr_K2, type = "b", pch = 19, lwd = 2,
        col = cols_K2, lty = 1, cex = 0.9,
        xlim = c(0, max(res$N) * 1.05), ylim = c(0, 100),
        xlab = "N (total samples)", ylab = "Power (K=2, FDR = 0.05)",
        main = "")
title(main = "Active KIM,  K = 2 harmonic F-test",
      adj = 0.5, font.main = 2, cex.main = 1.30, line = 1.0)
mtext("C", side = 3, line = 1.6, at = par("usr")[1],
      adj = 0, font = 2, cex = 1.20, col = "grey20")
for (k in seq_along(B_GRID_K2)) {
  arrows(res$N, pwr_K2[, k] - 1.96 * se_K2_[, k],
         res$N, pwr_K2[, k] + 1.96 * se_K2_[, k],
         code = 3, angle = 90, length = 0.04, col = cols_K2[k], lwd = 1.3)
}
abline(h = 80, lty = 3, col = "grey60", lwd = 1)
grid()
legend("bottomright", paste0("B = ", B_GRID_K2),
       col = cols_K2, lty = 1, pch = 19, lwd = 2,
       cex = 0.80, bty = "o", box.col = "grey70", box.lwd = 0.5,
       inset = 0.02, y.intersp = 0.90, bg = "white")

mtext("Active-design properties on Hughes 2009 mouse liver",
      outer = TRUE, side = 3, line = 0.2, font = 2, cex = 1.25)
dev.off()
cat(sprintf("\nWrote %s (%.1f KB)\n", FIG_PATH, file.size(FIG_PATH)/1024))
