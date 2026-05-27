#!/usr/bin/env Rscript
# =============================================================================
# Fig 6 v5: Identifiability cliff + B-invariance overlap (GTEx Liver K=2)
#
# Panel A: K=1 power vs N, B in {2,3,6,12,24}
#   B=2: rank-deficient for K=1 -> failure
#   B=3: identifiable but aliases K=2 content -> visible dip
#   B>=6: B-invariance overlap
#
# Panel B: K=2 power vs N, B in {3,4,5,8,12,24}
#   B=3,4: rank-deficient for K=2 -> failure
#   B=5: at identifiability boundary
#   B>=6: B-invariance overlap
#
# Panel C: K=2 phase MSE vs B at fixed N, r-tilde stratified
#   B=3,4: rank-deficient -> MSE explodes
#   B>=5: flat (information-optimal)
# =============================================================================
setwd("/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")
sink(stderr())
source("code/utils.R"); source("code/options.R"); source("code/simulation.R")
source("code/estimation.R"); source("code/bootstrap_sim.R"); source("code/detection.R")
sink()
suppressPackageStartupMessages(library(parallel))

psi <- readRDS("output/two_harmonic/results/pilot_2h_GTExLiver_topK300.rds")
cat(sprintf("Pilot: GTEx Liver K=2 (A2/A1 med=%.2f, prop_rhythmic=%.3f)\n",
            psi$diagnostics$A2_over_A1_med, psi$prop_rhythmic))

PERIOD <- 24; OMEGA0 <- 2*pi/PERIOD
NSIMS <- 100; NCORES <- 12; NGENES <- 2000; SEED_BASE <- 20260527

DISPLAY_N_AB <- c(20L, 40L, 60L, 80L, 100L, 120L, 150L, 200L, 250L)
B_K1 <- c(2L, 3L, 6L, 12L, 24L)
B_K2 <- c(3L, 4L, 5L, 8L, 12L, 24L)

N_FIXED_C <- 80L
B_C <- c(3L, 4L, 5L, 6L, 8L, 12L, 24L)
R_BINS <- c(0, 0.3, 0.5, Inf)
R_LABELS <- c("low (r<0.3)", "mid (0.3-0.5)", "high (>0.5)")

active_times <- function(N, B) {
  if (N < B) return(sort(sample(seq(0, PERIOD - PERIOD/B, length.out = B), N, replace = FALSE)))
  m <- N %/% B
  ts <- rep(seq(0, PERIOD - PERIOD/B, length.out = B), each = m)
  if (length(ts) < N) ts <- c(ts, sample(seq(0, PERIOD - PERIOD/B, length.out = B), N - length(ts), replace = TRUE))
  if (length(ts) > N) ts <- ts[seq_len(N)]
  ts
}
circ_se <- function(phi_hat, phi_true) {
  d <- (phi_hat - phi_true) %% PERIOD; pmin(d, PERIOD - d)
}

sim_one <- function(N, B, seed) {
  ts <- active_times(N, B); set.seed(seed)
  psi_local <- psi; psi_local$ngenes <- NGENES
  sim <- simCircadianSingleCohort2H(psi_local, ts, seed = seed)
  G <- nrow(sim$expr); n_unique <- length(unique(ts))
  X1 <- cbind(1, cos(OMEGA0*ts), sin(OMEGA0*ts))
  X2 <- cbind(X1, cos(2*OMEGA0*ts), sin(2*OMEGA0*ts))
  df1 <- length(ts) - 3L; df2 <- length(ts) - 5L
  pv1 <- pv2 <- rep(1, G); phi2 <- rep(NA_real_, G)
  k1_ok <- (n_unique >= 3 && df1 > 0)
  k2_ok <- (n_unique >= 5 && df2 > 0)
  if (k1_ok || k2_ok) {
    for (g in seq_len(G)) {
      y <- sim$expr[g, ]; R0 <- sum((y - mean(y))^2); if (R0 <= 0) next
      if (k1_ok) {
        f1 <- tryCatch(lm.fit(X1, y), error = function(e) NULL)
        if (!is.null(f1)) { R1 <- sum(f1$residuals^2)
          if (R1 > 0) pv1[g] <- pf(((R0-R1)/2)/(R1/df1), 2, df1, lower.tail=FALSE)
        }
      }
      if (k2_ok) {
        f2 <- tryCatch(lm.fit(X2, y), error = function(e) NULL)
        if (!is.null(f2) && !any(is.na(f2$coefficients[2:3]))) {
          R2 <- sum(f2$residuals^2)
          if (R2 > 0) pv2[g] <- pf(((R0-R2)/4)/(R2/df2), 4, df2, lower.tail=FALSE)
          phi2[g] <- (atan2(f2$coefficients[3], f2$coefficients[2]) / OMEGA0) %% PERIOD
        }
      }
    }
  }
  q1 <- p.adjust(pv1, "BH"); q2 <- p.adjust(pv2, "BH")
  rhy <- sim$is_rhythmic
  list(power_K1 = if (k1_ok) sum(q1 <= 0.05 & rhy) / max(1, sum(rhy)) else NA,
       power_K2 = if (k2_ok) sum(q2 <= 0.05 & rhy) / max(1, sum(rhy)) else NA,
       errs_K2  = if (k2_ok) data.frame(err2 = circ_se(phi2[rhy], sim$phi1[rhy])^2,
                                          r = sim$r_values[rhy]) else NULL)
}

# ----- Panel A: K=1 power -----
cat("\n== Panel A: K=1 power ==\n")
pwr_K1 <- matrix(NA_real_, length(DISPLAY_N_AB), length(B_K1),
                  dimnames = list(paste0("N=",DISPLAY_N_AB), paste0("B=",B_K1)))
for (j in seq_along(DISPLAY_N_AB)) {
  for (k in seq_along(B_K1)) {
    N <- DISPLAY_N_AB[j]; B <- B_K1[k]; if (N < B) next
    res <- mclapply(seq_len(NSIMS), function(i) sim_one(N, B, SEED_BASE + (N*1000L+B+i)*7919L),
                    mc.cores = NCORES, mc.preschedule = FALSE)
    res <- res[vapply(res, is.list, logical(1))]
    if (length(res) > 0) {
      v <- sapply(res, function(x) x$power_K1)
      pwr_K1[j, k] <- mean(v, na.rm = TRUE)
    }
  }
  cat(sprintf("  N=%3d: %s\n", DISPLAY_N_AB[j], paste(sprintf("%5.1f%%", 100*pwr_K1[j,]), collapse=" ")))
}

# ----- Panel B: K=2 power -----
cat("\n== Panel B: K=2 power ==\n")
pwr_K2 <- matrix(NA_real_, length(DISPLAY_N_AB), length(B_K2),
                  dimnames = list(paste0("N=",DISPLAY_N_AB), paste0("B=",B_K2)))
for (j in seq_along(DISPLAY_N_AB)) {
  for (k in seq_along(B_K2)) {
    N <- DISPLAY_N_AB[j]; B <- B_K2[k]; if (N < B) next
    res <- mclapply(seq_len(NSIMS), function(i) sim_one(N, B, SEED_BASE + (N*1000L+B+i)*7919L),
                    mc.cores = NCORES, mc.preschedule = FALSE)
    res <- res[vapply(res, is.list, logical(1))]
    if (length(res) > 0) {
      v <- sapply(res, function(x) x$power_K2)
      pwr_K2[j, k] <- mean(v, na.rm = TRUE)
    }
  }
  cat(sprintf("  N=%3d: %s\n", DISPLAY_N_AB[j], paste(sprintf("%5.1f%%", 100*pwr_K2[j,]), collapse=" ")))
}

# ----- Panel C: K=2 phase MSE vs B at fixed N -----
cat(sprintf("\n== Panel C: K=2 phase MSE vs B (N=%d) ==\n", N_FIXED_C))
M_phase <- matrix(NA_real_, length(B_C), length(R_LABELS),
                   dimnames = list(paste0("B=",B_C), R_LABELS))
for (k in seq_along(B_C)) {
  B <- B_C[k]
  res <- mclapply(seq_len(NSIMS * 2), function(i)
    sim_one(N_FIXED_C, B, SEED_BASE + (B*1000L+i)*7919L),
    mc.cores = NCORES, mc.preschedule = FALSE)
  res <- res[vapply(res, is.list, logical(1))]
  if (length(res) == 0) next
  errs <- do.call(rbind, lapply(res, function(x) x$errs_K2))
  if (is.null(errs) || nrow(errs) == 0) next
  errs <- errs[is.finite(errs$err2) & is.finite(errs$r) & errs$r > 0, ]
  bin_idx <- cut(errs$r, breaks = R_BINS, include.lowest = TRUE, labels = FALSE)
  for (b in seq_along(R_LABELS)) {
    sel <- bin_idx == b
    if (sum(sel, na.rm = TRUE) > 5)
      M_phase[k, b] <- median(errs$err2[sel], na.rm = TRUE)
  }
  cat(sprintf("  B=%2d: %s\n", B, paste(sprintf("%6.2f", M_phase[k,]), collapse=" ")))
}

saveRDS(list(pwr_K1 = pwr_K1, pwr_K2 = pwr_K2, M_phase = M_phase,
              N_AB = DISPLAY_N_AB, B_K1 = B_K1, B_K2 = B_K2, B_C = B_C,
              N_FIXED_C = N_FIXED_C, R_LABELS = R_LABELS),
        "output/two_harmonic/results/fig6_v5_data.rds")

cols_A <- c("#fde725","#5ec962","#21918c","#3b528b","#440154")[seq_along(B_K1)]
cols_B <- c("#fde725","#5ec962","#21918c","#3b528b","#440154","#7d028c")[seq_along(B_K2)]
cols_C <- c("#440154","#21918c","#fde725")[seq_along(R_LABELS)]

pdf("submission/figures/Fig6_v5_cliff_overlap.pdf", width = 11.0, height = 3.8)
par(mfrow = c(1, 3), mai = c(0.95, 0.85, 0.55, 0.10),
    mgp = c(2.6, 0.55, 0), oma = c(0, 0, 2.0, 0),
    cex.axis = 0.90, cex.lab = 1.05, cex.main = 1.05, font.main = 2)

# Panel A
matplot(DISPLAY_N_AB, 100 * pwr_K1, type = "o", pch = 19, lwd = 1.8,
        col = cols_A, lty = 1, cex = 0.55, ylim = c(0, 100),
        xlim = c(0, max(DISPLAY_N_AB) * 1.05),
        xlab = "Sample size (n)", ylab = "Power (%)", main = "")
title(main = "A   K=1 Power", adj = 0.5, font.main = 2, cex.main = 1.05, line = 0.3)
abline(h = 80, lty = 2, col = "grey50"); grid()
legend("bottomright", paste0("B=", B_K1), col = cols_A,
       lty = 1, pch = 19, lwd = 1.4, cex = 0.55,
       bty = "o", box.col = "grey70", inset = 0.01, y.intersp = 0.78, bg = "white")

# Panel B
matplot(DISPLAY_N_AB, 100 * pwr_K2, type = "o", pch = 19, lwd = 1.8,
        col = cols_B, lty = 1, cex = 0.55, ylim = c(0, 100),
        xlim = c(0, max(DISPLAY_N_AB) * 1.05),
        xlab = "Sample size (n)", ylab = "Power (%)", main = "")
title(main = "B   K=2 Power", adj = 0.5, font.main = 2, cex.main = 1.05, line = 0.3)
abline(h = 80, lty = 2, col = "grey50"); grid()
legend("bottomright", paste0("B=", B_K2), col = cols_B,
       lty = 1, pch = 19, lwd = 1.4, cex = 0.55,
       bty = "o", box.col = "grey70", inset = 0.01, y.intersp = 0.78, bg = "white")

# Panel C
ymax <- max(M_phase, na.rm = TRUE) * 1.05
matplot(B_C, M_phase, type = "o", pch = 19, lwd = 1.8,
        col = cols_C, lty = 1, cex = 0.6,
        ylim = c(0, ymax), xlim = c(2, 26),
        xlab = "B (number of distinct phases)",
        ylab = expression(median ~ (hat(phi)[1] - phi[1*","*true])^2 ~ "(h"^2*")"),
        main = "")
title(main = sprintf("C   Phase MSE vs B (N=%d)", N_FIXED_C),
      adj = 0.5, font.main = 2, cex.main = 1.05, line = 0.3)
abline(v = 5, lty = 3, col = "grey60", lwd = 1.2)
grid()
legend("topright", R_LABELS, col = cols_C, lty = 1, pch = 19, lwd = 1.6,
       cex = 0.65, bty = "o", box.col = "grey70", inset = 0.01,
       y.intersp = 0.85, bg = "white")

mtext("Identifiability cliff and B-invariance above threshold (GTEx Liver K=2)",
      outer = TRUE, side = 3, line = 0.2, font = 2, cex = 1.10)
dev.off()
cat("\nSaved: submission/figures/Fig6_v5_cliff_overlap.pdf\n")
