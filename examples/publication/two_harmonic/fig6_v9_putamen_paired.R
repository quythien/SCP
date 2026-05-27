#!/usr/bin/env Rscript
# =============================================================================
# Fig 6 v9 - Active B vs m on Putamen Control (n=59, brain)
#
# Panel A: K=1 detector + K=1 simulator + K=1 paired top-K=300 sampling
#   (matched detector + simulator on K=1 truth)
# Panel B: K=2 detector + K=2 simulator + K=2 paired top-K=300 sampling
#   (matched detector + simulator on K=2 truth, with synthetic A2 = 0.5*A1
#    matching GTEx Liver's empirical A2/A1 ratio).
#
# Sampling convention: joint draw of (A, sigma, phase) from top-K=300 paired set,
# matching the package's paired_sigma=TRUE convention.
# =============================================================================
setwd("/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")
sink(stderr())
source("code/utils.R"); source("code/options.R"); source("code/simulation.R")
source("code/estimation.R"); source("code/bootstrap_sim.R"); source("code/detection.R")
sink()
suppressPackageStartupMessages(library(parallel))

bio <- readRDS("data/gse160521_putamen_ctrl_pilot.rds")
cat(sprintf("Pilot: Putamen Control (n=%d, prop_rhy=%.3f)\n",
            length(bio$cts), bio$prop_rhythmic))
cat(sprintf("  amplitude length: %d, sigma length: %d, phase length: %d\n",
            length(bio$amplitude), length(bio$sigma_rhythmic), length(bio$phase)))

# Top-K paired set = first length(sigma_rhythmic) entries of amp/phase paired with sigma
TOP_K <- length(bio$sigma_rhythmic)
A_top   <- bio$amplitude[seq_len(TOP_K)]
sig_top <- bio$sigma_rhythmic
phi_top <- bio$phase[seq_len(TOP_K)]
cat(sprintf("  paired top-K=%d: median r̃=%.2f (A/sigma)\n",
            TOP_K, median(A_top / sig_top)))

# Background distributions for non-rhythmic genes
mu_dist  <- bio$lBaselineExpr
sig_bg   <- exp(bio$lOD)

PERIOD <- 24; OMEGA0 <- 2*pi/PERIOD
NSIMS <- 100; NCORES <- 12; NGENES <- 2000; SEED_BASE <- 20260527
DISPLAY_N <- c(6L, 8L, 10L, 12L, 16L, 20L, 25L, 30L, 40L, 50L, 60L, 80L, 100L)
B_GRID    <- c(3L, 4L, 6L, 8L, 12L, 24L)

# Synthetic K=2 content: A2 = 0.5 * A1 paired per gene (matches GTEx Liver's A2/A1 ~ 0.56)
phi2_top <- (phi_top + runif(TOP_K, -3, 3)) %% PERIOD  # phi2 near phi1 with jitter
A2_top   <- 0.5 * A_top

active_times <- function(N, B) {
  m  <- N %/% B
  ts <- rep(seq(0, PERIOD - PERIOD/B, length.out = B), each = m)
  if (length(ts) < N) ts <- c(ts, sample(seq(0, PERIOD - PERIOD/B, length.out = B),
                                          N - length(ts), replace = TRUE))
  if (length(ts) > N) ts <- ts[seq_len(N)]
  ts
}

# Panel A: K=1 simulator + K=1 detector (paired top-K=300 sampling)
sim_K1 <- function(N, B, seed) {
  if (length(unique(active_times(N, B))) < 3) return(NA_real_)
  set.seed(seed)
  ts <- active_times(N, B)
  n_rhy <- round(NGENES * bio$prop_rhythmic)
  ji <- sample.int(TOP_K, n_rhy, replace = TRUE)   # PAIRED draw
  A <- A_top[ji]; phi <- phi_top[ji]; s_rhy <- sig_top[ji]
  mu <- sample(mu_dist, NGENES, replace = TRUE)
  is_rhy <- c(rep(TRUE, n_rhy), rep(FALSE, NGENES - n_rhy))
  s_full <- c(s_rhy, sample(sig_bg, NGENES - n_rhy, replace = TRUE))
  expr <- matrix(rnorm(NGENES * length(ts), sd = 1) * s_full, NGENES, length(ts))
  for (g in seq_len(n_rhy))
    expr[g, ] <- expr[g, ] + mu[g] + A[g] * cos(OMEGA0 * (ts - phi[g]))
  for (g in (n_rhy+1):NGENES) expr[g, ] <- expr[g, ] + mu[g]
  X1 <- cbind(1, cos(OMEGA0*ts), sin(OMEGA0*ts))
  df_full <- length(ts) - 3L
  if (df_full <= 0) return(NA_real_)
  pv <- numeric(NGENES)
  for (g in seq_len(NGENES)) {
    y <- expr[g, ]; R0 <- sum((y - mean(y))^2)
    if (R0 <= 0) { pv[g] <- 1; next }
    f <- tryCatch(lm.fit(X1, y), error = function(e) NULL)
    if (is.null(f)) { pv[g] <- 1; next }
    R <- sum(f$residuals^2)
    if (R <= 0 || any(is.na(f$coefficients))) { pv[g] <- 1; next }
    pv[g] <- pf(((R0-R)/2) / (R/df_full), 2, df_full, lower.tail = FALSE)
  }
  q <- p.adjust(pv, "BH")
  sum(q <= 0.05 & is_rhy) / max(1, sum(is_rhy))
}

# Panel B: K=2 simulator (with A2 = 0.5*A1) + K=2 detector
sim_K2 <- function(N, B, seed) {
  if (B < 5 || length(unique(active_times(N, B))) < 5) return(NA_real_)
  set.seed(seed)
  ts <- active_times(N, B)
  n_rhy <- round(NGENES * bio$prop_rhythmic)
  ji <- sample.int(TOP_K, n_rhy, replace = TRUE)   # PAIRED draw
  A1 <- A_top[ji]; phi1 <- phi_top[ji]; s_rhy <- sig_top[ji]
  A2 <- A2_top[ji]; phi2 <- phi2_top[ji]
  mu <- sample(mu_dist, NGENES, replace = TRUE)
  is_rhy <- c(rep(TRUE, n_rhy), rep(FALSE, NGENES - n_rhy))
  s_full <- c(s_rhy, sample(sig_bg, NGENES - n_rhy, replace = TRUE))
  expr <- matrix(rnorm(NGENES * length(ts), sd = 1) * s_full, NGENES, length(ts))
  for (g in seq_len(n_rhy))
    expr[g, ] <- expr[g, ] + mu[g] +
                 A1[g] * cos(OMEGA0 * (ts - phi1[g])) +
                 A2[g] * cos(2 * OMEGA0 * (ts - phi2[g]))
  for (g in (n_rhy+1):NGENES) expr[g, ] <- expr[g, ] + mu[g]
  X2 <- cbind(1, cos(OMEGA0*ts), sin(OMEGA0*ts),
              cos(2*OMEGA0*ts), sin(2*OMEGA0*ts))
  df_full <- length(ts) - 5L
  if (df_full <= 0) return(NA_real_)
  pv <- numeric(NGENES)
  for (g in seq_len(NGENES)) {
    y <- expr[g, ]; R0 <- sum((y - mean(y))^2)
    if (R0 <= 0) { pv[g] <- 1; next }
    f <- tryCatch(lm.fit(X2, y), error = function(e) NULL)
    if (is.null(f)) { pv[g] <- 1; next }
    R <- sum(f$residuals^2)
    if (R <= 0 || any(is.na(f$coefficients))) { pv[g] <- 1; next }
    pv[g] <- pf(((R0-R)/4) / (R/df_full), 4, df_full, lower.tail = FALSE)
  }
  q <- p.adjust(pv, "BH")
  sum(q <= 0.05 & is_rhy) / max(1, sum(is_rhy))
}

run_grid <- function(sim_fn, B_grid, label) {
  cat(sprintf("\n== %s ==\n", label))
  M  <- matrix(NA_real_, length(DISPLAY_N), length(B_grid),
                dimnames = list(paste0("N=", DISPLAY_N), paste0("B=", B_grid)))
  SE <- M
  for (k in seq_along(B_grid)) {
    for (j in seq_along(DISPLAY_N)) {
      N <- DISPLAY_N[j]; B <- B_grid[k]
      v <- mclapply(seq_len(NSIMS), function(i)
        sim_fn(N, B, SEED_BASE + (N*1000L + B + i)*7919L),
        mc.cores = NCORES, mc.preschedule = FALSE)
      v <- unlist(v); v <- v[is.finite(v)]
      if (length(v) > 0) {
        M[j, k]  <- mean(v)
        SE[j, k] <- sd(v) / sqrt(length(v))
      }
    }
    cat(sprintf("  B=%2d: %s\n", B_grid[k],
                paste(sprintf("%5.1f%%", 100*M[, k]), collapse = " ")))
  }
  list(M = M, SE = SE)
}

res_K1 <- run_grid(sim_K1, B_GRID, "Panel A: K=1 detector + K=1 simulator")
res_K2 <- run_grid(sim_K2, B_GRID, "Panel B: K=2 detector + K=2 simulator (A2=0.5*A1)")

saveRDS(list(N = DISPLAY_N, B = B_GRID,
              pwr_K1 = res_K1$M, se_K1 = res_K1$SE,
              pwr_K2 = res_K2$M, se_K2 = res_K2$SE),
        "output/two_harmonic/results/fig6_v9_data.rds")

cols <- c("#440154","#3b528b","#21918c","#5ec962","#fde725","#7d028c")[seq_along(B_GRID)]

draw_panel <- function(N, M, SE, B_grid, cols, letter, main_text, jitter_step = 1.2) {
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
draw_panel(DISPLAY_N, res_K1$M, res_K1$SE, B_GRID, cols,
           "A", "K=1 Power, FDR 5%")
draw_panel(DISPLAY_N, res_K2$M, res_K2$SE, B_GRID, cols,
           "B", "K=2 Power, FDR 5%")
mtext("Active design B vs m trade-off (Putamen Control, n=59)",
      outer = TRUE, side = 3, line = 0.5, font = 2, cex = 1.20)
dev.off()
cat("\nSaved: submission/figures/Fig6_active_BvsM.pdf\n")
