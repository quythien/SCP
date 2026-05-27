#!/usr/bin/env Rscript
# =============================================================================
# Fig 6 v7 - Active B vs m, native pilots per detector
#   Panel A: K=1 detector + K=1 simulator + K=1 paired pilot (estCircadianParam)
#   Panel B: K=2 detector + K=2 simulator + K=2 paired pilot (estCircadianParam2H)
# Jittered curves, +/- 1 SE error bars, B values include sub-threshold to show cliff.
# =============================================================================
setwd("/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")
sink(stderr())
source("code/utils.R"); source("code/options.R"); source("code/simulation.R")
source("code/estimation.R"); source("code/bootstrap_sim.R"); source("code/detection.R")
sink()
suppressPackageStartupMessages(library(parallel))

# Pilots
bio_K1 <- readRDS("data/gtex_Liver_single_pilot.rds")
psi_K2 <- readRDS("output/two_harmonic/results/pilot_2h_GTExLiver_topK300.rds")
cat(sprintf("K=1 pilot: n=%d, prop_rhythmic=%.3f, n_amp=%d, n_sigma=%d\n",
            length(bio_K1$cts), bio_K1$prop_rhythmic,
            length(bio_K1$amplitude), length(bio_K1$sigma_rhythmic)))
cat(sprintf("K=2 pilot: n=%d, prop_rhythmic=%.3f, A2/A1 med=%.2f\n",
            attr(psi_K2, "n_pilot"), psi_K2$prop_rhythmic,
            psi_K2$diagnostics$A2_over_A1_med))

PERIOD <- 24; OMEGA0 <- 2*pi/PERIOD
NSIMS <- 100; NCORES <- 12; NGENES <- 2000; SEED_BASE <- 20260527
DISPLAY_N <- c(20L, 40L, 60L, 80L, 100L, 120L, 150L, 200L, 250L)
B_K1 <- c(2L, 3L, 6L, 8L, 12L, 24L)
B_K2 <- c(3L, 4L, 6L, 8L, 12L, 24L)

active_times <- function(N, B) {
  if (N < B) return(sort(sample(seq(0, PERIOD - PERIOD/B, length.out = B), N, replace = FALSE)))
  m <- N %/% B
  ts <- rep(seq(0, PERIOD - PERIOD/B, length.out = B), each = m)
  if (length(ts) < N) ts <- c(ts, sample(seq(0, PERIOD - PERIOD/B, length.out = B), N - length(ts), replace = TRUE))
  if (length(ts) > N) ts <- ts[seq_len(N)]
  ts
}

# K=1 simulator: cosinor from K=1 paired pilot, K=1 detection
sim_K1 <- function(N, B, seed) {
  ts <- active_times(N, B); set.seed(seed)
  n_rhy <- round(NGENES * bio_K1$prop_rhythmic)
  ji <- sample.int(length(bio_K1$amplitude), n_rhy, replace = TRUE)
  A   <- bio_K1$amplitude[ji]
  phi <- sample(bio_K1$phase, n_rhy, replace = TRUE)
  s_rhy <- sample(bio_K1$sigma_rhythmic, n_rhy, replace = TRUE)
  mu <- sample(bio_K1$lBaselineExpr, NGENES, replace = TRUE)
  is_rhy <- c(rep(TRUE, n_rhy), rep(FALSE, NGENES - n_rhy))
  s_full <- c(s_rhy, sample(exp(bio_K1$lOD), NGENES - n_rhy, replace = TRUE))
  expr <- matrix(rnorm(NGENES * length(ts), sd = 1) * s_full, NGENES, length(ts))
  for (g in seq_len(n_rhy))
    expr[g, ] <- expr[g, ] + mu[g] + A[g] * cos(OMEGA0 * (ts - phi[g]))
  for (g in (n_rhy + 1):NGENES) expr[g, ] <- expr[g, ] + mu[g]

  X1 <- cbind(1, cos(OMEGA0*ts), sin(OMEGA0*ts))
  df_full <- length(ts) - 3L
  if (df_full <= 0 || length(unique(ts)) < 3) return(NA_real_)
  pv <- numeric(NGENES)
  for (g in seq_len(NGENES)) {
    y <- expr[g, ]; R0 <- sum((y - mean(y))^2)
    if (R0 <= 0) { pv[g] <- 1; next }
    f <- tryCatch(lm.fit(X1, y), error = function(e) NULL)
    if (is.null(f)) { pv[g] <- 1; next }
    R <- sum(f$residuals^2)
    if (R <= 0 || any(is.na(f$coefficients))) { pv[g] <- 1; next }
    pv[g] <- pf(((R0 - R)/2) / (R/df_full), 2, df_full, lower.tail = FALSE)
  }
  q <- p.adjust(pv, "BH")
  sum(q <= 0.05 & is_rhy) / max(1, sum(is_rhy))
}

# K=2 simulator + K=2 detector, native K=2 pilot
sim_K2 <- function(N, B, seed) {
  ts <- active_times(N, B); set.seed(seed)
  psi_local <- psi_K2; psi_local$ngenes <- NGENES
  sim <- simCircadianSingleCohort2H(psi_local, ts, seed = seed)
  X2 <- cbind(1, cos(OMEGA0*ts), sin(OMEGA0*ts),
              cos(2*OMEGA0*ts), sin(2*OMEGA0*ts))
  df_full <- length(ts) - 5L
  if (df_full <= 0 || length(unique(ts)) < 5) return(NA_real_)
  G <- nrow(sim$expr); pv <- numeric(G)
  for (g in seq_len(G)) {
    y <- sim$expr[g, ]; R0 <- sum((y - mean(y))^2)
    if (R0 <= 0) { pv[g] <- 1; next }
    f <- tryCatch(lm.fit(X2, y), error = function(e) NULL)
    if (is.null(f)) { pv[g] <- 1; next }
    R <- sum(f$residuals^2)
    if (R <= 0 || any(is.na(f$coefficients))) { pv[g] <- 1; next }
    pv[g] <- pf(((R0 - R)/4) / (R/df_full), 4, df_full, lower.tail = FALSE)
  }
  q <- p.adjust(pv, "BH"); rhy <- sim$is_rhythmic
  sum(q <= 0.05 & rhy) / max(1, sum(rhy))
}

run_grid <- function(N_grid, B_grid, sim_fn) {
  M  <- matrix(NA_real_, length(N_grid), length(B_grid),
                dimnames = list(paste0("N=",N_grid), paste0("B=",B_grid)))
  SE <- M
  for (j in seq_along(N_grid)) {
    for (k in seq_along(B_grid)) {
      N <- N_grid[j]; B <- B_grid[k]; if (N < B) next
      v <- mclapply(seq_len(NSIMS), function(i)
        sim_fn(N, B, SEED_BASE + (N*1000L + B + i)*7919L),
        mc.cores = NCORES, mc.preschedule = FALSE)
      v <- unlist(v); v <- v[is.finite(v)]
      if (length(v) > 0) { M[j,k] <- mean(v); SE[j,k] <- sd(v) / sqrt(length(v)) }
    }
    cat(sprintf("  N=%3d: %s\n", N_grid[j],
                paste(sprintf("%5.1f%%", 100*M[j,]), collapse=" ")))
  }
  list(M = M, SE = SE)
}

cat("\n== Panel A: K=1 detector + K=1 paired pilot ==\n")
res_K1 <- run_grid(DISPLAY_N, B_K1, sim_K1)
cat("\n== Panel B: K=2 detector + K=2 paired pilot ==\n")
res_K2 <- run_grid(DISPLAY_N, B_K2, sim_K2)

saveRDS(list(N = DISPLAY_N, B_K1 = B_K1, B_K2 = B_K2,
              pwr_K1 = res_K1$M, se_K1 = res_K1$SE,
              pwr_K2 = res_K2$M, se_K2 = res_K2$SE),
        "output/two_harmonic/results/fig6_v7_data.rds")

cols_K1 <- c("#7d028c","#440154","#3b528b","#21918c","#5ec962","#fde725")
cols_K2 <- c("#7d028c","#440154","#3b528b","#21918c","#5ec962","#fde725")

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
    M_k <- M[, k]; M_k[!is.finite(M_k)] <- 0  # rank-deficient B -> 0%
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
draw_panel(DISPLAY_N, res_K1$M, res_K1$SE, B_K1, cols_K1,
           "A", "K=1 Power, FDR 5% (K=1 pilot)")
draw_panel(DISPLAY_N, res_K2$M, res_K2$SE, B_K2, cols_K2,
           "B", "K=2 Power, FDR 5% (K=2 pilot)")
mtext("Active design B vs m trade-off (GTEx Liver pilot)",
      outer = TRUE, side = 3, line = 0.5, font = 2, cex = 1.20)
dev.off()
cat("\nSaved: submission/figures/Fig6_active_BvsM.pdf\n")
