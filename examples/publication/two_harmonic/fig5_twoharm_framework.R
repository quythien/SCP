#!/usr/bin/env Rscript
# =============================================================================
# Fig 5 - Two-harmonic framework operating characteristics (GTEx Liver)
#
# 1x2 layout:
#   A: K=2 marginal power vs N, multiple FDR-target lines
#      (alpha in {0.01, 0.05, 0.10, 0.20})
#   B: K=2 power stratified by first-harmonic effect size r1 = A1 / sigma,
#      one curve per N in DISPLAY_N, at BH-FDR <= 0.05
#
# Pilot: GTEx Liver via pilot_2h_GTExLiver.rds
# Outputs:
#   - output/two_harmonic/results/fig5_panelAB_GTExLiver.rds
#   - submission/figures/Fig5_twoharm_framework.pdf
# =============================================================================

RUN_SIMS <- FALSE

setwd("/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")
# Only the exported simCircadianSingleCohort2H is needed (in the sim path, which
# is skipped when replotting from cache), so load the installed package rather
# than the removed code/*.R files.
suppressPackageStartupMessages({ library(SCP); library(parallel) })

PILOT_PATH <- "output/two_harmonic/results/pilot_2h_GTExLiver_topK300.rds"
RES_PATH   <- "output/two_harmonic/results/fig5_panelAB_GTExLiver_topK300_unpaired.rds"
RES_PATH_A <- "output/two_harmonic/results/fig5_panelA_GTExLiver_topK300_paired.rds"
FIG_PATH   <- "submission/figures/Fig5_twoharm_framework.pdf"

PERIOD     <- 24
DISPLAY_N  <- c(20, 40, 60, 80, 100, 120, 150, 200, 250, 300)
ALPHA_GRID <- c(0.01, 0.05, 0.10, 0.20)
NSIMS      <- 100
NGENES_SIM <- 2000
NCORES     <- min(12, parallel::detectCores() - 4)
# 0.25-step bins; cap at r_max = 3
R_BREAKS   <- c(seq(0, 3, by = 0.25), Inf)
R_LABELS   <- c(sprintf("(%g,%g]", head(R_BREAKS, -2), R_BREAKS[seq(2, length(R_BREAKS)-1)]),
                ">3")
SEED_BASE  <- 20260526
VLINE_FDR  <- 0.05
VLINE_POWER <- 0.80

# ---- Load pilot ----
psi <- readRDS(PILOT_PATH)
cat(sprintf("Pilot: n=%d, top_k=%d, A2/A1 med=%.2f, prop_rhythmic=%.3f\n",
            attr(psi, "n_pilot"),
            psi$diagnostics$top_k_used,
            psi$diagnostics$A2_over_A1_med,
            psi$prop_rhythmic))

# ---- Simulator: one K=2 replicate -> p-values + r1 per gene ----
sim_one_rep <- function(N, seed, psi, ngenes = NGENES_SIM,
                         paired_sigma = TRUE) {
  set.seed(seed)
  B <- min(24L, N)
  m <- max(1L, N %/% B)
  ts <- rep(seq(0, PERIOD - PERIOD/B, length.out = B), each = m)
  if (length(ts) > N) ts <- ts[seq_len(N)]
  if (length(ts) < N) ts <- c(ts, sample(ts, N - length(ts), replace = TRUE))

  psi_local <- psi; psi_local$ngenes <- ngenes

  # Optional UNPAIRED sigma: break the (A1,A2,phi1,phi2)<->sigma pairing by
  # independently resampling sigma from the empirical pilot sigma vector.
  # Matches the visualization-friendly behavior of older estCircadianParam
  # pilots (Fig 1A/1B) that used unpaired sigma resampling.
  if (!paired_sigma) {
    n_pilot <- length(psi_local$amplitude)
    ji_sigma <- sample.int(n_pilot, n_pilot, replace = TRUE)
    psi_local$sigma_rhythmic <- psi_local$sigma_rhythmic[ji_sigma]
  }

  sim <- simCircadianSingleCohort2H(psi_local, ts, seed = seed)

  Xs <- cbind(1,
              cos(2*pi/PERIOD*ts), sin(2*pi/PERIOD*ts),
              cos(2*2*pi/PERIOD*ts), sin(2*2*pi/PERIOD*ts))
  df_full <- length(ts) - 5L
  G <- nrow(sim$expr)
  pv <- numeric(G)
  for (g in seq_len(G)) {
    y <- sim$expr[g, ]
    RSS0 <- sum((y - mean(y))^2)
    if (RSS0 <= 0) { pv[g] <- 1; next }
    f <- lm.fit(Xs, y)
    R <- sum(f$residuals^2)
    if (R <= 0) { pv[g] <- 1; next }
    pv[g] <- pf(((RSS0 - R)/4) / (R/df_full), 4, df_full, lower.tail = FALSE)
  }
  r1 <- sim$r_values
  list(pvalue = pv, is_rhythmic = sim$is_rhythmic, r_tilde = r1)
}

# ---- Run sims ----
if (RUN_SIMS || !file.exists(RES_PATH)) {
  cat(sprintf("Running sims: %d N values x %d sims x %d cores\n",
              length(DISPLAY_N), NSIMS, NCORES))

  t0 <- Sys.time()
  res_by_N <- vector("list", length(DISPLAY_N))
  names(res_by_N) <- as.character(DISPLAY_N)
  for (j in seq_along(DISPLAY_N)) {
    N <- DISPLAY_N[j]
    cat(sprintf("  N=%d ", N))
    sims <- mclapply(seq_len(NSIMS), function(i) {
      sim_one_rep(N, seed = SEED_BASE + (N * 1000L + i) * 7919L,
                  psi = psi, paired_sigma = FALSE)
    }, mc.cores = NCORES, mc.preschedule = FALSE)
    res_by_N[[j]] <- sims
    cat(sprintf("(%.0f s)\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
  saveRDS(res_by_N, RES_PATH)
  cat(sprintf("Saved %s\n", RES_PATH))
} else {
  res_by_N <- readRDS(RES_PATH)
}

# ---- Compute Panel A: power vs N at varying alpha (uses panel-A source) ----
res_panelA   <- if (file.exists(RES_PATH_A)) readRDS(RES_PATH_A) else res_by_N
DISPLAY_N_A  <- if (!identical(res_panelA, res_by_N)) {
                  c(20, 40, 60, 80, 100, 120, 150, 200, 250, 300)
                } else {
                  c(20, 40, 60, 80, 100, 120, 150, 200)
                }
DISPLAY_N_B  <- c(20, 40, 60, 80, 100, 120, 150, 200)
if (!identical(res_panelA, res_by_N))
  cat(sprintf("Panel A source: %s (N grid: %s)\n",
              basename(RES_PATH_A), paste(DISPLAY_N_A, collapse = ",")))

n_N_A <- length(DISPLAY_N_A)
n_N_B <- length(DISPLAY_N_B)
n_a   <- length(ALPHA_GRID)
pwr_mean <- matrix(NA_real_, nrow = n_N_A, ncol = n_a)
pwr_se   <- matrix(NA_real_, nrow = n_N_A, ncol = n_a)

for (j in seq_along(DISPLAY_N_A)) {
  sims <- res_panelA[[j]]
  # Defensive: drop any non-list entries (mclapply error captures)
  sims <- sims[vapply(sims, is.list, logical(1))]
  if (length(sims) == 0) next
  per_rep <- vapply(sims, function(s) {
    n_tgt <- sum(s$is_rhythmic)
    sapply(ALPHA_GRID, function(a) {
      q <- p.adjust(s$pvalue, "BH")
      sum(q <= a & s$is_rhythmic) / max(1, n_tgt)
    })
  }, FUN.VALUE = numeric(n_a))
  pwr_mean[j, ] <- rowMeans(per_rep, na.rm = TRUE)
  # SE of the mean across NSIMS=100 replications. Small but proper.
  pwr_se[j, ]   <- apply(per_rep, 1, sd, na.rm = TRUE) / sqrt(length(sims))
}

# ---- Compute Panel B: r1-stratified power at alpha=0.05 (primary source) ----
pwr_strat_mean <- matrix(NA_real_, nrow = n_N_B, ncol = length(R_LABELS))
pwr_strat_se   <- matrix(NA_real_, nrow = n_N_B, ncol = length(R_LABELS))

for (j in seq_along(DISPLAY_N_B)) {
  sims <- res_by_N[[j]]
  sims <- sims[vapply(sims, is.list, logical(1))]
  if (length(sims) == 0) next
  per_rep <- vapply(sims, function(s) {
    q <- p.adjust(s$pvalue, "BH")
    disc <- q <= 0.05 & s$is_rhythmic
    bins <- cut(s$r_tilde, breaks = R_BREAKS, include.lowest = TRUE, labels = FALSE)
    sapply(seq_along(R_LABELS), function(k) {
      in_k  <- s$is_rhythmic & !is.na(bins) & bins == k
      n_tgt <- sum(in_k)
      if (n_tgt == 0) NA_real_ else sum(disc & in_k) / n_tgt
    })
  }, FUN.VALUE = numeric(length(R_LABELS)))
  pwr_strat_mean[j, ] <- rowMeans(per_rep, na.rm = TRUE)
  pwr_strat_se[j, ]   <- apply(per_rep, 1, sd, na.rm = TRUE) / sqrt(length(sims))
}

# ---- N80 interpolation for Panel A vline ----
interp_n80 <- function(N_vals, p_vals, target = VLINE_POWER) {
  ok <- !is.na(p_vals) & is.finite(N_vals)
  N_vals <- N_vals[ok]; p_vals <- p_vals[ok]
  if (length(N_vals) < 2) return(NA_real_)
  if (max(p_vals) < target) return(NA_real_)
  if (min(p_vals) >= target) return(N_vals[1])
  i <- min(which(p_vals >= target))
  if (i == 1) return(N_vals[1])
  N_vals[i-1] + (target - p_vals[i-1]) * (N_vals[i] - N_vals[i-1]) /
                (p_vals[i] - p_vals[i-1])
}
idx_vfdr <- which.min(abs(ALPHA_GRID - VLINE_FDR))
vline_n  <- interp_n80(DISPLAY_N_A, pwr_mean[, idx_vfdr])

# ---- Render: 2-panel layout ----
pdf(FIG_PATH, width = 7.5, height = 4.5)
par(mfrow = c(1, 2), mai = c(1.30, 0.80, 0.50, 0.15),
    mgp = c(2.4, 0.55, 0), oma = c(0, 0, 2.4, 0),
    cex.axis = 0.85, cex.lab = 1.00, cex.main = 1.00, font.main = 2)

# ---- Panel A ----
alpha_cols <- c("darkgreen", "steelblue", "orange", "red")
matplot(DISPLAY_N_A, 100 * pwr_mean, type = "o", pch = 19, lwd = 2,
        col = alpha_cols, lty = 1, cex = 0.55,
        xlim = c(0, max(DISPLAY_N_A) * 1.05), ylim = c(0, 100),
        xlab = "Sample size (n)", ylab = "Power (%)",
        main = "")
title(main = "A   Power vs Sample Size",
      adj = 0.5, font.main = 2, cex.main = 1.00, line = 0.3)
for (a in seq_along(ALPHA_GRID)) {
  lo <- 100 * (pwr_mean[, a] - pwr_se[, a])
  hi <- 100 * (pwr_mean[, a] + pwr_se[, a])
  ok <- is.finite(lo) & is.finite(hi) & pwr_se[, a] > 0
  if (any(ok)) {
    arrows(DISPLAY_N_A[ok], lo[ok], DISPLAY_N_A[ok], hi[ok],
           code = 3, angle = 90, length = 0.04, col = alpha_cols[a], lwd = 1.2)
  }
}
abline(h = 80, lty = 2, col = "grey50", lwd = 1.2)
if (!is.na(vline_n) && is.finite(vline_n)) {
  abline(v = vline_n, lty = 2, col = adjustcolor("steelblue", 0.7), lwd = 1.4)
  text(vline_n, 35, sprintf("n=%d", round(vline_n)),
       col = "steelblue", cex = 0.70, adj = c(-0.10, 0.5), font = 2)
}
grid()
legend("bottomright", paste0("FDR ", round(100*ALPHA_GRID), "%"),
       col = alpha_cols, lty = 1, pch = 19, lwd = 1.5,
       cex = 0.55, bty = "o", box.col = "grey70", box.lwd = 0.5,
       inset = 0.01, y.intersp = 0.85, bg = "white")

# ---- Panel B ----
par(mgp = c(3.6, 0.55, 0))
size_colors <- rainbow(length(DISPLAY_N_B), s = 0.6, v = 0.8)
matplot(seq_along(R_LABELS), 100 * t(pwr_strat_mean),
        type = "l", lwd = 2, col = size_colors, lty = 1,
        xlim = c(0.5, length(R_LABELS) + 0.5), ylim = c(0, 100), bty = "l",
        xlab = expression(tilde(r) == A/sigma),
        ylab = "Power (%)",
        main = "", xaxt = "n")
title(main = bquote(bold("B   ") * bold("Stratified Power by") ~
                     bold(tilde(r)) ~ bold("(FDR 5%)")),
      adj = 0.5, font.main = 2, cex.main = 1.00, line = 0.3)
axis(1, at = seq_along(R_LABELS), labels = R_LABELS, las = 2, cex.axis = 0.65)
for (j in seq_len(n_N_B)) {
  points(seq_along(R_LABELS), 100 * pwr_strat_mean[j, ],
         pch = 19, col = size_colors[j], cex = 0.65)
  arrows(seq_along(R_LABELS),
         100*(pwr_strat_mean[j, ] - pwr_strat_se[j, ]),
         seq_along(R_LABELS),
         100*(pwr_strat_mean[j, ] + pwr_strat_se[j, ]),
         code = 3, angle = 90, length = 0.04, col = size_colors[j], lwd = 1.2)
}
grid()
legend("bottomright", paste0("n=", DISPLAY_N_B),
       col = size_colors, lty = 1, lwd = 1.2,
       cex = 0.45, bty = "o", box.col = "grey70", box.lwd = 0.5,
       inset = 0.01, y.intersp = 0.78, bg = "white", seg.len = 1.2)
par(mgp = c(3.0, 0.6, 0))

mtext("Two-harmonic Single-Cohort Power Analysis (GTEx Liver)",
      outer = TRUE, side = 3, line = 0.2, font = 2, cex = 1.15)
dev.off()
cat(sprintf("Wrote %s (%.1f KB)\n", FIG_PATH, file.size(FIG_PATH)/1024))
cat(sprintf("N80 (FDR %g%%): %s\n", VLINE_FDR * 100,
            ifelse(is.na(vline_n), "out of range",
                   sprintf("%d", round(vline_n)))))
