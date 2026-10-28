#!/usr/bin/env Rscript
# =============================================================================
# Fig 6 A/B - Balanced ACTIVE design, GTEx Liver, B-vs-m trade-off
#
#   Panel A: single-harmonic (K=1) biomarker detection power vs total N,
#            swept over B (distinct collection time points).
#   Panel B: two-harmonic  (K=2) biomarker detection power vs total N,
#            swept over B.
#
# The two panels use TWO DIFFERENT Liver pilots / TWO DIFFERENT true-positive
# denominators (this is the whole point of the rebuild):
#
#   Panel A true set = the K=1 pilot's paired top-K single-harmonic rhythmic set
#                      (data/gtex_Liver_single_pilot.rds; amplitude_spec paired
#                       with sigma_rhythmic, phase_spec; median r1 ~ 0.48, tight).
#                      Detector: K=1 cosinor F(2, N-3), BH-FDR 0.05.
#
#   Panel B true set = the K=2 pilot's paired top-K two-harmonic rhythmic set
#                      (output/two_harmonic/results/pilot_2h_GTExLiver_topK300.rds;
#                       includes weak-first-harmonic / genuine-12h genes, so the
#                       first-harmonic effect size is weaker, median r1 ~ 0.35).
#                      Detector: K=2 two-harmonic F(4, N-5), BH-FDR 0.05.
#
# Because the K=2 true set has weaker first-harmonic effect sizes AND the K=2
# test spends two extra degrees of freedom, Panel B reaches 80% power at a
# LARGER N than Panel A (the "two-harmonic sample-size penalty").
#
# Balanced active design: B distinct equally spaced TODs t_b = (b-1)*24/B, each
# with m = N/B replicates. B is swept; N runs over multiples of 24 so every B
# in each sweep divides N exactly.
#
# Run:  MC_CORES=48 Rscript examples/publication/two_harmonic/fig6AB_liver_active.R
# Saves: output/two_harmonic/results/fig6AB_liver_active.rds
#        + a diagnostic PNG to scratchpad/ for eyeballing.
# =============================================================================

setwd(Sys.getenv("SCP_ROOT", unset =
  "/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim"))
suppressPackageStartupMessages({ library(SCP); library(parallel) })

# ---- Config ----
PERIOD  <- 24
OMEGA0  <- 2 * pi / PERIOD
NSIMS   <- as.integer(Sys.getenv("NSIMS",   unset = "100"))
NGENES  <- as.integer(Sys.getenv("NGENES",  unset = "2500"))
FDR     <- 0.05
MC_CORES <- as.integer(Sys.getenv("MC_CORES", unset = "48"))
MC_CORES <- max(1L, min(MC_CORES, parallel::detectCores()))
SEED_BASE <- 20260724L
# Integer-safe seed (avoid 32-bit overflow at large N: (N*1000)*7919 can exceed
# .Machine$integer.max, which would silently produce NA seeds).
mk_seed <- function(N, B, i)
  as.integer((SEED_BASE + (N * 1000 + B * 31 + i) * 7919) %% 2147483647)

# N grid: multiples of 24 so every swept B divides N exactly (balanced m = N/B).
N_GRID  <- as.integer(seq(24, 288, by = 24))
B_K1    <- c(3L, 4L, 6L, 8L, 12L, 24L)   # K=1 identifiability: B >= 3
B_K2    <- c(6L, 8L, 12L, 24L)           # K=2 identifiability: B >= 5

OUT_RDS <- "output/two_harmonic/results/fig6AB_liver_active.rds"
OUT_PNG <- file.path(
  "/tmp/claude-547180/-home-qtp1-Projects-Circadian-Kyle-Circadian-analysis-main-R-v1-PowerSim",
  "491a8454-ce4d-44cf-b7af-479d8b17b691/scratchpad",
  "fig6AB_liver_active_diag.png")

# ---- Pilots ----
p1_raw <- readRDS("data/gtex_Liver_single_pilot.rds")     # K=1 single-harmonic
p2     <- readRDS("output/two_harmonic/results/pilot_2h_GTExLiver_topK300.rds")  # K=2

# Panel A pilot: paired top-K single-harmonic set. The shipped K=1 pilot stores
# the full rhythmic amplitude/phase (length 839) but only a top-K sigma vector
# (length 300); pairing amplitude_spec/phase_spec (both length 300) with
# sigma_rhythmic makes length(amplitude) == length(sigma_rhythmic), which trips
# simCircadianSingleCohort's `has_joint` branch => genuine paired (A, sigma, phi)
# draws (median r1 ~ 0.48, tight), the narrow single-harmonic true set.
p1 <- p1_raw
p1$amplitude <- p1_raw$amplitude_spec
p1$phase     <- p1_raw$phase_spec
stopifnot(length(p1$amplitude) == length(p1$sigma_rhythmic),
          length(p1$phase)     == length(p1$sigma_rhythmic))

cat(sprintf("Panel A (K=1) pilot: ngenes=%d prop_rhy=%.4f  median r1=%.3f (paired top-K=%d)\n",
            p1$ngenes, p1$prop_rhythmic,
            median(p1$amplitude / p1$sigma_rhythmic), length(p1$amplitude)))
cat(sprintf("Panel B (K=2) pilot: ngenes=%d prop_rhy=%.4f  median r1=%.3f  A2/A1 med=%.2f (top-K=%d)\n",
            p2$ngenes, p2$prop_rhythmic,
            median(p2$amplitude / p2$sigma_rhythmic),
            median(p2$amplitude2 / p2$amplitude), length(p2$amplitude)))

# ---- Balanced active collection grid ----
active_times <- function(N, B) {
  m  <- N %/% B
  ts <- rep(seq(0, PERIOD - PERIOD / B, length.out = B), each = m)
  if (length(ts) < N)
    ts <- c(ts, sample(seq(0, PERIOD - PERIOD / B, length.out = B),
                       N - length(ts), replace = TRUE))
  ts[seq_len(min(N, length(ts)))]
}

# ---- Vectorized K-harmonic rhythmicity F-test -> genome-wide power ----
ftest_power <- function(expr, ts, is_rhy, K) {
  N <- length(ts)
  if (K == 1L) { X <- cbind(1, cos(OMEGA0 * ts), sin(OMEGA0 * ts)); dfn <- 2L }
  else         { X <- cbind(1, cos(OMEGA0 * ts), sin(OMEGA0 * ts),
                            cos(2 * OMEGA0 * ts), sin(2 * OMEGA0 * ts)); dfn <- 4L }
  dfd <- N - ncol(X)
  if (dfd <= 0L) return(NA_real_)
  gmean <- rowMeans(expr)
  R0 <- rowSums((expr - gmean)^2)                 # per-gene null RSS
  fit <- lm.fit(X, t(expr))                       # fit all genes at once
  R  <- colSums(fit$residuals^2)                  # per-gene full-model RSS
  Fstat <- ((R0 - R) / dfn) / (R / dfd)
  pv <- pf(Fstat, dfn, dfd, lower.tail = FALSE)
  pv[!is.finite(pv) | R0 <= 0 | R <= 0] <- 1
  q <- p.adjust(pv, "BH")
  sum(q <= FDR & is_rhy) / max(1L, sum(is_rhy))
}

# ---- One replicate: Panel A (K=1) ----
sim_K1 <- function(N, B, seed) {
  ts <- active_times(N, B)
  if (length(unique(ts)) < 3L) return(NA_real_)
  set.seed(seed)
  pl <- p1; pl$ngenes <- NGENES
  s  <- SCP:::simCircadianSingleCohort(pl, ts, seed = seed)   # K=1 paired sim
  ftest_power(s$expr, ts, s$is_rhythmic, K = 1L)
}

# ---- One replicate: Panel B (K=2) ----
sim_K2 <- function(N, B, seed) {
  ts <- active_times(N, B)
  if (length(unique(ts)) < 5L) return(NA_real_)
  set.seed(seed)
  pl <- p2; pl$ngenes <- NGENES
  s  <- simCircadianSingleCohort2H(pl, ts, seed = seed)       # K=2 paired sim
  ftest_power(s$expr, ts, s$is_rhythmic, K = 2L)
}

# ---- Sweep N x B ----
run_grid <- function(sim_fn, B_grid, label) {
  cat(sprintf("\n== %s ==\n", label))
  M  <- matrix(NA_real_, length(N_GRID), length(B_grid),
               dimnames = list(paste0("N=", N_GRID), paste0("B=", B_grid)))
  SE <- M
  t0 <- Sys.time()
  for (k in seq_along(B_grid)) {
    for (j in seq_along(N_GRID)) {
      N <- N_GRID[j]; B <- B_grid[k]
      v <- mclapply(seq_len(NSIMS), function(i)
             sim_fn(N, B, mk_seed(N, B, i)),
             mc.cores = MC_CORES, mc.preschedule = FALSE)
      v <- unlist(v); v <- v[is.finite(v)]
      if (length(v) > 0L) { M[j, k] <- mean(v); SE[j, k] <- sd(v) / sqrt(length(v)) }
    }
    cat(sprintf("  B=%2d: %s  (%.0fs)\n", B_grid[k],
                paste(sprintf("%4.0f", 100 * M[, k]), collapse = " "),
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
  list(M = M, SE = SE)
}

t_start <- Sys.time()
res_K1 <- run_grid(sim_K1, B_K1, "Panel A: K=1 cosinor, paired narrow top-K")
res_K2 <- run_grid(sim_K2, B_K2, "Panel B: K=2 two-harmonic, paired top-K")
runtime_min <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))

# ---- N80 interpolation helper ----
interp_n80 <- function(N_vals, p_vals, target = 0.80) {
  ok <- is.finite(p_vals) & is.finite(N_vals)
  N_vals <- N_vals[ok]; p_vals <- p_vals[ok]
  if (length(N_vals) < 2 || max(p_vals) < target) return(NA_real_)
  if (min(p_vals) >= target) return(N_vals[1])
  i <- min(which(p_vals >= target))
  if (i == 1) return(N_vals[1])
  N_vals[i-1] + (target - p_vals[i-1]) * (N_vals[i] - N_vals[i-1]) /
                (p_vals[i] - p_vals[i-1])
}
n80_K1 <- vapply(seq_along(B_K1), function(k) interp_n80(N_GRID, res_K1$M[, k]), numeric(1))
n80_K2 <- vapply(seq_along(B_K2), function(k) interp_n80(N_GRID, res_K2$M[, k]), numeric(1))
names(n80_K1) <- paste0("B=", B_K1); names(n80_K2) <- paste0("B=", B_K2)

# ---- Save cache ----
dir.create(dirname(OUT_RDS), recursive = TRUE, showWarnings = FALSE)
saveRDS(list(
  N       = N_GRID,
  B_K1    = B_K1, B_K2 = B_K2,
  pwr_K1  = res_K1$M, se_K1 = res_K1$SE,
  pwr_K2  = res_K2$M, se_K2 = res_K2$SE,
  n80_K1  = n80_K1,  n80_K2 = n80_K2,
  NSIMS   = NSIMS, NGENES = NGENES, FDR = FDR, PERIOD = PERIOD,
  pilot_A = "data/gtex_Liver_single_pilot.rds (paired top-K, K=1)",
  pilot_B = "output/two_harmonic/results/pilot_2h_GTExLiver_topK300.rds (K=2)",
  meta    = list(design = "balanced active t_b=(b-1)*24/B, m=N/B",
                 detector_A = "K=1 cosinor F(2,N-3), BH 0.05",
                 detector_B = "K=2 two-harmonic F(4,N-5), BH 0.05",
                 runtime_min = runtime_min)
), OUT_RDS)
cat(sprintf("\nSaved cache: %s\n", OUT_RDS))

# ---- Diagnostic PNG ----
dir.create(dirname(OUT_PNG), recursive = TRUE, showWarnings = FALSE)
tryCatch(png(OUT_PNG, width = 1200, height = 560, res = 110, type = "cairo"),
         error = function(e) png(OUT_PNG, width = 1200, height = 560, res = 110))
par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.6, 0.8))
cols1 <- c("#440154","#3b528b","#21918c","#5ec962","#addc30","#fde725")
cols2 <- c("#3b528b","#21918c","#5ec962","#fde725")
matplot(N_GRID, 100 * res_K1$M, type = "o", pch = 19, lwd = 2, lty = 1,
        col = cols1[seq_along(B_K1)], ylim = c(0, 100),
        xlab = "Total N", ylab = "Power (%)",
        main = "A  K=1 cosinor (Liver, active)")
abline(h = 80, lty = 2, col = "grey50"); grid()
legend("bottomright", paste0("B=", B_K1), col = cols1[seq_along(B_K1)],
       lty = 1, pch = 19, lwd = 2, cex = 0.8, bg = "white")
matplot(N_GRID, 100 * res_K2$M, type = "o", pch = 19, lwd = 2, lty = 1,
        col = cols2[seq_along(B_K2)], ylim = c(0, 100),
        xlab = "Total N", ylab = "Power (%)",
        main = "B  K=2 two-harmonic (Liver, active)")
abline(h = 80, lty = 2, col = "grey50"); grid()
legend("bottomright", paste0("B=", B_K2), col = cols2[seq_along(B_K2)],
       lty = 1, pch = 19, lwd = 2, cex = 0.8, bg = "white")
dev.off()
cat(sprintf("Saved diagnostic PNG: %s\n", OUT_PNG))

# ---- Console report ----
cat(sprintf("\n---- N80 (power first >= 80%%), by B ----\n"))
cat("Panel A (K=1):\n"); for (k in seq_along(B_K1))
  cat(sprintf("  B=%2d -> N80 = %s\n", B_K1[k],
              ifelse(is.na(n80_K1[k]), ">max", sprintf("%.0f", n80_K1[k]))))
cat("Panel B (K=2):\n"); for (k in seq_along(B_K2))
  cat(sprintf("  B=%2d -> N80 = %s\n", B_K2[k],
              ifelse(is.na(n80_K2[k]), ">max", sprintf("%.0f", n80_K2[k]))))
cat(sprintf("\nMedian N80  K=1 = %.0f   K=2 = %.0f   (penalty = K2 - K1 = %.0f)\n",
            median(n80_K1, na.rm = TRUE), median(n80_K2, na.rm = TRUE),
            median(n80_K2, na.rm = TRUE) - median(n80_K1, na.rm = TRUE)))
cat(sprintf("Runtime: %.1f min | NSIMS=%d NGENES=%d cores=%d\n",
            runtime_min, NSIMS, NGENES, MC_CORES))
