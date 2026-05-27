#!/usr/bin/env Rscript
# Fig 5 Adrenal sims (K=2 GTEx Adrenal Gland) - both paired and unpaired
setwd("/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")
sink(stderr())
source("code/utils.R"); source("code/options.R"); source("code/simulation.R")
source("code/estimation.R"); source("code/bootstrap_sim.R"); source("code/detection.R")
sink()
suppressPackageStartupMessages(library(parallel))

PILOT_PATH       <- "output/two_harmonic/results/pilot_2h_GTExAdrenal_topK300.rds"
RES_PATH_PAIRED  <- "output/two_harmonic/results/fig5_panelA_GTExAdrenal_paired.rds"
RES_PATH_UNPAIR  <- "output/two_harmonic/results/fig5_panelAB_GTExAdrenal_unpaired.rds"

PERIOD     <- 24
DISPLAY_N  <- c(20, 40, 60, 80, 100, 120, 150, 200, 250, 300)
NSIMS      <- 100
NGENES_SIM <- 2000
NCORES     <- min(12, parallel::detectCores() - 4)
SEED_BASE  <- 20260526

psi <- readRDS(PILOT_PATH)
cat(sprintf("Pilot: n=%d   top_k=%d   A2/A1 med=%.2f   prop_rhythmic=%.3f\n",
            attr(psi, "n_pilot"),
            psi$diagnostics$top_k_used,
            psi$diagnostics$A2_over_A1_med,
            psi$prop_rhythmic))

sim_one_rep <- function(N, seed, psi, ngenes = NGENES_SIM,
                         paired_sigma = TRUE) {
  set.seed(seed)
  B <- min(24L, N)
  m <- max(1L, N %/% B)
  ts <- rep(seq(0, PERIOD - PERIOD/B, length.out = B), each = m)
  if (length(ts) > N) ts <- ts[seq_len(N)]
  if (length(ts) < N) ts <- c(ts, sample(ts, N - length(ts), replace = TRUE))
  psi_local <- psi; psi_local$ngenes <- ngenes
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
  list(pvalue = pv, is_rhythmic = sim$is_rhythmic, r_tilde = sim$r_values)
}

run_grid <- function(paired_sigma, out_path) {
  cat(sprintf("\n== Adrenal K=2 sims (paired_sigma=%s) ==\n", paired_sigma))
  t0 <- Sys.time()
  res_by_N <- vector("list", length(DISPLAY_N))
  names(res_by_N) <- as.character(DISPLAY_N)
  for (j in seq_along(DISPLAY_N)) {
    N <- DISPLAY_N[j]
    cat(sprintf("  N=%d ", N))
    sims <- mclapply(seq_len(NSIMS), function(i) {
      sim_one_rep(N, seed = SEED_BASE + (N * 1000L + i) * 7919L,
                  psi = psi, paired_sigma = paired_sigma)
    }, mc.cores = NCORES, mc.preschedule = FALSE)
    res_by_N[[j]] <- sims
    cat(sprintf("(%.0f s)\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
  saveRDS(res_by_N, out_path)
  cat(sprintf("Saved: %s\n", out_path))
}
run_grid(TRUE, RES_PATH_PAIRED)
run_grid(FALSE, RES_PATH_UNPAIR)
cat("\n== Done ==\n")
