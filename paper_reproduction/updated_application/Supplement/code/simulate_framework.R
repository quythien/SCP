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

