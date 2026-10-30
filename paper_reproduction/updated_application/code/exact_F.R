# Exact unedited function blocks from the validated engine.
# Source SHA256: 0a0e204d2ea0027a4d4afebf6de1d9374a1c1e312273e3762110b43e78a7e0c6

scp_time_sampler <- function(times) {
  circular_times <- times %% 24
  density <- density(c(circular_times - 24, circular_times,
                       circular_times + 24), from = 0, to = 24)
  list(x = density$x, probability = pmax(density$y, 1e-10))
}

scp_sample_cached_times <- function(N, sampler) {
  (sample(sampler$x, size = N, replace = TRUE,
           prob = sampler$probability) + rnorm(N, 0, 0.1)) %% 24
}

scp_f_noncentrality <- function(amplitude, phase, sigma, times, period = 24) {
  stopifnot(length(times) > 3L, all(is.finite(times)),
            length(amplitude) == length(phase), length(phase) == length(sigma),
            all(is.finite(amplitude)), all(is.finite(phase)),
            all(is.finite(sigma)), all(sigma > 0))
  omega <- 2 * pi / period
  cosine <- cos(omega * times)
  sine <- sin(omega * times)
  cosine <- cosine - mean(cosine)
  sine <- sine - mean(sine)
  g11 <- sum(cosine^2)
  g12 <- sum(cosine * sine)
  g22 <- sum(sine^2)
  stopifnot(g11 * g22 - g12^2 > 1e-12 * max(g11, g22)^2)
  b1 <- amplitude * cos(omega * phase)
  b2 <- amplitude * sin(omega * phase)
  value <- (g11 * b1^2 + 2 * g12 * b1 * b2 + g22 * b2^2) / sigma^2
  stopifnot(all(is.finite(value)), all(value >= -1e-10))
  pmax(value, 0)
}

simulate_scp_f_replicate <- function(bio, times, N, seed,
                                     retain_parameters = FALSE,
                                     time_sampler = NULL) {
  stopifnot(!isTRUE(bio$paired_2h), is.null(bio$alpha2_dist),
            is.null(bio$alpha3_dist), is.null(bio$omega_dist),
            is.null(bio$alpha_dist), is.null(bio$omega_rhythmic),
            is.null(bio$alpha_rhythmic))
  set.seed(seed, kind = "L'Ecuyer-CMRG")
  time_stream <- parallel::nextRNGStream(.Random.seed)
  statistic_stream <- parallel::nextRNGStream(time_stream)
  G <- bio$ngenes
  n_targets <- round(G * bio$prop_rhythmic)
  target_id <- sample(G, n_targets)
  truth <- logical(G)
  truth[target_id] <- TRUE
  stopifnot(length(bio$amplitude) == length(bio$sigma_rhythmic))
  index <- sample(length(bio$amplitude), n_targets, replace = TRUE)
  amplitude <- pmax(bio$amplitude[index], 0.05)
  sigma <- pmax(bio$sigma_rhythmic[index], 1e-6)
  phase <- sample(bio$phase, n_targets, replace = TRUE)

  assign(".Random.seed", time_stream, envir = .GlobalEnv)
  sampled_times <- if (is.null(time_sampler)) sample_scp_times(N, times) else
    scp_sample_cached_times(N, time_sampler)
  lambda <- scp_f_noncentrality(amplitude, phase, sigma, sampled_times, bio$period)
  assign(".Random.seed", statistic_stream, envir = .GlobalEnv)
  pvalue <- numeric(G)
  pvalue[!truth] <- runif(G - n_targets)
  f_statistic <- rf(n_targets, df1 = 2, df2 = N - 3, ncp = lambda)
  pvalue[truth] <- pf(f_statistic, df1 = 2, df2 = N - 3, lower.tail = FALSE)
  stopifnot(all(is.finite(pvalue)), all(pvalue >= 0 & pvalue <= 1))
  discovered <- p.adjust(pvalue, "BH") <= 0.05
  TD <- sum(discovered & truth)
  FD <- sum(discovered & !truth)
  result <- list(power = if (n_targets > 0L) TD / n_targets else NA_real_,
                 TD = TD, FD = FD, n_targets = n_targets,
                 FDP = FD / max(TD + FD, 1L))
  if (retain_parameters) {
    result$amplitude <- amplitude
    result$phase <- phase
    result$sigma <- sigma
    result$times <- sampled_times
    result$truth <- truth
    result$lambda <- lambda
    result$pvalue <- pvalue
  }
  result
}

