# Compare sample-size requirements under the supplied pilot calibrations.
# Each point is an exact Gaussian F-statistic Monte Carlo experiment with BH
# at 5 percent. An expression-matrix engine is retained for validation.
# Source with SCP_SIM_DEFINE_ONLY=true to load the helpers without running it.

# 1. Load the local package and choose reproducible settings -------------------

stage_root <- Sys.getenv("SCP_GTEX_V10_STAGE", "output/gtex_v10_chiral_protocol")
simulation_root <- file.path(stage_root, "simulation")
simulation_library <- Sys.getenv("SCP_SIM_LIBRARY")
if (nzchar(simulation_library)) .libPaths(c(simulation_library, .libPaths()))
suppressPackageStartupMessages(library(SCP))
stopifnot(utils::packageVersion("SCP") >= "0.5.0")
.CPP_LOADED <- TRUE
stopifnot(exists("sim_cosinor_expr_cpp", asNamespace("SCP"), inherits = FALSE),
          exists("cosinor_pvals_cpp", asNamespace("SCP"), inherits = FALSE))
simulation_cpp <- getFromNamespace("sim_cosinor_expr_cpp", "SCP")
test_cpp <- getFromNamespace("cosinor_pvals_cpp", "SCP")
sample_scp_times <- getFromNamespace("sampleTimesFromDist", "SCP")


# 2. Simulate one complete expression experiment ------------------------------
# Draw parameters in the same order as runSimsSingleCohort(). The phase draw
# is independent of the paired amplitude/noise draw, as in the package runner.

simulate_scp_replicate <- function(bio, times, N, seed,
                                   retain_pvalues = FALSE,
                                   retain_expression = FALSE) {
  set.seed(seed)
  G <- bio$ngenes
  n_targets <- round(G * bio$prop_rhythmic)
  target_id <- sample(G, n_targets)
  truth <- logical(G)
  truth[target_id] <- TRUE

  mesor <- bio$lBaselineExpr
  sigma <- exp(bio$lOD)
  stopifnot(length(mesor) == length(sigma))
  if (length(mesor) != G) {
    index <- sample(length(mesor), G, replace = TRUE)
    mesor <- mesor[index]
    sigma <- sigma[index]
  }
  amplitude <- numeric(G)
  phase <- numeric(G)
  if (n_targets > 0L) {
    stopifnot(length(bio$amplitude) == length(bio$sigma_rhythmic))
    index <- sample(length(bio$amplitude), n_targets, replace = TRUE)
    amplitude[target_id] <- pmax(bio$amplitude[index], 0.05)
    sigma[target_id] <- pmax(bio$sigma_rhythmic[index], 1e-6)
    phase[target_id] <- sample(bio$phase, n_targets, replace = TRUE)
  }

  sampled_times <- sample_scp_times(N, times)
  expression <- simulation_cpp(mesor, amplitude, phase, sigma, sampled_times,
                                period = bio$period, alpha2 = 0, alpha3 = 0)
  pvalue <- test_cpp(expression, sampled_times, period = bio$period)
  stopifnot(length(pvalue) == G, all(is.finite(pvalue)),
            all(pvalue >= 0 & pvalue <= 1))
  discovered <- p.adjust(pvalue, "BH") <= 0.05
  TD <- sum(discovered & truth)
  FD <- sum(discovered & !truth)
  result <- list(power = if (n_targets > 0L) TD / n_targets else NA_real_,
                 TD = TD, FD = FD, n_targets = n_targets,
                 FDP = FD / max(TD + FD, 1L))
  if (retain_pvalues) result$pvalue <- pvalue
  if (retain_expression) {
    result$expression <- expression
    result$times <- sampled_times
    result$truth <- truth
  }
  result
}


# 3. Simulate the exact conditional Gaussian F-statistic distribution ----------
# For a pure Gaussian cosinor, the fitted F statistic has the noncentral F law
# below. Sampling this law preserves the complete joint BH experiment without
# allocating a genes-by-subjects expression matrix for every replicate.

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


# 4. Keep seeds stable when resuming points or changing worker counts ----------
# The same N and replicate use the same starting seed across calibrations.
# Within each case every replicate has a unique seed. Shared seeds across
# cases are intentional and do not make the different pilot laws identical.

scp_replicate_seed <- function(N, replicate_id, seed_base = 20261027L) {
  stopifnot(all(replicate_id > 0L & replicate_id < 100000L))
  seed <- as.double(seed_base) + as.double(N) * 100000 + replicate_id
  stopifnot(all(seed < .Machine$integer.max))
  as.integer(seed)
}

scp_case_id <- function(entry) {
  gsub("[^A-Za-z0-9]+", "_", paste(entry$tissue, entry$time_basis, sep = "__"))
}


# 5. Summarize variation between simulated experiments -----------------------
# MC intervals describe the estimated mean under a fixed fitted pilot.
# They do not include pilot-estimation uncertainty or interpolation-grid error.

summarize_scp_point <- function(replicates, entry, N, G) {
  power <- replicates$power
  fdp <- replicates$FDP
  stopifnot(all(is.finite(power)), all(is.finite(fdp)))
  se_power <- sd(power) / sqrt(length(power))
  data.frame(
    case_id = scp_case_id(entry), tissue = entry$tissue,
    time_basis = entry$time_basis, n_pilot = entry$n_pilot,
    genes_fitted = entry$genes_fitted, G_simulated = G,
    N = N, nsims = nrow(replicates), power = mean(power),
    power_mcse = se_power,
    power_mc95_low = max(0, mean(power) - 1.96 * se_power),
    power_mc95_high = min(1, mean(power) + 1.96 * se_power),
    power_q025 = unname(quantile(power, 0.025)),
    power_q975 = unname(quantile(power, 0.975)),
    unconditional_FDR = mean(fdp), FDR_mcse = sd(fdp) / sqrt(length(fdp)),
    mean_TD = mean(replicates$TD), mean_FD = mean(replicates$FD),
    n_targets = unique(replicates$n_targets),
    zero_discovery_fraction = mean(replicates$TD + replicates$FD == 0),
    stringsAsFactors = FALSE
  )
}


# 6. Run or extend one sample size and retain its replicate summaries ---------

run_scp_point <- function(entry, N, nsims, cores, G, output_root,
                           engine = "exact_gaussian_F") {
  case_id <- scp_case_id(entry)
  checkpoint <- file.path(output_root, "replicates", case_id,
                          sprintf("N_%05d.rds", N))
  bio <- entry$bio
  bio$ngenes <- G
  signature <- digest::digest(list(bio = bio, times = entry$times,
                                   engine = engine, recipe_version = 2L,
                                   seed_stride = 100000L), algo = "sha256")
  existing <- NULL
  if (file.exists(checkpoint)) {
    saved <- readRDS(checkpoint)
    stopifnot(identical(saved$signature, signature), saved$N == N,
              saved$G == G, identical(saved$engine, engine),
              identical(saved$replicates$seed,
                        scp_replicate_seed(N, saved$replicates$replicate_id)))
    existing <- saved$replicates
    if (nrow(existing) >= nsims) return(saved$summary)
  }
  previous_n <- if (is.null(existing)) 0L else nrow(existing)
  ids <- seq.int(previous_n + 1L, nsims)
  seeds <- scp_replicate_seed(N, ids)
  sampler <- scp_time_sampler(entry$times)

  started <- proc.time()[["elapsed"]]
  values <- parallel::mclapply(seq_along(ids), function(j) {
    value <- if (engine == "exact_gaussian_F") {
      simulate_scp_f_replicate(bio, entry$times, N, seeds[j], time_sampler = sampler)
    } else if (engine == "actual_expression_cpp") {
      simulate_scp_replicate(bio, entry$times, N, seeds[j])
    } else stop("Unknown simulation engine.")
    data.frame(replicate_id = ids[j], seed = seeds[j], power = value$power,
                TD = value$TD, FD = value$FD, n_targets = value$n_targets,
                FDP = value$FDP)
  }, mc.cores = cores, mc.set.seed = FALSE, mc.preschedule = TRUE)
  stopifnot(all(vapply(values, is.data.frame, logical(1))))
  replicates <- rbind(existing, do.call(rbind, values))
  stopifnot(!anyDuplicated(replicates$replicate_id),
            !anyDuplicated(replicates$seed))
  summary <- summarize_scp_point(replicates, entry, N, G)
  summary$elapsed_seconds <- proc.time()[["elapsed"]] - started +
    if (is.null(existing)) 0 else saved$summary$elapsed_seconds
  dir.create(dirname(checkpoint), recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(signature = signature, N = N, G = G, nsims = nrow(replicates),
               engine = engine, shared_seeds_across_cases = TRUE,
               replicates = replicates, summary = summary), checkpoint)
  cat(case_id, "N", N, "power", round(summary$power, 4),
      "FDR", round(summary$unconditional_FDR, 4),
      "seconds", round(summary$elapsed_seconds, 1), "\n")
  flush.console()
  summary
}


# 7. Find a tested bracket around 80 percent power ----------------------------

scp_n80_bracket <- function(points) {
  points <- points[order(points$N), , drop = FALSE]
  high <- which(points$power >= 0.80)
  if (!length(high)) return(list(low = tail(points, 1), high = NULL))
  first_high <- high[1]
  low <- if (first_high > 1L) points[first_high - 1L, , drop = FALSE] else NULL
  list(low = low, high = points[first_high, , drop = FALSE])
}

scp_n80_summary <- function(points, entry) {
  bracket <- scp_n80_bracket(points)
  lower <- if (is.null(bracket$low)) NA_integer_ else bracket$low$N
  upper <- if (is.null(bracket$high)) NA_integer_ else bracket$high$N
  estimate <- NA_real_
  status <- if (is.null(bracket$high)) "not_reached" else "below_minimum_tested"
  if (!is.null(bracket$low) && !is.null(bracket$high)) {
    low <- bracket$low
    high <- bracket$high
    estimate <- low$N + (0.80 - low$power) * (high$N - low$N) /
      (high$power - low$power)
    status <- "bracketed"
  }
  data.frame(case_id = scp_case_id(entry), tissue = entry$tissue,
             time_basis = entry$time_basis, n_pilot = entry$n_pilot,
             N80_low = lower, N80_high = upper, N80_interpolated = estimate,
             N80_bracket_width = upper - lower,
             N80_bracket_relative_width = (upper - lower) / estimate,
             status = status,
             stringsAsFactors = FALSE)
}


# 8. Estimate Monte Carlo uncertainty in the interpolated crossing ------------

scp_n80_bootstrap <- function(points, entry, output_root, B = 1000L) {
  power_vectors <- lapply(points$N, function(N) {
    file <- file.path(output_root, "replicates", scp_case_id(entry),
                       sprintf("N_%05d.rds", N))
    readRDS(file)$replicates$power
  })
  set.seed(20261027L)
  crossings <- replicate(B, {
    bootstrap <- points
    bootstrap$power <- vapply(power_vectors, function(x) {
      mean(x[sample.int(length(x), length(x), replace = TRUE)])
    }, numeric(1))
    scp_n80_summary(bootstrap, entry)$N80_interpolated
  })
  finite <- crossings[is.finite(crossings)]
  if (!length(finite)) return(c(low = NA_real_, high = NA_real_, valid_fraction = 0))
  c(low = unname(quantile(finite, 0.025)), high = unname(quantile(finite, 0.975)),
    valid_fraction = length(finite) / B)
}


# 9. Run the coarse grid and refine only the crossing -------------------------

if (!identical(Sys.getenv("SCP_SIM_DEFINE_ONLY"), "true")) {
  nsims <- as.integer(Sys.getenv("SCP_NSIM", "200"))
  nsims_final <- max(nsims, as.integer(Sys.getenv("SCP_NSIM_FINAL", "1000")))
  cores <- min(6L, as.integer(Sys.getenv("SCP_CORES", "4")))
  G <- as.integer(Sys.getenv("SCP_NGENES", "5000"))
  max_N <- as.integer(Sys.getenv("SCP_MAX_N", "10000"))
  coarse_N <- as.integer(strsplit(Sys.getenv("SCP_COARSE_N",
    "12,24,48,96,192,384,768,1536,3072,6144,10000"), ",", fixed = TRUE)[[1]])
  coarse_N <- sort(unique(coarse_N[coarse_N >= 6L & coarse_N <= max_N]))
  stopifnot(nsims >= 2L, nsims_final < 100000L, cores >= 1L, G >= 100L,
            length(coarse_N) > 0L, max_N <= 20000L)
  pilot_file <- Sys.getenv("SCP_SIM_PILOTS", file.path(stage_root, "private", "report_pilots.rds"))
  output_root <- Sys.getenv("SCP_SIM_OUTPUT", file.path(simulation_root, "results_exactF"))
  entries <- readRDS(pilot_file)
  if (!is.null(entries$pilots)) entries <- entries$pilots
  if (all(vapply(entries, function(x) !is.null(x$cohort), logical(1))))
    entries <- Filter(function(x) identical(x$cohort, "matched"), entries)
  case_ids <- vapply(entries, scp_case_id, character(1))
  stopifnot(!anyDuplicated(case_ids))
  case_filter <- Sys.getenv("SCP_CASE_FILTER")
  if (nzchar(case_filter)) entries <- entries[grepl(case_filter, case_ids)]
  stopifnot(length(entries) > 0L)
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  run_settings <- list(
    engine = "Exact Gaussian F-statistic Monte Carlo",
    expression_validation_engine = "actual_expression_cpp",
    G = G, nsims_initial = nsims, nsims_final = nsims_final, cores = cores,
    coarse_N = coarse_N, maximum_N = max_N, BH_threshold = 0.05,
    matched_cases = length(entries), shared_seeds_across_cases = TRUE,
    seed_rule = "20261027 + N * 100000 + replicate_id",
    RNG = "L'Ecuyer-CMRG with parameter, time, and statistic substreams",
    bootstrap_replicates = 1000L,
    pilot_file_sha256 = digest::digest(file = pilot_file, algo = "sha256"),
    engine_sha256 = digest::digest(file = "examples/pilots/gtex_v10/05_compare_sample_sizes.R", algo = "sha256"),
    MC_interval_scope = "Fixed pilots. Excludes pilot uncertainty and interpolation-grid error."
  )
  saveRDS(run_settings, file.path(output_root, "run_settings.rds"))
  jsonlite::write_json(run_settings, file.path(output_root, "run_settings.json"),
                        pretty = TRUE, auto_unbox = TRUE)
  writeLines(capture.output(sessionInfo()), file.path(output_root, "sessionInfo.txt"))
  all_points <- list()
  all_crossings <- list()

  for (entry in entries) {
    points <- list()
    for (N in coarse_N) {
      points[[length(points) + 1L]] <- run_scp_point(entry, N, nsims, cores, G, output_root)
      if (tail(points, 1)[[1]]$power >= 0.85) break
    }

    repeat {
      combined <- do.call(rbind, points)
      bracket <- scp_n80_bracket(combined)
      if (is.null(bracket$high)) {
        if (bracket$low$nsims < nsims_final) {
          N <- bracket$low$N
        } else {
          larger_N <- coarse_N[coarse_N > max(combined$N)]
          if (!length(larger_N)) break
          N <- larger_N[1]
        }
      } else if (is.null(bracket$low)) {
        N <- max(6L, as.integer(floor(bracket$high$N / 2)))
        if (N >= bracket$high$N) break
      } else {
        low <- bracket$low$N
        high <- bracket$high$N
        needs_final <- c(low, high)[c(bracket$low$nsims, bracket$high$nsims) < nsims_final]
        if (length(needs_final)) {
          N <- needs_final[1]
        } else {
        if (high - low <= max(2, 0.025 * mean(c(low, high)))) break
        N <- as.integer(round(mean(c(low, high))))
        if (N <= low || N >= high) break
        }
      }
      old_index <- which(vapply(points, function(x) x$N, integer(1)) == N)
      final_point <- run_scp_point(entry, N, nsims_final, cores, G, output_root)
      if (length(old_index)) points[[old_index]] <- final_point else
        points[[length(points) + 1L]] <- final_point
    }

    combined <- do.call(rbind, points)
    combined <- combined[order(combined$N), , drop = FALSE]
    crossing <- scp_n80_summary(combined, entry)
    crossing$nsims_initial <- nsims
    crossing$nsims_final <- nsims_final
    uncertainty <- scp_n80_bootstrap(combined, entry, output_root)
    crossing$N80_mc95_low <- uncertainty["low"]
    crossing$N80_mc95_high <- uncertainty["high"]
    crossing$N80_mc_bootstrap_valid_fraction <- uncertainty["valid_fraction"]
    upper_point <- combined[which(combined$N == crossing$N80_high), , drop = FALSE]
    crossing$N80_FDR_N <- crossing$N80_high
    crossing$N80_FDR <- if (nrow(upper_point)) upper_point$unconditional_FDR else NA_real_
    crossing$N80_FDR_mcse <- if (nrow(upper_point)) upper_point$FDR_mcse else NA_real_
    all_points[[scp_case_id(entry)]] <- combined
    all_crossings[[scp_case_id(entry)]] <- crossing
    write.csv(combined, file.path(output_root, paste0(scp_case_id(entry), "_power.csv")),
                row.names = FALSE)
    write.csv(do.call(rbind, all_points), file.path(output_root, "power_summary.csv"),
                row.names = FALSE)
    write.csv(do.call(rbind, all_crossings), file.path(output_root, "N80_summary.csv"),
                row.names = FALSE)
  }
  cat("Simulation complete. N80 estimates are conditional on the fitted pilots.\n")
  cat("N80 MC intervals exclude pilot uncertainty and interpolation-grid error.\n")
}
