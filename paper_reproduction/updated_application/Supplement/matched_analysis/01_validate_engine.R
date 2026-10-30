# Validate the frozen-calibration, matched-truth comparison locally.

# 1. Read the original calibration and check its provenance ----------------
root <- normalizePath(".")
out <- file.path(root, "dissertation_figures/figure_2_5/official_matched_20261029")
prior <- file.path(root, "dissertation_figures/figure_2_5/matched_v10")
stopifnot(!file.exists(file.path(out, "validation_receipt.json")))
source(file.path(out, "00_paired_engine.R"))
contract <- jsonlite::read_json(file.path(root, "dissertation_figures/analysis_contract.json"),
                              simplifyVector = TRUE)
.libPaths(c(contract$package_library, .libPaths()))
library(SCP)
stopifnot(utils::packageDescription("SCP")$RemoteSha == contract$merged_commit)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) RhpcBLASctl::blas_set_num_threads(1L)
old <- jsonlite::read_json(file.path(prior, "validation_receipt.json"), simplifyVector = TRUE)
source_paths <- c(source_pilot = old$source_pilot,
                 source_sampler = old$passive_sampler_source,
                 prior_engine = file.path(prior, "00_paired_engine.R"),
                 prior_validation = file.path(prior, "01_validate_engine.R"),
                 prior_calibration = file.path(prior, "calibration_and_sampler.rds"))
source_hashes <- vapply(source_paths, matched_sha256, character(1))
stopifnot(source_hashes[["source_pilot"]] == old$source_pilot_sha256,
          source_hashes[["source_sampler"]] == old$passive_sampler_source_sha256,
          source_hashes[["prior_engine"]] == old$engine_sha256,
          source_hashes[["prior_validation"]] == old$validation_script_sha256,
          matched_sha256(file.path(out, "calibration_and_sampler.rds")) ==
            source_hashes[["prior_calibration"]])
input <- readRDS(file.path(out, "calibration_and_sampler.rds"))
bio <- input$bio
sampler <- input$sampler
stopifnot(lengths(bio[c("amplitude", "amplitude2", "phase", "phase2", "sigma_rhythmic")]) == 300L,
          bio$ngenes == 8917L, bio$period == 24, all(bio$phase2 >= 0 & bio$phase2 < 12),
          all(is.finite(sampler$x)), all(is.finite(sampler$probability)),
          all(sampler$probability >= 0), sum(sampler$probability) > 0)
write.csv(data.frame(role = names(source_paths), path = unname(source_paths),
                    sha256 = unname(source_hashes)),
          file.path(out, "source_provenance.csv"), row.names = FALSE)

# 2. Compare actual expression with coordinate and installed tests --------
checks <- fixtures <- list()
case <- 0L
for (N in c(12L, 24L, 48L)) {
  seed <- 612900L + N
  pair <- simulate_matched_pair(bio, sampler, N, seed, G = 256L, keep_fixture = TRUE)
  stopifnot(identical(pair$fixture$designs$passive$parameters,
                      pair$fixture$designs$balanced6$parameters),
            !identical(pair$fixture$RNG_streams$passive_noise,
                       pair$fixture$RNG_streams$balanced6_noise))
  for (design in c("passive", "balanced6")) {
    case <- case + 1L
    f <- pair$fixture$designs[[design]]
    times <- f$times
    draws <- f$parameters
    projection <- matched_projection(times)
    assign(".Random.seed", pair$fixture$RNG_streams$parameters, envir = .GlobalEnv)
    oracle_bio <- bio
    oracle_bio$ngenes <- 256L
    oracle <- simCircadianSingleCohort2H(oracle_bio, times, seed = NULL)
    stopifnot(identical(draws$truth, oracle$is_rhythmic),
              identical(draws$A1, oracle$A1), identical(draws$A2, oracle$A2),
              identical(draws$phi1, oracle$phi1), identical(draws$phi2, oracle$phi2))
    y <- oracle$expr / draws$sigma
    coordinates <- y %*% projection$Q[, 2:5, drop = FALSE]
    full <- lm.fit(projection$X, t(y))
    rss2 <- colSums(full$residuals^2)
    from_coordinates <- matched_pvalues(coordinates, rss2, N)
    first <- lm.fit(projection$X[, 1:3], t(y))
    rss1 <- colSums(first$residuals^2)
    rss0 <- rowSums((y - rowMeans(y))^2)
    expression_p1 <- pf(((rss0 - rss1) / 2) / (rss1 / (N - 3)),
                        2, N - 3, lower.tail = FALSE)
    expression_p2 <- pf(((rss0 - rss2) / 4) / (rss2 / (N - 5)),
                        4, N - 5, lower.tail = FALSE)
    detector <- getFromNamespace("detect_cosinor", "SCP")
    installed_p1 <- detector(oracle$expr, times, K = 1L, ebayes = FALSE)
    installed_p2 <- detector(oracle$expr, times, K = 2L, ebayes = FALSE)
    omega <- 2 * pi / 24
    coefficients <- cbind(draws$A1 * cos(omega * draws$phi1),
      draws$A1 * sin(omega * draws$phi1), draws$A2 * cos(2 * omega * draws$phi2),
      draws$A2 * sin(2 * omega * draws$phi2))
    mean_expression <- coefficients %*% t(projection$X[, -1]) + draws$mesor
    means <- matched_mean_coordinates(draws, projection)
    mean_error <- max(abs((mean_expression / draws$sigma) %*%
                          projection$Q[, 2:5] - means))

    # Reconstruct a full Gaussian expression matrix from this engine's draws.
    # Its residual direction is independent of the engine's chi-square radius.
    set.seed(991000L + case, kind = "L'Ecuyer-CMRG")
    Qfull <- qr.Q(projection$decomposition, complete = TRUE)
    direction <- matrix(rnorm(256L * (N - 5L)), nrow = 256L)
    residual <- direction * sqrt(f$residual_ss / rowSums(direction^2))
    standardized <- mean_expression / draws$sigma +
      (f$coordinates - means) %*% t(Qfull[, 2:5]) +
      residual %*% t(Qfull[, 6:N, drop = FALSE]) +
      rnorm(256L) %o% Qfull[, 1]
    reconstructed <- standardized * draws$sigma
    reconstructed_p1 <- detector(reconstructed, times, K = 1L, ebayes = FALSE)
    reconstructed_p2 <- detector(reconstructed, times, K = 2L, ebayes = FALSE)
    errors <- c(max(abs(from_coordinates$p1 - expression_p1)),
      max(abs(from_coordinates$p2 - expression_p2)),
      max(abs(from_coordinates$p1 - installed_p1)),
      max(abs(from_coordinates$p2 - installed_p2)),
      max(abs(f$pvalue_K1 - reconstructed_p1)),
      max(abs(f$pvalue_K2 - reconstructed_p2)))
    stopifnot(max(errors) < 1e-8, mean_error < 1e-10)
    if (design == "balanced6") stopifnot(identical(sort(unique(times)), seq(0, 20, by = 4)),
      all(table(times) == N / 6L))
    checks[[case]] <- data.frame(design = design, N = N, G = 256L,
      coordinate_vs_expression_K1 = errors[1], coordinate_vs_expression_K2 = errors[2],
      coordinate_vs_installed_K1 = errors[3], coordinate_vs_installed_K2 = errors[4],
      engine_vs_reconstructed_expression_K1 = errors[5],
      engine_vs_reconstructed_expression_K2 = errors[6],
      projected_mean_max_error = mean_error, seed = seed)
    fixtures[[case]] <- list(design = design, N = N, seed = seed, times = times,
      source_generator_expression = oracle$expr, engine_expression = reconstructed,
      truth = draws$truth, engine_p1 = f$pvalue_K1, engine_p2 = f$pvalue_K2)
  }
}
checks <- do.call(rbind, checks)
write.csv(checks, file.path(out, "expression_oracle_checks.csv"), row.names = FALSE)
saveRDS(fixtures, file.path(out, "synthetic_expression_oracle_fixtures.rds"), compress = "xz")

# 3. Verify balanced-six energies, passive replay and deterministic replay -
old_api <- new.env(parent = globalenv())
sys.source(file.path(prior, "00_paired_engine.R"), envir = old_api)
energy_errors <- list()
for (N in c(12L, 24L, 48L, 240L, 1200L)) {
  times <- matched_schedule(N, "balanced6", sampler)
  projection <- matched_projection(times)
  draws <- list(A1 = bio$amplitude, A2 = bio$amplitude2, phi1 = bio$phase,
                phi2 = bio$phase2, sigma = bio$sigma_rhythmic)
  means <- matched_mean_coordinates(draws, projection)
  expected1 <- N / 2 * (draws$A1 / draws$sigma)^2
  expected2 <- N / 2 * (draws$A2 / draws$sigma)^2
  error <- max(abs(rowSums(means[, 1:2]^2) - expected1),
               abs(rowSums(means[, 3:4]^2) - expected2))
  stopifnot(error < 1e-10)
  energy_errors[[as.character(N)]] <- error
  set.seed(3200L + N, kind = "L'Ecuyer-CMRG")
  previous_times <- old_api$matched_schedule(N, "passive", sampler)
  set.seed(3200L + N, kind = "L'Ecuyer-CMRG")
  stopifnot(identical(previous_times, matched_schedule(N, "passive", sampler)))
}
stopifnot(inherits(try(matched_schedule(20L, "balanced6", sampler), silent = TRUE), "try-error"))
a <- simulate_matched_pair(bio, sampler, 24L, 20261029L, 256L, TRUE)
b <- simulate_matched_pair(bio, sampler, 24L, 20261029L, 256L, TRUE)
stopifnot(identical(a, b))

# 4. Freeze the requested first-pass job plan and benchmark small batches --
N_grid <- c(12L, 24L, 48L, 72L, 96L, 120L, 144L, 168L, 192L, 216L,
            240L, 264L, 288L, 312L, 336L, 360L, 408L, 480L, 600L, 900L, 1200L)
settings <- list(title = "GTEx Liver", application = "Appendix A matched comparison",
  generating_law = "Same source-calibrated two-harmonic Gaussian truth for all four design-detector combinations",
  not_Figure_2_5A_pure_K1_truth = TRUE, G = 5000L, nsims = 600L,
  N_grid = N_grid, FDR = c(.01, .05, .10, .20), primary_FDR = .05,
  designs = c("passive", "balanced6"), balanced_phases = seq(0, 20, by = 4),
  source_samples = 262L, source_fitted_genes = 8917L, selected_joint_tuples = 300L,
  source_prop_rhythmic = bio$prop_rhythmic,
  expected_targets_per_experiment = round(5000L * bio$prop_rhythmic),
  seed_base = 20261029L, seed_rule = "seed_base + N * 100000 + replicate",
  parameter_pairing = "Same truth, baseline draws and joint five-way rhythmic tuples across both designs",
  noise_pairing = "Independent Gaussian noise across designs, separate L'Ecuyer-CMRG streams",
  detector_pairing = "Same times, truth, mean, Gaussian coordinates and residual for K1 and K2 within each design",
  streams = c("parameters", "passive_times", "passive_noise", "balanced6_noise"),
  passive_law = "Original grid/weights draw with replacement, independent N(0,0.1^2) jitter, modulo 24",
  K1_denominator = "Second-harmonic coordinate sum of squares plus full-model residual sum of squares",
  uncertainty = "Fixed-pilot Monte Carlo error only", benchmark_reused_in_production = FALSE)
saveRDS(settings, file.path(out, "settings.rds"), compress = "xz")
jobs <- expand.grid(chunk = 1:10, N = N_grid)
jobs$job_id <- seq_len(nrow(jobs))
jobs$replicate_start <- (jobs$chunk - 1L) * 60L + 1L
jobs$replicate_end <- jobs$chunk * 60L
write.csv(jobs, file.path(out, "job_manifest.csv"), row.names = FALSE)
benchmark <- list()
timing <- list()
for (N in c(240L, 1200L)) {
  start <- proc.time()[["elapsed"]]
  benchmark[[as.character(N)]] <- lapply(1:20, function(i)
    simulate_matched_pair(bio, sampler, N, seed = 910000L + N * 100L + i, G = 5000L))
  elapsed <- proc.time()[["elapsed"]] - start
  timing[[as.character(N)]] <- data.frame(N = N, shared_parameter_pairs = 20L,
    datasets = 40L, elapsed_seconds = elapsed, seconds_per_pair = elapsed / 20)
}
saveRDS(benchmark, file.path(out, "validation_benchmark.rds"), compress = "xz")
write.csv(do.call(rbind, timing), file.path(out, "benchmark_timing.csv"), row.names = FALSE)

# 5. Record the evidence required for independent review ------------------
receipt <- list(completed_at_UTC = format(Sys.time(), tz = "UTC", usetz = TRUE),
  status = "VALIDATED_LOCALLY_REQUIRES_INDEPENDENT_GATE", oracle_cases = nrow(checks),
  oracle_N = c(12L, 24L, 48L), oracle_designs = c("passive", "balanced6"),
  maximum_pvalue_error = max(as.matrix(checks[, 4:9])),
  maximum_projected_mean_error = max(checks$projected_mean_max_error),
  balanced_six_energy_errors = energy_errors, passive_schedule_byte_replay = TRUE,
  deterministic_pair_replay = TRUE, same_truth_and_parameters_across_designs = TRUE,
  independent_noise_streams_across_designs = TRUE,
  K1_omitted_harmonic_in_denominator = TRUE, five_way_tuples_preserved = TRUE,
  A1_floor = .05, A2_floor = 0, sigma_floor = 1e-6,
  source_A1_below_floor = sum(bio$amplitude < .05),
  source_A2_below_floor = sum(bio$amplitude2 < 0),
  calibration_and_sampler_sha256 = matched_sha256(file.path(out, "calibration_and_sampler.rds")),
  source_hashes = as.list(source_hashes), installed_package_SHA = contract$merged_commit,
  code_hashes = as.list(vapply(file.path(out, c("00_paired_engine.R", "01_validate_engine.R",
    "02_run_jobs.R", "03_summarize.R")), matched_sha256, character(1))),
  planned_pairs = length(N_grid) * settings$nsims,
  planned_datasets = length(N_grid) * settings$nsims * 2L,
  planned_detector_evaluations = length(N_grid) * settings$nsims * 4L,
  benchmark_seconds_per_pair = vapply(timing, function(x) x$seconds_per_pair, numeric(1)),
  estimated_serial_seconds = length(N_grid) * settings$nsims *
    mean(vapply(timing, function(x) x$seconds_per_pair, numeric(1))),
  production_launched = FALSE, real_individual_TOD_or_donor_ids_written = FALSE,
  existing_paper_package_dissertation_files_modified = FALSE)
jsonlite::write_json(receipt, file.path(out, "validation_receipt.json"),
                     pretty = TRUE, auto_unbox = TRUE, digits = 17)
writeLines(capture.output(sessionInfo()), file.path(out, "validation_sessionInfo.txt"))
cat(sprintf("Six expression-oracle cases passed; maximum p-value error %.3g.\n", receipt$maximum_pvalue_error))
cat(sprintf("Estimated serial first-pass runtime %.1f seconds from 40 validation pairs.\n", receipt$estimated_serial_seconds))
