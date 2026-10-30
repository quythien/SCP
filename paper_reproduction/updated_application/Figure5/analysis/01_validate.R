# Freeze exact Figure 2.1 and two-harmonic calibrations, then validate all B.

# 1. Verify the existing scientific sources without modifying them ---------
root <- normalizePath(".")
out <- file.path(root, "dissertation_figures/figure_2_5/official_panels_AB_20261029")
source(file.path(out, "00_exact_active_engine.R"))
stopifnot(!file.exists(file.path(out, "validation_receipt.json")))
contract <- jsonlite::read_json("dissertation_figures/analysis_contract.json", simplifyVector = TRUE)
.libPaths(c(contract$package_library, .libPaths()))
library(SCP)
stopifnot(utils::packageDescription("SCP")$RemoteSha == contract$merged_commit)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) RhpcBLASctl::blas_set_num_threads(1L)
paths <- c(K1_calibration = paste0("dissertation_figures/figure_2_1/caudate_liver_candidate/",
                                  "v10_comparison/calibration_and_sampler.rds"),
           K1_settings = paste0("dissertation_figures/figure_2_1/caudate_liver_candidate/",
                               "v10_comparison/settings.json"),
           K1_original_input = file.path(contract$stage, "private/report_pilots.rds"),
           K1_helper = "examples/pilots/gtex_v10/05_compare_sample_sizes.R",
           K2_calibration = "dissertation_figures/figure_2_5/official_matched_20261029/calibration_and_sampler.rds",
           K2_pilot = "output/two_harmonic/results/pilot_2h_GTExLiver_topK300.rds",
           original_style = "examples/publication/onecol/fig6_onecol.R")
hashes <- vapply(paths, active_sha256, character(1))
K1_settings <- jsonlite::read_json(paths[["K1_settings"]], simplifyVector = TRUE)
stopifnot(hashes[["K1_helper"]] == K1_settings$input_sha256$simulation_helper,
          hashes[["K1_original_input"]] == K1_settings$input_sha256$prior_Liver_input,
          hashes[["K2_calibration"]] == "5d86ea5aef6b0a05a95ab0aebf8bd35fc9d0056a03e5afb0358d77c9cae328cd",
          hashes[["K2_pilot"]] == "1a7f2873a938e425585b3a900d49789d09b9db67eb9978e6ab2d32fef06734cb")
K1_case <- readRDS(paths[["K1_calibration"]])$Liver
K1_original <- readRDS(paths[["K1_original_input"]])$Liver_TOD$bio
K2 <- readRDS(paths[["K2_calibration"]])$bio
for (field in c("amplitude", "sigma_rhythmic", "phase", "prop_rhythmic", "period"))
  stopifnot(identical(K1_case$bio[[field]], K1_original[[field]]))
calibrations <- list(A = K1_case$bio, B = K2)
stopifnot(length(calibrations$A$amplitude) == 300L,
          all(lengths(K2[c("amplitude", "amplitude2", "phase", "phase2", "sigma_rhythmic")]) == 300L))
saveRDS(calibrations, file.path(out, "calibrations.rds"), compress = "xz")
write.csv(data.frame(role = names(paths), path = unname(paths), sha256 = unname(hashes)),
          file.path(out, "source_provenance.csv"), row.names = FALSE)

# 2. Freeze the balanced schedules, independent panel streams and job plan --
settings <- list(title = "GTEx Liver", panels = c("A", "B"),
  truth_A = "Pure first harmonic, exact Figure 2.1 Liver calibration",
  truth_B = "Two-harmonic, verified joint five-way Liver calibration",
  G = 5000L, nsims = 600L, N_grid = seq(24L, 600L, 24L),
  B = list(A = active_B("A"), B = active_B("B")),
  FDR = c(.01, .05, .10, .20), primary_FDR = .05,
  seed_base = c(A = 20261027L, B = 120261027L),
  seed_rule = "panel seed base + N * 100000 + replicate",
  CRN = "Parameters, canonical normal coordinates and residual chi-square shared across B within each panel/N/replicate",
  panels_independent = TRUE, A_phase_draw = "Independent empirical phase draws as in Figure 2.1",
  B_phase_draw = "Joint A1,A2,phase1,phase2,sigma source tuple",
  expected_information = "X'X = diag(N, N/2, ..., N/2) for every requested B",
  expected_curve_overlap = TRUE, horizontal_or_vertical_plot_jitter = FALSE,
  uncertainty = "Fixed-pilot Monte Carlo error; no source-pilot uncertainty")
saveRDS(settings, file.path(out, "settings.rds"), compress = "xz")
jobs <- expand.grid(chunk = 1:10, N = settings$N_grid, panel = settings$panels,
                    stringsAsFactors = FALSE)
jobs$job_id <- seq_len(nrow(jobs))
jobs$replicate_start <- (jobs$chunk - 1L) * 60L + 1L
jobs$replicate_end <- jobs$chunk * 60L
write.csv(jobs, file.path(out, "job_manifest.csv"), row.names = FALSE)
information <- list()
for (panel in settings$panels) for (N in settings$N_grid) for (B in active_B(panel)) {
  K <- if (panel == "A") 1L else 2L
  geometry <- active_geometry(N, B, K)
  information[[paste(panel, N, B)]] <- data.frame(panel = panel, N = N, B = B,
    m = N / B, K = K, rank = geometry$rank,
    information_max_error = geometry$information_max_error)
}
write.csv(do.call(rbind, information), file.path(out, "balanced_information_checks.csv"), row.names = FALSE)

# 3. Check the original parameter laws and actual expression for every B ---
Sys.setenv(SCP_SIM_LIBRARY = contract$package_library, SCP_SIM_DEFINE_ONLY = "true")
source(paths[["K1_helper"]])
Sys.unsetenv("SCP_SIM_DEFINE_ONLY")
detector <- getFromNamespace("detect_cosinor", "SCP")
checks <- fixtures <- list()
case <- 0L
for (panel in settings$panels) for (N in c(24L, 48L)) {
  seed <- as.integer(settings$seed_base[[panel]] + N * 100000L + 91L)
  value <- simulate_active_panel(calibrations[[panel]], panel, N, seed, G = 256L, keep_fixture = TRUE)
  fixture <- value$fixture
  p <- fixture$parameters
  K <- fixture$K
  if (panel == "A") {
    source_bio <- calibrations$A
    source_bio$ngenes <- 256L
    source_draw <- simulate_scp_f_replicate(source_bio, NULL, N, seed,
      retain_parameters = TRUE, time_sampler = K1_case$sampler)
    stopifnot(identical(p$truth, source_draw$truth),
              identical(p$A1[p$target_id], source_draw$amplitude),
              identical(p$sigma[p$target_id], source_draw$sigma),
              identical(p$phi1[p$target_id], source_draw$phase), all(p$A2 == 0))
  }
  for (B in active_B(panel)) {
    case <- case + 1L
    geometry <- active_geometry(N, B, K)
    X <- geometry$X
    if (panel == "B") {
      source_bio <- calibrations$B
      source_bio$ngenes <- 256L
      assign(".Random.seed", fixture$parameter_stream, envir = .GlobalEnv)
      oracle <- simCircadianSingleCohort2H(source_bio, geometry$times, seed = NULL)
      stopifnot(identical(p$truth, oracle$is_rhythmic), identical(p$A1, oracle$A1),
        identical(p$A2, oracle$A2), identical(p$phi1, oracle$phi1), identical(p$phi2, oracle$phi2))
    }
    H <- X[, -1, drop = FALSE] / sqrt(N / 2)
    mean_expression <- active_coefficients(p, K) %*% t(X[, -1, drop = FALSE]) + p$mesor
    full_Q <- qr.Q(geometry$decomposition, complete = TRUE)
    set.seed(930000L + case, kind = "L'Ecuyer-CMRG")
    df <- N - 2L * K - 1L
    direction <- matrix(rnorm(256L * df), nrow = 256L)
    residual <- direction * sqrt(fixture$residual_ss / rowSums(direction^2))
    standardized <- mean_expression / p$sigma +
      (fixture$coordinates - fixture$means) %*% t(H) +
      residual %*% t(full_Q[, (2L * K + 2L):N, drop = FALSE]) +
      rnorm(256L) %o% rep(1 / sqrt(N), N)
    expression <- standardized * p$sigma
    full_fit <- lm.fit(X, t(expression))
    rss <- colSums(full_fit$residuals^2)
    rss0 <- rowSums((expression - rowMeans(expression))^2)
    regression_p <- pf(((rss0 - rss) / (2L * K)) / (rss / df),
                        2L * K, df, lower.tail = FALSE)
    installed_p <- detector(expression, geometry$times, K = K, ebayes = FALSE)
    regression_error <- max(abs(fixture$pvalue - regression_p))
    installed_error <- max(abs(fixture$pvalue - installed_p))
    means_error <- max(abs((mean_expression / p$sigma) %*% H - fixture$means))
    stopifnot(regression_error < 1e-8, installed_error < 1e-8, means_error < 1e-10)
    expected_energy <- N / 2 * ((p$A1 / p$sigma)^2 +
      if (K == 2L) (p$A2 / p$sigma)^2 else 0)
    energy_error <- max(abs(rowSums(fixture$means^2) - expected_energy))
    stopifnot(energy_error < 1e-10)
    checks[[case]] <- data.frame(panel = panel, K = K, N = N, B = B, m = N / B,
      G = 256L, regression_pvalue_max_error = regression_error,
      installed_pvalue_max_error = installed_error, means_max_error = means_error,
      energy_max_error = energy_error)
    fixtures[[case]] <- list(panel = panel, K = K, N = N, B = B, seed = seed,
      times = geometry$times, expression = expression, truth = p$truth,
      p_exact = fixture$pvalue, p_regression = regression_p, p_installed = installed_p)
  }
  by_alpha <- split(value$metrics, value$metrics$FDR_target)
  stopifnot(all(vapply(by_alpha, function(x) length(unique(x$power)) == 1L &&
                       length(unique(x$FDP)) == 1L, logical(1))))
}
checks <- do.call(rbind, checks)
write.csv(checks, file.path(out, "expression_oracle_checks.csv"), row.names = FALSE)
saveRDS(fixtures, file.path(out, "synthetic_expression_oracles.rds"), compress = "xz")

# 4. Check failure cases, exact replay and small full-size timing batches --
stopifnot(inherits(try(active_geometry(24, 4, 2), silent = TRUE), "try-error"),
          inherits(try(active_geometry(25, 6, 1), silent = TRUE), "try-error"))
replay1 <- simulate_active_panel(calibrations$A, "A", 24L, 20261029L, 256L, TRUE)
replay2 <- simulate_active_panel(calibrations$A, "A", 24L, 20261029L, 256L, TRUE)
stopifnot(identical(replay1, replay2))
benchmark <- timing <- list()
for (panel in settings$panels) for (N in c(240L, 600L)) {
  start <- proc.time()[["elapsed"]]
  key <- paste(panel, N, sep = "_")
  benchmark[[key]] <- lapply(1:20, function(i) simulate_active_panel(calibrations[[panel]],
    panel, N, as.integer(settings$seed_base[[panel]] + N * 100L + i), 5000L))
  elapsed <- proc.time()[["elapsed"]] - start
  timing[[key]] <- data.frame(panel = panel, N = N, panel_replicates = 20L,
    B_cells = length(active_B(panel)), seconds = elapsed, seconds_per_panel_replicate = elapsed / 20)
}
saveRDS(benchmark, file.path(out, "validation_benchmark.rds"), compress = "xz")
timing <- do.call(rbind, timing)
write.csv(timing, file.path(out, "benchmark_timing.csv"), row.names = FALSE)

# 5. Record evidence for the independent gate and explain expected overlap -
receipt <- list(completed_UTC = format(Sys.time(), tz = "UTC", usetz = TRUE),
  status = "VALIDATED_REQUIRES_INDEPENDENT_ROOT_GATE", expression_oracle_cases = nrow(checks),
  oracle_N = c(24L, 48L), all_requested_B_values_tested = TRUE,
  maximum_pvalue_error = max(checks$regression_pvalue_max_error, checks$installed_pvalue_max_error),
  maximum_information_error = max(vapply(information, function(x) x$information_max_error, numeric(1))),
  K1_Figure2_1_parameters_replay_exactly = TRUE,
  K2_installed_generator_parameters_replay_exactly = TRUE,
  B_information_identical_mathematically = TRUE, B_results_identical_under_CRN = TRUE,
  curve_separation_expected = FALSE, no_plot_jitter = TRUE,
  K1_true_second_harmonic_amplitude_zero = TRUE, source_parameter_floors_preserved = TRUE,
  nonidentifiable_B4_K2_refused = TRUE, unbalanced_N_refused = TRUE,
  deterministic_replay_verified = TRUE, G = 5000L, nsims_per_cell = 600L,
  planned_N_values = length(settings$N_grid), planned_panel_replicates = 30000L,
  planned_B_cell_replicates = 150000L, planned_jobs = nrow(jobs),
  estimated_serial_compute_seconds = 15000 * sum(tapply(timing$seconds_per_panel_replicate,
                                                       timing$panel, mean)),
  calibration_sha256 = active_sha256(file.path(out, "calibrations.rds")),
  source_hashes = as.list(hashes), installed_package_SHA = contract$merged_commit,
  full_production_launched = FALSE, current_dissertation_modified = FALSE,
  individual_real_TOD_or_donor_identifiers_written = FALSE)
jsonlite::write_json(receipt, file.path(out, "validation_receipt.json"),
                     pretty = TRUE, auto_unbox = TRUE, digits = 17)
writeLines(capture.output(sessionInfo()), file.path(out, "validation_sessionInfo.txt"))
cat(sprintf("Validated %d actual-expression cases; max p-value error %.3g.\n",
            nrow(checks), receipt$maximum_pvalue_error))
cat(sprintf("Expected B overlap confirmed. Estimated serial compute: %.1f seconds.\n",
            receipt$estimated_serial_compute_seconds))
