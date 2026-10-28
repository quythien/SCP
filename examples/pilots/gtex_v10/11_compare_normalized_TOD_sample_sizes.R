# Estimate N80 from the normalized-expression GTEx v10 TOD pilots.

# 1. Load the validated simulation helpers without running the old analysis -
stage <- Sys.getenv("SCP_GTEX_TOD_STAGE", "output/gtex_v10_normalized_TOD")
publication <- Sys.getenv("SCP_GTEX_TOD_PUBLICATION", stage)
pilot_file <- file.path(stage, "private", "report_pilots.rds")
signal_file <- file.path(publication, "summary", "gtex_v10_signal_summary.csv")
manifest_file <- file.path(publication, "inst", "extdata", "pilots", "manifest.csv")
output_root <- file.path(stage, "simulation", "results_exactF")
report_root <- file.path(stage, "report")
stopifnot(file.exists(pilot_file), file.exists(signal_file), file.exists(manifest_file),
          requireNamespace("digest", quietly = TRUE), requireNamespace("jsonlite", quietly = TRUE))
engine_file <- "examples/pilots/gtex_v10/05_compare_sample_sizes.R"
engine_hash <- digest::digest(file = engine_file, algo = "sha256")
stopifnot(identical(engine_hash,
  "0a0e204d2ea0027a4d4afebf6de1d9374a1c1e312273e3762110b43e78a7e0c6"))
old_define_only <- Sys.getenv("SCP_SIM_DEFINE_ONLY", unset = NA_character_)
Sys.setenv(SCP_SIM_DEFINE_ONLY = "true", SCP_GTEX_V10_STAGE = stage,
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
source(engine_file)
if (is.na(old_define_only)) Sys.unsetenv("SCP_SIM_DEFINE_ONLY") else
  Sys.setenv(SCP_SIM_DEFINE_ONLY = old_define_only)
if (requireNamespace("RhpcBLASctl", quietly = TRUE))
  RhpcBLASctl::blas_set_num_threads(1L)
stopifnot(utils::packageVersion("SCP") >= "0.5.1")

# 2. Check the ten calibrations against the newly staged full-gene fits -----
entries <- readRDS(pilot_file)
if (!is.null(entries$pilots)) entries <- entries$pilots
expected_tissues <- c("Adrenal", "Liver", "Muscle", "Thyroid", "Lung", "Pancreas",
  "Heart_Left_Ventricle", "Artery_Tibial", "Adipose_Subcutaneous", "Brain_Cortex")
stopifnot(length(entries) == 10L,
          setequal(vapply(entries, function(x) x$tissue, character(1)), expected_tissues),
          all(vapply(entries, function(x) identical(x$time_basis, "TOD"), logical(1))))
entries <- entries[match(expected_tissues, vapply(entries, function(x) x$tissue, character(1)))]
manifest <- utils::read.csv(manifest_file, stringsAsFactors = FALSE)
manifest <- manifest[manifest$dataset == "GTEx", , drop = FALSE]
signal <- utils::read.csv(signal_file, stringsAsFactors = FALSE)
stopifnot(nrow(manifest) == 46L, nrow(signal) == 46L,
          !anyDuplicated(manifest$tissue), !anyDuplicated(signal$tissue))
near <- function(x, y) isTRUE(all.equal(as.numeric(x), as.numeric(y),
  tolerance = 1e-12, check.attributes = FALSE))
validation <- list()
for (entry in entries) {
  tissue <- if (entry$tissue == "Muscle") "Muscle_Skeletal" else entry$tissue
  row <- manifest[match(tissue, manifest$tissue), , drop = FALSE]
  summary <- signal[match(tissue, signal$tissue), , drop = FALSE]
  path <- file.path(publication, "inst", "extdata", "pilots", row$file)
  pilot <- readRDS(path)
  fit <- pilot$explorer_fit
  stopifnot(!is.null(fit), nrow(fit) == entry$genes_fitted,
            !anyDuplicated(fit$gene), all(is.finite(fit$pvalue)),
            near(stats::p.adjust(fit$pvalue, "BH"), fit$qvalue))
  passing <- which(fit$pvalue < 0.01 & is.finite(fit$A) & fit$A > 0 &
    is.finite(fit$sigma) & fit$sigma > 0 & is.finite(fit$phi))
  selected <- head(passing[order(fit$pvalue[passing])], 300L)
  ratio <- fit$A[selected] / fit$sigma[selected]
  bio <- entry$bio
  checks <- c(release = identical(pilot$gtex_release, "v10"),
    TOD = identical(pilot$time_basis, "TOD"),
    n = near(pilot$n_samples, entry$n_pilot) && near(summary$n, entry$n_pilot),
    genes = near(pilot$ngenes, entry$genes_fitted) && near(summary$ngenes, entry$genes_fitted),
    selection = near(length(selected), entry$n_selected) && near(summary$n_selected, entry$n_selected),
    amplitude = near(bio$amplitude, fit$A[selected]),
    noise = near(bio$sigma_rhythmic, fit$sigma[selected]),
    phase = near(bio$phase, fit$phi[selected]),
    median_r = near(stats::median(ratio), entry$medianr) && near(summary$r_median_top300, entry$medianr),
    screening_fraction = near(length(passing) / nrow(fit), entry$prop) &&
      near(bio$prop_rhythmic, entry$prop) && near(summary$prop_rhythmic, entry$prop),
    BH05 = near(sum(fit$qvalue < 0.05), summary$n_BH05),
    paired_sigma = isTRUE(bio$paired_sigma),
    nominal_threshold = near(pilot$alpha_pilot, 0.01),
    top_k = near(pilot$pilot_top_k, 300L),
    time_vector = length(entry$times) == entry$n_pilot && all(is.finite(entry$times)) &&
      all(entry$times >= 0 & entry$times < 24),
    normalized_baselines = stats::median(bio$lBaselineExpr) > 3)
  stopifnot(all(checks), length(selected) > 0L,
            length(bio$amplitude) == length(bio$sigma_rhythmic))
  validation[[entry$tissue]] <- data.frame(tissue = entry$tissue,
    n = entry$n_pilot, genes = entry$genes_fitted, n_selected = length(selected),
    median_r = stats::median(ratio), prop_rhythmic = entry$prop,
    n_BH05 = sum(fit$qvalue < 0.05), median_baseline = stats::median(bio$lBaselineExpr),
    checks_passed = sum(checks), pilot_sha256 = digest::digest(file = path, algo = "sha256"))
  if (entry$tissue == "Liver") {
    second_path <- file.path(publication, "inst", "extdata", "pilots", sub("[.]rds$", "_2H.rds", row$file))
    second <- readRDS(second_path)$explorer_fit
    second_index <- match(fit$gene, second$gene)
    stopifnot(!anyNA(second_index), entry$n_pilot == 262L, entry$genes_fitted == 8917L,
      sum(fit$qvalue < .05) == 210L, sum(second$qvalue < .05) == 725L,
      sum(fit$qvalue < .05 & second$qvalue[second_index] < .05) == 207L)
    reference_file <- Sys.getenv("SCP_ORIGINAL_LIVER_FITS",
      "output/two_harmonic/results/tissue_Liver_fits.rds")
    if (file.exists(reference_file)) {
      reference <- readRDS(reference_file)
      reference_index <- match(fit$symbol, reference$gene)
      stopifnot(!anyNA(reference_index),
        max(abs(fit$qvalue - reference$q1H[reference_index])) < 1e-9,
        max(abs(second$qvalue[second_index] - reference$q2H[reference_index])) < 1e-9)
    }
  }
}
validation <- do.call(rbind, validation)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
validation_path <- file.path(output_root, "normalized_fit_validation.csv")
if (file.exists(validation_path)) {
  previous <- utils::read.csv(validation_path, stringsAsFactors = FALSE)
  stopifnot(isTRUE(all.equal(previous, validation, check.attributes = FALSE)))
} else utils::write.csv(validation, validation_path, row.names = FALSE)
cat("Ten normalized TOD calibrations and the exact original Liver result passed validation.\n")

# 3. Fix the Monte Carlo settings and record the new input identities -------
nsims <- max(200L, as.integer(Sys.getenv("SCP_NSIM", "200")))
nsims_final <- max(nsims, 1000L, as.integer(Sys.getenv("SCP_NSIM_FINAL", "1000")))
cores <- min(2L, max(1L, as.integer(Sys.getenv("SCP_CORES", "2"))))
G <- 5000L
max_N <- 10000L
coarse_N <- c(12L, 24L, 48L, 96L, 192L, 384L, 768L, 1536L, 3072L, 6144L, 10000L)
stopifnot(nsims_final < 100000L)
run_settings <- list(engine = "Exact Gaussian F-statistic Monte Carlo",
  expression_validation_engine = "actual_expression_cpp", G = G,
  nsims_initial = nsims, nsims_final = nsims_final, cores = cores,
  coarse_N = coarse_N, maximum_N = max_N, BH_threshold = 0.05,
  matched_cases = length(entries), cohort = "All valid season-adjusted TOD samples",
  expression_preprocessing = "TMM-normalized log-CPM without covariate regression",
  gene_filter = "First valid annotation symbol, then mean normalized log-CPM>3",
  alpha_pilot = 0.01, top_k = 300L, paired_sigma = TRUE,
  shared_seeds_across_cases = TRUE,
  seed_rule = "20261027 + N * 100000 + replicate_id",
  RNG = "L'Ecuyer-CMRG with parameter, time, and statistic substreams",
  bootstrap_replicates = 1000L,
  pilot_file_sha256 = digest::digest(file = pilot_file, algo = "sha256"),
  signal_file_sha256 = digest::digest(file = signal_file, algo = "sha256"),
  manifest_file_sha256 = digest::digest(file = manifest_file, algo = "sha256"),
  engine_sha256 = engine_hash,
  builder_sha256 = digest::digest(file = "examples/pilots/gtex_v10/10_build_normalized_TOD.R", algo = "sha256"),
  generator_sha256 = digest::digest(file = "examples/pilots/gtex_v10/11_compare_normalized_TOD_sample_sizes.R", algo = "sha256"),
  MC_interval_scope = "Fixed pilots. Excludes pilot uncertainty and interpolation-grid error.")
settings_path <- file.path(output_root, "run_settings.rds")
if (file.exists(settings_path)) {
  stopifnot(identical(readRDS(settings_path), run_settings))
} else {
  saveRDS(run_settings, settings_path)
  jsonlite::write_json(run_settings, file.path(output_root, "run_settings.json"),
    pretty = TRUE, auto_unbox = TRUE)
  writeLines(capture.output(sessionInfo()), file.path(output_root, "sessionInfo.txt"))
}

# 4. Run the coarse sample-size grid and refine the 80 percent crossing -----
if (!identical(Sys.getenv("SCP_NORMALIZED_VALIDATE_ONLY"), "true")) {
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
    upper <- combined[which(combined$N == crossing$N80_high), , drop = FALSE]
    crossing$N80_FDR_N <- crossing$N80_high
    crossing$N80_FDR <- if (nrow(upper)) upper$unconditional_FDR else NA_real_
    crossing$N80_FDR_mcse <- if (nrow(upper)) upper$FDR_mcse else NA_real_
    crossing$gtex_release <- "v10"
    crossing$genes_fitted <- entry$genes_fitted
    crossing$n_selected <- entry$n_selected
    crossing$median_r <- entry$medianr
    crossing$prop_rhythmic <- entry$prop
    crossing$calibration_law <- entry$calibration_law
    all_points[[scp_case_id(entry)]] <- combined
    all_crossings[[scp_case_id(entry)]] <- crossing
    utils::write.csv(combined, file.path(output_root,
      paste0(scp_case_id(entry), "_power.csv")), row.names = FALSE)
    utils::write.csv(do.call(rbind, all_points), file.path(output_root,
      "power_summary.csv"), row.names = FALSE)
    utils::write.csv(do.call(rbind, all_crossings), file.path(output_root,
      "N80_summary.csv"), row.names = FALSE)
  }

# 5. Save the ten-case table for the new normalized TOD report --------------
  comparison <- do.call(rbind, all_crossings)
  stopifnot(nrow(comparison) == 10L, !anyDuplicated(comparison$case_id))
  dir.create(report_root, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(comparison, file.path(report_root, "N80_comparison_table.csv"), row.names = FALSE)
  cat("New normalized TOD simulations complete for ten tissues.\n")
  cat("N80 Monte Carlo intervals are conditional on the fitted pilots.\n")
}
