# Build a synthetic Caudate/Liver example from the reviewed aggregate pilot.
# Run: Rscript 01_build_caudate_liver.R CALIBRATION_RDS SCP_LIBRARY OUTPUT_RDS
# The output contains simulated expression and scheduled simulation times.

# 1. Read explicit paths and verify the frozen calibration and generator.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Supply CALIBRATION_RDS, SCP_LIBRARY and a new OUTPUT_RDS path.", call. = FALSE)
}
calibration_file <- normalizePath(args[1L], mustWork = TRUE)
library_path <- normalizePath(args[2L], mustWork = TRUE)
output_file <- args[3L]
stopifnot(!file.exists(output_file), dir.exists(dirname(output_file)))
.libPaths(c(library_path, .libPaths()))
suppressPackageStartupMessages(library(SCP))
stopifnot(normalizePath(find.package("SCP")) == normalizePath(file.path(library_path, "SCP")))
calibration_sha256 <- digest::digest(calibration_file, file = TRUE, algo = "sha256")
stopifnot(calibration_sha256 == "eeb5a55c2bae14105756c1a842877c563bf7af4dc4ad867205a33fcea4b9611c")
simulator <- getFromNamespace("simCircadianDiff", "SCP")
simulator_sha256 <- digest::digest(paste(deparse(body(simulator)), collapse = "\n"),
                                  algo = "sha256", serialize = FALSE)
stopifnot(simulator_sha256 == "9d5a7e635c5eba7a3ac95e28158a6e19997898c260d700752022eaf729ca4971")
compiled_library <- getLoadedDLLs()[["SCP"]][["path"]]
compiled_sha256 <- digest::digest(compiled_library, file = TRUE, algo = "sha256")
stopifnot("_SCP_sim_cosinor_expr_cpp" %in% names(getDLLRegisteredRoutines("SCP")$.Call))
# Select the same compiled backend used by the reviewed production generator.
.CPP_LOADED <- TRUE
stopifnot(isTRUE(get(".CPP_LOADED", envir = asNamespace("SCP"), inherits = TRUE)))

# 2. Check the aggregate cohort and selected amplitude/noise pairs.
pilot <- readRDS(calibration_file)
stopifnot(inherits(pilot, "CircadianBioOptions"), pilot$ngenes == 8831L,
          inherits(pilot$cts, "SCPTimeDistribution"), pilot$cts$n_observed == 59L,
          inherits(pilot$cts2, "SCPTimeDistribution"), pilot$cts2$n_observed == 262L,
          pilot$period == 24, isTRUE(pilot$paired_sigma2),
          length(pilot$amplitude) == length(pilot$sigma_rhythmic),
          length(pilot$amplitude2) == length(pilot$sigma_rhythmic2),
          all(is.finite(c(pilot$amplitude, pilot$amplitude2))),
          all(c(pilot$amplitude, pilot$amplitude2) >= 0),
          all(is.finite(c(pilot$sigma_rhythmic, pilot$sigma_rhythmic2))),
          all(c(pilot$sigma_rhythmic, pilot$sigma_rhythmic2) > 0))

# 3. Define a balanced schedule and a fixed seed for independent samples.
ngenes <- 1200L
n_per_group <- 48L
seed <- 20261029L
scheduled_times <- rep(seq(0, 20, by = 4), each = 8L)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1")
sim <- simulator(ngenes = ngenes, n1 = n_per_group, n2 = n_per_group,
  lBaselineExpr = pilot$lBaselineExpr, lBaselineExpr2 = pilot$lBaselineExpr2,
  lOD = pilot$lOD, lOD2 = pilot$lOD2,
  amplitude = pilot$amplitude, amplitude2 = pilot$amplitude2,
  sigma_rhythmic = pilot$sigma_rhythmic, sigma_rhythmic2 = pilot$sigma_rhythmic2,
  prop_rhythmic = pilot$prop_rhythmic, prop_DR = pilot$prop_DR,
  prop_DP = pilot$prop_DP, prop_DM = pilot$prop_DM,
  phase_diff = pilot$phase_diff, dp_shift_mode = pilot$dp_shift_mode,
  amp_diff = pilot$amp_diff, mesor_diff = pilot$mesor_diff,
  period = pilot$period, design = "active", cts = scheduled_times,
  cts2 = scheduled_times, harmonics = c(0, 0), sim.seed = seed)

# 4. Label every gene and sample as synthetic and verify the realized truth.
genes <- sprintf("SyntheticGene%04d", seq_len(ngenes))
rownames(sim$expr1) <- rownames(sim$expr2) <- genes
colnames(sim$expr1) <- sprintf("SyntheticCaudate_%03d", seq_len(n_per_group))
colnames(sim$expr2) <- sprintf("SyntheticLiver_%03d", seq_len(n_per_group))
truth <- sim$ground_truth
truth$gene_id <- genes
stopifnot(identical(dim(sim$expr1), c(ngenes, n_per_group)),
          identical(dim(sim$expr2), c(ngenes, n_per_group)),
          all(is.finite(sim$expr1)), all(is.finite(sim$expr2)),
          identical(sim$times1, scheduled_times), identical(sim$times2, scheduled_times),
          !anyDuplicated(genes),
          !anyDuplicated(c(colnames(sim$expr1), colnames(sim$expr2))),
          identical(truth$is_DR, truth$diff_type %in% c(2L, 3L)),
          identical(truth$is_DP, truth$diff_type == 4L),
          identical(truth$is_DM, truth$diff_type == 5L),
          all(truth$mesor1[!truth$is_DM] == truth$mesor2[!truth$is_DM]),
          all(truth$phase1 >= 0 & truth$phase1 < 24),
          all(truth$phase2 >= 0 & truth$phase2 < 24))

# 5. Save only synthetic data, truth and compact provenance without local paths.
source_note <- paste("Synthetic Caudate/Liver expression generated from the reviewed",
  "joint aggregate calibration. These are simulated samples and scheduled",
  "simulation times, not observed donors or a reconstruction of either cohort.")
fixture <- list(
  caudate = list(expr = sim$expr1, times = sim$times1,
                 synthetic = TRUE, time_units = "hours",
                 time_type = "scheduled simulation times"),
  liver = list(expr = sim$expr2, times = sim$times2,
               synthetic = TRUE, time_units = "hours",
               time_type = "scheduled simulation times"),
  source = source_note,
  metadata = list(
    synthetic = TRUE, schema_version = 1L,
    purpose = "Runnable examples of fitting synthetic two-group expression data.",
    interpretation = paste("Conditional calibration across studies. The synthetic arm",
      "labels do not establish a biological tissue effect or reproduce official power results."),
    expression_units = "Simulated log-scale expression; not observed expression or counts.",
    time_units = "hours", phase_units = "hours", period = 24,
    time_type = "scheduled simulation times", design = "balanced active",
    distinct_times = 6L, independent_samples_per_time = 8L,
    genes = ngenes, independent_samples_per_group = n_per_group,
    noise_model = "Independent Gaussian errors for every simulated gene and sample in both arms.",
    signal_model = "First-harmonic differential cosinor with the saved DR, DP and DM parameter law.",
    calibration = list(name = "Caudate control / GTEx v10 seasonal-TOD Liver joint aggregate",
      sha256 = calibration_sha256, observed_cohort_counts = c(Caudate_control = 59L, Liver = 262L),
      common_tested_genes = 8831L, raw_p_threshold = 0.01, top_k_per_arm = 300L),
    generator = list(package = "SCP", version = as.character(packageVersion("SCP")),
      function_name = "simCircadianDiff", function_body_sha256 = simulator_sha256,
      compiled_library_sha256 = compiled_sha256,
      reference_source_archive_sha256 = "0db6aa8f7ce813ee0cac719c12bac056133ff8285784403128a0f4433618e51a",
      backend = "compiled", R_version = as.character(getRversion()),
      platform = R.version$platform,
      seed = seed, RNGkind = RNGkind()),
    observed_individual_times_included = FALSE, donor_identifiers_included = FALSE,
    observed_expression_included = FALSE),
  ground_truth = truth,
  parameter_draws = sim$parameter_draws,
  simulation_options = sim$simOptions)
saveRDS(fixture, output_file, compress = "xz", version = 3L)
stopifnot(identical(readRDS(output_file), fixture),
          identical(digest::digest(calibration_file, file = TRUE, algo = "sha256"),
                    calibration_sha256))
cat("Saved synthetic 1200-by-48 Caudate and Liver matrices.\n")
cat("SHA-256:", digest::digest(output_file, file = TRUE, algo = "sha256"), "\n")
