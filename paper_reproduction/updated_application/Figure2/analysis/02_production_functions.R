# Helpers for one recorded differential experiment at a time.

# Verify the actual review receipt and every input it covers.
verify_reviewed_production <- function(job, job_file, cases_file, pilot_file,
                                       helper_file, runner_file) {
  if (job$stage != "production") return(invisible(TRUE))
  hash <- function(path) digest::digest(path, file = TRUE, algo = "sha256")
  resolve <- function(path) if (grepl("^/", path)) path else file.path(dirname(job_file), path)
  stopifnot(is.character(job$independent_review_file),
    is.character(job$independent_review_sha256))
  review_file <- resolve(job$independent_review_file)
  stopifnot(file.exists(review_file), identical(hash(review_file), job$independent_review_sha256))
  review <- jsonlite::read_json(review_file, simplifyVector = TRUE)
  archive <- resolve(job$source_archive_file)
  stopifnot(file.exists(archive), identical(hash(archive), job$source_archive_sha256))
  actual <- list(model_id = job$model_id, expression_backend = job$expression_backend,
    source_pilot_sha256 = hash(pilot_file), cases_sha256 = hash(cases_file),
    production_functions_sha256 = hash(helper_file), runner_script_sha256 = hash(runner_file),
    source_archive_sha256 = hash(archive))
  stopifnot(identical(review$status, "pass"),
    all(vapply(names(actual), function(name) identical(review[[name]], actual[[name]]), logical(1))),
    identical(review$RNGkind, job$RNGkind),
    identical(review$installed_function_body_sha256, job$installed_function_body_sha256))
  invisible(TRUE)
}

# 1. Map a recorded replicate ID without changing the installed model ------
recorded_replicate_analysis <- function(replicate_id) {
  stopifnot(length(replicate_id) == 1L, is.finite(replicate_id),
    replicate_id > 0, replicate_id == as.integer(replicate_id))
  namespace <- asNamespace("SCP")
  scope <- new.env(parent = namespace)
  scope$original_simulator <- get("simCircadianDiff", envir = namespace)
  scope$recorded_id <- as.integer(replicate_id)
  scope$last_truth <- scope$last_draws <- NULL
  scope$simCircadianDiff <- function(...) {
    args <- list(...)
    stopifnot(identical(as.numeric(args$sim.seed), 1000))
    args$sim.seed <- recorded_id * 1000L
    stopifnot(is.finite(args$sim.seed + args$n1 * 7919L + args$n2 * 104729L))
    sim <- do.call(original_simulator, args)
    last_truth <<- sim$ground_truth
    last_draws <<- sim$parameter_draws
    sim
  }
  environment(scope$simCircadianDiff) <- scope
  for (name in c("runSimsDiff", "runDifferentialPower")) {
    fun <- get(name, envir = namespace)
    environment(fun) <- scope
    assign(name, fun, envir = scope)
  }
  scope
}

# 2. Independently check the realized endpoint nulls and both source pools --
check_realized_truth <- function(gt, draws, bio) {
  DR <- gt$amplitude1^2 / gt$sigma^2 != gt$amplitude2^2 / gt$sigma2^2
  both <- gt$amplitude1 > 0 & gt$amplitude2 > 0
  delta_phase <- ((gt$phase2 - gt$phase1 + bio$period / 2) %% bio$period) - bio$period / 2
  DP <- both & delta_phase != 0
  DM <- gt$mesor1 != gt$mesor2
  stopifnot(identical(DR, gt$diff_type %in% c(2L, 3L)),
    identical(DP, gt$diff_type == 4L), identical(DM, gt$diff_type == 5L),
    all(gt$mesor1[!DM] == gt$mesor2[!DM]),
    all(gt$sigma[both] == gt$sigma2[both]))
  first <- which(gt$diff_type == 2L)
  second <- which(gt$diff_type == 3L)
  i1 <- draws$pair_index1[first]
  i2 <- draws$pair_index2[second]
  stopifnot(!anyNA(c(i1, i2)),
    identical(draws$raw_amplitude1[first], unname(bio$amplitude[i1])),
    identical(draws$raw_sigma1[first], unname(bio$sigma_rhythmic[i1])),
    identical(draws$raw_amplitude2[second], unname(bio$amplitude2[i2])),
    identical(draws$raw_sigma2[second], unname(bio$sigma_rhythmic2[i2])),
    all(gt$amplitude1[first] == pmax(bio$amplitude[i1], .05)),
    all(gt$sigma[first] == pmax(bio$sigma_rhythmic[i1], 1e-6)),
    all(gt$amplitude2[second] == pmax(bio$amplitude2[i2], .05)),
    all(gt$sigma2[second] == pmax(bio$sigma_rhythmic2[i2], 1e-6)))
  data.frame(DR = sum(DR), DP = sum(DP), DM = sum(DM),
    DP_at_least2h = sum(DP & abs(delta_phase) >= 2),
    shared_mean_null = TRUE, shared_SNR_null = TRUE,
    both_source_index_pairs = TRUE)
}

# 3. Summarize the full BH universe and the separate larger-shift power -----
compact_differential_result <- function(result, case, valid, error = "", gt = NULL) {
  if (is.null(gt) && !is.null(result)) gt <- result$ground_truth[[1]]
  rows <- list()
  for (endpoint in c("DR", "DP", "DM")) {
    target <- if (!is.null(gt)) gt[[paste0("is_", endpoint)]] else NULL
    if (valid) {
      pv <- result[[paste0("pval_", endpoint)]][, 1, 1]
      qv <- result[[paste0("fdr_", endpoint)]][, 1, 1]
      stopifnot(all(is.finite(pv)), all(is.finite(qv)), length(pv) == 5000L,
        max(abs(qv - p.adjust(pv, "BH"))) < 1e-12)
    }
    for (threshold in c(.01, .05, .10, .20)) {
      discovered <- if (valid) qv <= threshold else NULL
      for (population in if (endpoint == "DP") c("all_assigned", "phase_at_least2h") else "all_assigned") {
        selected <- target
        if (!is.null(gt) && population == "phase_at_least2h")
          selected <- target & abs(gt$phase_diff) >= 2
        TD <- if (valid) sum(discovered & selected) else NA_integer_
        FD <- if (valid) sum(discovered & !target) else NA_integer_
        discoveries <- if (valid) sum(discovered) else NA_integer_
        count <- if (!is.null(gt)) sum(selected) else NA_integer_
        rows[[paste(endpoint, threshold, population)]] <- data.frame(
          case_id = case$case_id, design = case$design, B = case$B,
          N_per_group = case$N_per_group, replicate_id = case$replicate_id,
          endpoint = endpoint, target_population = population, BH_threshold = threshold,
          DCP_gate = .20, valid = valid, true_count = count,
          TD = TD, FD_against_actual_endpoint_null = FD, discoveries = discoveries,
          power = if (valid && count > 0) TD / count else NA_real_,
          FDP = if (!valid) NA_real_ else if (discoveries) FD / discoveries else 0,
          error = error, stringsAsFactors = FALSE)
      }
    }
  }
  do.call(rbind, rows)
}

# 4. Run or verify one immutable case --------------------------------------
run_recorded_case <- function(case, pilot, job, output) {
  hash <- function(path) digest::digest(path, file = TRUE, algo = "sha256")
  stopifnot(identical(job$expression_backend, "compiled"),
    isTRUE(get(".CPP_LOADED", envir = asNamespace("SCP"), inherits = TRUE)),
    all(c("_SCP_sim_cosinor_expr_cpp", "_SCP_cosinor_pvals_cpp") %in%
      names(getDLLRegisteredRoutines("SCP")$.Call)))
  path <- file.path(output, "cases", case$case_id)
  receipt_file <- file.path(path, "receipt.json")
  case_hash <- digest::digest(case[c("case_id", "design", "B", "N_per_group", "replicate_id")], algo = "sha256")
  if (dir.exists(path)) {
    stopifnot(file.exists(receipt_file))
    receipt <- jsonlite::read_json(receipt_file, simplifyVector = TRUE)
    stopifnot(identical(receipt$case_id, case$case_id),
      identical(receipt$case_sha256, case_hash),
      receipt$status %in% c("complete", "failed"),
      identical(receipt$source_pilot_sha256, job$source_pilot_sha256),
      identical(receipt$model_id, job$model_id),
      identical(receipt$compact_sha256, hash(file.path(path, "per_replicate.csv"))))
    if (is.character(receipt$raw_sha256) && !is.na(receipt$raw_sha256))
      stopifnot(identical(receipt$raw_sha256, hash(file.path(path, "raw_result.rds"))))
    return(data.frame(case_id = case$case_id, status = receipt$status,
      resumed = TRUE, elapsed_seconds = receipt$elapsed_seconds))
  }
  dir.create(path, recursive = TRUE)
  bio <- pilot
  bio$ngenes <- 5000L
  if (case$design == "active") {
    stopifnot(is.finite(case$B), case$B >= 2, case$N_per_group %% case$B == 0)
    schedule <- rep(seq(0, bio$period, length.out = case$B + 1L)[seq_len(case$B)],
      each = case$N_per_group / case$B)
    bio$cts <- bio$cts2 <- schedule
  }
  stopifnot(length(bio$amplitude2) == length(bio$sigma_rhythmic2),
    isTRUE(bio$paired_sigma2), !is.null(bio$sigma_rhythmic2))
  analysis_scope <- recorded_replicate_analysis(case$replicate_id)
  design <- SCP::CircadianDesignOptions(sample_sizes = case$N_per_group, nsims = 1,
    design = case$design, cts = bio$cts, test_types = c("DR", "DP", "DM"))
  analysis <- SCP::CircadianAnalysisOptions(alpha = .20, p.adjust.method = "BH",
    fdr_thresholds = c(.01, .05, .10, .20), parallel.ncores = 1L)
  jsonlite::write_json(list(case = case, model_id = job$model_id,
    source_pilot_sha256 = job$source_pilot_sha256,
    mapped_seed = case$replicate_id * 1000L,
    expression_seed = case$replicate_id * 1000L + case$N_per_group * (7919L + 104729L),
    started_UTC = format(Sys.time(), tz = "UTC", usetz = TRUE)),
    file.path(path, "input.json"), pretty = TRUE, auto_unbox = TRUE, digits = 17)
  started <- proc.time()[["elapsed"]]
  warnings <- character()
  error <- ""
  result <- tryCatch(withCallingHandlers(
    analysis_scope$runDifferentialPower(bio, design, analysis, methods = "DCP",
      test_types = c("DR", "DP", "DM"), mc.cores = 1, plot = FALSE, verbose = FALSE),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) {
      error <<- conditionMessage(e)
      NULL
    })
  elapsed <- proc.time()[["elapsed"]] - started
  valid <- !is.null(result)
  truth_checks <- NULL
  if (!is.null(result)) {
    truth_checks <- tryCatch(check_realized_truth(result$ground_truth[[1]],
      result$parameter_draws[[1]], bio), error = function(e) {
        error <<- conditionMessage(e)
        NULL
      })
    valid <- !is.null(truth_checks)
  }
  if (valid && length(warnings)) {
    valid <- FALSE
    error <- paste("Warnings require review:", paste(warnings, collapse = " | "))
  }
  compact <- tryCatch(compact_differential_result(result, case, valid, error,
    gt = if (!is.null(result)) result$ground_truth[[1]] else analysis_scope$last_truth),
    error = function(e) {
      error <<- conditionMessage(e)
      valid <<- FALSE
      compact_differential_result(NULL, case, FALSE, error, gt = analysis_scope$last_truth)
    })
  utils::write.csv(compact, file.path(path, "per_replicate.csv"), row.names = FALSE)
  writeLines(warnings, file.path(path, "warnings.txt"))
  if (!is.null(result)) {
    saveRDS(result, file.path(path, "raw_result.rds"), compress = "gzip")
    if (!is.null(truth_checks))
      utils::write.csv(truth_checks, file.path(path, "truth_checks.csv"), row.names = FALSE)
  }
  if (!valid) {
    saveRDS(list(error = error, ground_truth = analysis_scope$last_truth,
      parameter_draws = analysis_scope$last_draws), file.path(path, "failure.rds"))
  }
  truth <- if (!is.null(result)) result$ground_truth[[1]] else analysis_scope$last_truth
  receipt <- list(status = if (valid) "complete" else "failed", case_id = case$case_id,
    case_sha256 = case_hash,
    model_id = job$model_id, expression_backend = job$expression_backend,
    DLL_sha256 = job$runtime_DLL_sha256, source_pilot_sha256 = job$source_pilot_sha256,
    elapsed_seconds = elapsed, warnings = length(warnings), error = error,
    ground_truth_sha256 = if (!is.null(truth)) digest::digest(truth, algo = "sha256") else NA_character_,
    compact_sha256 = hash(file.path(path, "per_replicate.csv")),
    raw_sha256 = if (!is.null(result)) hash(file.path(path, "raw_result.rds")) else NA_character_,
    completed_UTC = format(Sys.time(), tz = "UTC", usetz = TRUE))
  jsonlite::write_json(receipt, receipt_file, pretty = TRUE, auto_unbox = TRUE, digits = 17)
  data.frame(case_id = case$case_id, status = receipt$status,
    resumed = FALSE, elapsed_seconds = elapsed)
}
