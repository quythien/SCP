# Portable base-R collection: Rscript 03_summarize.R JOB_DIRECTORY NEW_REPORT_DIRECTORY

# 1. Verify every planned job before calculating final estimates -----------
arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 2L)
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
bundle <- dirname(normalizePath(script, mustWork = TRUE))
source(file.path(bundle, "00_paired_engine.R"))
gate <- matched_verify_gate(bundle)
settings <- readRDS(file.path(bundle, "settings.rds"))
jobs <- read.csv(file.path(bundle, "job_manifest.csv"))
input_dir <- normalizePath(arguments[1], mustWork = TRUE)
report <- normalizePath(arguments[2], mustWork = FALSE)
stopifnot(!dir.exists(report))
paths <- file.path(input_dir, sprintf("job_%04d.rds", jobs$job_id))
stopifnot(all(file.exists(paths)), all(file.exists(paste0(paths, ".sha256"))))
values <- vector("list", length(paths))
for (i in seq_along(paths)) {
  stopifnot(matched_sha256(paths[i]) == readLines(paste0(paths[i], ".sha256")))
  values[[i]] <- readRDS(paths[i])
  stopifnot(identical(values[[i]]$job, jobs[i, , drop = FALSE]),
            identical(values[[i]]$gate_sha256, gate$gate_sha256),
            values[[i]]$calibration_sha256 ==
              matched_sha256(file.path(bundle, "calibration_and_sampler.rds")))
}
metrics <- do.call(rbind, lapply(values, `[[`, "metrics"))
key <- c("design", "N", "FDR_target", "replicate")
stopifnot(nrow(metrics) == length(settings$N_grid) * settings$nsims * 8L,
          !anyDuplicated(metrics[key]),
          all(metrics$n_targets == settings$expected_targets_per_experiment))
groups <- split(metrics, interaction(metrics$design, metrics$N, metrics$FDR_target, drop = TRUE))
stopifnot(all(vapply(groups, nrow, integer(1)) == settings$nsims))
dir.create(report, recursive = TRUE)
saveRDS(metrics, file.path(report, "all_replicate_metrics.rds"), compress = "xz")
mcse <- function(x) stats::sd(x) / sqrt(length(x))

# 2. Summarize power, FDP, discoveries and matched detector differences ----
summary <- do.call(rbind, lapply(groups, function(x) data.frame(
  design = x$design[1], N = x$N[1], FDR_target = x$FDR_target[1], nsims = nrow(x),
  power_K1 = mean(x$power_K1), power_K2 = mean(x$power_K2),
  power_MCSE_K1 = mcse(x$power_K1), power_MCSE_K2 = mcse(x$power_K2),
  K2_minus_K1_power = mean(x$paired_power_difference),
  K2_minus_K1_power_MCSE = mcse(x$paired_power_difference),
  mean_FDP_K1 = mean(x$FDP_K1), mean_FDP_K2 = mean(x$FDP_K2),
  FDP_MCSE_K1 = mcse(x$FDP_K1), FDP_MCSE_K2 = mcse(x$FDP_K2),
  K2_minus_K1_FDP = mean(x$FDP_K2 - x$FDP_K1),
  K2_minus_K1_FDP_MCSE = mcse(x$FDP_K2 - x$FDP_K1),
  mean_TD_K1 = mean(x$TD_K1), mean_TD_K2 = mean(x$TD_K2),
  mean_FD_K1 = mean(x$FD_K1), mean_FD_K2 = mean(x$FD_K2))))
summary <- summary[order(summary$design, summary$FDR_target, summary$N), ]
write.csv(summary, file.path(report, "power_and_FDP_summary.csv"), row.names = FALSE)

# 3. Keep within-replicate design differences and their Monte Carlo error --
design_groups <- split(metrics, interaction(metrics$N, metrics$FDR_target, drop = TRUE))
design_comparison <- do.call(rbind, lapply(design_groups, function(x) {
  passive <- x[x$design == "passive", ]
  balanced <- x[x$design == "balanced6", ]
  passive <- passive[order(passive$replicate), ]
  balanced <- balanced[order(balanced$replicate), ]
  stopifnot(identical(passive$replicate, balanced$replicate),
            identical(passive$seed, balanced$seed),
            identical(passive$n_targets, balanced$n_targets))
  do.call(rbind, lapply(1:2, function(K) {
    power_difference <- balanced[[paste0("power_K", K)]] - passive[[paste0("power_K", K)]]
    FDP_difference <- balanced[[paste0("FDP_K", K)]] - passive[[paste0("FDP_K", K)]]
    data.frame(N = x$N[1], FDR_target = x$FDR_target[1], K = K,
      nsims = nrow(passive), balanced6_minus_passive_power = mean(power_difference),
      balanced6_minus_passive_power_MCSE = mcse(power_difference),
      balanced6_minus_passive_FDP = mean(FDP_difference),
      balanced6_minus_passive_FDP_MCSE = mcse(FDP_difference))
  }))
}))
design_comparison <- design_comparison[order(design_comparison$K,
  design_comparison$FDR_target, design_comparison$N), ]
write.csv(design_comparison, file.path(report, "paired_design_comparison.csv"), row.names = FALSE)

# 4. Report tested N80 brackets without imposing a monotone power curve ----
crossings <- list()
for (design in settings$designs) for (alpha in settings$FDR) for (K in 1:2) {
  x <- summary[summary$design == design & summary$FDR_target == alpha, ]
  x <- x[order(x$N), ]
  p <- x[[paste0("power_K", K)]]
  se <- x[[paste0("power_MCSE_K", K)]]
  above <- which(p >= .8)
  first <- if (length(above)) above[1] else NA_integer_
  bracketed <- !is.na(first) && first > 1L
  estimate <- if (bracketed) x$N[first - 1L] + (.8 - p[first - 1L]) *
    (x$N[first] - x$N[first - 1L]) / (p[first] - p[first - 1L]) else NA_real_
  confident_above <- which(p - 1.96 * se >= .8)
  ca <- if (length(confident_above)) confident_above[1] else NA_integer_
  cb <- which(p + 1.96 * se < .8 & (is.na(ca) | seq_along(p) < ca))
  crossings[[paste(design, alpha, K)]] <- data.frame(design = design, FDR_target = alpha,
    K = K, N_low = if (bracketed) x$N[first - 1L] else NA_integer_,
    N_high = if (!is.na(first)) x$N[first] else NA_integer_,
    power_low = if (bracketed) p[first - 1L] else NA_real_,
    power_high = if (!is.na(first)) p[first] else NA_real_,
    MCSE_low = if (bracketed) se[first - 1L] else NA_real_,
    MCSE_high = if (!is.na(first)) se[first] else NA_real_,
    N80_interpolated = estimate,
    last_N_pointwise_upper_below80 = if (length(cb)) max(x$N[cb]) else NA_integer_,
    first_N_pointwise_lower_at80 = if (!is.na(ca)) x$N[ca] else NA_integer_,
    decreasing_grid_intervals = sum(diff(p) < 0),
    extension_needed = is.na(first),
    crossing_status = if (is.na(first)) "Above maximum tested N" else if (first == 1L)
      "At or below minimum tested N" else "Bracketed on tested grid")
}
crossings <- do.call(rbind, crossings)
write.csv(crossings, file.path(report, "N80_brackets.csv"), row.names = FALSE)
write.csv(crossings[crossings$FDR_target == settings$primary_FDR, ],
          file.path(report, "N80_brackets_BH05.csv"), row.names = FALSE)

# 5. Save interpretation and completion evidence for independent review ---
writeLines(c(
  "Appendix A: matched K1/K2 comparison under the same two-harmonic truth.",
  "This calibration also compares balanced-six versus passive sampling for each detector.",
  "These K1 curves do not represent the pure single-harmonic truth used for Figure 2.5A.",
  "Each of 600 replicates per N shares parameters and truth across designs. Design noise is independent.",
  "Both detectors share the same dataset within each design. K1 retains omitted second-harmonic energy in its residual denominator.",
  "Balanced-six uses phases 0, 4, 8, 12, 16 and 20 hours, with N/6 samples per phase.",
  "Passive sampling retains the frozen source grid and weights, Gaussian jitter SD 0.1 hours and modulo-24 wrap.",
  "BH thresholds are 0.01, 0.05, 0.10 and 0.20. BH 0.05 is primary. FDP is zero when no genes are discovered.",
  "All MCSEs describe fixed-pilot Monte Carlo error and exclude uncertainty in the fitted pilot.",
  "N80 brackets are tested-grid crossings. Interpolation is descriptive. Pointwise 1.96-MCSE bounds are not simultaneous confidence bands or an N80 confidence interval.",
  "Rows with extension_needed=TRUE require a separately reviewed grid extension. No detector or design advantage is assumed."
), file.path(report, "Interpretation.txt"))
receipt <- list(completed_at_UTC = format(Sys.time(), tz = "UTC", usetz = TRUE),
  jobs = nrow(jobs), paired_parameter_experiments = length(settings$N_grid) * settings$nsims,
  datasets = length(settings$N_grid) * settings$nsims * 2L,
  detector_evaluations = length(settings$N_grid) * settings$nsims * 4L,
  gate_sha256 = gate$gate_sha256, source_bundle_hashes = gate$hashes,
  metric_rows = nrow(metrics), total_worker_seconds = sum(vapply(values, `[[`, numeric(1), "elapsed_seconds")),
  N80_extensions_required = crossings[crossings$extension_needed, ],
  existing_deliverables_replaced = FALSE)
saveRDS(receipt, file.path(report, "completion_receipt.rds"), compress = "xz")
cat("Complete matched summaries and N80 brackets saved to", report, "\n")
