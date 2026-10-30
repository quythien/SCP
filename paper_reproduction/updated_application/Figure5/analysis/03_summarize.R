# Collect completed A/B jobs: Rscript 03_summarize.R JOB_DIRECTORY NEW_REPORT_DIRECTORY

# 1. Verify all job hashes, frozen calibration and complete replicate counts -
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
bundle <- dirname(normalizePath(script, mustWork = TRUE))
source(file.path(bundle, "00_exact_active_engine.R"))
gate <- active_verify_gate(bundle)
settings <- readRDS(file.path(bundle, "settings.rds"))
jobs <- read.csv(file.path(bundle, "job_manifest.csv"), stringsAsFactors = FALSE)
input <- normalizePath(args[1], mustWork = TRUE)
report <- normalizePath(args[2], mustWork = FALSE)
stopifnot(!dir.exists(report))
paths <- file.path(input, sprintf("job_%04d.rds", jobs$job_id))
stopifnot(all(file.exists(paths)), all(file.exists(paste0(paths, ".sha256"))))
values <- vector("list", length(paths))
calibration_hash <- active_sha256(file.path(bundle, "calibrations.rds"))
for (i in seq_along(paths)) {
  stopifnot(active_sha256(paths[i]) == readLines(paste0(paths[i], ".sha256")))
  values[[i]] <- readRDS(paths[i])
  stopifnot(identical(values[[i]]$job, jobs[i, , drop = FALSE]),
            values[[i]]$gate_sha256 == gate$gate_sha256,
            values[[i]]$calibration_sha256 == calibration_hash)
}
metrics <- do.call(rbind, lapply(values, `[[`, "metrics"))
stopifnot(nrow(metrics) == 600000L,
          !anyDuplicated(metrics[c("panel", "N", "B", "FDR_target", "replicate")]))
groups <- split(metrics, interaction(metrics$panel, metrics$N, metrics$B,
                                    metrics$FDR_target, drop = TRUE))
stopifnot(length(groups) == 1000L, all(vapply(groups, nrow, integer(1)) == 600L))
dir.create(report, recursive = TRUE)
saveRDS(metrics, file.path(report, "all_replicate_metrics.rds"), compress = "xz")

# 2. Retain power and FDP means, MCSEs and exact within-replicate B overlap -
mcse <- function(x) sd(x) / sqrt(length(x))
power <- do.call(rbind, lapply(groups, function(x) data.frame(panel = x$panel[1],
  K = x$K[1], N = x$N[1], B = x$B[1], m = x$m[1], FDR_target = x$FDR_target[1],
  nsims = nrow(x), power = mean(x$power), power_MCSE = mcse(x$power),
  mean_FDP = mean(x$FDP), FDP_MCSE = mcse(x$FDP), mean_TD = mean(x$TD), mean_FD = mean(x$FD))))
power <- power[order(power$panel, power$FDR_target, power$B, power$N), ]
write.csv(power, file.path(report, "power_and_FDP_summary.csv"), row.names = FALSE)
coupled <- split(metrics, interaction(metrics$panel, metrics$N, metrics$FDR_target,
                                    metrics$replicate, drop = TRUE))
overlap <- vapply(coupled, function(x) max(x$power) - min(x$power), numeric(1))
FDP_overlap <- vapply(coupled, function(x) max(x$FDP) - min(x$FDP), numeric(1))
stopifnot(all(overlap == 0), all(FDP_overlap == 0))
write.csv(data.frame(coupled_panel_replicates = length(overlap),
  maximum_power_difference_across_B = max(overlap),
  maximum_FDP_difference_across_B = max(FDP_overlap),
  reason = "Identical balanced harmonic information and common canonical Gaussian draws"),
  file.path(report, "B_overlap_check.csv"), row.names = FALSE)

# 3. Report tested N80 brackets and descriptive interpolation for every B --
crossings <- list()
for (panel in settings$panels) for (B in active_B(panel)) for (alpha in settings$FDR) {
  x <- power[power$panel == panel & power$B == B & power$FDR_target == alpha, ]
  x <- x[order(x$N), ]
  hits <- which(x$power >= .8)
  first <- if (length(hits)) hits[1] else NA_integer_
  bracketed <- !is.na(first) && first > 1L
  crossings[[paste(panel, B, alpha)]] <- data.frame(panel = panel, B = B,
    FDR_target = alpha, N_low = if (bracketed) x$N[first - 1L] else NA_integer_,
    N_high = if (!is.na(first)) x$N[first] else NA_integer_,
    power_low = if (bracketed) x$power[first - 1L] else NA_real_,
    power_high = if (!is.na(first)) x$power[first] else NA_real_,
    MCSE_low = if (bracketed) x$power_MCSE[first - 1L] else NA_real_,
    MCSE_high = if (!is.na(first)) x$power_MCSE[first] else NA_real_,
    N80_interpolated = if (bracketed) x$N[first - 1L] + (.8 - x$power[first - 1L]) *
      (x$N[first] - x$N[first - 1L]) / (x$power[first] - x$power[first - 1L]) else NA_real_,
    extension_needed = is.na(first),
    status = if (is.na(first)) "Above maximum tested N" else if (first == 1L)
      "At or below minimum tested N" else "Bracketed on tested grid")
}
crossings <- do.call(rbind, crossings)
write.csv(crossings, file.path(report, "N80_brackets.csv"), row.names = FALSE)

# 4. Supply numeric panel matrices for the existing five-panel layout ------
matrices <- list()
for (panel in settings$panels) {
  B_values <- active_B(panel)
  M <- SE <- matrix(NA_real_, length(settings$N_grid), length(B_values),
                   dimnames = list(settings$N_grid, paste0("B=", B_values)))
  for (j in seq_along(B_values)) {
    x <- power[power$panel == panel & power$B == B_values[j] & power$FDR_target == .05, ]
    x <- x[match(settings$N_grid, x$N), ]
    M[, j] <- x$power
    SE[, j] <- x$power_MCSE
  }
  matrices[[panel]] <- list(power = M, MCSE = SE)
}
cache <- list(N = settings$N_grid, B_K1 = active_B("A"), B_K2 = active_B("B"),
  pwr_K1 = matrices$A$power, se_K1 = matrices$A$MCSE,
  pwr_K2 = matrices$B$power, se_K2 = matrices$B$MCSE,
  NSIMS = settings$nsims, NGENES = settings$G, FDR = .05, PERIOD = 24,
  calibrations_differ_between_panels = TRUE, balanced_B_curves_coincide = TRUE,
  nonidentifiable_designs_not_included = TRUE, gate_sha256 = gate$gate_sha256)
saveRDS(cache, file.path(report, "figure_2_5_AB_cache.rds"), compress = "xz")
saveRDS(list(status = "complete", gate_sha256 = gate$gate_sha256,
  jobs = nrow(jobs), panel_replicates = 30000L, B_cell_replicates = 150000L,
  metric_rows = nrow(metrics), maximum_B_power_difference = max(overlap),
  extensions_required = sum(crossings$extension_needed)),
  file.path(report, "completion_receipt.rds"), compress = "xz")
cat("Saved active A/B summaries and five-panel-layout cache to", report, "\n")
