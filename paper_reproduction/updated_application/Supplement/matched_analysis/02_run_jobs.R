# Portable base-R worker: Rscript 02_run_jobs.R JOB_ID OUTPUT_DIRECTORY
# One process runs one 60-replicate job. Root coordinates process placement.

# 1. Locate the portable bundle and require its independent hash gate ------
arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 2L)
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
bundle <- dirname(normalizePath(script, mustWork = TRUE))
source(file.path(bundle, "00_paired_engine.R"))
gate <- matched_verify_gate(bundle)
settings <- readRDS(file.path(bundle, "settings.rds"))
input <- readRDS(file.path(bundle, "calibration_and_sampler.rds"))
jobs <- read.csv(file.path(bundle, "job_manifest.csv"))
job_id <- as.integer(arguments[1])
stopifnot(length(job_id) == 1L, is.finite(job_id), job_id %in% jobs$job_id)
job <- jobs[jobs$job_id == job_id, , drop = FALSE]
stopifnot(nrow(job) == 1L, settings$G == 5000L, settings$nsims == 600L,
          job$N %% 6L == 0L, job$replicate_start >= 1L,
          job$replicate_end <= settings$nsims)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
output <- normalizePath(arguments[2], mustWork = FALSE)
dir.create(output, recursive = TRUE, showWarnings = FALSE)
destination <- file.path(output, sprintf("job_%04d.rds", job_id))
checksum <- paste0(destination, ".sha256")
replicates <- seq.int(job$replicate_start, job$replicate_end)

# 2. Resume completed jobs only when their saved receipt and hash agree ----
if (file.exists(destination)) {
  stopifnot(file.exists(checksum), matched_sha256(destination) == readLines(checksum),
            identical(readRDS(destination)$job, job),
            identical(readRDS(destination)$gate_sha256, gate$gate_sha256))
  cat(sprintf("Job %d already complete and verified.\n", job_id))
  quit(save = "no", status = 0L)
}
stopifnot(!file.exists(checksum))

# 3. Simulate the paired designs and both detectors without nested workers -
start <- proc.time()[["elapsed"]]
rows <- vector("list", length(replicates))
for (j in seq_along(replicates)) {
  i <- replicates[j]
  seed <- as.integer(settings$seed_base + job$N * 100000L + i)
  value <- simulate_matched_pair(input$bio, input$sampler, job$N, seed, settings$G)
  value$metrics$replicate <- i
  value$metrics$seed <- seed
  rows[[j]] <- value$metrics
}
metrics <- do.call(rbind, rows)
key <- c("design", "N", "FDR_target", "replicate")
stopifnot(nrow(metrics) == length(replicates) * 8L, !anyDuplicated(metrics[key]),
          all(metrics$n_targets == settings$expected_targets_per_experiment),
          all(is.finite(metrics$power_K1)), all(is.finite(metrics$power_K2)),
          all(metrics$FDP_K1 >= 0 & metrics$FDP_K1 <= 1),
          all(metrics$FDP_K2 >= 0 & metrics$FDP_K2 <= 1))
result <- list(job = job, gate_sha256 = gate$gate_sha256,
  calibration_sha256 = matched_sha256(file.path(bundle, "calibration_and_sampler.rds")),
  metrics = metrics, elapsed_seconds = proc.time()[["elapsed"]] - start,
  completed_at_UTC = format(Sys.time(), tz = "UTC", usetz = TRUE),
  host = Sys.info()[["nodename"]], R_version = R.version.string,
  RNG_kind = RNGkind(), within_process_workers = 1L)

# 4. Publish a complete new job artifact and its integrity checksum --------
temporary <- tempfile(pattern = sprintf("job_%04d_partial_", job_id), tmpdir = output)
saveRDS(result, temporary, compress = "xz")
stopifnot(identical(readRDS(temporary)$metrics, metrics), !file.exists(destination),
          file.rename(temporary, destination))
writeLines(matched_sha256(destination), checksum)
cat(sprintf("Job %d: N=%d, %d parameter pairs, %.3f seconds.\n",
            job_id, job$N, length(replicates), result$elapsed_seconds))
