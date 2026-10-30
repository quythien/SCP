# Portable base-R worker: Rscript 02_run_job.R JOB_ID NEW_OUTPUT_DIRECTORY

# 1. Require the frozen, independently approved bundle ---------------------
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
bundle <- dirname(normalizePath(script, mustWork = TRUE))
source(file.path(bundle, "00_exact_active_engine.R"))
gate <- active_verify_gate(bundle)
settings <- readRDS(file.path(bundle, "settings.rds"))
calibrations <- readRDS(file.path(bundle, "calibrations.rds"))
jobs <- read.csv(file.path(bundle, "job_manifest.csv"), stringsAsFactors = FALSE)
job <- jobs[jobs$job_id == as.integer(args[1]), , drop = FALSE]
stopifnot(nrow(job) == 1L, settings$G == 5000L, settings$nsims == 600L,
          job$N %% 24L == 0L)
panel <- job$panel
K <- if (panel == "A") 1L else 2L
for (B in active_B(panel)) active_geometry(job$N, B, K)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
output <- normalizePath(args[2], mustWork = FALSE)
dir.create(output, recursive = TRUE, showWarnings = FALSE)
path <- file.path(output, sprintf("job_%04d.rds", job$job_id))
checksum <- paste0(path, ".sha256")
if (file.exists(path)) {
  old <- readRDS(path)
  stopifnot(file.exists(checksum), active_sha256(path) == readLines(checksum),
            identical(old$job, job), identical(old$gate_sha256, gate$gate_sha256))
  cat("Verified completed job", job$job_id, "\n")
  quit(save = "no", status = 0L)
}
stopifnot(!file.exists(checksum))

# 2. Run sixty independent panel replicates, coupled across B --------------
replicates <- seq.int(job$replicate_start, job$replicate_end)
rows <- vector("list", length(replicates))
start <- proc.time()[["elapsed"]]
for (j in seq_along(replicates)) {
  i <- replicates[j]
  seed <- as.integer(settings$seed_base[[panel]] + job$N * 100000L + i)
  value <- simulate_active_panel(calibrations[[panel]], panel, job$N, seed, settings$G)
  value$metrics$replicate <- i
  value$metrics$seed <- seed
  rows[[j]] <- value$metrics
}
metrics <- do.call(rbind, rows)
stopifnot(nrow(metrics) == length(replicates) * length(active_B(panel)) * 4L,
          !anyDuplicated(metrics[c("panel", "N", "B", "FDR_target", "replicate")]),
          all(metrics$n_targets == round(settings$G * calibrations[[panel]]$prop_rhythmic)))
result <- list(job = job, metrics = metrics, gate_sha256 = gate$gate_sha256,
  calibration_sha256 = active_sha256(file.path(bundle, "calibrations.rds")),
  elapsed_seconds = proc.time()[["elapsed"]] - start,
  host = Sys.info()[["nodename"]], R_version = R.version.string,
  completed_UTC = format(Sys.time(), tz = "UTC", usetz = TRUE))

# 3. Publish new immutable job data and its checksum -----------------------
temporary <- tempfile(pattern = "active_partial_", tmpdir = output)
saveRDS(result, temporary, compress = "xz")
stopifnot(identical(readRDS(temporary)$metrics, metrics), !file.exists(path),
          file.rename(temporary, path))
writeLines(active_sha256(path), checksum)
cat(sprintf("Job %d: panel %s, N=%d, %d replicates, %.3f seconds.\n",
            job$job_id, panel, job$N, length(replicates), result$elapsed_seconds))
