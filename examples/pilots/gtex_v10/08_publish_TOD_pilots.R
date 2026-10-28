# Prepare the TOD-only catalogue from the validated paired analysis.

# 1. Read the staged fits and choose a new publication directory -----------
stopifnot(requireNamespace("digest", quietly = TRUE),
          requireNamespace("jsonlite", quietly = TRUE))
stage <- Sys.getenv("SCP_GTEX_V10_STAGE", "output/gtex_v10_chiral_protocol")
publication <- Sys.getenv("SCP_GTEX_TOD_PUBLICATION",
                          file.path(stage, "publication_TOD"))
source_pilots <- file.path(stage, "inst", "extdata", "pilots")
native_pilots <- Sys.getenv("SCP_PILOT_SOURCE_ROOT", "inst/extdata/pilots")
manifest_file <- file.path(source_pilots, "manifest.csv")
summary_file <- file.path(stage, "summary", "gtex_v10_signal_summary.csv")
stopifnot(file.exists(manifest_file), file.exists(summary_file),
          dir.exists(native_pilots))
if (dir.exists(publication) &&
    length(list.files(publication, all.files = TRUE, no.. = TRUE)))
  stop("The publication directory contains files. Preserve it or choose a new SCP_GTEX_TOD_PUBLICATION.")
manifest <- utils::read.csv(manifest_file, stringsAsFactors = FALSE)
signal <- utils::read.csv(summary_file, stringsAsFactors = FALSE)

# 2. Retain TOD for GTEx and preserve the other study rows ------------------
gtex <- manifest$dataset == "GTEx"
keep <- !gtex | (manifest$gtex_release == "v10" & manifest$time_basis == "TOD")
manifest <- manifest[which(keep), , drop = FALSE]
gtex <- manifest$dataset == "GTEx"
stopifnot(sum(gtex) == 46L,
          !anyDuplicated(manifest[gtex, c("tissue", "condition")]))
manifest$is_default[gtex] <- TRUE
manifest$note[gtex & manifest$tissue == "Whole_Blood"] <- ""
manifest$citation[gtex] <- paste("GTEx Consortium (2020).",
  "The GTEx Consortium atlas of genetic regulatory effects across human tissues.",
  "Science 369(6509):1318-1330. doi:10.1126/science.aaz1776.")
signal <- signal[signal$gtex_release == "v10" & signal$time_basis == "TOD",
                 , drop = FALSE]
stopifnot(nrow(signal) == 46L, setequal(signal$tissue, manifest$tissue[gtex]))
signal$is_default <- TRUE

# 3. Count the observed TOD samples in one-hour bins ---------------------
sampling_rows <- list()
for (i in which(gtex)) {
  row <- manifest[i, , drop = FALSE]
  cohort_file <- file.path(stage, "private", "paired_cohorts", paste0(row$tissue, ".rds"))
  stopifnot(file.exists(cohort_file))
  times <- readRDS(cohort_file)$TOD
  stopifnot(is.numeric(times), length(times) == row$n, all(is.finite(times)))
  bins <- hist(times %% 24, breaks = 0:24, right = FALSE, include.lowest = TRUE, plot = FALSE)
  stopifnot(length(bins$counts) == 24L, sum(bins$counts) == row$n)
  sampling_rows[[length(sampling_rows) + 1L]] <- data.frame(
    dataset = row$dataset, tissue = row$tissue, condition = row$condition,
    gtex_release = row$gtex_release, time_basis = "TOD",
    bin_start = 0:23, bin_end = 1:24, count = as.integer(bins$counts),
    n_samples = as.integer(row$n), stringsAsFactors = FALSE)
}
sampling_histograms <- do.call(rbind, sampling_rows)
stopifnot(nrow(sampling_histograms) == 46L * 24L)

# 4. Copy the original TOD fits and verify every file ----------------------
destination_pilots <- file.path(publication, "inst", "extdata", "pilots")
dir.create(destination_pilots, recursive = TRUE, showWarnings = FALSE)
copied <- list()
for (relative in manifest$file[gtex]) {
  for (K in 1:2) {
    file <- if (K == 1L) relative else sub("[.]rds$", "_2H.rds", relative)
    source_file <- file.path(source_pilots, file)
    stopifnot(file.exists(source_file))
    pilot <- readRDS(source_file)
    stopifnot(identical(pilot$gtex_release, "v10"),
              identical(pilot$time_basis, "TOD"),
              !length(pilot$times), !length(pilot$cts), !length(pilot$raw))
    destination <- file.path(destination_pilots, file)
    dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
    stopifnot(file.copy(source_file, destination, overwrite = FALSE))
    hash <- digest::digest(file = source_file, algo = "sha256")
    stopifnot(identical(hash, digest::digest(file = destination, algo = "sha256")))
    copied[[length(copied) + 1L]] <- data.frame(file = file, dataset = "GTEx",
      K = K, sha256 = hash, stringsAsFactors = FALSE)
  }
}
stopifnot(length(copied) == 92L)

# 5. Preserve all existing non-GTEx pilot files ----------------------------
native_files <- list.files(native_pilots, pattern = "[.]rds$", recursive = TRUE)
native_files <- native_files[!grepl("^GTEx_", basename(native_files))]
for (file in native_files) {
  source_file <- file.path(native_pilots, file)
  destination <- file.path(destination_pilots, file)
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  stopifnot(file.copy(source_file, destination, overwrite = FALSE))
  hash <- digest::digest(file = source_file, algo = "sha256")
  stopifnot(identical(hash, digest::digest(file = destination, algo = "sha256")))
  copied[[length(copied) + 1L]] <- data.frame(file = file, dataset = "native",
    K = if (grepl("_2H[.]rds$", file)) 2L else 1L,
    sha256 = hash, stringsAsFactors = FALSE)
}

# 6. Save the active catalogue, summary, and publication receipt -----------
dir.create(file.path(publication, "summary"), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(manifest, file.path(destination_pilots, "manifest.csv"), row.names = FALSE)
utils::write.csv(signal, file.path(publication, "summary", "gtex_v10_signal_summary.csv"),
                 row.names = FALSE)
sampling_file <- file.path(destination_pilots, "gtex_v10_TOD_sampling_histograms.csv")
utils::write.csv(sampling_histograms, sampling_file, row.names = FALSE)
files <- do.call(rbind, copied)
utils::write.csv(files, file.path(publication, "pilot_hashes.csv"), row.names = FALSE)
jsonlite::write_json(list(package_version = read.dcf("DESCRIPTION", fields = "Version")[[1]],
  input_time = "Season-adjusted TOD", gtex_tissues = sum(gtex), gtex_pilot_files = 92L,
  native_pilot_files = length(native_files), source_manifest = manifest_file,
  source_manifest_sha256 = digest::digest(file = manifest_file, algo = "sha256"),
  source_summary_sha256 = digest::digest(file = summary_file, algo = "sha256"),
  pilot_values_changed = FALSE, DIP_in_active_catalogue = FALSE,
  sampling_histogram_file = "inst/extdata/pilots/gtex_v10_TOD_sampling_histograms.csv",
  sampling_histogram_sha256 = digest::digest(file = sampling_file, algo = "sha256"),
  sampling_histogram_bin_hours = 1L, sampling_histogram_bins_per_tissue = 24L,
  sampling_histogram_rows = nrow(sampling_histograms),
  sampling_histogram_counts = "Observed TOD counts with left-closed, right-open one-hour bins",
  cohort = "Existing matched v10 cohort used by the paired research analysis"),
  file.path(publication, "publication_recipe.json"), auto_unbox = TRUE, pretty = TRUE)
cat("TOD publication prepared:", normalizePath(publication), "\n")
