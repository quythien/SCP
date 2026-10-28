# Stage a demo bundle using the canonical app and active default pilots.
# Run from the repository root after integrating the validated v10 pilots.

# 1. Choose the source database and a new output directory -----------------
stopifnot(requireNamespace("SCP", quietly = TRUE),
          requireNamespace("digest", quietly = TRUE),
          requireNamespace("jsonlite", quietly = TRUE))
stage <- Sys.getenv("SCP_GTEX_V10_STAGE", "output/gtex_v10_chiral_protocol")
demo <- Sys.getenv("SCP_GTEX_V10_DEMO_DIR", file.path(stage, "demo"))
pilot_root <- normalizePath(Sys.getenv("SCP_DEMO_PILOT_ROOT",
                                       "inst/extdata/pilots"), mustWork = TRUE)
if (dir.exists(demo) && length(list.files(demo, all.files = TRUE, no.. = TRUE)))
  stop("The demo destination contains files. Preserve it before rebuilding or choose SCP_GTEX_V10_DEMO_DIR.")
manifest_file <- file.path(pilot_root, "manifest.csv")
manifest <- utils::read.csv(manifest_file, stringsAsFactors = FALSE)
version <- read.dcf("DESCRIPTION", fields = "Version")[[1]]
if (as.character(utils::packageVersion("SCP")) != version)
  stop("Install the current source version of SCP before generating the demo.")

# 2. Keep DIP for GTEx and the native input time for other studies ----------
manifest <- manifest[manifest$status == "ingested", , drop = FALSE]
manifest <- getFromNamespace(".scp_filter_time_basis", "SCP")(
  manifest, time_basis = "DIP")
gtex <- manifest$dataset == "GTEx"
stopifnot(any(gtex), all(manifest$gtex_release[gtex] == "v10"),
          all(manifest$time_basis[gtex] == "DIP"),
          !anyDuplicated(manifest[c("species", "dataset", "tissue", "condition")]))
dir.create(file.path(demo, "pilots"), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(manifest, file.path(demo, "pilots", "manifest.csv"),
                 row.names = FALSE)

# 3. Copy each available harmonic variant and verify its hash --------------
copied <- list()
for (i in seq_len(nrow(manifest))) {
  for (K in 1:2) {
    file <- if (K == 1L) manifest$file[i] else
      sub("\\.rds$", "_2H.rds", manifest$file[i])
    source_file <- file.path(pilot_root, file)
    if (!file.exists(source_file)) {
      if (K == 1L) stop("Missing default pilot: ", file)
      next
    }
    if (gtex[i]) {
      pilot <- readRDS(source_file)
      getFromNamespace(".scp_validate_gtex_pilot", "SCP")(
        pilot, manifest[i, , drop = FALSE])
      if (length(pilot$times) || length(pilot$cts) || length(pilot$raw$times))
        stop("The GTEx demo pilot contains individual sampling times: ", file)
    }
    destination <- file.path(demo, "pilots", file)
    dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
    stopifnot(file.copy(source_file, destination, overwrite = FALSE))
    hash <- digest::digest(file = source_file, algo = "sha256")
    stopifnot(identical(hash, digest::digest(file = destination, algo = "sha256")))
    copied[[length(copied) + 1L]] <- data.frame(file = file, K = K,
      dataset = manifest$dataset[i], time_basis = manifest$time_basis[i],
      sha256 = hash, stringsAsFactors = FALSE)
  }
}
files <- do.call(rbind, copied)
utils::write.csv(files, file.path(demo, "pilot_hashes.csv"), row.names = FALSE)

# 4. Use one app implementation and record its exact source ----------------
wrapper <- "inst/shiny_demo/app.R"
canonical_app <- "inst/shiny/app.R"
stopifnot(file.exists(wrapper), file.exists(canonical_app),
          file.copy(wrapper, file.path(demo, "app.R"), overwrite = FALSE))
writeLines(version, file.path(demo, "SCP_VERSION"))
writeLines(unname(tools::md5sum(canonical_app)), file.path(demo, "SCP_APP_MD5"))
revision <- Sys.getenv("SCP_RELEASE_SHA", "uncommitted source checkout")
jsonlite::write_json(list(package = "SCP", package_version = version,
  source_revision = revision, generated_at_UTC = format(Sys.time(), tz = "UTC", usetz = TRUE),
  canonical_app_sha256 = digest::digest(file = canonical_app, algo = "sha256"),
  source_manifest_sha256 = digest::digest(file = manifest_file, algo = "sha256"),
  gtEx_time_basis = "DIP", manifest_rows = nrow(manifest), pilot_files = nrow(files)),
  file.path(demo, "build_info.json"), auto_unbox = TRUE, pretty = TRUE)

# 5. Check representative pilots without running a power simulation --------
smoke_rows <- unique(c(which(gtex)[1], which(manifest$dataset == "GTEx" &
                                              manifest$tissue == "Liver")[1],
                       which(!gtex)[1]))
smoke_rows <- smoke_rows[is.finite(smoke_rows)]
for (i in smoke_rows) {
  pilot <- readRDS(file.path(demo, "pilots", manifest$file[i]))
  stopifnot(inherits(pilot, "CircadianBioOptions"))
  if (gtex[i]) stopifnot(identical(pilot$time_basis, "DIP"),
                         identical(pilot$gtex_release, "v10"))
}
stopifnot(!any(files$dataset == "GTEx" & files$time_basis != "DIP"))
cat("Demo staged:", normalizePath(demo), "\n",
    nrow(manifest), "default pilots and", nrow(files), "harmonic files.\n")
