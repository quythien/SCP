# Shared inputs for the native-time pilot recipes.

# 1. Use the installed package and choose a staging directory -------------
stopifnot(requireNamespace("SCP", quietly = TRUE),
          requireNamespace("digest", quietly = TRUE))
native_stage <- Sys.getenv("SCP_NATIVE_STAGE", "output/pilot_rebuild/native")
native_cache <- Sys.getenv("SCP_NATIVE_CACHE", file.path(native_stage, "cache_geo"))
native_index_file <- Sys.getenv("SCP_NATIVE_INDEX",
                              "examples/pilots/builder_source_index.csv")
dir.create(native_stage, recursive = TRUE, showWarnings = FALSE)
native_stage <- normalizePath(native_stage)
installed_pilots <- system.file("extdata", "pilots", package = "SCP")
source_pilots <- normalizePath("inst/extdata/pilots", mustWork = FALSE)
for (protected in c(installed_pilots, source_pilots)) {
  if (nzchar(protected) &&
      (identical(native_stage, protected) ||
       startsWith(native_stage, paste0(protected, .Platform$file.sep)) ||
       identical(normalizePath(file.path(native_stage, "pilots"), mustWork = FALSE), protected)))
    stop("Choose a staging directory outside the package pilot database.")
}
native_index <- read.csv(native_index_file, stringsAsFactors = FALSE)
selected <- trimws(strsplit(Sys.getenv("SCP_NATIVE_DATASETS", ""), ",")[[1]])
if (length(selected) && any(nzchar(selected)))
  native_index <- native_index[native_index$dataset %in% selected, , drop = FALSE]

# 2. Locate explicit local inputs ----------------------------------------
native_input <- function(variable) {
  path <- Sys.getenv(variable)
  if (!nzchar(path) || !file.exists(path))
    stop("Set ", variable, " to the required local source file.")
  path
}

# 3. Save prepared expression and times without touching bundled pilots ---
native_save_input <- function(expr, times, base_file, sources, notes,
                              symbol_map = NULL) {
  expr <- as.matrix(expr)
  storage.mode(expr) <- "numeric"
  times <- as.numeric(times)
  stopifnot(length(times) == ncol(expr), ncol(expr) >= 6L,
            all(is.finite(times)), !is.null(rownames(expr)),
            !anyDuplicated(rownames(expr)))
  hashes <- vapply(sources, function(path)
    digest::digest(file = path, algo = "sha256"), character(1), USE.NAMES = FALSE)
  labels <- names(sources)
  if (is.null(labels)) labels <- rep("", length(sources))
  labels[!nzchar(labels)] <- paste0("input", which(!nzchar(labels)))
  names(hashes) <- labels
  out <- file.path(native_stage, "inputs", base_file)
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(expr = expr, times = times, symbol_map = symbol_map,
               source_hashes = hashes, notes = notes), out, compress = "xz")
  invisible(out)
}
