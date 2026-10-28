# Prepare GEO expression using a reviewed sample plan.
source("examples/pilots/native/00_setup.R")
source("examples/pilots/native/geo_helpers.R")
stopifnot(requireNamespace("GEOquery", quietly = TRUE),
          requireNamespace("Biobase", quietly = TRUE),
          requireNamespace("matrixStats", quietly = TRUE))

# 1. Read the exact sample/time choices ----------------------------------
# The plan makes the series platform, arm, and input time explicit.
plan_file <- native_input("SCP_NATIVE_SAMPLE_PLAN")
plan <- read.csv(plan_file, stringsAsFactors = FALSE)
stopifnot(all(c("base_file", "sample", "time", "eset_index") %in% names(plan)))
geo <- unique(native_index[native_index$source_group == "GEO",
                           c("accession", "base_file", "preprocessing")])
dir.create(native_cache, recursive = TRUE, showWarnings = FALSE)

# 2. Load the specified series and align columns by sample key ------------
for (accession in unique(geo$accession)) {
  cache_dir <- file.path(native_cache, accession)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  series <- GEOquery::getGEO(accession, GSEMatrix = TRUE,
                             destdir = cache_dir, AnnotGPL = FALSE)
  rows <- geo[geo$accession == accession, , drop = FALSE]
  for (i in seq_len(nrow(rows))) {
    selected <- plan[plan$base_file == rows$base_file[i], , drop = FALSE]
    if (!nrow(selected)) stop("The sample plan has no rows for ", rows$base_file[i])
    stopifnot(!anyDuplicated(selected$sample), all(is.finite(selected$time)),
              length(unique(selected$eset_index)) == 1L)
    eset <- series[[as.integer(selected$eset_index[1])]]
    expression <- Biobase::exprs(eset)
    if (!nrow(expression)) {
      source_file <- unique(selected$expression_file)
      if (length(source_file) != 1L || !file.exists(source_file))
        stop("Provide the exact supplementary expression_file for ", accession)
      expression <- native_read_expression(source_file)
    } else {
      source_file <- list.files(cache_dir, pattern = "series_matrix.*gz$",
                                full.names = TRUE)
    }
    matched <- match(selected$sample, colnames(expression))
    if (anyNA(matched)) stop("Expression columns do not match the reviewed sample plan.")
    expr <- expression[, matched, drop = FALSE]

    # 3. Retain the original GEO preprocessing convention -----------------
    keep <- rowSums(is.na(expr)) < ncol(expr) / 2 &
      matrixStats::rowVars(expr, na.rm = TRUE) > 0
    expr <- expr[which(keep), , drop = FALSE]
    if (max(expr, na.rm = TRUE) > 1e4 && all(expr[is.finite(expr)] >= 0)) {
      library_size <- colSums(expr, na.rm = TRUE)
      stopifnot(all(library_size > 0))
      expr <- log2(t(t(expr) / library_size) * 1e6 + 1)
    }
    symbols <- native_probe_symbols(Biobase::fData(eset))
    native_save_input(expr, selected$time, rows$base_file[i],
                      c(sample_plan = plan_file, source_file),
                      rows$preprocessing[i], symbols)
  }
}
