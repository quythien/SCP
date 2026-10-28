# Prepare human datasets that require an explicit local expression/time input.
source("examples/pilots/native/00_setup.R")
stopifnot(requireNamespace("matrixStats", quietly = TRUE))

# 1. GSE160521 control striatum ------------------------------------------
if (any(native_index$dataset == "GSE160521")) {
  counts_file <- native_input("SCP_STRIATUM_COUNTS")
  metadata_file <- native_input("SCP_STRIATUM_METADATA")
  counts <- as.matrix(read.csv(gzfile(counts_file), row.names = 1,
                               check.names = FALSE))
  metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)
  stopifnot(all(c("HU", "Diagnostic.Category", "CorrectedTOD") %in% names(metadata)),
            !anyDuplicated(metadata$HU))
  regions <- c(Caudate = "C", NAc = "N", Putamen = "P")
  rows <- unique(native_index[native_index$dataset == "GSE160521",
                              c("tissue", "base_file")])
  for (i in seq_len(nrow(rows))) {
    columns <- grep(paste0(regions[rows$tissue[i]], "$"), colnames(counts))
    expr <- counts[, columns, drop = FALSE]
    subject <- as.integer(gsub("[A-Za-z]+$", "", colnames(expr)))
    matched <- match(subject, metadata$HU)
    keep <- which(!is.na(matched) &
                    metadata$Diagnostic.Category[matched] == "CONTROL" &
                    is.finite(metadata$CorrectedTOD[matched]))
    expr <- expr[, keep, drop = FALSE]
    times <- metadata$CorrectedTOD[matched[keep]] %% 24
    keep_gene <- rowSums(is.na(expr)) < ncol(expr) / 2 &
      matrixStats::rowVars(expr, na.rm = TRUE) > 0
    expr <- expr[which(keep_gene), , drop = FALSE]
    library_size <- colSums(expr)
    stopifnot(all(is.finite(library_size)), all(library_size > 0))
    expr <- log2(t(t(expr) / library_size) * 1e6 + 1)
    native_save_input(expr, times, rows$base_file[i],
                      c(counts = counts_file, metadata = metadata_file),
                      "Control samples matched by supplied subject key. Library-size CPM then log2(CPM + 1).")
  }
}

# 2. GSE71620 prefrontal cortex ------------------------------------------
if (any(native_index$dataset == "GSE71620")) {
  source_file <- native_input("SCP_PFC_INPUT")
  combined <- readRDS(source_file)
  matched <- match(colnames(combined$expr), combined$pheno$sample_name)
  keep <- which(!is.na(matched))
  expr <- combined$expr[, keep, drop = FALSE]
  metadata <- combined$pheno[matched[keep], , drop = FALSE]
  times <- if ("TOD.x" %in% names(metadata)) metadata$TOD.x else metadata$TOD.y
  age_group <- if ("AgeGroup" %in% names(metadata)) metadata$AgeGroup else metadata$age_group
  rows <- unique(native_index[native_index$dataset == "GSE71620",
                              c("tissue", "condition", "base_file")])
  for (i in seq_len(nrow(rows))) {
    keep <- is.finite(times) & !is.na(metadata$region)
    if (rows$condition[i] %in% c("Younger", "Older")) {
      keep <- keep & !is.na(age_group) & age_group == tolower(rows$condition[i])
    } else {
      keep <- keep & metadata$region == rows$tissue[i]
    }
    native_save_input(expr[, which(keep), drop = FALSE], times[which(keep)],
                      rows$base_file[i], source_file,
                      "Prepared expression and matched phenotype table. Original TOD-column precedence and age/region selections retained.")
  }
}
