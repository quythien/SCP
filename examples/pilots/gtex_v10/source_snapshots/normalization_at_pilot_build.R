# Prepare GTEx v10 using the published CHIRAL preprocessing functions.
# The public subject file supplies age groups, sex, and death type.

# 1. Read the public inputs and the pinned author source -------------------
stage <- Sys.getenv("SCP_GTEX_V10_STAGE", "output/gtex_v10_chiral_protocol")
full_path <- Sys.getenv("SCP_GTEX_V10_FULL")
qc_path <- Sys.getenv("SCP_GTEX_V10_SAMPLE_QC")
subject_path <- Sys.getenv("SCP_GTEX_V10_SUBJECTS")
author_source <- Sys.getenv("SCP_CHIRAL_SOURCE")
stopifnot(nzchar(full_path), file.exists(qc_path), file.exists(subject_path),
          dir.exists(author_source))
prepared <- file.path(stage, "prepared_inputs")
dir.create(prepared, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(stage, "private", "normalization_checkpoints"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(stage, "qc"), showWarnings = FALSE)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) RhpcBLASctl::blas_set_num_threads(1L)
workers <- max(1L, as.integer(Sys.getenv("SCP_GTEX_PREPROCESS_WORKERS", "4")))
source_file <- file.path(author_source, "Paper", "To_CPM.R")
author <- new.env(parent = globalenv())
author$mclapply <- parallel::mclapply
for (expression in parse(source_file)) {
  if (is.call(expression) && is.symbol(expression[[1]]) &&
      as.character(expression[[1]]) %in% c("<-", "=") && is.symbol(expression[[2]]) &&
      as.character(expression[[2]]) %in% c("run_edgeR", "Norm.CPM", "remove_covariates", "spliti"))
    eval(expression, envir = author)
}

# Verify downloaded files against the public source index before using them.
public_inputs <- jsonlite::read_json("examples/pilots/gtex_v10/public_inputs.json", simplifyVector = TRUE)
verify_public_file <- function(path, record) {
  stopifnot(file.exists(path), !is.null(record$sha256), !is.null(record$bytes))
  hash <- digest::digest(file = path, algo = "sha256")
  stopifnot(identical(hash, record$sha256), file.info(path)$size == record$bytes)
  hash
}
verified_inputs <- list(
  sample_QC_sha256 = verify_public_file(qc_path, public_inputs$inputs$sample_QC),
  public_subjects_sha256 = verify_public_file(subject_path, public_inputs$inputs$public_subjects))

# 2. Use public sample QC and the original public age-group categories ------
fields <- c("SAMPID", "SMRIN", "SMTSISCH", "SMMAPRT", "SMUNMPRT", "SMMPPD",
            "SMNABTCH", "SMATSSCR", "SMTSPAX", "SMPTHNTS")
samp <- as.data.frame(data.table::fread(qc_path, select = fields, showProgress = FALSE))
subjects <- as.data.frame(data.table::fread(subject_path, showProgress = FALSE))
stopifnot(all(c("SUBJID", "AGE", "SEX", "DTHHRDY") %in% names(subjects)),
          is.character(subjects$AGE), !anyDuplicated(subjects$SUBJID))
# Reuse the cached author-scale CPM when it is available. A first build can
# instead supply the official public reads file and run the same run_edgeR.
if (file.exists(full_path)) {
  cache_sha <- digest::digest(file = full_path, algo = "sha256")
  cache_record <- public_inputs$normalized_CPM_cache
  if (!is.null(cache_record) && identical(cache_sha, cache_record$sha256)) {
    stopifnot(file.info(full_path)$size == cache_record$bytes,
      identical(cache_record$raw_gene_reads_sha256, public_inputs$inputs$gene_reads$sha256))
  } else {
    receipt_path <- paste0(full_path, ".chiral_cache_receipt.json")
    stopifnot(file.exists(receipt_path))
    cache_record <- jsonlite::read_json(receipt_path, simplifyVector = TRUE)
    stopifnot(identical(cache_sha, cache_record$sha256),
      identical(cache_record$raw_gene_reads_sha256, public_inputs$inputs$gene_reads$sha256),
      identical(cache_record$source_sample_QC_sha256, verified_inputs$sample_QC_sha256),
      identical(cache_record$author_To_CPM_sha256, digest::digest(file = source_file, algo = "sha256")))
  }
  source_env <- new.env()
  load(full_path, envir = source_env)
  full <- source_env$CPM.all
  rm(source_env)
} else {
  reads_path <- Sys.getenv("SCP_GTEX_V10_READS")
  stopifnot(file.exists(reads_path), requireNamespace("edgeR", quietly = TRUE))
  verified_inputs$gene_reads_sha256 <- verify_public_file(reads_path, public_inputs$inputs$gene_reads)
  author$DGEList <- edgeR::DGEList
  author$calcNormFactors <- edgeR::calcNormFactors
  author$cpm <- edgeR::cpm
  author$gtex <- as.data.frame(data.table::fread(reads_path, skip = 2L, showProgress = FALSE))
  names(author$gtex) <- gsub(".", "-", names(author$gtex), fixed = TRUE)
  tissue_metadata <- data.table::fread(qc_path, select = c("SAMPID", "SMTSD"), showProgress = FALSE)
  author$gtex.tissue <- tissue_metadata$SMTSD[match(names(author$gtex), tissue_metadata$SAMPID)]
  labels <- unique(author$gtex.tissue[!is.na(author$gtex.tissue)])
  cpm_values <- parallel::mclapply(labels, author$run_edgeR, mc.cores = workers,
    mc.preschedule = FALSE)
  stopifnot(!any(vapply(cpm_values, inherits, logical(1), "try-error")))
  full <- setNames(lapply(cpm_values, function(x) {
    colnames(x) <- gsub(".", "-", colnames(x), fixed = TRUE)
    list(counts = as.matrix(x))
  }), gsub("^_|_$", "", gsub("[^[:alnum:]]+", "_", labels)))
  CPM.all <- full
  dir.create(dirname(full_path), recursive = TRUE, showWarnings = FALSE)
  save(CPM.all, file = full_path, compress = "gzip")
  cache_record <- list(sha256 = digest::digest(file = full_path, algo = "sha256"),
    bytes = file.info(full_path)$size, raw_gene_reads_sha256 = verified_inputs$gene_reads_sha256,
    source_sample_QC_sha256 = verified_inputs$sample_QC_sha256,
    author_To_CPM_sha256 = digest::digest(file = source_file, algo = "sha256"),
    normalization = "Unchanged author run_edgeR with v10 public inputs",
    edgeR_version = as.character(utils::packageVersion("edgeR")),
    R_version = R.version.string, gene_mean_min_exclusive = 10,
    method = "TMM and edgeR::cpm(log=TRUE, prior.count=0.25)")
  jsonlite::write_json(cache_record, paste0(full_path, ".chiral_cache_receipt.json"),
    auto_unbox = TRUE, pretty = TRUE)
  rm(CPM.all, cpm_values)
  rm(list = c("gtex", "gtex.tissue"), envir = author)
}
tissues <- names(full)[!grepl("^Cells", names(full))]
source_hashes <- list(full_CPM_sha256 = digest::digest(file = full_path, algo = "sha256"),
  sample_QC_sha256 = digest::digest(file = qc_path, algo = "sha256"),
  public_subjects_sha256 = digest::digest(file = subject_path, algo = "sha256"),
  author_To_CPM_sha256 = digest::digest(file = source_file, algo = "sha256"))
signature <- digest::digest(source_hashes, algo = "sha256")
jsonlite::write_json(list(verified_inputs = verified_inputs, CPM_cache = cache_record),
  file.path(stage, "qc", "public_input_verification.json"), auto_unbox = TRUE, pretty = TRUE)

# 3. Run the unchanged author function and checkpoint each tissue -----------
prepare_tissue <- function(tissue, kind) {
  path <- file.path(stage, "private", "normalization_checkpoints", paste0(tissue, "_", kind, ".rds"))
  if (file.exists(path)) {
    old <- readRDS(path)
    stopifnot(identical(old$signature, signature))
    return(old)
  }
  original <- full[[tissue]]$counts
  frame <- as.data.frame(original)
  # The author's data.frame CPM has dotted sample names. Restore the v10
  # names after preprocessing so the later donor join uses exact IDs.
  colnames(frame) <- gsub("-", ".", colnames(frame), fixed = TRUE)
  threshold <- if (kind == "strict") 3 else 0
  frame <- frame[rowMeans(frame) > threshold, , drop = FALSE]
  tissue_samp <- samp[samp$SAMPID %in% colnames(original), , drop = FALSE]
  eligible <- if (kind == "strict") subset(tissue_samp, SMRIN > 6 & SMMPPD > 40000000 &
    SMMAPRT > 0.8 & SMTSISCH > 0 & SMATSSCR %in% 0:1) else subset(tissue_samp, SMRIN > 4)
  # The author's >48-sample guard also handles the v10 one-column tissues.
  # Check it before data.frame subsetting can drop a single column to a vector.
  blocks <- if (nrow(eligible) <= 48L) list() else
    split(seq_len(nrow(frame)), ceiling(seq_len(nrow(frame)) / 250L))
  began <- proc.time()[["elapsed"]]
  adjusted_blocks <- lapply(seq_along(blocks), function(block) {
    input <- setNames(list(frame[blocks[[block]], , drop = FALSE]), tissue)
    result <- author$Norm.CPM(input, high_filter = identical(kind, "strict"),
      ncores = 1L, samp = tissue_samp, samp.2 = subjects)[[1]]
    if (block == 1L || block %% 20L == 0L)
      cat(tissue, kind, "gene block", block, "of", length(blocks), "\n")
    flush.console()
    result
  })
  adjusted <- if (all(vapply(adjusted_blocks, is.null, logical(1)))) NULL else
    do.call(rbind, adjusted_blocks)
  if (is.null(adjusted)) result <- list(data = NULL, n_before_complete_cases = 0L,
    n_complete_samples = 0L, n_genes = 0L) else {
    good <- is.finite(as.numeric(adjusted[1, ]))
    expr <- as.matrix(adjusted[, good, drop = FALSE])
    sample_map <- setNames(colnames(original), gsub("-", ".", colnames(original), fixed = TRUE))
    colnames(expr) <- unname(sample_map[colnames(expr)])
    stopifnot(!anyNA(colnames(expr)), all(is.finite(expr)))
    metadata <- samp[match(colnames(expr), samp$SAMPID), , drop = FALSE]
    metadata$SUBJID <- sub("^(GTEX-[^-]+)-.*", "\\1", metadata$SAMPID)
    metadata <- data.frame(metadata, subjects[match(metadata$SUBJID, subjects$SUBJID),
      c("AGE", "SEX", "DTHHRDY"), drop = FALSE], row.names = NULL)
    result <- list(data = list(counts = expr, metadata = metadata,
       gene_annotation = data.frame(gene_id = rownames(expr),
         ensembl = sub("_.*$", "", rownames(expr)),
         symbol = sub("^.*_", "", rownames(expr)), stringsAsFactors = FALSE)),
      n_before_complete_cases = ncol(adjusted), n_complete_samples = ncol(expr), n_genes = nrow(expr))
  }
  result$tissue <- tissue
  result$kind <- kind
  result$signature <- signature
  result$elapsed_seconds <- proc.time()[["elapsed"]] - began
  saveRDS(result, path, compress = "gzip")
  cat(tissue, kind, result$n_complete_samples, "complete samples", result$n_genes,
      "genes", round(result$elapsed_seconds, 1), "seconds\n")
  flush.console()
  result
}

# 4. Save strict inputs first so phase inference can begin ------------------
inventory <- list()
for (kind in c("strict", "broad")) {
  results <- parallel::mclapply(tissues, prepare_tissue, kind = kind,
    mc.cores = workers, mc.preschedule = FALSE, mc.set.seed = FALSE)
  failed <- vapply(results, inherits, logical(1), "try-error")
  if (any(failed)) {
    for (i in which(failed)) cat(tissues[i], as.character(results[[i]]), "\n")
    stop("Published preprocessing failed for the listed tissues.")
  }
  values <- setNames(lapply(results, function(x) x$data), tissues)
  values <- Filter(Negate(is.null), values)
  if (kind == "strict") {
    CPM.all.norm <- values
    destination <- file.path(prepared, "v10_CPM.all.norm.chiral.RData")
    save(CPM.all.norm, file = destination, compress = "gzip")
  } else {
    CPM.all.norm_large <- values
    destination <- file.path(prepared, "v10_CPM.all.norm_large.chiral.RData")
    save(CPM.all.norm_large, file = destination, compress = "gzip")
  }
  inventory[[kind]] <- do.call(rbind, lapply(results, function(x) data.frame(
    tissue = x$tissue, kind = kind, n_before_complete_cases = x$n_before_complete_cases,
    n_complete_samples = x$n_complete_samples, ngenes = x$n_genes,
    elapsed_seconds = x$elapsed_seconds)))
  jsonlite::write_json(list(kind = kind, n_tissues = length(values),
    file = basename(destination), sha256 = digest::digest(file = destination, algo = "sha256")),
    file.path(stage, "qc", paste0(kind, "_ready.json")), auto_unbox = TRUE, pretty = TRUE)
  cat(kind, "author-protocol inputs ready\n")
  flush.console()
}

# 5. Record the published recipe and exact source hashes --------------------
utils::write.csv(do.call(rbind, inventory), file.path(stage, "qc", "normalization_inventory.csv"), row.names = FALSE)
jsonlite::write_json(c(source_hashes, list(chiral_revision = "1f78a9a60b6d18facfedf5923d43bae9fa79f0b0",
  upstream_functions = c("Norm.CPM", "remove_covariates"),
  raw_normalization = "Cached public TMM log CPM from edgeR cpm(log=TRUE, prior.count=0.25). Raw mean count >10.",
  strict_gene_mean_threshold = 3, broad_gene_mean_threshold = 0,
  gene_filter_before_sample_QC = TRUE,
  covariates = "isch=cut(SMTSISCH,c(0,200,500,800,2000)), public AGE groups, SEX factor, DTHHRDY factor",
  formula = "y ~ isch + age + sex + death.type. Sex omitted if it has one level, as in the author code.",
  minimum_pre_regression_samples = 49L,
  complete_case_rule = "The author lm residuals retain missing columns as NA. CPM_to_E removes columns missing in the first gene.",
  identifier_adaptation = "Dotted author CPM sample names are mapped back to the exact original v10 sample IDs.",
  gene_batch_size = 250L,
  batching = "The unchanged author gene-wise function is evaluated in blocks with the same tissue sample metadata.",
  private_exact_age_used = FALSE)), file.path(stage, "qc", "normalization_recipe.json"),
  auto_unbox = TRUE, pretty = TRUE, digits = 16)
