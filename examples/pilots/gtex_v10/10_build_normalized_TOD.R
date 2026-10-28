# Build GTEx v10 pilots from normalized log CPM and seasonal TOD.

# 1. Read the inputs and choose a new output directory --------------------
stopifnot(requireNamespace("digest", quietly = TRUE),
          requireNamespace("jsonlite", quietly = TRUE))
stage <- Sys.getenv("SCP_GTEX_TOD_STAGE", "output/gtex_v10_normalized_TOD")
full_path <- Sys.getenv("SCP_GTEX_V10_FULL")
phenotype_path <- Sys.getenv("SCP_GTEX_PHENO")
stopifnot(nzchar(phenotype_path), file.exists(phenotype_path))
if (dir.exists(stage) && length(list.files(stage, all.files = TRUE, no.. = TRUE)))
  stop("The output directory contains files. Choose a new SCP_GTEX_TOD_STAGE.")
for (directory in c("private/cohorts", "qc", "summary", "inst/extdata/pilots"))
  dir.create(file.path(stage, directory), recursive = TRUE, showWarnings = FALSE)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) RhpcBLASctl::blas_set_num_threads(1L)
for (file in c("R/options.R", "R/estimation.R")) sys.source(file, envir = globalenv())
file_hash <- function(path) digest::digest(file = path, algo = "sha256")
json_write <- function(value, path) jsonlite::write_json(value, path,
  auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null", digits = 17)
public_inputs <- jsonlite::read_json("examples/pilots/gtex_v10/public_inputs.json",
                                   simplifyVector = TRUE)
phenotype <- readRDS(phenotype_path)
stopifnot(all(c("SUBJID", "CT") %in% names(phenotype)),
          !anyDuplicated(phenotype$SUBJID), is.numeric(phenotype$CT))
seasonal_time <- setNames(phenotype$CT %% 24, phenotype$SUBJID)
source_manifest <- utils::read.csv("inst/extdata/pilots/manifest.csv",
                                   stringsAsFactors = FALSE)
targets <- source_manifest[source_manifest$dataset == "GTEx" &
  source_manifest$status == "ingested" & source_manifest$time_basis == "TOD", , drop = FALSE]
native <- source_manifest[source_manifest$dataset != "GTEx", , drop = FALSE]
stopifnot(nrow(targets) == 46L, nrow(native) == 116L)
targets$source_tissue <- ifelse(targets$tissue == "Adrenal", "Adrenal_Gland", targets$tissue)
report_tissues <- c("Adrenal_Gland", "Liver", "Muscle_Skeletal", "Thyroid", "Lung", "Pancreas",
                   "Heart_Left_Ventricle", "Artery_Tibial", "Adipose_Subcutaneous", "Brain_Cortex")
priority <- match(targets$source_tissue, c("Liver", report_tissues))
targets <- targets[order(is.na(priority), priority, targets$tissue), , drop = FALSE]
pilot_root <- file.path(stage, "inst", "extdata", "pilots")
recipe <- list(release = "v10", period = 24, pilot_seed = 2025L,
  alpha_pilot = 0.01, top_k = 300L, paired_sigma = TRUE,
  normalization = "Within-tissue TMM log CPM with edgeR prior.count=0.25",
  raw_gene_mean_min_exclusive = 10, normalized_gene_mean_min_exclusive = 3,
  gene_filter_order = "First valid occurrence of each gene symbol, then mean log CPM >3 on the finite-TOD cohort",
  input_time = "Season-adjusted time of death (TOD)",
  cohort = "All source samples with a finite donor-matched seasonal TOD")
source_hashes <- list(public_gene_reads_sha256 = public_inputs$inputs$gene_reads$sha256,
  public_sample_QC_sha256 = public_inputs$inputs$sample_QC$sha256,
  seasonal_phenotype_sha256 = file_hash(phenotype_path),
  pilot_recipe_sha256 = file_hash("examples/pilots/gtex_v10/10_build_normalized_TOD.R"),
  options_source_sha256 = file_hash("R/options.R"),
  estimation_source_sha256 = file_hash("R/estimation.R"))

# 2. Load verified normalized CPM or rebuild it from the public raw reads --
if (nzchar(full_path) && file.exists(full_path)) {
  expression_sha <- file_hash(full_path)
  if (identical(expression_sha, public_inputs$normalized_CPM_cache$sha256)) {
    stopifnot(file.info(full_path)$size == public_inputs$normalized_CPM_cache$bytes)
    cache_mode <- "Verified normalized cache"
  } else {
    receipt_path <- paste0(full_path, ".normalized_TOD_receipt.json")
    stopifnot(file.exists(receipt_path))
    cache_receipt <- jsonlite::read_json(receipt_path, simplifyVector = TRUE)
    stopifnot(identical(expression_sha, cache_receipt$sha256),
      file.info(full_path)$size == cache_receipt$bytes,
      identical(cache_receipt$public_gene_reads_sha256, public_inputs$inputs$gene_reads$sha256),
      identical(cache_receipt$public_sample_QC_sha256, public_inputs$inputs$sample_QC$sha256),
      identical(cache_receipt$builder_sha256, source_hashes$pilot_recipe_sha256),
      cache_receipt$raw_gene_mean_min_exclusive == 10,
      cache_receipt$prior_count == 0.25, identical(cache_receipt$method, "TMM log CPM"))
    cache_mode <- "Verified rebuilt normalized cache"
  }
  environment <- new.env()
  load(full_path, envir = environment)
  stopifnot(exists("CPM.all", environment, inherits = FALSE))
  normalized <- environment$CPM.all
  rm(environment)
} else {
  stopifnot(requireNamespace("edgeR", quietly = TRUE),
            requireNamespace("data.table", quietly = TRUE))
  reads_path <- Sys.getenv("SCP_GTEX_V10_READS")
  qc_path <- Sys.getenv("SCP_GTEX_V10_SAMPLE_QC")
  stopifnot(file.exists(reads_path), file.exists(qc_path),
    identical(file_hash(reads_path), public_inputs$inputs$gene_reads$sha256),
    identical(file_hash(qc_path), public_inputs$inputs$sample_QC$sha256))
  data.table::setDTthreads(1L)
  reads <- data.table::fread(reads_path, skip = 2L, showProgress = FALSE)
  sample_metadata <- data.table::fread(qc_path, select = c("SAMPID", "SMTSD"),
                                     showProgress = FALSE)
  sample_names <- names(reads)[-(1:2)]
  tissue_names <- sample_metadata$SMTSD[match(sample_names, sample_metadata$SAMPID)]
  normalized <- list()
  for (label in unique(tissue_names[!is.na(tissue_names)])) {
    samples <- sample_names[!is.na(tissue_names) & tissue_names == label]
    counts <- as.matrix(reads[, samples, with = FALSE])
    rownames(counts) <- paste(reads$Name, reads$Description, sep = "_")
    counts <- counts[rowMeans(counts) > 10, , drop = FALSE]
    dge <- edgeR::calcNormFactors(edgeR::DGEList(counts = counts))
    log_cpm <- as.matrix(edgeR::cpm(dge, log = TRUE, prior.count = 0.25))
    gene_id <- rownames(log_cpm)
    donor <- sub("^GTEX[.-]([A-Z0-9]+).*", "GTEX-\\1", colnames(log_cpm))
    key <- gsub("^_|_$", "", gsub("[^[:alnum:]]+", "_", label))
    normalized[[key]] <- list(counts = log_cpm,
      gene_annotation = data.frame(gene_id = gene_id,
        symbol = as.character(reads$Description[match(gene_id,
          paste(reads$Name, reads$Description, sep = "_"))]), stringsAsFactors = FALSE),
      tod = unname(seasonal_time[donor]))
  }
  normalized <- normalized[order(names(normalized))]
  CPM.all <- normalized
  full_path <- file.path(stage, "private", "v10_CPM_full.rebuilt.RData")
  stopifnot(!file.exists(full_path))
  save(CPM.all, file = full_path, compress = "gzip")
  expression_sha <- file_hash(full_path)
  json_write(list(sha256 = expression_sha, bytes = file.info(full_path)$size,
    public_gene_reads_sha256 = public_inputs$inputs$gene_reads$sha256,
    public_sample_QC_sha256 = public_inputs$inputs$sample_QC$sha256,
    builder_sha256 = source_hashes$pilot_recipe_sha256,
    raw_gene_mean_min_exclusive = 10, prior_count = 0.25, method = "TMM log CPM",
    edgeR_version = as.character(utils::packageVersion("edgeR")), R_version = R.version.string),
    paste0(full_path, ".normalized_TOD_receipt.json"))
  rm(CPM.all, reads, sample_metadata, counts, dge, log_cpm)
  cache_mode <- "Rebuilt from verified public raw reads using this script"
}
source_hashes$expression_sha256 <- expression_sha
stopifnot(all(targets$source_tissue %in% names(normalized)))
json_write(c(recipe, source_hashes, list(cache_mode = cache_mode,
  cache_path = normalizePath(full_path), phenotype_path = normalizePath(phenotype_path),
  active_manifest_sha256 = file_hash("inst/extdata/pilots/manifest.csv"),
  R_version = R.version.string)), file.path(stage, "private", "build_inputs.json"))
cat("Normalized CPM and seasonal TOD verified.\n")
flush.console()

# 3. Fit full-universe models independently for validation and exploration -
full_gene_tables <- function(expr, times, symbol) {
  n <- length(times)
  x1 <- cbind(1, cos(2 * pi * times / 24), sin(2 * pi * times / 24))
  x2 <- cbind(x1, cos(4 * pi * times / 24), sin(4 * pi * times / 24))
  fit1 <- stats::lm.fit(x1, t(expr))
  fit2 <- stats::lm.fit(x2, t(expr))
  stopifnot(n > 5L, fit1$rank == 3L, fit2$rank == 5L)
  rss1 <- colSums(fit1$residuals^2)
  rss2 <- colSums(fit2$residuals^2)
  tss <- rowSums((expr - rowMeans(expr))^2)
  p1 <- stats::pf((pmax(tss - rss1, 0) / 2) / (rss1 / (n - 3)), 2, n - 3,
                 lower.tail = FALSE)
  p2 <- stats::pf((pmax(tss - rss2, 0) / 4) / (rss2 / (n - 5)), 4, n - 5,
                 lower.tail = FALSE)
  ph <- stats::pf((pmax(rss1 - rss2, 0) / 2) / (rss2 / (n - 5)), 2, n - 5,
                 lower.tail = FALSE)
  first <- data.frame(gene = rownames(expr), symbol = symbol, pvalue = p1,
    qvalue = stats::p.adjust(p1, "BH", n = nrow(expr)), p_K1 = p1,
    q_K1 = stats::p.adjust(p1, "BH", n = nrow(expr)),
    A = sqrt(fit1$coefficients[2, ]^2 + fit1$coefficients[3, ]^2),
    phi = (atan2(fit1$coefficients[3, ], fit1$coefficients[2, ]) * 12 / pi) %% 24,
    sigma = sqrt(rss1 / (n - 3)), mesor = fit1$coefficients[1, ], n = n,
    stringsAsFactors = FALSE, row.names = NULL)
  second <- data.frame(gene = first$gene, symbol = symbol, pvalue = p2,
    qvalue = stats::p.adjust(p2, "BH", n = nrow(expr)), p_joint = p2,
    q_joint = stats::p.adjust(p2, "BH", n = nrow(expr)), p_K1 = p1, q_K1 = first$q_K1,
    p_2h = ph, q_2h = stats::p.adjust(ph, "BH", n = nrow(expr)),
    A = sqrt(fit2$coefficients[2, ]^2 + fit2$coefficients[3, ]^2),
    phi = (atan2(fit2$coefficients[3, ], fit2$coefficients[2, ]) * 12 / pi) %% 24,
    A2 = sqrt(fit2$coefficients[4, ]^2 + fit2$coefficients[5, ]^2),
    phi2 = (atan2(fit2$coefficients[5, ], fit2$coefficients[4, ]) * 6 / pi) %% 12,
    sigma = sqrt(rss2 / (n - 5)), mesor = fit2$coefficients[1, ], n = n,
    stringsAsFactors = FALSE, row.names = NULL)
  list(first, second)
}
sampling_distribution <- function(times) {
  den <- stats::density(c(times - 24, times, times + 24), from = 0, to = 24)
  structure(list(grid = unname(den$x), weights = unname(pmax(den$y, 1e-10)),
    period = 24, jitter_sd = 0.1, n_observed = length(times)), class = "SCPTimeDistribution")
}
manifest_rows <- summary_rows <- build_rows <- sampling_rows <- check_rows <- cohort_rows <- list()
report_cases <- list()

# 4. Select each original-style cohort and fit the SCP pilot models --------
for (i in seq_len(nrow(targets))) {
  row <- targets[i, , drop = FALSE]
  tissue <- row$tissue
  item <- normalized[[row$source_tissue]]
  expr_all <- as.matrix(item$counts)
  stopifnot(is.numeric(expr_all), !anyDuplicated(colnames(expr_all)),
            !anyDuplicated(rownames(expr_all)), all(is.finite(expr_all)))
  samples <- gsub(".", "-", colnames(expr_all), fixed = TRUE)
  stopifnot(all(grepl("^GTEX-[A-Z0-9]+-", samples)))
  donors <- sub("^GTEX-([A-Z0-9]+)-.*", "GTEX-\\1", samples)
  times_all <- unname(seasonal_time[donors])
  eligible <- is.finite(times_all)
  stopifnot(sum(eligible) > 5L, !anyDuplicated(donors[eligible]))
  times <- times_all[eligible] %% 24
  cached_time <- as.numeric(item$tod) %% 24
  time_overlap <- eligible & is.finite(cached_time)
  stopifnot(length(cached_time) == ncol(expr_all), any(time_overlap))
  time_error <- max(abs(((times_all[time_overlap] - cached_time[time_overlap] + 12) %% 24) - 12))
  stopifnot(time_error < 1e-10)
  ann <- item$gene_annotation
  stopifnot(all(c("gene_id", "symbol") %in% names(ann)))
  symbols <- as.character(ann$symbol[match(rownames(expr_all), ann$gene_id)])
  first_symbol <- !duplicated(symbols) & !is.na(symbols) & nzchar(symbols)
  mean_log_cpm <- rowMeans(expr_all[, eligible, drop = FALSE])
  keep_gene <- first_symbol & mean_log_cpm > 3
  gene_ids <- sub("^(ENS[A-Z]*G[0-9]+).*", "\\1", rownames(expr_all)[keep_gene])
  symbols_kept <- symbols[keep_gene]
  stopifnot(all(grepl("^ENSG[0-9]+$", gene_ids)), !anyDuplicated(gene_ids),
            !anyDuplicated(symbols_kept), length(gene_ids) > 0L)
  expr <- expr_all[keep_gene, eligible, drop = FALSE]
  rownames(expr) <- gene_ids
  cohort_rows[[length(cohort_rows) + 1L]] <- data.frame(tissue = tissue,
    n_source = ncol(expr_all), n_used = sum(eligible), n_excluded_TOD = sum(!eligible),
    duplicate_donors = anyDuplicated(donors[eligible]),
    n_cached_TOD_overlap = sum(time_overlap), max_cached_TOD_difference_hours = time_error,
    source_genes = nrow(expr_all), first_valid_symbols = sum(first_symbol),
    genes_after_mean_filter = nrow(expr), stringsAsFactors = FALSE)
  tables <- full_gene_tables(expr, times, symbols_kept)
  full_fit <- tables[[1]]
  candidates <- which(is.finite(full_fit$pvalue) & full_fit$pvalue < recipe$alpha_pilot &
    is.finite(full_fit$A) & full_fit$A > 0 & is.finite(full_fit$sigma) & full_fit$sigma > 0)
  selected <- head(candidates[order(full_fit$pvalue[candidates])], recipe$top_k)
  stopifnot(length(selected) > 0L)
  ratio <- full_fit$A[selected] / full_fit$sigma[selected]
  summary <- data.frame(species = "human", dataset = "GTEx", tissue = tissue, condition = "All",
    gtex_release = "v10", time_basis = "TOD", is_default = TRUE,
    n = length(times), n_source = ncol(expr_all), n_with_TOD = sum(eligible), ngenes = nrow(expr),
    n_selected = length(selected), r_median_top300 = unname(stats::median(ratio)),
    r_q25_top300 = unname(stats::quantile(ratio, 0.25)),
    r_q75_top300 = unname(stats::quantile(ratio, 0.75)),
    prop_rhythmic = length(candidates) / nrow(expr),
    n_BH05 = sum(full_fit$qvalue < 0.05, na.rm = TRUE),
    n_BH05_K2 = sum(tables[[2]]$qvalue < 0.05, na.rm = TRUE), stringsAsFactors = FALSE)
  utils::write.csv(summary, file.path(stage, "qc", paste0(tissue, "_signal.csv")), row.names = FALSE)
  cohort <- list(tissue = tissue, source_tissue = row$source_tissue,
    expression = expr, samples = samples[eligible], donors = donors[eligible],
    TOD = times, genes = gene_ids, gene_symbols = symbols_kept,
    n_source = ncol(expr_all), n_excluded = sum(!eligible),
    source_expression_sha256 = expression_sha, seasonal_phenotype_sha256 = source_hashes$seasonal_phenotype_sha256)
  saveRDS(cohort, file.path(stage, "private", "cohorts", paste0(tissue, ".rds")), compress = "gzip")

  if (tissue == "Liver") {
    liver_check <- data.frame(n_samples = ncol(expr), n_genes = nrow(expr),
      K1 = summary$n_BH05, K2 = summary$n_BH05_K2,
      shared = sum(tables[[1]]$qvalue < 0.05 & tables[[2]]$qvalue < 0.05),
      added_K2 = sum(tables[[1]]$qvalue >= 0.05 & tables[[2]]$qvalue < 0.05),
      lost_K1 = sum(tables[[1]]$qvalue < 0.05 & tables[[2]]$qvalue >= 0.05))
    stopifnot(identical(as.integer(unlist(liver_check[1, ])),
                        c(262L, 8917L, 210L, 725L, 207L, 518L, 3L)))
    utils::write.csv(liver_check, file.path(stage, "qc", "Liver_paper_detection_check.csv"), row.names = FALSE)
    reference_path <- Sys.getenv("SCP_GTEX_LIVER_REFERENCE")
    reference_check <- list(external_TOD_matches_cached_TOD = time_error < 1e-10)
    if (nzchar(reference_path)) {
      reference <- readRDS(reference_path)
      reference_check$expression_matches_reference <- identical(unname(expr), unname(reference$counts))
      reference_check$symbols_match_reference <- identical(symbols_kept, rownames(reference$counts))
      reference_check$times_match_reference <- identical(unname(times), unname(reference$tod))
      reference_check$samples_match_reference <- identical(colnames(expr), reference$sample_ids)
      stopifnot(all(unlist(reference_check)))
    }
    liver_expression <- expr
    colnames(liver_expression) <- sprintf("sample_%03d", seq_len(ncol(expr)))
    saveRDS(list(expression = liver_expression, gene_id = gene_ids, gene_symbol = symbols_kept,
      times = times, source_expression_sha256 = expression_sha,
      source_match_checks = reference_check,
      prefilter_genes = data.frame(gene_id = rownames(expr_all), symbol = symbols,
        original_source_mean_logCPM = mean_log_cpm, first_valid_symbol = first_symbol,
        stringsAsFactors = FALSE)), file.path(stage, "private", "source_Liver_exact.rds"), compress = "gzip")
  }

  row_out <- row[, setdiff(names(row), "source_tissue"), drop = FALSE]
  row_out$n <- length(times)
  row_out$ngenes <- nrow(expr)
  row_out$n_source <- ncol(expr_all)
  row_out$n_excluded <- sum(!eligible)
  row_out$gtex_release <- "v10"
  row_out$time_basis <- "TOD"
  row_out$time_source <- "Season-adjusted time of death (TOD)"
  row_out$is_default <- TRUE
  row_out$has_k2 <- TRUE
  row_out$accession <- "phs000424"
  row_out$citation <- paste("GTEx Consortium (2020).",
    "The GTEx Consortium atlas of genetic regulatory effects across human tissues.",
    "Science 369(6509):1318-1330. doi:10.1126/science.aaz1776.")
  row_out$tod_n_unique <- length(unique(times))
  row_out$tod_cycles <- NA_real_
  row_out$tod_step_hr <- NA_character_
  row_out$tod_phases <- NA_character_
  row_out$note <- ""
  row_out$file <- file.path("human", "gtex_v10", "tod", basename(row$file))

  for (harmonic in 1:2) {
    set.seed(recipe$pilot_seed)
    pilot <- suppressWarnings(if (harmonic == 1L)
      estCircadianParam(expr, times, period = 24, min_rhythm_pval = recipe$alpha_pilot,
        top_k = recipe$top_k, paired_sigma = TRUE, verbose = FALSE) else
      estCircadianParam2H(expr, times, period = 24, min_rhythm_pval = recipe$alpha_pilot,
        top_k = recipe$top_k, paired_sigma = TRUE, verbose = FALSE))
    fit <- pilot$rhythm_fit
    independent <- tables[[harmonic]]
    index <- match(fit$gene, independent$gene)
    stopifnot(!anyNA(index))
    p_reference <- if (harmonic == 1L) independent$pvalue[index] else independent$p_K1[index]
    errors <- c(pvalue = max(abs(fit$pvalue - p_reference)),
      A = max(abs(fit$A - independent$A[index])),
      sigma = max(abs(fit$sigma - independent$sigma[index])),
      mesor = max(abs(fit$mesor - independent$mesor[index])),
      phi = max(abs(((fit$phi - independent$phi[index] + 12) %% 24) - 12)))
    if (harmonic == 2L) errors <- c(errors,
      p_joint = max(abs(fit$p_joint - independent$p_joint[index])),
      A2 = max(abs(fit$A2 - independent$A2[index])),
      phi2 = max(abs(((fit$phi2 - independent$phi2[index] + 6) %% 12) - 6)))
    stopifnot(all(is.finite(errors)), max(errors) < 1e-8,
              abs(pilot$prop_rhythmic - summary$prop_rhythmic) < 1e-12)
    check_rows[[length(check_rows) + 1L]] <- data.frame(tissue = tissue, K = harmonic,
      field = names(errors), maximum_absolute_error = unname(errors),
      tested_genes = nrow(expr), stored_rhythm_genes = nrow(fit), stringsAsFactors = FALSE)
    fit$symbol <- symbols_kept[match(fit$gene, gene_ids)]
    pilot$rhythm_fit <- fit
    chosen <- head(which(fit$pvalue < recipe$alpha_pilot), recipe$top_k)
    stopifnot(length(chosen) == length(selected))
    pilot$amplitude <- pilot$amplitude_spec <- unname(fit$A[chosen])
    pilot$sigma_rhythmic <- unname(fit$sigma[chosen])
    pilot$phase <- pilot$phase_spec <- unname(fit$phi[chosen])
    pilot$paired_sigma <- TRUE
    if (harmonic == 2L) {
      pilot$amplitude2 <- unname(fit$A2[chosen])
      pilot$phase2 <- unname(fit$phi2[chosen])
      pilot$diagnostics$screen_method <- "K1_F2"
    }
    pilot$tested_genes <- gene_ids
    pilot$explorer_fit <- independent
    pilot$times <- pilot$cts <- pilot$raw <- NULL
    pilot$time_distribution <- sampling_distribution(times)
    pilot$gtex_release <- "v10"
    pilot$time_basis <- "TOD"
    pilot$time_source <- row_out$time_source
    pilot$n_samples <- length(times)
    pilot$n_distinct_times <- length(unique(times))
    x <- cbind(1, cos(2 * pi * times / 24), sin(2 * pi * times / 24))
    if (harmonic == 2L) x <- cbind(x, cos(4 * pi * times / 24), sin(4 * pi * times / 24))
    pilot$design_rank <- qr(x)$rank
    stopifnot(pilot$design_rank == 2L * harmonic + 1L)
    pilot$source_provenance <- c(source_hashes, list(dataset = "GTEx", release = "v10",
      accession = "phs000424", data_license = "GTEx Portal Data License",
      license_url = "https://gtexportal.org/home/license",
      acknowledgment = "The data used for these analyses were obtained from the GTEx Portal on 10/28/2026.",
      no_endorsement_implied = TRUE, time_basis = "TOD", time_source = row_out$time_source,
      preprocessing = "Within-tissue TMM log CPM with edgeR prior.count=0.25",
      raw_gene_mean_min_exclusive = 10, normalized_gene_mean_min_exclusive = 3,
      first_valid_symbol_before_gene_mean_filter = TRUE,
      covariate_adjustment = FALSE, tissue_mean_preserved = TRUE,
      cohort = "All source samples with finite donor-matched seasonal TOD",
      n_source = ncol(expr_all), n_used = ncol(expr), n_excluded = sum(!eligible),
      n_genes = nrow(expr), seed = recipe$pilot_seed, period = 24, harmonic_order = harmonic,
      alpha_pilot = recipe$alpha_pilot, top_k = recipe$top_k, paired_sigma = TRUE,
      pilot_screen = "K1_F2", explorer_test = if (harmonic == 1L) "K1_F2" else "K2_F4",
      recipe = "examples/pilots/gtex_v10/10_build_normalized_TOD.R"))
    stopifnot(!length(pilot$times), !length(pilot$cts), !length(pilot$raw),
      !any(grepl("GTEX[-.]", unlist(lapply(pilot[c("rhythm_fit", "explorer_fit")],
        function(table) unlist(table[vapply(table, is.character, logical(1))]))))))
    path <- if (harmonic == 1L) row_out$file else sub("[.]rds$", "_2H.rds", row_out$file)
    destination <- file.path(pilot_root, path)
    dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
    stopifnot(!file.exists(destination))
    saveRDS(pilot, destination, compress = "xz")
    build_rows[[length(build_rows) + 1L]] <- data.frame(file = path, dataset = "GTEx",
      tissue = tissue, K = harmonic, n = length(times), ngenes = nrow(expr),
      sha256 = file_hash(destination), stringsAsFactors = FALSE)
    if (harmonic == 1L && row$source_tissue %in% report_tissues) {
      key <- if (row$source_tissue == "Adrenal_Gland") "Adrenal" else
        if (row$source_tissue == "Muscle_Skeletal") "Muscle" else row$source_tissue
      report_pilot <- pilot
      report_pilot$explorer_fit <- NULL
      report_pilot$ngenes <- 5000L
      report_cases[[paste0(key, "_TOD")]] <- list(bio = report_pilot, times = unname(times),
        tissue = key, time_basis = "TOD", cohort = "all_valid_TOD", n_pilot = length(times),
        genes_fitted = nrow(expr), n_selected = length(selected), medianr = summary$r_median_top300,
        prop = pilot$prop_rhythmic, calibration_law = "installed_loader_top300_v1")
    }
  }
  bins <- hist(times, breaks = 0:24, right = FALSE, include.lowest = TRUE, plot = FALSE)$counts
  stopifnot(length(bins) == 24L, sum(bins) == length(times))
  sampling_rows[[length(sampling_rows) + 1L]] <- data.frame(dataset = "GTEx", tissue = tissue,
    condition = row$condition, gtex_release = "v10", time_basis = "TOD", bin_start = 0:23,
    bin_end = 1:24, count = as.integer(bins), n_samples = length(times), stringsAsFactors = FALSE)
  manifest_rows[[length(manifest_rows) + 1L]] <- row_out
  summary_rows[[length(summary_rows) + 1L]] <- summary
  utils::write.csv(do.call(rbind, tail(check_rows, 2L)),
    file.path(stage, "qc", paste0(tissue, "_fit_validation.csv")), row.names = FALSE)
  cat(tissue, "TOD", ncol(expr), "samples", nrow(expr), "genes", length(selected),
      "selected, median signal-to-noise ratio", format(summary$r_median_top300, digits = 5), "\n")
  flush.console()
  normalized[[row$source_tissue]] <- NULL
  rm(item, expr_all, expr, cohort, tables, full_fit, independent, pilot, fit)
  invisible(gc())
}

# 5. Preserve native studies and save the TOD catalogue and report pilots --
manifest <- rbind(native, do.call(rbind, manifest_rows))
summary <- do.call(rbind, summary_rows)
sampling <- do.call(rbind, sampling_rows)
stopifnot(nrow(manifest) == 162L, nrow(summary) == 46L, nrow(sampling) == 1104L,
          length(report_cases) == 10L, length(build_rows) == 92L)
native_root <- Sys.getenv("SCP_PILOT_SOURCE_ROOT", "inst/extdata/pilots")
native_files <- list.files(native_root, pattern = "[.]rds$", recursive = TRUE)
native_files <- native_files[!grepl("^GTEx_", basename(native_files))]
stopifnot(length(native_files) == 215L)
native_hashes <- list()
for (file in native_files) {
  source <- file.path(native_root, file)
  destination <- file.path(pilot_root, file)
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  stopifnot(file.copy(source, destination, overwrite = FALSE),
            identical(file_hash(source), file_hash(destination)))
  native_hashes[[length(native_hashes) + 1L]] <- data.frame(file = file, dataset = "native",
    K = if (grepl("_2H[.]rds$", file)) 2L else 1L, sha256 = file_hash(destination))
}
utils::write.csv(manifest, file.path(pilot_root, "manifest.csv"), row.names = FALSE)
utils::write.csv(summary, file.path(stage, "summary", "gtex_v10_signal_summary.csv"), row.names = FALSE)
utils::write.csv(sampling, file.path(pilot_root, "gtex_v10_TOD_sampling_histograms.csv"), row.names = FALSE)
inventory <- do.call(rbind, build_rows)
utils::write.csv(inventory, file.path(stage, "qc", "pilot_build_inventory.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, check_rows), file.path(stage, "qc", "fit_validation.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, cohort_rows), file.path(stage, "qc", "cohort_validation.csv"), row.names = FALSE)
hashes <- rbind(inventory[c("file", "dataset", "K", "sha256")], do.call(rbind, native_hashes))
utils::write.csv(hashes, file.path(stage, "pilot_hashes.csv"), row.names = FALSE)
keys <- ifelse(report_tissues == "Adrenal_Gland", "Adrenal",
               ifelse(report_tissues == "Muscle_Skeletal", "Muscle", report_tissues))
report_cases <- report_cases[paste0(keys, "_TOD")]
saveRDS(report_cases, file.path(stage, "private", "report_pilots.rds"), compress = "xz")

# 6. Record the complete recipe and verify the new output tree ------------
stopifnot(length(list.files(pilot_root, pattern = "[.]rds$", recursive = TRUE)) == 307L,
          identical(native, manifest[manifest$dataset != "GTEx", , drop = FALSE]))
json_write(c(recipe, source_hashes, list(cache_mode = cache_mode, gtex_tissues = 46L,
  gtex_pilot_files = 92L, native_rows_preserved = 116L, native_pilot_files = 215L,
  n_report_cases = 10L, full_explorer_tables = TRUE, donor_times_in_public_pilots = FALSE,
  source_manifest_sha256 = file_hash("inst/extdata/pilots/manifest.csv"),
  manifest_sha256 = file_hash(file.path(pilot_root, "manifest.csv")),
  signal_summary_sha256 = file_hash(file.path(stage, "summary", "gtex_v10_signal_summary.csv")),
  histogram_sha256 = file_hash(file.path(pilot_root, "gtex_v10_TOD_sampling_histograms.csv")),
  report_pilots_sha256 = file_hash(file.path(stage, "private", "report_pilots.rds")),
  R_version = R.version.string)), file.path(stage, "qc", "build_recipe.json"))
cat("Normalized TOD pilots prepared:", normalizePath(stage), "\n")
