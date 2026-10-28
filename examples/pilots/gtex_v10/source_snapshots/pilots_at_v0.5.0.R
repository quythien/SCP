# Build paired GTEx v10 pilots for CHIRAL DIP and seasonal TOD.
# This script writes a separate staging tree and leaves paper inputs alone.

# 1. Load the fixed recipe and the fitted donor input times ----------------
source("examples/pilots/gtex_v10/00_settings.R")
stopifnot(file.exists(broad_path), file.exists(strict_path))
phenotype_path <- required_input("SCP_GTEX_PHENO")
seasonal_time <- read_seasonal_time(phenotype_path)
for (file in c("R/options.R", "R/estimation.R")) sys.source(file, envir = globalenv())
broad <- read_expression(broad_path, "CPM.all.norm_large")
inference <- readRDS(file.path(stage, "private", "v10_DIP_inference.rds"))
dip <- setNames(inference$donors$DIP, inference$donors$donor)
original_manifest <- utils::read.csv("inst/extdata/pilots/manifest.csv",
                                     stringsAsFactors = FALSE)
targets <- original_manifest[original_manifest$dataset == "GTEx" &
                               original_manifest$status == "ingested" &
                               !grepl("^Cells", original_manifest$tissue) &
                               original_manifest$time_basis == "DIP", , drop = FALSE]
stopifnot(nrow(targets) == 46L)
targets$source_tissue <- ifelse(targets$tissue == "Adrenal", "Adrenal_Gland", targets$tissue)
stopifnot(all(targets$source_tissue %in% names(broad)))
# Fit the report tissues first so simulations can run during the other fits.
priority <- match(targets$source_tissue, report_tissues)
targets <- targets[order(is.na(priority), priority, seq_len(nrow(targets))), , drop = FALSE]
pilot_root <- file.path(stage, "inst", "extdata", "pilots")
dir.create(pilot_root, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(stage, "summary"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(stage, "private", "paired_cohorts"), recursive = TRUE, showWarnings = FALSE)
public_sources <- jsonlite::read_json("examples/pilots/gtex_v10/public_inputs.json", simplifyVector = TRUE)
source_hashes <- list(public_gene_reads_sha256 = public_sources$inputs$gene_reads$sha256,
                      public_sample_QC_sha256 = public_sources$inputs$sample_QC$sha256,
                      public_subjects_sha256 = public_sources$inputs$public_subjects$sha256,
                      author_To_CPM_sha256 = file_hash(file.path(chiral_source, "Paper", "To_CPM.R")),
                      author_DIP_inf_sha256 = file_hash(recipe_file),
                      preprocessing_recipe_sha256 = file_hash("examples/pilots/gtex_v10/00_normalize_inputs.R"),
                      inference_recipe_sha256 = file_hash("examples/pilots/gtex_v10/01_infer_DIP.R"),
                      pilot_recipe_sha256 = file_hash("examples/pilots/gtex_v10/02_build_pilots.R"),
                      expression_sha256 = file_hash(broad_path),
                      strict_expression_sha256 = file_hash(strict_path),
                      seasonal_phenotype_sha256 = file_hash(phenotype_path),
                      inference_sha256 = file_hash(file.path(stage, "private", "v10_DIP_inference.rds")))

# 2. Fit all genes independently for the signal summaries ------------------
fit_full_genes <- function(expr, times) {
  n <- length(times)
  design <- cbind(1, cos(2 * pi * times / 24), sin(2 * pi * times / 24))
  fit <- stats::lm.fit(design, t(expr))
  stopifnot(n > 5L, fit$rank == 3L)
  rss <- colSums(fit$residuals^2)
  centered <- expr - rowMeans(expr)
  tss <- rowSums(centered^2)
  data.frame(gene = gene_ensembl(rownames(expr)), symbol = gene_symbols(rownames(expr)),
             pvalue = stats::pf((pmax(tss - rss, 0) / 2) / (rss / (n - 3)),
                                2, n - 3, lower.tail = FALSE),
             A = sqrt(fit$coefficients[2, ]^2 + fit$coefficients[3, ]^2),
             phi = (atan2(fit$coefficients[3, ], fit$coefficients[2, ]) * 12 / pi) %% 24,
             sigma = sqrt(rss / (n - 3)), mesor = fit$coefficients[1, ],
             stringsAsFactors = FALSE, row.names = NULL)
}

signal_summary <- function(fit, tissue, basis, n, n_source, n_dip, n_tod) {
  candidates <- which(is.finite(fit$pvalue) & fit$pvalue < recipe$alpha_pilot &
                        is.finite(fit$A) & fit$A > 0 &
                        is.finite(fit$sigma) & fit$sigma > 0)
  selected <- head(candidates[order(fit$pvalue[candidates])], recipe$top_k)
  ratio <- fit$A[selected] / fit$sigma[selected]
  quant <- function(p) if (length(ratio)) unname(stats::quantile(ratio, p)) else NA_real_
  data.frame(species = "human", dataset = "GTEx", tissue = tissue,
             condition = "All", gtex_release = "v10", time_basis = basis,
             is_default = basis == "DIP", n = n, n_source = n_source,
             n_with_DIP = n_dip, n_with_TOD = n_tod, ngenes = nrow(fit),
             n_selected = length(selected), r_median_top300 = quant(0.5),
             r_q25_top300 = quant(0.25), r_q75_top300 = quant(0.75),
             prop_rhythmic = length(candidates) / nrow(fit),
             n_BH05 = sum(stats::p.adjust(fit$pvalue, "BH", n = nrow(fit)) < 0.05, na.rm = TRUE),
             stringsAsFactors = FALSE)
}

label_pilot_genes <- function(pilot, expr) {
  labels <- data.frame(original = rownames(expr), gene = gene_ensembl(rownames(expr)),
                       symbol = gene_symbols(rownames(expr)), stringsAsFactors = FALSE)
  stopifnot(!anyDuplicated(labels$gene))
  index <- match(pilot$rhythm_fit$gene, labels$original)
  stopifnot(!anyNA(index))
  pilot$rhythm_fit$gene <- labels$gene[index]
  pilot$rhythm_fit$symbol <- labels$symbol[index]
  pilot$tested_genes <- labels$gene
  pilot
}

# 3. Build each tissue from the same donors and genes for both time bases ---
manifest_rows <- list()
summary_rows <- list()
private_cases <- list()
build_rows <- list()
for (i in seq_len(nrow(targets))) {
  row <- targets[i, , drop = FALSE]
  tissue <- row$tissue
  source_tissue <- row$source_tissue
  expr_all <- tissue_counts(broad[[source_tissue]])
  donors <- sample_donor(colnames(expr_all))
  stopifnot(!anyDuplicated(donors))
  metadata <- broad[[source_tissue]]$metadata
  sample_index <- match(colnames(expr_all), metadata$SAMPID)
  stopifnot(!anyNA(sample_index),
            identical(as.character(metadata$SUBJID[sample_index]), donors))
  dip_time <- unname(dip[donors])
  tod_time <- unname(seasonal_time[donors])
  matched <- is.finite(dip_time) & is.finite(tod_time)
  expr <- expr_all[, matched, drop = FALSE]
  keep_gene <- apply(expr, 1, function(x) all(is.finite(x)) && stats::sd(x) > 0)
  expr <- expr[keep_gene, , drop = FALSE]
  stopifnot(ncol(expr) > 5L, nrow(expr) > 0L, all(is.finite(expr)))
  private_cohort <- list(tissue = tissue, source_tissue = source_tissue,
                         samples = colnames(expr), donors = donors[matched],
                         DIP = dip_time[matched], TOD = tod_time[matched],
                         genes = rownames(expr), n_source = ncol(expr_all),
                         gene_qc_excluded = sum(!keep_gene))
  saveRDS(private_cohort, file.path(stage, "private", "paired_cohorts", paste0(tissue, ".rds")),
          compress = "xz")

  for (basis in c("DIP", "TOD")) {
    times <- private_cohort[[basis]]
    full_tables <- full_gene_tables(expr, times)
    full_fit <- full_tables[[1]]
    summary <- signal_summary(full_fit, tissue, basis, length(times), ncol(expr_all),
                               sum(is.finite(dip_time)), sum(is.finite(tod_time)))
    if (summary$n_selected == 0L)
      stop("No rhythmic candidates for ", tissue, " ", basis, ". Do not substitute a default pilot.")
    summary_rows[[length(summary_rows) + 1L]] <- summary
    relative <- file.path("human", "gtex_v10", tolower(basis), basename(row$file))
    source_label <- if (basis == "DIP") "CHIRAL-estimated time (DIP)" else
      "Season-adjusted time of death (TOD)"
    row_out <- row[, setdiff(names(row), "source_tissue"), drop = FALSE]
    row_out$n <- length(times)
    row_out$ngenes <- nrow(expr)
    row_out$file <- relative
    row_out$gtex_release <- "v10"
    row_out$accession <- "phs000424"
    row_out$license <- "GTEx Portal Data License"
    row_out$time_basis <- basis
    row_out$time_source <- source_label
    row_out$citation <- paste("GTEx Consortium. GTEx v10. phs000424.",
      "CHIRAL: Talamanca, Gobet and Naef (2023). Science 379:478-483. doi:10.1126/science.add0846")
    row_out$is_default <- basis == "DIP"
    row_out$n_source <- ncol(expr_all)
    row_out$n_excluded <- ncol(expr_all) - ncol(expr)
    row_out$tod_n_unique <- length(unique(times))
    row_out$tod_cycles <- NA_real_
    row_out$tod_step_hr <- NA_character_
    row_out$tod_phases <- NA_character_
    row_out$note <- if (source_tissue == "Whole_Blood")
      "Donor DIP transferred from other organs. Whole Blood was not used for inference." else ""
    row_out$has_k2 <- TRUE

    for (harmonic in 1:2) {
      set.seed(recipe$pilot_seed)
      pilot <- suppressWarnings(if (harmonic == 1L)
        estCircadianParam(expr, times, period = 24, min_rhythm_pval = recipe$alpha_pilot,
                           top_k = recipe$top_k, paired_sigma = TRUE, verbose = FALSE) else
        estCircadianParam2H(expr, times, period = 24, min_rhythm_pval = recipe$alpha_pilot,
                             top_k = recipe$top_k, paired_sigma = TRUE, verbose = FALSE))
      pilot <- label_pilot_genes(pilot, expr)
      if (harmonic == 1L) {
        fit <- pilot$rhythm_fit
        selected <- which(is.finite(fit$pvalue) & fit$pvalue < recipe$alpha_pilot &
                            is.finite(fit$A) & fit$A > 0 &
                            is.finite(fit$sigma) & fit$sigma > 0 & is.finite(fit$phi))
        selected <- head(selected[order(fit$pvalue[selected])], recipe$top_k)
        pilot$amplitude <- pilot$amplitude_spec <- unname(fit$A[selected])
        pilot$sigma_rhythmic <- unname(fit$sigma[selected])
        pilot$phase <- pilot$phase_spec <- unname(fit$phi[selected])
        pilot$paired_sigma <- TRUE
      } else stopifnot(length(pilot$amplitude) == length(pilot$sigma_rhythmic),
                         length(pilot$amplitude) == length(pilot$phase),
                         length(pilot$amplitude) == length(pilot$amplitude2),
                         length(pilot$amplitude) == length(pilot$phase2))
      pilot$amplitude_spec <- pilot$amplitude
      pilot$phase_spec <- pilot$phase
      pilot$explorer_fit <- full_tables[[harmonic]]
      if (harmonic == 2L) pilot$diagnostics$screen_method <- "K1_F2"
      pilot$times <- NULL
      pilot$cts <- NULL
      pilot$raw <- NULL
      pilot$time_distribution <- sampling_distribution(times)
      pilot$gtex_release <- "v10"
      pilot$time_basis <- basis
      pilot$time_source <- source_label
      pilot$n_samples <- length(times)
      pilot$n_distinct_times <- length(unique(times))
      x <- cbind(1, cos(2 * pi * times / 24), sin(2 * pi * times / 24))
      if (harmonic == 2L) x <- cbind(x, cos(4 * pi * times / 24), sin(4 * pi * times / 24))
      pilot$design_rank <- qr(x)$rank
      stopifnot(pilot$design_rank == 2L * harmonic + 1L)
      pilot$source_provenance <- c(source_hashes, list(dataset = "GTEx", release = "v10",
        time_basis = basis, time_source = source_label, chiral_revision = chiral_revision,
        accession = "phs000424", data_license = "GTEx Portal Data License",
        license_url = "https://gtexportal.org/home/license",
        acknowledgment = "The data used for the analyses described in this release were obtained from the GTEx Portal on 10/28/2026.",
        no_endorsement_implied = TRUE,
        chiral_package_version = recipe$chiral_package_version,
        preprocessing = "Published CHIRAL TMM log CPM and covariate adjustment",
        covariates = "Public age groups, sex, ischemic time groups, and death type",
        preprocessing_recipe = "examples/pilots/gtex_v10/00_normalize_inputs.R",
        inference_recipe = "examples/pilots/gtex_v10/01_infer_DIP.R",
        n_source = ncol(expr_all), n_used = ncol(expr), n_excluded = sum(!matched),
        gene_qc_excluded = sum(!keep_gene), n_genes = nrow(expr),
        seed = recipe$pilot_seed, period = 24, harmonic_order = harmonic,
        alpha_pilot = recipe$alpha_pilot, top_k = recipe$top_k, paired_sigma = TRUE,
        donor_matching = "exact v10 source identifiers", cohort = "paired DIP and TOD",
        recipe = "examples/pilots/gtex_v10/02_build_pilots.R"))
      if (harmonic == 2L) pilot$source_provenance$pilot_screen <- "K1_F2"
      path <- if (harmonic == 1L) relative else sub("[.]rds$", "_2H.rds", relative)
      destination <- file.path(pilot_root, path)
      dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
      saveRDS(pilot, destination, compress = "xz")
      build_rows[[length(build_rows) + 1L]] <- data.frame(
        tissue = tissue, time_basis = basis, K = harmonic, file = path,
        n = length(times), ngenes = nrow(expr), output_sha256 = file_hash(destination),
        stringsAsFactors = FALSE)
      if (harmonic == 1L && source_tissue %in% report_tissues) {
        report_key <- if (source_tissue == "Adrenal_Gland") "Adrenal" else
          if (source_tissue == "Muscle_Skeletal") "Muscle" else source_tissue
        report_pilot <- pilot
        report_pilot$explorer_fit <- NULL
        fit <- report_pilot$rhythm_fit
        selected <- which(is.finite(fit$pvalue) & fit$pvalue < recipe$alpha_pilot &
                            is.finite(fit$A) & fit$A > 0 &
                            is.finite(fit$sigma) & fit$sigma > 0 & is.finite(fit$phi))
        selected <- head(selected[order(fit$pvalue[selected])], recipe$top_k)
        report_pilot$amplitude <- report_pilot$amplitude_spec <- unname(fit$A[selected])
        report_pilot$sigma_rhythmic <- unname(fit$sigma[selected])
        report_pilot$phase <- report_pilot$phase_spec <- unname(fit$phi[selected])
        report_pilot$paired_sigma <- TRUE
        report_pilot$ngenes <- 5000L
        private_cases[[paste(report_key, basis, sep = "_")]] <- list(
          bio = report_pilot, times = unname(times), tissue = report_key,
          time_basis = basis, cohort = "matched", n_pilot = length(times),
          genes_fitted = nrow(expr), n_selected = summary$n_selected,
          medianr = summary$r_median_top300, prop = pilot$prop_rhythmic,
          calibration_law = "installed_loader_top300_v1")
      }
    }
    manifest_rows[[length(manifest_rows) + 1L]] <- row_out
  }
  if (length(private_cases) == 2L * length(report_tissues)) {
    early_report_path <- file.path(stage, "private", "report_pilots.rds")
    if (!file.exists(early_report_path)) {
      saveRDS(private_cases, early_report_path, compress = "xz")
      cat("All 20 matched report cases are ready for independent simulations.\n")
      flush.console()
    }
  }
  cat("Staged", tissue, "paired samples", ncol(expr), "genes", nrow(expr), "\n")
  flush.console()
}

# 4. Retain non-GTEx rows and publish an explicit current v10 summary -------
native <- original_manifest[original_manifest$dataset != "GTEx", , drop = FALSE]
manifest <- rbind(native, do.call(rbind, manifest_rows))
stopifnot(!anyDuplicated(manifest[c("species", "dataset", "tissue", "condition", "time_basis")] ),
          sum(manifest$dataset == "GTEx") == 92L)
utils::write.csv(manifest, file.path(pilot_root, "manifest.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, summary_rows),
                 file.path(stage, "summary", "gtex_v10_signal_summary.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, build_rows), file.path(stage, "qc", "pilot_build_inventory.csv"),
                 row.names = FALSE)
report_path <- file.path(stage, "private", "report_pilots.rds")
if (file.exists(report_path)) {
  previous_cases <- readRDS(report_path)
  if (!identical(previous_cases, private_cases)) {
    revision <- Sys.getenv("SCP_GTEX_REPORT_CASE_REVISION")
    expected_old_hash <- Sys.getenv("SCP_GTEX_REPORT_CASE_PREVIOUS_SHA256")
    if (!identical(revision, "installed_loader_top300_v1") ||
        !nzchar(expected_old_hash) || !identical(file_hash(report_path), expected_old_hash))
      stop("Report cases differ. Preserve their prior archive and set the explicit revision and previous hash.")
    saveRDS(private_cases, report_path, compress = "xz")
  }
} else saveRDS(private_cases, report_path, compress = "xz")
json_write(c(recipe, source_hashes, list(n_tissues = nrow(targets),
                                        n_pilots = length(build_rows),
                                        source_normalization = "Published CHIRAL preprocessing using public v10 age groups, sex, ischemic time groups, and death type")),
           file.path(stage, "qc", "pilot_recipe.json"))
cat("Staging complete:", length(build_rows), "pilot files and",
    length(private_cases), "private report cases.\n")
