# Check the published donor consensus and the paired pilot fits.
# These checks verify the saved computation. They do not establish true time.

# 1. Load the saved consensus and verify the author's calculation ----------
source("examples/pilots/gtex_v10/00_settings.R")
inference <- readRDS(file.path(stage, "private", "v10_DIP_inference.rds"))
reference <- upstream$Phi.from.phi_mat(inference$main_tips, ct = 1.95)
index <- match(inference$donors$donor, names(reference))
stopifnot(!anyNA(index))
error <- circular_difference(inference$donors$DIP_radians, unname(reference[index]))
stopifnot(max(error) < 1e-12)
json_write(list(n_compared = length(error), maximum_error_radians = max(error),
  function_used = "Pinned author Phi.from.phi_mat(ct=1.95) with min.tix=1",
  extra_stability_filters = FALSE), file.path(stage, "qc", "consensus_oracle.json"))
clock_qc <- inference$tissue_qc[, c("tissue", "n_clock_genes", "available_clock_genes", "missing_clock_genes")]
utils::write.csv(clock_qc, file.path(stage, "qc", "clock_gene_availability.csv"), row.names = FALSE)

# 2. Check shared cohorts, full-gene tests, and stored parameter tuples ------
broad <- read_expression(broad_path, "CPM.all.norm_large")
manifest <- utils::read.csv(file.path(stage, "inst", "extdata", "pilots", "manifest.csv"),
  stringsAsFactors = FALSE)
tissues <- unique(manifest$tissue[manifest$dataset == "GTEx"])
fit_checks <- list()
clock_checks <- list()
for (tissue in tissues) {
  cohort <- readRDS(file.path(stage, "private", "paired_cohorts", paste0(tissue, ".rds")))
  expr <- tissue_counts(broad[[cohort$source_tissue]])[cohort$genes, cohort$samples, drop = FALSE]
  stopifnot(ncol(expr) == length(cohort$DIP), ncol(expr) == length(cohort$TOD),
    !anyDuplicated(cohort$donors), all(is.finite(expr)))
  for (basis in c("DIP", "TOD")) {
    row <- manifest[manifest$dataset == "GTEx" & manifest$tissue == tissue &
      manifest$time_basis == basis, , drop = FALSE]
    stopifnot(nrow(row) == 1L)
    tables <- full_gene_tables(expr, cohort[[basis]])
    for (K in 1:2) {
      relative <- if (K == 1L) row$file else sub("[.]rds$", "_2H.rds", row$file)
      pilot <- readRDS(file.path(stage, "inst", "extdata", "pilots", relative))
      stopifnot(identical(pilot$explorer_fit, tables[[K]]),
        is.null(pilot$times), is.null(pilot$cts), is.null(pilot$raw),
        identical(pilot$time_basis, basis), identical(pilot$gtex_release, "v10"),
        inherits(pilot$time_distribution, "SCPTimeDistribution"),
        pilot$n_samples == ncol(expr), pilot$design_rank == 2L * K + 1L)
      fit_checks[[length(fit_checks) + 1L]] <- data.frame(tissue = tissue,
        time_basis = basis, K = K, n = ncol(expr), ngenes = nrow(expr),
        full_table_exact = TRUE, individual_time_fields_absent = TRUE)
    }
    independent <- full_gene_tables(expr[!toupper(gene_symbols(rownames(expr))) %in% clock_genes,
      , drop = FALSE], cohort[[basis]])[[1]]
    good <- which(is.finite(independent$p_K1) & independent$p_K1 < recipe$alpha_pilot &
      is.finite(independent$A) & independent$A > 0 & is.finite(independent$sigma) & independent$sigma > 0)
    selected <- head(good[order(independent$p_K1[good])], recipe$top_k)
    clock_checks[[length(clock_checks) + 1L]] <- data.frame(tissue = tissue,
      time_basis = basis, n = ncol(expr), ngenes = nrow(independent), n_selected = length(selected),
      medianr = median(independent$A[selected] / independent$sigma[selected]),
      prop_rhythmic = length(good) / nrow(independent), n_BH05 = sum(independent$q_K1 < 0.05),
      excluded_genes = "CHIRAL inference clock genes", stringsAsFactors = FALSE)
  }
}
utils::write.csv(do.call(rbind, fit_checks), file.path(stage, "qc", "pilot_fit_checks.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, clock_checks), file.path(stage, "qc", "clock_genes_excluded.csv"), row.names = FALSE)
cat("Published consensus and all paired full-gene fit tables verified.\n")

# 3. Describe the planned target-tissue omission check on matched donors -----
leave_rows <- list()
leave_private <- list()
for (source_tissue in report_tissues) {
  tissue <- if (source_tissue == "Adrenal_Gland") "Adrenal" else source_tissue
  if (!source_tissue %in% colnames(inference$main_tips)) {
    leave_private[[source_tissue]] <- list(status = "not_used_in_phase_inference")
    next
  }
  tips <- inference$main_tips[, colnames(inference$main_tips) != source_tissue, drop = FALSE]
  alternate <- upstream$Phi.from.phi_mat(tips, ct = 1.95)
  cohort <- readRDS(file.path(stage, "private", "paired_cohorts", paste0(tissue, ".rds")))
  alt_hours <- unname(alternate[cohort$donors]) * 12 / pi
  matched <- is.finite(alt_hours)
  stopifnot(sum(matched) > 5L)
  leave_private[[source_tissue]] <- list(status = "target_removed_before_author_consensus",
    donors = data.frame(donor = names(alternate), DIP = unname(alternate * 12 / pi)))
  expr <- tissue_counts(broad[[source_tissue]])[cohort$genes, cohort$samples[matched], drop = FALSE]
  expr <- expr[!toupper(gene_symbols(rownames(expr))) %in% clock_genes, , drop = FALSE]
  for (comparison in c("main_DIP", "leave_target_out_DIP")) {
    times <- if (comparison == "main_DIP") cohort$DIP[matched] else alt_hours[matched]
    fit <- full_gene_tables(expr, times)[[1]]
    good <- which(is.finite(fit$p_K1) & fit$p_K1 < recipe$alpha_pilot &
      is.finite(fit$A) & fit$A > 0 & is.finite(fit$sigma) & fit$sigma > 0)
    selected <- head(good[order(fit$p_K1[good])], recipe$top_k)
    leave_rows[[length(leave_rows) + 1L]] <- data.frame(tissue = tissue,
      status = "target_removed_before_author_consensus", comparison = comparison,
      n_original_paired = length(cohort$donors), n_lost = sum(!matched), n = sum(matched),
      ngenes = nrow(expr), n_selected = length(selected),
      medianr = median(fit$A[selected] / fit$sigma[selected]), prop_rhythmic = length(good) / nrow(expr),
      n_BH05 = sum(fit$q_K1 < 0.05))
  }
}
utils::write.csv(do.call(rbind, leave_rows), file.path(stage, "qc", "leave_target_out_summary.csv"), row.names = FALSE)
saveRDS(leave_private, file.path(stage, "private", "leave_target_out_DIP.rds"), compress = "xz")
cat("Target-tissue omission checks complete without changing the default consensus.\n")
