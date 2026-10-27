# Rebuild the GTEx database pilots from local source data.
# Run from the SCP source root. Outputs go to a separate staging directory.
# The Adrenal, Liver and Muscle paper inputs are never overwritten here.

# 1. Choose the local inputs and output directory -------------------------
source("examples/pilots/gtex_ct_source.R")
source("examples/pilots/update_gtex_signal_csvs.R")
gtex_source_bootstrap()

source_file <- Sys.getenv("SCP_GTEX_CPM_NORM")
phenotype_file <- Sys.getenv("SCP_GTEX_PHENO")
alias_file <- Sys.getenv("SCP_GTEX_DONOR_ALIASES", "data/local/gtex_v8_donor_aliases.csv")
stage <- Sys.getenv("SCP_GTEX_STAGE", "output/pilot_rebuild/gtex_v8")
stopifnot(nzchar(source_file), file.exists(source_file),
          nzchar(phenotype_file), file.exists(phenotype_file), file.exists(alias_file))
dir.create(file.path(stage, "inst/extdata/pilots/human"), recursive = TRUE,
           showWarnings = FALSE)
source_hash <- digest::digest(file = source_file, algo = "sha256")
phenotype_hash <- digest::digest(file = phenotype_file, algo = "sha256")
alias_hash <- digest::digest(file = alias_file, algo = "sha256")

# 2. Load the original expression matrices and the donor CT table ---------
source_env <- new.env()
load(source_file, envir = source_env)
stopifnot(exists("CPM.all.norm", envir = source_env, inherits = FALSE))
cpm <- source_env$CPM.all.norm
phenotype <- gtex_v8_phenotype(phenotype_file)
manifest <- read.csv("inst/extdata/pilots/manifest.csv", stringsAsFactors = FALSE)

# 3. Keep the original gene identifiers and symbols aligned --------------
pilot_name <- function(tissue) {
  key <- gsub(" - ", "_", tissue, fixed = TRUE)
  key <- gsub("[()]", "", key)
  key <- gsub(" ", "_", key, fixed = TRUE)
  if (key == "Adrenal_Gland") key <- "Adrenal"
  paste0("human/GTEx_", key, "_All.rds")
}

gene_labels <- function(ids) {
  gene <- sub("^(ENS[A-Z]*G[0-9]+).*", "\\1", ids)
  symbol <- sub("^ENS[A-Z]*G[0-9.]+_", "", ids)
  symbol[symbol == ids] <- gene[symbol == ids]
  stopifnot(!anyDuplicated(gene), !anyNA(gene))
  data.frame(gene = gene, symbol = symbol, stringsAsFactors = FALSE)
}

add_gene_labels <- function(pilot, labels) {
  gene <- sub("^(ENS[A-Z]*G[0-9]+).*", "\\1", pilot$rhythm_fit$gene)
  index <- match(gene, labels$gene)
  stopifnot(!anyNA(index))
  pilot$rhythm_fit$gene <- gene
  pilot$rhythm_fit$symbol <- labels$symbol[index]
  pilot$tested_genes <- labels$gene
  pilot
}

# 4. Fit an independent full gene table for spreadsheet summaries --------
# These fits retain every tested gene so BH counts use the full universe.
full_gene_fit <- function(expr, times) {
  n <- length(times)
  design <- cbind(1, cos(2 * pi * times / 24), sin(2 * pi * times / 24))
  fit <- lm.fit(design, t(expr))
  stopifnot(fit$rank == 3L)
  rss <- colSums(fit$residuals^2)
  tss <- rowSums((expr - rowMeans(expr))^2)
  data.frame(gene = rownames(expr),
             pvalue = pf(((tss - rss) / 2) / (rss / (n - 3)),
                         2, n - 3, lower.tail = FALSE),
             A = sqrt(fit$coefficients[2, ]^2 + fit$coefficients[3, ]^2),
             phi = (atan2(fit$coefficients[3, ], fit$coefficients[2, ]) *
                      24 / (2 * pi)) %% 24,
             sigma = sqrt(rss / (n - 3)), mesor = fit$coefficients[1, ],
             row.names = NULL, stringsAsFactors = FALSE)
}

# 5. Rebuild the 31 eligible v8 tissues for both harmonic orders ----------
updates <- list()
build_rows <- list()
for (tissue in names(cpm)) {
  original <- cpm[[tissue]]
  if (is.null(dim(original)) || any(dim(original) == 0L)) next
  if (gtex_v8_protected(tissue)) next
  file <- pilot_name(tissue)
  manifest_row <- which(manifest$file == file & manifest$status == "ingested")
  stopifnot(length(manifest_row) == 1L)
  input <- gtex_v8_extract(original, phenotype)
  expr <- input$expr
  times <- input$times
  stopifnot(ncol(expr) >= 6L, all(is.finite(expr)), all(is.finite(times)))
  labels <- gene_labels(rownames(expr))
  full_fit <- full_gene_fit(expr, times)

  for (harmonic in 1:2) {
    relative <- if (harmonic == 1L) file else sub("[.]rds$", "_2H.rds", file)
    previous <- readRDS(file.path("inst/extdata/pilots", relative))
    top_k <- previous$pilot_top_k
    stopifnot(length(top_k) == 1L, top_k > 0)
    set.seed(2025)
    if (harmonic == 1L) {
      pilot <- estCircadianParam(expr, times, period = 24,
                                min_rhythm_pval = 0.01, top_k = top_k,
                                paired_sigma = TRUE, verbose = FALSE)
    } else {
      pilot <- estCircadianParam2H(expr, times, period = 24,
                                  min_rhythm_pval = 0.01, top_k = top_k,
                                  paired_sigma = TRUE, verbose = FALSE)
      pilot$diagnostics$screen_method <- "K1_F2"
    }
    pilot <- add_gene_labels(pilot, labels)
    pilot$times <- unname(times)
    pilot$cts <- unname(times)
    pilot$source_provenance <- list(
      dataset = "GTEx", release = "v8", time_basis = "season_adjusted_recorded_CT",
      expression_sha256 = source_hash, phenotype_sha256 = phenotype_hash,
      donor_aliases_sha256 = alias_hash,
      n_source = input$n_source, n_used = length(times), n_excluded = input$n_excluded,
      seed = 2025L, period = 24, harmonic_order = harmonic,
      alpha_pilot = 0.01, top_k = top_k, paired_sigma = TRUE,
      recipe = "examples/pilots/rebuild_gtex_database.R")
    saveRDS(pilot, file.path(stage, "inst/extdata/pilots", relative), compress = "xz")
    build_rows[[length(build_rows) + 1L]] <- data.frame(
      tissue = tissue, file = relative, K = harmonic, n = length(times),
      n_excluded = input$n_excluded, ngenes = nrow(expr), top_k = top_k,
      stringsAsFactors = FALSE)
    if (harmonic == 1L) {
      updates[[length(updates) + 1L]] <- list(
        tissue = tissue, pilot_file = file, times = times, ngenes = nrow(expr),
        full_fit = full_fit, prop_rhythmic = pilot$prop_rhythmic)
    }
  }
  cat("Staged", tissue, "with", length(times), "samples\n")
}
stopifnot(length(updates) == 31L, length(build_rows) == 62L)

# 6. Update only the corresponding rows in the three signal CSVs ---------
csv_report <- update_gtex_signal_csvs(updates, stage, package_root = ".")
saveRDS(updates, file.path(stage, "full_gene_summary_inputs.rds"), compress = "xz")
write.csv(do.call(rbind, build_rows), file.path(stage, "build_inventory.csv"),
          row.names = FALSE)

# 7. Keep manifest source information explicit ---------------------------
if (!"gtex_release" %in% names(manifest)) manifest$gtex_release <- NA_character_
if (!"time_basis" %in% names(manifest)) manifest$time_basis <- NA_character_
v8_files <- vapply(names(cpm)[vapply(cpm, function(x) {
  !is.null(dim(x)) && all(dim(x) > 0L)
}, logical(1))], pilot_name, character(1))
v8_rows <- which(manifest$dataset == "GTEx" & manifest$file %in% v8_files)
manifest$gtex_release[v8_rows] <- "v8"
manifest$time_basis[v8_rows] <- "sample_barcode_field"
v10_tissues <- c("Cells_Cultured_fibroblasts", "Whole_Blood",
                 "Brain_Frontal_Cortex_BA9", "Brain_Caudate_basal_ganglia",
                 "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Putamen_basal_ganglia",
                 "Brain_Hypothalamus", "Brain_Hippocampus",
                 "Brain_Anterior_cingulate_cortex_BA24", "Brain_Cerebellar_Hemisphere",
                 "Brain_Substantia_nigra", "Kidney_Cortex", "Brain_Amygdala")
v10_rows <- which(manifest$dataset == "GTEx" & manifest$tissue %in% v10_tissues &
                   manifest$status == "ingested")
stopifnot(length(v10_rows) == 13L)
manifest$gtex_release[v10_rows] <- "v10"
manifest$time_basis[v10_rows] <- "season_adjusted_recorded_CT"

for (update in updates) {
  i <- match(update$pilot_file, manifest$file)
  manifest$n[i] <- length(update$times)
  manifest$ngenes[i] <- update$ngenes
  manifest$time_basis[i] <- "season_adjusted_recorded_CT"
  td <- gtex_signal_tod_summary(update$times)
  for (column in intersect(names(td), names(manifest))) {
    manifest[i, column] <- td[[column]]
  }
}
stopifnot(!anyDuplicated(manifest[c("species", "dataset", "tissue", "condition")]))
write.csv(manifest, file.path(stage, "inst/extdata/pilots/manifest.csv"), row.names = FALSE)
cat("Staging complete. Paper inputs and canonical package files were not replaced.\n")
