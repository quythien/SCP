# Stage the full v10 analysis matrices with newly inferred DIP.
# This script creates a local candidate. It does not upload data.

# 1. Read expression and the newly inferred donor times --------------------
source("examples/pilots/gtex_v10/00_settings.R")
stopifnot(file.exists(broad_path))
broad <- read_expression(broad_path, "CPM.all.norm_large")
inference <- readRDS(file.path(stage, "private", "v10_DIP_inference.rds"))
public_sources <- jsonlite::read_json("examples/pilots/gtex_v10/public_inputs.json", simplifyVector = TRUE)
manifest <- utils::read.csv(file.path(stage, "inst", "extdata", "pilots", "manifest.csv"),
                            stringsAsFactors = FALSE)
targets <- unique(manifest$tissue[manifest$dataset == "GTEx"])
source_tissues <- ifelse(targets == "Adrenal", "Adrenal_Gland", targets)
stopifnot(length(source_tissues) == 46L, all(source_tissues %in% names(broad)),
          !any(grepl("^Cells", source_tissues)))
candidate <- file.path(stage, "full_data_candidate")
dir.create(file.path(candidate, "CPM"), recursive = TRUE, showWarnings = FALSE)

# 2. Create project keys and keep the private identity bridge separate ------
all_samples <- sort(unique(unlist(lapply(broad[source_tissues], function(x)
  colnames(tissue_counts(x))))))
all_donors <- sort(unique(sample_donor(all_samples)))
sample_key <- setNames(sprintf("SCPV10S%06d", seq_along(all_samples)), all_samples)
donor_key <- setNames(sprintf("SCPV10D%06d", seq_along(all_donors)), all_donors)
bridge <- data.frame(source_sample = all_samples, project_sample = unname(sample_key),
                      source_donor = sample_donor(all_samples),
                      project_donor = unname(donor_key[sample_donor(all_samples)]),
                      stringsAsFactors = FALSE)
saveRDS(bridge, file.path(stage, "private", "full_data_identity_bridge.rds"), compress = "xz")

# 3. Retain every source column and gene, including unavailable DIP ---------
coverage <- list()
sample_tables <- list()
for (i in seq_along(source_tissues)) {
  tissue <- source_tissues[i]
  expr <- tissue_counts(broad[[tissue]])
  donors <- sample_donor(colnames(expr))
  donor_index <- match(donors, inference$donors$donor)
  dip <- inference$donors$DIP[donor_index]
  status <- inference$donors$status[donor_index]
  status[is.na(donor_index)] <- "not_in_inference_cohort"
  pilot_tissue <- targets[i]
  cohort <- readRDS(file.path(stage, "private", "paired_cohorts", paste0(pilot_tissue, ".rds")))
  samples <- data.frame(sample_id = unname(sample_key[colnames(expr)]),
                         donor_id = unname(donor_key[donors]), tissue = tissue,
                         DIP = unname(dip), DIP_status = unname(status),
                         pilot_included = colnames(expr) %in% cohort$samples,
                         stringsAsFactors = FALSE)
  genes <- data.frame(gene_id = gene_ensembl(rownames(expr)),
                       symbol = gene_symbols(rownames(expr)), stringsAsFactors = FALSE)
  stopifnot(!anyDuplicated(genes$gene_id), !anyNA(samples$sample_id),
            !anyNA(samples$donor_id), !anyNA(samples$DIP_status),
            nrow(samples) == ncol(expr), nrow(genes) == nrow(expr))
  colnames(expr) <- samples$sample_id
  rownames(expr) <- genes$gene_id
  path <- file.path(candidate, "CPM", paste0("GTEx_v10_", tissue, ".rds"))
  saveRDS(list(expression = expr, gene_annotation = genes, sample_metadata = samples,
                gtex_release = "v10", input_time_basis = "DIP", time_units = "hours",
                expression_scale = "CHIRAL covariate-adjusted TMM log CPM residuals",
                data_license = "https://gtexportal.org/home/license"),
           path, compress = "gzip")
  roundtrip <- readRDS(path)
  stopifnot(identical(roundtrip$expression, expr),
            identical(roundtrip$sample_metadata, samples),
            identical(roundtrip$gene_annotation, genes))
  coverage[[tissue]] <- data.frame(tissue = tissue, n_genes = nrow(expr),
                                    n_samples = ncol(expr), n_with_DIP = sum(is.finite(dip)),
                                    n_without_DIP = sum(!is.finite(dip)),
                                    file = file.path("CPM", basename(path)),
                                    sha256 = file_hash(path), stringsAsFactors = FALSE)
  sample_tables[[tissue]] <- samples
  cat("Staged full matrix", tissue, nrow(expr), "genes", ncol(expr), "samples\n")
  flush.console()
}

# 4. Write the source index, minimal metadata, and concise instructions -----
utils::write.csv(do.call(rbind, coverage), file.path(candidate, "coverage.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, sample_tables), file.path(candidate, "sample_metadata.csv"),
                 row.names = FALSE, na = "")
json_write(list(dataset = "GTEx", release = "v10", input_time_basis = "DIP",
                expression_sha256 = file_hash(broad_path),
                inference_input_sha256 = file_hash(strict_path),
                chiral_revision = chiral_revision, chiral_package_version = recipe$chiral_package_version, recipe = recipe,
                preprocessing_recipe = "examples/pilots/gtex_v10/00_normalize_inputs.R",
                inference_recipe = "examples/pilots/gtex_v10/01_infer_DIP.R",
                public_subjects_SHA256 = "2313caeaaf3718ac8aead8c6b9eab92b3f613509cfe1a6f7b6286310a7200e73",
                source_permissions_status = "GTEx Portal open-access data, distributed under the GTEx Portal Data License",
                license_url = "https://gtexportal.org/home/license",
                acknowledgment = "The data used for the analyses described in this release were obtained from the GTEx Portal on 10/28/2026.",
                no_endorsement_implied = TRUE,
                study_accession = "phs000424",
                public_inputs = public_sources$inputs,
                build_context = "Expression and DIP use public inputs from 00/01. Pilot membership comes from the paired cohorts created by 02, which requires the local TOD table.",
                normalized_CPM_cache = public_sources$normalized_CPM_cache,
                public_metadata_fields = c("sample_id", "donor_id", "tissue", "DIP", "DIP_status", "pilot_included"),
                declared_universe = "All supplied genes and samples in the 46 selected in vivo analysis matrices",
                normalization = "Published CHIRAL preprocessing. TMM log CPM followed by adjustment for public age groups, sex, ischemic time groups, and death type.",
                excluded_fields = c("recorded TOD", "seasonal CT", "clinical metadata",
                                     "original sample identifiers", "original donor identifiers", "identity bridge"),
                DIP_citation = "Talamanca, Gobet and Naef (2023). Science 379:478-483. doi:10.1126/science.add0846"),
           file.path(candidate, "provenance.json"))
writeLines(c("# GTEx v10 analysis matrices with CHIRAL-estimated input time",
             "", "Each tissue file contains expression, gene annotation, and minimal sample metadata.",
             "Expression values are TMM log CPM residuals after the published CHIRAL covariate adjustment.",
             "The full retained analysis matrix is included before pilot gene selection.",
             "This collection does not contain raw sequencing or the complete official GTEx release.",
             "", "DIP is donor input time in hours on [0, 24).",
             "Missing DIP remains missing. DIP_status records the reason.",
             "pilot_included identifies the shared analysis cohort used for the paired pilots.",
             "Expression and DIP can be rebuilt from the public inputs through scripts 00 and 01.",
             "This archive builder reads paired-cohort membership from script 02, which requires the local TOD table.",
             "The phase origin uses the published human muscle reference with PER3 peaking at 9 hours.",
             "DIP was inferred from v10 expression. Recorded death time was not used to infer or orient DIP.",
             "", "Project identifiers are pseudonyms. They do not establish anonymity.",
             "Original identifiers and the private identity bridge are excluded.",
             "", "Method: Talamanca, Gobet and Naef (2023). Science 379:478-483.",
             "https://doi.org/10.1126/science.add0846",
             "", "Read provenance.json and coverage.csv for source hashes and the retained data universe.",
             "The data used for the analyses described in this release were obtained from the GTEx Portal on 10/28/2026.",
             "GTEx v10 data are distributed under the GTEx Portal Data License.",
             "https://gtexportal.org/home/license",
             "This release does not imply endorsement by GTEx or the CHIRAL authors."),
           file.path(candidate, "README.md"))
files <- sort(list.files(candidate, recursive = TRUE, full.names = TRUE))
files <- files[basename(files) != "SHA256SUMS"]
relative <- substring(files, nchar(candidate) + 2L)
checksums <- vapply(files, file_hash, character(1))
writeLines(paste(checksums, relative, sep = "  "), file.path(candidate, "SHA256SUMS"))

# 5. Package one archive without altering its gene or sample coverage -------
archive <- file.path(stage, "GTEx_v10_full_DIP.tar.xz")
old_xz <- Sys.getenv("XZ_OPT", unset = NA_character_)
Sys.setenv(XZ_OPT = "-1")
status <- system2("tar", c("-cJf", shQuote(archive), "-C", shQuote(candidate), "."))
if (is.na(old_xz)) Sys.unsetenv("XZ_OPT") else Sys.setenv(XZ_OPT = old_xz)
stopifnot(status == 0L, file.exists(archive))
members <- utils::untar(archive, list = TRUE)
stopifnot(length(members) >= length(files) + 1L)
bytes <- file.info(archive)$size
json_write(list(archive = basename(archive), bytes = bytes, sha256 = file_hash(archive),
                github_release_size_eligible = bytes < 2 * 1024^3,
                publication_status = "Staged for release validation",
                n_tissues = length(source_tissues), n_samples = length(all_samples),
                n_donors = length(all_donors)), file.path(stage, "qc", "full_archive_receipt.json"))
cat("Full data archive staged for release validation.\n")
