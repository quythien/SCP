# Three local v8 CPM alternatives with explicit CT metadata.
# Run from the SCP source root. These donor-linked files stay local.

# 1. Read the original source files --------------------------------------
source("examples/pilots/gtex_ct_source.R")
source_file <- Sys.getenv("SCP_GTEX_CPM_NORM")
phenotype_file <- Sys.getenv("SCP_GTEX_PHENO")
alias_file <- Sys.getenv("SCP_GTEX_DONOR_ALIASES", "data/local/gtex_v8_donor_aliases.csv")
phase_file <- Sys.getenv("SCP_GTEX_CHIRAL_PHASES")
destination <- Sys.getenv("SCP_GTEX_LOCAL_CPM", "data/local/gtex_v8_seasonal_CT")
stopifnot(file.exists(source_file), file.exists(phenotype_file), file.exists(phase_file),
          file.exists(alias_file))
expected_phase_hash <- "0aa5436fd3a461f37448e61decf9d274dbee64591de5302c88685e9ab0e5173d"
if (digest::digest(file = phase_file, algo = "sha256") != expected_phase_hash) {
  stop("Use the public CHIRAL DIPs.RData version identified in this analysis.")
}
expression_env <- new.env()
load(source_file, envir = expression_env)
phenotype <- gtex_v8_phenotype(phenotype_file)
phase_env <- new.env()
load(phase_file, envir = phase_env)
phase <- phase_env$phi
stopifnot(is.numeric(phase), length(phase) == 914L, !anyDuplicated(names(phase)),
          all(is.finite(phase)), all(phase >= 0 & phase < 2 * pi))

# 2. Align the CHIRAL clock using all matched phenotype donors ------------
phase_hours <- as.numeric(phase[match(sub("^GTEX-", "", phenotype$SUBJID),
                                    names(phase))]) * 24 / (2 * pi)
matched <- is.finite(phase_hours) & is.finite(phenotype$CT)
stopifnot(sum(matched) == 914L)
coherence <- vapply(c(1, -1), function(direction) {
  Mod(mean(exp(1i * 2 * pi *
                 (phenotype$CT[matched] - direction * phase_hours[matched]) / 24)))
}, numeric(1))
direction <- c(1, -1)[which.max(coherence)]
offset <- Arg(mean(exp(1i * 2 * pi *
                        (phenotype$CT[matched] - direction * phase_hours[matched]) / 24))) *
  24 / (2 * pi)

# 3. Describe the sources without modifying them -------------------------
provenance <- list(
  release = "v8", default_time_basis = "season_adjusted_recorded_CT",
  expression_sha256 = digest::digest(file = source_file, algo = "sha256"),
  phenotype_sha256 = digest::digest(file = phenotype_file, algo = "sha256"),
  donor_aliases_sha256 = digest::digest(file = alias_file, algo = "sha256"),
  chiral_sha256 = digest::digest(file = phase_file, algo = "sha256"),
  chiral_url = "https://github.com/naef-lab/CHIRAL/blob/master/Paper/paper_data/DIPs.RData",
  chiral_blob = "0c86fc3679d8cc300b7152580777d4caff39544a",
  chiral_alignment = list(n = sum(matched), direction = direction,
                          offset_hours = offset, coherence = max(coherence)),
  recipe = "examples/pilots/create_local_gtex_cpm.R")
dir.create(destination, recursive = TRUE, showWarnings = FALSE)

# 4. Save one alternative for each paper tissue ---------------------------
tissues <- c("Adrenal Gland", "Liver", "Muscle - Skeletal")
file_keys <- c("Adrenal_Gland", "Liver", "Muscle_Skeletal")
expected_n <- c(187L, 154L, 748L)
expected_eligible <- c(187L, 154L, 743L)
for (i in seq_along(tissues)) {
  tissue <- tissues[i]
  original <- expression_env$CPM.all.norm[[tissue]]
  sample_metadata <- gtex_v8_times(colnames(original), phenotype)
  chiral_index <- match(sub("^GTEX-", "", sample_metadata$donor_resolved), names(phase))
  sample_metadata$CHIRAL_radians <- unname(phase[chiral_index])
  sample_metadata$CHIRAL_hours <- sample_metadata$CHIRAL_radians * 24 / (2 * pi)
  sample_metadata$CHIRAL_aligned_hours <-
    (direction * sample_metadata$CHIRAL_hours + offset) %% 24
  CPM.all.norm <- setNames(list(original), tissue)
  CT.by.tissue <- setNames(list(setNames(sample_metadata$seasonal_CT,
                                        sample_metadata$sample)), tissue)
  provenance$tissue <- tissue
  provenance$expression_object_sha256 <- digest::digest(original, algo = "sha256")
  stopifnot(identical(CPM.all.norm[[tissue]], original),
            identical(names(CT.by.tissue[[tissue]]), colnames(original)),
            ncol(original) == expected_n[i],
            sum(sample_metadata$eligible) == expected_eligible[i])
  output <- file.path(destination, paste0("GTEx_v8_", file_keys[i], "_seasonal_CT.RData"))
  if (file.exists(output)) stop("Alternative already exists. Preserve it before rebuilding: ", output)
  save(CPM.all.norm, CT.by.tissue, sample_metadata, provenance, file = output, compress = "xz")
  cat("Saved", basename(output), "with", sum(sample_metadata$eligible), "verified CT values\n")
}
