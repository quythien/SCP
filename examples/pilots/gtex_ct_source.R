# GTEx source preparation
# Run these sections from the package root before fitting a database pilot.

# 1. Load the local fitting code -------------------------------------------
gtex_source_bootstrap <- function(envir = parent.frame()) {
  suppressPackageStartupMessages(requireNamespace("SCP", quietly = TRUE))
  for (file in c("R/options.R", "R/estimation.R")) {
    if (!file.exists(file)) stop("Run this analysis from the SCP source root.")
    sys.source(file, envir = envir)
  }
  invisible(envir)
}

# 2. Read the phenotype table and resolve documented donor aliases --------
gtex_v8_phenotype <- function(
    path = Sys.getenv("SCP_GTEX_PHENO"),
    alias_path = Sys.getenv("SCP_GTEX_DONOR_ALIASES", "data/local/gtex_v8_donor_aliases.csv")) {
  if (!nzchar(path) || !file.exists(path)) {
    stop("Set SCP_GTEX_PHENO to the local gtex_TOD_pheno.rds file.")
  }
  if (!nzchar(alias_path) || !file.exists(alias_path)) {
    stop("Set SCP_GTEX_DONOR_ALIASES to the local donor-alias CSV.")
  }
  aliases <- read.csv(alias_path, stringsAsFactors = FALSE)
  if (!all(c("donor_original", "donor_resolved") %in% names(aliases))) {
    stop("The donor-alias CSV needs donor_original and donor_resolved columns.")
  }
  phenotype <- readRDS(path)
  attr(phenotype, "gtex_donor_aliases") <- aliases[c("donor_original", "donor_resolved")]
  attr(phenotype, "donor_aliases_path") <- normalizePath(alias_path, mustWork = TRUE)
  attr(phenotype, "donor_aliases") <- setNames(aliases$donor_resolved, aliases$donor_original)
  phenotype
}

gtex_v8_times <- function(sample_ids, phenotype) {
  sample_ids <- as.character(sample_ids)
  if (!length(sample_ids) || anyNA(sample_ids) || any(!nzchar(sample_ids)) ||
      anyDuplicated(sample_ids)) {
    stop("Sample identifiers must be present and unique.")
  }
  if (!is.data.frame(phenotype) ||
      !all(c("SUBJID", "CT") %in% names(phenotype))) {
    stop("Phenotype data must contain SUBJID and CT.")
  }
  donor_ids <- as.character(phenotype$SUBJID)
  if (anyNA(donor_ids) || any(!nzchar(donor_ids)) || anyDuplicated(donor_ids)) {
    stop("Phenotype donor identifiers must be present and unique.")
  }
  if (!is.numeric(phenotype$CT) ||
      any(!is.na(phenotype$CT) & !is.finite(phenotype$CT)) ||
      any(phenotype$CT < 0 | phenotype$CT >= 24, na.rm = TRUE)) {
    stop("Phenotype CT must be numeric hours in [0, 24), or NA.")
  }
  valid_id <- grepl("^GTEX[.-][A-Z0-9]+([.-]|$)", sample_ids)
  if (!all(valid_id)) stop("Unrecognized GTEx sample identifier.")
  donor_original <- sub("^GTEX[.-]([A-Z0-9]+).*", "GTEX-\\1", sample_ids)
  donor_resolved <- donor_original
  alias_table <- attr(phenotype, "gtex_donor_aliases", exact = TRUE)
  if (is.null(alias_table) || !is.data.frame(alias_table) ||
      !all(c("donor_original", "donor_resolved") %in% names(alias_table))) {
    stop("Load phenotype with gtex_v8_phenotype() and the declared local donor aliases.")
  }
  original <- as.character(alias_table$donor_original)
  resolved <- as.character(alias_table$donor_resolved)
  if (anyNA(original) || anyNA(resolved) || anyDuplicated(original) ||
      any(!grepl("^GTEX-[A-Z0-9]+$", original)) ||
      any(!grepl("^GTEX-[A-Z0-9]+$", resolved))) {
    stop("Donor aliases must contain unique original keys and valid donor identifiers.")
  }
  if (any(!resolved %in% donor_ids)) stop("A resolved donor alias is absent from phenotype.")
  aliases <- setNames(resolved, original)
  hit <- donor_original %in% names(aliases)
  donor_resolved[hit] <- unname(aliases[donor_original[hit]])
  index <- match(donor_resolved, donor_ids)
  ct <- unname(phenotype$CT[index])
  eligible <- is.finite(ct)
  status <- ifelse(is.na(index), "missing_donor", ifelse(eligible, "matched", "missing_CT"))
  data.frame(sample = sample_ids, donor_original = donor_original,
             donor_resolved = donor_resolved, seasonal_CT = ct,
             eligible = eligible, match_status = status, stringsAsFactors = FALSE)
}

# 3. Keep expression columns and CT aligned -------------------------------
gtex_v8_extract <- function(df, phenotype) {
  if (is.null(df) || is.null(dim(df)) || any(dim(df) == 0L)) return(NULL)
  if (is.null(colnames(df))) stop("Expression columns need GTEx sample identifiers.")
  metadata <- gtex_v8_times(colnames(df), phenotype)
  expr <- as.matrix(df[, metadata$eligible, drop = FALSE])
  if (!is.numeric(expr)) stop("Expression data must be numeric.")
  list(expr = expr, times = unname(metadata$seasonal_CT[metadata$eligible]),
       metadata = metadata, n_source = ncol(df), n_excluded = sum(!metadata$eligible))
}

# 4. Preserve the three paper-linked tissue pilots ------------------------
gtex_v8_protected <- function(tissue_or_file) {
  label <- tolower(gsub("[^a-zA-Z0-9]", "", basename(tissue_or_file)))
  label <- sub("^gtex", "", label)
  label <- sub("(all2hrds|allrds|2hrds|rds)$", "", label)
  label %in% c("adrenal", "adrenalgland", "liver", "muscleskeletal", "skeletalmuscle")
}

# 5. Send legacy GTEx entrypoints through the same staged analysis ---------
gtex_stage_database <- function() {
  analysis <- "examples/pilots/rebuild_gtex_database.R"
  if (!file.exists(analysis)) stop("Run this analysis from the SCP source root.")
  staging_env <- new.env(parent = globalenv())
  sys.source(analysis, envir = staging_env)
  invisible(staging_env$stage)
}
