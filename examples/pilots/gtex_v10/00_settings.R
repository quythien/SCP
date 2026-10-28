# Inputs and helpers for the published CHIRAL workflow on GTEx v10.

# 1. Read local locations and load the pinned package ----------------------
required_input <- function(name) {
  value <- Sys.getenv(name)
  if (!nzchar(value) || !file.exists(value)) stop("Set ", name, " to an existing local input.")
  normalizePath(value, mustWork = TRUE)
}
stage <- Sys.getenv("SCP_GTEX_V10_STAGE", "output/gtex_v10_chiral_protocol")
strict_path <- Sys.getenv("SCP_GTEX_V10_STRICT", file.path(stage, "prepared_inputs", "v10_CPM.all.norm.chiral.RData"))
broad_path <- Sys.getenv("SCP_GTEX_V10_BROAD", file.path(stage, "prepared_inputs", "v10_CPM.all.norm_large.chiral.RData"))
phenotype_path <- Sys.getenv("SCP_GTEX_PHENO")
chiral_source <- required_input("SCP_CHIRAL_SOURCE")
chiral_library <- required_input("SCP_CHIRAL_LIBRARY")
.libPaths(c(chiral_library, .libPaths()))
stopifnot(requireNamespace("CHIRAL", quietly = TRUE),
          as.character(utils::packageVersion("CHIRAL")) == "0.1.0")
dir.create(file.path(stage, "private"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(stage, "qc"), recursive = TRUE, showWarnings = FALSE)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) RhpcBLASctl::blas_set_num_threads(1L)

# 2. Read the author's clock genes and fixed analysis settings --------------
chiral_revision <- "1f78a9a60b6d18facfedf5923d43bae9fa79f0b0"
clock_environment <- new.env()
load(file.path(chiral_source, "Paper", "paper_data", "CGRs.RData"), envir = clock_environment)
clock_genes <- get(ls(clock_environment)[1], envir = clock_environment)
report_tissues <- c("Adrenal_Gland", "Liver", "Muscle_Skeletal", "Thyroid", "Lung", "Pancreas",
                    "Heart_Left_Ventricle", "Artery_Tibial", "Adipose_Subcutaneous", "Brain_Cortex")
recipe <- list(release = "v10", period = 24, method = "Published CHIRAL protocol",
  chiral_revision = chiral_revision, chiral_package_version = "0.1.0",
  clock_genes = clock_genes, inference_seed = 20261028L, max_iterations = 500L,
  standardize = TRUE, GTEx_names = TRUE, TSM = TRUE, q = 0.1, update_q = FALSE,
  consensus_ct = 1.95, consensus_min_tissues = 1L,
  phase_origin = "Published human muscle reference with PER3 peak at 9 hours",
  pilot_seed = 2025L, alpha_pilot = 0.01, top_k = 300L, paired_sigma = TRUE)

# 3. Load the published preparation and consensus functions only ------------
upstream <- new.env(parent = globalenv())
upstream$CHIRAL <- CHIRAL::CHIRAL
upstream$mclapply <- parallel::mclapply
upstream$lrtest <- lmtest::lrtest
recipe_file <- file.path(chiral_source, "Paper", "DIP_inf.R")
functions <- c("CPM_to_E", "harm_reg", "Fit_OUT", "Set_OUT", "order.out.setgene.hc",
               "Unique_donors", "order.from.hc.ref", "Create_phi_matrix", "Phi.from.phi_mat",
               "Make_big_OUT", "spliti")
for (expression in parse(recipe_file)) {
  if (is.call(expression) && is.symbol(expression[[1]]) &&
      as.character(expression[[1]]) %in% c("<-", "=") && is.symbol(expression[[2]]) &&
      as.character(expression[[2]]) %in% functions) eval(expression, envir = upstream)
}
stopifnot(all(vapply(functions, exists, logical(1), envir = upstream, inherits = FALSE)))

# 4. Keep donor matching exact and read TOD only when it is requested -------
sample_donor <- function(ids) {
  ids <- as.character(ids)
  stopifnot(!anyNA(ids), all(grepl("^GTEX[.-][A-Z0-9]+([.-]|$)", ids)))
  sub("^GTEX[.-]([A-Z0-9]+).*", "GTEX-\\1", ids)
}
read_seasonal_time <- function(path) {
  phenotype <- readRDS(path)
  stopifnot(all(c("SUBJID", "CT") %in% names(phenotype)), !anyDuplicated(phenotype$SUBJID))
  times <- setNames(phenotype$CT, phenotype$SUBJID)
  stopifnot(all(is.na(times) | (is.finite(times) & times >= 0 & times < 24)))
  times
}

# 5. Small helpers preserve the supplied expression and identifiers ---------
read_expression <- function(path, object) {
  environment <- new.env()
  load(path, envir = environment)
  stopifnot(exists(object, environment, inherits = FALSE))
  get(object, environment, inherits = FALSE)
}
tissue_counts <- function(item) {
  expression <- if (is.list(item) && !is.null(item$counts)) item$counts else item
  expression <- as.matrix(expression)
  stopifnot(is.numeric(expression), !is.null(rownames(expression)), !is.null(colnames(expression)),
            !anyDuplicated(rownames(expression)), !anyDuplicated(colnames(expression)))
  expression
}
gene_symbols <- function(ids) sub("\\|.*$", "", sub("^.*_", "", ids))
gene_ensembl <- function(ids) sub("^(ENS[A-Z]*G[0-9]+).*", "\\1", ids)
circular_difference <- function(x, y) abs(Arg(exp(1i * (x - y))))
file_hash <- function(path) digest::digest(file = path, algo = "sha256")
json_write <- function(x, path) jsonlite::write_json(x, path, auto_unbox = TRUE,
  pretty = TRUE, null = "null", na = "null", digits = 16)

# 6. Package pilots store a sampling distribution rather than linked times -
sampling_distribution <- function(times) {
  times <- times %% 24
  den <- stats::density(c(times - 24, times, times + 24), from = 0, to = 24)
  structure(list(grid = unname(den$x), weights = unname(pmax(den$y, 1e-10)),
                 period = 24, jitter_sd = 0.1, n_observed = length(times)),
            class = "SCPTimeDistribution")
}

# 7. Keep full-universe model tables for the app and FDR summaries ----------
full_gene_tables <- function(expr, times) {
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
  first <- data.frame(gene = gene_ensembl(rownames(expr)),
                        symbol = gene_symbols(rownames(expr)), pvalue = p1,
                        qvalue = stats::p.adjust(p1, "BH", n = nrow(expr)),
                        p_K1 = p1, q_K1 = stats::p.adjust(p1, "BH", n = nrow(expr)),
                        A = sqrt(fit1$coefficients[2, ]^2 + fit1$coefficients[3, ]^2),
                        phi = (atan2(fit1$coefficients[3, ], fit1$coefficients[2, ]) * 12 / pi) %% 24,
                        sigma = sqrt(rss1 / (n - 3)), mesor = fit1$coefficients[1, ], n = n,
                        stringsAsFactors = FALSE, row.names = NULL)
  second <- data.frame(gene = first$gene, symbol = first$symbol,
                         pvalue = p2, qvalue = stats::p.adjust(p2, "BH", n = nrow(expr)),
                         p_joint = p2, q_joint = stats::p.adjust(p2, "BH", n = nrow(expr)),
                         p_K1 = p1, q_K1 = first$q_K1,
                         p_2h = ph, q_2h = stats::p.adjust(ph, "BH", n = nrow(expr)),
                         A = sqrt(fit2$coefficients[2, ]^2 + fit2$coefficients[3, ]^2),
                         phi = (atan2(fit2$coefficients[3, ], fit2$coefficients[2, ]) * 12 / pi) %% 24,
                         A2 = sqrt(fit2$coefficients[4, ]^2 + fit2$coefficients[5, ]^2),
                         phi2 = (atan2(fit2$coefficients[5, ], fit2$coefficients[4, ]) * 6 / pi) %% 12,
                         sigma = sqrt(rss2 / (n - 5)), mesor = fit2$coefficients[1, ], n = n,
                         stringsAsFactors = FALSE, row.names = NULL)
  list(first, second)
}
