#' =======================================================================
#' Estimate Pilot Parameters from BA11 Younger Brain Data
#' =======================================================================
#'
#' One-time script that fits cosinor models to all genes in the BA11
#' younger brain microarray dataset and saves the empirical parameter
#' distributions as data/ba11_ba47_younger.rds.
#'
#' The resulting .rds file is used by CircadianBioOptions() as the
#' default parameter source (similar to PROPER's built-in cheung dataset).
#'
#' USAGE:
#'   cd PowerSim/
#'   Rscript data/estimate_pilot_params.R
#'
#' @author Thien Pham

# =====================================================================
# Setup
# =====================================================================

project_root <- "/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim"
setwd(project_root)

# Source the codebase (setup.R requires setwd into code/)
source_dir <- file.path(project_root, "code")
old_wd <- setwd(source_dir)
source("setup.R")
setwd(old_wd)

cat("====================================================================\n")
cat("ESTIMATING PILOT PARAMETERS FROM BA11 YOUNGER BRAIN DATA\n")
cat("====================================================================\n\n")

# =====================================================================
# Load BA11 younger brain expression data
# =====================================================================

cat("Loading BA11 younger brain expression data...\n")

COMBINED <- readRDS("/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Collaborative/Paper/Congruence/PNAS_aging/data/combined_data.rds")

expr_sample_names <- colnames(COMBINED$expr)
pheno_order <- match(expr_sample_names, COMBINED$pheno$sample_name)
valid_samples <- !is.na(pheno_order)
expr_sample_names <- expr_sample_names[valid_samples]
COMBINED$expr <- COMBINED$expr[, valid_samples]
pheno_order <- pheno_order[valid_samples]
pheno_data <- COMBINED$pheno[pheno_order, ]
pheno_data$tod <- if ("TOD.x" %in% colnames(pheno_data)) pheno_data$TOD.x else pheno_data$TOD.y
pheno_data$age_group_final <- if ("AgeGroup" %in% colnames(pheno_data)) pheno_data$AgeGroup else pheno_data$age_group
complete_samples <- !is.na(pheno_data$age_group_final) & !is.na(pheno_data$tod) &
  pheno_data$age_group_final %in% c("younger", "older")
pheno_clean <- pheno_data[complete_samples, ]
younger_idx <- pheno_clean$age_group_final == "younger"

expr_younger <- COMBINED$expr[, complete_samples][, younger_idx]
times_young <- pheno_clean$tod[younger_idx]

cat(sprintf("  Expression matrix: %d genes x %d samples\n", nrow(expr_younger), ncol(expr_younger)))
cat(sprintf("  Time points: %d unique values\n\n", length(unique(times_young))))

rm(COMBINED, expr_sample_names, pheno_order, valid_samples, pheno_data,
   complete_samples, pheno_clean, younger_idx)

# =====================================================================
# Estimate circadian parameters for all genes
# =====================================================================

cat("Fitting cosinor models to all genes...\n\n")

params <- estimate_circadian_params(expr_younger, times_young,
                                     period = 24,
                                     min_rhythm_pval = 0.1,
                                     verbose = TRUE)

# =====================================================================
# Build pilot dataset
# =====================================================================

raw <- params$raw
rhythmic <- raw$is_rhythmic & !is.na(raw$A) & raw$A > 0

# lBaselineExpr: log-mesor for all genes (drop NAs)
lBaselineExpr <- raw$M[!is.na(raw$M)]

# lOD: log-sigma for all genes (drop NAs and non-positive)
valid_sigma <- !is.na(raw$sigma) & raw$sigma > 0
lOD <- log(raw$sigma[valid_sigma])

# amplitude: rhythmic genes only
amplitude <- raw$A[rhythmic]

# phase: rhythmic genes only
phase <- raw$phi[rhythmic & !is.na(raw$phi)]

# proportion rhythmic
prop_rhythmic <- params$prop_rhythmic

pilot_data <- list(
  lBaselineExpr = lBaselineExpr,
  lOD           = lOD,
  amplitude     = amplitude,
  phase         = phase,
  prop_rhythmic = prop_rhythmic,
  source        = "BA11+BA47 younger, Human Brain Microarray",
  n_genes       = nrow(expr_younger),
  n_samples     = ncol(expr_younger)
)

cat("\n====================================================================\n")
cat("PILOT DATASET SUMMARY\n")
cat("====================================================================\n\n")
cat(sprintf("  Source:          %s\n", pilot_data$source))
cat(sprintf("  Genes:           %d\n", pilot_data$n_genes))
cat(sprintf("  Samples:         %d\n", pilot_data$n_samples))
cat(sprintf("  prop_rhythmic:   %.1f%%\n", 100 * pilot_data$prop_rhythmic))
cat(sprintf("  lBaselineExpr:   n=%d, mean=%.2f, sd=%.2f\n",
            length(pilot_data$lBaselineExpr),
            mean(pilot_data$lBaselineExpr),
            sd(pilot_data$lBaselineExpr)))
cat(sprintf("  lOD:             n=%d, mean=%.2f, sd=%.2f\n",
            length(pilot_data$lOD),
            mean(pilot_data$lOD),
            sd(pilot_data$lOD)))
cat(sprintf("  amplitude:       n=%d, mean=%.3f, median=%.3f\n",
            length(pilot_data$amplitude),
            mean(pilot_data$amplitude),
            median(pilot_data$amplitude)))
cat(sprintf("  phase:           n=%d\n", length(pilot_data$phase)))

# =====================================================================
# Save
# =====================================================================

out_file <- file.path(project_root, "data", "ba11_ba47_younger.rds")
saveRDS(pilot_data, file = out_file)
cat(sprintf("\nSaved: %s\n", out_file))
cat(sprintf("Size:  %.1f KB\n", file.size(out_file) / 1024))
cat("\nDone.\n")
