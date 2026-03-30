#' =======================================================================
#' DP Only Power Analysis (n up to 160)
#' =======================================================================
#'
#' USAGE:
#'   Rscript examples/run_dp_only_panelC160.R

set.seed(12345)

setwd("/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")
source_dir <- file.path(getwd(), "code")
old_wd <- setwd(source_dir)
source("setup.R")
setwd(old_wd)

cat("Loading pilot expression data (PFC younger: BA11 + BA47)...\n")
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
cat(sprintf("  Reference TOD: n=%d younger subjects\n\n", length(times_young)))
rm(COMBINED, expr_sample_names, pheno_order, valid_samples, pheno_data,
   complete_samples, pheno_clean, younger_idx)

cat("Estimating circadian parameters from pilot data...\n\n")
opts_bio <- estCircadianParam(
  data          = expr_younger,
  times         = times_young,
  period        = 24,
  prop_DR       = 0.15,
  prop_DP       = 0.10,
  prop_DA       = 0.10,
  phase_diff    = c(-6, 6),
  amp_diff      = c(2, 4)
)
rm(expr_younger)
opts_bio <- updateBioOptions(opts_bio, ngenes = 5000)

opts_design <- CircadianDesignOptions(
  sample_sizes = c(20, 40, 60, 80, 100, 120, 140, 160),
  nsims        = 50,
  design       = "passive",
  cts          = times_young,
  test_types   = c("DP")
)

opts_analysis <- CircadianAnalysisOptions(
  alpha           = 0.05,
  p.adjust.method = "BH",
  parallel.ncores = 1,
  amp.cutoff      = 0,
  target_effect   = 0.1,
  fdr_thresholds  = c(0.01, 0.05, 0.10, 0.20),
  reference_n     = 60
)

base_out <- file.path("output", "run_panelC160")
dir.create(base_out, recursive = TRUE, showWarnings = FALSE)

cat("====================================================================\n")
cat("ANALYSIS: DIFFERENTIAL PHASE (DP)\n")
cat("====================================================================\n\n")

opts_bio_DP <- updateBioOptions(opts_bio,
  prop_DR    = 0.00,
  prop_DP    = 0.15,
  prop_DA    = 0.00,
  phase_diff = c(-6, 6),
  amp_diff   = c(0.5, 2)
)

dp_power_raw <- runPowerAnalysis(opts_bio_DP, opts_design, opts_analysis,
                                 test_type = "DP")

dp_results_file <- file.path(base_out, "dp_power_raw_pvalues.rds")
save(dp_power_raw, file = dp_results_file)
cat(sprintf("\nDP results saved: %s\n", dp_results_file))
