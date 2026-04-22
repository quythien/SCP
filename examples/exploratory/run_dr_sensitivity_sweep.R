# DR sensitivity sweep (amplitude scaling)
# Usage:
#   Rscript examples/run_dr_sensitivity_sweep.R

set.seed(12345)

setwd("/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")
source_dir <- file.path(getwd(), "code")
old_wd <- setwd(source_dir)
source("setup.R")
setwd(old_wd)

source("code/plot_with_se.R")

# Sweep config
amp_scales <- c(1.0, 1.25, 1.5, 2.0)

# Output dir
base_out <- file.path("output", "run_dr_sensitivity_20260223")
out_dir  <- file.path(base_out, "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Load pilot data
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

# Estimate parameters
cat("Estimating circadian parameters from pilot data...\n\n")
opts_bio_base <- estCircadianParam(
  data          = expr_younger,
  times         = times_young,
  period        = 24,
  prop_DR       = 0.15,
  prop_DP       = 0.10,
  phase_diff    = c(-6, 6),
  amp_diff      = c(2, 4)
)
rm(expr_younger)
opts_bio_base <- updateBioOptions(opts_bio_base, ngenes = 5000)

# Design & analysis
opts_design <- CircadianDesignOptions(
  sample_sizes = c(20, 40, 60, 80, 100, 120, 140, 160),
  nsims        = 50,
  design       = "passive",
  cts          = times_young,
  test_types   = c("DR")
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

summary_table <- data.frame()

for (s in amp_scales) {
  cat(sprintf("\n=== DR sensitivity: amplitude scale = %.2f ===\n", s))
  opts_bio <- opts_bio_base
  opts_bio$amplitude <- opts_bio$amplitude * s

  dr_power_raw <- runPowerAnalysis(opts_bio, opts_design, opts_analysis, test_type = "DR")

  out_rds <- file.path(base_out, sprintf("dr_power_raw_ampScale_%0.2f.rds", s))
  save(dr_power_raw, file = out_rds)
  cat(sprintf("Saved: %s\n", out_rds))

  out_fig <- file.path(out_dir, sprintf("dr_power_ampScale_%0.2f.pdf", s))
  plotWithSE(out_rds, out_fig, test_name = sprintf("DR (A scale %.2f)", s), analysis.opts = opts_analysis)

  mp <- rowMeans(dr_power_raw$marginal_power, na.rm = TRUE)
  summary_table <- rbind(summary_table, data.frame(
    amp_scale = s,
    n = dr_power_raw$sample_sizes,
    marginal_power = mp
  ))
}

summary_path <- file.path(base_out, "dr_sensitivity_summary.csv")
write.csv(summary_table, summary_path, row.names = FALSE)
cat(sprintf("\nSummary saved: %s\n", summary_path))
