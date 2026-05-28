#' Debug DP p-values - Check amplitude vs p-value relationship

setwd("/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")

source("code/options.R")
source("code/simulation.R")
source("code/detection.R")
source("code/runner.R")

# Load real data for time distribution
COMBINED <- readRDS("/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Collaborative/Paper/Congruence/PNAS_aging/data/combined_data.rds")

expr_sample_names <- colnames(COMBINED$expr)
pheno_order <- match(expr_sample_names, COMBINED$pheno$sample_name)
valid_samples <- !is.na(pheno_order)
expr_sample_names <- expr_sample_names[valid_samples]
COMBINED$expr <- COMBINED$expr[, valid_samples]
pheno_order <- pheno_order[valid_samples]
pheno_data <- COMBINED$pheno[pheno_order, ]
pheno_data$tod <- if("TOD.x" %in% colnames(pheno_data)) pheno_data$TOD.x else pheno_data$TOD.y
pheno_data$age_group_final <- if("AgeGroup" %in% colnames(pheno_data)) pheno_data$AgeGroup else pheno_data$age_group
complete_samples <- !is.na(pheno_data$age_group_final) & !is.na(pheno_data$tod) &
  pheno_data$age_group_final %in% c("younger", "older")
pheno_clean <- pheno_data[complete_samples, ]
younger_idx <- pheno_clean$age_group_final == "younger"
times_young <- pheno_clean$tod[younger_idx]

cat("Reference data: n=", length(times_young), "\n\n")

# Run a small simulation to debug
set.seed(123)
n <- 60
ngenes <- 500
nsims <- 1
phase_shift <- 6

cat("Running debug simulation...\n")

sim_out <- runSimsDiff(
  sample_sizes = c(n),
  nsims = nsims,
  ngenes = ngenes,
  prop_rhythmic = 0.30,
  prop_DR = 0.00,
  prop_DP = 0.15,
  prop_DA = 0.00,
  phase_diff = c(-phase_shift, phase_shift),
  amp_diff = c(0.5, 2),
  cts = times_young,
  design = "passive",
  test_types = c("DP"),
  verbose = FALSE
)

# Analyze the first simulation
pval_DP <- sim_out$pval_DP[, 1, 1]
diff_type <- sim_out$diff_type[[1]]
effectsize <- sim_out$effectsize[[1]]

cat("\n=== DP P-VALUE DIAGNOSTICS ===\n\n")

# Genes that were tested (p-value < 1)
tested <- pval_DP < 1
cat("Genes tested (p < 1):", sum(tested), "\n\n")

# Focus on null genes (diff_type = 1, rhythmic both same)
null_idx <- which(diff_type == 1 & tested)
cat("=== NULL GENES (diff_type=1, rhythmic both, same phase) ===\n")
cat("Total null genes tested:", length(null_idx), "\n\n")

# Check r values for nulls
r_null <- pmin(effectsize$DR1[null_idx], effectsize$DR2[null_idx])

cat("Signal-to-noise (r = A/sigma) summary for nulls:\n")
print(summary(r_null))

# Count tiny p-values
pvals_null <- pval_DP[null_idx]
cat("\nP < 0.01:", sum(pvals_null < 0.01), "(", round(100*sum(pvals_null < 0.01)/length(pvals_null), 1), "%)\n")
cat("P < 0.05:", sum(pvals_null < 0.05), "(", round(100*sum(pvals_null < 0.05)/length(pvals_null), 1), "%)\n")
cat("P = 0:", sum(pvals_null == 0), "\n")

# Now check: what's the minimum r for tested genes?
cat("\n=== SIGNAL-TO-NOISE FOR TESTED GENES ===\n")
cat("All tested genes (n=", sum(tested), "):\n")
r_tested <- pmin(effectsize$DR1[tested], effectsize$DR2[tested])
print(summary(r_tested))

# The key issue: genes with very low amplitude should NOT be tested for phase difference
# because the phase estimate is unreliable
cat("\n=== RECOMMENDATION ===\n")
cat("Genes with r < 0.5 should be EXCLUDED from DP testing\n")
cat("Current: all 'both rhythmic' genes are tested regardless of amplitude\n")
cat("This leads to unstable phase estimates and tiny p-values for noise\n")
