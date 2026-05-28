#' Debug DP p-values - Check why null genes have tiny p-values

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

# By diff_type
for (dt in 0:5) {
  n_tested <- sum(tested & diff_type == dt)
  cat(sprintf("diff_type=%d: %d tested\n", dt, n_tested))
}

# Focus on null genes (diff_type = 1, rhythmic both same)
null_idx <- which(diff_type == 1 & tested)
cat("\n=== NULL GENES (diff_type=1, rhythmic both, same phase) ===\n")
cat("Total null genes tested:", length(null_idx), "\n\n")

# Check p-value distribution for nulls
pvals_null <- pval_DP[null_idx]
cat("Null p-value summary:\n")
print(quantile(pvals_null, probs = c(0, 0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99, 1)))

# Count tiny p-values
cat("\nP < 0.01:", sum(pvals_null < 0.01), "(", round(100*sum(pvals_null < 0.01)/length(pvals_null), 1), "%)\n")
cat("P < 0.05:", sum(pvals_null < 0.05), "(", round(100*sum(pvals_null < 0.05)/length(pvals_null), 1), "%)\n")
cat("P = 0:", sum(pvals_null == 0), "\n")

# Manual LR_diff test on a few null genes
cat("\n=== MANUAL TEST ON FEW NULL GENES ===\n")

if (length(null_idx) > 0) {
  # Test first few null genes manually
  for (i in head(null_idx, 5)) {
    expr1 <- sim_out$expr1[[1]][i, ]
    expr2 <- sim_out$expr2[[1]][i, ]

    result <- LR_diff(times_young, expr1, times_young, expr2, 24, FN = TRUE, type = "phase")

    cat(sprintf("\nGene %d:\n", i))
    cat(sprintf("  diff_type: %d\n", diff_type[i]))
    cat(sprintf("  phase_diff (truth): %.4f\n", effectsize$phase[i]))
    cat(sprintf("  p-value: %.6e\n", result$pvalue))
    cat(sprintf("  phase1 (fit): %.4f\n", result$phase1))
    cat(sprintf("  phase2 (fit): %.4f\n", result$phase2))
    cat(sprintf("  delta_phase (fit): %.4f\n", result$delta_phase))
  }
}

# Now test true DP genes
target_idx <- which(diff_type == 4 & tested)
cat("\n=== DP TARGET GENES (diff_type=4) ===\n")
cat("Total target genes tested:", length(target_idx), "\n")

if (length(target_idx) > 0) {
  pvals_target <- pval_DP[target_idx]
  cat("Target p-value summary:\n")
  print(quantile(pvals_target, probs = c(0, 0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99, 1)))

  cat("\nP < 0.01:", sum(pvals_target < 0.01), "(", round(100*sum(pvals_target < 0.01)/length(pvals_target), 1), "%)\n")
  cat("P < 0.05:", sum(pvals_target < 0.05), "(", round(100*sum(pvals_target < 0.05)/length(pvals_target), 1), "%)\n")
}
