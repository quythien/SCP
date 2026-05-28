#' DA Power Analysis - PROPER-style Stratified Power
#' Differential Amplitude - Stratified by r = A/σ
#' Save results for easy replotting

setwd("/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")

source("code/options.R")
source("code/simulation.R")
source("code/LRTest_diff_amp.R")
source("code/LRTest_diff_phase.R")
source("code/detection.R")
source("code/runner.R")
source("code/power.R")

cat("====================================================================\n")
cat("DA POWER ANALYSIS - PROPER-STYLE STRATIFIED POWER\n")
cat("Stratifying by amplitude (signal strength)\n")
cat("====================================================================\n\n")

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

cat(sprintf("Reference data: n=%d\n\n", length(times_young)))

# =====================================================================
# SIMULATION PARAMETERS
# =====================================================================

sample_sizes <- c(10, 20, 40, 60, 80, 100)
nsims <- 50
target_effect <- 0.1  # Threshold for "interesting" DA genes

cat(sprintf("Sample sizes: %s\n", paste(sample_sizes, collapse=", ")))
cat(sprintf("Simulations: %d per scenario\n", nsims))
cat(sprintf("Target effect threshold: %.2f\n\n", target_effect))

# r = A/σ strata (signal-to-noise ratio, like CircaPower)
r_strata <- c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.5, 3, 3.5, 4, 4.5, 5, Inf)
strata_labels <- c("(0,0.25]", "(0.25,0.5]", "(0.5,0.75]", "(0.75,1]",
                  "(1,1.25]", "(1.25,1.5]", "(1.5,1.75]", "(1.75,2]",
                  "(2,2.5]", "(2.5,3]", "(3,3.5]", "(3.5,4]",
                  "(4,4.5]", "(4.5,5]", ">5")

# =====================================================================
# RUN SIMULATIONS - DA (Differential Amplitude)
# =====================================================================

cat("====================================================================\n")
cat("RUNNING SIMULATIONS...\n")
cat("====================================================================\n\n")

# Storage arrays: [sample_size, r_stratum, sim]
strat_power <- array(NA, dim = c(length(sample_sizes), length(r_strata)-1, nsims))
strat_TD <- array(NA, dim = c(length(sample_sizes), length(r_strata)-1, nsims))
strat_FD <- array(NA, dim = c(length(sample_sizes), length(r_strata)-1, nsims))
strat_n_targets <- array(NA, dim = c(length(sample_sizes), length(r_strata)-1, nsims))
strat_n_null <- array(NA, dim = c(length(sample_sizes), length(r_strata)-1, nsims))

# Store r values per gene for visualization
r_list <- vector("list", nsims)

for (j in seq_along(sample_sizes)) {
  n <- sample_sizes[j]
  cat(sprintf(">>> n = %d\n", n))

  sim_out <- runSimsDiff(
    sample_sizes = c(n),
    nsims = nsims,
    ngenes = 5000,
    prop_rhythmic = 0.30,
    prop_DR = 0.00,
    prop_DP = 0.00,
    prop_DA = 0.15,  # 15% DA genes
    phase_diff = c(0, 0),
    amp_diff = c(0.5, 2),  # Amplitude ratio between 0.5 and 2
    cts = times_young,
    design = "passive",
    test_types = c("DA"),
    verbose = FALSE
  )

  fdr_DA <- sim_out$fdr_DA[, 1, ]  # genes x nsims

  for (i in 1:nsims) {
    # Get ground truth
    diff_type <- sim_out$diff_type[[i]]  # 5 = DA genes
    effectsize_amp <- sim_out$effectsize[[i]]$amp
    effectsize_DR1 <- sim_out$effectsize[[i]]$DR1  # A1/sigma = r for group 1
    effectsize_DR2 <- sim_out$effectsize[[i]]$DR2  # A2/sigma = r for group 2

    # For DA: use r = min(A1/sigma, A2/sigma)
    r_values <- pmin(effectsize_DR1, effectsize_DR2)

    is_DA <- diff_type == 5  # Type 5 = DA genes

    # Store r distribution for first sample size
    if (j == 1) {
      r_list[[i]] <- r_values[is_DA]
    }

    # Zg: All DA genes
    Zg <- ifelse(is_DA, 1, 0)

    # Zg2: Target DA genes (effect >= target_effect)
    Zg2 <- ifelse(is_DA & abs(effectsize_amp) >= target_effect, 1, 0)

    # Stratify by r
    r_for_strat <- r_values
    r_for_strat[!is_DA] <- 0
    n_strata <- length(r_strata) - 1
    xgr <- cut(r_for_strat, breaks = r_strata, include.lowest = TRUE, labels = FALSE)
    xgr[!is_DA] <- NA

    # Get discoveries
    discoveries <- fdr_DA[, i] <= 0.05

    for (k in 1:n_strata) {
      in_stratum <- xgr == k

      TD <- sum(discoveries & Zg2 == 1 & in_stratum, na.rm = TRUE)
      FD <- sum(discoveries & Zg == 0 & in_stratum, na.rm = TRUE)
      n_targets <- sum(Zg2 == 1 & in_stratum, na.rm = TRUE)
      n_null <- sum(Zg == 0 & in_stratum, na.rm = TRUE)

      if (n_targets > 0) {
        strat_power[j, k, i] <- TD / n_targets
      } else {
        strat_power[j, k, i] <- NA
      }

      strat_TD[j, k, i] <- TD
      strat_FD[j, k, i] <- FD
      strat_n_targets[j, k, i] <- n_targets
      strat_n_null[j, k, i] <- n_null
    }
  }

  cat(sprintf("  n=%d: ", n))
  for (k in 1:n_strata) {
    mean_p <- mean(strat_power[j, k, ], na.rm = TRUE)
    cat(sprintf("%s=%.0f%% ", strata_labels[k], 100 * mean_p))
  }
  cat("\n\n")
}

# =====================================================================
# SAVE RESULTS
# =====================================================================

cat("====================================================================\n")
cat("SAVING RESULTS\n")
cat("====================================================================\n\n")

da_power_results <- list(
  sample_sizes = sample_sizes,
  nsims = nsims,
  target_effect = target_effect,
  r_strata = r_strata,
  strata_labels = strata_labels,
  strat_power = strat_power,
  strat_TD = strat_TD,
  strat_FD = strat_FD,
  strat_n_targets = strat_n_targets,
  strat_n_null = strat_n_null,
  r_list = r_list
)

save(da_power_results, file = "output/da_power_stratified/da_power_signal_stratified_results.rds", ascii = TRUE)
cat("Results saved: output/da_power_stratified/da_power_signal_stratified_results.rds\n\n")

# =====================================================================
# QUICK SUMMARY
# =====================================================================

cat("MARGINAL POWER SUMMARY:\n")
cat("========================\n\n")

cat("Sample Size | Marginal Power\n")
cat("---------------------------\n")
for (j in seq_along(sample_sizes)) {
  marginal_TD <- sum(strat_TD[j, , ], na.rm = TRUE)
  marginal_targets <- sum(strat_n_targets[j, , ], na.rm = TRUE)
  marginal_power <- if (marginal_targets > 0) marginal_TD / marginal_targets else NA
  cat(sprintf("n = %-7d | %.1f%%\n", sample_sizes[j], 100 * marginal_power))
}

cat("\nAll simulations complete!\n")
