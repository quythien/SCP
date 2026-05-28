#' DP Power Analysis - Phase Shift × r Stratification
#' Varying both phase shift magnitude and signal-to-noise ratio
#' Save results for easy replotting

setwd("/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")

source("code/options.R")
source("code/simulation.R")
source("code/detection.R")
source("code/runner.R")
source("code/power.R")

cat("====================================================================\n")
cat("DP POWER ANALYSIS - PHASE SHIFT × r STRATIFICATION\n")
cat("Varying phase shift magnitude and signal-to-noise ratio\n")
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
phase_shifts <- c(0, 0.5, 1, 2, 4, 6, 8, 10, 12)
nsims <- 50
target_effect <- 0.1  # Threshold for "interesting" DP genes

cat(sprintf("Sample sizes: %s\n", paste(sample_sizes, collapse=", ")))
cat(sprintf("Phase shifts: %s hours\n", paste(phase_shifts, collapse=", ")))
cat(sprintf("Simulations: %d per scenario\n", nsims))
cat(sprintf("Target effect threshold: %.2f\n\n", target_effect))

# r = A/σ strata (signal-to-noise ratio)
r_strata <- c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.5, 3, 3.5, 4, 4.5, 5, Inf)
strata_labels <- c("(0,0.25]", "(0.25,0.5]", "(0.5,0.75]", "(0.75,1]",
                  "(1,1.25]", "(1.25,1.5]", "(1.5,1.75]", "(1.75,2]",
                  "(2,2.5]", "(2.5,3]", "(3,3.5]", "(3.5,4]",
                  "(4,4.5]", "(4.5,5]", ">5")

n_r_strata <- length(r_strata) - 1

# =====================================================================
# RUN SIMULATIONS
# =====================================================================

cat("====================================================================\n")
cat("RUNNING SIMULATIONS...\n")
cat("====================================================================\n\n")

# Storage arrays: [phase_shift, sample_size, r_stratum, sim]
strat_power <- array(NA, dim = c(length(phase_shifts), length(sample_sizes), n_r_strata, nsims))
strat_TD <- array(NA, dim = c(length(phase_shifts), length(sample_sizes), n_r_strata, nsims))
strat_FD <- array(NA, dim = c(length(phase_shifts), length(sample_sizes), n_r_strata, nsims))
strat_n_targets <- array(NA, dim = c(length(phase_shifts), length(sample_sizes), n_r_strata, nsims))

# Total scenarios
total_scenarios <- length(phase_shifts) * length(sample_sizes) * nsims
current_scenario <- 0

for (p in seq_along(phase_shifts)) {
  phase_shift <- phase_shifts[p]
  cat(sprintf(">>> Phase shift: %.1f hours\n", phase_shift))

  for (j in seq_along(sample_sizes)) {
    n <- sample_sizes[j]
    cat(sprintf("  >>> n = %d\n", n))

    sim_out <- runSimsDiff(
      sample_sizes = c(n),
      nsims = nsims,
      ngenes = 5000,
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

    fdr_DP <- sim_out$fdr_DP[, 1, ]  # genes x nsims

    for (i in 1:nsims) {
      current_scenario <- current_scenario + 1
      if (current_scenario %% 50 == 0) {
        cat(sprintf("    Progress: %d/%d scenarios\n", current_scenario, total_scenarios))
      }

      # Get ground truth
      diff_type <- sim_out$diff_type[[i]]  # 4 = DP genes
      effectsize_phase <- sim_out$effectsize[[i]]$phase
      effectsize_DR1 <- sim_out$effectsize[[i]]$DR1  # A1/sigma = r for group 1
      effectsize_DR2 <- sim_out$effectsize[[i]]$DR2  # A2/sigma = r for group 2

      # For DP: use r = min(A1/sigma, A2/sigma)
      r_values <- pmin(effectsize_DR1, effectsize_DR2)
      is_DP <- diff_type == 4

      # Zg: All DP genes
      Zg <- ifelse(is_DP, 1, 0)

      # Zg2: Target DP genes (effect >= target_effect)
      Zg2 <- ifelse(is_DP & effectsize_phase >= target_effect, 1, 0)

      # Stratify ALL genes by r
      r_for_strat <- r_values
      r_for_strat[!is_DP] <- 0  # Non-DP genes go to lowest bin
      xgr <- cut(r_for_strat, breaks = r_strata, include.lowest = TRUE, labels = FALSE)
      xgr[!is_DP] <- NA  # Remove non-DP from strata

      # Get discoveries
      discoveries <- fdr_DP[, i] <= 0.05

      for (k in 1:n_r_strata) {
        in_stratum <- xgr == k

        # True Discoveries
        TD <- sum(discoveries & Zg2 == 1 & in_stratum, na.rm = TRUE)

        # False Discoveries
        FD <- sum(discoveries & Zg == 0 & in_stratum, na.rm = TRUE)

        # Target genes
        n_targets <- sum(Zg2 == 1 & in_stratum, na.rm = TRUE)

        # Store results
        strat_power[p, j, k, i] <- if (n_targets > 0) TD / n_targets else NA
        strat_TD[p, j, k, i] <- TD
        strat_FD[p, j, k, i] <- FD
        strat_n_targets[p, j, k, i] <- n_targets
      }
    }
  }
}

cat("\n====================================================================\n")
cat("SIMULATIONS COMPLETE\n")
cat("====================================================================\n\n")

# =====================================================================
# SAVE RESULTS (before plotting, for debugging)
# =====================================================================

dp_phase_results <- list(
  sample_sizes = sample_sizes,
  phase_shifts = phase_shifts,
  nsims = nsims,
  target_effect = target_effect,
  r_strata = r_strata,
  strata_labels = strata_labels,
  strat_power = strat_power,
  strat_TD = strat_TD,
  strat_FD = strat_FD,
  strat_n_targets = strat_n_targets
)

save(dp_phase_results, file = "output/dp_power_phase_shift/dp_power_phase_shift_results.rds", ascii = TRUE)
cat("Results saved: output/dp_power_phase_shift/dp_power_phase_shift_results.rds\n\n")

# =====================================================================
# QUICK SUMMARY CHECK
# =====================================================================

cat("MARGINAL POWER SUMMARY (across all r strata):\n")
cat("==============================================\n\n")

cat("Phase Shift | ")
for (n in sample_sizes) {
  cat(sprintf(" n=%d |", n))
}
cat("\n")
cat(paste0(rep("-", 20 + length(sample_sizes)*7), collapse = ""), "\n")

for (p in seq_along(phase_shifts)) {
  cat(sprintf("%11s | ", phase_shifts[p]))
  for (j in seq_along(sample_sizes)) {
    marginal_TD <- sum(strat_TD[p, j, , ], na.rm = TRUE)
    marginal_targets <- sum(strat_n_targets[p, j, , ], na.rm = TRUE)
    marginal_power <- if (marginal_targets > 0) marginal_TD / marginal_targets else NA
    cat(sprintf(" %4.1f%% |", 100 * marginal_power))
  }
  cat("\n")
}

cat("\nAll simulations complete!\n")
