#' DP Power Analysis - Varying Phase Shift
#' Test power for different phase shift magnitudes: 2, 4, 6, 8, 10, 12 hours
#' This validates that the DP implementation correctly detects larger shifts

setwd("/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")

source("code/options.R")
source("code/simulation.R")
source("code/LRTest_diff_phase.R")
source("code/detection.R")
source("code/runner.R")
source("code/power.R")

cat("====================================================================\n")
cat("DP POWER ANALYSIS - VARYING PHASE SHIFT\n")
cat("Testing power for phase shifts: 2, 4, 6, 8, 10, 12 hours\n")
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

# Test different phase shifts
phase_shifts <- c(2, 4, 6, 8, 10, 12)  # hours
sample_sizes <- c(20, 40, 60, 80, 100)   # fewer sizes for speed
nsims <- 30  # fewer sims for speed

cat(sprintf("Phase shifts: %s hours\n", paste(phase_shifts, collapse=", ")))
cat(sprintf("Sample sizes: %s\n", paste(sample_sizes, collapse=", ")))
cat(sprintf("Simulations: %d per scenario\n\n", nsims))

# Storage for results
results_by_shift <- list()

# =====================================================================
# RUN SIMULATIONS FOR EACH PHASE SHIFT
# =====================================================================

for (shift_idx in seq_along(phase_shifts)) {
  phase_shift <- phase_shifts[shift_idx]

  cat(sprintf("====================================================================\n"))
  cat(sprintf("PHASE SHIFT: %d hours\n", phase_shift))
  cat(sprintf("====================================================================\n\n"))

  # Storage for this phase shift
  pvalues_DP <- array(NA, dim = c(length(sample_sizes), 5000, nsims))
  is_target_list <- vector("list", nsims)
  is_null_list <- vector("list", nsims)

  for (j in seq_along(sample_sizes)) {
    n <- sample_sizes[j]
    cat(sprintf(">>> n = %d\n", n))

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

    pval_DP <- sim_out$pval_DP[, 1, ]
    pvalues_DP[j, , ] <- pval_DP

    # Store ground truth
    diff_type <- sim_out$diff_type[[1]]
    effectsize_phase <- sim_out$effectsize[[1]]$phase

    for (i in 1:nsims) {
      if (j == 1) {
        is_target_list[[i]] <- diff_type == 4  # DP genes
        is_null_list[[i]] <- diff_type == 1    # Rhythmic both, same phase
      }
    }

    # Calculate marginal power at FDR 0.05
    total_targets <- sum(is_target_list[[1]]) * nsims

    power_by_n <- sapply(1:nsims, function(i) {
      qvals <- p.adjust(pvalues_DP[j, , i], method = "BH")
      discoveries <- qvals <= 0.05
      TD <- sum(discoveries & is_target_list[[i]])
      TD / total_targets
    })

    cat(sprintf("  Mean power at FDR 5%%: %.1f%%\n", 100 * mean(power_by_n)))
  }

  # Store results
  results_by_shift[[shift_idx]] <- list(
    phase_shift = phase_shift,
    sample_sizes = sample_sizes,
    nsims = nsims,
    pvalues = pvalues_DP,
    is_target_list = is_target_list,
    is_null_list = is_null_list
  )
}

# =====================================================================
# CALCULATE POWER SUMMARY BY PHASE SHIFT
# =====================================================================

cat("\n====================================================================\n")
cat("POWER SUMMARY BY PHASE SHIFT\n")
cat("====================================================================\n\n")

power_summary <- matrix(NA, nrow = length(sample_sizes), ncol = length(phase_shifts))
colnames(power_summary) <- paste0(phase_shifts, "hr")
rownames(power_summary) <- paste0("n=", sample_sizes)

for (shift_idx in seq_along(phase_shifts)) {
  res <- results_by_shift[[shift_idx]]
  total_targets <- sum(res$is_target_list[[1]]) * res$nsims

  for (j in seq_along(sample_sizes)) {
    power_vals <- sapply(1:res$nsims, function(i) {
      qvals <- p.adjust(res$pvalues[j, , i], method = "BH")
      discoveries <- qvals <= 0.05
      TD <- sum(discoveries & res$is_target_list[[i]])
      TD / total_targets
    })
    power_summary[j, shift_idx] <- mean(power_vals)
  }
}

cat("Marginal Power at FDR 5%%:\n")
cat("(Percentage of true DP genes detected)\n\n")
print(round(100 * power_summary, 1))

# Also check null calibration (type I error)
cat("\n====================================================================\n")
cat("NULL CALIBRATION (Type I Error at FDR 5%%)\n")
cat("====================================================================\n\n")

null_summary <- matrix(NA, nrow = length(sample_sizes), ncol = length(phase_shifts))
colnames(null_summary) <- paste0(phase_shifts, "hr")
rownames(null_summary) <- paste0("n=", sample_sizes)

for (shift_idx in seq_along(phase_shifts)) {
  res <- results_by_shift[[shift_idx]]

  for (j in seq_along(sample_sizes)) {
    null_rates <- sapply(1:res$nsims, function(i) {
      qvals <- p.adjust(res$pvalues[j, , i], method = "BH")
      discoveries <- qvals <= 0.05
      FD <- sum(discoveries & res$is_null_list[[i]])
      total_nulls <- sum(res$is_null_list[[i]])
      FD / total_nulls
    })
    null_summary[j, shift_idx] <- mean(null_rates)
  }
}

cat("False Discovery Rate on Null Genes:\n")
cat("(Should be close to 5%% if well-calibrated)\n\n")
print(round(100 * null_summary, 1))

# =====================================================================
# SAVE RESULTS
# =====================================================================

cat("\n====================================================================\n")
cat("SAVING RESULTS\n")
cat("====================================================================\n\n")

if (!dir.exists("output/dp_power_by_shift")) {
  dir.create("output/dp_power_by_shift", recursive = TRUE)
}

dp_varying_shift <- list(
  phase_shifts = phase_shifts,
  sample_sizes = sample_sizes,
  nsims = nsims,
  results_by_shift = results_by_shift,
  power_summary = power_summary,
  null_summary = null_summary
)

save(dp_varying_shift, file = "output/dp_power_by_shift/dp_varying_phase_shift.rds")
cat("Results saved: output/dp_power_by_shift/dp_varying_phase_shift.rds\n\n")

cat("\nAll simulations complete!\n")
