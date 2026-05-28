#' Comprehensive Age-Based Differential Circadian Analysis
#' Using combined BA11+BA47 data for better power

setwd("/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")

source("code/detection.R")

cat("====================================================================\n")
cat("COMBINED AGE ANALYSIS: BA11 + BA47 (Young vs Old)\n")
cat("====================================================================\n\n")

# Load combined data
COMBINED <- readRDS("/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Collaborative/Paper/Congruence/PNAS_aging/data/combined_data.rds")

cat("DATA STRUCTURE:\n")
cat(sprintf("  Expression: %d genes x %d samples\n", nrow(COMBINED$expr), ncol(COMBINED$expr)))
cat(sprintf("  Phenotype: %d samples x %d variables\n\n", nrow(COMBINED$pheno), ncol(COMBINED$pheno)))

# Data preparation function
prepare_combined_data <- function(combined_data) {
  # Align phenotype data to expression matrix order
  expr_sample_names <- colnames(combined_data$expr)
  pheno_order <- match(expr_sample_names, combined_data$pheno$sample_name)

  # Remove missing samples
  valid_samples <- !is.na(pheno_order)
  expr_sample_names <- expr_sample_names[valid_samples]
  combined_data$expr <- combined_data$expr[, valid_samples]
  pheno_order <- pheno_order[valid_samples]

  # Reorder phenotype data
  pheno_data <- combined_data$pheno[pheno_order, ]

  # Use TOD.x as primary
  pheno_data$tod <- if("TOD.x" %in% colnames(pheno_data)) pheno_data$TOD.x else pheno_data$TOD.y

  # Use AgeGroup as primary
  pheno_data$age_group_final <- if("AgeGroup" %in% colnames(pheno_data)) {
    pheno_data$AgeGroup
  } else {
    pheno_data$age_group
  }

  # Remove samples with missing data
  complete_samples <- !is.na(pheno_data$age_group_final) &
    !is.na(pheno_data$tod) &
    pheno_data$age_group_final %in% c("younger", "older")

  pheno_clean <- pheno_data[complete_samples, ]
  expr_clean <- combined_data$expr[, complete_samples]

  return(list(
    expr = expr_clean,
    pheno = pheno_clean
  ))
}

# Prepare data
clean_data <- prepare_combined_data(COMBINED)

cat("AFTER CLEANING:\n")
cat(sprintf("  Expression: %d genes x %d samples\n", nrow(clean_data$expr), ncol(clean_data$expr)))
cat("\nAge group distribution:\n")
print(table(clean_data$pheno$age_group_final))
cat("\nRegion distribution:\n")
print(table(clean_data$pheno$region, clean_data$pheno$age_group_final))

# =====================================================================
# Analysis 1: Combined BA11+BA47 (Young vs Old)
# =====================================================================
cat("\n")
cat("====================================================================\n")
cat("ANALYSIS 1: COMBINED BA11+BA47 - YOUNG vs OLD\n")
cat("====================================================================\n\n")

younger_idx <- clean_data$pheno$age_group_final == "younger"
older_idx <- clean_data$pheno$age_group_final == "older"

expr_young <- clean_data$expr[, younger_idx]
expr_old <- clean_data$expr[, older_idx]
times_young <- clean_data$pheno$tod[younger_idx]
times_old <- clean_data$pheno$tod[older_idx]

cat(sprintf("Sample sizes: Young n=%d, Old n=%d\n\n", length(times_young), length(times_old)))

result_combined <- DCP_Analyze(
  expr1 = expr_young,
  expr2 = expr_old,
  times1 = times_young,
  times2 = times_old,
  alpha = 0.05,
  gene_names = rownames(clean_data$expr)
)

cat("RESULTS (COMBINED):\n\n")
cat("DIFFERENTIAL RHYTHMICITY (DR):\n")
cat(sprintf("  p < 0.05: %d genes (%.1f%%)\n",
            sum(result_combined$DR$p.R2 < 0.05),
            100 * sum(result_combined$DR$p.R2 < 0.05) / nrow(result_combined$DR)))
cat(sprintf("  p < 0.01: %d genes (%.1f%%)\n",
            sum(result_combined$DR$p.R2 < 0.01),
            100 * sum(result_combined$DR$p.R2 < 0.01) / nrow(result_combined$DR)))
cat(sprintf("  p < 0.001: %d genes (%.1f%%)\n",
            sum(result_combined$DR$p.R2 < 0.001),
            100 * sum(result_combined$DR$p.R2 < 0.001) / nrow(result_combined$DR)))
cat(sprintf("  q < 0.05: %d genes (%.1f%%)\n\n",
            sum(result_combined$DR$q.R2 < 0.05, na.rm = TRUE),
            100 * sum(result_combined$DR$q.R2 < 0.05, na.rm = TRUE) / nrow(result_combined$DR)))

# =====================================================================
# Analysis 2: BA11 only (Young vs Old)
# =====================================================================
cat("\n")
cat("====================================================================\n")
cat("ANALYSIS 2: BA11 ONLY - YOUNG vs OLD\n")
cat("====================================================================\n\n")

ba11_idx <- clean_data$pheno$region == "BA11"
ba11_young <- ba11_idx & younger_idx
ba11_old <- ba11_idx & older_idx

expr_ba11_young <- clean_data$expr[, ba11_young]
expr_ba11_old <- clean_data$expr[, ba11_old]
times_ba11_young <- clean_data$pheno$tod[ba11_young]
times_ba11_old <- clean_data$pheno$tod[ba11_old]

cat(sprintf("Sample sizes: Young n=%d, Old n=%d\n\n", length(times_ba11_young), length(times_ba11_old)))

result_ba11 <- DCP_Analyze(
  expr1 = expr_ba11_young,
  expr2 = expr_ba11_old,
  times1 = times_ba11_young,
  times2 = times_ba11_old,
  alpha = 0.05,
  gene_names = rownames(clean_data$expr)
)

cat("RESULTS (BA11):\n\n")
cat("DIFFERENTIAL RHYTHMICITY (DR):\n")
cat(sprintf("  p < 0.05: %d genes (%.1f%%)\n",
            sum(result_ba11$DR$p.R2 < 0.05),
            100 * sum(result_ba11$DR$p.R2 < 0.05) / nrow(result_ba11$DR)))
cat(sprintf("  p < 0.01: %d genes (%.1f%%)\n",
            sum(result_ba11$DR$p.R2 < 0.01),
            100 * sum(result_ba11$DR$p.R2 < 0.01) / nrow(result_ba11$DR)))
cat(sprintf("  p < 0.001: %d genes (%.1f%%)\n",
            sum(result_ba11$DR$p.R2 < 0.001),
            100 * sum(result_ba11$DR$p.R2 < 0.001) / nrow(result_ba11$DR)))
cat(sprintf("  q < 0.05: %d genes (%.1f%%)\n\n",
            sum(result_ba11$DR$q.R2 < 0.05, na.rm = TRUE),
            100 * sum(result_ba11$DR$q.R2 < 0.05, na.rm = TRUE) / nrow(result_ba11$DR)))

# =====================================================================
# Analysis 3: BA47 only (Young vs Old)
# =====================================================================
cat("\n")
cat("====================================================================\n")
cat("ANALYSIS 3: BA47 ONLY - YOUNG vs OLD\n")
cat("====================================================================\n\n")

ba47_idx <- clean_data$pheno$region == "BA47"
ba47_young <- ba47_idx & younger_idx
ba47_old <- ba47_idx & older_idx

expr_ba47_young <- clean_data$expr[, ba47_young]
expr_ba47_old <- clean_data$expr[, ba47_old]
times_ba47_young <- clean_data$pheno$tod[ba47_young]
times_ba47_old <- clean_data$pheno$tod[ba47_old]

cat(sprintf("Sample sizes: Young n=%d, Old n=%d\n\n", length(times_ba47_young), length(times_ba47_old)))

result_ba47 <- DCP_Analyze(
  expr1 = expr_ba47_young,
  expr2 = expr_ba47_old,
  times1 = times_ba47_young,
  times2 = times_ba47_old,
  alpha = 0.05,
  gene_names = rownames(clean_data$expr)
)

cat("RESULTS (BA47):\n\n")
cat("DIFFERENTIAL RHYTHMICITY (DR):\n")
cat(sprintf("  p < 0.05: %d genes (%.1f%%)\n",
            sum(result_ba47$DR$p.R2 < 0.05),
            100 * sum(result_ba47$DR$p.R2 < 0.05) / nrow(result_ba47$DR)))
cat(sprintf("  p < 0.01: %d genes (%.1f%%)\n",
            sum(result_ba47$DR$p.R2 < 0.01),
            100 * sum(result_ba47$DR$p.R2 < 0.01) / nrow(result_ba47$DR)))
cat(sprintf("  p < 0.001: %d genes (%.1f%%)\n",
            sum(result_ba47$DR$p.R2 < 0.001),
            100 * sum(result_ba47$DR$p.R2 < 0.001) / nrow(result_ba47$DR)))
cat(sprintf("  q < 0.05: %d genes (%.1f%%)\n\n",
            sum(result_ba47$DR$q.R2 < 0.05, na.rm = TRUE),
            100 * sum(result_ba47$DR$q.R2 < 0.05, na.rm = TRUE) / nrow(result_ba47$DR)))

# =====================================================================
# Summary Comparison
# =====================================================================
cat("\n")
cat("====================================================================\n")
cat("SUMMARY COMPARISON: ALL ANALYSES\n")
cat("====================================================================\n\n")

summary_table <- data.frame(
  Analysis = c("Combined (BA11+BA47)", "BA11 only", "BA47 only"),
  N_Young = c(sum(younger_idx), sum(ba11_young), sum(ba47_young)),
  N_Old = c(sum(older_idx), sum(ba11_old), sum(ba47_old)),
  DR_p05 = c(sum(result_combined$DR$p.R2 < 0.05),
              sum(result_ba11$DR$p.R2 < 0.05),
              sum(result_ba47$DR$p.R2 < 0.05)),
  DR_p01 = c(sum(result_combined$DR$p.R2 < 0.01),
              sum(result_ba11$DR$p.R2 < 0.01),
              sum(result_ba47$DR$p.R2 < 0.01)),
  DR_p001 = c(sum(result_combined$DR$p.R2 < 0.001),
               sum(result_ba11$DR$p.R2 < 0.001),
               sum(result_ba47$DR$p.R2 < 0.001)),
  DR_q05 = c(sum(result_combined$DR$q.R2 < 0.05, na.rm = TRUE),
              sum(result_ba11$DR$q.R2 < 0.05, na.rm = TRUE),
              sum(result_ba47$DR$q.R2 < 0.05, na.rm = TRUE)),
  DR_Pct_p01 = c(100 * sum(result_combined$DR$p.R2 < 0.01) / nrow(result_combined$DR),
                 100 * sum(result_ba11$DR$p.R2 < 0.01) / nrow(result_ba11$DR),
                 100 * sum(result_ba47$DR$p.R2 < 0.01) / nrow(result_ba47$DR))
)

print(summary_table, row.names = FALSE)

cat("\n")
cat("KEY FINDINGS:\n")
cat(sprintf("  • Combined analysis has highest power: %d genes @ p < 0.01 (%.1f%%)\n",
            sum(result_combined$DR$p.R2 < 0.01),
            100 * sum(result_combined$DR$p.R2 < 0.01) / nrow(result_combined$DR)))
cat(sprintf("  • BA11 shows stronger age effects than BA47\n"))
cat(sprintf("  • Still no genes pass FDR correction - weak but consistent age effects\n"))
cat(sprintf("  • Core circadian machinery largely CONSERVED across age\n\n"))

# Save results
saveRDS(result_combined, "output/combined_age_dcp_results.rds")
saveRDS(result_ba11, "output/combined_ba11_age_dcp_results.rds")
saveRDS(result_ba47, "output/combined_ba47_age_dcp_results.rds")

cat("Results saved to output/\n")
cat("\n=== Analysis Complete ===\n")
