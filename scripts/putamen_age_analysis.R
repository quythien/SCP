#' Age-Based Differential Circadian Analysis on Putamen Data
#' Compare Young vs Old controls for differential rhythmicity

setwd("/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")

source("code/detection.R")

cat("====================================================================\n")
cat("AGE-BASED DIFFERENTIAL CIRCADIAN ANALYSIS: PUTAMEN CONTROLS\n")
cat("====================================================================\n\n")

# Load data
put_expr <- read.csv("/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Collaborative/NPC/data/Putamen_CPMfiltered_logCPM_1215_rm97_rm231.csv", row.names = 1)
put_clinical <- read.csv("/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Collaborative/NPC/data/DS_clinical_1221_rm97_rm231_matchIndex34.csv", row.names = 1)

# Filter for controls only
ctrl_idx <- which(put_clinical$Diagnostic.Category == "CONTROL")
ctrl_clinical <- put_clinical[ctrl_idx, ]
ctrl_expr <- put_expr[, ctrl_idx]

cat("CONTROL SAMPLES:\n")
cat(sprintf("  Total: %d\n", length(ctrl_idx)))
cat(sprintf("  Age range: %.1f - %.1f years\n", min(ctrl_clinical$Age), max(ctrl_clinical$Age)))
cat(sprintf("  Mean age: %.1f years\n", mean(ctrl_clinical$Age)))
cat(sprintf("  Median age: %.1f years\n\n", median(ctrl_clinical$Age)))

# Define age groups (median split)
age_median <- median(ctrl_clinical$Age)
young_idx <- which(ctrl_clinical$Age <= age_median)
old_idx <- which(ctrl_clinical$Age > age_median)

cat("AGE GROUPS:\n")
cat(sprintf("  Young (age <= %.1f): n = %d\n", age_median, length(young_idx)))
cat(sprintf("    Age range: %.1f - %.1f years\n", min(ctrl_clinical$Age[young_idx]), max(ctrl_clinical$Age[young_idx])))
cat(sprintf("    Mean age: %.1f years\n", mean(ctrl_clinical$Age[young_idx])))
cat(sprintf("  Old (age > %.1f): n = %d\n", age_median, length(old_idx)))
cat(sprintf("    Age range: %.1f - %.1f years\n", min(ctrl_clinical$Age[old_idx]), max(ctrl_clinical$Age[old_idx])))
cat(sprintf("    Mean age: %.1f years\n\n", mean(ctrl_clinical$Age[old_idx])))

# Get expression and time for each group
expr_young <- ctrl_expr[, young_idx]
expr_old <- ctrl_expr[, old_idx]
times_young <- ctrl_clinical$CorrectedTOD[young_idx]
times_old <- ctrl_clinical$CorrectedTOD[old_idx]

cat("TIME OF DEATH DISTRIBUTIONS:\n")
cat(sprintf("  Young: %.1f - %.1f hours\n", min(times_young), max(times_young)))
cat(sprintf("  Old: %.1f - %.1f hours\n\n", min(times_old), max(times_old)))

# Run DCP_Analyze
cat("====================================================================\n")
cat("RUNNING DIFFERENTIAL CIRCADIAN ANALYSIS\n")
cat("====================================================================\n\n")

result <- DCP_Analyze(
  expr1 = expr_young,
  expr2 = expr_old,
  times1 = times_young,
  times2 = times_old,
  alpha = 0.05,
  gene_names = rownames(put_expr)
)

# Summary results
cat("====================================================================\n")
cat("RESULTS: YOUNG vs OLD CONTROLS\n")
cat("====================================================================\n\n")

# DR results
cat("DIFFERENTIAL RHYTHMICITY (DR):\n")
cat(sprintf("  p < 0.05: %d genes (%.1f%%)\n",
            sum(result$DR$p.R2 < 0.05),
            100 * sum(result$DR$p.R2 < 0.05) / nrow(result$DR)))
cat(sprintf("  p < 0.01: %d genes (%.1f%%)\n",
            sum(result$DR$p.R2 < 0.01),
            100 * sum(result$DR$p.R2 < 0.01) / nrow(result$DR)))
cat(sprintf("  q < 0.05: %d genes (%.1f%%)\n\n",
            sum(result$DR$q.R2 < 0.05, na.rm = TRUE),
            100 * sum(result$DR$q.R2 < 0.05, na.rm = TRUE) / nrow(result$DR)))

# DP results
cat("DIFFERENTIAL PHASE (DP):\n")
cat(sprintf("  p < 0.05: %d genes (%.1f%%)\n",
            sum(result$DP$p.phase < 0.05),
            100 * sum(result$DP$p.phase < 0.05) / nrow(result$DP)))
cat(sprintf("  p < 0.01: %d genes (%.1f%%)\n",
            sum(result$DP$p.phase < 0.01),
            100 * sum(result$DP$p.phase < 0.01) / nrow(result$DP)))
cat(sprintf("  q < 0.05: %d genes (%.1f%%)\n\n",
            sum(result$DP$q.phase < 0.05, na.rm = TRUE),
            100 * sum(result$DP$q.phase < 0.05, na.rm = TRUE) / nrow(result$DP)))

# DA results
cat("DIFFERENTIAL AMPLITUDE (DA):\n")
cat(sprintf("  p < 0.05: %d genes (%.1f%%)\n",
            sum(result$DA$p.amp < 0.05),
            100 * sum(result$DA$p.amp < 0.05) / nrow(result$DA)))
cat(sprintf("  p < 0.01: %d genes (%.1f%%)\n",
            sum(result$DA$p.amp < 0.01),
            100 * sum(result$DA$p.amp < 0.01) / nrow(result$DA)))
cat(sprintf("  q < 0.05: %d genes (%.1f%%)\n\n",
            sum(result$DA$q.amp < 0.05, na.rm = TRUE),
            100 * sum(result$DA$q.amp < 0.05, na.rm = TRUE) / nrow(result$DA)))

# Joint rhythmicity classification
cat("JOINT RHYTHMICITY CLASSIFICATION:\n")
joint_table <- table(result$classification)
print(joint_table)

cat("\nINTERPRETATION:\n")
dr_p05 <- sum(result$DR$p.R2 < 0.05)
dr_q05 <- sum(result$DR$q.R2 < 0.05, na.rm = TRUE)

if (dr_q05 > 0) {
  cat(sprintf("  • %d genes show differential rhythmicity between young and old\n", dr_q05))
  cat(sprintf("  • This suggests AGE affects circadian rhythms in Putamen\n"))
} else if (dr_p05 > 0) {
  cat(sprintf("  • %d genes show suggestive differential rhythmicity (p < 0.05)\n", dr_p05))
  cat(sprintf("  • No genes pass FDR correction - weak evidence for age effects\n"))
} else {
  cat(sprintf("  • No evidence for age-related differences in circadian rhythms\n"))
  cat(sprintf("  • Core circadian machinery appears stable with age\n"))
}

# Save results
result$age_median <- age_median
result$young_n <- length(young_idx)
result$old_n <- length(old_idx)
result$young_age_range <- c(min(ctrl_clinical$Age[young_idx]), max(ctrl_clinical$Age[young_idx]))
result$old_age_range <- c(min(ctrl_clinical$Age[old_idx]), max(ctrl_clinical$Age[old_idx]))

saveRDS(result, "output/putamen_age_dcp_results.rds")
cat("\nResults saved to: output/putamen_age_dcp_results.rds\n")

cat("\n=== Analysis Complete ===\n")
