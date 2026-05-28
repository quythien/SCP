#' Quick test of LRTest_diff_phase function

setwd("/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")

source("code/options.R")
source("code/simulation.R")
source("code/LRTest_diff_phase.R")
source("code/detection.R")

cat("Testing LRTest_diff_phase function...\n\n")

# Create simple test data with same phase (null case)
set.seed(123)
n <- 50
tt1 <- runif(n, 0, 24)
tt2 <- runif(n, 0, 24)

# Same phase - should have large p-value
Amp1 <- 2
Phase1 <- 6
Offset1 <- 3
yy1 <- Amp1 * sin(2*pi/24 * (tt1 + Phase1)) + Offset1 + rnorm(n, 0, 0.5)

# Same parameters
yy2 <- Amp1 * sin(2*pi/24 * (tt2 + Phase1)) + Offset1 + rnorm(n, 0, 0.5)

cat("Test 1: Same phase (NULL case)\n")
result1 <- LRTest_diff_phase(tt1, yy1, tt2, yy2, period = 24, FN = TRUE)
cat(sprintf("  p-value: %.6f\n", result1$pvalue))
cat(sprintf("  phase1: %.4f, phase2: %.4f\n", result1$phase_1, result1$phase_2))
cat(sprintf("  Expected: p-value should be LARGE (close to 1)\n\n"))

# Different phase - should have small p-value
Phase2 <- 12  # 6 hour difference
yy3 <- Amp1 * sin(2*pi/24 * (tt2 + Phase2)) + Offset1 + rnorm(n, 0, 0.5)

cat("Test 2: Different phase (6 hour shift)\n")
result2 <- LRTest_diff_phase(tt1, yy1, tt2, yy3, period = 24, FN = TRUE)
cat(sprintf("  p-value: %.6e\n", result2$pvalue))
cat(sprintf("  phase1: %.4f, phase2: %.4f\n", result2$phase_1, result2$phase_2))
cat(sprintf("  Expected: p-value should be SMALL\n\n"))

# Test LR_diff wrapper
cat("Test 3: LR_diff wrapper function\n")
result3 <- LR_diff(tt1, yy1, tt2, yy3, 24, FN = TRUE, type = "phase")
cat(sprintf("  p-value: %.6e\n", result3$pvalue))
cat(sprintf("  delta_phase: %.4f\n", result3$delta_phase))
