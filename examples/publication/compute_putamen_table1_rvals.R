#' Compute "Median r (7 core clock)" and "Median r (top 100 rhythmic)"
#' for the Ketchesin GSE160521 Putamen control cohort (n=59), matching
#' the metric used elsewhere in Table 1.
setwd("/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")

obs <- read.csv("/home/qtp1/Projects/Circadian/Kyle/Kyle_multiBrainRegion/Putamen_observed_para_control_rm97_symbols.csv",
                row.names = 1, stringsAsFactors = FALSE)
n  <- 59  # control sample size

# Intrinsic effect size r = A/sigma where sigma is the per-gene
# residual noise (independent of n). In the cosinor regression with
# uniform-on-period sampling, Var(y) = A^2/2 + sigma^2 and
# R^2 = (A^2/2)/(A^2/2 + sigma^2) = r^2/(r^2 + 2), so
# r = sqrt(2 R^2 / (1 - R^2)).
obs$r_inferred <- with(obs, sqrt(2 * R2 / pmax(1 - R2, 1e-9)))

# --- 7 core clock genes (Zong 2023 / CircaPower core panel) ---
core_clock <- c("ARNTL", "CLOCK", "NR1D1", "NR1D2", "PER1", "PER2", "PER3",
                "CRY1", "CRY2", "DBP")
present    <- intersect(core_clock, obs$symbols)
cat("Core clock genes available:\n"); print(present)

core_rows  <- obs[obs$symbols %in% present, ]
cat("\nCore clock r values:\n")
print(core_rows[, c("symbols","A","R2","pvalue","r_inferred")])

# Pick the canonical 7 used in Zong 2023 if all present; otherwise top by R^2
pick7 <- head(core_clock[core_clock %in% present], 7)
r_7core <- median(obs$r_inferred[obs$symbols %in% pick7])
cat(sprintf("\nMedian r over %d core clock genes: %.3f\n", length(pick7), r_7core))

# --- Top 100 rhythmic by p-value ---
top100  <- obs[order(obs$pvalue), ][seq_len(100), ]
r_top100 <- median(top100$r_inferred)
cat(sprintf("Median r over top 100 rhythmic: %.3f\n", r_top100))
cat(sprintf("Median r over top 300 rhythmic: %.3f\n",
            median(obs[order(obs$pvalue), ][seq_len(300), "r_inferred"])))

cat(sprintf("\n=> Table 1 row: Putamen   59   %.2f   %.2f\n", r_7core, r_top100))
