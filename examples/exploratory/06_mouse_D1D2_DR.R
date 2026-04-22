#' =======================================================================
#' 06_mouse_D1D2_DR.R — Mouse D1 vs D2 Cell Type DR Power Analysis
#' =======================================================================
#'
#' PURPOSE
#'   Explore differential rhythmicity (DR) power for mouse D1 vs D2 cell
#'   types (likely striatal dopamine receptor-expressing neurons). This
#'   dataset has a realistic multi-replicate active design:
#'     - B = 6 ZT time points (2, 6, 10, 14, 18, 22h), every 4h
#'     - m ≈ 7-8 replicates per time point per group
#'     - N ≈ 45 subjects per group (total pilot n=91)
#'
#'   Effect size (r ≈ 0.65) is moderate — between baboon (r ≈ 1.7) and
#'   human aging (r ≈ 0.56). Unique because this dataset has true
#'   biological replicates at each time point, making the B vs m tradeoff
#'   directly relevant.
#'
#' DATA
#'   data/mouse_clinicalinfo_03082021_rmOutliers.csv  — sample metadata
#'   data/mouse_D1D2_logCPMfiltered_counts.csv        — raw counts (ENSMUSG)
#'   Note: expression values are raw counts (not log-transformed despite
#'   the filename); apply log2(CPM+1) normalization before fitting.
#'
#' DESIGN
#'   - Active, controlled ZT (6 time points: 2,6,10,14,18,22h)
#'   - Groups: D1 cell type vs D2 cell type
#'   - n per group: D1=45, D2=46
#'   - Expected n80 for DR: moderate (~40–80 based on r≈0.65, prop_DR≈22%)
#'
#' OUTPUTS
#'   output/06_mouse_D1D2_DR_<timestamp>/
#'     s0_signal_summary.txt              (prop_rhy, r, prop_DR)
#'     s1_d1d2_boot_grid.rds
#'     figures/s1_d1d2_bootstrap.pdf
#'     figures/s2_bm_tradeoff.pdf
#'
#' USAGE
#'   Rscript examples/exploratory/06_mouse_D1D2_DR.R
#'   POWERSIM_SMOKE=1 Rscript examples/exploratory/06_mouse_D1D2_DR.R
#'
#' @author Thien Pham

# =====================================================================
# SETTINGS
# =====================================================================

SMOKE       <- nzchar(Sys.getenv("POWERSIM_SMOKE"))

NGENES      <- if (SMOKE) 500L  else 5000L
NBOOT       <- if (SMOKE) 3L    else 100L
NSIMS_INNER <- if (SMOKE) 3L    else 20L

# B sweep: 3 or 6 ZT time points
# B=3 (every 8h: e.g., ZT2, ZT10, ZT18) vs B=6 (full coverage every 4h)
B_VALS      <- c(3L, 6L)

# N_values: all must be divisible by max(B_VALS)=6
# Smoke: small grid; Production: wider grid
N_GRID      <- if (SMOKE) c(24L, 36L, 48L)  else c(24L, 36L, 48L, 60L, 72L, 96L, 120L)

RHYTHM_PVAL <- 0.05

DATA_PHENO  <- "data/mouse_clinicalinfo_03082021_rmOutliers.csv"   # local data/
DATA_EXPR   <- "data/mouse_D1D2_logCPMfiltered_counts.csv"         # local data/

cat(sprintf("Mode: %s | NGENES=%d | NBOOT=%d | NSIMS_INNER=%d\n",
            if (SMOKE) "SMOKE" else "PRODUCTION",
            NGENES, NBOOT, NSIMS_INNER))

# =====================================================================
# PATH CONFIGURATION
# =====================================================================
# Set POWERSIM_ROOT as an environment variable for portability:
#   Linux/Mac shell:  export POWERSIM_ROOT=/path/to/PowerSim
#   R console:        Sys.setenv(POWERSIM_ROOT="/path/to/PowerSim")

POWERSIM_ROOT <- {
  env <- Sys.getenv("POWERSIM_ROOT", unset = "")
  if (nzchar(env)) env else
    "/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim"
}

# =====================================================================
# SETUP
# =====================================================================

POWERSIM_DIR <- POWERSIM_ROOT
setwd(POWERSIM_DIR)
source_dir <- file.path(getwd(), "code")
old_wd <- setwd(source_dir)
source("setup.R")
setwd(old_wd)

out_dir <- file.path("output", sprintf("06_mouse_D1D2_DR_%s", format(Sys.time(), "%Y%m%d_%H%M")))
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("Output: %s/\n\n", out_dir))

opts_analysis <- CircadianAnalysisOptions(
  alpha           = 0.05,
  p.adjust.method = "BH",
  fdr_thresholds  = c(0.01, 0.05, 0.10),
  reference_n     = 60
)

# =====================================================================
# SECTION 0: LOAD AND CHARACTERIZE DATA
# =====================================================================

cat("====================================================================\n")
cat("SECTION 0: Load and characterize D1/D2 mouse data\n")
cat("====================================================================\n\n")

pheno <- read.csv(DATA_PHENO, row.names = 1)
cat(sprintf("Pheno: %d samples, columns: %s\n", nrow(pheno),
            paste(names(pheno), collapse=", ")))

# Raw counts (despite filename); prepCircadianData normalizes to log2(CPM+1)
# and aligns sample order via pheno$sample column
prep    <- prepCircadianData(DATA_EXPR, times = "time", input_type = "counts",
                             pheno = pheno, sample_col = "sample")
log_mat <- prep$data
cat(sprintf("log2(CPM+1) range: %.2f to %.2f  median: %.2f\n",
            min(log_mat), max(log_mat), median(log_mat)))

# Split by cell type
d1_samp <- pheno$sample[pheno$cell == "D1"]
d2_samp <- pheno$sample[pheno$cell == "D2"]

mat_d1 <- log_mat[, colnames(log_mat) %in% d1_samp, drop=FALSE]
mat_d2 <- log_mat[, colnames(log_mat) %in% d2_samp, drop=FALSE]

tod_d1 <- pheno$time[match(colnames(mat_d1), pheno$sample)]
tod_d2 <- pheno$time[match(colnames(mat_d2), pheno$sample)]

cat(sprintf("\nD1: n=%d, ZT=%s\n", ncol(mat_d1), paste(sort(unique(tod_d1)), collapse=",")))
cat(sprintf("D2: n=%d, ZT=%s\n", ncol(mat_d2), paste(sort(unique(tod_d2)), collapse=",")))
cat("Replicates per ZT (D1):\n")
print(table(ZT=tod_d1))
cat("Replicates per ZT (D2):\n")
print(table(ZT=tod_d2))

# Filter low-expression genes (keep genes expressed in ≥ half of D1 samples)
keep <- rowSums(mat_d1 > 1) >= floor(ncol(mat_d1) / 2)
mat_d1 <- mat_d1[keep, , drop=FALSE]
mat_d2 <- mat_d2[keep, , drop=FALSE]
cat(sprintf("\nGenes after expression filter (logCPM>1 in ≥50%% of D1): %d\n", nrow(mat_d1)))

# Subsample for cosinor estimation
set.seed(42)
g_sub <- sample(nrow(mat_d1), min(NGENES, nrow(mat_d1)))
mat_d1_s <- mat_d1[g_sub, , drop=FALSE]
mat_d2_s <- mat_d2[g_sub, , drop=FALSE]

# Cosinor fits
fit_d1 <- fitCosinorAll(mat_d1_s, times=tod_d1, period=24)
fit_d2 <- fitCosinorAll(mat_d2_s, times=tod_d2, period=24)

rhy_d1 <- !is.na(fit_d1$pvalue) & fit_d1$pvalue < RHYTHM_PVAL
rhy_d2 <- !is.na(fit_d2$pvalue) & fit_d2$pvalue < RHYTHM_PVAL
dr_mask <- xor(rhy_d1, rhy_d2)

r_d1 <- as.numeric(fit_d1$A[rhy_d1]) / as.numeric(fit_d1$sigma[rhy_d1])
r_d2 <- as.numeric(fit_d2$A[rhy_d2]) / as.numeric(fit_d2$sigma[rhy_d2])

prop_DR_d1d2 <- mean(dr_mask)
r_d1_med     <- median(r_d1, na.rm=TRUE)

cat(sprintf("\n=== D1 vs D2 Circadian Signal (p < %.2f) ===\n", RHYTHM_PVAL))
cat(sprintf("D1: prop_rhythmic=%.1f%%  r_median=%.2f\n", mean(rhy_d1)*100, r_d1_med))
cat(sprintf("D2: prop_rhythmic=%.1f%%  r_median=%.2f\n", mean(rhy_d2)*100,
            median(r_d2, na.rm=TRUE)))
cat(sprintf("prop_DR (D1 vs D2): %.1f%%\n", prop_DR_d1d2*100))

# Save signal summary
s0_lines <- c(
  "Mouse D1 vs D2 cell type — circadian signal summary",
  sprintf("Data: %d genes (log2 CPM+1, filtered) x D1(n=%d) D2(n=%d)", nrow(mat_d1), ncol(mat_d1), ncol(mat_d2)),
  sprintf("Design: active, B=6 ZT (2,6,10,14,18,22h), ~3-4 repl/ZT/group"),
  sprintf("D1: prop_rhythmic=%.1f%%  r_median=%.2f  (p<%.2f)", mean(rhy_d1)*100, r_d1_med, RHYTHM_PVAL),
  sprintf("D2: prop_rhythmic=%.1f%%  r_median=%.2f", mean(rhy_d2)*100, median(r_d2,na.rm=TRUE)),
  sprintf("prop_DR (D1 vs D2): %.1f%%", prop_DR_d1d2*100),
  "",
  "Comparison with other datasets:",
  "  Mouse GSE54651 LIV vs CER: r~2.9, prop_DR~25%, n80 < 24",
  "  Baboon LUN vs CER:         r~1.7, prop_DR~41%, n80 ~24-36",
  "  Mouse D1 vs D2 (this):     r~0.65, prop_DR~22%, n80 TBD",
  "  Human aging young vs old:  r~0.5, prop_DR~14%, n80 > 300"
)
writeLines(s0_lines, file.path(out_dir, "s0_signal_summary.txt"))
cat(sprintf("\nSaved: %s/s0_signal_summary.txt\n", out_dir))

# =====================================================================
# SECTION 1: BOOTSTRAP DESIGN GRID (B vs m tradeoff)
# =====================================================================

cat("\n====================================================================\n")
cat("SECTION 1: Bootstrap design grid — D1 vs D2 (sweep B and N)\n")
cat("====================================================================\n\n")

cat(sprintf("B_values: %s\n", paste(B_VALS, collapse=", ")))
cat(sprintf("N_values: %s\n", paste(N_GRID, collapse=", ")))
cat(sprintf("  (B=3: ZT every 8h; B=6: full 4h coverage)\n\n"))

# ZT time points for the design vector (full B=6 coverage = pilot TOD)
d1_zt_full <- sort(unique(tod_d1))   # 2, 6, 10, 14, 18, 22
boot_opts <- CircadianBootstrapOptions(
  design        = "active",
  design_vector = d1_zt_full,
  B_values      = B_VALS,
  N_values      = N_GRID,
  nboot         = NBOOT,
  nsims_inner   = NSIMS_INNER,
  seed          = 42
)

bio_d1d2 <- estCircadianParam(
  data            = mat_d1_s,
  times           = tod_d1,
  period          = 24,
  prop_DR         = prop_DR_d1d2,
  prop_DP         = 0,
  min_rhythm_pval = RHYTHM_PVAL,
  verbose         = FALSE
)

cat("Running bootstrap design grid (D1 vs D2)...\n")
s1_boot <- runBootstrapDesignGrid(
  pilot_data    = mat_d1_s,
  pilot_times   = tod_d1,
  boot.opts     = boot_opts,
  analysis.opts = opts_analysis,
  bio_diff.opts = bio_d1d2,
  verbose       = TRUE
)
cat("\n--- Bootstrap Summary (D1 vs D2 DR, FDR 5%) ---\n")
summaryBootstrapDesignGrid(s1_boot, test_type = "DR", fdr_threshold = 0.05)

plotBootstrapDesignGrid(s1_boot, test_type = "DR", fdr_threshold = 0.05,
                        output_file = file.path(fig_dir, "s1_d1d2_bootstrap.pdf"))
saveRDS(s1_boot, file.path(out_dir, "s1_d1d2_boot_grid.rds"))
cat(sprintf("Figure: %s/figures/s1_d1d2_bootstrap.pdf\n", out_dir))

# =====================================================================
# SECTION 2: B vs m TRADEOFF PLOT
# =====================================================================

cat("\n====================================================================\n")
cat("SECTION 2: B vs m tradeoff — D1 vs D2\n")
cat("====================================================================\n\n")

fig_bm <- file.path(fig_dir, "s2_bm_tradeoff.pdf")
pdf(fig_bm, width=7, height=5)
par(mar=c(4,4,3,1))

fdr_idx <- which.min(abs(s1_boot$fdr_thresholds - 0.05))
n_N  <- length(s1_boot$N_values)
n_B  <- length(s1_boot$B_values)

# power_mean is [n_N x n_B x n_tests]; extract DR (fdr_idx)
pm <- matrix(s1_boot$power_mean[, , fdr_idx], nrow=n_N, ncol=n_B)
lo <- matrix(s1_boot$power_ci_lo[, , fdr_idx], nrow=n_N, ncol=n_B)
hi <- matrix(s1_boot$power_ci_hi[, , fdr_idx], nrow=n_N, ncol=n_B)

cols_b <- c("steelblue", "darkorange")
if (n_B > 2) cols_b <- rainbow(n_B, s=0.7, v=0.85)

plot(NA, xlim=range(s1_boot$N_values), ylim=c(0,1),
     xlab="N per group", ylab="Power (FDR 5%)",
     main=sprintf("Mouse D1 vs D2: B vs m tradeoff\nprop_DR=%.0f%%, r~%.2f",
                  prop_DR_d1d2*100, r_d1_med))
abline(h=0.80, lty=2, col="gray40")
abline(h=c(0.2,0.4,0.6), lty=3, col="gray85")

for (b in seq_len(n_B)) {
  Bv <- s1_boot$B_values[b]
  col_b <- cols_b[b]
  polygon(c(s1_boot$N_values, rev(s1_boot$N_values)),
          c(lo[,b], rev(hi[,b])),
          col=adjustcolor(col_b, 0.15), border=NA)
  lines(s1_boot$N_values, pm[,b], col=col_b, lwd=2, lty=b)
  points(s1_boot$N_values, pm[,b], col=col_b, pch=15+b, cex=0.9)
}

# Approximate m for each B at minimum N
m_labels <- sprintf("B=%d (m~%d repl/ZT)", s1_boot$B_values,
                    pmax(1L, round(min(N_GRID)/s1_boot$B_values)))
legend("bottomright", legend=m_labels,
       col=cols_b, lwd=2, lty=seq_len(n_B), pch=15+seq_len(n_B),
       bty="n", cex=0.85)

dev.off()
cat(sprintf("Figure: %s\n", fig_bm))

# =====================================================================
# WRAP-UP
# =====================================================================

# Extract n80 (bootstrap median)
fdr_idx <- which.min(abs(s1_boot$fdr_thresholds - 0.05))
pm_margin <- rowMeans(matrix(s1_boot$power_mean[,,fdr_idx], nrow=n_N, ncol=n_B), na.rm=TRUE)
lo_margin <- rowMeans(matrix(s1_boot$power_ci_lo[,,fdr_idx], nrow=n_N, ncol=n_B), na.rm=TRUE)
hi_margin <- rowMeans(matrix(s1_boot$power_ci_hi[,,fdr_idx], nrow=n_N, ncol=n_B), na.rm=TRUE)

n80_med <- N_GRID[which(pm_margin >= 0.80)[1]]
n80_lo  <- N_GRID[which(hi_margin >= 0.80)[1]]
n80_hi  <- N_GRID[which(lo_margin >= 0.80)[1]]

cat("\n====================================================================\n")
cat("06_mouse_D1D2_DR COMPLETE\n")
cat("====================================================================\n\n")

cat("--- Summary ---\n")
cat(sprintf("D1 vs D2: prop_DR=%.1f%%  r_median(D1)=%.2f\n", prop_DR_d1d2*100, r_d1_med))
cat(sprintf("n80 (bootstrap median): %s\n", ifelse(is.na(n80_med), ">max(N)", n80_med)))
cat(sprintf("n80 95%% CI: [%s, %s]\n",
            ifelse(is.na(n80_lo), ">max(N)", n80_lo),
            ifelse(is.na(n80_hi), ">max(N)", n80_hi)))
cat(sprintf("\nOutput: %s/\n", out_dir))
