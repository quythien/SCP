#' =======================================================================
#' dr_camo_power.R — DR Power Study: Cross-Species and Cross-Tissue
#' =======================================================================
#'
#' PURPOSE
#'   Compare differential rhythmicity (DR) power across three settings:
#'     1. BA11 vs BA47 (age groups, same brain tissue)   [reference]
#'     2. Baboon LUN vs GTEx Human LUN (cross-species, same tissue)
#'     3. Baboon LUN vs Baboon ILE (cross-tissue, same species)
#'
#'   Hypothesis: DR effect size is larger in cross-tissue (LUN vs ILE)
#'   than in age-group comparison (BA11 vs BA47), because tissue-specific
#'   circadian programs differ more than age-related changes.
#'
#' DATA
#'   CAMO.bab.hum.RData:
#'     - gtex$CPM.large.clean: Human GTEx, 26 tissues, log-scaled CPM
#'     - baboon_withTOD$baboon: Baboon expression, 26 tissues, raw CPM
#'     - baboon_withTOD$tod: ZT times (0,2,4,...,22h, n=12 per tissue)
#'     - gtex$tod: Real post-mortem TOD (continuous 0-24h)
#'     - All datasets share 5066 human-ortholog ENSEMBL gene IDs
#'
#' OUTPUTS
#'   output/dr_camo_<timestamp>/
#'     figures/dr_power_comparison.pdf   (3-way DR power curves)
#'     figures/rhythmicity_summary.pdf   (prop rhythmic per group)
#'     dr_camo_summary.txt
#'     results.rds
#'
#' USAGE
#'   Rscript examples/exploratory/dr_camo_power.R
#'   POWERSIM_SMOKE=1 Rscript examples/exploratory/dr_camo_power.R
#'
#' @author Thien Pham

# =====================================================================
# SETTINGS
# =====================================================================

SMOKE <- nzchar(Sys.getenv("POWERSIM_SMOKE"))

NGENES   <- if (SMOKE) 500L  else 5000L
NSIMS    <- if (SMOKE) 5L    else 50L
N_GRID   <- if (SMOKE) c(20L, 40L, 60L) else c(20L, 40L, 60L, 80L, 100L, 120L, 160L)

# Threshold for classifying rhythmic genes when estimating prop_DR.
# Using raw p < 0.01 (not BH-adjusted) since the goal is to estimate
# the biological proportion, not to control FDR.
# BH FDR is too conservative for small gene sets (e.g., 500 subsampled).
RHYTHM_PVAL <- 0.01   # raw p-value threshold

CAMO_PATH <- "/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Collaborative/GTEXdata/data/CAMO.bab.hum.RData"

# =====================================================================
# SETUP
# =====================================================================

setwd("/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")
old_wd <- setwd(file.path(getwd(), "code"))
source("setup.R")
setwd(old_wd)

run_tag  <- format(Sys.time(), "%Y%m%d_%H%M")
base_out <- file.path("output", paste0("dr_camo_", run_tag))
fig_dir  <- file.path(base_out, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

cat("=======================================================================\n")
cat("DR Power Study: Cross-Species and Cross-Tissue Comparison\n")
cat("=======================================================================\n\n")
cat(sprintf("Settings: ngenes=%d, nsims=%d, N grid = %s\n",
            NGENES, NSIMS, paste(N_GRID, collapse=", ")))
cat(sprintf("Output: %s/\n\n", base_out))

t_start <- proc.time()

# =====================================================================
# LOAD CAMO DATA
# =====================================================================

cat("Loading CAMO data...\n")
load(CAMO_PATH)

# --- GTEx LUN (human lung, log-scaled CPM, n=578) ---
gtex_lun_mat  <- as.matrix(gtex$CPM.large.clean[["LUN"]])
class(gtex_lun_mat) <- "numeric"
gtex_lun_tod  <- gtex$tod[["LUN"]]
cat(sprintf("  GTEx LUN:    %d genes x %d samples | TOD range [%.1f, %.1f]\n",
            nrow(gtex_lun_mat), ncol(gtex_lun_mat),
            min(gtex_lun_tod), max(gtex_lun_tod)))

# --- Baboon LUN (raw CPM, n=12, ZT 0-22h) ---
bab_lun_df   <- baboon_withTOD$baboon[["LUN"]]
bab_lun_cols <- grep("^LUN\\.ZT", colnames(bab_lun_df))
bab_lun_mat  <- log2(as.matrix(bab_lun_df[, bab_lun_cols]) + 1)
rownames(bab_lun_mat) <- rownames(bab_lun_df)
bab_lun_tod  <- baboon_withTOD$tod[["LUN"]]
cat(sprintf("  Baboon LUN:  %d genes x %d samples | ZT [%.0f, %.0f]h\n",
            nrow(bab_lun_mat), ncol(bab_lun_mat),
            min(bab_lun_tod), max(bab_lun_tod)))

# --- Baboon ILE (raw CPM, n=12, ZT 0-22h) ---
bab_ile_df   <- baboon_withTOD$baboon[["ILE"]]
bab_ile_cols <- grep("^ILE\\.ZT", colnames(bab_ile_df))
bab_ile_mat  <- log2(as.matrix(bab_ile_df[, bab_ile_cols]) + 1)
rownames(bab_ile_mat) <- rownames(bab_ile_df)
bab_ile_tod  <- baboon_withTOD$tod[["ILE"]]
cat(sprintf("  Baboon ILE:  %d genes x %d samples | ZT [%.0f, %.0f]h\n\n",
            nrow(bab_ile_mat), ncol(bab_ile_mat),
            min(bab_ile_tod), max(bab_ile_tod)))

# Align genes across all three (should already match, but be safe)
common_genes <- Reduce(intersect, list(rownames(gtex_lun_mat),
                                       rownames(bab_lun_mat),
                                       rownames(bab_ile_mat)))
cat(sprintf("Common genes across all 3 groups: %d\n\n", length(common_genes)))
gtex_lun_mat <- gtex_lun_mat[common_genes, ]
bab_lun_mat  <- bab_lun_mat[common_genes, ]
bab_ile_mat  <- bab_ile_mat[common_genes, ]

# Subsample genes for speed if needed
if (NGENES < nrow(gtex_lun_mat)) {
  set.seed(42)
  gene_idx     <- sample(nrow(gtex_lun_mat), NGENES)
  gtex_lun_mat <- gtex_lun_mat[gene_idx, ]
  bab_lun_mat  <- bab_lun_mat[gene_idx, ]
  bab_ile_mat  <- bab_ile_mat[gene_idx, ]
  cat(sprintf("Subsampled to %d genes for speed.\n\n", NGENES))
}

# =====================================================================
# FIT COSINOR TO EACH GROUP — ESTIMATE RHYTHMICITY
# =====================================================================

cat(sprintf("Fitting cosinor to each group (raw p < %.2f)...\n", RHYTHM_PVAL))

.fit_group <- function(expr_mat, times, label) {
  result <- fitCosinorAll(expr_mat, times)
  n_rhythmic <- sum(result$pvalue < RHYTHM_PVAL, na.rm = TRUE)
  prop <- n_rhythmic / nrow(result)
  cat(sprintf("  %-20s: %d / %d rhythmic = %.1f%%  (raw p < %.2f)\n",
              label, n_rhythmic, nrow(result), 100 * prop, RHYTHM_PVAL))
  list(fit = result,
       rhythmic = result$pvalue < RHYTHM_PVAL,
       prop_rhythmic = prop)
}

fit_gtex_lun <- .fit_group(gtex_lun_mat, gtex_lun_tod, "GTEx LUN (human)")
fit_bab_lun  <- .fit_group(bab_lun_mat,  bab_lun_tod,  "Baboon LUN")
fit_bab_ile  <- .fit_group(bab_ile_mat,  bab_ile_tod,  "Baboon ILE")

# =====================================================================
# ESTIMATE EMPIRICAL DR PROPORTIONS
# =====================================================================

cat("\nEstimating empirical DR proportions...\n")

.compute_DR_prop <- function(rhythmic1, rhythmic2, label1, label2) {
  r1 <- rhythmic1 & !rhythmic2   # rhythmic in 1 only
  r2 <- rhythmic2 & !rhythmic1   # rhythmic in 2 only
  both <- rhythmic1 & rhythmic2
  neither <- !rhythmic1 & !rhythmic2
  n <- length(rhythmic1)
  cat(sprintf("  %s vs %s:\n", label1, label2))
  cat(sprintf("    Rhythmic in %s only:  %d (%.1f%%)\n", label1, sum(r1), 100*mean(r1)))
  cat(sprintf("    Rhythmic in %s only: %d (%.1f%%)\n", label2, sum(r2), 100*mean(r2)))
  cat(sprintf("    Rhythmic in both:     %d (%.1f%%)\n", sum(both), 100*mean(both)))
  cat(sprintf("    DR total (either):    %d (%.1f%%)\n", sum(r1|r2), 100*mean(r1|r2)))
  mean(r1 | r2)  # prop_DR = fraction with discordant rhythmicity
}

prop_DR_cross_species <- .compute_DR_prop(
  fit_bab_lun$rhythmic, fit_gtex_lun$rhythmic, "Baboon LUN", "GTEx LUN")
prop_DR_cross_tissue  <- .compute_DR_prop(
  fit_bab_lun$rhythmic, fit_bab_ile$rhythmic,  "Baboon LUN", "Baboon ILE")

cat(sprintf("\nSummary:\n"))
cat(sprintf("  Cross-species prop_DR (Bab LUN vs GTEx LUN): %.3f (%.1f%%)\n",
            prop_DR_cross_species, 100 * prop_DR_cross_species))
cat(sprintf("  Cross-tissue  prop_DR (Bab LUN vs Bab ILE):  %.3f (%.1f%%)\n\n",
            prop_DR_cross_tissue, 100 * prop_DR_cross_tissue))

# =====================================================================
# PARAMETER ESTIMATION VIA estCircadianParam
# =====================================================================

cat("Estimating circadian parameters for power simulation...\n\n")

opts_analysis <- CircadianAnalysisOptions(reference_n = 60L)

# Cross-species: use GTEx LUN as pilot (larger, more reliable)
cat("--- Cross-species pilot: GTEx LUN ---\n")
opts_bio_cross_species <- estCircadianParam(
  data         = gtex_lun_mat,
  times        = gtex_lun_tod,
  prop_DR      = prop_DR_cross_species,
  prop_DP      = 0.00,
  phase_diff   = c(0, 0),
  amp_diff     = c(1, 1),
  verbose      = TRUE
)

# Cross-tissue: use baboon LUN as pilot (clean, even spacing)
cat("\n--- Cross-tissue pilot: Baboon LUN ---\n")
opts_bio_cross_tissue <- estCircadianParam(
  data         = bab_lun_mat,
  times        = bab_lun_tod,
  prop_DR      = prop_DR_cross_tissue,
  prop_DP      = 0.00,
  phase_diff   = c(0, 0),
  amp_diff     = c(1, 1),
  verbose      = TRUE
)

# =====================================================================
# DR POWER ANALYSIS — CROSS-SPECIES (GTEx LUN passive design)
# =====================================================================

cat("\n====================================================================\n")
cat("DR POWER: Cross-Species (Baboon LUN vs GTEx Human LUN)\n")
cat("====================================================================\n")
cat("Design: passive (real GTEx TOD distribution)\n\n")

opts_design_cs <- CircadianDesignOptions(
  sample_sizes = N_GRID,
  nsims        = NSIMS,
  design       = "passive",
  cts          = gtex_lun_tod,
  test_types   = "DR"
)

dr_cross_species <- runPowerAnalysis(opts_bio_cross_species, opts_design_cs, opts_analysis,
                                     test_type = "DR")
saveRDS(dr_cross_species, file.path(base_out, "dr_cross_species.rds"))

cs_summary <- summaryRunPower(dr_cross_species)
cat("\nCross-species DR power (FDR 5%):\n")
print(cs_summary[, c("n", "Power", "TypeI_Error", "FDR")])

# =====================================================================
# DR POWER ANALYSIS — CROSS-TISSUE (Baboon LUN vs ILE, active design)
# =====================================================================

cat("\n====================================================================\n")
cat("DR POWER: Cross-Tissue (Baboon LUN vs Baboon ILE)\n")
cat("====================================================================\n")
cat("Design: active (ZT 0-22h, even spacing)\n\n")

opts_design_ct <- CircadianDesignOptions(
  sample_sizes = N_GRID,
  nsims        = NSIMS,
  design       = "active",
  test_types   = "DR"
)

dr_cross_tissue <- runPowerAnalysis(opts_bio_cross_tissue, opts_design_ct, opts_analysis,
                                    test_type = "DR")
saveRDS(dr_cross_tissue, file.path(base_out, "dr_cross_tissue.rds"))

ct_summary <- summaryRunPower(dr_cross_tissue)
cat("\nCross-tissue DR power (FDR 5%):\n")
print(ct_summary[, c("n", "Power", "TypeI_Error", "FDR")])

# =====================================================================
# COMPARISON FIGURE: 3-WAY DR POWER CURVES
# =====================================================================

cat("\nGenerating comparison figure...\n")

pdf(file.path(fig_dir, "dr_power_comparison.pdf"), width = 8, height = 5)
par(mar = c(4.5, 4.5, 3, 8), xpd = TRUE)

cols <- c("steelblue", "darkorange", "forestgreen")
ltys <- c(1, 2, 3)
pwrs <- list(cs_summary, ct_summary)
labels <- c(
  sprintf("Bab LUN vs GTEx LUN (prop_DR=%.0f%%)", 100 * prop_DR_cross_species),
  sprintf("Bab LUN vs Bab ILE (prop_DR=%.0f%%)",  100 * prop_DR_cross_tissue)
)

plot(N_GRID, rep(NA, length(N_GRID)), ylim = c(0, 1),
     xlab = "Sample size per group (N)",
     ylab = "DR Power (FDR 5%)",
     main = "DR Power: Cross-Species vs Cross-Tissue",
     type = "n", las = 1)
abline(h = 0.80, lty = 2, col = "grey50")
abline(h = seq(0, 1, 0.2), lty = 3, col = "grey90")
text(max(N_GRID) * 0.98, 0.82, "80%", col = "grey50", cex = 0.8, adj = 1)

for (i in seq_along(pwrs)) {
  pwr <- pwrs[[i]]
  ns  <- pwr$n
  pow <- pwr$Power
  lines(ns, pow, col = cols[i], lty = ltys[i], lwd = 2)
  points(ns, pow, col = cols[i], pch = 16, cex = 1.2)
}

legend("right", inset = c(-0.42, 0),
       legend = labels,
       col = cols[1:2], lty = ltys[1:2], lwd = 2, pch = 16,
       bty = "n", cex = 0.8)
dev.off()
cat(sprintf("Figure: %s\n", file.path(fig_dir, "dr_power_comparison.pdf")))

# =====================================================================
# RHYTHMICITY SUMMARY FIGURE
# =====================================================================

pdf(file.path(fig_dir, "rhythmicity_summary.pdf"), width = 7, height = 5)
par(mar = c(5, 5, 3, 2))

groups <- c("GTEx LUN\n(human)", "Baboon LUN", "Baboon ILE")
props  <- c(fit_gtex_lun$prop_rhythmic,
            fit_bab_lun$prop_rhythmic,
            fit_bab_ile$prop_rhythmic)
bp <- barplot(props * 100,
              names.arg = groups,
              col = c("steelblue", "darkorange", "forestgreen"),
              ylab = "% Rhythmic genes (BH FDR 5%)",
              main = "Rhythmic Gene Proportion by Group",
              ylim = c(0, max(props * 100) * 1.3),
              las = 1)
text(bp, props * 100 + 1, sprintf("%.1f%%", props * 100), cex = 0.9, font = 2)
dev.off()
cat(sprintf("Figure: %s\n", file.path(fig_dir, "rhythmicity_summary.pdf")))

# =====================================================================
# WRAP-UP
# =====================================================================

t_elapsed <- proc.time() - t_start

summary_lines <- c(
  "DR Power Study: Cross-Species and Cross-Tissue Comparison",
  sprintf("Run: %s", run_tag),
  sprintf("Settings: ngenes=%d, nsims=%d, N grid=%s",
          NGENES, NSIMS, paste(N_GRID, collapse=",")),
  "",
  "=== Rhythmicity ===",
  sprintf("GTEx LUN (human):   prop_rhythmic = %.1f%%", 100 * fit_gtex_lun$prop_rhythmic),
  sprintf("Baboon LUN:         prop_rhythmic = %.1f%%", 100 * fit_bab_lun$prop_rhythmic),
  sprintf("Baboon ILE:         prop_rhythmic = %.1f%%", 100 * fit_bab_ile$prop_rhythmic),
  "",
  "=== Empirical DR Proportions ===",
  sprintf("Cross-species (Bab LUN vs GTEx LUN): prop_DR = %.1f%%",
          100 * prop_DR_cross_species),
  sprintf("Cross-tissue  (Bab LUN vs Bab ILE):  prop_DR = %.1f%%",
          100 * prop_DR_cross_tissue),
  "",
  "=== DR Power at FDR 5% ===",
  "Cross-species:",
  paste(sprintf("  n=%d: %.1f%%", cs_summary$n, 100 * cs_summary$Power), collapse="\n"),
  "Cross-tissue:",
  paste(sprintf("  n=%d: %.1f%%", ct_summary$n, 100 * ct_summary$Power), collapse="\n"),
  "",
  sprintf("Runtime: %.1f min", t_elapsed[3] / 60),
  sprintf("Output:  %s/", base_out)
)

writeLines(summary_lines, file.path(base_out, "dr_camo_summary.txt"))

cat("\n=======================================================================\n")
cat("COMPLETE\n")
cat("=======================================================================\n\n")
writeLines(summary_lines)
cat("\n")
