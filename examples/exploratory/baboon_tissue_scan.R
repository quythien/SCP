#' =======================================================================
#' baboon_tissue_scan.R — Baboon Tissue DR Power Scan
#' =======================================================================
#'
#' PURPOSE
#'   1. Scan all 61 baboon tissues from baboon_withTOD.RData for circadian
#'      signal strength (r = A/σ, prop_rhythmic at p<0.05).
#'   2. Select top tissue pairs with high DR potential (large rhythmicity
#'      difference between tissues).
#'   3. Run DR power analysis on selected pairs.
#'
#' DATA
#'   baboon_withTOD.RData (readRDS):
#'     - Object: list with $baboon_data (list of tissues), $tod (ZT times)
#'     - 61 tissues, 6396 genes, raw CPM
#'     - ZT design: 0,2,4,...,22h (12 samples per tissue, n=1/time-point)
#'
#' OUTPUTS
#'   output/baboon_scan_<timestamp>/
#'     tissue_signal.pdf        (ranked bar chart of r per tissue)
#'     dr_power_pairs.pdf       (power curves for selected tissue pairs)
#'     tissue_signal.csv        (summary table)
#'     results.rds
#'
#' USAGE
#'   Rscript examples/exploratory/baboon_tissue_scan.R
#'   POWERSIM_SMOKE=1 Rscript examples/exploratory/baboon_tissue_scan.R
#'
#' @author Thien Pham

# =====================================================================
# SETTINGS
# =====================================================================

SMOKE <- nzchar(Sys.getenv("POWERSIM_SMOKE"))

NGENES       <- if (SMOKE) 500L  else 5000L
NSIMS        <- if (SMOKE) 5L    else 50L
N_GRID       <- if (SMOKE) c(12L, 24L, 36L) else c(12L, 24L, 36L, 48L, 60L, 72L, 96L)
N_TOP_PAIRS  <- if (SMOKE) 2L    else 5L     # number of tissue pairs to power-analyze

RHYTHM_PVAL  <- 0.05   # p-value threshold for rhythmic classification
DATA_PATH    <- "/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Collaborative/GTEXdata/data/CAMO_PRC_hmb.RData"
# Data contains: baboon_withTOD$baboon (list of 61 tissues), baboon_withTOD$tod (list of TOD per tissue)
# baboon: active ZT design (0,2,...,22h), 12 samples/tissue, raw CPM, 4938 human-ortholog genes

cat(sprintf("Mode: %s | NGENES=%d | NSIMS=%d | N_GRID=%s | N_TOP_PAIRS=%d\n",
            if (SMOKE) "SMOKE" else "PRODUCTION",
            NGENES, NSIMS, paste(N_GRID, collapse=","), N_TOP_PAIRS))

# =====================================================================
# SETUP
# =====================================================================

POWERSIM_DIR <- "/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim"
setwd(POWERSIM_DIR)
source_dir <- file.path(getwd(), "code")
old_wd <- setwd(source_dir)
source("setup.R")
setwd(old_wd)

out_dir <- file.path("output", sprintf("baboon_scan_%s", format(Sys.time(), "%Y%m%d_%H%M")))
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("Output: %s/\n", out_dir))

# =====================================================================
# SECTION 1: LOAD DATA
# =====================================================================

cat("\n====================================================================\n")
cat("SECTION 1: LOAD CAMO_PRC_hmb.RData (baboon_withTOD object)\n")
cat("====================================================================\n\n")

load(DATA_PATH)   # loads: gtex, baboon_withTOD, mice
cat(sprintf("Objects loaded: gtex, baboon_withTOD, mice\n"))

# baboon_withTOD$baboon: named list of 61 tissues, each genes x 12 samples (raw CPM)
# baboon_withTOD$tod:    named list of TOD vectors per tissue (ZT 0,2,...,22h)
bab_expr_list <- baboon_withTOD$baboon
bab_tod_list  <- baboon_withTOD$tod    # list, one vector per tissue

tissue_names <- names(bab_expr_list)
n_tissues    <- length(tissue_names)
n_genes_all  <- nrow(bab_expr_list[[1]])

cat(sprintf("Tissues: %d\n", n_tissues))
cat(sprintf("Genes: %d\n", n_genes_all))
cat(sprintf("First tissue ZT: %s\n", paste(bab_tod_list[[1]], collapse=",")))
cat(sprintf("First 10 tissues: %s\n", paste(tissue_names[1:10], collapse=", ")))

# =====================================================================
# SECTION 2: SCAN ALL TISSUES FOR CIRCADIAN SIGNAL
# =====================================================================

cat("\n====================================================================\n")
cat("SECTION 2: SCAN ALL TISSUES — circadian signal strength\n")
cat("====================================================================\n\n")

cat(sprintf("Fitting cosinor to all %d tissues (%d genes each)...\n", n_tissues, n_genes_all))

tissue_signal <- data.frame(
  tissue        = tissue_names,
  prop_rhythmic = NA_real_,
  amp_mean      = NA_real_,
  amp_median    = NA_real_,
  sigma_mean    = NA_real_,
  r_mean        = NA_real_,   # effect size A/sigma
  r_median      = NA_real_,
  stringsAsFactors = FALSE
)

for (i in seq_along(tissue_names)) {
  tis <- tissue_names[i]
  mat <- bab_expr_list[[tis]]

  # log2(CPM+1) transform (raw CPM in data); convert data.frame → matrix
  mat <- data.matrix(log2(mat + 1))
  mat <- mat[rowSums(mat > 0) >= 6, , drop = FALSE]   # keep expressed genes

  # Subsample genes for speed (use full set for scan, NGENES for power)
  g_scan <- min(nrow(mat), 2000L)
  if (g_scan < nrow(mat)) {
    set.seed(42 + i)
    mat <- mat[sample(nrow(mat), g_scan), , drop = FALSE]
  }

  # Per-tissue TOD
  tod_tis <- bab_tod_list[[tis]]

  # Fit cosinor
  fit <- fitCosinorAll(mat, times = tod_tis, period = 24)

  rhythmic    <- !is.na(fit$pvalue) & fit$pvalue < RHYTHM_PVAL
  prop_rhy    <- mean(rhythmic, na.rm=TRUE)
  amp_rhy     <- as.numeric(fit$A[rhythmic])
  sigma_all   <- as.numeric(fit$sigma)

  tissue_signal$prop_rhythmic[i] <- prop_rhy
  tissue_signal$amp_mean[i]      <- if (any(rhythmic)) mean(amp_rhy, na.rm=TRUE) else NA_real_
  tissue_signal$amp_median[i]    <- if (any(rhythmic)) median(amp_rhy, na.rm=TRUE) else NA_real_
  tissue_signal$sigma_mean[i]    <- mean(sigma_all, na.rm=TRUE)
  tissue_signal$r_mean[i]        <- if (any(rhythmic)) mean(amp_rhy / sigma_all[rhythmic], na.rm=TRUE) else NA_real_
  tissue_signal$r_median[i]      <- if (any(rhythmic)) median(amp_rhy / sigma_all[rhythmic], na.rm=TRUE) else NA_real_

  cat(sprintf("  [%2d/%2d] %-8s  prop_rhythmic=%.1f%%  r_mean=%.2f\n",
              i, n_tissues, tis,
              prop_rhy * 100,
              ifelse(is.na(tissue_signal$r_mean[i]), 0, tissue_signal$r_mean[i])))
}

# Sort by r_median descending
tissue_signal <- tissue_signal[order(-tissue_signal$r_median, na.last=TRUE), ]

cat("\n--- Top 15 tissues by effect size (r = A/σ) ---\n")
print(head(tissue_signal[, c("tissue","prop_rhythmic","r_mean","r_median","amp_median","sigma_mean")], 15),
      row.names=FALSE, digits=3)

# Save CSV
write.csv(tissue_signal, file.path(out_dir, "tissue_signal.csv"), row.names=FALSE)
cat(sprintf("\nSaved: %s/tissue_signal.csv\n", out_dir))

# =====================================================================
# SECTION 3: PLOT TISSUE SIGNAL RANKING
# =====================================================================

cat("\n====================================================================\n")
cat("SECTION 3: PLOT tissue signal ranking\n")
cat("====================================================================\n\n")

fig_signal <- file.path(fig_dir, "tissue_signal.pdf")
pdf(fig_signal, width=14, height=8)

# Sort for plotting
ts_plot <- tissue_signal[!is.na(tissue_signal$r_median), ]
ts_plot <- ts_plot[order(ts_plot$r_median), ]   # ascending for horizontal bar

par(mar=c(4,7,3,2))
barplot(
  ts_plot$r_median,
  names.arg = ts_plot$tissue,
  horiz     = TRUE,
  las       = 1,
  col       = ifelse(ts_plot$r_median > median(ts_plot$r_median, na.rm=TRUE),
                     "steelblue", "lightblue"),
  xlab      = "Effect size r = A/σ (median over rhythmic genes, p<0.05)",
  main      = sprintf("Baboon Circadian Signal Strength — %d tissues", n_tissues),
  cex.names = 0.65,
  cex.axis  = 0.8,
  cex.lab   = 0.9
)
abline(v = c(0.5, 1.0, 1.5), lty=2, col="gray50")
legend("bottomright", legend=c("Above median","Below median"),
       fill=c("steelblue","lightblue"), bty="n", cex=0.8)

dev.off()
cat(sprintf("Figure: %s\n", fig_signal))

# =====================================================================
# SECTION 4: SELECT TOP TISSUE PAIRS FOR DR ANALYSIS
# =====================================================================

cat("\n====================================================================\n")
cat("SECTION 4: SELECT TOP TISSUE PAIRS\n")
cat("====================================================================\n\n")

# Strategy: pick pairs with:
#   - Both tissues rhythmic (prop_rhythmic > 0.10)
#   - Large |difference| in prop_rhythmic → high expected prop_DR
#   - Both have sufficient signal (r > 0.5)

good_tis <- tissue_signal$tissue[
  !is.na(tissue_signal$prop_rhythmic) &
  tissue_signal$prop_rhythmic > 0.10 &
  !is.na(tissue_signal$r_median) &
  tissue_signal$r_median > 0.5
]

cat(sprintf("Tissues meeting criteria (prop_rhy>10%%, r>0.5): %d\n", length(good_tis)))
cat(paste(good_tis, collapse=", "), "\n\n")

# Build all pairs and rank by |diff in prop_rhythmic|
if (length(good_tis) >= 2) {
  pairs_df <- do.call(rbind, lapply(seq_along(good_tis), function(i) {
    lapply(seq_along(good_tis), function(j) {
      if (j <= i) return(NULL)
      ti <- good_tis[i]; tj <- good_tis[j]
      ri <- tissue_signal$prop_rhythmic[tissue_signal$tissue == ti]
      rj <- tissue_signal$prop_rhythmic[tissue_signal$tissue == tj]
      data.frame(tis1=ti, tis2=tj, prop_rhy1=ri, prop_rhy2=rj,
                 diff_prop=abs(ri - rj),
                 r1=tissue_signal$r_median[tissue_signal$tissue == ti],
                 r2=tissue_signal$r_median[tissue_signal$tissue == tj],
                 stringsAsFactors=FALSE)
    })
  }))
  pairs_df <- do.call(rbind, pairs_df[!sapply(pairs_df, is.null)])
  pairs_df <- pairs_df[order(-pairs_df$diff_prop), ]

  cat(sprintf("All candidate pairs: %d\n", nrow(pairs_df)))
  cat("Top pairs by |Δprop_rhythmic|:\n")
  print(head(pairs_df, min(10, nrow(pairs_df))), row.names=FALSE, digits=3)

  top_pairs <- head(pairs_df, N_TOP_PAIRS)
} else {
  cat("WARNING: fewer than 2 tissues meet criteria — relaxing r threshold to 0.3\n")
  good_tis <- tissue_signal$tissue[
    !is.na(tissue_signal$prop_rhythmic) &
    tissue_signal$prop_rhythmic > 0.05 &
    !is.na(tissue_signal$r_median) &
    tissue_signal$r_median > 0.3
  ]
  # Build pairs with relaxed criteria
  pairs_df <- do.call(rbind, lapply(seq_along(good_tis), function(i) {
    lapply(seq_along(good_tis), function(j) {
      if (j <= i) return(NULL)
      ti <- good_tis[i]; tj <- good_tis[j]
      ri <- tissue_signal$prop_rhythmic[tissue_signal$tissue == ti]
      rj <- tissue_signal$prop_rhythmic[tissue_signal$tissue == tj]
      data.frame(tis1=ti, tis2=tj, prop_rhy1=ri, prop_rhy2=rj,
                 diff_prop=abs(ri - rj),
                 r1=tissue_signal$r_median[tissue_signal$tissue == ti],
                 r2=tissue_signal$r_median[tissue_signal$tissue == tj],
                 stringsAsFactors=FALSE)
    })
  }))
  pairs_df <- do.call(rbind, pairs_df[!sapply(pairs_df, is.null)])
  pairs_df <- pairs_df[order(-pairs_df$diff_prop), ]
  top_pairs <- head(pairs_df, N_TOP_PAIRS)
}

cat(sprintf("\nSelected %d pairs for DR power analysis:\n", nrow(top_pairs)))
print(top_pairs[, c("tis1","tis2","prop_rhy1","prop_rhy2","diff_prop","r1","r2")],
      row.names=FALSE, digits=3)

# =====================================================================
# SECTION 5: DR POWER ANALYSIS FOR SELECTED PAIRS
# =====================================================================

cat("\n====================================================================\n")
cat("SECTION 5: DR POWER ANALYSIS — selected tissue pairs\n")
cat("====================================================================\n\n")

# Helper: fit cosinor and estimate circadian parameters for one tissue
.fit_tissue <- function(tissue_name, n_genes_use) {
  mat <- bab_expr_list[[tissue_name]]
  mat <- data.matrix(log2(mat + 1))    # data.frame → numeric matrix
  mat <- mat[rowSums(mat > 0) >= 6, , drop=FALSE]

  # Subsample
  if (nrow(mat) > n_genes_use) {
    set.seed(99)
    mat <- mat[sample(nrow(mat), n_genes_use), , drop=FALSE]
  }
  list(data=mat, times=bab_tod_list[[tissue_name]])
}

# Helper: compute prop_DR between two tissues (discordant rhythmicity)
.compute_DR_prop <- function(fit1, fit2, pval_thresh=RHYTHM_PVAL) {
  genes_common <- intersect(rownames(fit1$data), rownames(fit2$data))
  # fitCosinorAll returns rows in same order as input matrix
  # We need to work with named results
  # fit1 and fit2 are lists with $data and $times
  # Re-fit on common genes only
  mat1 <- fit1$data[genes_common, , drop=FALSE]
  mat2 <- fit2$data[genes_common, , drop=FALSE]
  res1 <- fitCosinorAll(mat1, times=fit1$times, period=24)
  res2 <- fitCosinorAll(mat2, times=fit2$times, period=24)
  rhythmic1 <- res1$pvalue < pval_thresh
  rhythmic2 <- res2$pvalue < pval_thresh
  # DR = rhythmic in exactly one group
  dr_mask <- xor(rhythmic1, rhythmic2)
  list(
    prop_DR      = mean(dr_mask),
    prop_rhy1    = mean(rhythmic1),
    prop_rhy2    = mean(rhythmic2),
    n_genes      = length(genes_common)
  )
}

# Analysis options (shared across all pairs)
opts_analysis <- CircadianAnalysisOptions(
  alpha           = 0.05,
  p.adjust.method = "BH",
  fdr_thresholds  = c(0.01, 0.05, 0.10),
  reference_n     = max(N_GRID)
)

# Active design, cts=NULL -> auto evenly-spaced over [0, 24) for each N.
# Baboon pilot used ZT 0,2,...,22 (12 evenly-spaced points); evenly-spaced
# is equivalent and generalizes to any N in the sweep.
opts_design <- CircadianDesignOptions(
  sample_sizes = N_GRID,
  nsims        = NSIMS,
  design       = "active",
  cts          = NULL,
  test_types   = "DR"
)

# Storage for power results
pair_power_list <- list()

for (p_idx in seq_len(nrow(top_pairs))) {
  tis1 <- top_pairs$tis1[p_idx]
  tis2 <- top_pairs$tis2[p_idx]

  cat(sprintf("\n--- Pair %d/%d: %s vs %s ---\n", p_idx, nrow(top_pairs), tis1, tis2))

  # Load and subsample tissue data
  t1 <- .fit_tissue(tis1, NGENES)
  t2 <- .fit_tissue(tis2, NGENES)

  # Find common genes
  genes_common <- intersect(rownames(t1$data), rownames(t2$data))
  n_common     <- length(genes_common)
  cat(sprintf("  Common genes: %d\n", n_common))

  mat1 <- t1$data[genes_common, , drop=FALSE]
  mat2 <- t2$data[genes_common, , drop=FALSE]

  # Fit cosinor and compute prop_DR (use per-tissue TOD)
  res1 <- fitCosinorAll(mat1, times=t1$times, period=24)
  res2 <- fitCosinorAll(mat2, times=t2$times, period=24)
  rhy1 <- !is.na(res1$pvalue) & res1$pvalue < RHYTHM_PVAL
  rhy2 <- !is.na(res2$pvalue) & res2$pvalue < RHYTHM_PVAL
  dr_mask <- xor(rhy1, rhy2)

  prop_DR   <- mean(dr_mask)
  prop_rhy1 <- mean(rhy1)
  prop_rhy2 <- mean(rhy2)
  cat(sprintf("  prop_rhythmic: %s=%.1f%%, %s=%.1f%%  |  prop_DR=%.1f%%\n",
              tis1, prop_rhy1*100, tis2, prop_rhy2*100, prop_DR*100))

  # Estimate circadian params using tissue with larger rhythmic fraction as pilot
  pilot_tis  <- if (prop_rhy1 >= prop_rhy2) tis1 else tis2
  pilot_mat  <- if (pilot_tis == tis1) mat1 else mat2

  cat(sprintf("  Using %s as pilot for parameter estimation\n", pilot_tis))

  pilot_tod  <- bab_tod_list[[pilot_tis]]
  bio_pilot <- estCircadianParam(
    data        = pilot_mat,
    times       = pilot_tod,
    period      = 24,
    prop_DR     = prop_DR,
    prop_DP     = 0,
    prop_DA     = 0,
    min_rhythm_pval = RHYTHM_PVAL
  )

  cat(sprintf("  Bio params: prop_rhythmic=%.1f%%, r_mean=%.2f, prop_DR=%.1f%%\n",
              bio_pilot$prop_rhythmic * 100,
              mean(bio_pilot$amplitude / bio_pilot$sigma, na.rm=TRUE),   # estCircadianParam uses $amplitude
              bio_pilot$prop_DR * 100))

  # Run DR power analysis
  cat(sprintf("  Running power analysis (nsims=%d, N_GRID=%s)...\n",
              NSIMS, paste(N_GRID, collapse=",")))
  power_out <- runPowerAnalysis(
    bio.opts      = bio_pilot,
    design.opts   = opts_design,
    analysis.opts = opts_analysis,
    test_type     = "DR"
  )

  pwr_summary <- summaryRunPower(power_out)
  cat(sprintf("  Power at N=%d: %.1f%%\n",
              N_GRID[length(N_GRID)],
              pwr_summary$Power[nrow(pwr_summary)] * 100))

  pair_power_list[[p_idx]] <- list(
    tis1          = tis1,
    tis2          = tis2,
    prop_rhy1     = prop_rhy1,
    prop_rhy2     = prop_rhy2,
    prop_DR       = prop_DR,
    pilot_tis     = pilot_tis,
    bio_pilot     = bio_pilot,
    power_out     = power_out,
    pwr_summary   = pwr_summary
  )
}

# =====================================================================
# SECTION 6: PLOT DR POWER CURVES
# =====================================================================

cat("\n====================================================================\n")
cat("SECTION 6: PLOT DR power curves\n")
cat("====================================================================\n\n")

fig_power <- file.path(fig_dir, "dr_power_pairs.pdf")
n_pairs   <- length(pair_power_list)

pdf(fig_power, width=7, height=5)

cols <- rainbow(n_pairs, s=0.7, v=0.8)
ltys <- rep(c(1,2,3,4,5), length.out=n_pairs)

# Plot power curves
plot(NA, xlim=range(N_GRID), ylim=c(0,100),
     xlab="N per group", ylab="Power (%)",
     main="Baboon DR Power — Cross-Tissue Pairs")
abline(h=80, lty=2, col="gray40", lwd=1.2)
abline(h=c(20,40,60), lty=3, col="gray80")

legend_labels <- character(n_pairs)
for (p_idx in seq_len(n_pairs)) {
  pr <- pair_power_list[[p_idx]]
  pwr <- pr$pwr_summary$Power * 100

  lines(N_GRID, pwr, col=cols[p_idx], lty=ltys[p_idx], lwd=2)
  points(N_GRID, pwr, col=cols[p_idx], pch=16, cex=0.8)

  legend_labels[p_idx] <- sprintf("%s vs %s (DR=%.0f%%)",
                                  pr$tis1, pr$tis2, pr$prop_DR * 100)
}

legend("bottomright", legend=legend_labels, col=cols, lty=ltys, lwd=2,
       cex=0.7, bty="n")

dev.off()
cat(sprintf("Figure: %s\n", fig_power))

# =====================================================================
# WRAP-UP
# =====================================================================

cat("\n====================================================================\n")
cat("BABOON TISSUE SCAN COMPLETE\n")
cat("====================================================================\n\n")

cat("--- Tissue signal ranking (top 20) ---\n")
print(head(tissue_signal[, c("tissue","prop_rhythmic","r_mean","r_median")], 20),
      row.names=FALSE, digits=3)

cat("\n--- DR power summary ---\n")
for (p_idx in seq_len(n_pairs)) {
  pr  <- pair_power_list[[p_idx]]
  pwr <- pr$pwr_summary$Power * 100
  n80 <- N_GRID[which(pwr >= 80)[1]]
  cat(sprintf("  %s vs %s: prop_DR=%.1f%%, n80=%s\n",
              pr$tis1, pr$tis2, pr$prop_DR * 100,
              ifelse(is.na(n80), ">max(N)", as.character(n80))))
}

# Save results
results <- list(
  tissue_signal   = tissue_signal,
  top_pairs       = top_pairs,
  pair_power_list = pair_power_list
)
saveRDS(results, file.path(out_dir, "results.rds"))
cat(sprintf("\nSaved: %s/results.rds\n", out_dir))
cat(sprintf("Figures: %s/\n", fig_dir))
