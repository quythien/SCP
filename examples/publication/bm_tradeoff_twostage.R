#' =======================================================================
#' bm_tradeoff_twostage.R — B vs m Tradeoff (Two-Stage DR Framework)
#' =======================================================================
#'
#' Regenerates paper/PowerSim/figures/bm_tradeoff.pdf using the two-stage
#' (plug-in) DR framework. Three active-design dataset pairs, unified
#' B = {3,4,6,8,12} and N = {12,24,48,72,96,120,144,192,240,288}.
#'
#' Panels ordered left to right by decreasing r (SNR):
#'   A. Mouse LIV vs CER  (GSE54651, r~2.9, prop_DR~27%)
#'   B. Baboon LUN vs CER (CAMO,     r~1.7, prop_DR~41%)
#'   C. Mouse D1 vs D2    (D1D2,     r~0.65, prop_DR~20%)
#'
#' USAGE:
#'   Rscript examples/publication/bm_tradeoff_twostage.R
#'   SMOKE_TEST=true Rscript examples/publication/bm_tradeoff_twostage.R
#'
#' @author Thien Pham

# =====================================================================
# 0. Setup
# =====================================================================
SMOKE_TEST  <- identical(Sys.getenv("SMOKE_TEST"), "true")
GLOBAL_SEED <- 2025L
set.seed(GLOBAL_SEED)

old_wd <- setwd("code"); source("setup.R"); setwd(old_wd)

B_VALS      <- c(3L, 4L, 6L, 8L, 12L)
N_GRID      <- if (SMOKE_TEST) c(12L, 24L, 48L) else
                 c(12L, 24L, 48L, 72L, 96L, 120L, 144L, 192L, 240L, 288L)
NSIMS       <- if (SMOKE_TEST) 5L   else 200L
NGENES_CORE <- if (SMOKE_TEST) 500L else 5000L
RHYTHM_PVAL <- 0.05
N_CORES     <- as.integer(Sys.getenv("MC_CORES", unset = "4"))

cat(sprintf("Mode     : %s\n", if (SMOKE_TEST) "SMOKE" else "PRODUCTION"))
cat(sprintf("B values : %s\n", paste(B_VALS, collapse = ", ")))
cat(sprintf("N grid   : %s\n", paste(N_GRID,  collapse = ", ")))
cat(sprintf("nsims    : %d\n", NSIMS))
cat(sprintf("ngenes   : %d\n", NGENES_CORE))
cat(sprintf("mc.cores : %d\n", N_CORES))

out_dir   <- "output/bm_tradeoff_twostage"
cache_dir <- file.path(out_dir, "cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

opts_analysis <- CircadianAnalysisOptions(
  alpha           = 0.05,
  p.adjust.method = "BH",
  parallel.ncores = 1L
)

# =====================================================================
# Helper: evenly-spaced active design vector for B time points
# =====================================================================
make_design_vec <- function(B) seq(0, 24 * (1 - 1/B), length.out = B)

# =====================================================================
# 1. Load datasets and estimate two-group pilots
# =====================================================================
cat("\n--- Loading datasets ---\n")

# ----- A: Mouse LIV vs CER (GSE54651) --------------------------------
dat_mouse <- readRDS("data/mice_GSE54651_CPM.RData")
prep_liv  <- prepCircadianData(dat_mouse$count_clean[["LIV"]],
                               times = dat_mouse$tod[["LIV"]], input_type = "log2")
prep_cer_liv <- prepCircadianData(dat_mouse$count_clean[["CER"]],
                                  times = dat_mouse$tod[["CER"]], input_type = "log2")
rm(dat_mouse)
keep_liv <- rowSums(prep_liv$data > 0) >= 4
keep_cer_liv <- rowSums(prep_cer_liv$data > 0) >= 4
common_liv <- intersect(rownames(prep_liv$data)[keep_liv],
                         rownames(prep_cer_liv$data)[keep_cer_liv])
set.seed(7L)
g_liv    <- sample(common_liv, min(NGENES_CORE, length(common_liv)))
mat_liv  <- prep_liv$data[g_liv, ]
mat_cer_liv <- prep_cer_liv$data[g_liv, ]
tod_liv  <- prep_liv$times
tod_cer_liv <- prep_cer_liv$times
rm(prep_liv, prep_cer_liv)

# ----- B: Baboon LUN vs CER ------------------------------------------
load("data/CAMO_PRC_hmb.RData")
prep_lun <- prepCircadianData(baboon_withTOD$baboon[["LUN"]],
                               times = baboon_withTOD$tod[["LUN"]], input_type = "cpm")
prep_cer_bab <- prepCircadianData(baboon_withTOD$baboon[["CER"]],
                                   times = baboon_withTOD$tod[["CER"]], input_type = "cpm")
rm(baboon_withTOD, gtex, mice)
keep_lun <- rowSums(prep_lun$data > 0) >= 6
keep_cer_bab <- rowSums(prep_cer_bab$data > 0) >= 6
common_bab <- intersect(rownames(prep_lun$data)[keep_lun],
                         rownames(prep_cer_bab$data)[keep_cer_bab])
set.seed(7L)
g_bab    <- sample(common_bab, min(NGENES_CORE, length(common_bab)))
mat_lun  <- prep_lun$data[g_bab, ]
mat_cer_bab <- prep_cer_bab$data[g_bab, ]
tod_lun  <- prep_lun$times
tod_cer_bab <- prep_cer_bab$times
rm(prep_lun, prep_cer_bab)

# ----- C: Mouse D1 vs D2 ---------------------------------------------
pheno  <- read.csv("data/mouse_clinicalinfo_03082021_rmOutliers.csv", row.names = 1)
prep_d1d2 <- prepCircadianData("data/mouse_D1D2_logCPMfiltered_counts.csv",
                                times = "time", input_type = "counts",
                                pheno = pheno, sample_col = "sample")
log_mat <- prep_d1d2$data
d1_samp <- pheno$sample[pheno$cell == "D1"]
d2_samp <- pheno$sample[pheno$cell == "D2"]
mat_d1  <- log_mat[, colnames(log_mat) %in% d1_samp, drop = FALSE]
mat_d2  <- log_mat[, colnames(log_mat) %in% d2_samp, drop = FALSE]
tod_d1  <- pheno$time[match(colnames(mat_d1), pheno$sample)]
tod_d2  <- pheno$time[match(colnames(mat_d2), pheno$sample)]
keep_d1 <- rowSums(mat_d1 > 1) >= floor(ncol(mat_d1) / 2)
mat_d1  <- mat_d1[keep_d1, ]; mat_d2 <- mat_d2[keep_d1, ]
set.seed(42L)
g_d1    <- sample(nrow(mat_d1), min(NGENES_CORE, nrow(mat_d1)))
mat_d1  <- mat_d1[g_d1, ]; mat_d2 <- mat_d2[g_d1, ]
rm(pheno, prep_d1d2, log_mat)

cat(sprintf("  LIV: %d genes x %d samples | CER: %d samples\n",
            nrow(mat_liv), ncol(mat_liv), ncol(mat_cer_liv)))
cat(sprintf("  LUN: %d genes x %d samples | CER: %d samples\n",
            nrow(mat_lun), ncol(mat_lun), ncol(mat_cer_bab)))
cat(sprintf("  D1:  %d genes x %d samples | D2:  %d samples\n",
            nrow(mat_d1), ncol(mat_d1), ncol(mat_d2)))

# =====================================================================
# 2. Estimate two-group pilots (with caching)
# =====================================================================
cat("\n--- Estimating two-group pilots ---\n")

estimate_or_cache <- function(tag, data1, times1, data2, times2) {
  rds <- file.path(cache_dir, sprintf("%s_pilot.rds", tag))
  if (file.exists(rds)) {
    cat(sprintf("  [cache] pilot %s\n", tag)); return(readRDS(rds))
  }
  cat(sprintf("  [fit]   pilot %s\n", tag))
  bio <- estCircadianParamTwoGroup(
    data_1 = data1, data_2 = data2,
    times_1 = times1, times_2 = times2,
    period = 24, min_rhythm_pval = RHYTHM_PVAL, verbose = FALSE
  )
  bio <- updateBioOptions(bio, ngenes = NGENES_CORE)
  bio <- updateBioOptions(bio,
    prop_DR    = bio$prop_DR,
    prop_DP    = 0.00,
    phase_diff = c(0, 0),
    amp_diff   = c(1, 1)
  )
  saveRDS(bio, rds); bio
}

pilot_liv <- estimate_or_cache("liv", mat_liv, tod_liv, mat_cer_liv, tod_cer_liv)
pilot_bab <- estimate_or_cache("bab", mat_lun, tod_lun, mat_cer_bab, tod_cer_bab)
pilot_d1  <- estimate_or_cache("d1",  mat_d1,  tod_d1,  mat_d2,      tod_d2)
rm(mat_liv, mat_cer_liv, mat_lun, mat_cer_bab, mat_d1, mat_d2)

# Compute r_med for each pilot (for labels)
# Two-group pilots from estCircadianParamTwoGroup do not populate sigma_rhythmic;
# fall back to exp(median lOD) as the noise estimate.
r_med <- function(bio) {
  if (!is.null(bio$sigma_rhythmic) && length(bio$sigma_rhythmic) > 0) {
    median(bio$amplitude / bio$sigma_rhythmic, na.rm = TRUE)
  } else {
    sigma_med <- exp(median(bio$lOD, na.rm = TRUE))
    median(bio$amplitude / sigma_med, na.rm = TRUE)
  }
}

datasets <- list(
  liv = list(bio = pilot_liv, label = "Mouse LIV vs CER\n(GSE54651, active)",
             color = "steelblue",   tod = tod_liv),
  bab = list(bio = pilot_bab, label = "Baboon LUN vs CER\n(CAMO, active)",
             color = "darkorange",  tod = tod_lun),
  d1  = list(bio = pilot_d1,  label = "Mouse D1 vs D2\n(D1D2, active)",
             color = "forestgreen", tod = tod_d1)
)

for (nm in names(datasets)) {
  bio <- datasets[[nm]]$bio
  cat(sprintf("  [%s] prop_DR=%.1f%%  r_med=%.2f  n_pilot=%d\n",
              nm, 100 * bio$prop_DR, r_med(bio), length(datasets[[nm]]$tod)))
}

# =====================================================================
# 3. Run two-stage DR power for all (dataset x B) combinations
# =====================================================================
cat("\n--- Running two-stage DR power grid (parallel over dataset x B) ---\n")

jobs <- expand.grid(nm = names(datasets), B = B_VALS,
                    stringsAsFactors = FALSE)
cat(sprintf("Total jobs: %d  (mc.cores=%d)\n", nrow(jobs), N_CORES))

all_results <- parallel::mclapply(seq_len(nrow(jobs)), function(i) {
  nm <- jobs$nm[i]
  B  <- jobs$B[i]
  rds <- file.path(cache_dir, sprintf("%s_B%02d.rds", nm, B))
  if (file.exists(rds)) {
    cat(sprintf("  [cache] %s B=%d\n", nm, B))
    return(readRDS(rds))
  }
  cat(sprintf("  [run]   %s B=%d\n", nm, B))
  bio <- datasets[[nm]]$bio
  bio$cts2 <- NULL  # both groups use the same active B-point design
  opts_design <- CircadianDesignOptions(
    sample_sizes = N_GRID,
    nsims        = NSIMS,
    design       = "active",
    cts          = make_design_vec(B),
    test_types   = "DR"
  )
  set.seed(GLOBAL_SEED + B)
  res <- runPowerAnalysis(bio, opts_design, opts_analysis, test_type = "DR")
  saveRDS(res, rds)
  res
}, mc.cores = N_CORES)

# Reshape into named list: results[[nm]][[paste0("B", B)]]
results <- lapply(names(datasets), function(nm) {
  idx <- which(jobs$nm == nm)
  r   <- all_results[idx]
  names(r) <- paste0("B", B_VALS)
  r
})
names(results) <- names(datasets)

# =====================================================================
# 4. Generate figure
# =====================================================================
fig_path <- "paper/PowerSim/figures/bm_tradeoff.pdf"
dir.create(dirname(fig_path), recursive = TRUE, showWarnings = FALSE)

pdf(fig_path, width = 17, height = 5)
par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3.5, 1.5))

for (nm in names(datasets)) {
  d    <- datasets[[nm]]
  bio  <- d$bio
  cols <- colorRampPalette(c("gray75", d$color))(length(B_VALS))
  pchs <- c(16, 17, 15, 18, 8)[seq_along(B_VALS)]

  plot(NA, xlim = range(N_GRID), ylim = c(0, 100), las = 1, xaxt = "n",
       xlab = "N per group", ylab = "Power (FDR 5%)",
       main = sprintf("%s\nprop_DR=%.0f%%  r_med=%.2f",
                      d$label, 100 * bio$prop_DR, r_med(bio)))
  axis(1, at = N_GRID, las = 2, cex.axis = 0.8)
  abline(h = c(20, 40, 60, 80), lty = 3, col = "gray85")

  for (bi in seq_along(B_VALS)) {
    B   <- B_VALS[bi]
    res <- results[[nm]][[paste0("B", B)]]
    pwr <- rowMeans(res$marginal_power, na.rm = TRUE) * 100
    lines(N_GRID,  pwr, col = cols[bi], lwd = 2, lty = bi)
    points(N_GRID, pwr, col = cols[bi], pch = pchs[bi], cex = 0.9)
  }

  m_at_min <- pmax(1L, floor(min(N_GRID) / B_VALS))
  legend("bottomright",
         legend = sprintf("B=%2d (m~%d at N=%d)", B_VALS, m_at_min, min(N_GRID)),
         col = cols, lwd = 2, lty = seq_along(B_VALS), pch = pchs,
         bty = "n", cex = 0.82)
}

dev.off()
cat(sprintf("\nFigure saved: %s\n", fig_path))

# =====================================================================
# 5. Print n80 summary
# =====================================================================
cat("\n--- n80 summary (DR power >= 80%) ---\n")
for (nm in names(datasets)) {
  cat(sprintf("\n%s:\n", datasets[[nm]]$label))
  for (bi in seq_along(B_VALS)) {
    B   <- B_VALS[bi]
    pwr <- rowMeans(results[[nm]][[paste0("B", B)]]$marginal_power, na.rm = TRUE)
    n80 <- N_GRID[which(pwr >= 0.80)[1]]
    cat(sprintf("  B=%2d: n80=%s\n", B,
                ifelse(is.na(n80), ">max", as.character(n80))))
  }
}
cat("\n=== Done ===\n")
