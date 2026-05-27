#' =======================================================================
#' fig5_mh_framework.R — Full MH framework on Baboon LUN active pilot
#' =======================================================================
#'
#' Mirror of fig5_fmm_framework.R, but with multi-harmonic (MH) truth.
#' Because the MH (K = 2) data and the K-harmonic regression at K = 2
#' are the same parametric family, the K = 2 detector is correctly
#' specified here (no truncation / no geometric-decay constraint to
#' relax). This sets the upper benchmark on power vs DCP.
#'
#' Pipeline:
#'   cosinor pilot params (estCircadianParam)
#'   -> MH-generated data via runSimsSingleCohort(harmonics = c(alpha2, 0))
#'   -> K = 2 harmonic F-test (method = "FMM", K = 2; equivalent to MH)
#'   -> BH-FDR power calculation.
#'
#' Three panels:
#'   Panel A — marginal power vs N, one line per BH-FDR alpha
#'   Panel B — power stratified by r-tilde at FDR = 0.05, one line per N
#'   Panel C — alpha2 sweep: power vs N at varying alpha2
#'
#' OUTPUT:
#'   output/mh_alternative/results/fig5_mh_framework_<ts>.rds
#'
#' USAGE:
#'   SMOKE_TEST=true Rscript examples/publication/mh_alternative/fig5_mh_framework.R
#'   Rscript examples/publication/mh_alternative/fig5_mh_framework.R

SMOKE_TEST  <- identical(Sys.getenv("SMOKE_TEST"), "true")
NGENES      <- if (SMOKE_TEST) 500L  else 2000L
NSIMS       <- if (SMOKE_TEST) 5L    else 20L
N_CORES     <- as.integer(Sys.getenv("MC_CORES", unset = "20"))
RHYTHM_PVAL <- 0.01
GLOBAL_SEED <- 2025L

N_GRID    <- if (SMOKE_TEST) c(24L, 48L) else
             c(12L, 24L, 36L, 48L, 60L, 72L, 96L, 120L, 144L)
FDR_GRID  <- c(0.01, 0.05, 0.10, 0.20)
# Panel C: per-gene alpha2 ~ Beta(1, eta_a2) sweep. eta -> 0 = mass near 1
# (heavy non-cosinor); large eta = mass near 0 (cosinor-like).
ETA_A2_GRID  <- if (SMOKE_TEST) c(0.5, 2, 20) else c(0.5, 1, 2, 5, 20)
ETA_A2_ANCHOR <- 2.33   # mean(alpha2) approx 0.3, matches scalar anchor

cat(sprintf("Mode         : %s\n", if (SMOKE_TEST) "SMOKE" else "PRODUCTION"))
cat(sprintf("NGENES       : %d\n", NGENES))
cat(sprintf("NSIMS        : %d\n", NSIMS))
cat(sprintf("N_GRID       : %s\n", paste(N_GRID,      collapse = ", ")))
cat(sprintf("ETA_A2_GRID  : %s\n", paste(ETA_A2_GRID, collapse = ", ")))
cat(sprintf("ETA_A2_ANCHOR: %g  (mean alpha2 = %.3f)\n",
            ETA_A2_ANCHOR, 1 / (1 + ETA_A2_ANCHOR)))
cat(sprintf("FDR_GRID     : %s\n", paste(FDR_GRID,    collapse = ", ")))

setwd("/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")
old_wd <- setwd("code"); source("setup.R"); setwd(old_wd)

out_dir_res <- "output/mh_alternative/results"
dir.create(out_dir_res, recursive = TRUE, showWarnings = FALSE)

# ====================================================================
# 1. Baboon LUN pilot (cosinor estimation)
# ====================================================================
cat("\n=== Step 1: load + cosinor-fit Baboon LUN ===\n")
load("data/CAMO_PRC_hmb.RData")
prep_lun <- prepCircadianData(baboon_withTOD$baboon[["LUN"]],
                               times = baboon_withTOD$tod[["LUN"]] %% 24,
                               input_type = "cpm")
mat_lun  <- prep_lun$data[rowSums(prep_lun$data > 0) >= 6, , drop = FALSE]
tod_lun  <- prep_lun$times
rm(baboon_withTOD, gtex, mice, prep_lun)
set.seed(GLOBAL_SEED)
g_idx   <- sample(nrow(mat_lun), min(NGENES, nrow(mat_lun)))
mat_sub <- mat_lun[g_idx, ]

bio_lun <- estCircadianParam(mat_sub, tod_lun,
                              min_rhythm_pval = RHYTHM_PVAL,
                              verbose = TRUE)
bio_lun$ngenes <- NGENES
bio_lun$omega_rhythmic <- NULL; bio_lun$alpha_rhythmic <- NULL
bio_lun$omega_dist <- NULL;     bio_lun$alpha_dist <- NULL

# Active design B = 12
B_VAL <- 12L
cts_active <- seq(0, 24 * (1 - 1/B_VAL), length.out = B_VAL)

# r-strata for Panel B
r_strata    <- makeAdaptiveRStrata(bio_lun, bin_width = 0.25)
strata_lbls <- paste0("(",
                      round(head(r_strata, -1), 2), ", ",
                      round(tail(r_strata, -1), 2), "]")
cat(sprintf("r_strata: %s\n",
            paste(round(r_strata[is.finite(r_strata)], 2), collapse = ", ")))

# ====================================================================
# 2. Panels A and B: shared sim at alpha2 = ALPHA2_ANCHOR
#    Save pvalues + r_values_list for offline re-thresholding / stratification.
# ====================================================================
cat(sprintf("\n=== Step 2: N sweep at alpha2_g ~ Beta(1, %g) (K = 2 detector) ===\n",
            ETA_A2_ANCHOR))
analysis_AB <- CircadianAnalysisOptions(alpha           = 0.05,
                                         p.adjust.method = "BH",
                                         fdr_thresholds  = FDR_GRID,
                                         r_strata        = r_strata,
                                         strata_labels   = strata_lbls)

G       <- bio_lun$ngenes
n_N     <- length(N_GRID)
n_sims  <- NSIMS
n_strat <- length(strata_lbls)
n_fdr   <- length(FDR_GRID)

power_marg  <- array(NA_real_, dim = c(n_N, n_fdr, n_sims))
power_strat <- array(NA_real_, dim = c(n_N, n_strat, n_fdr, n_sims))

# Anchor bio.opts for Panels A and B: per-gene alpha2 ~ Beta(1, ETA_A2_ANCHOR)
bio_anchor <- bio_lun
bio_anchor$alpha2_dist <- list(family = "beta", a = 1, b = ETA_A2_ANCHOR)
bio_anchor$alpha3_dist <- NULL

for (j in seq_along(N_GRID)) {
  N <- N_GRID[j]
  if (N %% B_VAL != 0L) next
  m <- N %/% B_VAL
  cts_design <- rep(cts_active, each = m)
  design <- CircadianDesignOptions(sample_sizes = N, nsims = n_sims,
                                    design = "active",
                                    cts = cts_design, B_values = B_VAL)
  set.seed(GLOBAL_SEED)
  res <- runSimsSingleCohort(bio_anchor, design, analysis_AB,
                              method   = "FMM", K = 2L,
                              mc.cores = N_CORES, verbose = FALSE)
  pvals      <- res$pvalues[1, , ]
  rvals_list <- res$r_values_list[[1]]

  for (s in seq_len(n_sims)) {
    p_s <- pvals[, s]
    r_s <- rvals_list[[s]]
    is_rh <- r_s > 0
    r_for_strat <- ifelse(is_rh, r_s, 0)
    xgr <- cut(r_for_strat, breaks = r_strata, include.lowest = TRUE, labels = FALSE)
    xgr[!is_rh] <- NA

    for (a_i in seq_len(n_fdr)) {
      fdr_g <- p.adjust(p_s, method = "BH")
      disc  <- fdr_g <= FDR_GRID[a_i]
      n_rh  <- sum(is_rh)
      TD    <- sum(disc & is_rh, na.rm = TRUE)
      power_marg[j, a_i, s] <- if (n_rh > 0) TD / n_rh else NA_real_
      for (k in seq_len(n_strat)) {
        in_k <- !is.na(xgr) & xgr == k
        n_tgt <- sum(is_rh & in_k)
        TD_k  <- sum(disc & is_rh & in_k, na.rm = TRUE)
        power_strat[j, k, a_i, s] <- if (n_tgt > 0) TD_k / n_tgt else NA_real_
      }
    }
  }
  cat(sprintf("  N=%3d: power @ FDR 0.05 = %.3f\n", N,
              mean(power_marg[j, which(FDR_GRID == 0.05), ], na.rm = TRUE)))
}

# ====================================================================
# 3. Panel C: eta_alpha2 sweep at FDR = 0.05 (K = 2 detector)
# ====================================================================
cat("\n=== Step 3: eta_alpha2 sweep (K = 2 detector, FDR = 0.05) ===\n")
analysis_C <- CircadianAnalysisOptions(alpha           = 0.05,
                                        p.adjust.method = "BH",
                                        fdr_thresholds  = 0.05)
power_a2 <- matrix(NA_real_, nrow = length(N_GRID), ncol = length(ETA_A2_GRID),
                   dimnames = list(paste0("N=", N_GRID),
                                   paste0("eta_a2=", ETA_A2_GRID)))
for (k in seq_along(ETA_A2_GRID)) {
  bio_k <- bio_lun
  bio_k$alpha2_dist <- list(family = "beta", a = 1, b = ETA_A2_GRID[k])
  bio_k$alpha3_dist <- NULL
  pwr <- vapply(N_GRID, function(N) {
    if (N %% B_VAL != 0L) return(NA_real_)
    m <- N %/% B_VAL
    cts_design <- rep(cts_active, each = m)
    design <- CircadianDesignOptions(sample_sizes = N, nsims = n_sims,
                                      design = "active",
                                      cts = cts_design, B_values = B_VAL)
    set.seed(GLOBAL_SEED)
    res <- runSimsSingleCohort(bio_k, design, analysis_C,
                                method   = "FMM", K = 2L,
                                mc.cores = N_CORES, verbose = FALSE)
    mean(res$marginal_power, na.rm = TRUE)
  }, numeric(1))
  power_a2[, k] <- pwr
  cat(sprintf("  eta_a2=%.2f: %s\n", ETA_A2_GRID[k],
              paste(sprintf("%.2f", pwr), collapse = " ")))
}

# ====================================================================
# 4. Save
# ====================================================================
out <- list(
  framework       = "MH alternative -- full MH framework with per-gene Beta(1, eta) alpha2",
  detector        = "K = 2 harmonic F-test (LimoRhyde / limma)",
  pilot           = "baboon LUN (active, GSE98965)",
  GLOBAL_SEED     = GLOBAL_SEED,
  NGENES          = NGENES, NSIMS = NSIMS,
  N_GRID          = N_GRID,
  ETA_A2_GRID     = ETA_A2_GRID,
  FDR_GRID        = FDR_GRID,
  eta_a2_anchor   = ETA_A2_ANCHOR,
  r_strata        = r_strata,
  strata_labels   = strata_lbls,
  power_marg      = power_marg,
  power_strat     = power_strat,
  power_eta       = power_a2     # naming aligned with FMM-version replot
)
ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
rds <- file.path(out_dir_res, sprintf("fig5_mh_framework_%s.rds", ts))
saveRDS(out, rds)
cat(sprintf("\nSaved: %s\n", rds))
