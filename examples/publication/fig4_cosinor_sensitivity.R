#' =======================================================================
#' fig4_cosinor_sensitivity.R — Sensitivity of cosinor-based SCP (DCP)
#'                              to FMM-shaped data
#' =======================================================================
#'
#' Companion to fig4_sensitivity.R (which evaluates the full FMM
#' framework's K=2 detector). This script evaluates the COSINOR-BASED
#' analysis pipeline (cosinor pilot estimation + DCP detection) under
#' FMM-simulated data, quantifying the empirical cost of analysing
#' non-cosinor waveforms with a cosinor detector.
#'
#' Pipeline (cosinor framework, §2.4 of paper):
#'   (i)   pilot params via estCircadianParam (cosinor fits, not FMM)
#'   (ii)  data generated under FMM truth via omega_dist/alpha_dist
#'   (iii) detection via DCP (method="DCP" in runSimsSingleCohort)
#'   (iv)  power evaluated at BH-FDR alpha = 0.05
#'
#' Panel A (omega/eta sweep): omega_g ~ Beta(1, eta), eta in {0, 0.5, 1, 2, 5, 20};
#'   alpha distribution fixed at empirical sigma_alpha.
#' Panel B (alpha sweep): omega distribution fixed at empirical Beta(1, eta_hat);
#'   alpha_g ~ vonMises(0, kappa), kappa = (24/(2*pi*sigma_alpha))^2,
#'   sigma_alpha in {0, 0.5, 1, 2, 4} h.
#'   Includes eta = 0 baseline (pure cosinor truth) per Sections 2.4 plan.
#'
#' OUTPUT:
#'   output/sensitivity/results/fig4_cosinor_sensitivity_<ts>.rds
#'
#' USAGE:
#'   SMOKE_TEST=true  Rscript examples/publication/fig4_cosinor_sensitivity.R
#'   Rscript examples/publication/fig4_cosinor_sensitivity.R
#'
#' REPRODUCIBILITY: see fig4_sensitivity.R header — identical conventions
#' (GLOBAL_SEED, NGENES, NSIMS, MC_CORES, anchor pilot baboon LUN at B=12).

SMOKE_TEST  <- identical(Sys.getenv("SMOKE_TEST"), "true")
NGENES      <- if (SMOKE_TEST) 500L  else 2000L
NSIMS       <- if (SMOKE_TEST) 5L    else 20L
N_CORES     <- as.integer(Sys.getenv("MC_CORES", unset = "20"))
TOP_K_FMM   <- if (SMOKE_TEST) 50L   else 300L
RHYTHM_PVAL <- 0.01
GLOBAL_SEED <- 2025L

N_GRID    <- if (SMOKE_TEST) c(24L, 36L) else c(12L, 24L, 36L, 48L, 60L, 72L, 96L, 120L, 144L)
# Panel A: eta sweep (eta = 0 represents pure cosinor truth, full power baseline)
BETA_GRID <- if (SMOKE_TEST) c(0, 1, 5)  else c(0, 0.5, 1, 2, 5, 20)
# Panel B: sigma_alpha sweep
SDHR_GRID <- if (SMOKE_TEST) c(0, 1, 4)  else c(0, 0.5, 1, 2, 4)

cat(sprintf("Mode      : %s\n", if (SMOKE_TEST) "SMOKE" else "PRODUCTION"))
cat(sprintf("NGENES    : %d\n", NGENES))
cat(sprintf("NSIMS     : %d\n", NSIMS))
cat(sprintf("N_GRID    : %s\n", paste(N_GRID,    collapse = ", ")))
cat(sprintf("BETA_GRID : %s\n", paste(BETA_GRID, collapse = ", ")))
cat(sprintf("SDHR_GRID : %s\n", paste(SDHR_GRID, collapse = ", ")))
cat(sprintf("MC_CORES  : %d\n", N_CORES))

old_wd <- setwd("code"); source("setup.R"); setwd(old_wd)

out_dir_res <- "output/sensitivity/results"
dir.create(out_dir_res, recursive = TRUE, showWarnings = FALSE)

analysis <- CircadianAnalysisOptions(alpha           = 0.05,
                                      p.adjust.method = "BH",
                                      fdr_thresholds  = 0.05)

# ====================================================================
# 1. Anchor pilot: Baboon LUN — cosinor estimation
# ====================================================================
# We also load the cached FMM fit (bio_lun_fmm) to recover sigma_alpha_hat
# and eta_hat anchors for the sweep — these characterise the *data-generating*
# regime even though the *analysis* uses cosinor estimation.
cat("\n=== Step 1: load Baboon LUN ===\n")
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

# Cosinor-based pilot estimate (the §2.4 analysis framework)
cat("Step 1a: cosinor-based pilot estimation (estCircadianParam)...\n")
bio_lun_cos <- estCircadianParam(mat_sub, tod_lun,
                                  min_rhythm_pval = RHYTHM_PVAL,
                                  verbose = TRUE)
bio_lun_cos$ngenes <- NGENES

# Sweep anchors from the FMM-fit diagnostics (data-generating regime)
fmm_anchor_rds <- file.path(out_dir_res,
                            sprintf("bio_lun_fmm_NG%d_K%d.rds",
                                    NGENES, TOP_K_FMM))
if (file.exists(fmm_anchor_rds)) {
  bio_lun_fmm <- readRDS(fmm_anchor_rds)
  beta_hat        <- bio_lun_fmm$diagnostics$beta_hat
  sigma_alpha_hat <- bio_lun_fmm$diagnostics$sigma_alpha_hat
} else {
  cat("FMM diagnostic anchors not cached; running estCircadianParamFMM ...\n")
  bio_lun_fmm <- estCircadianParamFMM(mat_sub, tod_lun,
                                       min_rhythm_pval = RHYTHM_PVAL,
                                       top_k    = TOP_K_FMM,
                                       mc.cores = N_CORES,
                                       verbose  = TRUE)
  bio_lun_fmm$ngenes <- NGENES
  saveRDS(bio_lun_fmm, fmm_anchor_rds)
  beta_hat        <- bio_lun_fmm$diagnostics$beta_hat
  sigma_alpha_hat <- bio_lun_fmm$diagnostics$sigma_alpha_hat
}
cat(sprintf("Sweep anchors: eta_hat=%.3f  sigma_alpha_hat=%.3fh\n",
            beta_hat, sigma_alpha_hat))

# Active design: equispaced B = 12, m = N/12 reps each
B_VAL <- 12L
cts_active <- seq(0, 24 * (1 - 1/B_VAL), length.out = B_VAL)

run_sweep_curve <- function(bio_k, label) {
  pwr <- vapply(N_GRID, function(N) {
    if (N %% B_VAL != 0L) return(NA_real_)
    m <- N %/% B_VAL
    cts_design <- rep(cts_active, each = m)
    design <- CircadianDesignOptions(sample_sizes = N, nsims = NSIMS,
                                      design = "active",
                                      cts = cts_design, B_values = B_VAL)
    set.seed(GLOBAL_SEED)
    res <- runSimsSingleCohort(bio_k, design, analysis,
                                method   = "DCP",   # cosinor F-test, K = 1
                                mc.cores = N_CORES, verbose = FALSE)
    mean(res$marginal_power, na.rm = TRUE)
  }, numeric(1))
  cat(sprintf("  %-22s : %s\n", label,
              paste(sprintf("%.2f", pwr), collapse = " ")))
  pwr
}

# Helper: set FMM-truth sweep parameters on the cosinor pilot bio.opts
set_fmm_truth <- function(bio_base, eta, sd_hours) {
  bio_k <- bio_base
  bio_k$omega_rhythmic <- NULL
  bio_k$alpha_rhythmic <- NULL
  if (eta == 0) {
    # Cosinor truth: omega == 1 for all genes
    bio_k$omega_dist <- list(family = "fixed", value = 1.0)
  } else {
    bio_k$omega_dist <- list(family = "beta", a = 1, b = eta)
  }
  bio_k$alpha_dist <- list(family = "normal", mean = 0, sd_hours = sd_hours)
  bio_k
}

# ====================================================================
# 2. Panel A — eta sweep (sigma_alpha fixed at empirical)
# ====================================================================
cat("\n=== Step 2: eta sweep (DCP under FMM truth) ===\n")
power_A <- matrix(NA_real_, nrow = length(N_GRID), ncol = length(BETA_GRID),
                  dimnames = list(paste0("N=", N_GRID),
                                  paste0("eta=", BETA_GRID)))
for (k in seq_along(BETA_GRID)) {
  bio_k <- set_fmm_truth(bio_lun_cos, eta = BETA_GRID[k],
                          sd_hours = sigma_alpha_hat)
  power_A[, k] <- run_sweep_curve(bio_k, sprintf("eta=%.2f", BETA_GRID[k]))
}

# ====================================================================
# 3. Panel B — sigma_alpha sweep (omega-dist fixed at empirical Beta(1, eta_hat))
# ====================================================================
cat("\n=== Step 3: sigma_alpha sweep (DCP under FMM truth) ===\n")
power_B <- matrix(NA_real_, nrow = length(N_GRID), ncol = length(SDHR_GRID),
                  dimnames = list(paste0("N=", N_GRID),
                                  paste0("sdh=", SDHR_GRID)))
for (k in seq_along(SDHR_GRID)) {
  bio_k <- set_fmm_truth(bio_lun_cos, eta = beta_hat,
                          sd_hours = SDHR_GRID[k])
  power_B[, k] <- run_sweep_curve(bio_k, sprintf("sd_alpha=%.1fh", SDHR_GRID[k]))
}

# ====================================================================
# 4. Save
# ====================================================================
out <- list(
  bio_lun_cos     = bio_lun_cos,
  beta_hat        = beta_hat,
  sigma_alpha_hat = sigma_alpha_hat,
  N_grid          = N_GRID,
  beta_grid       = BETA_GRID,
  sdhr_grid       = SDHR_GRID,
  power_A         = power_A,
  power_B         = power_B,
  detector        = "DCP (cosinor F-test, K = 1)",
  framework       = "cosinor-based SCP (Section 2.4 of paper)"
)
ts       <- format(Sys.time(), "%Y%m%d_%H%M%S")
rds_path <- file.path(out_dir_res,
                      sprintf("fig4_cosinor_sensitivity_%s.rds", ts))
saveRDS(out, rds_path)
cat(sprintf("\nSaved: %s\n", rds_path))
