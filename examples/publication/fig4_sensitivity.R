#' =======================================================================
#' fig4_sensitivity.R — Sensitivity to cosinor violation
#' =======================================================================
#'
#' Two orthogonal sensitivity sweeps anchored on Baboon LUN (CAMO):
#'
#'   Panel A — omega sweep:   omega_g ~ Beta(1, beta), sweep beta in
#'             {0.01, 0.1, 0.5, 1, 2, 5, 20, 50}.
#'             beta -> 0+ : pure cosinor (omega ~ 1)
#'             beta = 1   : Uniform(0, 1)
#'             beta -> Inf: arrhythmic limit (omega ~ 0)
#'             Empirical beta_hat from Baboon LUN annotated as a vertical
#'             reference line.
#'
#'   Panel B — alpha sweep:   alpha_g ~ vonMises(alpha_hat_g, kappa = 1/sd_rad^2),
#'             sweep sd_hours in {0, 0.5, 1, 2, 4}.
#'             sd_hours = 0: uses empirical alpha_hat exactly
#'             sd_hours large: phase becomes circular-uniform around alpha_hat
#'             Empirical sigma_alpha_hat (std dev of alpha_hat in hours)
#'             annotated as a vertical reference line.
#'
#' Layout: 2x2 grid
#'   (top-left)  empirical omega_hat histogram + Beta(1, beta_hat) overlay
#'   (top-right) empirical alpha_hat histogram + von Mises overlay
#'   (bottom-left)  power vs beta (omega sweep)
#'   (bottom-right) power vs sd_hours (alpha sweep)
#'
#' USAGE:
#'   Rscript examples/publication/fig4_sensitivity.R
#'   SMOKE_TEST=true Rscript examples/publication/fig4_sensitivity.R

SMOKE_TEST  <- identical(Sys.getenv("SMOKE_TEST"), "true")
NGENES      <- if (SMOKE_TEST) 500L  else 5000L
NSIMS       <- if (SMOKE_TEST) 10L   else 50L
# N_FIXED is the sample size where power is evaluated. Choose N near the
# elbow of the cosinor power curve so the omega sweep traverses 0% -> 100%.
# For Baboon LUN (active design, B=12), N=20 produces ~0.84 cosinor power.
N_FIXED     <- as.integer(Sys.getenv("N_FIXED", unset = "20"))
N_CORES     <- as.integer(Sys.getenv("MC_CORES", unset = "20"))
TOP_K_FMM   <- if (SMOKE_TEST) 50L  else 200L
RHYTHM_PVAL <- 0.01
GLOBAL_SEED <- 2025L

# Sweep grids
BETA_GRID  <- c(0.01, 0.1, 0.5, 1, 2, 5, 20, 50)
SDHR_GRID  <- c(0, 0.5, 1, 2, 4)

cat(sprintf("Mode      : %s\n", if (SMOKE_TEST) "SMOKE" else "PRODUCTION"))
cat(sprintf("NGENES    : %d\n", NGENES))
cat(sprintf("NSIMS     : %d\n", NSIMS))
cat(sprintf("N_FIXED   : %d\n", N_FIXED))
cat(sprintf("MC_CORES  : %d\n", N_CORES))
cat(sprintf("TOP_K_FMM : %d\n", TOP_K_FMM))
cat(sprintf("BETA_GRID : %s\n", paste(BETA_GRID,  collapse = ", ")))
cat(sprintf("SDHR_GRID : %s\n", paste(SDHR_GRID,  collapse = ", ")))

old_wd <- setwd("code"); source("setup.R"); setwd(old_wd)

out_dir_fig <- "output/sensitivity/figures"
out_dir_res <- "output/sensitivity/results"
dir.create(out_dir_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_res, recursive = TRUE, showWarnings = FALSE)

analysis <- CircadianAnalysisOptions(alpha = 0.05, p.adjust.method = "BH",
                                      fdr_thresholds = 0.05)

# =====================================================================
# 1. Load Baboon LUN and fit FMM per top-K rhythmic gene
# =====================================================================
cat("\n=== Step 1: anchor pilot — Baboon LUN ===\n")
load("data/CAMO_PRC_hmb.RData")
prep_lun <- prepCircadianData(baboon_withTOD$baboon[["LUN"]],
                               times = baboon_withTOD$tod[["LUN"]] %% 24,
                               input_type = "cpm")
mat_lun <- prep_lun$data[rowSums(prep_lun$data > 0) >= 6, , drop = FALSE]
tod_lun <- prep_lun$times
rm(baboon_withTOD, gtex, mice, prep_lun)

set.seed(GLOBAL_SEED)
g_idx   <- sample(nrow(mat_lun), min(NGENES, nrow(mat_lun)))
mat_sub <- mat_lun[g_idx, ]

# Cache the FMM-fitted bio object — expensive to recompute
bio_rds <- file.path(out_dir_res, "bio_lun_fmm.rds")
if (file.exists(bio_rds) && !identical(Sys.getenv("REFIT"), "true")) {
  cat(sprintf("Loading cached FMM fit: %s\n", bio_rds))
  bio_lun <- readRDS(bio_rds)
} else {
  cat("Fitting FMM per top-K gene...\n")
  bio_lun <- estCircadianParamFMM(mat_sub, tod_lun,
                                   min_rhythm_pval = RHYTHM_PVAL,
                                   top_k    = TOP_K_FMM,
                                   mc.cores = N_CORES,
                                   verbose  = TRUE)
  bio_lun$ngenes <- NGENES
  saveRDS(bio_lun, bio_rds)
}

beta_hat        <- bio_lun$diagnostics$beta_hat
sigma_alpha_hat <- bio_lun$diagnostics$sigma_alpha_hat
R2_median       <- bio_lun$diagnostics$R2_median
omega_emp       <- bio_lun$diagnostics$omega_emp
alpha_emp       <- bio_lun$diagnostics$alpha_emp

cat(sprintf("Anchor diagnostics: beta_hat=%.3f  sigma_alpha_hat=%.3fh  R2_med=%.3f  n_fitted=%d\n",
            beta_hat, sigma_alpha_hat, R2_median, bio_lun$diagnostics$n_fitted))

# =====================================================================
# 2. Power sweep over Beta(1, beta) — omega only
# =====================================================================
cat("\n=== Step 2: omega sensitivity sweep ===\n")
power_beta <- numeric(length(BETA_GRID))

for (k in seq_along(BETA_GRID)) {
  bv <- BETA_GRID[k]
  cat(sprintf("  beta = %g ... ", bv))
  bio_k <- bio_lun
  # Override pairing logic: omega_dist samples per simulation,
  # alpha stays at empirical (paired_alpha=TRUE in bio)
  bio_k$omega_rhythmic <- NULL  # let omega_dist take over
  bio_k$omega_dist     <- list(family = "beta", a = 1, b = bv)
  bio_k$alpha_dist     <- NULL   # use empirical alpha (no jitter)

  design_k <- CircadianDesignOptions(sample_sizes = N_FIXED,
                                      nsims        = NSIMS,
                                      design       = "active",
                                      cts          = seq(0, 24 * (1 - 1/12), length.out = 12),
                                      B_values     = 12L)

  set.seed(GLOBAL_SEED)
  res <- runSimsSingleCohort(bio_k, design_k, analysis,
                              mc.cores = N_CORES, verbose = FALSE)
  power_beta[k] <- mean(res$marginal_power, na.rm = TRUE)
  cat(sprintf("power = %.3f\n", power_beta[k]))
}

# =====================================================================
# 3. Power sweep over alpha noise sd_hours — alpha only
# =====================================================================
cat("\n=== Step 3: alpha sensitivity sweep ===\n")
power_sdhr <- numeric(length(SDHR_GRID))

for (k in seq_along(SDHR_GRID)) {
  sv <- SDHR_GRID[k]
  cat(sprintf("  sd_hours = %g ... ", sv))
  bio_k <- bio_lun
  bio_k$omega_rhythmic <- bio_lun$omega_rhythmic  # keep empirical omega
  bio_k$omega_dist     <- NULL
  bio_k$alpha_dist     <- list(family = "normal", mean = 0, sd_hours = sv)

  design_k <- CircadianDesignOptions(sample_sizes = N_FIXED,
                                      nsims        = NSIMS,
                                      design       = "active",
                                      cts          = seq(0, 24 * (1 - 1/12), length.out = 12),
                                      B_values     = 12L)

  set.seed(GLOBAL_SEED)
  res <- runSimsSingleCohort(bio_k, design_k, analysis,
                              mc.cores = N_CORES, verbose = FALSE)
  power_sdhr[k] <- mean(res$marginal_power, na.rm = TRUE)
  cat(sprintf("power = %.3f\n", power_sdhr[k]))
}

# =====================================================================
# 4. Save results
# =====================================================================
out <- list(
  bio_lun         = bio_lun,
  beta_grid       = BETA_GRID,
  power_beta      = power_beta,
  sdhr_grid       = SDHR_GRID,
  power_sdhr      = power_sdhr,
  N_fixed         = N_FIXED,
  beta_hat        = beta_hat,
  sigma_alpha_hat = sigma_alpha_hat,
  omega_emp       = omega_emp,
  alpha_emp       = alpha_emp,
  R2_median       = R2_median
)
ts       <- format(Sys.time(), "%Y%m%d_%H%M%S")
rds_path <- file.path(out_dir_res, sprintf("fig4_sensitivity_%s.rds", ts))
saveRDS(out, rds_path)
cat(sprintf("\nResults saved -> %s\n", rds_path))

# =====================================================================
# 5. Plot — 2x2 layout
# =====================================================================
cat("\n=== Step 4: render Fig 4 ===\n")

fig_paths <- c(
  file.path(out_dir_fig, "fig4_sensitivity.pdf"),
  "output/main_figures/Fig4_sensitivity.pdf"
)
dir.create("output/main_figures", recursive = TRUE, showWarnings = FALSE)

for (out_pdf in fig_paths) {
  pdf(out_pdf, width = 12, height = 10)
  par(mfrow = c(2, 2), mar = c(4.4, 4.4, 3, 1.2), las = 1,
      cex.lab = 1.15, cex.axis = 1.0, cex.main = 1.1, font.main = 2)

  # --- (a) empirical omega histogram + Beta(1, beta_hat) overlay ---
  hist(omega_emp, breaks = seq(0, 1, by = 0.05), freq = FALSE,
       main = "(A) Empirical omega_hat from FMM fits",
       xlab = expression(hat(omega)), col = "lightsteelblue",
       border = "white", xlim = c(0, 1))
  curve(dbeta(x, 1, beta_hat), add = TRUE, col = "darkred", lwd = 2)
  abline(v = 1, lty = 3, col = "grey50")
  legend("topright", bty = "n", lwd = 2, col = "darkred",
         legend = sprintf("Beta(1, %.2f) MoM fit", beta_hat))

  # --- (b) empirical alpha histogram + von Mises overlay (kappa = 1/sd_rad^2) ---
  sd_rad <- sigma_alpha_hat * (2 * pi / 24)
  kappa  <- 1 / sd_rad^2
  hist(alpha_emp, breaks = 18, freq = FALSE,
       main = "(B) Empirical alpha_hat from FMM fits",
       xlab = expression(hat(alpha) ~ "(rad)"),
       col = "lightsalmon", border = "white", xlim = c(0, 2 * pi))
  # von Mises density at the mean of empirical alpha
  mu_alpha <- mean(alpha_emp)
  vm_dens  <- function(x) exp(kappa * cos(x - mu_alpha)) /
                          (2 * pi * besselI(kappa, 0))
  curve(vm_dens, add = TRUE, col = "darkred", lwd = 2)
  legend("topright", bty = "n", lwd = 2, col = "darkred",
         legend = sprintf("vonMises(mu=%.2f, sd_h=%.2f)", mu_alpha, sigma_alpha_hat))

  # --- (c) power vs beta (omega sweep) ---
  plot(BETA_GRID, power_beta, type = "o", pch = 19, lwd = 2.2,
       col = "steelblue", log = "x",
       ylim = c(0, 1),
       xlab = expression(beta ~ "of Beta(1," ~ beta ~ ")"),
       ylab = sprintf("Power at N=%d", N_FIXED),
       main = "(C) omega sensitivity")
  abline(v = beta_hat, lty = 2, col = "darkred", lwd = 1.5)
  abline(h = 0.80, lty = 3, col = "grey50")
  text(beta_hat, 0.05, sprintf("hat(beta)=%.2f", beta_hat),
       pos = 4, cex = 0.9, col = "darkred")

  # --- (d) power vs sd_hours (alpha sweep) ---
  plot(SDHR_GRID, power_sdhr, type = "o", pch = 19, lwd = 2.2,
       col = "steelblue",
       ylim = c(0, 1),
       xlab = "Phase noise sd (hours)",
       ylab = sprintf("Power at N=%d", N_FIXED),
       main = "(D) alpha sensitivity")
  abline(v = sigma_alpha_hat, lty = 2, col = "darkred", lwd = 1.5)
  abline(h = 0.80, lty = 3, col = "grey50")
  text(sigma_alpha_hat, 0.05, sprintf("hat(sigma)[alpha]=%.2fh", sigma_alpha_hat),
       pos = 4, cex = 0.9, col = "darkred")

  mtext("Fig 4: Sensitivity to cosinor violation (Baboon LUN, n=12, active)",
        outer = TRUE, line = -1.5, cex = 1.05, font = 2)

  dev.off()
  cat(sprintf("Saved: %s\n", out_pdf))
}

cat("\n=== Done ===\n")
