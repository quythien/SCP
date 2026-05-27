#' =======================================================================
#' fig4_mh_sensitivity.R — Sensitivity of cosinor-based SCP (DCP)
#'                          to MH-shaped data (multi-harmonic truth)
#' =======================================================================
#'
#' Same scientific question as fig4_cosinor_sensitivity.R, but the
#' data-generating model is the multi-harmonic (MH) cosine series
#' instead of the Frequency-Modulated Möbius (FMM):
#'
#'   y_g(t) = mu_g + A_g [ cos(omega0 (t - phi_g))
#'                       + alpha2 * cos(2 omega0 (t - phi_g)) ]
#'                       + eps_g(t)
#'
#' Pipeline:
#'   (i)   pilot params via estCircadianParam (cosinor fits, NOT FMM)
#'   (ii)  data generated under MH truth via the harmonics = c(alpha2, 0)
#'         argument to runSimsSingleCohort (alpha3 fixed at 0 per
#'         "K up to 2" scope)
#'   (iii) detection via DCP (cosinor F-test, method = "DCP")
#'   (iv)  power evaluated at BH-FDR alpha = 0.05
#'
#' Panel A — alpha2 sweep (phase distribution fixed at empirical):
#'   alpha2 in {0, 0.1, 0.2, 0.3, 0.5}; alpha3 = 0; phase_g empirical.
#'   One curve per alpha2 value.
#'
#' Panel B — phase-dispersion sweep (alpha2 fixed at anchor):
#'   phi_g ~ vonMises(0, kappa), kappa = (24 / (2*pi*sigma_phi))^2,
#'   sigma_phi in {0, 0.5, 1, 2, 4} h; alpha2 = ALPHA2_ANCHOR.
#'   Includes sigma_phi = 0 (all genes peak at same time) baseline.
#'
#' ALPHA2_ANCHOR represents a "moderately non-cosinor" empirical anchor.
#' Default 0.3 (second-harmonic amplitude is 30% of first), comparable
#' to the LUN FMM eta_hat = 1.80 regime.
#'
#' OUTPUT:
#'   output/mh_alternative/results/fig4_mh_sensitivity_<ts>.rds
#'
#' USAGE:
#'   SMOKE_TEST=true Rscript examples/publication/mh_alternative/fig4_mh_sensitivity.R
#'   Rscript examples/publication/mh_alternative/fig4_mh_sensitivity.R

SMOKE_TEST  <- identical(Sys.getenv("SMOKE_TEST"), "true")
NGENES      <- if (SMOKE_TEST) 500L  else 2000L
NSIMS       <- if (SMOKE_TEST) 5L    else 20L
N_CORES     <- as.integer(Sys.getenv("MC_CORES", unset = "20"))
RHYTHM_PVAL <- 0.01
GLOBAL_SEED <- 2025L

N_GRID    <- if (SMOKE_TEST) c(24L, 36L) else
             c(12L, 24L, 36L, 48L, 60L, 72L, 96L, 120L, 144L)

# Three sweep axes, parallel to the FMM Fig 4 sweeps but with per-gene
# alpha2 / alpha3 drawn from Beta(1, eta) distributions (FMM-eta-style).
# Beta(1, eta): small eta concentrates mass near 1 (heavy non-cosinor);
# large eta concentrates mass near 0 (cosinor-like). E.g.,
#   eta = 0.5 -> mean(alpha) = 0.67 (mostly violated)
#   eta = 1   -> Uniform(0, 1)
#   eta = 20  -> mean(alpha) = 0.05 (mostly cosinor)
ETA_A2_GRID <- if (SMOKE_TEST) c(0.5, 2, 20) else c(0.5, 1, 2, 5, 20)
ETA_A3_GRID <- if (SMOKE_TEST) c(0.5, 2, 20) else c(0.5, 1, 2, 5, 20)
SDHR_GRID   <- if (SMOKE_TEST) c(0, 1, 4)    else c(0, 0.5, 1, 2, 4)

# Anchor eta corresponds to mean(alpha) approx 0.3, matching the prior
# scalar anchor: 0.3 = 1/(1 + eta) -> eta_anchor approx 2.33.
ETA_A2_ANCHOR <- 2.33
ETA_A3_ANCHOR <- 20      # near-cosinor for the 3rd harmonic by default

cat(sprintf("Mode          : %s\n", if (SMOKE_TEST) "SMOKE" else "PRODUCTION"))
cat(sprintf("NGENES        : %d\n", NGENES))
cat(sprintf("NSIMS         : %d\n", NSIMS))
cat(sprintf("N_GRID        : %s\n", paste(N_GRID,        collapse = ", ")))
cat(sprintf("ETA_A2_GRID   : %s\n", paste(ETA_A2_GRID,   collapse = ", ")))
cat(sprintf("ETA_A3_GRID   : %s\n", paste(ETA_A3_GRID,   collapse = ", ")))
cat(sprintf("SDHR_GRID     : %s\n", paste(SDHR_GRID,     collapse = ", ")))
cat(sprintf("ETA_A2_ANCHOR : %g  (mean alpha2 = %.3f)\n",
            ETA_A2_ANCHOR, 1 / (1 + ETA_A2_ANCHOR)))
cat(sprintf("ETA_A3_ANCHOR : %g  (mean alpha3 = %.3f)\n",
            ETA_A3_ANCHOR, 1 / (1 + ETA_A3_ANCHOR)))

setwd("/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")
old_wd <- setwd("code"); source("setup.R"); setwd(old_wd)

out_dir_res <- "output/mh_alternative/results"
dir.create(out_dir_res, recursive = TRUE, showWarnings = FALSE)

analysis <- CircadianAnalysisOptions(alpha           = 0.05,
                                      p.adjust.method = "BH",
                                      fdr_thresholds  = 0.05)

# ====================================================================
# 1. Anchor pilot: Baboon LUN — cosinor estimation
# ====================================================================
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

cat("Step 1a: cosinor-based pilot estimation (estCircadianParam)...\n")
bio_lun <- estCircadianParam(mat_sub, tod_lun,
                              min_rhythm_pval = RHYTHM_PVAL,
                              verbose = TRUE)
bio_lun$ngenes <- NGENES
# Ensure NO per-gene FMM params -> cosinor + harmonics path will trigger.
bio_lun$omega_rhythmic <- NULL
bio_lun$alpha_rhythmic <- NULL
bio_lun$omega_dist     <- NULL
bio_lun$alpha_dist     <- NULL

# Active design B = 12
B_VAL <- 12L
cts_active <- seq(0, 24 * (1 - 1/B_VAL), length.out = B_VAL)

# ====================================================================
# 2. Power-sweep helper (returns marginal power per N for a given harmonic vector)
# ====================================================================
run_sweep_curve <- function(bio_k, label, harmonics = c(0, 0)) {
  pwr <- vapply(N_GRID, function(N) {
    if (N %% B_VAL != 0L) return(NA_real_)
    m <- N %/% B_VAL
    cts_design <- rep(cts_active, each = m)
    design <- CircadianDesignOptions(sample_sizes = N, nsims = NSIMS,
                                      design = "active",
                                      cts = cts_design, B_values = B_VAL)
    set.seed(GLOBAL_SEED)
    res <- runSimsSingleCohort(bio_k, design, analysis,
                                method   = "DCP",
                                harmonics = harmonics,
                                mc.cores = N_CORES, verbose = FALSE)
    mean(res$marginal_power, na.rm = TRUE)
  }, numeric(1))
  cat(sprintf("  %-30s : %s\n", label,
              paste(sprintf("%.2f", pwr), collapse = " ")))
  pwr
}

# Helper to set per-gene alpha2/alpha3 Beta(1, eta) distributions on bio.opts.
set_alpha_dist <- function(bio_base, eta_a2 = NULL, eta_a3 = NULL) {
  bio_k <- bio_base
  bio_k$alpha2_dist <- if (!is.null(eta_a2))
    list(family = "beta", a = 1, b = eta_a2) else NULL
  bio_k$alpha3_dist <- if (!is.null(eta_a3))
    list(family = "beta", a = 1, b = eta_a3) else NULL
  bio_k
}

# Helper: set a vonMises-like phase distribution on bio.opts.
# The simulator draws phi_g from bio.opts$phase (in hours). For
# sigma_phi = 0 all rhythmic genes share the same phase (mean of pilot);
# otherwise we sample around the mean with the requested sd.
set_phase_dist <- function(bio_base, sd_hours, n_draw = 5000L) {
  bio_k <- bio_base
  phi_bar <- mean(bio_base$phase, na.rm = TRUE)
  if (sd_hours <= 0) {
    bio_k$phase <- rep(phi_bar, n_draw)
  } else {
    sd_rad <- sd_hours * 2 * pi / 24
    kappa  <- 1 / sd_rad^2
    .rvonmises <- function(n, mu, kappa) {
      # Best-rejection sampler (Wood 1994); same as simulation.R's internal.
      a <- 1 + sqrt(1 + 4 * kappa^2)
      b <- (a - sqrt(2 * a)) / (2 * kappa)
      r <- (1 + b^2) / (2 * b)
      out <- numeric(n)
      i <- 1
      while (i <= n) {
        u1 <- runif(1); z <- cos(pi * u1)
        f  <- (1 + r * z) / (r + z)
        c  <- kappa * (r - f)
        u2 <- runif(1)
        if (c * (2 - c) > u2 || log(c / u2) + 1 - c >= 0) {
          u3 <- runif(1)
          out[i] <- (mu + sign(u3 - 0.5) * acos(f)) %% (2 * pi)
          i <- i + 1
        }
      }
      out
    }
    set.seed(GLOBAL_SEED + 17L)
    phi_rad <- .rvonmises(n_draw, mu = phi_bar * 2 * pi / 24, kappa = kappa)
    bio_k$phase <- (phi_rad * 24 / (2 * pi)) %% 24
  }
  bio_k
}

# ====================================================================
# 3. Panel D -- eta_alpha2 sweep:  alpha2_g ~ Beta(1, eta_a2),  alpha3 = 0
# ====================================================================
cat("\n=== Step 2: eta_alpha2 sweep (DCP under MH truth, per-gene alpha2 ~ Beta(1, eta_a2)) ===\n")
power_A <- matrix(NA_real_, nrow = length(N_GRID), ncol = length(ETA_A2_GRID),
                  dimnames = list(paste0("N=", N_GRID),
                                  paste0("eta_a2=", ETA_A2_GRID)))
for (k in seq_along(ETA_A2_GRID)) {
  bio_k <- set_alpha_dist(bio_lun, eta_a2 = ETA_A2_GRID[k], eta_a3 = NULL)
  power_A[, k] <- run_sweep_curve(bio_k,
                                  label = sprintf("eta_a2=%.2f", ETA_A2_GRID[k]))
}

# ====================================================================
# 4. Panel E -- eta_alpha3 sweep:  alpha3_g ~ Beta(1, eta_a3),  alpha2_g ~ anchor dist
# ====================================================================
cat(sprintf("\n=== Step 3: eta_alpha3 sweep (DCP under MH truth, alpha2_g ~ Beta(1, %g)) ===\n",
            ETA_A2_ANCHOR))
power_B <- matrix(NA_real_, nrow = length(N_GRID), ncol = length(ETA_A3_GRID),
                  dimnames = list(paste0("N=", N_GRID),
                                  paste0("eta_a3=", ETA_A3_GRID)))
for (k in seq_along(ETA_A3_GRID)) {
  bio_k <- set_alpha_dist(bio_lun, eta_a2 = ETA_A2_ANCHOR,
                                    eta_a3 = ETA_A3_GRID[k])
  power_B[, k] <- run_sweep_curve(bio_k,
                                  label = sprintf("eta_a3=%.2f", ETA_A3_GRID[k]))
}

# ====================================================================
# 5. Panel F -- phase-dispersion sweep (alpha2_g ~ anchor dist, alpha3 = 0)
# ====================================================================
cat(sprintf("\n=== Step 4: phase-dispersion sweep (DCP under MH truth, alpha2_g ~ Beta(1, %g)) ===\n",
            ETA_A2_ANCHOR))
power_C <- matrix(NA_real_, nrow = length(N_GRID), ncol = length(SDHR_GRID),
                  dimnames = list(paste0("N=", N_GRID),
                                  paste0("sdh=", SDHR_GRID)))
for (k in seq_along(SDHR_GRID)) {
  bio_k <- set_phase_dist(bio_lun, sd_hours = SDHR_GRID[k])
  bio_k <- set_alpha_dist(bio_k, eta_a2 = ETA_A2_ANCHOR, eta_a3 = NULL)
  power_C[, k] <- run_sweep_curve(bio_k,
                                  label = sprintf("sd_phi=%.1fh", SDHR_GRID[k]))
}

# ====================================================================
# 5. Save
# ====================================================================
out <- list(
  framework       = "MH alternative -- per-gene Beta(1, eta) alpha2/alpha3 sweeps",
  detector        = "DCP (cosinor F-test, K = 1)",
  pilot           = "baboon LUN (active, GSE98965)",
  GLOBAL_SEED     = GLOBAL_SEED,
  NGENES          = NGENES, NSIMS = NSIMS,
  N_grid          = N_GRID,
  eta_a2_grid     = ETA_A2_GRID,
  eta_a3_grid     = ETA_A3_GRID,
  sdhr_grid       = SDHR_GRID,
  eta_a2_anchor   = ETA_A2_ANCHOR,
  eta_a3_anchor   = ETA_A3_ANCHOR,
  power_A         = power_A,         # eta_alpha2 sweep
  power_B         = power_B,         # eta_alpha3 sweep
  power_C         = power_C          # phi-dispersion sweep
)
ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
rds <- file.path(out_dir_res, sprintf("fig4_mh_sensitivity_%s.rds", ts))
saveRDS(out, rds)
cat(sprintf("\nSaved: %s\n", rds))
