#' =======================================================================
#' fig4_sensitivity.R — Sensitivity to cosinor violation (FMM (K=2))
#' =======================================================================
#'
#' Two orthogonal sensitivity sweeps anchored on Baboon LUN. Both panels
#' use FMM (K=2) detection and parametric distributions for ω and α.
#'
#' Panel A — ω sweep (α distribution fixed at empirical):
#'   ω_g ~ Beta(1, β),    sweep β over a log grid
#'   x-axis: N (total samples), one curve per β value
#'   β̂ from MoM fit on empirical ω̂_g marked as a curve label
#'
#' Panel B — α sweep (ω distribution fixed at empirical Beta(1, β̂)):
#'   α_g ~ vonMises(0, κ),  κ = 1/sd_rad²,  sd_rad = sd_hours · 2π/24
#'   x-axis: N (total samples), one curve per σ_α (hours) value
#'   σ̂_α from sd of empirical α̂_g marked as a curve label
#'   Expected: curves overlap (rotation-invariance of cosinor / FMM (K=2))
#'
#' Detector: detect_FMM (K=2 harmonic LRT) only.  DCP behavior is closed-form
#' (NCP = (r·c(ω))²·N/2) and would just shift the whole sweep down with
#' attenuation c(ω)² — uninformative. The interesting question is whether
#' the K-harmonic LRT preserves power under non-cosinor truth.
#'
#' USAGE:
#'   Rscript examples/publication/fig4_sensitivity.R
#'   SMOKE_TEST=true Rscript examples/publication/fig4_sensitivity.R

SMOKE_TEST  <- identical(Sys.getenv("SMOKE_TEST"), "true")
NGENES      <- if (SMOKE_TEST) 500L  else 2000L
NSIMS       <- if (SMOKE_TEST) 5L    else 20L
N_CORES     <- as.integer(Sys.getenv("MC_CORES", unset = "20"))
TOP_K_FMM   <- if (SMOKE_TEST) 50L  else 300L  # codebase default min(300, n_cand)
RHYTHM_PVAL <- 0.01
GLOBAL_SEED <- 2025L

# x-axis: sample size grid
# N must be divisible by B_VAL=12 (active design with 12 timepoints).
N_GRID  <- if (SMOKE_TEST) c(24L, 36L) else c(12L, 24L, 36L, 48L, 60L, 72L, 96L, 120L, 144L)

# Sweep grids
BETA_GRID <- if (SMOKE_TEST) c(0.5, 1, 5)        else c(0.5, 1, 2, 5, 20)
SDHR_GRID <- if (SMOKE_TEST) c(0, 1, 4)          else c(0, 0.5, 1, 2, 4)

cat(sprintf("Mode      : %s\n", if (SMOKE_TEST) "SMOKE" else "PRODUCTION"))
cat(sprintf("NGENES    : %d\n", NGENES))
cat(sprintf("NSIMS     : %d\n", NSIMS))
cat(sprintf("N_GRID    : %s\n", paste(N_GRID,    collapse = ", ")))
cat(sprintf("BETA_GRID : %s\n", paste(BETA_GRID, collapse = ", ")))
cat(sprintf("SDHR_GRID : %s\n", paste(SDHR_GRID, collapse = ", ")))
cat(sprintf("TOP_K_FMM : %d\n", TOP_K_FMM))
cat(sprintf("MC_CORES  : %d\n", N_CORES))

old_wd <- setwd("code"); source("setup.R"); setwd(old_wd)

out_dir_fig <- "output/sensitivity/figures"
out_dir_res <- "output/sensitivity/results"
dir.create(out_dir_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_res, recursive = TRUE, showWarnings = FALSE)

analysis <- CircadianAnalysisOptions(alpha           = 0.05,
                                      p.adjust.method = "BH",
                                      fdr_thresholds  = 0.05)

# ====================================================================
# 1. Anchor pilot: Baboon LUN
# ====================================================================
cat("\n=== Step 1: anchor pilot — Baboon LUN ===\n")
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

bio_rds <- file.path(out_dir_res,
                     sprintf("bio_lun_fmm_NG%d_K%d.rds", NGENES, TOP_K_FMM))
if (file.exists(bio_rds) && !identical(Sys.getenv("REFIT"), "true")) {
  cat(sprintf("Loading cached FMM fit: %s\n", bio_rds))
  bio_lun <- readRDS(bio_rds)
} else {
  cat("Fitting FMM per top-K rhythmic gene...\n")
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
cat(sprintf("Anchor: β̂=%.3f  σ̂_α=%.3fh  R²=%.3f  n_fitted=%d\n",
            beta_hat, sigma_alpha_hat, bio_lun$diagnostics$R2_median,
            bio_lun$diagnostics$n_fitted))

# Active design: equispaced 12 timepoints — m = N/12 reps each
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
                                method   = "FMM", K = 2L,
                                mc.cores = N_CORES, verbose = FALSE)
    mean(res$marginal_power, na.rm = TRUE)
  }, numeric(1))
  cat(sprintf("  %-22s : %s\n", label,
              paste(sprintf("%.2f", pwr), collapse = " ")))
  pwr
}

# ====================================================================
# 2. Panel A — ω sweep (α distribution fixed at empirical)
# ====================================================================
cat("\n=== Step 2: ω sweep (FMM (K=2)) ===\n")
power_A <- matrix(NA_real_, nrow = length(N_GRID), ncol = length(BETA_GRID),
                  dimnames = list(paste0("N=", N_GRID),
                                  paste0("beta=", BETA_GRID)))
for (k in seq_along(BETA_GRID)) {
  bio_k <- bio_lun
  bio_k$omega_rhythmic <- NULL
  bio_k$omega_dist     <- list(family = "beta", a = 1, b = BETA_GRID[k])
  bio_k$alpha_rhythmic <- NULL
  bio_k$alpha_dist     <- list(family = "normal", mean = 0,
                                sd_hours = sigma_alpha_hat)  # fixed at empirical
  power_A[, k] <- run_sweep_curve(bio_k, sprintf("β=%.2f", BETA_GRID[k]))
}

# ====================================================================
# 3. Panel B — α sweep (ω distribution fixed at empirical Beta(1, β̂))
# ====================================================================
cat("\n=== Step 3: α sweep (FMM (K=2)) ===\n")
power_B <- matrix(NA_real_, nrow = length(N_GRID), ncol = length(SDHR_GRID),
                  dimnames = list(paste0("N=", N_GRID),
                                  paste0("sdh=", SDHR_GRID)))
for (k in seq_along(SDHR_GRID)) {
  bio_k <- bio_lun
  bio_k$omega_rhythmic <- NULL
  bio_k$omega_dist     <- list(family = "beta", a = 1, b = beta_hat)  # fixed at empirical
  bio_k$alpha_rhythmic <- NULL
  bio_k$alpha_dist     <- list(family = "normal", mean = 0,
                                sd_hours = SDHR_GRID[k])
  power_B[, k] <- run_sweep_curve(bio_k, sprintf("σ_α=%.1fh", SDHR_GRID[k]))
}

# ====================================================================
# 4. Save results
# ====================================================================
out <- list(
  bio_lun         = bio_lun,
  beta_hat        = beta_hat,
  sigma_alpha_hat = sigma_alpha_hat,
  N_grid          = N_GRID,
  beta_grid       = BETA_GRID,
  sdhr_grid       = SDHR_GRID,
  power_A         = power_A,
  power_B         = power_B
)
ts       <- format(Sys.time(), "%Y%m%d_%H%M%S")
rds_path <- file.path(out_dir_res, sprintf("fig4_sensitivity_%s.rds", ts))
saveRDS(out, rds_path)
cat(sprintf("\nSaved: %s\n", rds_path))

# ====================================================================
# 5. Plot — 1×2 layout
# ====================================================================
cat("\n=== Step 4: render Fig 4 ===\n")
source("examples/publication/_pub_style.R")

fig_paths <- c(
  file.path(out_dir_fig, "fig4_sensitivity.pdf"),
  "output/main_figures/Fig4_sensitivity.pdf"
)
dir.create("output/main_figures", recursive = TRUE, showWarnings = FALSE)

for (out_pdf in fig_paths) {
  cairo_pdf(out_pdf, width = 7.2, height = 3.4)
  pub_par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.4, 1.0))

  pal_A <- pub_palette_sequential(length(BETA_GRID))
  pal_B <- pub_palette_sequential(length(SDHR_GRID))

  # --- Panel a: omega sweep ---
  plot(NA, xlim = range(N_GRID), ylim = c(0, 1),
       xlab = "N (total samples)", ylab = "Power",
       main = expression(omega ~ "sweep"))
  panel_label("a")
  abline_80pct()
  for (k in seq_along(BETA_GRID)) {
    lty_k <- if (abs(BETA_GRID[k] - beta_hat) < 0.05) 1 else 2
    lwd_k <- if (abs(BETA_GRID[k] - beta_hat) < 0.05) 2.2 else 1.5
    lines(N_GRID, power_A[, k], col = pal_A[k], lwd = lwd_k, lty = lty_k)
    points(N_GRID, power_A[, k], pch = 19, col = pal_A[k], cex = 0.7)
  }
  pub_legend("bottomright",
             legend = sprintf("%g", BETA_GRID),
             col = pal_A, lwd = 1.6,
             title = expression(beta))

  # --- Panel b: alpha sweep ---
  plot(NA, xlim = range(N_GRID), ylim = c(0, 1),
       xlab = "N (total samples)", ylab = "Power",
       main = expression(alpha ~ "sweep"))
  panel_label("b")
  abline_80pct()
  for (k in seq_along(SDHR_GRID)) {
    lines(N_GRID, power_B[, k], col = pal_B[k], lwd = 1.6)
    points(N_GRID, power_B[, k], pch = 19, col = pal_B[k], cex = 0.7)
  }
  pub_legend("bottomright",
             legend = sprintf("%g", SDHR_GRID),
             col = pal_B, lwd = 1.6,
             title = expression(sigma[alpha] ~ "(h)"))

  dev.off()
  cat(sprintf("Saved: %s\n", out_pdf))
}

cat("\n=== Done ===\n")
