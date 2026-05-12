#' =======================================================================
#' supp_fmm_diagnostic.R — Q1 evidence: does data violate cosinor?
#' =======================================================================
#'
#' Multi-tissue diagnostic figure showing empirical FMM-fit distributions
#' across active-design tissues. Establishes that real data does NOT obey
#' the cosinor assumption (ω = 1) — motivating the §2.4 sensitivity
#' analysis and §2.5 active-vs-passive comparison.
#'
#' For each tissue:
#'   - Histogram of empirical ω̂_g (per-gene FMM shape parameter, top-K)
#'   - Beta(1, β̂) MoM overlay
#'   - Histogram of empirical α̂_g (peak location, radians)
#'   - vonMises(α̂_med, kappa) overlay
#'
#' Plus a summary table: tissue, n, n_fitted, β̂, σ̂_α (hours), R²-median.
#'
#' USAGE:
#'   Rscript examples/publication/supp_fmm_diagnostic.R
#'   SMOKE_TEST=true Rscript examples/publication/supp_fmm_diagnostic.R

SMOKE_TEST  <- identical(Sys.getenv("SMOKE_TEST"), "true")
NGENES_FIT  <- if (SMOKE_TEST) 500L  else 5000L
TOP_K_FMM   <- if (SMOKE_TEST) 50L   else 200L
N_CORES     <- as.integer(Sys.getenv("MC_CORES", unset = "20"))
RHYTHM_PVAL <- 0.01
GLOBAL_SEED <- 2025L

# Six tissues spanning a range of cosinor violation severity.
# Per the prior FMM scan: KIC strongest violation (median ω≈0.33),
# ADC mildest (median ω≈0.58). All Baboon CAMO active-design.
tissues <- c("KIC", "LUN", "ADM", "LIV", "ADC", "PAN")

cat(sprintf("Mode      : %s\n", if (SMOKE_TEST) "SMOKE" else "PRODUCTION"))
cat(sprintf("NGENES_FIT: %d\n", NGENES_FIT))
cat(sprintf("TOP_K_FMM : %d\n", TOP_K_FMM))
cat(sprintf("MC_CORES  : %d\n", N_CORES))
cat(sprintf("TISSUES   : %s\n", paste(tissues, collapse = ", ")))

old_wd <- setwd("code"); source("setup.R"); setwd(old_wd)

out_dir_fig <- "output/supplementary/figures"
out_dir_res <- "output/supplementary/results"
dir.create(out_dir_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_res, recursive = TRUE, showWarnings = FALSE)

# =====================================================================
# 1. Fit FMM per tissue (cache-aware)
# =====================================================================
cat("\n=== Step 1: per-tissue FMM fits ===\n")
load("data/CAMO_PRC_hmb.RData")

bio_list <- list()

for (t in tissues) {
  cat(sprintf("\n--- %s ---\n", t))
  m_full <- baboon_withTOD$baboon[[t]]
  if (is.null(m_full)) { cat(sprintf("  SKIP %s — not in CAMO\n", t)); next }
  prep <- prepCircadianData(m_full,
                             times = baboon_withTOD$tod[[t]] %% 24,
                             input_type = "cpm")
  mat  <- prep$data[rowSums(prep$data > 0) >= 6, , drop = FALSE]
  tod  <- prep$times
  set.seed(GLOBAL_SEED)
  g_idx <- sample(nrow(mat), min(NGENES_FIT, nrow(mat)))
  mat_sub <- mat[g_idx, ]

  # Reuse cached fit from Fig 4/5 if available, otherwise refit
  cached_path <- file.path("output", "active_vs_passive", "results",
                            sprintf("bio_baboon_%s_fmm.rds", t))
  if (file.exists(cached_path) && !identical(Sys.getenv("REFIT"), "true")) {
    cat(sprintf("Loading cached: %s\n", cached_path))
    bio_t <- readRDS(cached_path)
  } else {
    bio_t <- estCircadianParamFMM(mat_sub, tod,
                                   min_rhythm_pval = RHYTHM_PVAL,
                                   top_k    = TOP_K_FMM,
                                   mc.cores = N_CORES,
                                   verbose  = TRUE)
    bio_t$ngenes  <- NGENES_FIT
    bio_t$.t_name <- t
    bio_t$.n_pilot <- ncol(prep$data)
    saveRDS(bio_t, file.path(out_dir_res, sprintf("bio_baboon_%s_supp.rds", t)))
  }
  bio_t$.t_name  <- t
  bio_t$.n_pilot <- ncol(prep$data)
  bio_list[[t]] <- bio_t
}
rm(baboon_withTOD, gtex, mice)

# =====================================================================
# 2. Summary table
# =====================================================================
cat("\n=== Diagnostic summary ===\n")
summary_df <- data.frame(
  tissue          = character(0),
  n_pilot         = integer(0),
  n_fitted        = integer(0),
  prop_rhythmic   = numeric(0),
  beta_hat        = numeric(0),
  omega_median    = numeric(0),
  sigma_alpha_hat = numeric(0),
  R2_median       = numeric(0),
  stringsAsFactors = FALSE
)
for (t in names(bio_list)) {
  b <- bio_list[[t]]; d <- b$diagnostics
  summary_df <- rbind(summary_df, data.frame(
    tissue          = t,
    n_pilot         = b$.n_pilot,
    n_fitted        = d$n_fitted,
    prop_rhythmic   = b$prop_rhythmic,
    beta_hat        = d$beta_hat,
    omega_median    = median(d$omega_emp),
    sigma_alpha_hat = d$sigma_alpha_hat,
    R2_median       = d$R2_median
  ))
}
print(summary_df, row.names = FALSE, digits = 3)

ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
saveRDS(list(bio_list = bio_list, summary = summary_df),
        file.path(out_dir_res, sprintf("supp_fmm_diagnostic_%s.rds", ts)))
write.csv(summary_df,
          file.path(out_dir_res, sprintf("supp_fmm_diagnostic_summary_%s.csv", ts)),
          row.names = FALSE)

# =====================================================================
# 3. Plot — one row per tissue, two columns (omega, alpha)
# =====================================================================
cat("\n=== Step 2: render diagnostic figure ===\n")
source("examples/publication/_pub_style.R")

n_t <- length(bio_list)
fig_paths <- c(
  file.path(out_dir_fig, "supp_fmm_diagnostic.pdf"),
  "output/main_figures/SuppFig_FMM_diagnostic.pdf"
)
dir.create("output/main_figures", recursive = TRUE, showWarnings = FALSE)

col_omega   <- "#B0C4DE"   # lightsteelblue, soft
col_alpha   <- "#FAA582"   # soft salmon
col_fit     <- "#9E2A2B"   # muted dark red
col_median  <- "#0072B2"   # Wong blue, matches detector palette

tissue_names <- names(bio_list)
letters_seq  <- letters[seq_len(2 * length(tissue_names))]

for (out_pdf in fig_paths) {
  cairo_pdf(out_pdf, width = 7.2, height = 1.9 * n_t + 0.4)
  pub_par(mfrow = c(n_t, 2), mar = c(3.6, 4.0, 2.4, 0.8))

  for (ti in seq_along(tissue_names)) {
    t <- tissue_names[ti]
    b <- bio_list[[t]]; d <- b$diagnostics
    omega_emp <- d$omega_emp
    alpha_emp <- d$alpha_emp
    beta_hat  <- d$beta_hat
    sa_hat    <- d$sigma_alpha_hat

    # --- omega panel ---
    hist(omega_emp, breaks = seq(0, 1, by = 0.05), freq = FALSE,
         main = sprintf("%s, omega", t),
         xlab = expression(hat(omega)),
         col = col_omega, border = "white", xlim = c(0, 1))
    panel_label(letters_seq[2L * ti - 1L])
    curve(dbeta(x, 1, beta_hat), add = TRUE, col = col_fit, lwd = 1.8)
    abline(v = 1, lty = 3, col = "grey60")
    abline(v = median(omega_emp), lty = 2, col = col_median, lwd = 1.2)
    pub_legend("topright",
               legend = c(sprintf("median = %.2f", median(omega_emp)),
                          expression(Beta(1, hat(beta)))),
               col = c(col_median, col_fit),
               lty = c(2, 1), lwd = c(1.2, 1.8))

    # --- alpha panel ---
    sd_rad <- sa_hat * (2 * pi / 24)
    kappa  <- if (sd_rad > 0) 1 / sd_rad^2 else Inf
    mu_a   <- mean(alpha_emp)

    hist(alpha_emp, breaks = 18, freq = FALSE,
         main = sprintf("%s, alpha", t),
         xlab = expression(hat(alpha) ~ "(rad)"),
         col = col_alpha, border = "white", xlim = c(0, 2 * pi))
    panel_label(letters_seq[2L * ti])
    if (is.finite(kappa)) {
      vm_dens <- function(x) exp(kappa * cos(x - mu_a)) /
                              (2 * pi * besselI(kappa, 0))
      curve(vm_dens, add = TRUE, col = col_fit, lwd = 1.8)
    }
    pub_legend("topright",
               legend = c(sprintf("mean = %.2f rad", mu_a),
                          expression(vonMises(hat(mu), hat(kappa)))),
               col = c(col_fit, col_fit),
               lty = c(NA, 1), lwd = c(NA, 1.8),
               pch = c(NA, NA))
  }

  dev.off()
  cat(sprintf("Saved: %s\n", out_pdf))
}

cat("\n=== Done ===\n")
