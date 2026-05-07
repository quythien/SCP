#' =======================================================================
#' fig5_active_vs_passive.R — B vs m trade-off under cosinor vs FMM truth
#' =======================================================================
#'
#' Two-row figure comparing the (B, m) sampling design trade-off between
#' two truth conditions on three baboon CAMO active-design tissues:
#'
#'   Row 1 — Cosinor truth: omega = 1, alpha = 0 forced. Power should be
#'           approximately invariant to B at fixed N (current paper claim).
#'
#'   Row 2 — FMM truth: per-gene empirical (omega_hat, alpha_hat) from
#'           FMM fits to the pilot data. Same (A, sigma, phi) as Row 1.
#'           Hypothesis: under non-sinusoidal waveforms, more time-points
#'           (higher B) better captures the asymmetric peak, so power
#'           should INCREASE with B.
#'
#' Layout: 2 rows × 3 columns
#'   Cols: Baboon LIV (R²=0.77), Baboon LUN (R²=0.85), Baboon KIC (R²=0.90)
#'         All three confirmed cosinor-violating (median omega 0.33-0.53).
#'
#' USAGE:
#'   Rscript examples/publication/fig5_active_vs_passive.R
#'   SMOKE_TEST=true Rscript examples/publication/fig5_active_vs_passive.R

SMOKE_TEST  <- identical(Sys.getenv("SMOKE_TEST"), "true")
NGENES      <- if (SMOKE_TEST) 500L  else 5000L
NSIMS       <- if (SMOKE_TEST) 10L   else 50L
N_CORES     <- as.integer(Sys.getenv("MC_CORES", unset = "20"))
TOP_K_FMM   <- if (SMOKE_TEST) 50L   else 200L
RHYTHM_PVAL <- 0.01
GLOBAL_SEED <- 2025L

# B-vs-m grid — fix N = B*m at canonical values
N_TOTAL     <- if (SMOKE_TEST) c(24L) else c(24L, 36L, 48L)
B_GRID      <- c(2L, 3L, 4L, 6L, 8L, 12L)

# Tissue list (Baboon CAMO active-design)
tissues <- c("LIV", "LUN", "KIC")

cat(sprintf("Mode      : %s\n", if (SMOKE_TEST) "SMOKE" else "PRODUCTION"))
cat(sprintf("NGENES    : %d\n", NGENES))
cat(sprintf("NSIMS     : %d\n", NSIMS))
cat(sprintf("MC_CORES  : %d\n", N_CORES))
cat(sprintf("TOP_K_FMM : %d\n", TOP_K_FMM))
cat(sprintf("N_TOTAL   : %s\n", paste(N_TOTAL, collapse = ", ")))
cat(sprintf("B_GRID    : %s\n", paste(B_GRID,  collapse = ", ")))
cat(sprintf("TISSUES   : %s\n", paste(tissues, collapse = ", ")))

old_wd <- setwd("code"); source("setup.R"); setwd(old_wd)

out_dir_fig <- "output/active_vs_passive/figures"
out_dir_res <- "output/active_vs_passive/results"
dir.create(out_dir_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_res, recursive = TRUE, showWarnings = FALSE)

analysis <- CircadianAnalysisOptions(alpha = 0.05, p.adjust.method = "BH",
                                      fdr_thresholds = 0.05)

# =====================================================================
# 1. Load all 3 baboon tissues and fit FMM per top-K rhythmic gene
# =====================================================================
cat("\n=== Step 1: load tissues + fit FMM ===\n")
load("data/CAMO_PRC_hmb.RData")

bio_list <- list()
mat_list <- list()
tod_list <- list()

for (t in tissues) {
  cat(sprintf("\n--- %s ---\n", t))
  prep <- prepCircadianData(baboon_withTOD$baboon[[t]],
                             times = baboon_withTOD$tod[[t]] %% 24,
                             input_type = "cpm")
  mat  <- prep$data[rowSums(prep$data > 0) >= 6, , drop = FALSE]
  tod  <- prep$times
  set.seed(GLOBAL_SEED)
  g_idx <- sample(nrow(mat), min(NGENES, nrow(mat)))
  mat_sub <- mat[g_idx, ]

  bio_rds <- file.path(out_dir_res, sprintf("bio_baboon_%s_fmm.rds", t))
  if (file.exists(bio_rds) && !identical(Sys.getenv("REFIT"), "true")) {
    cat(sprintf("Loading cached: %s\n", bio_rds))
    bio_t <- readRDS(bio_rds)
  } else {
    bio_t <- estCircadianParamFMM(mat_sub, tod,
                                   min_rhythm_pval = RHYTHM_PVAL,
                                   top_k    = TOP_K_FMM,
                                   mc.cores = N_CORES,
                                   verbose  = TRUE)
    bio_t$ngenes <- NGENES
    saveRDS(bio_t, bio_rds)
  }

  bio_list[[t]] <- bio_t
  mat_list[[t]] <- mat_sub
  tod_list[[t]] <- tod
}
rm(baboon_withTOD, gtex, mice)

# Print diagnostic summary
cat("\n=== Diagnostics summary ===\n")
for (t in tissues) {
  d <- bio_list[[t]]$diagnostics
  cat(sprintf("  %s: beta_hat=%.3f  sigma_alpha_hat=%.3fh  R2_med=%.3f  n_fitted=%d\n",
              t, d$beta_hat, d$sigma_alpha_hat, d$R2_median, d$n_fitted))
}

# =====================================================================
# 2. B-vs-m grid simulation per tissue, per truth condition
# =====================================================================
# For each (tissue, truth, N, B) cell, simulate active design at exactly
# B equispaced timepoints, m = N/B reps per timepoint. Skip cells where
# B does not divide N evenly.

run_truth <- function(bio_t, tod, truth) {
  # truth = "cosinor" -> force omega=1, alpha_dist=NULL
  # truth = "fmm"     -> use empirical paired (omega, alpha) from bio_t
  bio_use <- bio_t
  if (truth == "cosinor") {
    bio_use$omega_rhythmic <- NULL
    bio_use$omega_dist     <- list(family = "fixed", value = 1.0)
    bio_use$alpha_rhythmic <- NULL
    bio_use$alpha_dist     <- NULL
  }
  # else: keep paired empirical

  result_grid <- expand.grid(N = N_TOTAL, B = B_GRID, stringsAsFactors = FALSE)
  result_grid$power <- NA_real_

  for (i in seq_len(nrow(result_grid))) {
    N <- result_grid$N[i]; B <- result_grid$B[i]
    if (N %% B != 0L) next
    m <- N %/% B
    tps <- seq(0, 24 * (1 - 1/B), length.out = B)
    cts_design <- rep(tps, each = m)

    design_i <- CircadianDesignOptions(sample_sizes = N,
                                        nsims        = NSIMS,
                                        design       = "active",
                                        cts          = cts_design,
                                        B_values     = B)
    set.seed(GLOBAL_SEED + i)
    res <- runSimsSingleCohort(bio_use, design_i, analysis,
                                mc.cores = N_CORES, verbose = FALSE)
    result_grid$power[i] <- mean(res$marginal_power, na.rm = TRUE)
  }
  result_grid
}

cat("\n=== Step 2: simulate B-vs-m grid ===\n")
results <- list()
for (t in tissues) {
  cat(sprintf("\n--- Tissue: %s ---\n", t))
  results[[t]] <- list(
    cosinor = run_truth(bio_list[[t]], tod_list[[t]], "cosinor"),
    fmm     = run_truth(bio_list[[t]], tod_list[[t]], "fmm")
  )
  cat(sprintf("  cosinor truth power range: [%.3f, %.3f]\n",
              min(results[[t]]$cosinor$power, na.rm = TRUE),
              max(results[[t]]$cosinor$power, na.rm = TRUE)))
  cat(sprintf("  FMM     truth power range: [%.3f, %.3f]\n",
              min(results[[t]]$fmm$power, na.rm = TRUE),
              max(results[[t]]$fmm$power, na.rm = TRUE)))
}

# =====================================================================
# 3. Save results
# =====================================================================
out <- list(
  bio_list = bio_list,
  results  = results,
  N_TOTAL  = N_TOTAL,
  B_GRID   = B_GRID,
  tissues  = tissues
)
ts       <- format(Sys.time(), "%Y%m%d_%H%M%S")
rds_path <- file.path(out_dir_res, sprintf("fig5_active_vs_passive_%s.rds", ts))
saveRDS(out, rds_path)
cat(sprintf("\nResults saved -> %s\n", rds_path))

# =====================================================================
# 4. Plot — 2 rows × 3 cols
# =====================================================================
cat("\n=== Step 3: render Fig 5 ===\n")

fig_paths <- c(
  file.path(out_dir_fig, "fig5_active_vs_passive.pdf"),
  "output/main_figures/Fig5_active_vs_passive.pdf"
)
dir.create("output/main_figures", recursive = TRUE, showWarnings = FALSE)

# Distinct color per N value
n_pal <- c("#1f77b4", "#ff7f0e", "#2ca02c")[seq_along(N_TOTAL)]

for (out_pdf in fig_paths) {
  pdf(out_pdf, width = 14, height = 8)
  par(mfrow = c(2, length(tissues)), mar = c(4.5, 4.5, 3, 1.2), las = 1,
      cex.lab = 1.1, cex.axis = 1.0, cex.main = 1.05, font.main = 2)

  panel_letters <- LETTERS[seq_len(2 * length(tissues))]
  pi <- 1L

  for (truth in c("cosinor", "fmm")) {
    for (t in tissues) {
      df <- results[[t]][[truth]]
      d  <- bio_list[[t]]$diagnostics

      title_main <- sprintf("(%s) %s — %s truth%s",
        panel_letters[pi], t,
        if (truth == "cosinor") "cosinor" else "FMM",
        if (truth == "fmm") sprintf(" (β̂=%.2f)", d$beta_hat) else "")
      pi <- pi + 1L

      plot(NA, xlim = range(B_GRID), ylim = c(0, 1),
           xlab = "B (number of distinct timepoints)",
           ylab = "Power",
           main = title_main, log = "x")
      abline(h = 0.80, lty = 3, col = "grey50")
      for (j in seq_along(N_TOTAL)) {
        sub <- df[df$N == N_TOTAL[j], ]
        sub <- sub[!is.na(sub$power), ]
        if (nrow(sub) == 0) next
        sub <- sub[order(sub$B), ]
        lines(sub$B, sub$power, col = n_pal[j], lwd = 2.2)
        points(sub$B, sub$power, pch = 19, col = n_pal[j], cex = 1.0)
      }
      if (t == tissues[1] && truth == "cosinor") {
        legend("bottomright", bty = "n",
               legend = sprintf("N=%d", N_TOTAL),
               lwd = 2.2, col = n_pal, cex = 0.95)
      }
    }
  }

  mtext("Fig 5: B-vs-m trade-off under cosinor vs FMM truth (Baboon CAMO, active)",
        outer = TRUE, line = -1.3, cex = 1.05, font = 2)

  dev.off()
  cat(sprintf("Saved: %s\n", out_pdf))
}

cat("\n=== Done ===\n")
