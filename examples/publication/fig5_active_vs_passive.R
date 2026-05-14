#' =======================================================================
#' fig5_active_vs_passive.R — B vs m trade-off under empirical FMM truth,
#'                            comparing DCP and FMM (K=2 harmonic LRT) detectors
#' =======================================================================
#'
#' Two-row figure on three baboon CAMO active-design tissues. Both rows
#' use the SAME empirical FMM truth (paired (A, σ, ω, α) per gene from
#' estCircadianParamFMM); only the detector differs.
#'
#'   Row 1 — DCP detection
#'     Empirical confirmation that DCP power is B-INVARIANT even under
#'     non-cosinor (FMM) truth. Lines per panel (one per B value) should
#'     overlap because the cosinor F-test projects onto the first harmonic
#'     and its NCP = (A·c(ω))²·N/(2σ²) is independent of B.
#'
#'   Row 2 — FMM (K=2 harmonic LRT) detection
#'     Demonstrates the practical advantage of FMM: lines per panel
#'     (one per B value) SPREAD with B because higher B → better resolves
#'     the asymmetric peak shape and improves likelihood-ratio statistic.
#'
#' Layout: 2 rows × 3 columns
#'   Cols (3 baboon tissues):  LIV, LUN, KIC
#'   All have R² > 0.77 FMM fits and clear cosinor violation
#'   (omega-median 0.40-0.56).
#'
#' USAGE:
#'   Rscript examples/publication/fig5_active_vs_passive.R
#'   SMOKE_TEST=true Rscript examples/publication/fig5_active_vs_passive.R

SMOKE_TEST  <- identical(Sys.getenv("SMOKE_TEST"), "true")
NGENES      <- if (SMOKE_TEST) 500L  else 2000L
NSIMS       <- if (SMOKE_TEST) 5L    else 20L
N_CORES     <- as.integer(Sys.getenv("MC_CORES", unset = "20"))
TOP_K_FMM   <- if (SMOKE_TEST) 50L   else 300L  # codebase default min(300, n_cand)
RHYTHM_PVAL <- 0.01
GLOBAL_SEED <- 2025L

# Sample size grid (small enough to escape saturation, divisible by all B values)
N_GRID  <- if (SMOKE_TEST) c(12L, 24L) else c(24L, 48L, 72L, 96L, 120L, 144L)
B_GRID  <- if (SMOKE_TEST) c(3L, 12L)  else c(4L, 6L, 8L, 12L, 24L)

# Three most extreme baboon tissues from the May 2026 omega-median scan:
# KIC ω-med=0.40, KIM ω-med=0.35, SUN ω-med=0.31. All R² > 0.82 FMM fits.
# These are the tissues most likely to show FMM > DCP at the population level.
tissues <- c("KIC", "KIM", "SUN")

cat(sprintf("Mode      : %s\n", if (SMOKE_TEST) "SMOKE" else "PRODUCTION"))
cat(sprintf("NGENES    : %d\n", NGENES))
cat(sprintf("NSIMS     : %d\n", NSIMS))
cat(sprintf("MC_CORES  : %d\n", N_CORES))
cat(sprintf("TOP_K_FMM : %d\n", TOP_K_FMM))
cat(sprintf("N_GRID    : %s\n", paste(N_GRID, collapse = ", ")))
cat(sprintf("B_GRID    : %s\n", paste(B_GRID, collapse = ", ")))
cat(sprintf("TISSUES   : %s\n", paste(tissues, collapse = ", ")))

old_wd <- setwd("code"); source("setup.R"); setwd(old_wd)

out_dir_fig <- "output/active_vs_passive/figures"
out_dir_res <- "output/active_vs_passive/results"
dir.create(out_dir_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_res, recursive = TRUE, showWarnings = FALSE)

analysis <- CircadianAnalysisOptions(alpha           = 0.05,
                                      p.adjust.method = "BH",
                                      fdr_thresholds  = 0.05)

# ====================================================================
# 1. Load tissues + fit FMM (cached)
# ====================================================================
cat("\n=== Step 1: load tissues + fit FMM ===\n")
load("data/CAMO_PRC_hmb.RData")

bio_list <- list()
for (t in tissues) {
  cat(sprintf("\n--- %s ---\n", t))
  prep <- prepCircadianData(baboon_withTOD$baboon[[t]],
                             times = baboon_withTOD$tod[[t]] %% 24,
                             input_type = "cpm")
  mat  <- prep$data[rowSums(prep$data > 0) >= 6, , drop = FALSE]
  set.seed(GLOBAL_SEED)
  g    <- sample(nrow(mat), min(NGENES, nrow(mat)))
  mat_sub <- mat[g, ]

  bio_rds <- file.path(out_dir_res,
                        sprintf("bio_baboon_%s_NG%d_K%d.rds",
                                t, NGENES, TOP_K_FMM))
  if (file.exists(bio_rds) && !identical(Sys.getenv("REFIT"), "true")) {
    cat(sprintf("Loading cached: %s\n", bio_rds))
    bio_t <- readRDS(bio_rds)
  } else {
    bio_t <- estCircadianParamFMM(mat_sub, prep$times,
                                   min_rhythm_pval = RHYTHM_PVAL,
                                   top_k    = TOP_K_FMM,
                                   mc.cores = N_CORES,
                                   verbose  = TRUE)
    bio_t$ngenes <- NGENES
    saveRDS(bio_t, bio_rds)
  }
  bio_list[[t]] <- bio_t
}
rm(baboon_withTOD, gtex, mice)

# Diagnostic summary
cat("\n=== Diagnostics ===\n")
for (t in tissues) {
  d <- bio_list[[t]]$diagnostics
  cat(sprintf("  %s: β̂=%.3f  σ̂_α=%.3fh  R²=%.3f  n_fitted=%d\n",
              t, d$beta_hat, d$sigma_alpha_hat, d$R2_median, d$n_fitted))
}

# ====================================================================
# 2. B-vs-m simulation per (tissue, detector, B)
# ====================================================================
# Truth: paired empirical (A, σ, ω, α) — bio_list already has these.
# Detector varies via runSimsSingleCohort(method=...).

run_curve <- function(bio_t, B, method, K = 2L) {
  N_valid <- N_GRID[N_GRID %% B == 0L]
  pwr <- vapply(N_valid, function(N) {
    m <- N %/% B
    cts_design <- rep(seq(0, 24*(1 - 1/B), length.out = B), each = m)
    design_i <- CircadianDesignOptions(sample_sizes = N, nsims = NSIMS,
                                        design = "active",
                                        cts = cts_design, B_values = B)
    set.seed(GLOBAL_SEED + N)
    res <- runSimsSingleCohort(bio_t, design_i, analysis,
                                method = method, K = K,
                                mc.cores = N_CORES, verbose = FALSE)
    mean(res$marginal_power, na.rm = TRUE)
  }, numeric(1))
  data.frame(N = N_valid, power = pwr, B = B, method = method,
             stringsAsFactors = FALSE)
}

cat("\n=== Step 2: B-vs-m simulation ===\n")
results <- list()
for (t in tissues) {
  cat(sprintf("\n--- Tissue: %s ---\n", t))
  results[[t]] <- list()
  for (m in c("DCP", "FMM")) {
    cat(sprintf("  %s detector\n", m))
    rows <- list()
    for (B in B_GRID) {
      df <- run_curve(bio_list[[t]], B, m)
      rows[[as.character(B)]] <- df
      cat(sprintf("    B=%-2d : %s\n", B,
                  paste(sprintf("N=%d:%.2f", df$N, df$power),
                        collapse = "  ")))
    }
    results[[t]][[m]] <- do.call(rbind, rows)
  }
}

# ====================================================================
# 3. Save
# ====================================================================
out <- list(
  bio_list = bio_list,
  results  = results,
  N_GRID   = N_GRID,
  B_GRID   = B_GRID,
  tissues  = tissues
)
ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
rds_path <- file.path(out_dir_res, sprintf("fig5_active_vs_passive_%s.rds", ts))
saveRDS(out, rds_path)
cat(sprintf("\nSaved: %s\n", rds_path))

# ====================================================================
# 4. Plot — 2 rows × 3 cols
# ====================================================================
cat("\n=== Step 3: render Fig 5 ===\n")
fig_paths <- c(
  file.path(out_dir_fig, "fig5_active_vs_passive.pdf"),
  "output/main_figures/Fig5_active_vs_passive.pdf"
)
dir.create("output/main_figures", recursive = TRUE, showWarnings = FALSE)

# Color per B value (viridis-like)
.pal <- function(n) {
  cols <- c("#440154", "#3b528b", "#21918c", "#5ec962", "#fde725", "#bb3754", "#ff7e6d")
  cols[seq_len(min(n, length(cols)))]
}
pal_B <- .pal(length(B_GRID))

for (out_pdf in fig_paths) {
  cairo_pdf(out_pdf, width = 14, height = 8)
  par(mfrow = c(2, length(tissues)), mar = c(4.4, 4.4, 3, 1.2), las = 1,
      cex.lab = 1.05, cex.axis = 0.95, cex.main = 1.05, font.main = 2)

  panel_idx <- 1L
  for (m in c("DCP", "FMM")) {
    for (t in tissues) {
      df <- results[[t]][[m]]
      d  <- bio_list[[t]]$diagnostics

      plot(NA, xlim = range(N_GRID), ylim = c(0, 1),
           xlab = "N (total samples)",
           ylab = sprintf("Power (%s)", m),
           main = sprintf("%s — %s   (β̂=%.2f, R²=%.2f)",
                           t, m, d$beta_hat, d$R2_median))
      abline(h = 0.80, lty = 3, col = "grey50")
      for (k in seq_along(B_GRID)) {
        sub <- df[df$B == B_GRID[k] & !is.na(df$power), ]
        sub <- sub[order(sub$N), ]
        if (nrow(sub) == 0) next
        lines(sub$N, sub$power, col = pal_B[k], lwd = 2.0)
        points(sub$N, sub$power, pch = 19, col = pal_B[k], cex = 0.9)
      }
      if (panel_idx == 1) {
        legend("bottomright", bty = "n", title = "B",
               legend = sprintf("%d", B_GRID),
               lwd = 2, col = pal_B, cex = 0.9)
      }
      panel_idx <- panel_idx + 1L
    }
  }

  mtext("Fig 5: B-vs-m trade-off under empirical FMM truth — DCP (top) vs FMM K=2 (bottom)",
        outer = TRUE, line = -1.4, cex = 1.0, font = 2)

  dev.off()
  cat(sprintf("Saved: %s\n", out_pdf))
}

cat("\n=== Done ===\n")
