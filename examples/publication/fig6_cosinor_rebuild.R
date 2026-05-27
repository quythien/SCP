#' =======================================================================
#' fig6_cosinor_rebuild.R — Cosinor-only Fig 6 (active KIM design)
#' =======================================================================
#'
#' Replaces the FMM-based Fig 6. NO FMM, NO omega_rhythmic. Two frameworks,
#' each simulated with its OWN matched generator and detected with the
#' matching test (the "correct generator for each" principle):
#'
#'   Panel A — Cosinor framework
#'     Truth : estCircadianParam(KIM)            (1-harmonic, paired A,sigma)
#'     Detect: method = "DCP" (cosinor F-test)
#'     B-vs-m power vs N at B in {6,8,12,24}.
#'
#'   Panel B — Phase estimation (cosinor framework only)
#'     Truth : cosinor; estimate phi_hat^cos via fitCosinorAll_fast.
#'     Median circular phase MSE (h^2) vs N at B in {4,6,8,12,24}.
#'     With pure cosinor truth there is NO first-harmonic-projection bias
#'     floor, so MSE declines ~1/(N r^2) instead of saturating early.
#'
#'   Panel C — Two-harmonic framework
#'     Truth : estCircadianParam2H(KIM)          (paired A1,phi1,A2,phi2,sigma)
#'     Detect: method = "FMM", K = 2 (K-harmonic F-test)
#'     B-vs-m power vs N at B in {6,8,12,24}.
#'
#' OUTPUT:
#'   output/phase_mse/results/fig6_cosinor_rebuild_<ts>.rds
#'   output/main_figures/Fig6_active_design.pdf
#'   submission/figures/Fig6_active_design.pdf
#'
#' USAGE:
#'   SMOKE_TEST=true Rscript examples/publication/fig6_cosinor_rebuild.R
#'   Rscript examples/publication/fig6_cosinor_rebuild.R

setwd("/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")

SMOKE_TEST  <- identical(Sys.getenv("SMOKE_TEST"), "true")
NGENES      <- if (SMOKE_TEST) 500L  else 2000L
NSIMS       <- if (SMOKE_TEST) 5L    else 100L
N_CORES     <- as.integer(Sys.getenv("MC_CORES", unset = "20"))
TOP_K       <- if (SMOKE_TEST) 50L   else 300L
RHYTHM_PVAL <- 0.01
GLOBAL_SEED <- 2025L
PERIOD      <- 24
R_FLOOR     <- 0.0    # phase MSE over ALL rhythmic genes (GTEx Liver median r=0.35; r>=1 subset empty)

# Sweeps. Power panels start at B=6 (K=2 Nyquist B>=2K+1=5); phase panel
# adds B=4 (cosinor only needs B>=3).
B_GRID_PWR  <- if (SMOKE_TEST) c(6L, 12L)        else c(6L, 8L, 12L, 24L)
B_GRID_MSE  <- if (SMOKE_TEST) c(4L, 12L)        else c(4L, 6L, 8L, 12L, 24L)
# Low-SNR passive pilot (median r=0.35): extend N upward so power spans
# its full rise (~0 to ~0.85) and the phase-MSE shrink is visible.
N_GRID      <- if (SMOKE_TEST) c(48L, 96L)       else c(24L, 48L, 96L, 144L, 192L, 288L)

cat(sprintf("Mode   : %s\n", if (SMOKE_TEST) "SMOKE" else "PRODUCTION"))
cat(sprintf("NGENES : %d | NSIMS : %d | TOP_K : %d | CORES : %d\n",
            NGENES, NSIMS, TOP_K, N_CORES))
cat(sprintf("N_GRID : %s\n", paste(N_GRID, collapse = ", ")))

source("examples/publication/_pub_style.R")
old_wd <- setwd("code"); source("setup.R"); setwd(old_wd)

out_dir_res <- "output/phase_mse/results"
dir.create(out_dir_res, recursive = TRUE, showWarnings = FALSE)

analysis <- CircadianAnalysisOptions(alpha = 0.05, p.adjust.method = "BH",
                                     fdr_thresholds = 0.05)

# ====================================================================
# 1. Load GTEx Liver two-harmonic pilot; derive matched cosinor truth
# ====================================================================
# Passive pilot (n=262), well-powered 2nd harmonic. We estimate the
# (A1, phi1, A2, phi2, sigma) truth here, then SIMULATE under an active
# design below. GTEx Liver has median A2/A1=0.56 (real 2nd harmonic ->
# K=2 advantage) and median SNR r=0.35 (low -> phase MSE shrinks with N).
# The cosinor framework's truth is the first-harmonic-only restriction of
# the SAME pilot, so Panel C shows exactly what adding K=2 buys.
DATASET  <- "GTEx Liver"
cat(sprintf("\n=== Step 1: load 2H pilot (%s) + derive cosinor truth ===\n", DATASET))
PILOT_RDS <- "output/two_harmonic/results/pilot_2h_GTExLiver_topK300.rds"
bio_2h <- readRDS(PILOT_RDS)
bio_2h$ngenes <- NGENES
stopifnot(isTRUE(bio_2h$paired_2h), !is.null(bio_2h$amplitude2))

# Matched cosinor truth: drop the 2nd-harmonic block from the same pilot.
bio_cos <- bio_2h
bio_cos$amplitude2 <- NULL
bio_cos$phase2     <- NULL
bio_cos$paired_2h  <- FALSE

stopifnot(is.null(bio_cos$omega_rhythmic), is.null(bio_cos$amplitude2))
r_cos  <- bio_2h$amplitude / bio_2h$sigma_rhythmic
ratio  <- bio_2h$amplitude2 / bio_2h$amplitude
cat(sprintf("%s 2H pilot: n=%d  median r=A1/sigma=%.2f  median A2/A1=%.2f  frac A2/A1>0.5=%.2f\n",
            DATASET, length(bio_2h$amplitude), median(r_cos, na.rm = TRUE),
            median(ratio, na.rm = TRUE), mean(ratio > 0.5, na.rm = TRUE)))

# ====================================================================
# 2. B-vs-m power (matched generator + detector per framework)
# ====================================================================
run_power_curve <- function(bio, B, method, K = 2L) {
  N_valid <- N_GRID[N_GRID %% B == 0L]
  pwr <- vapply(N_valid, function(N) {
    m   <- N %/% B
    cts <- rep(seq(0, PERIOD * (1 - 1/B), length.out = B), each = m)
    des <- CircadianDesignOptions(sample_sizes = N, nsims = NSIMS,
                                  design = "active", cts = cts, B_values = B)
    set.seed(GLOBAL_SEED + N)
    res <- runSimsSingleCohort(bio, des, analysis, method = method, K = K,
                               mc.cores = N_CORES, verbose = FALSE)
    mean(res$marginal_power, na.rm = TRUE)
  }, numeric(1))
  data.frame(N = N_valid, power = pwr, B = B)
}

cat("\n=== Step 2: B-vs-m power ===\n")
res_dcp <- list(); res_k2 <- list()
for (B in B_GRID_PWR) {
  cat(sprintf("  cosinor(DCP)  B=%d\n", B))
  res_dcp[[as.character(B)]] <- run_power_curve(bio_cos, B, "DCP")
  cat(sprintf("  twoharm(K=2)  B=%d\n", B))
  res_k2[[as.character(B)]]  <- run_power_curve(bio_2h,  B, "FMM", K = 2L)
}

# ====================================================================
# 3. Phase MSE (cosinor framework only, cosinor truth)
# ====================================================================
circular_diff <- function(a, b, p = 24) { d <- abs(a - b) %% p; pmin(d, p - d) }

run_mse_cell <- function(B, N) {
  if (N %% B != 0L) return(NULL)
  m   <- N %/% B
  cts <- rep(seq(0, PERIOD * (1 - 1/B), length.out = B), each = m)
  med <- numeric(NSIMS)
  for (s in seq_len(NSIMS)) {
    set.seed(GLOBAL_SEED + 1000L * B + 7L * N + s)
    sim <- simCircadianSingleCohort(bio_cos, cts)   # cosinor truth (no omega)
    is_rh <- as.logical(sim$is_rhythmic)
    # cosinor truth: acrophase in hours is phase_g[rhythmic]; sim returns it
    phi_truth <- sim$phase_g[is_rh] %% PERIOD
    r_truth   <- sim$r_values[is_rh]
    hi <- !is.na(r_truth) & r_truth >= R_FLOOR
    fit <- fitCosinorAll_fast(sim$expr, cts, period = PERIOD)
    d   <- circular_diff(fit$phi[is_rh], phi_truth, PERIOD)
    med[s] <- median(d[hi]^2, na.rm = TRUE)
  }
  list(B = B, N = N, med_raw = med,
       med_mean = mean(med, na.rm = TRUE),
       med_lo = mean(med, na.rm = TRUE) - 1.96 * sd(med, na.rm = TRUE) / sqrt(sum(is.finite(med))),
       med_hi = mean(med, na.rm = TRUE) + 1.96 * sd(med, na.rm = TRUE) / sqrt(sum(is.finite(med))))
}

cat("\n=== Step 3: phase MSE (cosinor) ===\n")
mse <- list()
for (N in N_GRID) for (B in B_GRID_MSE) {
  if (N %% B != 0L) next
  cell <- run_mse_cell(B, N)
  if (!is.null(cell)) {
    mse[[sprintf("B%d_N%d", B, N)]] <- cell
    cat(sprintf("  B=%2d N=%3d  medMSE=%.3f h^2\n", B, N, cell$med_mean))
  }
}

# ====================================================================
# 4. Save
# ====================================================================
out <- list(GLOBAL_SEED = GLOBAL_SEED, NGENES = NGENES, NSIMS = NSIMS,
            TOP_K = TOP_K, PERIOD = PERIOD, pilot = "GTEx Liver 2H (passive pilot, active-design sim)",
            B_GRID_PWR = B_GRID_PWR, B_GRID_MSE = B_GRID_MSE, N_GRID = N_GRID,
            res_dcp = res_dcp, res_k2 = res_k2, mse = mse,
            median_r_cosinor = median(r_cos, na.rm = TRUE))
ts  <- format(Sys.time(), "%Y%m%d_%H%M%S")
rds <- file.path(out_dir_res, sprintf("fig6_cosinor_rebuild_%s.rds", ts))
saveRDS(out, rds)
cat(sprintf("\nSaved cache: %s\n", rds))

# ====================================================================
# 5. Render 3-panel PDF
# ====================================================================
draw_fig6 <- function(out_pdf) {
  cairo_pdf(out_pdf, width = 11.5, height = 4.5)
  par(mfrow = c(1, 3), mai = c(1.30, 1.15, 0.85, 0.18), mgp = c(3.6, 0.75, 0),
      oma = c(0, 0, 0.4, 0), family = "Helvetica", las = 1, tcl = -0.3,
      cex.axis = 1.35, cex.lab = 1.55, cex.main = 1.55, font.main = 2)
  pal_pwr <- pub_palette_sequential(length(B_GRID_PWR))
  pal_mse <- pub_palette_sequential(length(B_GRID_MSE))

  # --- Panel A: cosinor DCP power ---
  plot(NA, xlim = range(N_GRID), ylim = c(0, 1),
       xlab = expression(N ~ "(total samples)"),
       ylab = "Power (DCP, FDR = 0.05)",
       main = expression(bold("GTEx Liver,  cosinor F-test")))
  panel_label("A"); abline_80pct()
  for (k in seq_along(B_GRID_PWR)) {
    df <- res_dcp[[as.character(B_GRID_PWR[k])]]
    lines(df$N, df$power, col = pal_pwr[k], lwd = 1.8)
    points(df$N, df$power, col = pal_pwr[k], pch = 19, cex = 0.55)
  }
  pub_legend("bottomright",
             legend = sapply(B_GRID_PWR, function(b) as.expression(bquote(B == .(b)))),
             col = pal_pwr, lwd = 1.8, title = NULL)

  # --- Panel B: cosinor phase MSE ---
  med_m <- matrix(NA_real_, length(B_GRID_MSE), length(N_GRID),
                  dimnames = list(B_GRID_MSE, N_GRID))
  med_l <- med_m; med_h <- med_m
  for (i in seq_along(B_GRID_MSE)) for (j in seq_along(N_GRID)) {
    cell <- mse[[sprintf("B%d_N%d", B_GRID_MSE[i], N_GRID[j])]]
    if (!is.null(cell)) { med_m[i,j] <- cell$med_mean; med_l[i,j] <- cell$med_lo; med_h[i,j] <- cell$med_hi }
  }
  y_max <- max(med_h, med_m, na.rm = TRUE) * 1.05
  plot(NA, xlim = range(N_GRID), ylim = c(0, y_max),
       xlab = expression(N ~ "(total samples)"),
       ylab = expression("Median phase MSE (h" ^ 2 * ")"),
       main = expression(bold("Phase estimation  (" * hat(phi)[g]^{cos} * ")")))
  panel_label("B")
  for (i in seq_along(B_GRID_MSE)) {
    yi <- med_m[i, ]; if (all(is.na(yi))) next
    lines(N_GRID, yi, col = pal_mse[i], lwd = 1.8)
    points(N_GRID, yi, col = pal_mse[i], pch = 19, cex = 0.7)
    keep <- is.finite(yi) & is.finite(med_l[i,]) & is.finite(med_h[i,])
    if (any(keep)) arrows(N_GRID[keep], med_l[i,keep], N_GRID[keep], med_h[i,keep],
                          angle = 90, code = 3, length = 0.03, col = pal_mse[i], lwd = 1.0)
  }
  legend("topright",
         legend = sapply(B_GRID_MSE, function(b) as.expression(bquote(B == .(b)))),
         col = pal_mse, lwd = 1.8, pch = 19, bty = "o", bg = "white", cex = 0.85)

  # --- Panel C: two-harmonic K=2 power ---
  plot(NA, xlim = range(N_GRID), ylim = c(0, 1),
       xlab = expression(N ~ "(total samples)"),
       ylab = "Power (K = 2, FDR = 0.05)",
       main = expression(bold("GTEx Liver,  K = 2 harmonic F-test")))
  panel_label("C"); abline_80pct()
  for (k in seq_along(B_GRID_PWR)) {
    df <- res_k2[[as.character(B_GRID_PWR[k])]]
    lines(df$N, df$power, col = pal_pwr[k], lwd = 1.8)
    points(df$N, df$power, col = pal_pwr[k], pch = 19, cex = 0.55)
  }
  pub_legend("bottomright",
             legend = sapply(B_GRID_PWR, function(b) as.expression(bquote(B == .(b)))),
             col = pal_pwr, lwd = 1.8, title = NULL)
  dev.off()
  cat(sprintf("Saved: %s\n", out_pdf))
}

out_pdfs <- c("output/main_figures/Fig6_active_design.pdf",
              "submission/figures/Fig6_active_design.pdf")
for (q in out_pdfs) dir.create(dirname(q), recursive = TRUE, showWarnings = FALSE)
for (q in out_pdfs) draw_fig6(q)
