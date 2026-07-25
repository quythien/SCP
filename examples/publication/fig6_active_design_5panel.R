#' Fig 6 generator: B-invariance of power under a balanced active design.
#' Five panels: (A,B) single-cohort biomarker detection, K=1 and K=2, on the
#' GSE160521 Putamen control pilot; (C,D,E) two-group differential power for
#' DR, DP and DM on the Putamen control vs Caudate control comparison.
#'
#' The figure combines two calibrations:
#'   * A/B biomarker panels are read from a pre-computed cache,
#'     output/two_harmonic/results/fig6_v9_data.rds, produced by the shipped
#'     two-harmonic active B-sweep on the Putamen control pilot (see FIGURES.md;
#'     that generator is git-recoverable from commit a1ac962).
#'   * C/D/E differential panels are produced here in two stages: a two-group
#'     pilot calibration (Stage 1) and a DR/DP/DM active B-sweep (Stage 2).
#' Stage 3 renders all five panels into submission/figures.
#'
#' Reproducibility notes:
#'   * Requires the controlled-access GSE160521 striatal CPM matrices and the
#'     matched-donor clinical file (CAMO server path below); runs on this server
#'     only. The shipped per-region rhythm_fit summaries are de-identified
#'     aggregates and are NOT redistributed.
#'   * Activate the fast C++ cosinor path before running (see FIGURES.md).
#'   * Stages 1-2 are compute-heavy and cache to RDS; on a re-run the script
#'     skips them if the caches already exist and goes straight to the render.
#'   * All randomness is seeded (calibration seed 2025; per-B sweep seed
#'     GLOBAL_SEED + B), so the figure is reproducible on this server.

setwd(Sys.getenv("SCP_ROOT", unset = "."))  # run from the package root
suppressWarnings(suppressMessages({ library(Rcpp); library(SCP) }))
dyn.load(system.file("libs", "SCP.so", package = "SCP")); .CPP_LOADED <- TRUE

GLOBAL_SEED <- 2025L
PERIOD      <- 24
FDR         <- 0.05
N_CORES     <- as.integer(Sys.getenv("MC_CORES", unset = "48"))
B_GRID      <- c(4L, 6L, 8L, 12L, 24L)       # identifiable differential designs
N_GRID      <- c(24L, 48L, 96L, 144L, 192L, 240L)
NSIMS       <- 60L
NGENES      <- 3000L

# GSE160521 striatal CPM + the matched-donor clinical file are controlled-access
# (CAMO server); not redistributed. Point SCP_GSE160521_DIR at your local copy to
# recalibrate, or run from the cached RDS below.
KD          <- Sys.getenv("SCP_GSE160521_DIR", unset = "data/gse160521")
AB_CACHE    <- "output/two_harmonic/results/fig6_v9_data.rds"          # A/B biomarker (Putamen ctrl)
PILOT_RDS   <- "output/diagnostics/pilot_put_vs_cau_control.rds"       # Stage 1 output
SWEEP_RDS   <- "output/diagnostics/Fig6_diff_Bsweep_PutCau.rds"        # Stage 2 output
OUT_PDF     <- "submission/figures/Fig6_active_design_5panel.pdf"      # Stage 3 output

# ============================================================================
# Stage 1 -- Two-group pilot calibration: Putamen control (g1) vs Caudate control (g2)
# ============================================================================
if (!file.exists(PILOT_RDS)) {
  clin <- read.csv(file.path(KD, "DS_clinical_1221_rm97_rm231_matchIndex34.csv"), row.names = 1)
  ctl  <- clin[clin$Diagnostic.Category == "CONTROL", ]

  load_region <- function(fpat) {
    f    <- list.files(KD, pattern = fpat, full.names = TRUE)
    ex   <- as.matrix(read.csv(f[1], row.names = 1, check.names = FALSE))
    cols <- intersect(ctl$pair, colnames(ex))
    tod  <- ctl$CorrectedTOD[match(cols, ctl$pair)] %% 24
    list(ex = ex[, cols, drop = FALSE], tod = tod)
  }
  put <- load_region("^Putamen_CPMfiltered_logCPM.*csv$")
  cau <- load_region("^Caudate_CPMfiltered_logCPM.*csv$")

  g <- intersect(rownames(put$ex), rownames(cau$ex))               # align genes
  P <- put$ex[g, , drop = FALSE]; C <- cau$ex[g, , drop = FALSE]

  set.seed(GLOBAL_SEED)
  bio <- estCircadianParamTwoGroup(
    data_1 = P, data_2 = C,
    times_1 = put$tod, times_2 = cau$tod,
    period = PERIOD, paired_sigma = TRUE, verbose = TRUE)
  saveRDS(bio, PILOT_RDS)
  cat(sprintf("Stage 1: prop_DR=%.4f prop_DP=%.4f prop_DM=%.4f -> %s\n",
              bio$prop_DR, bio$prop_DP, bio$prop_DM, PILOT_RDS))
}

# ============================================================================
# Stage 2 -- Differential active B-sweep (DR, DP, DM) over the (N, B) grid
# ============================================================================
if (!file.exists(SWEEP_RDS)) {
  bio <- readRDS(PILOT_RDS); bio$ngenes <- NGENES
  analysis  <- CircadianAnalysisOptions(alpha = FDR, p.adjust.method = "BH", fdr_thresholds = FDR)
  ENDPOINTS <- c("DR", "DP", "DM")
  prepFn    <- get(".prepDiffStratified", envir = asNamespace("SCP"))
  r_breaks  <- c(0, Inf)

  pw <- lapply(ENDPOINTS, function(e) matrix(NA_real_, length(N_GRID), length(B_GRID),
                                             dimnames = list(N_GRID, B_GRID)))
  se <- lapply(ENDPOINTS, function(e) matrix(NA_real_, length(N_GRID), length(B_GRID),
                                             dimnames = list(N_GRID, B_GRID)))
  names(pw) <- names(se) <- ENDPOINTS

  for (k in seq_along(B_GRID)) {
    B   <- B_GRID[k]
    cts <- seq(0, PERIOD * (1 - 1/B), length.out = B)          # balanced, equally spaced
    design <- CircadianDesignOptions(sample_sizes = N_GRID, nsims = NSIMS, design = "active",
                                     cts = cts, B_values = B, test_types = ENDPOINTS)
    set.seed(GLOBAL_SEED + B)
    res <- runDifferentialPower(bio, design, analysis, methods = "DCP",
                                test_types = ENDPOINTS, plot = FALSE, verbose = FALSE,
                                mc.cores = N_CORES)
    for (e in ENDPOINTS) {
      prep <- prepFn(res, e, fdr_thresholds = FDR, r_breaks = r_breaks)
      ms   <- prep$marginal_sim[, 1, , drop = TRUE]
      if (is.null(dim(ms))) ms <- matrix(ms, nrow = length(N_GRID))
      pw[[e]][, k] <- rowMeans(ms, na.rm = TRUE)
      se[[e]][, k] <- apply(ms, 1, function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
    }
    cat(sprintf("Stage 2: B=%d done\n", B))
  }
  saveRDS(list(pilot = "GSE160521 Putamen vs Caudate control (paired)",
               B_GRID = B_GRID, N_GRID = N_GRID, NSIMS = NSIMS, NGENES = NGENES, FDR = FDR,
               endpoints = ENDPOINTS, power = pw, se = se,
               prop_DR = bio$prop_DR, prop_DP = bio$prop_DP, prop_DM = bio$prop_DM),
          SWEEP_RDS)
  cat(sprintf("Stage 2: saved %s\n", SWEEP_RDS))
}

# ============================================================================
# Stage 3 -- Render the 5-panel figure
# ============================================================================
ab   <- readRDS(AB_CACHE)     # A/B biomarker panels (Putamen control)
diff <- readRDS(SWEEP_RDS)    # C/D/E differential panels (Putamen vs Caudate)

# B -> colour map, consistent across every panel (viridis; last = magenta as shipped)
B_ALL   <- c(3L, 4L, 6L, 8L, 12L, 24L)
COL_ALL <- c("#440154", "#3b528b", "#21918c", "#5ec962", "#fde725", "#7d028c")
names(COL_ALL) <- B_ALL

draw_ab <- function(N, M, SE, B_grid, letter, main_text, jitter_step = 1.2,
                    force_zero_cols = integer(0)) {
  cols <- COL_ALL[as.character(B_grid)]
  matplot(N, 100 * M, type = "n", ylim = c(0, 100), xlim = c(0, max(N) * 1.05),
          xlab = "Sample size (N)", ylab = "Power (%)", main = "")
  title(main = sprintf("%s   %s", letter, main_text),
        adj = 0.5, font.main = 2, cex.main = 1.10, line = 0.3)
  abline(h = 80, lty = 2, col = "grey50", lwd = 1.2); grid()
  for (k in seq_along(B_grid)) {
    xj <- N + (k - (length(B_grid) + 1) / 2) * jitter_step
    lo <- 100 * pmax(0, M[, k] - SE[, k]); hi <- 100 * pmin(1, M[, k] + SE[, k])
    ok <- is.finite(lo) & is.finite(hi) & SE[, k] > 0
    if (any(ok))
      arrows(xj[ok], lo[ok], xj[ok], hi[ok], code = 3, angle = 90,
             length = 0.025, col = cols[k], lwd = 1.0)
    M_k <- M[, k]
    if (k %in% force_zero_cols) M_k[!is.finite(M_k)] <- 0 else M_k[!is.finite(M_k)] <- NA
    lines(xj, 100 * M_k, type = "o", pch = 19, lwd = 1.6, col = cols[k], cex = 0.55)
  }
  legend("bottomright", paste0("B=", B_grid), col = cols, lty = 1, pch = 19,
         lwd = 1.4, cex = 0.62, bty = "o", box.col = "grey70", box.lwd = 0.5,
         inset = 0.01, y.intersp = 0.82, bg = "white", seg.len = 1.4)
}

draw_cde <- function(N, M, SE, B_grid, letter, main_text, jitter_step = 2.2,
                     show_legend = FALSE) {
  cols <- COL_ALL[as.character(B_grid)]
  plot(NA, xlim = c(0, max(N) * 1.03), ylim = c(0, 100),
       xlab = "Sample size (N)", ylab = "Power (%)", main = "")
  title(main = sprintf("%s   %s", letter, main_text),
        adj = 0.5, font.main = 2, cex.main = 1.10, line = 0.3)
  abline(h = 80, lty = 2, col = "grey50", lwd = 1.2); grid()
  for (k in seq_along(B_grid)) {
    xj <- N + (k - (length(B_grid) + 1) / 2) * jitter_step
    y  <- 100 * M[, k]; s <- 100 * SE[, k]
    lo <- pmax(0, y - s); hi <- pmin(100, y + s)
    ok <- is.finite(lo) & is.finite(hi) & s > 0
    if (any(ok))
      arrows(xj[ok], lo[ok], xj[ok], hi[ok], code = 3, angle = 90,
             length = 0.025, col = cols[k], lwd = 1.0)
    lines(xj, y, type = "o", pch = 19, lwd = 1.6, col = cols[k], cex = 0.55)
  }
  if (show_legend)
    legend("bottomright", paste0("B=", B_grid), col = cols, lty = 1, pch = 19,
           lwd = 1.4, cex = 0.62, bty = "o", box.col = "grey70", box.lwd = 0.5,
           inset = 0.01, y.intersp = 0.82, bg = "white", seg.len = 1.4)
}

dir.create(dirname(OUT_PDF), recursive = TRUE, showWarnings = FALSE)
pdf(OUT_PDF, width = 11.0, height = 7.6)
layout(matrix(c(1, 1, 1, 2, 2, 2,
                3, 3, 4, 4, 5, 5), nrow = 2, byrow = TRUE))
par(mai = c(0.82, 0.80, 0.42, 0.15), mgp = c(2.3, 0.55, 0), oma = c(0, 0, 3.3, 0),
    cex.axis = 0.95, cex.lab = 1.10, font.main = 2)

b34 <- which(ab$B %in% c(3L, 4L))   # non-identifiable K=2 designs, drawn flat at 0
draw_ab(ab$N, ab$pwr_K1, ab$se_K1, ab$B, "A", "Biomarker detection, K=1 (cosinor)")
draw_ab(ab$N, ab$pwr_K2, ab$se_K2, ab$B, "B", "Biomarker detection, K=2 (two-harmonic)",
        force_zero_cols = b34)
draw_cde(diff$N_GRID, diff$power$DR, diff$se$DR, diff$B_GRID, "C",
         "Differential rhythmicity (DR)", show_legend = TRUE)
draw_cde(diff$N_GRID, diff$power$DP, diff$se$DP, diff$B_GRID, "D",
         "Differential phase (DP)")
draw_cde(diff$N_GRID, diff$power$DM, diff$se$DM, diff$B_GRID, "E",
         "Differential mesor (DM)")
mtext(bquote("Active design B vs m trade-off (Putamen control, n = 59, " *
             tilde(r) * " = 1.06)"),
      outer = TRUE, side = 3, line = 1.2, font = 2, cex = 1.18)
mtext("Biomarker detection (A, B) and differential endpoints (C, D, E)",
      outer = TRUE, side = 3, line = 0.1, font = 1, cex = 0.95)
dev.off()
cat(sprintf("Stage 3: saved %s\n", OUT_PDF))
