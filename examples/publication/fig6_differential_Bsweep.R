#' =======================================================================
#' fig6_differential_Bsweep.R — Differential B-invariance (advisor addition)
#' =======================================================================
#'
#' Question (advisor): B-invariance of biomarker detection is Zong theory,
#' but does the SAME B-invariance hold for the DIFFERENTIAL endpoints
#' (DR = differential rhythmicity, DP = differential phase,
#'  DM = differential mesor) under a balanced active design?
#'
#' Design: balanced, equally spaced ACTIVE design with B distinct TODs and
#' m = N/B replicates per TOD, sweeping B in {4,6,8,12,24} at matched N.
#' Two-group pilot = GTEx Adrenal Gland vs Liver (same as paper Fig 2),
#' cached at data/gtex_adr_vs_liv_pilot_paired.rds.
#'
#' Output:
#'   output/diagnostics/Fig6_differential_Bsweep.pdf
#'   output/diagnostics/Fig6_differential_Bsweep.rds
#'
#' USAGE:
#'   SMOKE_TEST=true Rscript examples/publication/fig6_differential_Bsweep.R
#'   MC_CORES=48 Rscript examples/publication/fig6_differential_Bsweep.R

setwd("/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")

SMOKE_TEST <- identical(Sys.getenv("SMOKE_TEST"), "true")

## ---- fast C++ path (195x on cosinor fits) --------------------------------
suppressWarnings(suppressMessages({library(Rcpp); library(SCP)}))
dyn.load(system.file("libs", "SCP.so", package = "SCP")); .CPP_LOADED <- TRUE

GLOBAL_SEED <- 2025L
PERIOD      <- 24
FDR         <- 0.05
N_CORES     <- as.integer(Sys.getenv("MC_CORES", unset = "48"))

# N grid divisible by every B in {4,6,8,12,24}  (LCM = 24)
B_GRID <- if (SMOKE_TEST) c(4L, 24L)             else c(4L, 6L, 8L, 12L, 24L)
N_GRID <- if (SMOKE_TEST) c(24L, 48L)            else c(24L, 48L, 96L, 144L, 192L, 240L)
NSIMS  <- if (SMOKE_TEST) 5L                     else 80L
NGENES <- if (SMOKE_TEST) 800L                   else 3000L

cat(sprintf("Mode   : %s\n", if (SMOKE_TEST) "SMOKE" else "PRODUCTION"))
cat(sprintf("B grid : %s\n", paste(B_GRID, collapse = ", ")))
cat(sprintf("N grid : %s\n", paste(N_GRID, collapse = ", ")))
cat(sprintf("NSIMS  : %d | NGENES : %d | CORES : %d\n", NSIMS, NGENES, N_CORES))

out_dir <- "output/diagnostics"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ====================================================================
# 1. Two-group pilot (GTEx Adrenal vs Liver), same as Fig 2
# ====================================================================
bio <- readRDS("data/gtex_adr_vs_liv_pilot_paired.rds")
bio$ngenes <- NGENES
cat(sprintf("\nPilot: ADR vs LIV  n1=%d n2=%d  prop_DR=%.3f prop_DP=%.3f prop_DM=%.3f\n",
            length(bio$cts), length(bio$cts2),
            bio$prop_DR, bio$prop_DP, bio$prop_DM))

analysis <- CircadianAnalysisOptions(alpha = FDR, p.adjust.method = "BH",
                                     fdr_thresholds = FDR)

ENDPOINTS <- c("DR", "DP", "DM")

# ====================================================================
# 2. Sweep B; extract marginal power per endpoint per N
# ====================================================================
# marginal power = (# true-endpoint genes discovered at FDR<=0.05) /
#                  (# true-endpoint genes), averaged over sims.
# Computed by the package's own .prepDiffStratified() (R/plot_diff.R),
# field marginal_sim[j, t, s]; r_breaks irrelevant to the marginal.
prepFn <- get(".prepDiffStratified", envir = asNamespace("SCP"))
r_breaks <- c(0, Inf)

# storage: results[[endpoint]] = matrix [N x B] of mean power (+ SE)
pw  <- lapply(ENDPOINTS, function(e) matrix(NA_real_, length(N_GRID), length(B_GRID),
                                            dimnames = list(N_GRID, B_GRID)))
se  <- lapply(ENDPOINTS, function(e) matrix(NA_real_, length(N_GRID), length(B_GRID),
                                            dimnames = list(N_GRID, B_GRID)))
names(pw) <- names(se) <- ENDPOINTS

for (k in seq_along(B_GRID)) {
  B   <- B_GRID[k]
  cts <- seq(0, PERIOD * (1 - 1/B), length.out = B)  # B distinct balanced TODs
  cat(sprintf("\n=== B = %d  (TODs: %s) ===\n", B,
              paste(round(cts, 1), collapse = ",")))

  design <- CircadianDesignOptions(
    sample_sizes = N_GRID, nsims = NSIMS, design = "active",
    cts = cts, B_values = B, test_types = ENDPOINTS)

  set.seed(GLOBAL_SEED + B)
  res <- runDifferentialPower(bio, design, analysis,
                              methods    = "DCP",
                              test_types = ENDPOINTS,
                              plot       = FALSE,
                              verbose    = FALSE,
                              mc.cores   = N_CORES)

  for (e in ENDPOINTS) {
    prep <- prepFn(res, e, fdr_thresholds = FDR, r_breaks = r_breaks)
    # marginal_sim: [n_sizes, n_thresh=1, nsims]
    ms <- prep$marginal_sim[, 1, , drop = TRUE]   # [n_sizes x nsims]
    if (is.null(dim(ms))) ms <- matrix(ms, nrow = length(N_GRID))
    pw[[e]][, k] <- rowMeans(ms, na.rm = TRUE)
    se[[e]][, k] <- apply(ms, 1, function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
    cat(sprintf("  %s power: %s\n", e,
                paste(sprintf("%.2f", pw[[e]][, k]), collapse = " ")))
  }
}

# ====================================================================
# 3. Save results
# ====================================================================
out <- list(pilot = "GTEx Adrenal vs Liver (paired)", B_GRID = B_GRID,
            N_GRID = N_GRID, NSIMS = NSIMS, NGENES = NGENES, FDR = FDR,
            endpoints = ENDPOINTS, power = pw, se = se,
            prop_DR = bio$prop_DR, prop_DP = bio$prop_DP, prop_DM = bio$prop_DM)
rds <- file.path(out_dir, "Fig6_differential_Bsweep.rds")
saveRDS(out, rds)
cat(sprintf("\nSaved: %s\n", rds))

# ====================================================================
# 4. Figure: 1x3 (DR | DP | DM), one curve per B
# ====================================================================
cols <- c("#440154", "#3b528b", "#21918c", "#5ec962", "#fde725")[seq_along(B_GRID)]
ep_titles <- c(DR = "Differential rhythmicity (DR)",
               DP = "Differential phase (DP)",
               DM = "Differential mesor (DM)")

pdf(file.path(out_dir, "Fig6_differential_Bsweep.pdf"), width = 12.5, height = 4.4)
par(mfrow = c(1, 3), mai = c(0.95, 0.9, 0.55, 0.18), mgp = c(2.5, 0.6, 0),
    oma = c(0, 0, 2.0, 0), cex.axis = 0.95, cex.lab = 1.15, cex.main = 1.15,
    font.main = 2)
letters_abc <- c("A", "B", "C")
for (ei in seq_along(ENDPOINTS)) {
  e <- ENDPOINTS[ei]
  ymax <- 100
  plot(NA, xlim = c(0, max(N_GRID) * 1.03), ylim = c(0, ymax),
       xlab = "Sample size (n)", ylab = "Power (%)", main = "")
  title(main = sprintf("%s   %s", letters_abc[ei], ep_titles[[e]]),
        adj = 0.5, line = 0.3, font.main = 2)
  abline(h = 80, lty = 2, col = "grey50", lwd = 1.1); grid()
  for (k in seq_along(B_GRID)) {
    y  <- 100 * pw[[e]][, k]
    s  <- 100 * se[[e]][, k]
    lo <- pmax(0, y - s); hi <- pmin(100, y + s)
    ok <- is.finite(lo) & is.finite(hi) & s > 0
    if (any(ok)) arrows(N_GRID[ok], lo[ok], N_GRID[ok], hi[ok],
                        code = 3, angle = 90, length = 0.02, col = cols[k], lwd = 0.9)
    lines(N_GRID, y, type = "o", pch = 19, lwd = 1.7, col = cols[k], cex = 0.6)
  }
  legend("bottomright", paste0("B=", B_GRID), col = cols, lty = 1, pch = 19,
         lwd = 1.5, cex = 0.7, bty = "o", box.col = "grey70", inset = 0.01,
         bg = "white", y.intersp = 0.85)
}
mtext(sprintf("Active-design B vs m trade-off for differential endpoints (GTEx Adrenal vs Liver; N_sim=%d, BH-FDR %.2f)",
              NSIMS, FDR),
      outer = TRUE, side = 3, line = 0.4, font = 2, cex = 1.05)
dev.off()
cat(sprintf("Saved: %s\n", file.path(out_dir, "Fig6_differential_Bsweep.pdf")))
cat("\n=== Done ===\n")
