#' ============================================================
#' fig2_panelA_only.R - Fig 2 trimmed to 1x3 (Panel A only)
#'
#' Replots the differential ADR vs LIV power figure as a single
#' row of three panels (DR, DP, DM), each showing marginal power
#' vs sample size at four FDR thresholds. Reuses the cached
#' differential simulation results - no new sims needed.
#'
#' Output:
#'   submission/figures/Fig2_differential_ADR_vs_LIV.pdf
#'   output/main_figures/Fig2_differential_ADR_vs_LIV.pdf
#' ============================================================

setwd("/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")
old_wd <- setwd("code"); source("setup.R"); setwd(old_wd)

pick_latest <- function(dir, pattern) {
  cands <- list.files(dir, pattern, full.names = TRUE)
  if (length(cands) == 0) stop("No match for ", pattern, " in ", dir)
  cands[order(file.mtime(cands), decreasing = TRUE)[1]]
}

rds_path <- pick_latest("output/differential/results",
                        "^diff_power_ADR_vs_LIV.*\\.rds$")
cat(sprintf("Differential RDS: %s\n", basename(rds_path)))
res <- readRDS(rds_path)

DISPLAY_N      <- c(20, 40, 60, 80, 100, 120, 150, 200)
FDR_THRESHOLDS <- c(0.01, 0.05, 0.10, 0.20)
THRESH_COLS    <- c("darkgreen", "steelblue", "orange", "red")
THRESH_LABELS  <- paste0("FDR ", round(100 * FDR_THRESHOLDS), "%")
ENDPOINTS      <- c("DR", "DP", "DM")
VLINE_POWER    <- 0.80
VLINE_FDR      <- 0.20

# ---- Compute marginal power per endpoint ----
n_sizes <- length(res$sample_sizes)
nsims   <- res$nsims
ngenes  <- res$ngenes

compute_marginal <- function(ep) {
  target_types <- switch(ep, DR = c(2L, 3L), DP = 4L, DM = 5L)
  fdr_mat <- res[[paste0("fdr_", ep)]]
  marginal_sim <- array(NA_real_, dim = c(n_sizes, length(FDR_THRESHOLDS), nsims))
  for (s in seq_len(nsims)) {
    dt <- res$diff_type[[s]]
    is_tgt <- dt %in% target_types
    n_tgt  <- sum(is_tgt)
    for (j in seq_len(n_sizes)) {
      fdr_g <- fdr_mat[, j, s]
      for (t in seq_along(FDR_THRESHOLDS)) {
        disc <- !is.na(fdr_g) & fdr_g <= FDR_THRESHOLDS[t]
        marginal_sim[j, t, s] <- if (n_tgt > 0) sum(disc & is_tgt) / n_tgt else NA_real_
      }
    }
  }
  list(
    mean = apply(marginal_sim, c(1, 2), mean, na.rm = TRUE),
    se   = apply(marginal_sim, c(1, 2),
                 function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
  )
}

cat("Computing marginal power for DR, DP, DM...\n")
marg <- lapply(ENDPOINTS, compute_marginal)
names(marg) <- ENDPOINTS

# ---- N80 per endpoint ----
n80 <- sapply(ENDPOINTS, function(ep) {
  np <- npower(res, target_power = VLINE_POWER, fdr = VLINE_FDR, endpoint = ep)
  if (is.na(np$n)) NA_real_ else np$n
})
cat(sprintf("N80 (@ FDR %g%%): DR=%s, DP=%s, DM=%s\n",
            VLINE_FDR * 100,
            ifelse(is.na(n80["DR"]), "NA", n80["DR"]),
            ifelse(is.na(n80["DP"]), "NA", n80["DP"]),
            ifelse(is.na(n80["DM"]), "NA", n80["DM"])))

disp_idx <- which(res$sample_sizes %in% DISPLAY_N)
ss_disp  <- res$sample_sizes[disp_idx]

# ---- Render 1x3 ----
out_pdf <- tempfile(fileext = ".pdf")
pdf(out_pdf, width = 13.5, height = 4.5)
par(mfrow = c(1, 3), mai = c(1.05, 1.05, 0.85, 0.20),
    mgp = c(3.2, 0.65, 0), oma = c(0, 0, 2.4, 0),
    cex.axis = 1.25, cex.lab = 1.45, font.main = 2, cex.main = 1.40)

for (ei in seq_along(ENDPOINTS)) {
  ep   <- ENDPOINTS[ei]
  M    <- marg[[ep]]$mean
  SE   <- marg[[ep]]$se
  letter <- LETTERS[ei]

  matplot(ss_disp, 100 * M[disp_idx, , drop = FALSE],
          type = "b", pch = 19, lwd = 2.2,
          col = THRESH_COLS, lty = 1,
          xlim = c(0, max(ss_disp) * 1.05), ylim = c(0, 100),
          xlab = "Sample size (n)", ylab = sprintf("%s Power (%%)", ep),
          main = "")
  title(main = sprintf("%s   %s - Power vs Sample Size", letter, ep),
        adj = 0, font.main = 2, cex.main = 1.30, line = 0.7)

  for (t in seq_along(FDR_THRESHOLDS)) {
    arrows(ss_disp, 100 * (M[disp_idx, t] - 1.96 * SE[disp_idx, t]),
           ss_disp, 100 * (M[disp_idx, t] + 1.96 * SE[disp_idx, t]),
           code = 3, angle = 90, length = 0.03,
           col = THRESH_COLS[t], lwd = 1.2)
  }

  abline(h = 80, lty = 2, col = "grey50", lwd = 1.3)
  vline_n <- n80[ep]
  if (!is.na(vline_n) && is.finite(vline_n)) {
    abline(v = vline_n, lty = 2, col = adjustcolor("steelblue", 0.7), lwd = 1.8)
    text(vline_n, 95, sprintf("n = %d", vline_n),
         col = "steelblue", cex = 0.95, adj = c(1.05, 0.5), font = 2)
  }
  grid()

  if (ei == 1L) {
    legend("bottomright", THRESH_LABELS,
           col = THRESH_COLS, lty = 1, pch = 19, lwd = 2.2,
           cex = 0.66, bty = "o", box.col = "grey70", box.lwd = 0.5,
           inset = 0.01, y.intersp = 0.88, bg = "white")
  }
}

mtext("Differential Circadian Power Analysis (GTEx Adrenal Gland vs Liver)",
      outer = TRUE, side = 3, line = 0.6, font = 2, cex = 1.4)
dev.off()

for (d in c("submission/figures/Fig2_differential_ADR_vs_LIV.pdf",
            "output/main_figures/Fig2_differential_ADR_vs_LIV.pdf")) {
  dir.create(dirname(d), recursive = TRUE, showWarnings = FALSE)
  file.copy(out_pdf, d, overwrite = TRUE)
  cat(sprintf("Saved: %s\n", d))
}

cat("Done.\n")
