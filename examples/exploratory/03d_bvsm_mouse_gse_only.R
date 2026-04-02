#' B vs m bootstrap for Mouse LIV vs CER (GSE54651) — B vs m section only.
#' Saves bm_boot_grid.rds and updates paper/PowerSim/figures/bm_tradeoff.pdf.

SMOKE <- nzchar(Sys.getenv("POWERSIM_SMOKE"))
NBOOT       <- if (SMOKE) 5L  else 50L
NSIMS_INNER <- if (SMOKE) 5L  else 20L
NCORES      <- if (SMOKE) 2L  else 20L
NGENES      <- if (SMOKE) 500L else 5000L
RHYTHM_PVAL <- 0.05

B_VALS    <- c(4L, 6L, 8L)
N_GRID_BM <- if (SMOKE) c(24L, 48L) else c(24L, 48L, 72L, 96L)

POWERSIM_ROOT <- {
  env <- Sys.getenv("POWERSIM_ROOT", unset = "")
  if (nzchar(env)) env else
    "/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim"
}
setwd(POWERSIM_ROOT)
source_dir <- file.path(getwd(), "code")
old_wd <- setwd(source_dir); source("setup.R"); setwd(old_wd)

cat(sprintf("Mode: %s | NBOOT=%d | NSIMS_INNER=%d | NCORES=%d\n",
            if (SMOKE) "SMOKE" else "PRODUCTION", NBOOT, NSIMS_INNER, NCORES))
cat(sprintf("B values: %s\nN grid:   %s\n\n",
            paste(B_VALS, collapse=", "), paste(N_GRID_BM, collapse=", ")))

opts_analysis <- CircadianAnalysisOptions(
  alpha = 0.05, p.adjust.method = "BH",
  fdr_thresholds = c(0.01, 0.05, 0.10, 0.20),
  reference_n = 24L
)

# Load data
cat("Loading Mouse LIV vs CER (GSE54651)...\n")
dat_mouse <- readRDS("data/mice_GSE54651_CPM.RData")
prep_liv  <- prepCircadianData(dat_mouse$count_clean[["LIV"]],
                                times = dat_mouse$tod[["LIV"]], input_type = "log2")
prep_cer  <- prepCircadianData(dat_mouse$count_clean[["CER"]],
                                times = dat_mouse$tod[["CER"]], input_type = "log2")
rm(dat_mouse)

tod_liv  <- prep_liv$times
tod_cer  <- prep_cer$times
keep_liv <- rowSums(prep_liv$data > 0) >= 4
keep_cer <- rowSums(prep_cer$data > 0) >= 4
common_g <- intersect(rownames(prep_liv$data)[keep_liv],
                       rownames(prep_cer$data)[keep_cer])
set.seed(7)
g_idx   <- sample(common_g, min(NGENES, length(common_g)))
mat_liv <- prep_liv$data[g_idx, , drop = FALSE]
mat_cer <- prep_cer$data[g_idx, , drop = FALSE]
rm(prep_liv, prep_cer)
cat(sprintf("  LIV: %d genes x %d samples | ZT: %s\n",
            nrow(mat_liv), ncol(mat_liv),
            paste(sort(unique(tod_liv)), collapse=", ")))

# Estimate parameters
cat("Estimating parameters...\n")
opts_bio <- estCircadianParamTwoGroup(
  data_1 = mat_liv, data_2 = mat_cer,
  times_1 = tod_liv, times_2 = tod_cer,
  period = 24, min_rhythm_pval = RHYTHM_PVAL, verbose = FALSE
)
opts_bio <- updateBioOptions(opts_bio, ngenes = NGENES)

fit_liv   <- fitCosinorAll(mat_liv, times = tod_liv, period = 24)
rhy_liv   <- !is.na(fit_liv$pvalue) & fit_liv$pvalue < RHYTHM_PVAL
r_liv_med <- median(as.numeric(fit_liv$A[rhy_liv]) /
                    as.numeric(fit_liv$sigma[rhy_liv]), na.rm = TRUE)
prop_DR   <- opts_bio$prop_DR
cat(sprintf("  r_median=%.2f  prop_DR=%.1f%%\n\n", r_liv_med, prop_DR * 100))

# Bootstrap B vs m grid
design_vec <- sort(unique(tod_liv))
boot_opts  <- CircadianBootstrapOptions(
  design        = "active",
  design_vector = design_vec,
  B_values      = B_VALS,
  N_values      = N_GRID_BM,
  nboot         = NBOOT,
  nsims_inner   = NSIMS_INNER,
  seed          = 42L
)
bio_bm <- updateBioOptions(opts_bio,
  prop_DR = prop_DR, prop_DP = 0, prop_DA = 0,
  phase_diff = c(0, 0), amp_diff = c(1, 1)
)

cat("Running bootstrap B vs m grid...\n")
bm_boot <- runBootstrapDesignGrid(
  pilot_data    = mat_liv,
  pilot_times   = tod_liv,
  boot.opts     = boot_opts,
  analysis.opts = opts_analysis,
  bio_diff.opts = bio_bm,
  verbose       = TRUE,
  mc.cores      = NCORES
)

# Save RDS
rds_out <- "output/03d_mouse_gse_bm_boot.rds"
saveRDS(bm_boot, rds_out)
cat(sprintf("\nSaved: %s\n", rds_out))

summaryBootstrapDesignGrid(bm_boot, test_type = "DR", fdr_threshold = 0.05)

# n80 per B
n_N <- length(bm_boot$N_values); n_B <- length(bm_boot$B_values)
pm  <- matrix(bm_boot$power_mean[,,1], nrow = n_N, ncol = n_B)
lo  <- matrix(bm_boot$power_ci_lo[,,1], nrow = n_N, ncol = n_B)
hi  <- matrix(bm_boot$power_ci_hi[,,1], nrow = n_N, ncol = n_B)
cat("\n--- n80 per B (bootstrap, FDR 5%) ---\n")
for (b in seq_len(n_B)) {
  n80_b  <- N_GRID_BM[which(pm[,b] >= 0.80)[1]]
  n80_lo <- N_GRID_BM[which(hi[,b] >= 0.80)[1]]
  n80_hi <- N_GRID_BM[which(lo[,b] >= 0.80)[1]]
  cat(sprintf("  B=%d: n80=%s  [95%% CI: %s, %s]\n",
              bm_boot$B_values[b],
              ifelse(is.na(n80_b), ">max", n80_b),
              ifelse(is.na(n80_lo), ">max", n80_lo),
              ifelse(is.na(n80_hi), ">max", n80_hi)))
}

# Regenerate combined bm_tradeoff.pdf
r3b <- readRDS("output/03b_power_core_active_20260331_1159/bm_boot_grid.rds")
r3c <- readRDS("output/03c_power_core_mouse_20260331_1215/bm_boot_grid.rds")

fig_out <- "paper/PowerSim/figures/bm_tradeoff.pdf"
pdf(fig_out, width = 17, height = 5)
par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3.5, 1.5))

plot_bm_panel <- function(bm, title, col_hi, n80_line = TRUE) {
  N   <- bm$N_values; B <- bm$B_values
  pm  <- matrix(bm$power_mean[,,1], nrow=length(N), ncol=length(B))
  lo  <- matrix(bm$power_ci_lo[,,1], nrow=length(N), ncol=length(B))
  hi  <- matrix(bm$power_ci_hi[,,1], nrow=length(N), ncol=length(B))
  cols <- colorRampPalette(c("gray70", col_hi))(length(B))
  pchs <- c(16, 17, 15)[seq_along(B)]
  plot(NA, xlim=range(N), ylim=c(0,1), las=1, xaxt="n",
       xlab="N per group", ylab="DR power (FDR 5%)", main=title)
  axis(1, at=N)
  abline(h=0.80, lty=2, col="gray40")
  abline(h=c(0.2, 0.4, 0.6), lty=3, col="gray85")
  for (bi in seq_along(B)) {
    polygon(c(N,rev(N)), c(lo[,bi],rev(hi[,bi])),
            col=adjustcolor(cols[bi],0.15), border=NA)
    lines(N,  pm[,bi], col=cols[bi], lwd=2, lty=bi)
    points(N, pm[,bi], col=cols[bi], pch=pchs[bi], cex=0.9)
  }
  m_approx <- pmax(1L, round(min(N)/B))
  legend("bottomright", legend=sprintf("B=%d (m~%d/ZT)", B, m_approx),
         col=cols, lwd=2, lty=seq_along(B), pch=pchs, bty="n", cex=0.85)
  mtext("(bootstrap 95% CI)", side=3, line=0, cex=0.7, col="gray50")
}

plot_bm_panel(bm_boot,
  title = sprintf("Mouse LIV vs CER (r=%.1f)\nprop_DR=%.0f%%, B=4/6/8", r_liv_med, prop_DR*100),
  col_hi = "steelblue")
plot_bm_panel(r3b,
  title = "Baboon LUN vs CER (r=1.72)\nprop_DR=41%, B=4/6/12",
  col_hi = "darkorange")
plot_bm_panel(r3c,
  title = "Mouse D1 vs D2 neurons (r=0.65)\nprop_DR=20%, B=3/6",
  col_hi = "forestgreen")

dev.off()
cat(sprintf("\nUpdated figure: %s\n", fig_out))
cat("Done.\n")
