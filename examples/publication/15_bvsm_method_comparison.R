#' =======================================================================
#' 15_bvsm_method_comparison.R — B vs m Tradeoff: Method Comparison
#' =======================================================================
#'
#' Single-cohort rhythmicity detection power across B and N for three
#' detection methods:
#'   1. DCP  — cosinor K=1 F-test (B-invariant under sinusoidal truth)
#'   2. JTK  — JTK_CYCLE via MetaCycle (favors replication depth)
#'   3. MH   — adaptive multi-harmonic, K=floor((B-1)/2) (favors B when alpha2>0)
#'
#' Three pilot datasets (single group):
#'   A. Mouse LIV   (GSE54651)  r~2.88  strong
#'   B. Baboon LUN  (CAMO)      r~1.72  moderate
#'   C. Mouse D1    (D1D2)      r~0.65  weak / brain
#'
#' Outputs (Figure 3: method comparison under cosinor truth):
#'   output/bvsm_method_comparison/results/results_<dataset>_<method>_a2_0.rds
#'   output/bvsm_method_comparison/results/results_all.rds
#'   output/bvsm_method_comparison/figures/fig3_method_comparison.pdf
#'
#' USAGE:
#'   Rscript examples/publication/15_bvsm_method_comparison.R
#'   SMOKE_TEST=true Rscript examples/publication/15_bvsm_method_comparison.R

# =====================================================================
# 0. Setup
# =====================================================================
SMOKE_TEST  <- identical(Sys.getenv("SMOKE_TEST"), "true")
GLOBAL_SEED <- 2025L
set.seed(GLOBAL_SEED)

old_wd <- setwd("code"); source("setup.R"); setwd(old_wd)

B_VALS      <- c(3L, 4L, 6L, 8L, 12L)
N_GRID      <- if (SMOKE_TEST) c(12L, 24L, 48L) else seq(12L, 144L, by = 12L)
ALPHA2_VALS <- if (SMOKE_TEST) c(0, 1.0) else c(0, 0.5, 0.75, 1.0)
NSIMS       <- if (SMOKE_TEST) 5L  else 100L
NGENES      <- if (SMOKE_TEST) 200L else 5000L
FDR_THRESH  <- 0.05
N_CORES     <- as.integer(Sys.getenv("MC_CORES", unset = "60"))

cat(sprintf("Mode     : %s\n", if (SMOKE_TEST) "SMOKE" else "PRODUCTION"))
cat(sprintf("B vals   : %s\n", paste(B_VALS, collapse = ", ")))
cat(sprintf("N grid   : %s\n", paste(N_GRID,  collapse = ", ")))
cat(sprintf("alpha2   : %s\n", paste(ALPHA2_VALS, collapse = ", ")))
cat(sprintf("nsims    : %d\n", NSIMS))
cat(sprintf("ngenes   : %d\n", NGENES))
cat(sprintf("mc.cores : %d\n\n", N_CORES))

out_dir <- "output/bvsm_method_comparison"
dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "results"),  recursive = TRUE, showWarnings = FALSE)

# =====================================================================
# 1. Load pilot datasets — single group from each pair
# =====================================================================
cat("--- Loading datasets ---\n")

## A: Mouse LIV (GSE54651) -------------------------------------------
dat_mouse <- readRDS("data/mice_GSE54651_CPM.RData")
prep_liv  <- prepCircadianData(dat_mouse$count_clean[["LIV"]],
                               times      = dat_mouse$tod[["LIV"]],
                               input_type = "log2")
mat_liv  <- prep_liv$data[rowSums(prep_liv$data > 0) >= 4, , drop = FALSE]
tod_liv  <- prep_liv$times
rm(dat_mouse, prep_liv)

## B: Baboon LUN (CAMO) -----------------------------------------------
load("data/CAMO_PRC_hmb.RData")
prep_lun <- prepCircadianData(baboon_withTOD$baboon[["LUN"]],
                               times      = baboon_withTOD$tod[["LUN"]],
                               input_type = "cpm")
mat_lun  <- prep_lun$data[rowSums(prep_lun$data > 0) >= 6, , drop = FALSE]
tod_lun  <- prep_lun$times
rm(baboon_withTOD, gtex, mice, prep_lun)

## C: Mouse D1 (D1D2 striatum) ----------------------------------------
pheno    <- read.csv("data/mouse_clinicalinfo_03082021_rmOutliers.csv", row.names = 1)
prep_d1  <- prepCircadianData("data/mouse_D1D2_logCPMfiltered_counts.csv",
                               times      = "time",
                               input_type = "counts",
                               pheno      = pheno,
                               sample_col = "sample")
d1_samp  <- pheno$sample[pheno$cell == "D1"]
mat_d1   <- prep_d1$data[, colnames(prep_d1$data) %in% d1_samp, drop = FALSE]
tod_d1   <- pheno$time[match(colnames(mat_d1), pheno$sample)]
mat_d1   <- mat_d1[rowSums(mat_d1 > 1) >= floor(ncol(mat_d1) / 2), , drop = FALSE]
rm(pheno, prep_d1)

datasets <- list(
  LIV = list(mat = mat_liv, tod = tod_liv, label = "Mouse LIV (r~2.88)"),
  LUN = list(mat = mat_lun, tod = tod_lun, label = "Baboon LUN (r~1.72)"),
  D1  = list(mat = mat_d1,  tod = tod_d1,  label = "Mouse D1 (r~0.65)")
)
rm(mat_liv, mat_lun, mat_d1, tod_liv, tod_lun, tod_d1)

# =====================================================================
# 2. Estimate pilot parameters and run power via unified API
# =====================================================================
cat("\n--- Estimating pilot parameters and running power ---\n")
printMethodGuidance(methods = c("DCP", "JTK", "MH"), verbose = TRUE)

all_results <- list()

for (ds_name in names(datasets)) {
  ds <- datasets[[ds_name]]
  cat(sprintf("\n==============================\n"))
  cat(sprintf("Dataset: %s\n", ds$label))
  cat(sprintf("==============================\n"))

  bio <- estCircadianParam(ds$mat, times = ds$tod, period = 24, verbose = TRUE)
  bio$ngenes <- NGENES   # use full matrix for estimation, cap simulation gene count

  design <- CircadianDesignOptions(
    sample_sizes = N_GRID,
    nsims        = NSIMS,
    design       = "active",
    cts          = seq(0, 24 * (1 - 1/B_VALS[1]), length.out = B_VALS[1]),
    B_values     = B_VALS
  )

  analysis <- CircadianAnalysisOptions(
    alpha           = FDR_THRESH,
    p.adjust.method = "BH",
    fdr_thresholds  = FDR_THRESH
  )

  # JTK: alpha2=0 only (direction is alpha2-invariant under mean-collapse)
  # DCP and MH: sweep full ALPHA2_VALS
  for (meth_name in c("DCP", "JTK", "MH")) {
    a2_vals <- if (meth_name == "JTK") 0 else ALPHA2_VALS
    cat(sprintf("\n  Method: %s  alpha2: %s\n", meth_name, paste(a2_vals, collapse = ", ")))

    set.seed(GLOBAL_SEED)
    res <- runSingleCohortPower(bio, design, analysis,
                                 methods     = meth_name,
                                 alpha2      = a2_vals,
                                 mc.cores    = N_CORES,
                                 plot        = FALSE,
                                 verbose     = FALSE)

    key <- sprintf("%s_%s", ds_name, meth_name)
    all_results[[key]] <- res$power_df

    for (a2 in a2_vals) {
      saveRDS(res$power_df[res$power_df$alpha2 == a2, ],
              file.path(out_dir, "results",
                        sprintf("results_%s_%s_a2_%.2f.rds", ds_name, meth_name, a2)))
    }
  }
}

# Combined tidy frame
res_df <- do.call(rbind, all_results)
saveRDS(res_df, file.path(out_dir, "results", "results_all.rds"))
cat("\nSaved: results_all.rds\n")

# =====================================================================
# 3. Figure 3: Method comparison under cosinor truth (alpha2=0)
# =====================================================================
library(ggplot2)

b_colors <- c("3"="#E41A1C","4"="#FF7F00","6"="#4DAF4A","8"="#377EB8","12"="#984EA3")
b_labels <- setNames(paste0("B=", B_VALS), as.character(B_VALS))

cat("\nGenerating Figure 3...\n")

fig3_dat <- res_df[res_df$alpha2 == 0, ]
fig3_dat$method_label <- factor(fig3_dat$method,
  levels = c("DCP", "JTK", "MH"),
  labels = c("DCP (K=1 cosinor)", "JTK_CYCLE", "Multi-harmonic (adaptive K)"))
fig3_dat$dataset_label <- factor(fig3_dat$dataset %||% {
  # recover dataset from key if not stored
  gsub("_(DCP|JTK|MH)$", "", names(all_results)[match(fig3_dat$method, fig3_dat$method)])
}, levels = c("LIV", "LUN", "D1"),
  labels = c("Mouse LIV\n(r~2.88, strong)",
             "Baboon LUN\n(r~1.72, moderate)",
             "Mouse D1\n(r~0.65, weak)"))
fig3_dat$B_fac <- factor(fig3_dat$B, levels = B_VALS)

p3 <- ggplot(fig3_dat, aes(x = N, y = 100 * power, colour = B_fac, group = B_fac)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  geom_hline(yintercept = 80, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  facet_grid(method_label ~ dataset_label) +
  scale_colour_manual(values = b_colors, labels = b_labels, name = NULL) +
  scale_x_continuous(breaks = seq(12, 144, by = 24)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  labs(
    x     = "Total sample size N",
    y     = "Detection power (%)",
    title = "B vs m tradeoff: method comparison under cosinor truth"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey92"),
    legend.position  = "bottom",
    panel.grid.minor = element_blank()
  )

fig3_path <- file.path(out_dir, "figures", "fig3_method_comparison.pdf")
ggsave(fig3_path, p3, width = 11, height = 9)
cat(sprintf("Saved: %s\n", fig3_path))

cat("\n=== Done ===\n")
