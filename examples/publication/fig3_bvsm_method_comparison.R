#' =======================================================================
#' fig3_bvsm_method_comparison.R — B vs m Tradeoff: DCP B-invariance
#' =======================================================================
#'
#' Single-cohort rhythmicity detection power across B and N for DCP only,
#' verifying B-invariance under sinusoidal truth (alpha2=0).
#'   DCP — cosinor K=1 F-test (B-invariant under sinusoidal truth)
#'
#' Three pilot datasets (single group):
#'   A. Mouse LIV   (GSE54651)  r~2.88  strong
#'   B. Baboon LUN  (CAMO)      r~1.72  moderate
#'   C. Mouse D1    (D1D2)      r~0.65  weak / brain
#'
#' Outputs (Figure 3: DCP B-invariance across three datasets):
#'   output/bvsm_method_comparison/results/results_<dataset>_DCP_a2_<val>.rds
#'   output/bvsm_method_comparison/results/results_all.rds
#'   output/bvsm_method_comparison/figures/fig3_method_comparison.pdf
#'
#' USAGE:
#'   Rscript examples/publication/fig3_bvsm_method_comparison.R
#'   SMOKE_TEST=true Rscript examples/publication/fig3_bvsm_method_comparison.R
#'   DATASET=LUN METHOD=DCP MC_CORES=40 Rscript examples/publication/fig3_bvsm_method_comparison.R

# =====================================================================
# 0. Setup
# =====================================================================
SMOKE_TEST  <- identical(Sys.getenv("SMOKE_TEST"), "true")
GLOBAL_SEED <- 2025L
set.seed(GLOBAL_SEED)

old_wd <- setwd("code"); source("setup.R"); setwd(old_wd)

B_VALS      <- c(3L, 4L, 6L, 8L, 12L)
N_GRID      <- if (SMOKE_TEST) c(12L, 24L, 48L) else seq(12L, 96L, by = 12L)
ALPHA2_VALS <- if (SMOKE_TEST) c(0, 0.5, 1.0) else c(0, 0.25, 0.5, 0.75, 1.0)
NSIMS       <- if (SMOKE_TEST) 5L  else 30L
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

DATASET_FILTER <- Sys.getenv("DATASET", unset = "ALL")
METHOD_FILTER  <- Sys.getenv("METHOD",  unset = "ALL")
if (DATASET_FILTER != "ALL") {
  datasets <- datasets[names(datasets) == DATASET_FILTER]
  cat(sprintf("Dataset filter : %s\n", DATASET_FILTER))
}
if (METHOD_FILTER != "ALL") cat(sprintf("Method filter  : %s\n\n", METHOD_FILTER))

# =====================================================================
# 2. Estimate pilot parameters and run power via unified API
# =====================================================================
cat("\n--- Estimating pilot parameters and running power ---\n")
printMethodGuidance(methods = "DCP", verbose = TRUE)

all_results <- list()

for (ds_name in names(datasets)) {
  ds <- datasets[[ds_name]]
  cat(sprintf("\n==============================\n"))
  cat(sprintf("Dataset: %s\n", ds$label))
  cat(sprintf("==============================\n"))

  bio <- estCircadianParam(ds$mat, times = ds$tod, period = 24, verbose = TRUE)
  bio$ngenes <- NGENES   # use full matrix for estimation, cap simulation gene count

  analysis <- CircadianAnalysisOptions(
    alpha           = FDR_THRESH,
    p.adjust.method = "BH",
    fdr_thresholds  = FDR_THRESH
  )

  # DCP: sweep full ALPHA2_VALS
  design <- CircadianDesignOptions(
    sample_sizes = N_GRID, nsims = NSIMS, design = "active",
    cts = seq(0, 24 * (1 - 1/B_VALS[1]), length.out = B_VALS[1]),
    B_values = B_VALS
  )

  for (meth_name in "DCP") {
    if (METHOD_FILTER != "ALL" && meth_name != METHOD_FILTER) next
    a2_vals <- ALPHA2_VALS
    cat(sprintf("\n  Method: %s  alpha2: %s\n", meth_name, paste(a2_vals, collapse = ", ")))

    set.seed(GLOBAL_SEED)
    res <- runSingleCohortGrid(bio, design, analysis,
                                methods  = meth_name,
                                alpha2   = a2_vals,
                                mc.cores = N_CORES,
                                verbose  = FALSE)
    res$power_df$dataset <- ds_name

    key <- sprintf("%s_%s", ds_name, meth_name)
    all_results[[key]] <- res$power_df

    for (a2 in a2_vals) {
      saveRDS(res$power_df[res$power_df$alpha2 == a2, ],
              file.path(out_dir, "results",
                        sprintf("results_%s_%s_a2_%.2f.rds", ds_name, meth_name, a2)))
    }
  }
}

# =====================================================================
# 2b. Skip figure generation for partial runs
# =====================================================================
if (DATASET_FILTER != "ALL" || METHOD_FILTER != "ALL") {
  cat("\nPartial run — skipping figure generation.\n")
  cat("\n=== Done ===\n")
  quit(save = "no")
}

# Combined tidy frame
res_df <- do.call(rbind, all_results)
saveRDS(res_df, file.path(out_dir, "results", "results_all.rds"))
cat("\nSaved: results_all.rds\n")

# =====================================================================
# 3. Figure 3: DCP B-invariance across three datasets (alpha2=0)
# =====================================================================
library(ggplot2)

b_colors <- c("3"="#E41A1C","4"="#FF7F00","6"="#4DAF4A","8"="#377EB8","12"="#984EA3")
b_labels <- setNames(paste0("B=", B_VALS), as.character(B_VALS))

cat("\nGenerating Figure 3...\n")

fig3_dat <- res_df[res_df$alpha2 == 0, ]
fig3_dat$dataset_label <- factor(fig3_dat$dataset,
  levels = c("LIV", "LUN", "D1"),
  labels = c("Mouse LIV\n(r~2.88, strong)",
             "Baboon LUN\n(r~1.72, moderate)",
             "Mouse D1\n(r~0.65, weak)"))
fig3_dat$B_fac <- factor(fig3_dat$B, levels = B_VALS)

p3 <- ggplot(fig3_dat, aes(x = N, y = 100 * power, colour = B_fac, group = B_fac)) +
  geom_errorbar(aes(ymin = 100*(power - power_se), ymax = 100*(power + power_se)),
                width = 1.5, linewidth = 0.4) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  geom_hline(yintercept = 80, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  facet_wrap(~ dataset_label, nrow = 1) +
  scale_colour_manual(values = b_colors, labels = b_labels, name = "Time bins (B)") +
  scale_fill_manual(values   = b_colors, labels = b_labels, name = "Time bins (B)") +
  scale_x_continuous(breaks = seq(12, 96, by = 12)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  labs(
    x     = "Total sample size N",
    y     = "Power (%)",
    title = "B vs m trade-off: DCP power under cosinor truth"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey92"),
    strip.text       = element_text(size = 10),
    legend.position  = "bottom",
    panel.grid.minor = element_blank(),
    plot.title       = element_text(hjust = 0.5, face = "bold")
  )

fig3_path <- file.path(out_dir, "figures", "fig3_method_comparison.pdf")
ggsave(fig3_path, p3, width = 13, height = 5)
cat(sprintf("Saved: %s\n", fig3_path))

cat("\n=== Done ===\n")
