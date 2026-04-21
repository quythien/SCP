#' =======================================================================
#' Single-Cohort Circadian Power Analysis — GTEx Adrenal Gland & Liver
#' =======================================================================
#'
#' WHAT THIS DOES:
#'   Figure 1 equivalent for GTEx Adrenal Gland and Liver, showing
#'   single-cohort rhythmicity power under passive design with n up to 200.
#'
#' OUTPUTS:
#'   output/single_cohort/figures/single_cohort_power_GTEx_AdrenalGland.pdf
#'   output/single_cohort/figures/single_cohort_power_GTEx_Liver.pdf
#'   output/single_cohort/results/single_cohort_power_GTEx_AdrenalGland_<ts>.rds
#'   output/single_cohort/results/single_cohort_power_GTEx_Liver_<ts>.rds
#'
#' USAGE:
#'   cd /path/to/PowerSim
#'   Rscript examples/publication/14_single_cohort_gtex_ADR_LIV.R
#'
#' @author Thien Pham

# =====================================================================
# 0. Setup
# =====================================================================
SMOKE_TEST  <- identical(Sys.getenv("SMOKE_TEST"), "true")
GLOBAL_SEED <- 2025L
set.seed(GLOBAL_SEED)

old_wd <- setwd("code")
source("setup.R")
setwd(old_wd)

out_dir_fig <- "output/single_cohort/figures"
out_dir_res <- "output/single_cohort/results"
dir.create(out_dir_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_res, recursive = TRUE, showWarnings = FALSE)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
cat(sprintf("\n=== GTEx Single-Cohort Power [%s] ===\n", timestamp))

# =====================================================================
# 1. Simulation grid
# =====================================================================
if (SMOKE_TEST) {
  sample_sizes <- c(20, 40, 60, 80, 100)
  nsims        <- 5L
  n_cores      <- 1L
} else {
  sample_sizes <- c(20, 30, 40, 50, 60, 70, 80, 100, 120, 150, 200)
  nsims        <- 200L
  n_cores      <- as.integer(Sys.getenv("MC_CORES", unset = "4"))
}

cat(sprintf("Sample sizes : %s\n", paste(sample_sizes, collapse = ", ")))
cat(sprintf("nsims        : %d\n", nsims))
cat(sprintf("mc.cores     : %d\n", n_cores))

# =====================================================================
# 2. Load GTEx CPM
# =====================================================================
gtex_path <- "/home/qtp1/Projects/Collaborative/GTEXdata/TOD Xiangning/data/CPM/CPM.all.norm.RData"
cat(sprintf("\nLoading GTEx: %s\n", gtex_path))
load(gtex_path)

extract_tissue <- function(tissue_name) {
  df   <- CPM.all.norm[[tissue_name]]
  ids  <- as.character(colnames(df))
  hhmm <- sapply(strsplit(ids, "\\."), function(x) if (length(x) >= 3) x[3] else NA)
  hrs  <- as.numeric(substr(hhmm, 1, 2)) + as.numeric(substr(hhmm, 3, 4)) / 60
  ok   <- !is.na(hrs)
  list(expr = as.matrix(df[, ok]), times = hrs[ok], n = sum(ok))
}

tissues <- list(
  AdrenalGland = list(gtex_name = "Adrenal Gland",  label = "GTEx Adrenal Gland"),
  Liver        = list(gtex_name = "Liver",           label = "GTEx Liver")
)

tissue_filter <- Sys.getenv("TISSUE", unset = "")
if (nchar(tissue_filter) > 0 && tissue_filter %in% names(tissues))
  tissues <- tissues[tissue_filter]

# =====================================================================
# 3. Loop over tissues
# =====================================================================
for (tname in names(tissues)) {
  ti <- tissues[[tname]]
  cat(sprintf("\n====================================================\n"))
  cat(sprintf("  %s\n", ti$label))
  cat(sprintf("====================================================\n"))

  t <- extract_tissue(ti$gtex_name)
  cat(sprintf("n=%d samples, %d genes\n", t$n, nrow(t$expr)))

  # Estimate or load pilot
  pilot_rds <- file.path("data", sprintf("gtex_%s_single_pilot.rds", tname))
  if (file.exists(pilot_rds)) {
    cat(sprintf("Loading cached pilot: %s\n", pilot_rds))
    bio <- readRDS(pilot_rds)
  } else {
    cat("Estimating pilot parameters...\n")
    set.seed(GLOBAL_SEED)
    bio <- estCircadianParam(t$expr, t$times, verbose = TRUE)
    saveRDS(bio, pilot_rds)
    cat(sprintf("Pilot saved: %s\n", pilot_rds))
  }

  cat(sprintf("Pilot: n_genes=%d  prop_rhythmic=%.1f%%  n_pilot=%d\n",
              bio$ngenes, 100 * bio$prop_rhythmic, length(bio$cts)))

  design <- CircadianDesignOptions(
    sample_sizes = sample_sizes,
    nsims        = nsims,
    design       = "passive",
    cts          = bio$cts
  )

  analysis <- CircadianAnalysisOptions(
    alpha           = 0.05,
    p.adjust.method = "BH",
    r_strata        = makeAdaptiveRStrata(bio, bin_width = 0.25)
  )

  fig_path <- file.path(out_dir_fig,
    sprintf("single_cohort_power_GTEx_%s.pdf", tname))
  rds_path <- file.path(out_dir_res,
    sprintf("single_cohort_power_GTEx_%s_%s.rds", tname, timestamp))

  set.seed(GLOBAL_SEED)
  res <- runSingleCohortPower(bio, design, analysis,
                               methods     = "DCP",
                               mc.cores    = n_cores,
                               plot        = TRUE,
                               output_file = fig_path,
                               verbose     = TRUE)

  saveRDS(res, rds_path)
  cat(sprintf("Results saved -> %s\n", rds_path))
}

cat("\n=== Done ===\n")
