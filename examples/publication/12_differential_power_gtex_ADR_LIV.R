#' =======================================================================
#' Differential Circadian Power Analysis — GTEx Adrenal Gland vs Liver
#' =======================================================================
#'
#' WHAT THIS DOES:
#'   Simulation-based power analysis for transcriptome-wide differential
#'   circadian detection (DR, DP, DM) for:
#'     Adrenal Gland vs Liver (GTEx v8): strong DR + DP + DM, n=187/154
#'
#' OUTPUTS:
#'   output/differential/results/diff_power_ADR_vs_LIV_<ts>.rds
#'   output/differential/figures/diff_power_fig_ADR_vs_LIV_<ts>.pdf
#'
#' USAGE:
#'   cd /path/to/PowerSim
#'   Rscript examples/publication/12_differential_power_gtex_ADR_LIV.R
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

out_dir_fig <- "output/differential/figures"
out_dir_res <- "output/differential/results"
dir.create(out_dir_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_res, recursive = TRUE, showWarnings = FALSE)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
cat(sprintf("\n=== GTEx ADR vs LIV Differential Power [%s] ===\n", timestamp))

# =====================================================================
# 1. Simulation grid
# =====================================================================
if (SMOKE_TEST) {
  sample_sizes <- c(30, 60, 100)
  nsims        <- 5L
  n_cores      <- 1L
} else {
  sample_sizes <- c(20, 30, 40, 50, 60, 80, 100, 120, 150, 200)
  nsims        <- 200L
  n_cores      <- as.integer(Sys.getenv("MC_CORES", unset = "4"))
}

cat(sprintf("Global seed  : %d\n", GLOBAL_SEED))
cat(sprintf("Sample sizes : %s\n", paste(sample_sizes, collapse = ", ")))
cat(sprintf("nsims        : %d\n", nsims))
cat(sprintf("mc.cores     : %d\n", n_cores))

# =====================================================================
# 2. Load GTEx CPM and extract Adrenal / Liver
# =====================================================================
gtex_path <- "/home/qtp1/Projects/Collaborative/GTEXdata/TOD Xiangning/data/CPM/CPM.all.norm.RData"
cat(sprintf("\nLoading GTEx CPM: %s\n", gtex_path))
load(gtex_path)

extract_tissue <- function(tissue_name) {
  df   <- CPM.all.norm[[tissue_name]]
  ids  <- as.character(colnames(df))
  hhmm <- sapply(strsplit(ids, "\\."), function(x) if (length(x) >= 3) x[3] else NA)
  hrs  <- as.numeric(substr(hhmm, 1, 2)) + as.numeric(substr(hhmm, 3, 4)) / 60
  ok   <- !is.na(hrs)
  list(expr = as.matrix(df[, ok]), times = hrs[ok], n = sum(ok))
}

adr <- extract_tissue("Adrenal Gland")
liv <- extract_tissue("Liver")
common_genes <- intersect(rownames(adr$expr), rownames(liv$expr))

cat(sprintf("Adrenal Gland : n=%d samples, %d genes\n", adr$n, nrow(adr$expr)))
cat(sprintf("Liver         : n=%d samples, %d genes\n", liv$n, nrow(liv$expr)))
cat(sprintf("Common genes  : %d\n", length(common_genes)))

expr_adr <- adr$expr[common_genes, ]
expr_liv <- liv$expr[common_genes, ]

# =====================================================================
# 3. Estimate pilot parameters
# =====================================================================
pilot_rds <- file.path("data", "gtex_adr_vs_liv_pilot.rds")

if (file.exists(pilot_rds)) {
  cat(sprintf("\nLoading existing pilot: %s\n", pilot_rds))
  bio <- readRDS(pilot_rds)
} else {
  cat("\n--- Estimating two-group pilot parameters ---\n")
  set.seed(GLOBAL_SEED)
  bio <- estCircadianParamTwoGroup(
    data_1  = expr_adr, data_2  = expr_liv,
    times_1 = adr$times, times_2 = liv$times,
    verbose = TRUE
  )
  saveRDS(bio, pilot_rds)
  cat(sprintf("Pilot saved -> %s\n", pilot_rds))
}

cat(sprintf("\nPilot: n1=%d  n2=%d  DR=%.3f  DP=%.3f  DM=%.3f\n",
            length(bio$cts), length(bio$cts2),
            bio$prop_DR, bio$prop_DP, bio$prop_DM))

# =====================================================================
# 4. Run simulation
# =====================================================================
design <- CircadianDesignOptions(
  sample_sizes = sample_sizes,
  nsims        = nsims,
  design       = "passive",
  cts          = bio$cts,
  test_types   = c("DR", "DP", "DM")
)

analysis <- CircadianAnalysisOptions(
  alpha           = 0.05,
  p.adjust.method = "BH"
)

rds_path <- file.path(out_dir_res,
  sprintf("diff_power_ADR_vs_LIV_%s.rds", timestamp))
fig_path <- file.path(out_dir_fig,
  sprintf("diff_power_fig_ADR_vs_LIV_%s.pdf", timestamp))

set.seed(GLOBAL_SEED)
res <- runDifferentialPower(bio, design, analysis,
                             methods     = "DCP",
                             test_types  = c("DR", "DP", "DM"),
                             mc.cores    = n_cores,
                             plot        = TRUE,
                             output_file = fig_path,
                             verbose     = TRUE)

saveRDS(res, rds_path)
cat(sprintf("Results saved -> %s\n", rds_path))
cat(sprintf("Figure saved  -> %s\n", fig_path))

cat("\n=== Done ===\n")
cat(sprintf("Results : %s\n", rds_path))
cat(sprintf("Figure  : %s\n", fig_path))
