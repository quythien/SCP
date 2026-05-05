#' =======================================================================
#' fig4_fmm_violation.R — FMM Waveform Robustness: Active vs Passive
#' =======================================================================
#'
#' 6-row x 3-column layout (rows 1-2: omega sweep, 3-4: beta sweep, 5-6: alpha sweep)
#' x 2 detection methods (DCP and FMM-LRT).
#'
#' Row 1/3/5 (Active, B=12, every 2h):
#'   A. Mouse LIV (GSE54651)  r~2.88  strong
#'   B. Baboon LUN (CAMO)     r~1.72  moderate
#'   C. Mouse D1 (D1D2)       r~0.65  weak
#'
#' Row 2/4/6 (Passive — KDE-sampled TOD from pilot):
#'   D. GTEx Adrenal Gland       r~1.03  strong passive
#'   E. Putamen (Kyle/GSE160521) r~0.66  moderate passive
#'   F. NAc (Kyle/GSE160521)     r~0.69  weak passive
#'
#' Detection: DCP (correct for cosinor truth) and FMM-LRT (correct for FMM truth).
#' FMM-LRT uses a pre-computed empirical null table for calibrated p-values.
#'
#' USAGE:
#'   Rscript examples/publication/fig4_fmm_violation.R
#'   SMOKE_TEST=true Rscript examples/publication/fig4_fmm_violation.R

SMOKE_TEST  <- identical(Sys.getenv("SMOKE_TEST"), "true")
GLOBAL_SEED <- 2025L
set.seed(GLOBAL_SEED)

old_wd <- setwd("code"); source("setup.R"); setwd(old_wd)
library(FMM)
suppressPackageStartupMessages(library(readxl))

OMEGA_VALS <- if (SMOKE_TEST) c(0.0, 0.5, 1.0) else c(0.0, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0)
# All N divisible by B_ACTIVE=12 so active design uses every cell
N_GRID     <- if (SMOKE_TEST) c(12L, 24L, 36L) else as.integer(seq(12L, 192L, by = 12L))
B_ACTIVE   <- 12L    # every 2h — Hughes 2017 >=2h recommendation
NSIMS           <- if (SMOKE_TEST) 5L   else 30L
NGENES          <- if (SMOKE_TEST) 200L else 5000L
NGENES_FMM_LRT  <- if (SMOKE_TEST) 100L else 1000L   # reduced G for FMM-LRT speed
FDR_THRESH <- 0.05
N_CORES      <- as.integer(Sys.getenv("MC_CORES",      unset = "8"))
GENE_CORES   <- as.integer(Sys.getenv("GENE_CORES",   unset = "4"))  # for FMM-LRT gene-level

cat(sprintf("Mode     : %s\n", if (SMOKE_TEST) "SMOKE" else "PRODUCTION"))
cat(sprintf("omega    : %s\n", paste(OMEGA_VALS, collapse = ", ")))
cat(sprintf("N grid   : %s\n", paste(N_GRID, collapse = ", ")))
cat(sprintf("B_ACTIVE : %d (every %.0fh)\n", B_ACTIVE, 24/B_ACTIVE))
cat(sprintf("nsims    : %d\n", NSIMS))
cat(sprintf("ngenes   : %d\n", NGENES))
cat(sprintf("mc.cores : %d\n\n", N_CORES))

out_dir <- "output/bvsm_method_comparison"
dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "results"),  recursive = TRUE, showWarnings = FALSE)

# =====================================================================
# 1. Load datasets
# =====================================================================
cat("--- Loading datasets ---\n")

## ---- ACTIVE ROW --------------------------------------------------- ##

## A: Mouse LIV (GSE54651)
dat_mouse <- readRDS("data/mice_GSE54651_CPM.RData")
prep_liv  <- prepCircadianData(dat_mouse$count_clean[["LIV"]],
                               times = dat_mouse$tod[["LIV"]], input_type = "log2")
mat_liv  <- prep_liv$data[rowSums(prep_liv$data > 0) >= 4, , drop = FALSE]
tod_liv  <- prep_liv$times
rm(dat_mouse, prep_liv)

## B: Baboon LUN (CAMO)
load("data/CAMO_PRC_hmb.RData")
prep_lun <- prepCircadianData(baboon_withTOD$baboon[["LUN"]],
                               times = baboon_withTOD$tod[["LUN"]], input_type = "cpm")
mat_lun  <- prep_lun$data[rowSums(prep_lun$data > 0) >= 6, , drop = FALSE]
tod_lun  <- prep_lun$times
rm(baboon_withTOD, gtex, mice, prep_lun)

## C: Mouse D1 (D1D2 striatum)
pheno   <- read.csv("data/mouse_clinicalinfo_03082021_rmOutliers.csv", row.names = 1)
prep_d1 <- prepCircadianData("data/mouse_D1D2_logCPMfiltered_counts.csv",
                              times = "time", input_type = "counts",
                              pheno = pheno, sample_col = "sample")
d1_samp <- pheno$sample[pheno$cell == "D1"]
mat_d1  <- prep_d1$data[, colnames(prep_d1$data) %in% d1_samp, drop = FALSE]
tod_d1  <- pheno$time[match(colnames(mat_d1), pheno$sample)]
mat_d1  <- mat_d1[rowSums(mat_d1 > 1) >= floor(ncol(mat_d1)/2), , drop = FALSE]
rm(pheno, prep_d1)

## ---- PASSIVE ROW --------------------------------------------------- ##

## Helper: extract GTEx tissue (CPM format) and parse TOD from column names
extract_gtex <- function(tissue_name, cpm_list) {
  df   <- cpm_list[[tissue_name]]
  ids  <- as.character(colnames(df))
  hhmm <- sapply(strsplit(ids, "\\."), function(x) if (length(x) >= 3) x[3] else NA)
  hrs  <- as.numeric(substr(hhmm, 1, 2)) + as.numeric(substr(hhmm, 3, 4)) / 60
  ok   <- !is.na(hrs)
  expr <- log2(as.matrix(df[, ok]) + 1)
  list(expr = expr, times = hrs[ok])
}

gtex_path <- "/home/qtp1/Projects/Collaborative/GTEXdata/TOD Xiangning/data/CPM/CPM.all.norm.RData"
load(gtex_path)

## D: GTEx Adrenal Gland (passive, strong)
adr_raw <- extract_gtex("Adrenal Gland", CPM.all.norm)
mat_adr <- adr_raw$expr[rowSums(adr_raw$expr > 0) >= 10, , drop = FALSE]
tod_adr <- adr_raw$times

rm(CPM.all.norm, adr_raw)

## F: NAc (Kyle GSE160521, passive, weak) — introduced in Fig 1
bio_nac_pilot <- readRDS("data/gse160521_nac_ctrl_pilot.rds")
stopifnot(!is.null(bio_nac_pilot$cts), length(bio_nac_pilot$cts) > 0)
tod_nac <- bio_nac_pilot$cts

## E: Putamen CTL (passive, moderate) — Kyle GSE160521, introduced in Fig 1 + Fig 2
## Load pre-estimated pilot directly (avoids reloading raw counts)
bio_put_pilot <- readRDS("data/gse160521_putamen_ctrl_pilot.rds")
tod_put <- bio_put_pilot$cts   # 59 sample TOD values
## Dummy mat/tod for estCircadianParam — we will use bio_put_pilot directly
mat_put <- NULL   # signal to use pre-estimated pilot below

# =====================================================================
# 2. Dataset registry
# =====================================================================
active_datasets <- list(
  LIV = list(mat = mat_liv, tod = tod_liv,
             label = "Mouse Liver\nr̃ ≈ 2.88",  snr = "Strong"),
  LUN = list(mat = mat_lun, tod = tod_lun,
             label = "Baboon Lung\nr̃ ≈ 1.72", snr = "Moderate"),
  D1  = list(mat = mat_d1,  tod = tod_d1,
             label = "Mouse D1 Striatum\nr̃ ≈ 0.65", snr = "Weak")
)
passive_datasets <- list(
  ADR  = list(mat = mat_adr, tod = tod_adr, bio_pre = NULL,
              label = "Human Adrenal Gland\nr̃ ≈ 1.03", snr = "Strong"),
  PUT  = list(mat = NULL,    tod = tod_put, bio_pre = bio_put_pilot,
              label = "Human Putamen\nr̃ ≈ 0.66",        snr = "Moderate"),
  NAC  = list(mat = NULL,    tod = tod_nac, bio_pre = bio_nac_pilot,
              label = "Human NAc\nr̃ ≈ 0.69",            snr = "Weak")
)
rm(mat_liv, mat_lun, mat_d1, tod_liv, tod_lun, tod_d1,
   mat_adr, tod_adr, bio_put_pilot, bio_nac_pilot)

# =====================================================================
# 3. Power simulation helper
# =====================================================================
run_fmm_power <- function(bio, N_vals, sweep_vals, sweep_param = "omega",
                           fixed_omega = 1.0, fixed_beta = pi, fixed_alpha = NULL,
                           B_val, design_type, pilot_tod = NULL,
                           method = "DCP", null_table = NULL,
                           nsims, ngenes, fdr_thresh, n_cores, gene_cores = 1L,
                           ds_name, seed = GLOBAL_SEED) {
  bio$ngenes <- ngenes
  cts_template <- if (design_type == "active")
    seq(0, 24 * (1 - 1/B_val), length.out = B_val) else NULL

  rows <- list()
  for (sv in sweep_vals) {
    omega_val <- if (sweep_param == "omega") sv else fixed_omega
    beta_val  <- if (sweep_param == "beta")  sv else fixed_beta
    alpha_val <- if (sweep_param == "alpha") sv else fixed_alpha  # NULL = draw from pilot
    cat(sprintf("    %s=%.3f  [%s]\n", sweep_param, sv, method))
    N_valid <- if (design_type == "active") N_vals[N_vals %% B_val == 0L] else N_vals

    run_one_N <- function(N_val) {
      sims <- vapply(seq_len(nsims), function(s) {
        cts <- if (design_type == "active") {
          rep(cts_template, each = N_val / B_val)
        } else {
          sampleTimesFromDist(N_val, pilot_tod)
        }
        dat <- simCircadianFMM(bio, cts,
                               omega       = omega_val,
                               beta        = beta_val,
                               alpha_fixed = alpha_val,
                               seed        = seed + N_val * 1000L + s)
        pv <- if (method == "FMM_LRT") {
          detect_FMM_LRT(dat$expr, cts, null_table = null_table,
                         mc.cores = gene_cores)
        } else {
          detect_DCP(dat$expr, cts)
        }
        adj <- p.adjust(pv, method = "BH")
        if (!any(dat$is_rhythmic)) return(NA_real_)
        sum(adj[dat$is_rhythmic] <= fdr_thresh, na.rm = TRUE) /
          sum(dat$is_rhythmic)
      }, numeric(1))
      data.frame(N = N_val, sweep_val = sv, sweep_param = sweep_param,
                 omega = omega_val, beta = beta_val,
                 alpha = if (!is.null(alpha_val)) alpha_val else NA_real_,
                 method = method,
                 power    = mean(sims, na.rm = TRUE),
                 power_se = sd(sims,   na.rm = TRUE) / sqrt(sum(!is.na(sims))),
                 stringsAsFactors = FALSE)
    }

    set.seed(seed)
    res_list <- parallel::mclapply(N_valid, run_one_N, mc.cores = n_cores)
    rows <- c(rows, res_list)
  }
  df         <- do.call(rbind, rows)
  df$dataset <- ds_name
  df
}

# =====================================================================
# 4. Run simulations
# =====================================================================
BETA_VALS   <- if (SMOKE_TEST) c(0, pi/2, pi) else
               c(0, pi/4, pi/2, 3*pi/4, pi, 3*pi/2)
ALPHA_VALS  <- if (SMOKE_TEST) c(0, 8, 16) else
               c(0, 4, 8, 12, 16, 20)        # acrophase sweep (hours)
OMEGA_FIXED <- 0.5   # fixed omega for beta/alpha sweeps

# Detection methods to run
# DCP sweeps (omega, beta, alpha) are already done — only run FMM_LRT now
METHODS <- c("FMM_LRT")

# Load FMM-LRT null calibration table (built by build_FMM_LRT_null_table)
fmm_null_table <- NULL
null_tbl_path  <- "output/fmm_lrt_null_table.rds"
if (file.exists(null_tbl_path)) {
  fmm_null_table <- readRDS(null_tbl_path)
  cat(sprintf("Loaded FMM-LRT null table (%d n values)\n", length(fmm_null_table$n)))
} else {
  cat("WARNING: FMM-LRT null table not found — FMM_LRT will use chi-squared(4) fallback.\n")
  cat("Run build_FMM_LRT_null_table() first for calibrated p-values.\n")
}

all_results <- list()

# Helper: run one dataset for a given sweep + design + detection method
run_dataset <- function(ds_name, ds, bio, design_type, sweep_param,
                         sweep_vals, fixed_omega, fixed_beta, fixed_alpha = NULL,
                         method = "DCP", null_table = NULL, tag) {
  pilot_tod <- if (design_type == "passive") ds$tod else NULL
  df <- run_fmm_power(bio, N_GRID, sweep_vals,
                       sweep_param  = sweep_param,
                       fixed_omega  = fixed_omega,
                       fixed_beta   = fixed_beta,
                       fixed_alpha  = fixed_alpha,
                       B_val        = if (design_type=="active") B_ACTIVE else NULL,
                       design_type  = design_type,
                       pilot_tod    = pilot_tod,
                       method       = method,
                       null_table   = null_table,
                       nsims = NSIMS,
                       ngenes = if (method == "FMM_LRT") NGENES_FMM_LRT else NGENES,
                       fdr_thresh = FDR_THRESH, n_cores = N_CORES,
                       gene_cores = if (method == "FMM_LRT") GENE_CORES else 1L,
                       ds_name = ds_name)
  design_label <- if (design_type=="active") "Active (B=12, every 2h)" else "Passive"
  sweep_label  <- switch(sweep_param,
    omega = sprintf("%s — omega sweep (beta=pi)", design_label),
    beta  = sprintf("%s — beta sweep (omega=%.1f)", design_label, fixed_omega),
    alpha = sprintf("%s — alpha sweep (omega=%.1f)", design_label, fixed_omega))
  df$design      <- if (design_type=="active") "Active (B=12, every 2h)" else "Passive (KDE TOD)"
  df$design_row  <- sweep_label
  df$panel_label <- ds$label
  df$snr_col     <- ds$snr
  saveRDS(df, file.path(out_dir, "results",
          sprintf("results_%s_FMM_%s_%s.rds", ds_name, tag, method)))
  df
}

# ─── Helper: estimate pilot bio once per dataset ─────────────────────
get_bio <- function(ds) {
  if (!is.null(ds$bio_pre)) ds$bio_pre else
    estCircadianParam(ds$mat, times = ds$tod, period = 24, verbose = TRUE)
}

# ─── Run all 6 rows × 2 methods ──────────────────────────────────────
sweep_specs <- list(
  list(row_a = 1, row_p = 2, param = "omega",
       vals  = OMEGA_VALS, fw = 1.0,        fb = pi,       fa = NULL,
       tag   = "omega"),
  list(row_a = 3, row_p = 4, param = "beta",
       vals  = BETA_VALS,  fw = OMEGA_FIXED, fb = pi,       fa = NULL,
       tag   = "beta"),
  list(row_a = 5, row_p = 6, param = "alpha",
       vals  = ALPHA_VALS, fw = OMEGA_FIXED, fb = pi,       fa = 0,
       tag   = "alpha")
)

for (meth in METHODS) {
  cat(sprintf("\n========== METHOD: %s ==========\n", meth))
  nt <- if (meth == "FMM_LRT") fmm_null_table else NULL

  for (sp in sweep_specs) {
    cat(sprintf("\n--- %s sweep  (rows %d/%d) ---\n",
                toupper(sp$param), sp$row_a, sp$row_p))

    # Active datasets
    cat(sprintf("  == Active (rows %d) ==\n", sp$row_a))
    for (ds_name in names(active_datasets)) {
      ds  <- active_datasets[[ds_name]]
      cat(sprintf("    %s\n", gsub("\n.*","",ds$label)))
      bio <- get_bio(ds)
      key <- sprintf("r%d_%s_%s", sp$row_a, ds_name, meth)
      all_results[[key]] <- run_dataset(
        ds_name, ds, bio, "active", sp$param, sp$vals,
        fixed_omega = sp$fw, fixed_beta = sp$fb, fixed_alpha = sp$fa,
        method = meth, null_table = nt,
        tag = sprintf("active_%s", sp$tag))
    }

    # Passive datasets
    cat(sprintf("  == Passive (rows %d) ==\n", sp$row_p))
    for (ds_name in names(passive_datasets)) {
      ds  <- passive_datasets[[ds_name]]
      cat(sprintf("    %s\n", gsub("\n.*","",ds$label)))
      bio <- get_bio(ds)
      key <- sprintf("r%d_%s_%s", sp$row_p, ds_name, meth)
      all_results[[key]] <- run_dataset(
        ds_name, ds, bio, "passive", sp$param, sp$vals,
        fixed_omega = sp$fw, fixed_beta = sp$fb, fixed_alpha = sp$fa,
        method = meth, null_table = nt,
        tag = sprintf("passive_%s", sp$tag))
    }
  }
}

# =====================================================================
# 5. Combine and plot
# =====================================================================
library(ggplot2)

full_df <- do.call(rbind, all_results)
# Convert SE → SD for visible error bars (SD = SE * sqrt(nsims))
full_df$power_sd <- full_df$power_se * sqrt(NSIMS)
saveRDS(full_df, file.path(out_dir, "results", "results_FMM_all.rds"))

# ── Color schemes ──────────────────────────────────────────────────
# ω sweep: dark→light blue
omega_levels <- as.character(sort(OMEGA_VALS))
n_omega      <- length(omega_levels)
omega_colors <- setNames(
  colorRampPalette(c("#08306B", "#4292C6", "#C6DBEF"))(n_omega),
  omega_levels
)
omega_labels <- setNames(paste0("ω=", omega_levels), omega_levels)

# β sweep: dark→light red/orange
beta_levels  <- as.character(round(BETA_VALS, 4))
n_beta       <- length(beta_levels)
beta_colors  <- setNames(
  colorRampPalette(c("#7F0000", "#EF6548", "#FEE8C8"))(n_beta),
  beta_levels
)
beta_denom <- c("0", "π/4", "π/2", "3π/4", "π", "5π/4", "3π/2", "7π/4", "2π")
beta_labels <- setNames(
  beta_denom[seq_len(n_beta)],
  beta_levels
)
omega_labels <- setNames(paste0("omega=", omega_levels), omega_levels)

active_labels  <- sapply(active_datasets,  function(d) d$label)
passive_labels <- sapply(passive_datasets, function(d) d$label)

full_df$omega_fac   <- factor(as.character(full_df$omega), levels = omega_levels)
# design_row levels are re-derived internally by plotFMMViolation via grepl("^Active", ...)
# Do not re-factor here to avoid level mismatch with actual sweep_label strings
full_df$snr_col <- factor(full_df$snr_col, levels = c("Strong", "Moderate", "Weak"))

# Build column headers that include both SNR tier and dataset names
snr_col_labels <- c(
  Strong   = sprintf("Strong\nActive: %s\nPassive: %s",
    sub("\n.*","", active_datasets$LIV$label),
    sub("\n.*","", passive_datasets$ADR$label)),
  Moderate = sprintf("Moderate\nActive: %s\nPassive: %s",
    sub("\n.*","", active_datasets$LUN$label),
    sub("\n.*","", passive_datasets$PUT$label)),
  Weak     = sprintf("Weak\nActive: %s\nPassive: %s",
    sub("\n.*","", active_datasets$D1$label),
    sub("\n.*","", passive_datasets$NAC$label))
)

x_breaks <- N_GRID[N_GRID %% 48 == 0 | N_GRID <= 48]

theme_fig4 <- theme_bw(base_size = 11) + theme(
  strip.background  = element_rect(fill = "grey92"),
  strip.text.y      = element_text(size = 8),
  strip.text.x      = element_text(size = 10, face = "bold"),
  legend.position   = "bottom",
  panel.grid.minor  = element_blank(),
  axis.text.x       = element_text(angle = 45, hjust = 1, size = 8),
  plot.title        = element_text(hjust = 0.5, face = "bold", size = 11)
)

col_labeller <- as_labeller(snr_col_labels)

# ── Generate one figure per detection method ─────────────────────
for (meth in METHODS) {
  df_meth <- full_df[full_df$method == meth, ]
  fig_path <- file.path(out_dir, "figures",
                        sprintf("fig4_fmm_%s.pdf", tolower(meth)))
  plotFMMViolation(df_meth, nsims = NSIMS, omega_fixed = OMEGA_FIXED,
                   output_file = fig_path, width = 16, height = 21)
  # also save to main_figures
  main_path <- sprintf("output/main_figures/Fig4_FMM_%s.pdf", meth)
  file.copy(fig_path, main_path, overwrite = TRUE)
  cat(sprintf("Saved: %s\n", fig_path))
}
cat("\n=== Done ===\n")
