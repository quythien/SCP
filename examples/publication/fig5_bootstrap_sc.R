#' =======================================================================
#' fig5_bootstrap_sc.R — Single-Cohort Bootstrap: Plug-in vs Bootstrap CI
#' =======================================================================
#'
#' Three passive-design panels showing how bootstrap CI width depends on
#' pilot size and where the planned study sits on the power curve. All three
#' panels use paired_sigma=TRUE so plug-in and bootstrap operate on a
#' consistent (A, σ) sampling scheme.
#'
#'   A. Seney ACC ctrl       n_pilot=60   small  pilot, post-mortem brain
#'      → CI peaks ~10pp at N=80 (elbow), collapses by N=120
#'
#'   B. GTEx Pancreas        n_pilot=249  moderate pilot
#'      → CI ~6-10pp across N=100-200; needs N~200 for 80% power
#'
#'   C. GTEx Thyroid         n_pilot=416  large  pilot
#'      → CI persists ~6-7pp across N=80-160 (long elbow); needs N~200 for 80%
#'
#' Story: bootstrap quantifies genuine uncertainty wherever the power curve
#' is in the elbow region (30-75%) — even with large pilots. Plug-in and
#' bootstrap means agree (gap <2pp); the bootstrap value is the CI width.
#'
#' USAGE:
#'   Rscript examples/publication/fig5_bootstrap_sc.R
#'   SMOKE_TEST=true Rscript examples/publication/fig5_bootstrap_sc.R

SMOKE_TEST  <- identical(Sys.getenv("SMOKE_TEST"), "true")
NGENES      <- if (SMOKE_TEST) 300L  else 5000L
NBOOT       <- if (SMOKE_TEST) 5L    else 50L
NSIMS       <- if (SMOKE_TEST) 5L    else 30L
NSIMS_INNER <- if (SMOKE_TEST) 5L    else 25L
N_CORES     <- as.integer(Sys.getenv("MC_CORES", unset = "6"))
RHYTHM_PVAL <- 0.01
B_VAL       <- 1L                      # passive: B is not a power-relevant parameter; placeholder for API compatibility
GLOBAL_SEED <- 2025L

cat(sprintf("Mode        : %s\n", if (SMOKE_TEST) "SMOKE" else "PRODUCTION"))
cat(sprintf("NGENES      : %d\n", NGENES))
cat(sprintf("NBOOT       : %d\n", NBOOT))
cat(sprintf("NSIMS       : %d\n", NSIMS))
cat(sprintf("NSIMS_INNER : %d\n", NSIMS_INNER))
cat(sprintf("MC_CORES    : %d\n\n", N_CORES))

old_wd <- setwd("code"); source("setup.R"); setwd(old_wd)
suppressPackageStartupMessages(library(readxl))

out_dir <- "output/bootstrap_sc"
dir.create(file.path(out_dir, "results"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)

analysis <- CircadianAnalysisOptions(
  alpha = 0.05, p.adjust.method = "BH", fdr_thresholds = 0.05
)

# -----------------------------------------------------------------------
# Helper: plug-in + bootstrap for one PASSIVE dataset
# -----------------------------------------------------------------------
.run_passive_comparison <- function(mat, tod, N_grid, label, out_prefix) {
  cat(sprintf("\n--- %s ---\n", label))

  # Estimate pilot fresh with paired_sigma=TRUE so plug-in matches bootstrap
  # subsample only for cosinor fit speed; bootstrap gets the FULL gene matrix
  # so prop_rhythmic and top-K pool match the plug-in estimate
  set.seed(GLOBAL_SEED)
  if (nrow(mat) > NGENES) {
    g_idx <- sample(nrow(mat), NGENES)
    mat_s <- mat[g_idx, , drop = FALSE]
  } else {
    mat_s <- mat
  }

  bio_sc <- estCircadianParam(mat_s, tod,
                               min_rhythm_pval = RHYTHM_PVAL,
                               paired_sigma    = TRUE,
                               verbose         = FALSE)
  bio_sc$ngenes <- NGENES
  cat(sprintf("    n_pilot=%d  prop_rhy=%.1f%%  n_rhythmic=%d\n",
              ncol(mat), 100 * bio_sc$prop_rhythmic, length(bio_sc$amplitude)))

  # Plug-in (passive)
  design_plug <- CircadianDesignOptions(
    sample_sizes = N_grid, nsims = NSIMS, design = "passive", cts = tod
  )
  cat("  Running plug-in (passive)...\n")
  set.seed(GLOBAL_SEED)
  plugin_res <- tryCatch(
    runSimsSingleCohort(bio_sc, design_plug, analysis,
                        mc.cores = N_CORES, verbose = FALSE),
    error = function(e) { warning(sprintf("Plug-in failed: %s", e$message)); NULL }
  )

  # Bootstrap (passive) — B_VAL is a placeholder; passive power is N-only
  boot_opts <- CircadianBootstrapOptions(
    design_vector = tod, B_values = B_VAL, N_values = N_grid,
    nboot = NBOOT, nsims_inner = NSIMS_INNER,
    design = "passive", seed = GLOBAL_SEED
  )
  cat(sprintf("  Running bootstrap (%d draws, mc.cores=%d)...\n", NBOOT, N_CORES))
  set.seed(GLOBAL_SEED)
  bt_res <- tryCatch(
    runBootstrapDesignGrid(
      pilot_data      = mat_s,
      pilot_times     = tod,
      bio_diff.opts   = bio_sc,
      boot.opts       = boot_opts,
      analysis.opts   = analysis,
      mode            = "single",
      min_rhythm_pval = RHYTHM_PVAL,
      mc.cores        = N_CORES,
      verbose         = FALSE
    ),
    error = function(e) { warning(sprintf("Bootstrap failed: %s", e$message)); NULL }
  )

  if (is.null(plugin_res) || is.null(bt_res)) {
    cat("  FAILED — skipping\n"); return(NULL)
  }

  plugin_mean <- apply(plugin_res$marginal_power, 1, mean, na.rm = TRUE)
  cat(sprintf("  %-5s  plug-in  bt-mean  CI-width\n", "N"))
  for (i in seq_along(N_grid)) {
    tp  <- plugin_mean[i]
    bm  <- bt_res$power_mean[i, 1, 1]
    lo  <- bt_res$power_ci_lo[i, 1, 1]
    hi  <- bt_res$power_ci_hi[i, 1, 1]
    cat(sprintf("  N=%-4d  %5.1f%%   %5.1f%%   %5.1fpp\n",
        N_grid[i], 100*tp, 100*bm, 100*(hi-lo)))
  }

  out <- list(label = label, n_pilot = ncol(mat), N_grid = N_grid,
              plugin = plugin_res, boot = bt_res)
  saveRDS(out, paste0(out_prefix, ".rds"))
  cat(sprintf("  Saved: %s.rds\n", basename(out_prefix)))
  out
}


# =======================================================================
# PANEL A: Putamen SCZ (n=28, passive) — Kyle GSE160521 multiBrain
# =======================================================================
cat("================================================================\n")
cat("PANEL A: Putamen SCZ (n=28, passive)\n")
cat("================================================================\n")

KYLE_DIR <- Sys.getenv("KYLE_MULTIBRAINREGION_DIR",
                        unset = "/home/qtp1/Projects/Circadian/Kyle/Kyle_multiBrainRegion")
N_GRID_SCZ <- if (SMOKE_TEST) c(20L, 40L, 60L) else
              c(20L, 40L, 60L, 80L, 100L, 120L, 140L, 160L)

clin_scz <- read.csv(file.path(KYLE_DIR, "DS_clinical_1221_rm97_rm231_matchIndex34.csv"),
                      row.names = 1)
scz_meta <- clin_scz[clin_scz$Diagnostic.Category == "SCZ", ]
expr_scz_raw <- as.matrix(read.csv(
  file.path(KYLE_DIR, "Putamen_CPMfiltered_logCPM_1215_rm97_rm231.csv"),
  row.names = 1, check.names = FALSE))
scz_cols <- intersect(scz_meta$pair, colnames(expr_scz_raw))
mat_scz  <- expr_scz_raw[, scz_cols, drop = FALSE]
tod_scz  <- scz_meta$CorrectedTOD[match(scz_cols, scz_meta$pair)] %% 24
prep_scz <- prepCircadianData(mat_scz, times = tod_scz, input_type = "log2")
mat_scz  <- prep_scz$data; tod_scz <- prep_scz$times
rm(clin_scz, expr_scz_raw, prep_scz)
cat(sprintf("Loaded Putamen-SCZ: %d genes x %d samples; TOD [%.2f, %.2f]\n",
            nrow(mat_scz), ncol(mat_scz), min(tod_scz), max(tod_scz)))

pSCZ <- .run_passive_comparison(mat_scz, tod_scz, N_GRID_SCZ,
                                 sprintf("Putamen SCZ (n=%d, passive)", ncol(mat_scz)),
                                 file.path(out_dir, "results", "panelA_putamen_scz"))


# =======================================================================
# PANEL B: Seney CTL ACC (n=60, passive)
# =======================================================================
cat("================================================================\n")
cat("PANEL B: Seney ACC ctrl (n=60, passive)\n")
cat("================================================================\n")

N_GRID_SEN <- if (SMOKE_TEST) c(40L, 80L, 120L) else
              c(40L, 60L, 80L, 100L, 120L, 140L, 160L, 200L, 240L, 280L)

meta_s   <- read_excel("data/MD5_MetaData_1-15-25.xlsx")
tod_s    <- read_excel("data/TOD.xlsx")
expr_raw <- as.matrix(read.csv("data/ACC_RNA_filtered_normalized.csv",
                                row.names = 1, check.names = FALSE))
col_ids  <- gsub("[A-Za-z]+$", "", colnames(expr_raw))
meta_idx <- match(col_ids, as.character(meta_s$HU_NUM))
tod_idx  <- match(col_ids, as.character(tod_s$HU_NUM))
tod_hour <- as.numeric(format(tod_s$OFFC_TIME[tod_idx], "%H")) +
            as.numeric(format(tod_s$OFFC_TIME[tod_idx], "%M")) / 60
disease  <- meta_s$Disease[meta_idx]
ok       <- !is.na(disease) & !is.na(tod_hour)
prep_sen <- prepCircadianData(expr_raw[, ok], times = tod_hour[ok], input_type = "log2")
ctrl_idx <- disease[ok] == 1
mat_sen  <- prep_sen$data[, ctrl_idx, drop = FALSE]
tod_sen  <- prep_sen$times[ctrl_idx]
rm(meta_s, tod_s, expr_raw, prep_sen)

pA <- .run_passive_comparison(mat_sen, tod_sen, N_GRID_SEN,
                               sprintf("Seney ACC ctrl (n=%d, passive)", ncol(mat_sen)),
                               file.path(out_dir, "results", "panelA_seney"))


# =======================================================================
# PANEL B: GTEx Pancreas (n=249, passive)
# =======================================================================
cat("\n================================================================\n")
cat("PANEL B: GTEx Pancreas (n=249, passive)\n")
cat("================================================================\n")

N_GRID_PAN <- if (SMOKE_TEST) c(40L, 80L, 120L) else
              c(40L, 80L, 120L, 160L, 200L, 240L, 280L, 320L)

GTEX_CPM_PATH <- Sys.getenv("GTEX_CPM_PATH",
                  unset = "/home/qtp1/Projects/Collaborative/GTEXdata/TOD Xiangning/data/CPM/CPM.all.norm.RData")
load(GTEX_CPM_PATH)
.extract_gtex <- function(tissue) {
  df   <- CPM.all.norm[[tissue]]
  ids  <- as.character(colnames(df))
  hhmm <- sapply(strsplit(ids, "\\."),
                 function(x) if (length(x) >= 3) x[3] else NA)
  hrs  <- suppressWarnings(as.numeric(substr(hhmm, 1, 2)) +
                           as.numeric(substr(hhmm, 3, 4)) / 60)
  ok   <- !is.na(hrs)
  list(mat = as.matrix(df[, ok]), tod = hrs[ok])
}
pan_d <- .extract_gtex("Pancreas")
mat_pan <- pan_d$mat; tod_pan <- pan_d$tod
rm(pan_d)

pB <- .run_passive_comparison(mat_pan, tod_pan, N_GRID_PAN,
                               sprintf("GTEx Pancreas (n=%d, passive)", ncol(mat_pan)),
                               file.path(out_dir, "results", "panelB_pancreas"))


# =======================================================================
# PANEL C: GTEx Thyroid (n=407, passive)
# =======================================================================
cat("\n================================================================\n")
cat("PANEL C: GTEx Thyroid (n=407, passive)\n")
cat("================================================================\n")

N_GRID_THY <- if (SMOKE_TEST) c(80L, 160L, 240L) else
              c(40L, 80L, 120L, 160L, 200L, 240L, 280L, 320L, 400L)

thy_d <- .extract_gtex("Thyroid")
mat_thy <- thy_d$mat; tod_thy <- thy_d$tod
rm(thy_d, CPM.all.norm)

pC <- .run_passive_comparison(mat_thy, tod_thy, N_GRID_THY,
                               sprintf("GTEx Thyroid (n=%d, passive)", ncol(mat_thy)),
                               file.path(out_dir, "results", "panelC_thyroid"))


# =======================================================================
# Combined figure (PDF only — no PNG)
# =======================================================================
cat("\n=== Generating combined figure ===\n")
panels <- Filter(Negate(is.null), list(pSCZ, pA, pB, pC))

if (length(panels) > 0) {
  fig_path_local <- file.path(out_dir, "figures", "fig_bootstrap_sc.pdf")
  fig_path_main  <- "output/main_figures/Fig3_bootstrap_singlecohort.pdf"
  dir.create("output/main_figures", recursive = TRUE, showWarnings = FALSE)

  for (out_pdf in c(fig_path_local, fig_path_main)) {
    n_panels <- length(panels)
    layout_dims <- if (n_panels >= 4) c(2L, 2L)
                    else if (n_panels == 3) c(1L, 3L)
                    else c(1L, n_panels)
    pdf_h <- if (layout_dims[1] == 2) 9 else 5
    cairo_pdf(out_pdf, width = 14, height = pdf_h)
    par(mfrow = layout_dims, mar = c(4.2, 4.2, 3, 1), las = 1,
        cex.lab = 1.1, cex.axis = 1.0, cex.main = 1.05, font.main = 2)

    panel_letters <- LETTERS[seq_along(panels)]
    for (pi in seq_along(panels)) {
      p      <- panels[[pi]]
      N_grid <- p$N_grid
      tp     <- apply(p$plugin$marginal_power, 1, mean, na.rm = TRUE)
      bt_mn  <- p$boot$power_mean[, 1, 1]
      bt_lo  <- p$boot$power_ci_lo[, 1, 1]
      bt_hi  <- p$boot$power_ci_hi[, 1, 1]

      plot(N_grid, tp, type = "l", lwd = 2.2, col = "steelblue",
           ylim = c(0, 1), xlab = "N (total samples)", ylab = "Power",
           main = sprintf("(%s) %s", panel_letters[pi], p$label))
      # Bootstrap mean trend line (thin) + vertical 95% CI error bars at each tested N
      lines(N_grid, bt_mn, lwd = 1, col = "tomato", lty = 2)
      arrows(N_grid, bt_lo, N_grid, bt_hi,
             code = 3, angle = 90, length = 0.05, lwd = 2, col = "tomato")
      points(N_grid, bt_mn, pch = 19, col = "tomato", cex = 0.9)
      abline(h = 0.80, lty = 2, col = "grey50")
      if (pi == 1) {
        legend("bottomright", bty = "n",
               legend = c("Plug-in (point estimate)",
                          "Bootstrap mean + 95% CI"),
               lwd    = c(2.2, 2),
               lty    = c(1, NA),
               pch    = c(NA, 19),
               col    = c("steelblue", "tomato"))
      }
    }
    dev.off()
    cat(sprintf("Saved: %s\n", out_pdf))
  }
}

saveRDS(panels, file.path(out_dir, "results", "all_panels.rds"))
cat("\n=== Done ===\n")
