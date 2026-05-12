#' =======================================================================
#' fig5_active_vs_passive_v4.R — §2.5 active vs passive design comparison
#' =======================================================================
#'
#' 2x2 layout: rows = detector (DCP, FMM K=2); cols = design (active baboon,
#' passive Seney). Demonstrates that the K-harmonic LRT generalizes from
#' controlled (active) to uncontrolled (passive) sampling.
#'
#' STORY:
#'   Top-left  (active, DCP):     B-invariance under FMM truth (theory check)
#'   Top-right (passive, DCP):    DCP power on Seney ACC (cosinor-like truth)
#'   Bot-left  (active, K=2):     K=2 advantage over DCP at moderate-violation tissues
#'   Bot-right (passive, K=2):    K=2 ≈ DCP on cosinor-like truth (no df cost dominance)

setwd("/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")
old_wd <- setwd("code"); source("setup.R"); setwd(old_wd)
suppressPackageStartupMessages(library(readxl))

SMOKE_TEST  <- identical(Sys.getenv("SMOKE_TEST"), "true")
NGENES      <- if (SMOKE_TEST) 500L  else 2000L
NSIMS       <- if (SMOKE_TEST) 5L    else 25L
N_CORES     <- as.integer(Sys.getenv("MC_CORES", unset = "30"))
TOP_K_FMM   <- if (SMOKE_TEST) 50L   else 300L
RHYTHM_PVAL <- 0.01
GLOBAL_SEED <- 2025L
K_HARM      <- 2L

# Active design (KIM single tissue, B sweep)
ACTIVE_TISSUE <- "KIM"
N_GRID_ACTIVE <- if (SMOKE_TEST) c(24L, 48L) else c(24L, 48L, 72L, 96L, 120L, 144L)
B_GRID_ACTIVE <- if (SMOKE_TEST) c(6L, 12L)  else c(6L, 8L, 12L, 24L)
# All B values satisfy the K=2 identifiability requirement B >= 2K+1 = 5
# (Nyquist condition). Active designs with B < 5 should use K=1 (cosinor).

# Passive design (Seney ACC)
N_GRID_PASSIVE <- if (SMOKE_TEST) c(40L, 80L) else c(40L, 60L, 80L, 100L, 120L, 160L, 200L, 240L)

cat(sprintf("Mode      : %s\n", if (SMOKE_TEST) "SMOKE" else "PRODUCTION"))
cat(sprintf("NGENES    : %d\n", NGENES))
cat(sprintf("NSIMS     : %d\n", NSIMS))

out_dir_fig <- "output/active_vs_passive/figures"
out_dir_res <- "output/active_vs_passive/results"
dir.create(out_dir_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_res, recursive = TRUE, showWarnings = FALSE)

analysis <- CircadianAnalysisOptions(alpha = 0.05, p.adjust.method = "BH",
                                     fdr_thresholds = 0.05)

# ====================================================================
# 1. Active pilot: KIM (use cached if available)
# ====================================================================
cat("\n=== Step 1: active pilot — KIM ===\n")
load("data/CAMO_PRC_hmb.RData")
prep_kim <- prepCircadianData(baboon_withTOD$baboon[[ACTIVE_TISSUE]],
                               times = baboon_withTOD$tod[[ACTIVE_TISSUE]] %% 24,
                               input_type = "cpm")
mat_kim  <- prep_kim$data[rowSums(prep_kim$data > 0) >= 6, , drop = FALSE]
tod_kim  <- prep_kim$times
rm(baboon_withTOD, gtex, mice, prep_kim)

set.seed(GLOBAL_SEED)
g_kim <- sample(nrow(mat_kim), min(NGENES, nrow(mat_kim)))
mat_kim_sub <- mat_kim[g_kim, ]

active_bio_path <- file.path(out_dir_res,
                             sprintf("bio_baboon_%s_NG%d_K%d.rds",
                                     ACTIVE_TISSUE, NGENES, TOP_K_FMM))
if (file.exists(active_bio_path) && !identical(Sys.getenv("REFIT"), "true")) {
  cat(sprintf("Loading cached active bio: %s\n", active_bio_path))
  bio_active <- readRDS(active_bio_path)
} else {
  bio_active <- estCircadianParamFMM(mat_kim_sub, tod_kim,
                                     min_rhythm_pval = RHYTHM_PVAL,
                                     top_k = TOP_K_FMM,
                                     mc.cores = N_CORES,
                                     verbose = TRUE)
  bio_active$ngenes <- NGENES
  saveRDS(bio_active, active_bio_path)
}

# ====================================================================
# 2. Passive pilot: Seney ACC ctrl
# ====================================================================
cat("\n=== Step 2: passive pilot — Seney ACC ctrl ===\n")
meta_s <- read_excel("data/MD5_MetaData_1-15-25.xlsx")
tod_s  <- read_excel("data/TOD.xlsx")
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

set.seed(GLOBAL_SEED)
g_sen <- sample(nrow(mat_sen), min(NGENES, nrow(mat_sen)))
mat_sen_sub <- mat_sen[g_sen, ]

passive_bio_path <- file.path(out_dir_res,
                              sprintf("bio_seney_passive_NG%d.rds", NGENES))
if (file.exists(passive_bio_path) && !identical(Sys.getenv("REFIT"), "true")) {
  cat(sprintf("Loading cached passive bio: %s\n", passive_bio_path))
  bio_passive <- readRDS(passive_bio_path)
} else {
  bio_passive <- estCircadianParam(mat_sen_sub, tod_sen,
                                   min_rhythm_pval = RHYTHM_PVAL,
                                   verbose = TRUE)
  bio_passive$ngenes <- NGENES
  saveRDS(bio_passive, passive_bio_path)
}
cat(sprintf("Passive bio: prop_rhythmic=%.3f  median_r=%.3f\n",
            bio_passive$prop_rhythmic,
            median(bio_passive$amplitude / bio_passive$sigma_rhythmic, na.rm = TRUE)))

# ====================================================================
# 3. Active row: B-vs-N for both detectors on KIM
# ====================================================================
cat("\n=== Step 3: active simulations (KIM) ===\n")
results_active <- list(DCP = list(), FMM = list())
for (mname in c("DCP", "FMM")) {
  for (B in B_GRID_ACTIVE) {
    N_valid <- N_GRID_ACTIVE[N_GRID_ACTIVE %% B == 0L]
    pwr <- vapply(N_valid, function(N) {
      m <- N %/% B
      cts_design <- rep(seq(0, 24*(1 - 1/B), length.out = B), each = m)
      design_i <- CircadianDesignOptions(sample_sizes = N, nsims = NSIMS,
                                         design = "active", cts = cts_design,
                                         B_values = B)
      r <- runSimsSingleCohort(bio_active, design_i, analysis,
                               method = mname, K = K_HARM,
                               mc.cores = N_CORES, verbose = FALSE)
      mean(r$marginal_power, na.rm = TRUE)
    }, numeric(1))
    results_active[[mname]][[as.character(B)]] <-
      data.frame(N = N_valid, power = pwr, B = B, method = mname)
    cat(sprintf("  active %s B=%-2d : %s\n", mname, B,
                paste(sprintf("N=%d:%.2f", N_valid, pwr), collapse = "  ")))
  }
}

# ====================================================================
# 4. Passive row: power-vs-N for both detectors on Seney
# ====================================================================
cat("\n=== Step 4: passive simulations (Seney) ===\n")
results_passive <- list(DCP = NULL, FMM = NULL)
for (mname in c("DCP", "FMM")) {
  pwr <- vapply(N_GRID_PASSIVE, function(N) {
    design_i <- CircadianDesignOptions(sample_sizes = N, nsims = NSIMS,
                                       design = "passive", cts = tod_sen)
    r <- runSimsSingleCohort(bio_passive, design_i, analysis,
                             method = mname, K = K_HARM,
                             mc.cores = N_CORES, verbose = FALSE)
    mean(r$marginal_power, na.rm = TRUE)
  }, numeric(1))
  results_passive[[mname]] <- data.frame(N = N_GRID_PASSIVE, power = pwr,
                                         method = mname)
  cat(sprintf("  passive %s : %s\n", mname,
              paste(sprintf("N=%d:%.2f", N_GRID_PASSIVE, pwr), collapse = "  ")))
}

# ====================================================================
# 5. Save
# ====================================================================
out_data <- list(
  bio_active     = bio_active,
  bio_passive    = bio_passive,
  results_active = results_active,
  results_passive = results_passive,
  N_GRID_ACTIVE  = N_GRID_ACTIVE,
  B_GRID_ACTIVE  = B_GRID_ACTIVE,
  N_GRID_PASSIVE = N_GRID_PASSIVE,
  K_HARM         = K_HARM
)
ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
saveRDS(out_data, file.path(out_dir_res, sprintf("fig5_v4_active_passive_%s.rds", ts)))

# ====================================================================
# 6. Plot — 2x2
# ====================================================================
cat("\n=== Step 5: render Fig 5 v4 (2x2) ===\n")
source("examples/publication/_pub_style.R")

fig_paths <- c(
  file.path(out_dir_fig, "fig5_v4_active_vs_passive.pdf"),
  "output/main_figures/Fig5_active_vs_passive.pdf"
)
dir.create("output/main_figures", recursive = TRUE, showWarnings = FALSE)

pal_B  <- pub_palette_sequential(length(B_GRID_ACTIVE))
pal_DT <- pub_palette_detector()
col_DCP <- unname(pal_DT["DCP"])
col_FMM <- unname(pal_DT["FMM"])

for (out_pdf in fig_paths) {
  cairo_pdf(out_pdf, width = 7.2, height = 6.4)
  pub_par(mfrow = c(2, 2), mar = c(4.0, 4.2, 2.4, 1.0))

  # --- Panel a: Active, DCP ---
  plot(NA, xlim = range(N_GRID_ACTIVE), ylim = c(0, 1),
       xlab = "N (total samples)", ylab = "Power",
       main = "Active, DCP")
  panel_label("a")
  abline_80pct()
  for (k in seq_along(B_GRID_ACTIVE)) {
    df <- results_active$DCP[[as.character(B_GRID_ACTIVE[k])]]
    lines(df$N, df$power, col = pal_B[k], lwd = 1.8)
    points(df$N, df$power, pch = 19, col = pal_B[k], cex = 0.7)
  }
  pub_legend("bottomright", legend = sprintf("%d", B_GRID_ACTIVE),
             col = pal_B, lwd = 1.6, title = "B (timepoints)")

  # --- Panel b: Passive, DCP ---
  pasD <- results_passive$DCP
  plot(NA, xlim = range(N_GRID_PASSIVE), ylim = c(0, 1),
       xlab = "N (total samples)", ylab = "Power",
       main = "Passive, DCP")
  panel_label("b")
  abline_80pct()
  lines(pasD$N, pasD$power, col = col_DCP, lwd = 2.0)
  points(pasD$N, pasD$power, pch = 19, col = col_DCP, cex = 0.7)

  # --- Panel c: Active, FMM K=2 ---
  plot(NA, xlim = range(N_GRID_ACTIVE), ylim = c(0, 1),
       xlab = "N (total samples)", ylab = "Power",
       main = sprintf("Active, K = %d", K_HARM))
  panel_label("c")
  abline_80pct()
  for (k in seq_along(B_GRID_ACTIVE)) {
    df <- results_active$FMM[[as.character(B_GRID_ACTIVE[k])]]
    lines(df$N, df$power, col = pal_B[k], lwd = 1.8)
    points(df$N, df$power, pch = 19, col = pal_B[k], cex = 0.7)
  }

  # --- Panel d: Passive, FMM K=2 with DCP overlay ---
  pasF <- results_passive$FMM
  plot(NA, xlim = range(N_GRID_PASSIVE), ylim = c(0, 1),
       xlab = "N (total samples)", ylab = "Power",
       main = sprintf("Passive, K = %d", K_HARM))
  panel_label("d")
  abline_80pct()
  lines(pasF$N, pasF$power, col = col_FMM, lwd = 2.0)
  points(pasF$N, pasF$power, pch = 19, col = col_FMM, cex = 0.7)
  lines(pasD$N, pasD$power, col = col_DCP, lwd = 1.4, lty = 2)
  pub_legend("bottomright",
             legend = c(sprintf("K = %d", K_HARM), "DCP (overlay)"),
             col = c(col_FMM, col_DCP), lwd = c(2.0, 1.4),
             lty = c(1, 2))

  dev.off()
  cat(sprintf("Saved: %s\n", out_pdf))
}

cat("\n=== Fig 5 v4 done ===\n")
