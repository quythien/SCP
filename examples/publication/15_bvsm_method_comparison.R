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
#' Outputs:
#'   Figure 1 — method comparison under cosinor truth (alpha2=0)
#'   Figure 2 — cosinor violation, DCP vs MH, D1 only
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
omega       <- 2 * pi / 24

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
# 1. Multi-harmonic helpers (adaptive K)
# =====================================================================
max_K <- function(B) floor((B - 1L) / 2L)

make_X <- function(times, K) {
  cols <- lapply(seq_len(K), function(k)
    cbind(cos(k * omega * times), sin(k * omega * times)))
  cbind(1, do.call(cbind, cols))
}

# =====================================================================
# 2. Load pilot datasets — single group from each pair
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
# 3. Estimate pilot parameters per dataset
# =====================================================================
cat("\n--- Estimating pilot parameters ---\n")

pilot_params <- lapply(names(datasets), function(ds_name) {
  ds    <- datasets[[ds_name]]
  fits  <- fitCosinorAll(ds$mat, times = ds$tod, period = 24)

  rhy_idx  <- !is.na(fits$pvalue) & fits$pvalue < 0.05
  A_all    <- as.numeric(fits$A)
  sig_all  <- as.numeric(fits$sigma)

  rhy_pool <- which(rhy_idx  & is.finite(A_all)   & A_all > 0 & is.finite(sig_all))
  arr_pool <- which(!rhy_idx & is.finite(sig_all) & sig_all > 0)

  # Sample following pilot proportions
  pi_R  <- mean(rhy_idx, na.rm = TRUE)
  n_rhy <- max(10L, round(NGENES * pi_R))
  n_arr <- NGENES - n_rhy

  set.seed(7L)
  g_rhy <- sample(rhy_pool, min(n_rhy, length(rhy_pool)))
  g_arr <- sample(arr_pool, min(n_arr, length(arr_pool)))

  r_med <- median(A_all[g_rhy] / sig_all[g_rhy], na.rm = TRUE)
  cat(sprintf("  %-5s  pi_R=%.0f%%  rhy=%d  arr=%d  median_r=%.2f\n",
              ds_name, 100*pi_R, length(g_rhy), length(g_arr), r_med))

  list(
    A_rhy   = A_all[g_rhy],
    sig_rhy = sig_all[g_rhy],
    sig_arr = sig_all[g_arr],
    n_rhy   = length(g_rhy),
    n_arr   = length(g_arr),
    pi_R    = pi_R,
    r_med   = r_med
  )
})
names(pilot_params) <- names(datasets)

# =====================================================================
# 4. Power cell functions
# =====================================================================

## 4a. DCP (K=1 cosinor F-test) --------------------------------------
power_cell_dcp <- function(A_rhy, sig_rhy, sig_arr, B, N, alpha2) {
  if (N %% B != 0) return(NA_real_)
  m      <- N / B
  cts    <- seq(0, 24 * (1 - 1/B), length.out = B)
  cts_n  <- rep(cts, each = m)

  X   <- cbind(1, cos(omega * cts_n), sin(omega * cts_n))
  XtX <- solve(t(X) %*% X)
  P   <- X %*% XtX %*% t(X)
  P0  <- matrix(1/N, N, N)
  df1 <- 2L; df2 <- N - 3L

  ngr <- length(A_rhy); nga <- length(sig_arr); ng <- ngr + nga

  sim_power <- replicate(NSIMS, {
    pvals <- numeric(ng)
    for (g in seq_len(ngr)) {
      phi  <- runif(1, 0, 2*pi)
      mu   <- A_rhy[g] * (cos(omega*cts_n - phi) +
                           alpha2 * cos(2*omega*cts_n - phi))
      y    <- rnorm(N, mu, sig_rhy[g])
      SSR  <- sum((P%*%y  - P0%*%y)^2)
      SSE  <- sum(y^2) - sum((P%*%y)^2)
      pvals[g] <- pf((SSR/df1)/(SSE/df2), df1, df2, lower.tail = FALSE)
    }
    for (g in seq_len(nga)) {
      y    <- rnorm(N, 0, sig_arr[g])
      SSR  <- sum((P%*%y  - P0%*%y)^2)
      SSE  <- sum(y^2) - sum((P%*%y)^2)
      pvals[ngr+g] <- pf((SSR/df1)/(SSE/df2), df1, df2, lower.tail = FALSE)
    }
    adj <- p.adjust(pvals, method = "BH")
    sum(adj[seq_len(ngr)] <= FDR_THRESH) / ngr
  })
  mean(sim_power, na.rm = TRUE)
}

## 4b. JTK_CYCLE ------------------------------------------------------
power_cell_jtk <- function(A_rhy, sig_rhy, sig_arr, B, N, alpha2) {
  if (N %% B != 0) return(NA_real_)
  m    <- N / B
  cts  <- seq(0, 24 * (1 - 1/B), length.out = B)
  cts_n <- rep(cts, each = m)
  ngr  <- length(A_rhy); nga <- length(sig_arr); ng <- ngr + nga
  gnames <- c(paste0("R", seq_len(ngr)), paste0("A", seq_len(nga)))

  sim_power <- replicate(NSIMS, {
    expr <- matrix(NA_real_, ng, N)
    for (g in seq_len(ngr)) {
      phi  <- runif(1, 0, 2*pi)
      mu   <- A_rhy[g] * (cos(omega*cts_n - phi) +
                           alpha2 * cos(2*omega*cts_n - phi))
      expr[g, ] <- rnorm(N, mu, sig_rhy[g])
    }
    for (g in seq_len(nga))
      expr[ngr+g, ] <- rnorm(N, 0, sig_arr[g])
    rownames(expr) <- gnames

    pvals <- tryCatch(
      detect_JTK(expr, cts_n, gene_names = gnames, period = 24),
      error = function(e) rep(1, ng)
    )
    adj <- p.adjust(pvals, method = "BH")
    sum(adj[seq_len(ngr)] <= FDR_THRESH, na.rm = TRUE) / ngr
  })
  mean(sim_power, na.rm = TRUE)
}

## 4c. Multi-harmonic (adaptive K) ------------------------------------
power_cell_mh <- function(A_rhy, sig_rhy, sig_arr, B, N, alpha2) {
  if (N %% B != 0) return(NA_real_)
  m    <- N / B
  cts  <- seq(0, 24 * (1 - 1/B), length.out = B)
  cts_n <- rep(cts, each = m)
  K    <- max_K(B)
  Xf   <- make_X(cts_n, K)
  df1  <- 2L * K
  df2  <- N - ncol(Xf)
  XtXi <- solve(t(Xf) %*% Xf)
  Pf   <- Xf %*% XtXi %*% t(Xf)
  P0   <- matrix(1/N, N, N)
  ngr  <- length(A_rhy); nga <- length(sig_arr); ng <- ngr + nga

  sim_power <- replicate(NSIMS, {
    pvals <- numeric(ng)
    for (g in seq_len(ngr)) {
      phi  <- runif(1, 0, 2*pi)
      mu   <- A_rhy[g] * (cos(omega*cts_n - phi) +
                           alpha2 * cos(2*omega*cts_n - phi))
      y    <- rnorm(N, mu, sig_rhy[g])
      SSR  <- sum((Pf%*%y - P0%*%y)^2)
      SSE  <- sum(y^2) - sum((Pf%*%y)^2)
      pvals[g] <- pf((SSR/df1)/(SSE/df2), df1, df2, lower.tail = FALSE)
    }
    for (g in seq_len(nga)) {
      y    <- rnorm(N, 0, sig_arr[g])
      SSR  <- sum((Pf%*%y - P0%*%y)^2)
      SSE  <- sum(y^2) - sum((Pf%*%y)^2)
      pvals[ngr+g] <- pf((SSR/df1)/(SSE/df2), df1, df2, lower.tail = FALSE)
    }
    adj <- p.adjust(pvals, method = "BH")
    sum(adj[seq_len(ngr)] <= FDR_THRESH) / ngr
  })
  mean(sim_power, na.rm = TRUE)
}

# =====================================================================
# 5. Run grid
# =====================================================================
methods <- list(
  DCP = power_cell_dcp,
  JTK = power_cell_jtk,
  MH  = power_cell_mh
)

all_results <- list()

for (ds_name in names(datasets)) {
  pp <- pilot_params[[ds_name]]
  cat(sprintf("\n==============================\n"))
  cat(sprintf("Dataset: %s\n", datasets[[ds_name]]$label))
  cat(sprintf("==============================\n"))

  ds_results <- list()

  for (meth_name in names(methods)) {
    meth_fn <- methods[[meth_name]]
    cat(sprintf("\n  Method: %s\n", meth_name))

    # Determine alpha2 values: JTK only runs alpha2=0 (direction is alpha2-invariant)
    a2_vals <- if (meth_name == "JTK") 0 else ALPHA2_VALS

    for (alpha2 in a2_vals) {
      cat(sprintf("    alpha2=%.2f ...", alpha2))
      t0 <- proc.time()[["elapsed"]]

      # Build grid of (N, B) cells
      grid <- expand.grid(N = N_GRID, B = B_VALS)
      grid <- grid[grid$N %% grid$B == 0, ]  # keep valid cells only

      pwr_list <- parallel::mclapply(
        seq_len(nrow(grid)),
        function(i) {
          set.seed(GLOBAL_SEED + i)
          meth_fn(pp$A_rhy, pp$sig_rhy, pp$sig_arr,
                  B      = grid$B[i],
                  N      = grid$N[i],
                  alpha2 = alpha2)
        },
        mc.cores = N_CORES
      )

      grid$power <- unlist(pwr_list)
      grid$alpha2 <- alpha2
      grid$method <- meth_name
      grid$dataset <- ds_name

      key <- sprintf("%s_%s_a%.2f", ds_name, meth_name, alpha2)
      ds_results[[key]] <- grid

      elapsed <- proc.time()[["elapsed"]] - t0
      cat(sprintf(" done (%.0fs)\n", elapsed))
    }
  }

  all_results[[ds_name]] <- ds_results
  saveRDS(ds_results,
          file.path(out_dir, "results", sprintf("results_%s.rds", ds_name)))
  cat(sprintf("\nSaved: results_%s.rds\n", ds_name))
}

# Combined results table
res_df <- do.call(rbind, unlist(all_results, recursive = FALSE))
saveRDS(res_df, file.path(out_dir, "results", "results_all.rds"))

# =====================================================================
# 6. Figures
# =====================================================================
library(ggplot2)
library(dplyr)

b_colors <- c("3"="#E41A1C","4"="#FF7F00","6"="#4DAF4A","8"="#377EB8","12"="#984EA3")
b_labels <- setNames(paste0("B=", B_VALS), as.character(B_VALS))

## Figure 1: Method comparison under cosinor truth (alpha2=0) ---------
cat("\nGenerating Figure 1...\n")

fig1_dat <- res_df %>%
  filter(alpha2 == 0) %>%
  mutate(
    method  = factor(method,  levels = c("DCP","JTK","MH"),
                     labels = c("DCP (K=1 cosinor)", "JTK_CYCLE", "Multi-harmonic (adaptive K)")),
    dataset = factor(dataset, levels = c("LIV","LUN","D1"),
                     labels = c("Mouse LIV\n(r~2.88, strong)",
                                "Baboon LUN\n(r~1.72, moderate)",
                                "Mouse D1\n(r~0.65, weak)")),
    B = factor(B, levels = B_VALS)
  )

p1 <- ggplot(fig1_dat, aes(x = N, y = 100*power, colour = B, group = B)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  geom_hline(yintercept = 80, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  facet_grid(method ~ dataset) +
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

fig1_path <- file.path(out_dir, "figures", "fig1_method_comparison.pdf")
ggsave(fig1_path, p1, width = 11, height = 9)
cat(sprintf("Saved: %s\n", fig1_path))

## Figure 2: Cosinor violation — DCP vs MH, D1 only ------------------
cat("Generating Figure 2...\n")

fig2_dat <- res_df %>%
  filter(dataset == "D1", method %in% c("DCP","MH"), N %in% c(48L, 96L)) %>%
  mutate(
    method  = factor(method, levels = c("DCP","MH"),
                     labels = c("DCP (K=1 cosinor)", "Multi-harmonic (adaptive K)")),
    N_label = paste0("N=", N),
    B       = factor(B, levels = B_VALS)
  )

p2 <- ggplot(fig2_dat, aes(x = alpha2, y = 100*power, colour = B, group = B)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_hline(yintercept = 80, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  facet_grid(N_label ~ method) +
  scale_colour_manual(values = b_colors, labels = b_labels, name = NULL) +
  scale_x_continuous(breaks = ALPHA2_VALS) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  labs(
    x     = expression(alpha[2] ~ "(second harmonic ratio A"[2]*"/A"[1]*")"),
    y     = "Detection power (%)",
    title = "Cosinor violation: DCP vs multi-harmonic (Mouse D1, r~0.65)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey92"),
    legend.position  = "bottom",
    panel.grid.minor = element_blank()
  )

fig2_path <- file.path(out_dir, "figures", "fig2_violation.pdf")
ggsave(fig2_path, p2, width = 8, height = 7)
cat(sprintf("Saved: %s\n", fig2_path))

cat("\n=== Done ===\n")
