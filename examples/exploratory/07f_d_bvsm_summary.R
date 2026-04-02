#' =======================================================================
#' 07f_d_bvsm_summary.R — B vs m Tradeoff: Cross-Dataset Summary Figure
#' =======================================================================
#' Loads RDS outputs from 07f_a/b/c and produces:
#'   1. s4_bvsm_grid.pdf  — 3-row (dataset) x 4-col (N) panels,
#'                           power vs alpha2, one line per B value.
#'                           Directly shows where B-invariance holds/breaks.
#'   2. s4_bvsm_alpha0.pdf — Power vs N at alpha2=0 across B (Prop 1 check).
#'   3. s4_bvsm_drop.pdf   — Power drop (alpha2=0.5 minus alpha2=0) vs N,
#'                            colored by B. Shows where distortion costs most.
#'
#' Run after 07f_a, 07f_b, 07f_c complete:
#'   Rscript examples/exploratory/07f_d_bvsm_summary.R

POWERSIM_ROOT <- {
  env <- Sys.getenv("POWERSIM_ROOT", unset = "")
  if (nzchar(env)) env else
    "/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim"
}
setwd(POWERSIM_ROOT)

RUN_TAG <- Sys.getenv("RUN_TAG", format(Sys.Date(), "%Y%m%d"))
out_dir <- file.path("output", paste0("07_bvsm_", RUN_TAG))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------
# Load datasets
# -----------------------------------------------------------------------
rds_paths <- list(
  mouse  = file.path(out_dir, "s1_mouse_bvsm.rds"),
  baboon = file.path(out_dir, "s2_baboon_bvsm.rds"),
  d1d2   = file.path(out_dir, "s3_d1d2_bvsm.rds")
)

missing <- names(rds_paths)[!sapply(rds_paths, file.exists)]
if (length(missing) > 0) {
  stop(sprintf(
    "Missing RDS files (jobs still running?): %s\n  Expected in: %s",
    paste(missing, collapse = ", "), out_dir
  ))
}

ds_list <- lapply(rds_paths, readRDS)
cat("Loaded datasets:\n")
for (nm in names(ds_list)) {
  d <- ds_list[[nm]]
  cat(sprintf("  %-8s  label=%s  pilot_n=%d  r=%.2f  B=%s  N=%s\n",
              nm, d$label, d$pilot_n,
              NA,   # r stored in label string
              paste(d$B_vals, collapse = "/"),
              paste(d$N_grid, collapse = "/")))
}

# Order: descending r (Mouse > Baboon > D1D2)
# Labels already encode r in d$label (e.g. "GSE54651 (r~2.9)")
ds_order <- c("mouse", "baboon", "d1d2")
ds_cols  <- c(mouse = "steelblue", baboon = "darkorange", d1d2 = "forestgreen")

# -----------------------------------------------------------------------
# Figure 1: Full grid — 3 rows (dataset) x 4 cols (N), power vs alpha2
# -----------------------------------------------------------------------
cat("\nGenerating s4_bvsm_grid.pdf ...\n")

# Determine grid dimensions from first dataset
n_N    <- length(ds_list[[1]]$N_grid)
n_rows <- length(ds_order)

pdf(file.path(out_dir, "s4_bvsm_grid.pdf"),
    width = 4 * n_N, height = 4.5 * n_rows)
par(mfrow = c(n_rows, n_N), mar = c(4.5, 4.5, 3.2, 1.2))

for (ds_nm in ds_order) {
  d      <- ds_list[[ds_nm]]
  a2_vals <- d$harm_grid$alpha2
  B_vals  <- d$B_vals
  N_grid  <- d$N_grid
  b_cols  <- colorRampPalette(c("gray70", ds_cols[ds_nm]))(length(B_vals))
  b_pchs  <- c(16, 17, 15, 18)[seq_along(B_vals)]

  for (ni in seq_along(N_grid)) {
    N <- N_grid[ni]
    plot(NA, xlim = range(a2_vals), ylim = c(0, 1), las = 1,
         xlab = expression(alpha[2]),
         ylab = "Power (DR, FDR 5%)",
         main = sprintf("%s\nN = %d", d$label, N))
    abline(h = 0.80, lty = 3, col = "gray70")
    abline(v = 0.50, lty = 3, col = "gray70")

    for (bi in seq_along(B_vals)) {
      pm <- d$power_mean[bi, ni, ]
      se <- d$power_se[bi,   ni, ]
      polygon(c(a2_vals, rev(a2_vals)),
              c(pm - se, rev(pm + se)),
              col = adjustcolor(b_cols[bi], 0.12), border = NA)
      lines(a2_vals,  pm, col = b_cols[bi], lwd = 2)
      points(a2_vals, pm, col = b_cols[bi], pch = b_pchs[bi], cex = 1.1)
    }

    if (ni == 1)
      legend("topright",
             legend = sprintf("B=%d (m=%d)", B_vals, N / B_vals),
             col = b_cols, lwd = 2, pch = b_pchs, bty = "n", cex = 0.75)
  }
}
dev.off()
cat(sprintf("  -> %s\n", file.path(out_dir, "s4_bvsm_grid.pdf")))

# -----------------------------------------------------------------------
# Figure 2: alpha2=0 — Power vs N for all B values (Proposition 1 check)
# -----------------------------------------------------------------------
cat("Generating s4_bvsm_alpha0.pdf ...\n")

pdf(file.path(out_dir, "s4_bvsm_alpha0.pdf"),
    width = 5 * n_rows, height = 5)
par(mfrow = c(1, n_rows), mar = c(4.5, 4.5, 3.2, 1.2))

for (ds_nm in ds_order) {
  d       <- ds_list[[ds_nm]]
  h0      <- which(d$harm_grid$alpha2 == 0)[1]
  B_vals  <- d$B_vals
  N_grid  <- d$N_grid
  b_cols  <- colorRampPalette(c("gray70", ds_cols[ds_nm]))(length(B_vals))
  b_pchs  <- c(16, 17, 15, 18)[seq_along(B_vals)]

  plot(NA, xlim = range(N_grid), ylim = c(0, 1), las = 1,
       xlab = "N", ylab = "Power (alpha2 = 0)",
       main = sprintf("%s\nProp. 1 check (pure cosinor)", d$label),
       xaxt = "n")
  axis(1, at = N_grid)
  abline(h = 0.80, lty = 3, col = "gray70")

  for (bi in seq_along(B_vals)) {
    pm <- d$power_mean[bi, , h0]
    lines(N_grid,  pm, col = b_cols[bi], lwd = 2)
    points(N_grid, pm, col = b_cols[bi], pch = b_pchs[bi], cex = 1.1)
  }
  legend("bottomright",
         legend = paste0("B=", B_vals),
         col = b_cols, lwd = 2, pch = b_pchs, bty = "n", cex = 0.8)
}
dev.off()
cat(sprintf("  -> %s\n", file.path(out_dir, "s4_bvsm_alpha0.pdf")))

# -----------------------------------------------------------------------
# Figure 3: Power drop at alpha2=0.5 vs alpha2=0, by B and N
# -----------------------------------------------------------------------
cat("Generating s4_bvsm_drop.pdf ...\n")

pdf(file.path(out_dir, "s4_bvsm_drop.pdf"),
    width = 5 * n_rows, height = 5)
par(mfrow = c(1, n_rows), mar = c(4.5, 4.5, 3.2, 1.2))

for (ds_nm in ds_order) {
  d       <- ds_list[[ds_nm]]
  h0      <- which(d$harm_grid$alpha2 == 0)[1]
  h05     <- which(abs(d$harm_grid$alpha2 - 0.5) < 1e-6)[1]
  if (is.na(h05)) {
    cat(sprintf("  %s: no alpha2=0.5 in harm_grid, skipping drop panel\n", ds_nm))
    plot.new(); title(main = sprintf("%s\n(no alpha2=0.5)", d$label))
    next
  }
  B_vals  <- d$B_vals
  N_grid  <- d$N_grid
  b_cols  <- colorRampPalette(c("gray70", ds_cols[ds_nm]))(length(B_vals))
  b_pchs  <- c(16, 17, 15, 18)[seq_along(B_vals)]

  drop_mat <- d$power_mean[, , h05] - d$power_mean[, , h0]  # [B, N]

  ylim <- range(drop_mat, na.rm = TRUE)
  ylim[1] <- min(ylim[1], -0.5)
  ylim[2] <- max(ylim[2],  0.05)

  plot(NA, xlim = range(N_grid), ylim = ylim, las = 1,
       xlab = "N", ylab = "Power drop  (alpha2=0.5 minus alpha2=0)",
       main = sprintf("%s\nPower cost of 2nd harmonic", d$label),
       xaxt = "n")
  axis(1, at = N_grid)
  abline(h = 0, lty = 1, col = "gray40")
  abline(h = c(-0.10, -0.20), lty = 3, col = "gray70")

  for (bi in seq_along(B_vals)) {
    dr <- drop_mat[bi, ]
    lines(N_grid,  dr, col = b_cols[bi], lwd = 2)
    points(N_grid, dr, col = b_cols[bi], pch = b_pchs[bi], cex = 1.1)
  }
  legend("bottomright",
         legend = paste0("B=", B_vals),
         col = b_cols, lwd = 2, pch = b_pchs, bty = "n", cex = 0.8)
}
dev.off()
cat(sprintf("  -> %s\n", file.path(out_dir, "s4_bvsm_drop.pdf")))

# -----------------------------------------------------------------------
# Print summary table: power at alpha2=0 and alpha2=0.5 for each dataset
# -----------------------------------------------------------------------
cat("\n--- Summary: power at alpha2=0 and alpha2=0.5 ---\n")
for (ds_nm in ds_order) {
  d    <- ds_list[[ds_nm]]
  h0   <- which(d$harm_grid$alpha2 == 0)[1]
  h05  <- which(abs(d$harm_grid$alpha2 - 0.5) < 1e-6)[1]
  cat(sprintf("\n%s (pilot_n=%d):\n", d$label, d$pilot_n))
  cat(sprintf("  %-10s", ""))
  cat(sprintf("  N=%-4d", d$N_grid), "\n")
  for (bi in seq_along(d$B_vals)) {
    cat(sprintf("  B=%-3d a2=0:   ", d$B_vals[bi]))
    cat(sprintf("  %.2f", d$power_mean[bi, , h0]), "\n")
    if (!is.na(h05)) {
      cat(sprintf("  B=%-3d a2=0.5: ", d$B_vals[bi]))
      cat(sprintf("  %.2f", d$power_mean[bi, , h05]), "\n")
    }
  }
}

cat(sprintf("\nDone. Output: %s/\n", out_dir))
