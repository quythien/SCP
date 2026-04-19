# Single-cohort Figure 1: 3-panel power summary
# Panel A: marginal power vs n (sim + CircaPower theoretical)
# Panel B: expected true discoveries vs n (sim mean ± SE + loess)
# Panel C: r-stratified power vs n (one line per stratum)

#' Plot single-cohort Figure 1 (3 panels)
#'
#' @param res      Output from \code{runSimsSingleCohort()}.
#' @param bio.opts \code{CircadianBioOptions} used to generate \code{res}.
#' @param out_pdf  Path for PDF output. If NULL, plots to current device.
#' @param title    Overall figure title (default empty).
#' @param alpha    Nominal alpha for CircaPower theoretical line (default 0.05).
#' @param strata_to_show Integer indices of r strata to plot in panel C (NULL = all).
#' @param width    PDF width in inches (default 10).
#' @param height   PDF height in inches (default 4).
#' @return Invisibly returns a list of the data frames used for each panel.
#' @export
plotSingleCohortFig1 <- function(res, bio.opts, out_pdf = NULL, title = "",
                                 alpha = 0.05, strata_to_show = NULL,
                                 width = 10, height = 4) {

  sample_sizes  <- res$sample_sizes
  nsims         <- res$nsims
  r_strata      <- res$r_strata
  strata_labels <- res$strata_labels
  n0_cp         <- res$n0_circapower

  # ------------------------------------------------------------------
  # Panel A data: marginal power
  # ------------------------------------------------------------------
  pow_mean <- rowMeans(res$marginal_power, na.rm = TRUE)
  pow_se   <- apply(res$marginal_power, 1, function(x) sd(x, na.rm=TRUE) / sqrt(sum(!is.na(x))))

  # CircaPower theoretical curve
  r_med <- if (!is.null(bio.opts$sigma_rhythmic))
    median(bio.opts$amplitude / bio.opts$sigma_rhythmic, na.rm = TRUE)
  else
    median(bio.opts$amplitude / exp(median(bio.opts$lOD, na.rm = TRUE)), na.rm = TRUE)

  n_grid  <- seq(min(sample_sizes) - 2, max(sample_sizes) + 10, by = 2)
  n_grid  <- n_grid[n_grid >= 4]
  cp_pow  <- vapply(n_grid, function(n) {
    lam <- r_med^2 * n / 2
    q   <- qf(1 - alpha, 2, n - 3)
    1 - pf(q, 2, n - 3, ncp = lam)
  }, numeric(1))

  df_a <- data.frame(n = sample_sizes, power = pow_mean, se = pow_se)
  df_cp <- data.frame(n = n_grid, power = cp_pow)

  # ------------------------------------------------------------------
  # Panel B data: expected true discoveries
  # ------------------------------------------------------------------
  td_mean <- rowMeans(res$marginal_TD, na.rm = TRUE)
  td_se   <- apply(res$marginal_TD, 1, function(x) sd(x, na.rm=TRUE) / sqrt(sum(!is.na(x))))
  df_b    <- data.frame(n = sample_sizes, td = td_mean, se = td_se)

  # ------------------------------------------------------------------
  # Panel C data: r-stratified power
  # ------------------------------------------------------------------
  n_strata <- length(strata_labels)
  if (is.null(strata_to_show)) strata_to_show <- seq_len(n_strata)

  df_c_list <- lapply(strata_to_show, function(k) {
    pm <- res$strat_power[, k, ]
    data.frame(
      n      = sample_sizes,
      power  = rowMeans(pm, na.rm = TRUE),
      se     = apply(pm, 1, function(x) sd(x, na.rm=TRUE) / sqrt(sum(!is.na(x)))),
      stratum = strata_labels[k]
    )
  })
  df_c <- do.call(rbind, df_c_list)

  # ------------------------------------------------------------------
  # Colours for panel C strata
  # ------------------------------------------------------------------
  pal <- colorRampPalette(c("#2b83ba", "#abdda4", "#fdae61", "#d7191c"))
  strat_cols <- setNames(pal(length(strata_to_show)), strata_labels[strata_to_show])

  # ------------------------------------------------------------------
  # Open device
  # ------------------------------------------------------------------
  if (!is.null(out_pdf)) pdf(out_pdf, width = width, height = height)
  old_par <- par(mfrow = c(1, 3), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
  on.exit({ par(old_par); if (!is.null(out_pdf)) dev.off() }, add = TRUE)

  # ---- Panel A ----
  ylim_a <- c(0, 1)
  plot(df_a$n, df_a$power, type = "n", ylim = ylim_a,
       xlab = "Sample size (n)", ylab = "Power",
       main = "A  Marginal power")
  abline(h = 0.8, lty = 2, col = "grey60")
  lines(df_cp$n, df_cp$power, col = "#999999", lwd = 1.5, lty = 3)
  .add_ribbon(df_a$n, df_a$power, df_a$se)
  lines(df_a$n, df_a$power, col = "#2b83ba", lwd = 2)
  points(df_a$n, df_a$power, pch = 19, col = "#2b83ba", cex = 0.8)
  if (!is.na(n0_cp)) abline(v = n0_cp, lty = 2, col = "#d7191c", lwd = 1.2)
  legend("bottomright", bty = "n", cex = 0.75,
         legend = c("Simulation", "CircaPower (median r)", "80% target", "n0 (CP)"),
         lty = c(1, 3, 2, 2), lwd = c(2, 1.5, 1, 1.2),
         col = c("#2b83ba", "#999999", "grey60", "#d7191c"))

  # ---- Panel B ----
  plot(df_b$n, df_b$td, type = "n",
       xlab = "Sample size (n)", ylab = "Expected true discoveries",
       main = "B  Expected discoveries")
  .add_ribbon(df_b$n, df_b$td, df_b$se)
  lines(df_b$n, df_b$td, col = "#2b83ba", lwd = 2)
  points(df_b$n, df_b$td, pch = 19, col = "#2b83ba", cex = 0.8)
  # Loess smoother
  lo <- tryCatch(loess(td ~ n, data = df_b, span = 0.75), error = function(e) NULL)
  if (!is.null(lo)) {
    n_pred <- seq(min(df_b$n), max(df_b$n), length.out = 100)
    lo_pred <- predict(lo, newdata = data.frame(n = n_pred))
    lines(n_pred, lo_pred, col = "#d7191c", lwd = 1.5, lty = 2)
  }
  legend("topleft", bty = "n", cex = 0.75,
         legend = c("Simulation mean \u00b1SE", "Loess"),
         lty = c(1, 2), lwd = c(2, 1.5), col = c("#2b83ba", "#d7191c"))

  # ---- Panel C ----
  strat_means <- lapply(strata_to_show, function(k) rowMeans(res$strat_power[, k, ], na.rm=TRUE))
  all_vals    <- unlist(strat_means)
  ylim_c      <- c(0, min(1, max(all_vals, na.rm=TRUE) * 1.1 + 0.05))
  plot(range(sample_sizes), ylim_c, type = "n",
       xlab = "Sample size (n)", ylab = "Power",
       main = "C  Power by r-stratum")
  abline(h = 0.8, lty = 2, col = "grey60")
  for (ki in seq_along(strata_to_show)) {
    k   <- strata_to_show[ki]
    col <- strat_cols[strata_labels[k]]
    pm  <- rowMeans(res$strat_power[, k, ], na.rm=TRUE)
    lines(sample_sizes, pm, col = col, lwd = 1.8)
    points(sample_sizes, pm, pch = 19, col = col, cex = 0.7)
  }
  legend("bottomright", bty = "n", cex = 0.68,
         legend = strata_labels[strata_to_show],
         lty = 1, lwd = 1.8,
         col = strat_cols[strata_labels[strata_to_show]])

  if (nchar(title) > 0) mtext(title, outer = TRUE, cex = 1.0, font = 2)

  invisible(list(panel_a = df_a, panel_b = df_b, panel_c = df_c,
                 circapower = df_cp, n0_circapower = n0_cp))
}


.add_ribbon <- function(x, y, se, col = "#2b83ba33") {
  polygon(c(x, rev(x)), c(y + se, rev(y - se)),
          col = col, border = NA)
}
