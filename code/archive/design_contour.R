# =============================================================================
# design_contour.R
# Optimal design analysis: TOD coverage (B time bins) vs. total sample size (n)
#
# Key idea: Under passive design, power depends on the effective sample size
#   n_eff = n * d,  where d in (0, 0.5] is the design efficiency.
# For fixed total n, spreading samples across more distinct time bins (larger B)
# increases d toward 0.5 (better temporal coverage), while fewer bins means more
# replicates per bin (higher n per bin but lower d). The contour plots show
# iso-power lines in (n, B) space, letting the user identify the cost-optimal
# design for a target power level.
#
# Main functions:
#   computeDesignEfficiency(times, period)
#   powerDesignGrid(bio.opts, n_grid, B_grid, nsims, cts, ...)
#   plotDesignContour(grid_results, test, target_powers, ...)
#
# Requires (sourced before use):
#   simulation.R  (simCircadianDiff, sampleTimesFromDist)
#   detection.R   (DCP_Rhythmicity, DCP_DiffR2, DCP_DiffPar)
#   options.R     (%||%)
# =============================================================================


#' Compute Design Efficiency for a Set of Collection Times
#'
#' @description Evaluates how uniformly a set of time points covers the 24-hour
#'   cycle for circadian power analysis. The design efficiency at gene phase
#'   \eqn{\phi} is
#'   \deqn{d(\phi) = \frac{1}{n}\sum_{i=1}^n \cos^2(\omega t_i - \omega\phi),}
#'   which equals 0.5 for perfectly equally-spaced times (all \eqn{\phi}) and
#'   varies across \eqn{\phi} for clustered times.
#'
#' @param times Numeric vector of collection times in [0, 24).
#' @param period Period in hours (default 24).
#' @param n_phi Number of gene-phase grid points (default 360).
#'
#' @return Named list:
#'   \item{d_mean}{Mean of d(phi); always \eqn{\approx 0.5}.}
#'   \item{d_min}{Worst-case d across all gene phases. Equals 0.5 for
#'     equally-spaced active design; smaller values indicate coverage gaps.}
#'   \item{d_sd}{SD of d(phi) across phases. Zero for equally-spaced design;
#'     larger values indicate unequal power across genes.}
#'   \item{d_eff}{d_mean - d_sd: penalises both low mean and high variability.}
computeDesignEfficiency <- function(times, period = 24, n_phi = 360) {
  omega    <- 2 * pi / period
  phi_grid <- seq(0, period, length.out = n_phi + 1)[-1]
  d_by_phi <- vapply(phi_grid, function(phi) {
    mean(cos(omega * times - omega * phi)^2)
  }, numeric(1))
  list(
    d_mean   = mean(d_by_phi),
    d_min    = min(d_by_phi),
    d_sd     = sd(d_by_phi),
    d_eff    = mean(d_by_phi) - sd(d_by_phi)
  )
}


# --- Internal time-generation helpers ----------------------------------------

# Create exact fixed times for an active design with B bins and n samples.
.makeFixedTimes <- function(B, n, period = 24) {
  bin_centers <- seq(0, period - period / B, by = period / B)
  m <- floor(n / B)
  r <- n %% B
  times <- rep(bin_centers, each = m)
  if (r > 0) times <- c(times, bin_centers[1:r])
  times
}


#' Run Power Simulations Over a (n, B) Design Grid
#'
#' @description Evaluates DR and DP power across a grid of total sample sizes
#'   \code{n} and number of distinct time bins \code{B}. For each cell, times
#'   are fixed at \code{B} equally-spaced bin centers with replicates per bin
#'   (active design). If \code{cts} is provided, an uncontrolled passive
#'   reference row is also run at each \code{n}.
#'
#' @param bio.opts A \code{CircadianBioOptions} object.
#' @param n_grid Integer vector of total sample sizes per group
#'   (default \code{c(20, 40, 60, 80, 100)}).
#' @param B_grid Integer vector of time bins per group
#'   (default \code{c(4, 6, 8, 12, 24)}). Cells with \code{B > n} are skipped.
#' @param nsims Number of simulation replicates per grid cell (default 100;
#'   use 200-500 for final figures).
#' @param cts Optional numeric vector of pilot collection times for an
#'   uncontrolled passive reference row.
#' @param test_types Tests to run: \code{"DR"}, \code{"DP"}, or both.
#' @param alpha BH FDR threshold (default 0.05).
#' @param verbose Print progress (default TRUE).
#'
#' @return A \code{PowerDesignGrid} data frame with columns:
#'   \code{n, B, design_type, d_min, d_sd, EDR_DR, EDR_DP, pFDR_DR, pFDR_DP,
#'   EDR_DR_se, EDR_DP_se}. (EDR is targeted power; pFDR is empirical FDR.)
powerDesignGrid <- function(bio.opts,
                             n_grid     = c(20, 40, 60, 80, 100),
                             B_grid     = c(4, 6, 8, 12, 24),
                             nsims      = 100,
                             cts        = NULL,
                             test_types = c("DR", "DP"),
                             alpha      = 0.05,
                             verbose    = TRUE) {

  stopifnot(inherits(bio.opts, "CircadianBioOptions"))
  period <- bio.opts$period %||% 24

  # Build (n, B) grid; drop cells where B > n
  grid <- expand.grid(n = n_grid, B = B_grid, stringsAsFactors = FALSE)
  grid <- grid[grid$B <= grid$n, ]
  grid$design_type <- "active"

  # Add uncontrolled passive reference if cts provided
  if (!is.null(cts)) {
    passive_rows <- data.frame(n = n_grid, B = NA_integer_,
                               design_type = "passive",
                               stringsAsFactors = FALSE)
    grid <- rbind(grid, passive_rows)
  }

  grid <- grid[order(grid$design_type, grid$n, grid$B), ]
  rownames(grid) <- NULL

  grid$d_min     <- NA_real_;  grid$d_sd      <- NA_real_
  grid$EDR_DR    <- NA_real_;  grid$EDR_DR_se <- NA_real_
  grid$EDR_DP    <- NA_real_;  grid$EDR_DP_se <- NA_real_
  grid$pFDR_DR   <- NA_real_;  grid$pFDR_DP   <- NA_real_

  n_cells <- nrow(grid)

  for (idx in seq_len(n_cells)) {
    n_i   <- grid$n[idx]
    B_i   <- grid$B[idx]
    dtype <- grid$design_type[idx]

    if (verbose) {
      if (is.na(B_i)) {
        cat(sprintf("[%d/%d] passive (pilot TOD),  n = %d\n", idx, n_cells, n_i))
      } else {
        cat(sprintf("[%d/%d] active B = %2d,  n = %3d\n", idx, n_cells, B_i, n_i))
      }
    }

    # Choose design and times for this cell
    if (dtype == "passive") {
      cts_i  <- cts
      cts2_i <- if (!is.null(bio.opts$cts2)) bio.opts$cts2 else cts_i
      design_i <- "passive"
    } else {
      # Active: fixed equally-spaced times with replicates per bin
      fixed_times <- .makeFixedTimes(B_i, n_i, period)
      attr(fixed_times, "fixed_times") <- TRUE
      cts_i  <- fixed_times
      cts2_i <- fixed_times
      design_i <- "active"
    }

    # Accumulators — use NA so that untested / undefined metrics stay NA,
    # not 0 (a numeric(nsims) init would silently average zeros for missing sims)
    edr_dr_v  <- rep(NA_real_, nsims);  pfdr_dr_v <- rep(NA_real_, nsims)
    edr_dp_v  <- rep(NA_real_, nsims);  pfdr_dp_v <- rep(NA_real_, nsims)
    d_min_v   <- rep(NA_real_, nsims);  d_sd_v    <- rep(NA_real_, nsims)

    for (s in seq_len(nsims)) {

      # Simulate two-group data (active or passive)
      sim_data <- simCircadianDiff(
        ngenes        = bio.opts$ngenes,
        n1            = n_i,
        n2            = n_i,
        lBaselineExpr = bio.opts$lBaselineExpr,
        lOD           = bio.opts$lOD,
        lOD2          = bio.opts$lOD2,         # F̂_σ2 for group 2 (NULL = share with group 1)
        amplitude     = bio.opts$amplitude,
        amplitude2    = bio.opts$amplitude2,   # F̂_A2 for g2-only DR (NULL = fall back to amplitude)
        prop_rhythmic = bio.opts$prop_rhythmic,
        prop_DR       = bio.opts$prop_DR,
        prop_DP       = bio.opts$prop_DP,
        prop_DA       = 0,
        phase_diff    = bio.opts$phase_diff,
        dp_shift_mode = bio.opts$dp_shift_mode %||% "uniform",
        period        = period,
        design        = design_i,
        cts           = cts_i,
        cts2          = cts2_i,
        sim.seed      = s * 1000L
      )

      # Design efficiency from group-1 realized times
      deff         <- computeDesignEfficiency(sim_data$times1, period)
      d_min_v[s]   <- deff$d_min
      d_sd_v[s]    <- deff$d_sd

      # Run DCP detection
      det <- .dcpDetect(sim_data, alpha = alpha, test_types = test_types)

      # DR metrics
      if ("DR" %in% test_types) {
        z_star_dr <- .zStarDR(sim_data)
        Z_dr      <- sim_data$ground_truth$diff_type %in% c(2, 3)
        res       <- .edrpFDR(det$DR_disc, Z_dr, z_star_dr)
        edr_dr_v[s]  <- res$EDR
        pfdr_dr_v[s] <- res$pFDR
      }

      # DP metrics
      if ("DP" %in% test_types) {
        z_star_dp <- .zStarDP(sim_data)
        Z_dp      <- sim_data$ground_truth$diff_type == 4
        res       <- .edrpFDR(det$DP_disc, Z_dp, z_star_dp)
        edr_dp_v[s]  <- res$EDR
        pfdr_dp_v[s] <- res$pFDR
      }
    }

    grid$d_min[idx]     <- mean(d_min_v)
    grid$d_sd[idx]      <- mean(d_sd_v)
    grid$EDR_DR[idx]    <- mean(edr_dr_v,  na.rm = TRUE)
    grid$EDR_DR_se[idx] <- sd(edr_dr_v,    na.rm = TRUE) / sqrt(sum(!is.na(edr_dr_v)))
    grid$EDR_DP[idx]    <- mean(edr_dp_v,  na.rm = TRUE)
    grid$EDR_DP_se[idx] <- sd(edr_dp_v,    na.rm = TRUE) / sqrt(sum(!is.na(edr_dp_v)))
    grid$pFDR_DR[idx]   <- mean(pfdr_dr_v, na.rm = TRUE)
    grid$pFDR_DP[idx]   <- mean(pfdr_dp_v, na.rm = TRUE)
  }

  class(grid) <- c("PowerDesignGrid", "data.frame")
  grid
}


# --- Detection wrapper -------------------------------------------------------

.dcpDetect <- function(sim_data, alpha, test_types) {
  n_genes    <- nrow(sim_data$expr1)
  gene_names <- rownames(sim_data$expr1)
  if (is.null(gene_names)) gene_names <- paste0("Gene", seq_len(n_genes))

  DR_disc <- rep(FALSE, n_genes)
  DP_disc <- rep(FALSE, n_genes)

  tryCatch({
    x1 <- format_for_DCP(sim_data$expr1, sim_data$times1, gene_names)
    x2 <- format_for_DCP(sim_data$expr2, sim_data$times2, gene_names)

    rhythm_res <- DCP_Rhythmicity(x1 = x1, x2 = x2,
                                  method = "Sidak_FS", period = 24,
                                  alpha = alpha, amp.cutoff = 0,
                                  CI = FALSE, p.adjust.method = "BH",
                                  parallel.ncores = 1)

    if ("DR" %in% test_types) {
      n_testable <- sum(rhythm_res$rhythm.joint$TOJR != "arrhy")
      if (n_testable > 0) {
        dr_res    <- DCP_DiffR2(rhythm_res, method = "LR", alpha = alpha)
        match_idx <- match(gene_names, dr_res$gname)
        pval_dr   <- rep(1, n_genes)
        pval_dr[!is.na(match_idx)] <- dr_res$p.R2[match_idx[!is.na(match_idx)]]
        DR_disc <- p.adjust(pval_dr, method = "BH") <= alpha
      }
    }

    if ("DP" %in% test_types) {
      n_testable <- sum(rhythm_res$rhythm.joint$TOJR == "both")
      if (n_testable > 0) {
        dp_res    <- DCP_DiffPar(rhythm_res, Par = "A&phase", alpha = alpha)
        match_idx <- match(gene_names, dp_res$gname)
        pval_dp   <- rep(1, n_genes)
        pval_dp[!is.na(match_idx)] <- dp_res$p.delta.peak[match_idx[!is.na(match_idx)]]
        DP_disc <- p.adjust(pval_dp, method = "BH") <= alpha
      }
    }
  }, error = function(e) {
    warning("DCP detection failed (sim skipped): ", conditionMessage(e))
  })

  list(DR_disc = DR_disc, DP_disc = DP_disc)
}


# --- Target indicators (mirror runner.R logic) --------------------------------

# DR target: diff_type in {2,3} AND max(r_g1, r_g2) >= 0.1
.zStarDR <- function(sim_data, threshold = 0.1) {
  dt <- sim_data$ground_truth$diff_type
  (dt %in% c(2, 3)) &
    (pmax(sim_data$effectsize_DR1, sim_data$effectsize_DR2) >= threshold)
}

# DP target: diff_type == 4 AND effective DP SNR >= 0.1
.zStarDP <- function(sim_data, threshold = 0.1) {
  dt <- sim_data$ground_truth$diff_type
  (dt == 4) & (sim_data$effectsize_phase >= threshold)
}

# Targeted power (EDR) and FDR from discovery + ground-truth vectors
.edrpFDR <- function(discoveries, Z_diff, z_star) {
  n_star <- sum(z_star, na.rm = TRUE)
  n_disc <- sum(discoveries, na.rm = TRUE)
  EDR  <- if (n_star > 0) sum(discoveries & z_star,  na.rm = TRUE) / n_star else NA_real_
  pFDR <- if (n_disc > 0) sum(discoveries & !Z_diff, na.rm = TRUE) / n_disc else NA_real_
  list(EDR = EDR, pFDR = pFDR)
}


# =============================================================================
# Plotting
# =============================================================================

#' Plot Design Contour: Iso-Power Curves in (n, B) Space
#'
#' @description Produces a filled-contour plot of targeted power (EDR) as a function of
#'   total sample size \code{n} (x-axis) and number of distinct time bins
#'   \code{B} (y-axis). Iso-power curves at \code{target_powers} are overlaid.
#'   An optional horizontal dashed line marks an uncontrolled passive reference.
#'
#' @param grid_results \code{PowerDesignGrid} from \code{powerDesignGrid()}.
#' @param test \code{"DR"} or \code{"DP"}.
#' @param target_powers Contour levels (default \code{c(0.2, 0.4, 0.6, 0.8)}).
#' @param show_pFDR If \code{TRUE}, add a second panel with the FDR surface.
#' @param passive_B Numeric: draw a horizontal reference line at this B value
#'   (e.g., the effective B of the uncontrolled pilot design).
#' @param title Plot title (auto-generated if \code{NULL}).
#'
#' @return Invisibly returns the targeted-power (EDR) matrix used for contouring
#'   (rows = B values, columns = n values).
plotDesignContour <- function(grid_results,
                               test          = c("DR", "DP"),
                               target_powers = c(0.2, 0.4, 0.6, 0.8),
                               show_pFDR     = FALSE,
                               passive_B     = NULL,
                               title         = NULL,
                               show_iso_nd   = FALSE) {

  test     <- match.arg(test)
  edr_col  <- paste0("EDR_",  test)
  pfdr_col <- paste0("pFDR_", test)

  # Quota-B cells only
  df <- grid_results[!is.na(grid_results$B) &
                       grid_results$design_type == "active", ]
  df <- df[order(df$n, df$B), ]

  n_vals <- sort(unique(df$n))
  B_vals <- sort(unique(df$B))

  # Build matrices: rows=B, cols=n
  edr_mat  <- matrix(NA_real_, nrow = length(B_vals), ncol = length(n_vals),
                     dimnames = list(B_vals, n_vals))
  pfdr_mat <- matrix(NA_real_, nrow = length(B_vals), ncol = length(n_vals),
                     dimnames = list(B_vals, n_vals))

  for (i in seq_len(nrow(df))) {
    ri <- as.character(df$B[i])
    ci <- as.character(df$n[i])
    edr_mat[ri, ci]  <- df[[edr_col]][i]
    pfdr_mat[ri, ci] <- df[[pfdr_col]][i]
  }

  if (is.null(title)) {
    title <- sprintf("%s Targeted power: iso-power curves in (sample size, time bins) space", test)
  }

  # Color palette (high-contrast)
  blue_to_red <- colorRampPalette(c("#2c7bb6","#00a6ca","#00ccbc","#90eb9d",
                                     "#ffff8c","#f9d057","#f29e2e","#e76818",
                                     "#d7191c"))(25)

  n_panels <- if (show_pFDR) 2L else 1L
  old_par  <- par(mfrow = c(1, n_panels), mar = c(5, 5, 4.5, 2.5) + 0.1)
  on.exit(par(old_par), add = TRUE)

  # Helper to draw one contour panel
  .draw_panel <- function(mat, col_pal, levels_fill, contour_levels,
                           main_title, ylab_str = "Number of time bins (B)",
                           key_title = "Targeted power") {
    # filled.contour: x=rows, y=cols of z, so z = t(mat) (rows=n, cols=B)
    filled.contour(
      x = n_vals,
      y = B_vals,
      z = t(mat),
      levels    = levels_fill,
      col       = col_pal,
      xlab      = "Sample size per group (n)",
      ylab      = ylab_str,
      main      = main_title,
      key.title = title(key_title, cex.main = 0.9),
      key.axes  = axis(4, las = 1, cex.axis = 0.8),
      plot.axes = {
        axis(1, at = n_vals, cex.axis = 0.9)
        axis(2, at = B_vals, cex.axis = 0.9)
        # Iso-power contours
        contour(x = n_vals, y = B_vals, z = t(mat),
                levels = contour_levels,
                add = TRUE, lwd = 2.2, col = "black",
                labcex = 0.9, drawlabels = TRUE)
        # Optional n * d ≈ constant reference lines (theoretical, d = 0.5 for active)
        if (show_iso_nd) {
          for (nd in c(10, 20, 40, 60, 80)) {
            abline(b = -1, a = nd / 0.5 / 10, lty = 3, col = "grey60", lwd = 0.8)
          }
        }
        # Passive reference
        if (!is.null(passive_B)) {
          abline(h = passive_B, lty = 2, col = "steelblue", lwd = 2)
          mtext(sprintf("passive (B \u2248 %g)", passive_B),
                side = 4, at = passive_B, las = 1, cex = 0.75,
                col = "steelblue", line = 0.3)
        }
      }
    )
  }

  # Panel 1: targeted power (EDR)
  z_min <- min(edr_mat, na.rm = TRUE)
  z_max <- max(edr_mat, na.rm = TRUE)
  if (!is.finite(z_min) || !is.finite(z_max) || z_min == z_max) {
    z_min <- 0; z_max <- 1
  }
  .draw_panel(
    mat           = edr_mat,
    col_pal       = blue_to_red,
    levels_fill   = seq(z_min, z_max, length.out = 21),
    contour_levels = target_powers,
    main_title    = title,
    key_title     = "Targeted power"
  )

  # Panel 2 (optional): FDR
  if (show_pFDR) {
    z_min_pf <- min(pfdr_mat, na.rm = TRUE)
    z_max_pf <- max(pfdr_mat, na.rm = TRUE)
    if (!is.finite(z_min_pf) || !is.finite(z_max_pf) || z_min_pf == z_max_pf) {
      z_min_pf <- 0; z_max_pf <- 0.5
    }
    .draw_panel(
      mat           = pfdr_mat,
      col_pal       = blue_to_red,
      levels_fill   = seq(z_min_pf, z_max_pf, length.out = 21),
      contour_levels = c(0.05, 0.10, 0.20),
      main_title    = sprintf("%s FDR surface\n(dashed = nominal 5%%)", test),
      key_title     = "FDR"
    )
    abline(h = 0.05, lty = 3, col = "red", lwd = 1.5)
  }

  invisible(edr_mat)
}


#' Plot Design Factor d(phi) vs Gene Phase
#'
#' @description Shows WHY B doesn't drive power and WHY active beats passive.
#'   For each design (active (equally-spaced) or passive pilot), plots the per-gene design
#'   factor d(phi) = (1/n) sum cos^2(omega*t_i - omega*phi) as a function of gene
#'   peak phase phi.  Active quota-B designs all give d(phi) = 0.5 (flat, by the
#'   geometric series identity), while the passive pilot has valleys (underpowered
#'   genes) and peaks (overpowered genes).
#'
#'   Since power = Pr[F > F_crit | lambda = n * r^2 * d(phi)], a flat d = 0.5
#'   means power depends only on n and r, not on phi or B.
#'
#' @param pilot_times Numeric vector of pilot TOD values (for passive curve).
#' @param B_active Integer vector of bin counts for active designs to overlay.
#' @param period Period in hours (default 24).
#' @param n_phi Number of phase grid points (default 360).
#' @param title Plot title.
#'
#' @return Invisibly returns a list with d_by_phi for each design.
#' @export
plotDesignFactor <- function(pilot_times,
                              B_active = c(4, 6, 12, 24),
                              period   = 24,
                              n_phi    = 360,
                              title    = "Design factor d(\u03c6) vs gene peak phase") {

  omega    <- 2 * pi / period
  phi_grid <- seq(0, period, length.out = n_phi + 1)[-1]

  # Passive: compute d(phi) from pilot times
  d_passive <- vapply(phi_grid, function(phi)
    mean(cos(omega * pilot_times - omega * phi)^2), numeric(1))

  # Active: d(phi) = 0.5 exactly for all phi and all B >= 2
  # (geometric series identity: sum_{k=0}^{B-1} e^{4*pi*i*k/B} = 0)
  d_active_const <- 0.5

  # Color palette
  B_cols <- setNames(
    colorRampPalette(c("#2166ac","#4dac26","#d6604d","#8073ac","#e08214"))(length(B_active)),
    as.character(B_active)
  )

  d_eff_val <- mean(d_passive) - sd(d_passive)
  d_min_val <- min(d_passive)
  d_max_val <- max(d_passive)

  # Full y-range: show the complete passive oscillation with margin
  ylo <- max(0,  d_min_val - 0.06)
  yhi <- min(1,  d_max_val + 0.06)
  yrange <- c(ylo, yhi)

  plot(phi_grid, d_passive,
       type = "n",
       xlim = c(0, period), ylim = yrange,
       xlab = "Gene peak phase (hours)",
       ylab = "Design factor  d",
       main = title,
       xaxt = "n", las = 1)
  axis(1, at = seq(0, period, by = 4),
       labels = paste0(seq(0, period, by = 4), "h"))

  # Horizontal grid
  abline(h = seq(0, 1, by = 0.1), col = "grey92", lwd = 0.6)

  # Shade regions where passive < 0.5 (power loss) and > 0.5 (power gain)
  # Power loss (below 0.5): red shading
  loss_idx <- d_passive < 0.5
  if (any(loss_idx)) {
    polygon(c(phi_grid, rev(phi_grid)),
            c(pmin(d_passive, 0.5), rep(0.5, length(phi_grid))),
            col = adjustcolor("#d73027", alpha.f = 0.18), border = NA)
  }
  # Power gain (above 0.5): blue shading
  gain_idx <- d_passive > 0.5
  if (any(gain_idx)) {
    polygon(c(phi_grid, rev(phi_grid)),
            c(pmax(d_passive, 0.5), rep(0.5, length(phi_grid))),
            col = adjustcolor("#4575b4", alpha.f = 0.18), border = NA)
  }

  # Active reference line at d = 0.5 (one line, labeled)
  abline(h = 0.5, col = "#2166ac", lwd = 2.5, lty = 1)

  # Passive curve — drawn last so it is always on top
  lines(phi_grid, d_passive, col = "black", lwd = 2.2, lty = 1)

  # Dotted lines for d_min and d_max of passive
  abline(h = d_min_val, col = "#d73027", lty = 2, lwd = 1.2)
  abline(h = d_max_val, col = "#4575b4", lty = 2, lwd = 1.2)

  # Annotations (right-side)
  x_ann <- period * 0.99
  text(x_ann, 0.5 + (yhi - 0.5) * 0.25,
       sprintf("Active (any B): d = 0.50"),
       adj = 1, cex = 0.75, col = "#2166ac", font = 2)
  text(x_ann, d_min_val - (d_min_val - ylo) * 0.4,
       sprintf("Passive d_min = %.2f", d_min_val),
       adj = 1, cex = 0.72, col = "#d73027")
  text(x_ann, d_max_val + (yhi - d_max_val) * 0.4,
       sprintf("Passive d_max = %.2f", d_max_val),
       adj = 1, cex = 0.72, col = "#4575b4")

  legend("bottomleft", bty = "n", cex = 0.82,
         legend = c("Passive pilot TOD: d varies by gene phase",
                    "Active (equally-spaced, any B >= 2): d = 0.5 always",
                    "Power loss vs active  (passive < 0.5)",
                    "Power gain vs active  (passive > 0.5)"),
         col    = c("black", "#2166ac", "#d73027",      "#4575b4"),
         lty    = c(1,       1,         NA,              NA),
         lwd    = c(2.2,     2.5,       NA,              NA),
         fill   = c(NA,      NA,
                    adjustcolor("#d73027", 0.3),
                    adjustcolor("#4575b4", 0.3)),
         border = c(NA, NA, NA, NA))

  invisible(list(phi_grid  = phi_grid,
                 d_passive = d_passive,
                 d_active  = 0.5,
                 d_eff     = d_eff_val,
                 d_min     = d_min_val,
                 d_max     = d_max_val))
}


#' Plot Active vs Passive Power as Line Chart
#'
#' @description Line plot of targeted power (EDR) vs sample size n, with one line per B value
#'   (active (equally-spaced) designs) and a separate line for the passive pilot reference.
#'   More interpretable than a filled contour when B has little effect on power.
#'
#' @param grid_results \code{PowerDesignGrid} from \code{powerDesignGrid()}.
#' @param test \code{"DR"} or \code{"DP"}.
#' @param title Plot title (auto-generated if \code{NULL}).
#' @param target_power Optional horizontal reference line for a target power level.
#'
#' @return Invisibly returns a data.frame with the plotted values.
#' @export
plotDesignLines <- function(grid_results,
                             test         = c("DR", "DP"),
                             title        = NULL,
                             target_power = NULL) {

  test    <- match.arg(test)
  edr_col <- paste0("EDR_", test)
  se_col  <- paste0("EDR_", test, "_se")

  df      <- as.data.frame(grid_results)
  active  <- df[!is.na(df$B) & df$design_type == "active", ]
  passive <- df[!is.na(df$design_type) & df$design_type == "passive", ]

  n_vals  <- sort(unique(df$n))
  B_vals  <- sort(unique(active$B))

  if (is.null(title))
    title <- sprintf("%s Targeted power: active (equally-spaced) vs passive pilot", test)

  # Color palette for B lines
  B_cols <- setNames(
    colorRampPalette(c("#2166ac", "#4dac26", "#d6604d", "#8073ac", "#e08214"))(length(B_vals)),
    as.character(B_vals)
  )

  # y-axis range
  all_edr <- c(active[[edr_col]], passive[[edr_col]])
  edr_max  <- suppressWarnings(max(all_edr, na.rm = TRUE))
  if (!is.finite(edr_max)) {
    message(sprintf("plotDesignLines: no finite targeted-power values for %s — skipping plot.", test))
    return(invisible(NULL))
  }
  ymax    <- min(1, edr_max * 1.10 + 0.05)

  plot(NA, xlim = range(n_vals), ylim = c(0, ymax),
       xlab = "Sample size per group (n)",
       ylab = sprintf("Targeted power (%s)", test),
       main = title,
       xaxt = "n", las = 1)
  axis(1, at = n_vals)
  abline(h = seq(0, 1, by = 0.1), col = "grey90", lwd = 0.5)
  abline(v = n_vals,              col = "grey90", lwd = 0.5)

  # Active lines: one per B
  for (b in B_vals) {
    sub  <- active[active$B == b, ]
    sub  <- sub[order(sub$n), ]
    bchr <- as.character(b)
    lines(sub$n, sub[[edr_col]], col = B_cols[bchr], lwd = 2)
    points(sub$n, sub[[edr_col]], col = B_cols[bchr], pch = 16, cex = 0.9)
    # Error ribbon (±1 SE)
    if (se_col %in% names(sub)) {
      polygon(c(sub$n, rev(sub$n)),
              c(sub[[edr_col]] + sub[[se_col]], rev(sub[[edr_col]] - sub[[se_col]])),
              col = adjustcolor(B_cols[bchr], alpha.f = 0.15), border = NA)
    }
  }

  # Passive line
  if (nrow(passive) > 0) {
    passive <- passive[order(passive$n), ]
    lines(passive$n, passive[[edr_col]], col = "black", lwd = 2.2, lty = 2)
    points(passive$n, passive[[edr_col]], col = "black", pch = 4, cex = 1)
  }

  # Target power reference
  if (!is.null(target_power))
    abline(h = target_power, lty = 3, col = "red", lwd = 1.5)

  # Legend
  legend("topleft", bty = "n", cex = 0.82,
         legend = c(paste0("Active B=", B_vals), "Passive (pilot TOD)"),
         col    = c(B_cols[as.character(B_vals)], "black"),
         lty    = c(rep(1, length(B_vals)), 2),
         lwd    = c(rep(2, length(B_vals)), 2.2),
         pch    = c(rep(16, length(B_vals)), 4))

  invisible(rbind(
    cbind(active[, c("n","B","design_type", edr_col)]),
    cbind(passive[, c("n","B","design_type", edr_col)])
  ))
}


#' Plot Power in (N, m) Space: Sample Size vs Replicates per Time Point
#'
#' @description Plots active-design power with total sample size N on the
#'   x-axis and replicates per time point m = N/B on the y-axis.  Each data
#'   point from the simulation grid appears as a labelled bubble.  Diagonal
#'   lines show iso-B constraints (m = N/B for each B value).  Vertical colour
#'   bands, derived by averaging power across all B values at each N, visually
#'   confirm that power is driven by N not by m or B.
#'
#' @param grid_results \code{PowerDesignGrid} from \code{powerDesignGrid()}.
#' @param test \code{"DR"} or \code{"DP"}.
#' @param title Plot title (auto-generated if \code{NULL}).
#' @param passive Logical: overlay passive reference points (default \code{TRUE}).
#'
#' @return Invisibly returns the data frame used for plotting.
#' @export
plotDesignContourNM <- function(grid_results,
                                 test    = c("DR", "DP"),
                                 title   = NULL,
                                 passive = TRUE) {

  test    <- match.arg(test)
  edr_col <- paste0("EDR_", test)

  # Active cells only; compute m = N / B
  active <- as.data.frame(grid_results[
    !is.na(grid_results$B) & grid_results$design_type == "active", ])
  active$m <- active$n / active$B

  n_vals <- sort(unique(active$n))
  B_vals <- sort(unique(active$B))

  if (is.null(title))
    title <- sprintf(
      "%s Power: sample size (N) vs replicates per time point (m = N/B)\n[diagonal lines = iso-B; vertical bands = power driven by N]", test)

  # ---- Colour scale for power ----
  edr_range <- range(active[[edr_col]], na.rm = TRUE)
  if (diff(edr_range) < 1e-6) edr_range <- c(0, 1)
  n_cols   <- 100
  col_pal  <- colorRampPalette(c("#313695","#4575b4","#74add1","#abd9e9",
                                  "#ffffbf","#fee090","#f46d43","#a50026"))(n_cols)
  edr_to_col <- function(edr) {
    idx <- round((edr - edr_range[1]) / diff(edr_range) * (n_cols - 1)) + 1
    col_pal[pmax(1, pmin(n_cols, idx))]
  }

  # ---- Vertical band colours: mean targeted power (EDR) per N ----
  n_mean_edr <- tapply(active[[edr_col]], active$n, mean, na.rm = TRUE)

  # Y-axis range: from 0 to max(m) + padding
  m_max <- max(active$m) * 1.10

  # ---- Base plot ----
  old_par <- par(mar = c(5, 5, 4, 6))
  on.exit(par(old_par), add = TRUE)

  plot(NA, xlim = range(n_vals) + c(-8, 8), ylim = c(0, m_max),
       xlab = "Sample size per group (N)",
       ylab = "Replicates per time point  (m = N / B)",
       main = title, xaxt = "n", las = 1)
  axis(1, at = n_vals)

  # ---- Draw vertical colour bands (width = gap between n values) ----
  n_sorted <- sort(n_vals)
  half_gap <- min(diff(n_sorted)) / 2
  for (nv in n_sorted) {
    rect(nv - half_gap, 0, nv + half_gap, m_max,
         col = adjustcolor(edr_to_col(n_mean_edr[as.character(nv)]), 0.35),
         border = NA)
  }

  # ---- Grid lines ----
  abline(v = n_vals, col = "white", lwd = 0.5)
  abline(h = pretty(c(0, m_max), 8), col = "grey85", lwd = 0.5)

  # ---- Iso-B diagonal lines: m = N / B ----
  n_ext <- c(min(n_vals) - half_gap, max(n_vals) + half_gap)
  B_cols <- setNames(
    colorRampPalette(c("#1b7837","#762a83","#e08214","#2166ac","#d73027"))(length(B_vals)),
    as.character(B_vals))
  for (b in B_vals) {
    m_ends <- n_ext / b
    lines(n_ext, m_ends, col = B_cols[as.character(b)], lwd = 1.6, lty = 2)
    # Label at right end, just inside
    text(max(n_vals) + half_gap * 0.8, max(n_vals) / b,
         sprintf("B=%d", b), adj = 0, cex = 0.72,
         col = B_cols[as.character(b)], font = 2, xpd = TRUE)
  }

  # ---- Passive reference points ----
  if (passive) {
    pass_df <- as.data.frame(grid_results[
      !is.na(grid_results$design_type) & grid_results$design_type == "passive", ])
    if (nrow(pass_df) > 0) {
      # Passive has no B; plot along m=0 baseline with a rug
      for (i in seq_len(nrow(pass_df))) {
        nv <- pass_df$n[i]
        points(nv, 0.6, pch = 23, bg = edr_to_col(pass_df[[edr_col]][i]),
               col = "black", cex = 2.2, lwd = 1.2)
        text(nv, 0.6, sprintf("%.0f%%", 100 * pass_df[[edr_col]][i]),
             cex = 0.52, col = "white", font = 2)
      }
      text(min(n_vals) - half_gap * 0.8, 0.6, "Passive",
           adj = 1, cex = 0.72, col = "grey30", xpd = TRUE)
    }
  }

  # ---- Data bubbles: one per (N, B) cell ----
  for (i in seq_len(nrow(active))) {
    nv  <- active$n[i]
    mv  <- active$m[i]
    edr <- active[[edr_col]][i]
    se  <- active[[paste0(edr_col, "_se")]][i]
    points(nv, mv, pch = 21, bg = edr_to_col(edr),
           col = "grey20", cex = 2.8, lwd = 0.8)
    text(nv, mv, sprintf("%.0f%%", 100 * edr),
         cex = 0.52, col = "white", font = 2)
  }

  # ---- Colour-scale legend (right side) ----
  legend_x <- grconvertX(1.01, from = "npc", to = "user")
  legend_y0 <- grconvertY(0.15, from = "npc", to = "user")
  legend_y1 <- grconvertY(0.85, from = "npc", to = "user")
  n_leg <- 20
  y_leg <- seq(legend_y0, legend_y1, length.out = n_leg + 1)
  leg_cols <- colorRampPalette(c("#313695","#4575b4","#74add1","#abd9e9",
                                  "#ffffbf","#fee090","#f46d43","#a50026"))(n_leg)
  rect_w <- half_gap * 0.6
  for (k in seq_len(n_leg)) {
    rect(legend_x, y_leg[k], legend_x + rect_w, y_leg[k + 1],
         col = leg_cols[k], border = NA, xpd = TRUE)
  }
  text(legend_x + rect_w + 0.5, legend_y0,
       sprintf("%.0f%%", 100 * edr_range[1]), adj = 0, cex = 0.65, xpd = TRUE)
  text(legend_x + rect_w + 0.5, legend_y1,
       sprintf("%.0f%%", 100 * edr_range[2]), adj = 0, cex = 0.65, xpd = TRUE)
  text(legend_x + rect_w / 2, legend_y1 + diff(range(y_leg)) * 0.06,
       "Targeted power", adj = 0.5, cex = 0.75, font = 2, xpd = TRUE)

  invisible(active)
}


#' Print Method for PowerDesignGrid
#' @export
print.PowerDesignGrid <- function(x, ...) {
  cat("PowerDesignGrid\n")
  cat(sprintf("  Total cells:  %d\n", nrow(x)))
  cat(sprintf("  n range:      [%g, %g]\n", min(x$n), max(x$n)))
  B_vals <- x$B[!is.na(x$B)]
  if (length(B_vals) > 0)
    cat(sprintf("  B range:      [%g, %g]\n", min(B_vals), max(B_vals)))
  if ("passive" %in% x$design_type)
    cat("  Passive ref:  included\n")
  if (any(!is.na(x$EDR_DR)))
    cat(sprintf("  Targeted power (DR) range: [%.3f, %.3f]\n",
                min(x$EDR_DR, na.rm = TRUE), max(x$EDR_DR, na.rm = TRUE)))
  if (any(!is.na(x$EDR_DP)))
    cat(sprintf("  Targeted power (DP) range: [%.3f, %.3f]\n",
                min(x$EDR_DP, na.rm = TRUE), max(x$EDR_DP, na.rm = TRUE)))
  invisible(x)
}


#' Compute Simulation-Based Design Utility and Return Optimal Design
#'
#' @description
#' Evaluates the utility function
#'   U(D) = w_DR * EDR_DR + w_DP * EDR_DP - lambda * 1[FDR > alpha]
#' for each row of a \code{PowerDesignGrid} result and returns the
#' optimal design D* (Eq. 13 in the paper).  Under equally-spaced active
#' designs, Proposition 1 guarantees that U depends only on n, so D*
#' is the minimum n that meets the target power.
#'
#' @param grid_results  A \code{PowerDesignGrid} data frame from
#'   \code{powerDesignGrid()}.
#' @param w_DR   Weight on targeted power (EDR_DR) (default 1).
#' @param w_DP   Weight on targeted power (EDR_DP) (default 0).
#' @param lambda Penalty for FDR exceeding \code{alpha} (default 1).
#' @param alpha  Nominal FDR threshold (default 0.05).
#' @param target_power  If provided, D* is additionally constrained to
#'   cells with U >= target_power (printed separately).
#' @param verbose Print summary table (default TRUE).
#'
#' @return Invisibly returns the input data frame with an added column
#'   \code{utility}, sorted descending by utility.
#' @export
computeDesignUtility <- function(grid_results,
                                  w_DR         = 1,
                                  w_DP         = 0,
                                  lambda       = 1,
                                  alpha        = 0.05,
                                  target_power = NULL,
                                  verbose      = TRUE) {

  df <- as.data.frame(grid_results)

  # FDR penalty: use DR pFDR if w_DR > 0, else DP pFDR
  pfdr_col <- if (w_DR >= w_DP && "pFDR_DR" %in% names(df)) "pFDR_DR" else "pFDR_DP"
  pfdr_val <- if (pfdr_col %in% names(df)) df[[pfdr_col]] else rep(0, nrow(df))

  edr_dr <- if ("EDR_DR" %in% names(df)) df$EDR_DR else rep(0, nrow(df))
  edr_dp <- if ("EDR_DP" %in% names(df)) df$EDR_DP else rep(0, nrow(df))

  # Replace NaN/NA with 0 when the corresponding weight is 0,
  # so that 0 * NaN = 0 rather than NaN
  dr_contrib <- if (w_DR == 0) rep(0, nrow(df)) else w_DR * replace(edr_dr, is.na(edr_dr), 0)
  dp_contrib <- if (w_DP == 0) rep(0, nrow(df)) else w_DP * replace(edr_dp, is.na(edr_dp), 0)

  df$utility <- dr_contrib +
                dp_contrib -
                lambda * as.numeric(!is.na(pfdr_val) & pfdr_val > alpha)

  df <- df[order(-df$utility, df$n), ]

  if (verbose) {
    cat("=== Design Utility U(D) ===\n")
    cat(sprintf("  U(D) = %.2f * TargetedPower_DR + %.2f * TargetedPower_DP - %.2f * 1[FDR > %.2f]\n\n",
                w_DR, w_DP, lambda, alpha))

    print_cols <- intersect(c("design_type","n","B","EDR_DR","EDR_DP",
                               pfdr_col,"utility"), names(df))
    top <- head(df[, print_cols], 10)
    for (col in c("EDR_DR","EDR_DP","utility")) {
      if (col %in% names(top)) top[[col]] <- round(top[[col]], 3)
    }
    print(top, row.names = FALSE)

    dstar <- df[1, ]
    cat(sprintf("\nOptimal design D*:\n"))
    cat(sprintf("  design_type = %s\n", dstar$design_type))
    cat(sprintf("  n           = %g\n", dstar$n))
    if (!is.na(dstar$B)) cat(sprintf("  B           = %g\n", dstar$B))
    cat(sprintf("  U(D*)       = %.3f\n", dstar$utility))

    if (!is.null(target_power)) {
      feasible <- df[!is.na(df$utility) & df$utility >= target_power, ]
      if (nrow(feasible) == 0) {
        cat(sprintf("\n  No design meets target power %.2f in the grid.\n",
                    target_power))
      } else {
        # minimum n among feasible
        dstar_n <- feasible[which.min(feasible$n), ]
        cat(sprintf("\nSmallest n meeting U >= %.2f:\n", target_power))
        cat(sprintf("  design_type = %s,  n = %g",
                    dstar_n$design_type, dstar_n$n))
        if (!is.na(dstar_n$B)) cat(sprintf(",  B = %g", dstar_n$B))
        cat(sprintf(",  U = %.3f\n", dstar_n$utility))
      }
    }
  }

  invisible(df)
}
