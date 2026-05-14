#' =======================================================================
#' plot_fmm.R — FMM Waveform Robustness Plots (Figures 4, 4b)
#' =======================================================================
#'
#' Plotting functions for FMM-based cosinor violation analyses.
#' Two main functions:
#'   plotFMMViolation()     — single-cohort 4-row figure (Fig 4)
#'   plotFMMDifferential()  — differential 3-panel figure (Fig 4b)

# =====================================================================
# Single-cohort FMM plot (Figure 4)
# =====================================================================

#' Plot FMM Waveform Robustness: Single-Cohort
#'
#' @description
#' Produces a 2-row figure (one panel per sweep type) showing DCP power under
#' FMM cosinor violation: row 1 sweeps omega at fixed beta; row 2 sweeps beta
#' at fixed omega. SD bars are drawn at each N. Each panel is faceted by
#' design type (active / passive) and SNR column (strong / moderate / weak).
#'
#' The FMM model (Rueda, Rodríguez-Collado & Peddada, 2019, Sci Rep,
#' doi:10.1038/s41598-019-54569-1) parameterises waveform shape via
#' \code{omega} and orientation via \code{beta}. \code{omega = 1} recovers
#' pure cosinor, \code{omega = 0} gives a flat (arrhythmic) signal.
#'
#' @param x          Data frame with columns \code{N, sweep_val, sweep_param,
#'                   omega, beta, power, power_se, dataset, design_row, snr_col}.
#'                   Produced by \code{runFMMViolationAnalysis()} or the
#'                   \code{fig4_fmm_violation.R} script.
#' @param nsims      Number of simulation replicates (for SE → SD conversion).
#' @param omega_fixed Omega value used in the beta sweep rows (for subtitle).
#' @param tissue_annotations Optional \code{data.frame} with columns
#'   \code{design_type}, \code{snr_col}, \code{label}, \code{N}, \code{y}
#'   providing per-panel text annotations (e.g. dataset name and median
#'   r-tilde). If \code{NULL} (default), no annotations are drawn.  The
#'   caller is responsible for populating this from their pilot data — hard-coded
#'   tissue names must not appear inside this package function.
#' @param output_file PDF path. NULL = current device.
#' @param width PDF width in inches (default 16).
#' @param height PDF height in inches (default 14).
#'
#' @return Invisibly returns the combined patchwork object.
#' @seealso \code{\link{simCircadianFMM}}, \code{\link{plotFMMDifferential}}
#' @export
plotFMMViolation <- function(x,
                              nsims             = 30L,
                              omega_fixed       = 0.5,
                              tissue_annotations = NULL,
                              output_file       = NULL,
                              width             = 16,
                              height            = 14) {

  for (pkg in c("ggplot2","patchwork"))
    if (!requireNamespace(pkg, quietly = TRUE))
      stop(sprintf("Package '%s' is required for plotFMMViolation().", pkg))

  full_df <- x
  if (!"power_sd" %in% names(full_df))
    full_df$power_sd <- full_df$power_se * sqrt(nsims)

  full_df$snr_col <- factor(full_df$snr_col, levels = c("Strong","Moderate","Weak"))

  # ── Color schemes ─────────────────────────────────────────────
  omega_levels <- as.character(sort(unique(full_df$omega[full_df$sweep_param=="omega"])))
  n_omega      <- length(omega_levels)
  omega_colors <- stats::setNames(
    grDevices::colorRampPalette(c("#08306B","#4292C6","#C6DBEF"))(n_omega),
    omega_levels)
  # Parseable labels: "omega==0" renders as ω=0 via plotmath
  omega_labels <- stats::setNames(
    paste0("omega==", omega_levels),
    omega_levels)

  beta_levels  <- as.character(round(sort(unique(full_df$sweep_val[full_df$sweep_param=="beta"])),4))
  n_beta       <- length(beta_levels)
  beta_colors  <- stats::setNames(
    grDevices::colorRampPalette(c("#7F0000","#EF6548","#FEE8C8"))(n_beta),
    beta_levels)
  # Parseable beta labels: "pi/4" renders as π/4
  beta_denom  <- c("0","pi/4","pi/2","3*pi/4","pi","5*pi/4","3*pi/2")
  beta_labels <- stats::setNames(beta_denom[seq_len(n_beta)], beta_levels)

  # Alpha sweep color scheme (green palette: 0h dark → 20h light)
  has_alpha    <- "alpha" %in% full_df$sweep_param
  alpha_levels <- if (has_alpha)
    as.character(sort(unique(full_df$sweep_val[full_df$sweep_param=="alpha"]))) else character(0)
  n_alpha      <- length(alpha_levels)
  alpha_colors <- stats::setNames(
    grDevices::colorRampPalette(c("#004529","#41AB5D","#C7E9C0"))(max(n_alpha,1)),
    alpha_levels)
  alpha_labels <- stats::setNames(paste0("alpha==", alpha_levels, "*'h'"), alpha_levels)

  N_all    <- sort(unique(full_df$N))
  x_breaks <- N_all   # show all N values on x-axis

  theme_base <- ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      strip.background  = ggplot2::element_rect(fill = "grey92"),
      strip.text.y      = ggplot2::element_text(size = 8),
      strip.text.x      = ggplot2::element_text(size = 10, face = "bold"),
      legend.position   = "bottom",
      panel.grid.minor  = ggplot2::element_blank(),
      axis.text.x       = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
      plot.title        = ggplot2::element_text(hjust = 0.5, face = "bold", size = 11),
      plot.caption      = ggplot2::element_text(hjust = 0, size = 8, colour = "grey40")
    )

  # Per-panel text annotation (dataset name, median r-tilde, etc.)
  # Supplied by the caller via tissue_annotations; NULL disables annotations.
  tissue_ann <- tissue_annotations
  if (!is.null(tissue_ann)) {
    tissue_ann$design_type <- factor(tissue_ann$design_type,
      levels = c("Active (B=12, 2h)", "Passive"))
    tissue_ann$snr_col <- factor(tissue_ann$snr_col,
      levels = c("Strong", "Moderate", "Weak"))
  }

  make_panel <- function(sub_df, fac_levels, colors, labels, legend_name,
                          x_lab, title_str, legend_ncol = 1, jitter_x = 0) {
    sub_df$sweep_fac   <- factor(as.character(round(sub_df$sweep_val, 4L)),
                                  levels = fac_levels)
    sub_df$design_type <- ifelse(grepl("^Active", sub_df$design_row),
                                  "Active (B=12, 2h)", "Passive")
    sub_df$design_type <- factor(sub_df$design_type,
                                  levels = c("Active (B=12, 2h)","Passive"))
    sub_df$snr_col     <- factor(sub_df$snr_col, levels = c("Strong","Moderate","Weak"))

    # Apply x-jitter by offsetting N per sweep level (helps when lines overlap)
    if (jitter_x > 0) {
      lev_idx <- as.integer(sub_df$sweep_fac)
      n_lev   <- length(fac_levels)
      offsets <- seq(-(n_lev-1)/2, (n_lev-1)/2, length.out = n_lev) * jitter_x
      sub_df$N_jit <- sub_df$N + offsets[lev_idx]
    } else {
      sub_df$N_jit <- sub_df$N
    }

    # Build facet_grid plot with legend inside the empty panel
    p <- ggplot2::ggplot(sub_df, ggplot2::aes(
        x = N_jit, y = 100 * power, colour = sweep_fac,
        group = interaction(sweep_fac, dataset))) +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = 100*(power - power_sd),
                     ymax = 100*(power + power_sd)),
        width = diff(range(N_all)) * 0.012, linewidth = 0.45, alpha = 0.8) +
      ggplot2::geom_point(size = 1.5) +
      ggplot2::geom_hline(yintercept = 80, linetype = "dashed",
                          colour = "grey50", linewidth = 0.4) +
      ggplot2::facet_grid(design_type ~ snr_col) +
      ggplot2::scale_colour_manual(
        values = colors,
        labels = parse(text = unname(labels)),
        name   = legend_name) +
      ggplot2::guides(colour = ggplot2::guide_legend(
        ncol = legend_ncol,
        override.aes = list(size = 3)
      )) +
      ggplot2::scale_x_continuous(breaks = x_breaks) +
      ggplot2::scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
      ggplot2::labs(x = x_lab, y = "Power (%)", title = title_str) +
      theme_base +
      ggplot2::theme(
        legend.position          = "inside",
        legend.position.inside   = c(0.175, 0.26),
        legend.justification     = c(0.5, 0.5),
        legend.direction         = if (legend_ncol > 1) "horizontal" else "vertical",
        legend.background        = ggplot2::element_rect(fill="white", colour="grey70", linewidth=0.4),
        legend.key.size          = grid::unit(0.7, "lines"),
        legend.text              = ggplot2::element_text(size = 10),
        legend.title             = ggplot2::element_text(size = 11, face = "bold"),
        legend.margin            = ggplot2::margin(5, 8, 5, 8)
      )

    # Add optional per-panel dataset label annotations supplied by the caller
    if (!is.null(tissue_ann)) {
      p <- p + ggplot2::geom_text(
        data = tissue_ann,
        ggplot2::aes(x = N, y = y, label = label),
        inherit.aes = FALSE, hjust = 0, vjust = 1,
        size = 2.4, colour = "grey20", lineheight = 0.9)
    }

    # Post-process gtable: blank all elements of the empty (Passive/Strong) cell
    g <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(p))

    # 1. Blank: panel content, x-axis below Strong column, panel border
    for (nm in c("panel-1-2", "axis-b-1")) {
      idx <- which(g$layout$name == nm)
      for (i in idx) g$grobs[[i]] <- grid::nullGrob()
    }

    # 2. Move y-axis from left of empty panel (l=6) to left of Passive/Moderate (l=8)
    #    axis-l-2: y-axis for the entire Passive row, currently at l=6
    #    panel-2-2: Passive/Moderate panel, at l=9 → its left gap is l=8
    al2_idx <- which(g$layout$name == "axis-l-2")
    if (length(al2_idx) > 0) {
      # Save the y-axis grob, then blank it at its original position
      al2_grob  <- g$grobs[[al2_idx[1L]]]
      al2_row_t <- g$layout$t[al2_idx[1L]]
      al2_row_b <- g$layout$b[al2_idx[1L]]
      for (i in al2_idx) g$grobs[[i]] <- grid::nullGrob()

      # Place at l=8 (spacer column just left of Passive/Moderate panel at l=9)
      g <- gtable::gtable_add_grob(g, al2_grob,
        t = al2_row_t, b = al2_row_b,
        l = 8L, r = 8L,
        name = "axis-l-passive-moved")
    }

    cowplot::ggdraw(g)
  }

  # make_panel now handles legend positioning + gtable blanking internally
  p_omega <- make_panel(
    full_df[full_df$sweep_param=="omega",],
    omega_levels, omega_colors, omega_labels,
    legend_name = expression(omega),
    x_lab       = "",
    title_str   = expression(paste("(A) Waveform shape: varying ", omega,
                                    "   (fixed ", beta, " = ", pi, ")")),
    legend_ncol = 4L   # horizontal: 7 omega values → 2 rows of 4+3
  )

  p_beta <- make_panel(
    full_df[full_df$sweep_param=="beta",],
    beta_levels, beta_colors, beta_labels,
    legend_name = expression(beta),
    x_lab       = if (has_alpha) "" else "Total sample size N",
    title_str   = bquote("(B) Waveform orientation: varying" ~ beta ~
                          "  (fixed" ~ omega ~ "=" ~ .(omega_fixed) ~ ")"),
    legend_ncol = 3L
  )

  if (has_alpha && n_alpha > 0) {
    p_alpha <- make_panel(
      full_df[full_df$sweep_param=="alpha",],
      alpha_levels, alpha_colors, alpha_labels,
      legend_name = expression(alpha ~ "(hours)"),
      x_lab       = "Total sample size N",
      title_str   = expression("(C) Location parameter " * alpha *
                                ":  active = flat  |  passive = varies with TOD"),
      legend_ncol = 3L,
      jitter_x    = 2.5   # separate overlapping lines in passive panels
    )
    combined <- p_omega / p_beta / p_alpha +
      patchwork::plot_layout(heights = c(1,1,1)) +
      patchwork::plot_annotation(
        title = "FMM K-harmonic LRT waveform robustness: shape, orientation, and acrophase",
        theme = ggplot2::theme(
          plot.title = ggplot2::element_text(hjust=0.5, face="bold", size=13))
      )
  } else {
    combined <- p_omega / p_beta +
      patchwork::plot_layout(heights = c(1,1)) +
      patchwork::plot_annotation(
        title = "FMM K-harmonic LRT waveform robustness",
        theme = ggplot2::theme(
          plot.title = ggplot2::element_text(hjust=0.5, face="bold", size=13))
      )
  }

  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, combined, width = width, height = height)
    message(sprintf("Saved: %s", output_file))
  } else {
    print(combined)
  }
  invisible(combined)
}


# =====================================================================
# Differential FMM plot (Figure 4b)
# =====================================================================

#' Plot FMM Waveform Robustness: Differential Endpoints
#'
#' @description
#' Produces a 1-row x 3-column figure showing DCP differential power (DR, DP, DM)
#' as a function of N for multiple omega values under FMM cosinor violation.
#'
#' The FMM model (Rueda, Rodríguez-Collado & Peddada, 2019, Sci Rep,
#' doi:10.1038/s41598-019-54569-1) parameterises waveform shape via
#' \code{omega}: \code{omega = 1} recovers pure cosinor; \code{omega = 0}
#' gives a flat (arrhythmic) signal; intermediate values produce increasingly
#' peaked waveforms.  \code{beta} controls orientation (peak location) and
#' does not affect power (beta-invariance of DCP amplitude testing).
#'
#' @param x          Data frame with columns \code{N, omega, DR, DP, DM,
#'                   SE_DR, SE_DP, SE_DM}. Produced by \code{fig4b_fmm_differential.R}
#'                   or any equivalent script calling \code{simCircadianDiffFMM()}.
#' @param nsims      Number of simulation replicates (for SE to SD conversion).
#' @param dataset_label Character label for the plot subtitle describing the
#'                   pilot dataset used (e.g. supplied by the caller from their
#'                   pilot data metadata). If \code{NULL}, an auto-generated
#'                   label based on the N and omega ranges is used.
#' @param output_file PDF path. \code{NULL} = current device.
#' @param width PDF width in inches (default 14).
#' @param height PDF height in inches (default 5).
#'
#' @return Invisibly returns the ggplot object.
#' @seealso \code{\link{simCircadianDiffFMM}}, \code{\link{plotFMMViolation}}
#' @export
plotFMMDifferential <- function(x,
                                 nsims         = 20L,
                                 dataset_label = NULL,
                                 output_file   = NULL,
                                 width         = 14,
                                 height        = 5) {

  for (pkg in c("ggplot2","tidyr"))
    if (!requireNamespace(pkg, quietly = TRUE))
      stop(sprintf("Package '%s' is required for plotFMMDifferential().", pkg))

  df <- x
  # Auto-generate dataset label if not provided
  if (is.null(dataset_label)) {
    N_max  <- max(df$N, na.rm = TRUE)
    N_min  <- min(df$N, na.rm = TRUE)
    dataset_label <- sprintf("N per group: %d-%d | omega: %s",
      N_min, N_max,
      paste(sort(unique(round(df$omega, 2))), collapse=", "))
  }
  # SE → SD
  df$SD_DR <- df$SE_DR * sqrt(nsims)
  df$SD_DP <- df$SE_DP * sqrt(nsims)
  df$SD_DM <- df$SE_DM * sqrt(nsims)

  # Reshape to long
  long_df <- tidyr::pivot_longer(df, cols = c(DR, DP, DM),
                                  names_to = "endpoint", values_to = "power")
  # Match on both N and omega so each row gets the SD from its own (N, omega) cell,
  # not just the first matching N (which was wrong when multiple omega levels share a N).
  row_key     <- paste(long_df$N, long_df$omega, sep = "_")
  df_key      <- paste(df$N,      df$omega,      sep = "_")
  match_idx   <- match(row_key, df_key)
  long_df$power_sd <- ifelse(long_df$endpoint=="DR", df$SD_DR[match_idx],
                      ifelse(long_df$endpoint=="DP", df$SD_DP[match_idx],
                                                     df$SD_DM[match_idx]))
  long_df$endpoint <- factor(long_df$endpoint,
    levels = c("DR","DP","DM"),
    labels = c("Differential Rhythmicity (DR)",
               "Differential Phase (DP)",
               "Differential Mesor (DM)"))

  omega_levels <- as.character(sort(unique(long_df$omega)))
  n_omega      <- length(omega_levels)
  omega_colors <- stats::setNames(
    grDevices::colorRampPalette(c("#08306B","#4292C6","#C6DBEF"))(n_omega),
    omega_levels)
  omega_labels <- stats::setNames(paste0("omega==", omega_levels), omega_levels)
  long_df$omega_fac <- factor(as.character(long_df$omega), levels = omega_levels)

  N_breaks <- sort(unique(long_df$N))
  x_breaks <- N_breaks  # show all N values

  p <- ggplot2::ggplot(long_df, ggplot2::aes(
      x = N, y = 100 * power, colour = omega_fac, group = omega_fac)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = 100*(power - power_sd),
                   ymax = 100*(power + power_sd)),
      width = diff(range(N_breaks)) * 0.015, linewidth = 0.45, alpha = 0.8) +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::geom_hline(yintercept = 80, linetype = "dashed",
                        colour = "grey50", linewidth = 0.4) +
    ggplot2::facet_wrap(~ endpoint, nrow = 1) +
    ggplot2::scale_colour_manual(values = omega_colors,
                                  labels = parse(text = unname(omega_labels)),
                                  name = expression(omega)) +
    ggplot2::scale_x_continuous(breaks = x_breaks) +
    ggplot2::scale_y_continuous(limits = c(0,100), breaks = seq(0,100,20)) +
    ggplot2::labs(
      x        = "Total sample size N per group",
      y        = "Power (%)",
      title    = "Power analysis under cosinor violation: differential endpoints",
      subtitle = dataset_label,
      caption  = sprintf("+/-1 SD across %d replicates", nsims)
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "grey92"),
      strip.text       = ggplot2::element_text(size = 10, face = "bold"),
      legend.position  = "bottom",
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
      plot.title       = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.subtitle    = ggplot2::element_text(hjust = 0.5),
      plot.caption     = ggplot2::element_text(hjust = 0, size = 8, colour = "grey40")
    )

  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, p, width = width, height = height)
    message(sprintf("Saved: %s", output_file))
  } else {
    print(p)
  }
  invisible(p)
}
