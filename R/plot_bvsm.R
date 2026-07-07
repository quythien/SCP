#' =======================================================================
#' plot_bvsm.R, B vs m Trade-off Plot (Figure 3)
#' =======================================================================
#'
#' Plots DCP power vs N for multiple B values, demonstrating B-invariance.
#' Input: SCPSingleResult from runSingleCohortGrid() or a tidy data.frame
#' with columns N, B, alpha2, power, power_se.

# =====================================================================
# Main plot function
# =====================================================================

#' Plot B vs m Trade-off Power Grid
#'
#' @description
#' Produces a faceted power-vs-N figure for each dataset in the result,
#' with one line per B value. SD bars (+/-1 SD across simulation replicates)
#' are drawn at each N. B-invariance appears as overlapping lines.
#'
#' @param x         An \code{SCPSingleResult} from \code{runSingleCohortGrid()},
#'                  or a tidy \code{data.frame} with columns
#'                  \code{N, B, power, power_se, dataset} (and optionally
#'                  \code{alpha2}).
#' @param nsims     Number of simulation replicates used (to convert SE -> SD).
#'                  Default 30.
#' @param alpha2    Which \code{alpha2} value to plot (default 0 = pure cosinor).
#' @param fdr_power  Reference power level for dashed line (default 0.80).
#' @param dataset_labels Named character vector mapping dataset codes to display
#'                  labels. If NULL, uses dataset column values directly.
#' @param output_file Path to save PDF. NULL = current graphics device.
#' @param width     PDF width in inches (default 14).
#' @param height    PDF height in inches (default 5).
#'
#' @return Invisibly returns the ggplot object.
#' @export
plotBvsMPower <- function(x,
                           nsims         = 30L,
                           alpha2        = 0,
                           fdr_power     = 0.80,
                           dataset_labels = NULL,
                           output_file   = NULL,
                           width         = 14,
                           height        = 5) {

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required for plotBvsMPower().")

  # -- Extract tidy data frame -------------------------------------
  if (inherits(x, "SCPSingleResult")) {
    df <- x$power_df
  } else if (is.list(x) && !is.data.frame(x) &&
             all(vapply(x, inherits, FALSE, "SCPSingleResult"))) {
    df <- do.call(rbind, lapply(x, `[[`, "power_df"))
  } else if (is.data.frame(x)) {
    df <- x
  } else {
    stop("x must be an SCPSingleResult, a list of SCPSingleResult, or a data.frame with columns N, B, power, power_se")
  }

  required_cols <- c("N", "B", "power", "power_se")
  missing <- setdiff(required_cols, names(df))
  if (length(missing) > 0)
    stop("Missing columns: ", paste(missing, collapse = ", "))

  # Filter alpha2
  if ("alpha2" %in% names(df))
    df <- df[df$alpha2 == alpha2, , drop = FALSE]

  # Convert SE -> SD
  df$power_sd <- df$power_se * sqrt(nsims)

  # Dataset labels, auto-generate from dataset column if not provided
  if (is.null(dataset_labels) && "dataset" %in% names(df)) {
    ds_vals <- sort(unique(as.character(df$dataset)))
    dataset_labels <- stats::setNames(ds_vals, ds_vals)
  }
  if (!is.null(dataset_labels) && "dataset" %in% names(df)) {
    df$dataset_label <- factor(
      dataset_labels[as.character(df$dataset)],
      levels = dataset_labels
    )
  } else if ("dataset" %in% names(df)) {
    df$dataset_label <- factor(df$dataset)
  } else {
    df$dataset_label <- "Dataset"
  }

  B_vals   <- sort(unique(df$B))
  b_colors <- setNames(
    c("#E41A1C","#FF7F00","#4DAF4A","#377EB8","#984EA3",
      "#A65628","#F781BF","#999999")[seq_len(length(B_vals))],
    as.character(B_vals)
  )
  b_labels <- setNames(paste0("B=", B_vals), as.character(B_vals))
  df$B_fac <- factor(df$B, levels = B_vals)

  N_breaks <- sort(unique(df$N))

  p <- ggplot2::ggplot(df, ggplot2::aes(
      x = N, y = 100 * power, colour = B_fac, group = B_fac)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = 100*(power - power_sd),
                   ymax = 100*(power + power_sd)),
      width = diff(range(N_breaks)) * 0.012,
      linewidth = 0.5) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_hline(yintercept = 100 * fdr_power,
                        linetype = "dashed", colour = "grey50", linewidth = 0.4) +
    ggplot2::facet_wrap(~ dataset_label, nrow = 1) +
    ggplot2::scale_colour_manual(values = b_colors, labels = b_labels,
                                  name = "Time bins (B)") +
    ggplot2::scale_x_continuous(breaks = N_breaks) +
    ggplot2::scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    ggplot2::labs(
      x       = "Total sample size N",
      y       = "Power (%)",
      title   = "B vs m trade-off: DCP power is B-invariant under cosinor truth",
      caption = sprintf("+/-1 SD across %d simulation replicates", nsims)
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "grey92"),
      strip.text       = ggplot2::element_text(size = 10),
      legend.position  = "bottom",
      panel.grid.minor = ggplot2::element_blank(),
      plot.title       = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.caption     = ggplot2::element_text(hjust = 0, size = 8, colour = "grey40"),
      axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1, size = 8)
    )

  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, p, width = width, height = height)
    message(sprintf("Saved: %s", output_file))
  } else {
    print(p)
  }
  invisible(p)
}
