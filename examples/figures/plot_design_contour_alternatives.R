#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript examples/plot_design_contour_alternatives.R <design_contour_grid.rds>")
}

grid_path <- args[1]
grid <- readRDS(grid_path)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("ggplot2 is required for plotting")
}

suppressPackageStartupMessages(library(ggplot2))

out_dir <- file.path(dirname(grid_path), "figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

df <- grid[!is.na(grid$B) & grid$design_type == "active", ]
if (nrow(df) == 0) stop("No active-design rows found in grid.")

df$n <- as.integer(df$n)
df$B <- as.integer(df$B)

plot_power_vs_B <- function(df, test = "DR") {
  edr_col <- paste0("EDR_", test)
  if (!edr_col %in% names(df)) return(NULL)
  ggplot(df, aes(x = B, y = .data[[edr_col]], color = factor(n), group = n)) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.5) +
    scale_x_continuous(breaks = sort(unique(df$B))) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(x = "Number of time bins (B)", y = "Targeted power",
         color = "n per group",
         title = sprintf("%s power vs B (fixed n)", test))
}

plot_power_vs_n <- function(df, test = "DR") {
  edr_col <- paste0("EDR_", test)
  if (!edr_col %in% names(df)) return(NULL)
  ggplot(df, aes(x = n, y = .data[[edr_col]], color = factor(B), group = B)) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.5) +
    scale_x_continuous(breaks = sort(unique(df$n))) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(x = "Sample size per group (n)", y = "Targeted power",
         color = "B bins",
         title = sprintf("%s power vs n (fixed B)", test))
}

plot_delta_heatmap <- function(df, test = "DR", baseline_B = 4) {
  edr_col <- paste0("EDR_", test)
  if (!edr_col %in% names(df)) return(NULL)
  base <- df[df$B == baseline_B, c("n", edr_col)]
  names(base)[2] <- "EDR_base"
  df2 <- merge(df, base, by = "n", all.x = TRUE)
  df2$delta <- df2[[edr_col]] - df2$EDR_base
  ggplot(df2, aes(x = n, y = B, fill = delta)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_x_continuous(breaks = sort(unique(df2$n))) +
    scale_y_continuous(breaks = sort(unique(df2$B))) +
    scale_fill_gradient2(low = "#2c7bb6", mid = "white", high = "#d7191c",
                         midpoint = 0, labels = scales::percent_format(accuracy = 0.1)) +
    labs(x = "Sample size per group (n)", y = "Number of time bins (B)",
         fill = sprintf("Targeted power - (B=%d)", baseline_B),
         title = sprintf("%s: gain from increasing B (baseline B=%d)", test, baseline_B))
}

plot_frontier <- function(df, test = "DR") {
  edr_col <- paste0("EDR_", test)
  if (!edr_col %in% names(df)) return(NULL)
  best <- do.call(rbind, lapply(split(df, df$n), function(d) {
    d <- d[order(d[[edr_col]], d$B), ]
    d[nrow(d), , drop = FALSE]
  }))
  ggplot(best, aes(x = n, y = B, color = .data[[edr_col]])) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2) +
    scale_x_continuous(breaks = sort(unique(df$n))) +
    scale_y_continuous(breaks = sort(unique(df$B))) +
    scale_color_gradient(low = "#2c7bb6", high = "#d7191c",
                         labels = scales::percent_format(accuracy = 1)) +
    labs(x = "Sample size per group (n)", y = "Optimal B (max targeted power)",
         color = "Targeted power",
         title = sprintf("%s: optimal B by n (frontier)", test))
}

plot_neff <- function(df, test = "DR") {
  edr_col <- paste0("EDR_", test)
  if (!edr_col %in% names(df)) return(NULL)
  df$n_eff <- df$n * pmax(0.0, 0.5 - df$d_sd)  # approx d_eff for active designs
  ggplot(df, aes(x = n_eff, y = .data[[edr_col]], color = factor(B))) +
    geom_point(size = 2) +
    geom_smooth(se = FALSE, linewidth = 0.6) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(x = "Approx. effective n (n × (0.5 − d_sd))",
         y = "Targeted power", color = "B bins",
         title = sprintf("%s: power vs effective n (active designs)", test))
}

plot_set <- function(test) {
  p1 <- plot_power_vs_B(df, test)
  p2 <- plot_power_vs_n(df, test)
  p3 <- plot_delta_heatmap(df, test)
  p4 <- plot_frontier(df, test)
  p5 <- plot_neff(df, test)

  if (!is.null(p1)) ggsave(file.path(out_dir, sprintf("design_alt_%s_vs_B.pdf", test)),
                           p1, width = 6.5, height = 4.2)
  if (!is.null(p2)) ggsave(file.path(out_dir, sprintf("design_alt_%s_vs_n.pdf", test)),
                           p2, width = 6.5, height = 4.2)
  if (!is.null(p3)) ggsave(file.path(out_dir, sprintf("design_alt_%s_deltaB.pdf", test)),
                           p3, width = 6.5, height = 4.2)
  if (!is.null(p4)) ggsave(file.path(out_dir, sprintf("design_alt_%s_frontier.pdf", test)),
                           p4, width = 6.5, height = 4.2)
  if (!is.null(p5)) ggsave(file.path(out_dir, sprintf("design_alt_%s_neff.pdf", test)),
                           p5, width = 6.5, height = 4.2)
}

plot_set("DR")
plot_set("DP")

cat("Wrote alternative plots to:\n", out_dir, "\n")
