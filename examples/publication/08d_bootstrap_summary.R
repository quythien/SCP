#' =======================================================================
#' 08d_bootstrap_summary.R — Two-Stage vs Bootstrap: Summary Figure
#' =======================================================================
#' Loads RDS outputs from 08a/08b/08c and generates:
#'   - s4_ci_width_summary.pdf  (Panel A: all power curves; Panel B: CI width bar)
#'   - comparison_summary.txt
#'
#' Run AFTER 08a, 08b, 08c have all completed:
#'   POWERSIM_ROOT=$ROOT RUN_TAG=20260401 Rscript examples/publication/08d_bootstrap_summary.R

POWERSIM_ROOT <- {
  env <- Sys.getenv("POWERSIM_ROOT", unset = "")
  if (nzchar(env)) env else
    "/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim"
}
setwd(POWERSIM_ROOT)

RUN_TAG <- Sys.getenv("RUN_TAG", format(Sys.Date(), "%Y%m%d"))
out_dir <- file.path("output", paste0("08_two_stage_vs_bootstrap_", RUN_TAG))

if (!dir.exists(out_dir)) {
  stop(sprintf("Output dir not found: %s\nDid you set RUN_TAG correctly?", out_dir))
}

source_dir <- file.path(getwd(), "code")
old_wd <- setwd(source_dir); source("setup.R"); setwd(old_wd)

cat(sprintf("Loading results from: %s/\n\n", out_dir))

# -----------------------------------------------------------------------
# Load per-dataset RDS files
# -----------------------------------------------------------------------
rds_map <- list(
  baboon = file.path(out_dir, "s1_baboon_comparison.rds"),
  d1d2   = file.path(out_dir, "s2_d1d2_comparison.rds"),
  seney  = file.path(out_dir, "s3_seney_comparison.rds")
)

# Metadata needed for the summary plot (matches what .run_comparison returns)
meta_map <- list(
  baboon = list(label = "Baboon LUN",   color = "darkorange"),
  d1d2   = list(label = "D1D2 D1",      color = "forestgreen"),
  seney  = list(label = "Seney CTL",    color = "firebrick")
)

results_all   <- list()
summary_lines <- c("Two-Stage vs Bootstrap Real-Data Comparison",
                   sprintf("RUN_TAG: %s", RUN_TAG), "")

for (nm in names(rds_map)) {
  f <- rds_map[[nm]]
  if (!file.exists(f)) {
    warning(sprintf("Missing RDS: %s — skipping %s", f, nm))
    next
  }
  dat <- readRDS(f)
  comp_df   <- dat$comparison$comparison
  ci_widths <- comp_df$boot_ci_hi - comp_df$boot_ci_lo

  results_all[[nm]] <- list(
    label       = meta_map[[nm]]$label,
    color       = meta_map[[nm]]$color,
    n_pilot     = ncol(dat$two_stage$pilot_data %||% matrix(nrow=1, ncol=1)),
    N_grid      = comp_df$n,
    comparison  = dat$comparison,
    ci_widths   = ci_widths,
    ts_result   = dat$two_stage,
    boot_result = dat$boot
  )

  summary_lines <- c(summary_lines,
    sprintf("%s: two-stage n80=%s  boot_median n80=%s  mean_CI=%.0f pp",
            toupper(nm),
            ifelse(is.na(dat$comparison$n80_two_stage), ">max", dat$comparison$n80_two_stage),
            ifelse(is.na(dat$comparison$n80_boot_median), ">max", round(dat$comparison$n80_boot_median)),
            100 * mean(ci_widths, na.rm = TRUE)))
  cat(sprintf("  Loaded: %s\n", nm))
}

if (length(results_all) == 0) stop("No datasets loaded.")

# =======================================================================
# SECTION 4: SUMMARY — CI width vs pilot size
# =======================================================================
cat("\n====================================================================\n")
cat("SECTION 4: CI width summary across pilot sizes\n")
cat("====================================================================\n\n")

fig_s4 <- file.path(out_dir, "s4_ci_width_summary.pdf")
pdf(fig_s4, width = 12, height = 5)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3.5, 1.5))

# Panel A: all power curves overlaid
all_N <- sort(unique(unlist(lapply(results_all, `[[`, "N_grid"))))
plot(NA, xlim = range(all_N), ylim = c(0, 1),
     xlab = "N per group", ylab = "Power",
     main = "Two-stage (solid) vs Bootstrap (dashed +/- CI)\nby dataset", las = 1)
abline(h = 0.80, lty = 3, col = "gray50")

for (nm in names(results_all)) {
  sc  <- results_all[[nm]]
  cmp <- sc$comparison$comparison
  col <- sc$color

  polygon(c(cmp$n, rev(cmp$n)),
          c(cmp$boot_ci_lo, rev(cmp$boot_ci_hi)),
          col = adjustcolor(col, 0.12), border = NA)
  lines(cmp$n, cmp$boot_power_mean,  col = col, lwd = 2, lty = 2)
  lines(cmp$n, cmp$two_stage_power,  col = col, lwd = 2.5, lty = 1)
  points(cmp$n, cmp$two_stage_power, col = col, pch = 16, cex = 0.8)
  arrows(cmp$n, cmp$two_stage_power - cmp$two_stage_se,
         cmp$n, cmp$two_stage_power + cmp$two_stage_se,
         length = 0.05, angle = 90, code = 3, col = col, lwd = 1.2)
}

legend("bottomright",
       legend = sapply(results_all, `[[`, "label"),
       col    = sapply(results_all, `[[`, "color"),
       lwd = 2, bty = "n", cex = 0.78)

# Panel B: mean CI width bar chart
mean_widths <- sapply(results_all, function(sc) 100 * mean(sc$ci_widths, na.rm = TRUE))
pilot_ns    <- sapply(results_all, `[[`, "n_pilot")
cols_bar    <- sapply(results_all, `[[`, "color")
labels_bar  <- sapply(results_all, `[[`, "label")

ord <- order(pilot_ns)
bp  <- barplot(mean_widths[ord],
               names.arg = labels_bar[ord],
               col       = cols_bar[ord],
               ylim      = c(0, max(mean_widths) * 1.25),
               ylab      = "Mean bootstrap CI width (pp)",
               main      = "CI width vs pilot size\n(wider = more uncertainty)",
               las = 2, cex.names = 0.75)
text(bp, mean_widths[ord] + max(mean_widths) * 0.04,
     sprintf("n=%d", pilot_ns[ord]), cex = 0.85, col = "gray20")

dev.off()
cat(sprintf("Figure: %s\n", fig_s4))

summary_lines <- c(summary_lines, "",
  "Key interpretation:",
  "  - Two-stage n80 and bootstrap median n80 should be close",
  "  - Bootstrap CI width decreases as pilot size increases",
  "  - Small pilots (Baboon n=12) give wide CI — two-stage has false precision there"
)

writeLines(summary_lines, file.path(out_dir, "comparison_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")
cat(sprintf("\nOutput: %s/\nDone.\n", out_dir))
