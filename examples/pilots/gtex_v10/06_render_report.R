# Render the paired v10 sample-size comparison after validation.
# This report is separate from the manuscript and does not rebuild its figures.

# 1. Read the fixed pilots and completed simulation results ----------------
stage <- Sys.getenv("SCP_GTEX_V10_STAGE", "output/gtex_v10_chiral_protocol")
results <- file.path(stage, "simulation", "results_exactF")
report <- file.path(stage, "report")
dir.create(report, recursive = TRUE, showWarnings = FALSE)
required <- c("N80_summary.csv", "power_summary.csv", "run_settings.json")
stopifnot(all(file.exists(file.path(results, required))))
crossings <- utils::read.csv(file.path(results, "N80_summary.csv"))
power <- utils::read.csv(file.path(results, "power_summary.csv"))
pilots <- readRDS(file.path(stage, "private", "report_pilots.rds"))
if (!is.null(pilots$pilots)) pilots <- pilots$pilots
pilots <- Filter(function(x) identical(x$cohort, "matched"), pilots)
settings <- jsonlite::read_json(file.path(results, "run_settings.json"), simplifyVector = TRUE)
coverage <- jsonlite::read_json(file.path(stage, "qc", "donor_coverage.json"), simplifyVector = TRUE)
stopifnot(nrow(crossings) == 20L, length(unique(crossings$tissue)) == 10L,
          all(crossings$time_basis %in% c("DIP", "TOD")),
          all(table(crossings$tissue) == 2L), settings$G == 5000L,
          settings$nsims_initial >= 200L, settings$nsims_final >= 1000L)

# 2. Match every result to its saved pilot calibration ---------------------
pilot_table <- do.call(rbind, lapply(pilots, function(x) {
  data.frame(tissue = x$tissue, time_basis = x$time_basis,
             n_pilot = x$n_pilot, genes_fitted = x$genes_fitted,
             n_selected = x$n_selected, median_r = x$medianr,
             prop_rhythmic = x$prop, stringsAsFactors = FALSE)
}))
key <- function(x) paste(x$tissue, x$time_basis, sep = "__")
stopifnot(!anyDuplicated(key(pilot_table)), !anyDuplicated(key(crossings)))
index <- match(key(crossings), key(pilot_table))
stopifnot(!anyNA(index), all(crossings$n_pilot == pilot_table$n_pilot[index]))
comparison <- cbind(crossings, pilot_table[index,
  c("genes_fitted", "n_selected", "median_r", "prop_rhythmic"), drop = FALSE])
comparison$input_time <- ifelse(comparison$time_basis == "DIP",
  "CHIRAL-estimated time (DIP)", "Season-adjusted time of death (TOD)")
comparison$gtex_release <- "v10"
comparison$interval_scope <- "Monte Carlo uncertainty conditional on fixed pilots"
for (tissue in unique(comparison$tissue)) {
  rows <- pilot_table[pilot_table$tissue == tissue, ]
  stopifnot(length(unique(rows$n_pilot)) == 1L,
            length(unique(rows$genes_fitted)) == 1L)
}

# 3. Keep one current set of report outputs -------------------------------
names_out <- c("N80_comparison_table.csv", "N80_ten_tissue_summary.csv",
               "n80_ten_tissues.png", "n80_ten_tissues.pdf",
               "power_curves_ten_tissues.png", "power_curves_ten_tissues.pdf",
               "N80_comparison_report.html")
if (any(file.exists(file.path(report, names_out))))
  stop("Archive and verify the current report outputs before replacing them.")
utils::write.csv(comparison, file.path(report, names_out[1]), row.names = FALSE)
wide_columns <- c("tissue", "n_pilot", "median_r", "prop_rhythmic",
                  "N80_interpolated", "N80_mc95_low", "N80_mc95_high", "status")
wide <- merge(comparison[comparison$time_basis == "DIP", wide_columns],
              comparison[comparison$time_basis == "TOD", wide_columns],
              by = "tissue", suffixes = c("_DIP", "_TOD"), sort = FALSE)
utils::write.csv(wide, file.path(report, names_out[2]), row.names = FALSE)

# 4. Draw N80 and the full power curves with standard plotting tools -------
library(ggplot2)
colors <- c(DIP = "#D55E00", TOD = "#0072B2")
tissue_order <- unique(comparison$tissue)
labels <- function(x) {
  lookup <- c(Adrenal = "Adrenal gland", Liver = "Liver", Muscle = "Skeletal muscle",
    Adipose_Subcutaneous = "Subcutaneous adipose", Artery_Tibial = "Tibial artery",
    Brain_Cortex = "Cerebral cortex", Thyroid = "Thyroid", Lung = "Lung",
    Pancreas = "Pancreas", Heart_Left_Ventricle = "Heart, left ventricle")
  result <- unname(lookup[x])
  result[is.na(result)] <- gsub("_", " ", x[is.na(result)], fixed = TRUE)
  result
}
comparison$tissue_label <- factor(labels(comparison$tissue),
                                   levels = rev(labels(tissue_order)))
comparison$plot_N80 <- ifelse(is.finite(comparison$N80_interpolated),
  comparison$N80_interpolated, ifelse(is.finite(comparison$N80_high),
    comparison$N80_high, settings$maximum_N))
position <- position_dodge(width = 0.5)
n80_plot <- ggplot(comparison, aes(x = plot_N80, y = tissue_label,
                                    color = time_basis, group = time_basis)) +
  geom_errorbar(aes(xmin = N80_mc95_low, xmax = N80_mc95_high),
                orientation = "y", width = 0.2, position = position, na.rm = TRUE) +
  geom_point(aes(shape = status), size = 2.7, position = position) +
  scale_x_log10(breaks = c(25, 50, 100, 200, 500, 1000, 2000, 5000, 10000)) +
  scale_color_manual(values = colors, name = "Input time") +
  labs(title = "GTEx v10: sample size for 80% detection power",
       x = "N80 (samples, logarithmic scale)", y = NULL,
       caption = paste("Intervals show Monte Carlo uncertainty with fixed pilots.",
                       "Pilot and inferred-time uncertainty are excluded.", sep = "\n")) +
  theme_bw(base_size = 12) + theme(legend.position = "top", panel.grid.minor = element_blank())
if (all(comparison$status == "bracketed")) n80_plot <- n80_plot + guides(shape = "none")
power$tissue_label <- factor(labels(power$tissue), levels = labels(tissue_order))
curve_plot <- ggplot(power, aes(x = N, y = power, color = time_basis)) +
  geom_hline(yintercept = 0.8, linetype = 3, color = "grey50") +
  geom_ribbon(aes(ymin = power_mc95_low, ymax = power_mc95_high,
                  fill = time_basis), alpha = 0.13, color = NA) +
  geom_line(linewidth = 0.65) + geom_point(size = 1.0) +
  facet_wrap(~tissue_label, ncol = 2) + scale_x_log10() +
  scale_color_manual(values = colors, name = "Input time") +
  scale_fill_manual(values = colors, name = "Input time") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = "GTEx v10: paired input-time examples",
       x = "Samples (logarithmic scale)", y = "Detection power",
       caption = "K = 1. BH at 0.05. Each tissue uses matched donors and genes.") +
  theme_bw(base_size = 11) + theme(legend.position = "top", panel.grid.minor = element_blank())
png_device <- if (identical(Sys.info()[["sysname"]], "Darwin")) "quartz" else "cairo"
ggsave(file.path(report, names_out[3]), n80_plot, width = 11, height = 7.5,
       dpi = 170, device = grDevices::png, type = png_device)
ggsave(file.path(report, names_out[4]), n80_plot, width = 11, height = 7.5,
       device = grDevices::pdf)
ggsave(file.path(report, names_out[5]), curve_plot, width = 11, height = 15,
       dpi = 170, device = grDevices::png, type = png_device)
ggsave(file.path(report, names_out[6]), curve_plot, width = 11, height = 15,
       device = grDevices::pdf)

# 5. Embed the figures and tables in a standalone report ------------------
escape_html <- function(x) {
  x <- gsub("&", "&amp;", as.character(x), fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  gsub(">", "&gt;", x, fixed = TRUE)
}
shown <- comparison[, c("tissue", "time_basis", "n_pilot", "genes_fitted", "n_selected", "median_r",
                        "prop_rhythmic", "N80_interpolated",
                        "N80_mc95_low", "N80_mc95_high", "status")]
for (column in names(shown)[vapply(shown, is.numeric, logical(1))])
  shown[[column]] <- ifelse(is.finite(shown[[column]]),
                             format(round(shown[[column]], 3), trim = TRUE), "Unavailable")
table_rows <- apply(shown, 1, function(row)
  paste0("<tr><td>", paste(escape_html(row), collapse = "</td><td>"), "</td></tr>"))
image_tag <- function(name, alt) paste0('<img alt="', alt,
  '" src="data:image/png;base64,', base64enc::base64encode(file.path(report, name)), '">')
diagnostic_table <- function(name, columns) {
  path <- file.path(stage, "qc", name)
  stopifnot(file.exists(path))
  x <- utils::read.csv(path, stringsAsFactors = FALSE)[, columns, drop = FALSE]
  for (column in names(x)[vapply(x, is.numeric, logical(1))])
    x[[column]] <- ifelse(is.finite(x[[column]]), format(round(x[[column]], 3), trim = TRUE), "Unavailable")
  rows <- apply(x, 1, function(row) paste0("<tr><td>",
    paste(escape_html(row), collapse = "</td><td>"), "</td></tr>"))
  c('<div class="table"><table><thead><tr><th>',
    paste(escape_html(names(x)), collapse = '</th><th>'),
    '</th></tr></thead><tbody>', rows, '</tbody></table></div>')
}
html <- c('<!doctype html><html lang="en"><meta charset="utf-8">',
  '<title>GTEx v10 input-time comparison</title>',
  '<style>body{max-width:1120px;margin:32px auto;padding:0 20px;font:16px/1.55 system-ui;color:#17212b}img{width:100%;height:auto}table{border-collapse:collapse;font-size:13px;width:100%}th,td{padding:7px;border-bottom:1px solid #ddd;text-align:left}thead{background:#edf2f7}.table{overflow:auto}p{max-width:900px}</style>',
  '<h1>GTEx v10 input-time examples</h1>',
  '<p>SCP pairs expression with donor time information to estimate power under a chosen model and sampling design. This report compares CHIRAL-estimated time (DIP) with season-adjusted time of death (TOD) in ten tissues. Each comparison uses the same eligible expression matrix, donors, and genes.</p>',
  '<p>The median signal-to-noise ratio describes up to 300 selected genes passing the nominal pilot screen at p &lt; 0.01. Each input time selects its own gene set. The table shows the selected counts. The rhythmic fraction is the fraction passing that screening rule, rather than a measured biological prevalence.</p>',
  paste0('<p>The published CHIRAL procedure produced DIP for ', coverage$n_with_DIP, ' of ',
    coverage$n_source_donors, ' source donors. ', coverage$native_converged_tissues, ' of ',
    coverage$n_tissues, ' tissue fits reached native convergence. The pilot tables report the matched cohorts.</p>'),
  paste0('<p>Detection uses a single-harmonic cosinor model and BH at 0.05 with ',
    settings$G, ' simulated genes, ', settings$nsims_initial,
    ' initial experiments per sample size, and at least ', settings$nsims_final,
    ' experiments near the N80 bracket. The calculation uses the exact conditional Gaussian F distribution and retains a full-expression validation engine.</p>'),
  '<p>N80 is the sample size at which the mean proportion of simulated true rhythmic genes detected reaches 80%. The displayed intervals describe Monte Carlo uncertainty conditional on fixed pilots. They exclude pilot, inferred-time, cohort, and grid-interpolation uncertainty. Stronger fits or lower N80 do not independently validate the inferred times.</p>',
  image_tag(names_out[3], 'Paired N80 estimates for ten v10 tissues'),
  '<h2>Paired pilot and N80 summary</h2><div class="table"><table><thead><tr><th>',
  paste(escape_html(names(shown)), collapse = '</th><th>'),
  '</th></tr></thead><tbody>', table_rows, '</tbody></table></div>',
  image_tag(names_out[5], 'Detection power curves for the paired v10 examples'),
  '<h2>Clock-gene exclusion</h2>',
  '<p>These fits exclude the 12 reference clock genes from the tested gene universe. They describe sensitivity to using those genes in the pilot summary.</p>',
  diagnostic_table("clock_genes_excluded.csv", c("tissue", "time_basis", "n", "ngenes", "n_selected", "medianr", "prop_rhythmic", "n_BH05")),
  '<h2>Target-tissue exclusion</h2>',
  '<p>The target tissue is removed before rebuilding donor consensus with the published defaults. The two fits use the same retained donors and genes. Coverage loss is reported. These checks do not independently validate phase accuracy.</p>',
  diagnostic_table("leave_target_out_summary.csv", c("tissue", "comparison", "n_original_paired", "n_lost", "n", "ngenes", "n_selected", "medianr", "prop_rhythmic", "n_BH05")),
  '<p>Expression preprocessing follows the published CHIRAL covariate adjustment. The manuscript inputs and figures were not regenerated.</p>',
  '<p>CHIRAL: <a href="https://doi.org/10.1126/science.add0846">Talamanca, Gobet and Naef (2023)</a>. GTEx project: <a href="https://doi.org/10.1126/science.aaz1776">GTEx Consortium (2020)</a>. Expression release: v10.</p>',
  '</html>')
writeLines(html, file.path(report, names_out[7]))
jsonlite::write_json(setNames(lapply(file.path(report, names_out), function(path)
  digest::digest(file = path, algo = "sha256")), names_out),
  file.path(report, "artifact_sha256.json"), auto_unbox = TRUE, pretty = TRUE)
cat("Report written to", normalizePath(report), "\n")
