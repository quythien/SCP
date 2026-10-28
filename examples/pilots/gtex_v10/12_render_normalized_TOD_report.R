# Render the normalized-expression GTEx v10 TOD report from saved results.

# 1. Read the new normalized TOD results and staged pilot catalogue -------------
stopifnot(requireNamespace("ggplot2", quietly = TRUE),
          requireNamespace("digest", quietly = TRUE),
          requireNamespace("jsonlite", quietly = TRUE),
          requireNamespace("base64enc", quietly = TRUE))
stage <- Sys.getenv("SCP_GTEX_TOD_STAGE", "output/gtex_v10_normalized_TOD")
publication <- Sys.getenv("SCP_GTEX_TOD_PUBLICATION", stage)
report <- Sys.getenv("SCP_GTEX_TOD_REPORT",
                     file.path(stage, "report"))
results <- file.path(stage, "simulation", "results_exactF")
comparison_file <- file.path(stage, "report", "N80_comparison_table.csv")
signal_file <- file.path(publication, "summary", "gtex_v10_signal_summary.csv")
pilot_root <- file.path(publication, "inst", "extdata", "pilots")
manifest_file <- file.path(pilot_root, "manifest.csv")
source_files <- c(comparison = comparison_file, signal = signal_file,
  manifest = manifest_file, N80 = file.path(results, "N80_summary.csv"),
  power = file.path(results, "power_summary.csv"),
  settings = file.path(results, "run_settings.json"))
stopifnot(all(file.exists(source_files)))
comparison_raw <- utils::read.csv(comparison_file, colClasses = "character",
                                   check.names = FALSE)
comparison_raw <- comparison_raw[comparison_raw$time_basis == "TOD", , drop = FALSE]
power_raw <- utils::read.csv(source_files[["power"]], colClasses = "character",
                            check.names = FALSE)
power_raw <- power_raw[power_raw$time_basis == "TOD", , drop = FALSE]
signal_raw <- utils::read.csv(signal_file, colClasses = "character",
                             check.names = FALSE)
rownames(comparison_raw) <- NULL
rownames(power_raw) <- NULL
comparison <- utils::type.convert(comparison_raw, as.is = TRUE)
power <- utils::type.convert(power_raw, as.is = TRUE)
signal <- utils::type.convert(signal_raw, as.is = TRUE)
manifest <- utils::read.csv(manifest_file, stringsAsFactors = FALSE)
settings <- jsonlite::read_json(source_files[["settings"]], simplifyVector = TRUE)
stopifnot(nrow(comparison) == 10L, length(unique(comparison$tissue)) == 10L,
  !anyDuplicated(comparison$case_id), all(comparison$time_basis == "TOD"),
  all(comparison$gtex_release == "v10"), all(comparison$status == "bracketed"),
  nrow(signal) == 46L, !anyDuplicated(signal$tissue),
  all(signal$time_basis == "TOD"), all(signal$gtex_release == "v10"),
  settings$G == 5000L, settings$nsims_initial >= 200L,
  settings$nsims_final >= 1000L, settings$BH_threshold == 0.05)

# 2. Match the summaries to the staged normalized TOD pilot objects ----------------
near <- function(x, y) isTRUE(all.equal(as.numeric(x), as.numeric(y),
                                       tolerance = 1e-12, check.attributes = FALSE))
gtex <- manifest[manifest$dataset == "GTEx", , drop = FALSE]
stopifnot(nrow(gtex) == 46L, all(gtex$gtex_release == "v10"),
          all(gtex$time_basis == "TOD"), !anyDuplicated(gtex$tissue),
          setequal(gtex$tissue, signal$tissue))
validation <- list()
for (i in seq_len(nrow(signal))) {
  row <- signal[i, , drop = FALSE]
  entry <- gtex[match(row$tissue, gtex$tissue), , drop = FALSE]
  path <- file.path(pilot_root, entry$file)
  pilot <- readRDS(path)
  ratio <- pilot$amplitude / pilot$sigma_rhythmic
  qvalue <- pilot$explorer_fit$qvalue
  checks <- c(n = near(pilot$n_samples, row$n),
    genes = near(pilot$ngenes, row$ngenes),
    selected = near(length(ratio), row$n_selected),
    median_r = near(stats::median(ratio), row$r_median_top300),
    prop_rhythmic = near(pilot$prop_rhythmic, row$prop_rhythmic),
    n_BH05 = near(sum(is.finite(qvalue) & qvalue < 0.05), row$n_BH05),
    time_basis = identical(pilot$time_basis, "TOD"),
    release = identical(pilot$gtex_release, "v10"),
    nominal_screen = near(pilot$alpha_pilot, 0.01),
    selection_cap = near(pilot$pilot_top_k, 300L))
  stopifnot(all(checks), length(qvalue) == row$ngenes)
  validation[[i]] <- data.frame(tissue = row$tissue, n = row$n,
    ngenes = row$ngenes, n_selected = length(ratio),
    median_r = stats::median(ratio), prop_rhythmic = pilot$prop_rhythmic,
    n_BH05 = row$n_BH05, checks_passed = sum(checks),
    pilot_sha256 = digest::digest(file = path, algo = "sha256"),
    stringsAsFactors = FALSE)
}
validation <- do.call(rbind, validation)
canonical_tissue <- ifelse(comparison$tissue == "Muscle", "Muscle_Skeletal",
                           comparison$tissue)
index <- match(canonical_tissue, signal$tissue)
stopifnot(!anyNA(index), near(comparison$n_pilot, signal$n[index]),
  near(comparison$genes_fitted, signal$ngenes[index]),
  near(comparison$n_selected, signal$n_selected[index]),
  near(comparison$median_r, signal$r_median_top300[index]),
  near(comparison$prop_rhythmic, signal$prop_rhythmic[index]))
crossings <- utils::read.csv(source_files[["N80"]], stringsAsFactors = FALSE)
crossings <- crossings[crossings$time_basis == "TOD", , drop = FALSE]
index_crossings <- match(comparison$case_id, crossings$case_id)
stopifnot(nrow(crossings) == 10L, !anyNA(index_crossings))
for (column in intersect(names(crossings), names(comparison))) {
  left <- comparison[[column]]
  right <- crossings[[column]][index_crossings]
  stopifnot(if (is.numeric(left)) near(left, right) else identical(left, right))
}
index_power <- match(power$case_id, comparison$case_id)
stopifnot(!anyNA(index_power), setequal(power$case_id, comparison$case_id),
  !anyDuplicated(power[c("case_id", "N")]),
  near(power$n_pilot, comparison$n_pilot[index_power]),
  near(power$genes_fitted, comparison$genes_fitted[index_power]),
  all(power$G_simulated == settings$G), all(power$nsims >= settings$nsims_initial),
  all(is.finite(power$power)), all(power$power >= 0 & power$power <= 1))

# 3. Write exact saved values to a new set of report files -----------------
outputs <- c("GTEx_v10_TOD_N80.csv", "GTEx_v10_TOD_power_curves.csv",
  "GTEx_v10_TOD_signal_summary.csv", "pilot_validation.csv",
  "n80_TOD_ten_tissues.png", "n80_TOD_ten_tissues.pdf",
  "power_curves_TOD_ten_tissues.png", "power_curves_TOD_ten_tissues.pdf",
  "GTEx_v10_TOD_report.html", "report_recipe.json", "artifact_sha256.json")
if (any(file.exists(file.path(report, outputs))))
  stop("The TOD report outputs already exist. Preserve them or choose a new SCP_GTEX_TOD_REPORT.")
dir.create(report, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(comparison_raw, file.path(report, outputs[1]), row.names = FALSE)
utils::write.csv(power_raw, file.path(report, outputs[2]), row.names = FALSE)
stopifnot(file.copy(signal_file, file.path(report, outputs[3]), overwrite = FALSE))
utils::write.csv(validation, file.path(report, outputs[4]), row.names = FALSE)
stopifnot(identical(comparison_raw, utils::read.csv(file.path(report, outputs[1]),
  colClasses = "character", check.names = FALSE)),
  identical(power_raw, utils::read.csv(file.path(report, outputs[2]),
    colClasses = "character", check.names = FALSE)),
  identical(digest::digest(file = signal_file, algo = "sha256"),
    digest::digest(file = file.path(report, outputs[3]), algo = "sha256")))

# 4. Draw the N80 estimates and the saved TOD power curves -----------------
library(ggplot2)
tissue_names <- c(Adrenal = "Adrenal gland", Liver = "Liver",
  Muscle = "Skeletal muscle", Muscle_Skeletal = "Skeletal muscle",
  Adipose_Subcutaneous = "Subcutaneous adipose", Artery_Tibial = "Tibial artery",
  Brain_Cortex = "Cerebral cortex", Thyroid = "Thyroid", Lung = "Lung",
  Pancreas = "Pancreas", Heart_Left_Ventricle = "Heart, left ventricle")
labels <- function(x) {
  values <- unname(tissue_names[x])
  values[is.na(values)] <- gsub("_", " ", x[is.na(values)], fixed = TRUE)
  values
}
tissue_order <- comparison$tissue
comparison$tissue_label <- factor(labels(comparison$tissue),
                                   levels = rev(labels(tissue_order)))
color <- "#227C83"
n80_plot <- ggplot(comparison, aes(x = N80_interpolated, y = tissue_label)) +
  geom_errorbar(aes(xmin = N80_mc95_low, xmax = N80_mc95_high),
                orientation = "y", width = 0.22, color = color) +
  geom_point(size = 3, color = color) +
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0.01, 0.05))) +
  labs(title = "GTEx v10: sample size for 80% detection power",
    subtitle = "Season-adjusted TOD. Passive sampling. Single-harmonic cosinor.",
    x = "N80 (total samples)", y = NULL,
    caption = "Intervals show Monte Carlo uncertainty with fixed pilots.") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.caption = element_text(hjust = 0))
power$tissue_label <- factor(labels(power$tissue), levels = labels(tissue_order))
curve_plot <- ggplot(power, aes(x = N, y = power)) +
  geom_hline(yintercept = 0.8, linetype = 3, color = "grey50") +
  geom_ribbon(aes(ymin = power_mc95_low, ymax = power_mc95_high),
              fill = color, alpha = 0.14) +
  geom_line(color = color, linewidth = 0.65) +
  geom_point(color = color, size = 1) +
  facet_wrap(~tissue_label, ncol = 2) +
  scale_x_log10(breaks = c(12, 48, 192, 768, 3072, 10000)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = "GTEx v10: detection power with season-adjusted TOD",
    x = "Total samples (logarithmic scale)", y = "Detection power",
    caption = "K = 1. BH at 0.05. Shaded bands show fixed-pilot Monte Carlo intervals.") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), plot.caption = element_text(hjust = 0))
png_device <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else grDevices::png
ggsave(file.path(report, outputs[5]), n80_plot, width = 11, height = 7.5,
       dpi = 170, device = png_device)
ggsave(file.path(report, outputs[6]), n80_plot, width = 11, height = 7.5,
       device = grDevices::pdf, useDingbats = FALSE)
ggsave(file.path(report, outputs[7]), curve_plot, width = 11, height = 14,
       dpi = 170, device = png_device)
ggsave(file.path(report, outputs[8]), curve_plot, width = 11, height = 14,
       device = grDevices::pdf, useDingbats = FALSE)

# 5. Embed the figures and all tissue summaries in an offline HTML report --
escape_html <- function(x) {
  x <- gsub("&", "&amp;", as.character(x), fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  gsub(">", "&gt;", x, fixed = TRUE)
}
render_table <- function(x, id) {
  rows <- apply(x, 1, function(row) paste0("<tr><td>",
    paste(escape_html(row), collapse = "</td><td>"), "</td></tr>"))
  c(paste0('<div class="table"><table id="', id, '"><thead><tr><th>'),
    paste(escape_html(names(x)), collapse = "</th><th>"),
    "</th></tr></thead><tbody>", rows, "</tbody></table></div>")
}
image_tag <- function(name, alt) paste0('<img alt="', alt,
  '" src="data:image/png;base64,', base64enc::base64encode(file.path(report, name)), '">')
index <- match(canonical_tissue, signal$tissue)
shown_n80 <- data.frame(Tissue = labels(comparison$tissue),
  n_pilot = comparison$n_pilot, genes_fitted = comparison$genes_fitted,
  n_selected = comparison$n_selected,
  median_r = sprintf("%.3f", comparison$median_r),
  prop_rhythmic = sprintf("%.4f", comparison$prop_rhythmic),
  n_BH05 = signal$n_BH05[index], N80 = sprintf("%.1f", comparison$N80_interpolated),
  MC95_low = sprintf("%.1f", comparison$N80_mc95_low),
  MC95_high = sprintf("%.1f", comparison$N80_mc95_high), check.names = FALSE)
shown_signal <- data.frame(Tissue = labels(signal$tissue), n_pilot = signal$n,
  genes_fitted = signal$ngenes, n_selected = signal$n_selected,
  median_r = sprintf("%.3f", signal$r_median_top300),
  prop_rhythmic = sprintf("%.4f", signal$prop_rhythmic),
  n_BH05 = signal$n_BH05, check.names = FALSE)
html <- c('<!doctype html><html lang="en"><meta charset="utf-8">',
  '<meta name="viewport" content="width=device-width, initial-scale=1">',
  '<title>GTEx v10: season-adjusted TOD</title>',
  '<style>body{max-width:1120px;margin:32px auto;padding:0 20px;font:16px/1.55 system-ui;color:#17212b}h1,h2{line-height:1.25}h2{margin-top:32px}img{width:100%;height:auto}table{border-collapse:collapse;font-size:13px;width:100%}th,td{padding:7px;border-bottom:1px solid #d9e1e5;text-align:left;white-space:nowrap}thead{background:#edf4f4}.table{overflow:auto}p{max-width:940px}a{color:#176a72}.note{color:#46545d}</style>',
  '<h1>GTEx v10 with season-adjusted TOD</h1>',
  '<p>SCP uses pilot expression data and sample times to estimate power under a specified circadian model and sampling design. This report uses GTEx v10 TMM-normalized log-CPM and season-adjusted time of death (TOD). Tissue expression means are retained. The fits and power calculations come from the new normalized-expression build.</p>',
  '<p>Pilot selection uses nominal p &lt; 0.01 and retains up to 300 genes ranked by p-value. The median signal-to-noise ratio describes those selected genes. The n_selected column gives their number. The prop_rhythmic column is the fraction of all fitted genes passing the nominal pilot screen before the 300-gene cap. The n_BH05 column separately counts discoveries at BH FDR 5% across all fitted genes.</p>',
  '<h2>Sample size for 80% detection power</h2>',
  paste0('<p>The ten-tissue calculation uses passive sampling, a single-harmonic cosinor test (K = 1), BH at 0.05, ',
    settings$G, ' simulated genes, ', settings$nsims_initial,
    ' initial experiments per sample size, and at least ', settings$nsims_final,
    ' experiments near each N80 bracket. The saved calculation uses the exact conditional Gaussian F distribution.</p>'),
  '<p>N80 is the interpolated sample size at which the mean proportion of simulated true rhythmic genes detected reaches 80%. Intervals describe Monte Carlo uncertainty conditional on the fitted pilots. They do not include uncertainty in the pilot estimates or grid interpolation. The calculation describes the fitted model and sampling design.</p>',
  image_tag(outputs[5], 'N80 estimates for ten tissues using season-adjusted TOD'),
  render_table(shown_n80, 'n80-table'),
  '<p class="note">The tables round values for display. The accompanying CSV files retain the saved values without rounding.</p>',
  image_tag(outputs[7], 'Saved TOD detection power curves for ten tissues'),
  '<h2>Signal summary for all 46 tissues</h2>',
  render_table(shown_signal, 'signal-table'),
  '<h2>Sources and reproducibility</h2>',
  paste0('<p>The report matches all 46 signal-summary rows to the staged normalized TOD pilot objects and verifies the ten sample-size calibrations against those same objects. It retains ',
    nrow(power), ' saved sample-size points across the ten tissues. This rendering step does not rerun pilot fits or Monte Carlo experiments.</p>'),
  '<p>Expression is TMM-normalized log-CPM without covariate regression. The build retains samples with valid season-adjusted TOD, maps the annotation symbols, keeps their first occurrence, and then requires mean log-CPM above 3. The original Liver reference retains 262 samples, 8,917 genes, 210 single-harmonic discoveries and 725 joint two-harmonic discoveries at BH FDR 5%.</p>',
  '<p>Normalization source: <a href="https://doi.org/10.1126/science.add0846">Talamanca, Gobet and Naef (2023)</a>. GTEx project: <a href="https://doi.org/10.1126/science.aaz1776">GTEx Consortium (2020)</a>. Expression release: v10.</p>',
  '</html>')
writeLines(html, file.path(report, outputs[9]))

# 6. Save source identities, validation results, and output checksums -------
recipe <- list(created_at_UTC = format(Sys.time(), tz = "UTC", usetz = TRUE),
  expression_release = "v10", input_time = "Season-adjusted TOD", K = 1L,
  expression_preprocessing = "TMM-normalized log-CPM without covariate regression",
  gene_filter = "First valid annotation symbol, then mean normalized log-CPM>3",
  builder_sha256 = digest::digest(file = "examples/pilots/gtex_v10/10_build_normalized_TOD.R", algo = "sha256"),
  simulation_generator_sha256 = digest::digest(file = "examples/pilots/gtex_v10/11_compare_normalized_TOD_sample_sizes.R", algo = "sha256"),
  source_files = as.list(normalizePath(source_files)),
  source_sha256 = as.list(setNames(vapply(source_files, function(path)
    digest::digest(file = path, algo = "sha256"), character(1)), names(source_files))),
  generator_sha256 = digest::digest(file = "examples/pilots/gtex_v10/12_render_normalized_TOD_report.R", algo = "sha256"),
  publication = normalizePath(publication), report = normalizePath(report),
  nominal_pilot_threshold = 0.01, selected_gene_cap = 300L, BH_threshold = 0.05,
  signal_tissues = nrow(signal), N80_cases = nrow(comparison),
  saved_power_points = nrow(power), saved_MC_experiments = sum(power$nsims),
  pilot_validation_checks = sum(validation$checks_passed),
  source_N80_cells_retained_exactly = TRUE, source_power_cells_retained_exactly = TRUE,
  signal_CSV_copied_exactly = TRUE, refit_performed = FALSE,
  new_simulation_performed = FALSE,
  MC_interval_scope = "Fixed pilots. Excludes pilot, timing, cohort, and interpolation uncertainty.")
jsonlite::write_json(recipe, file.path(report, outputs[10]), auto_unbox = TRUE,
                     pretty = TRUE, digits = 16)
checksums <- setNames(vapply(file.path(report, outputs[1:10]), function(path)
  digest::digest(file = path, algo = "sha256"), character(1)), outputs[1:10])
jsonlite::write_json(as.list(checksums), file.path(report, outputs[11]),
                     auto_unbox = TRUE, pretty = TRUE)
cat("Normalized TOD report written to", normalizePath(report), "\n",
    nrow(comparison), "N80 cases and", nrow(signal), "tissue summaries verified.\n")
