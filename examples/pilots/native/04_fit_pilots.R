# Fit prepared inputs into the separate staging directory.
source("examples/pilots/native/00_setup.R")

# 1. Match each physical pilot to its preparation and stored settings ------
audit <- list()
for (i in seq_len(nrow(native_index))) {
  row <- native_index[i, , drop = FALSE]
  input_file <- file.path(native_stage, "inputs", row$base_file)
  if (!file.exists(input_file)) next
  input <- readRDS(input_file)
  set.seed(row$seed)

  # 2. Fit the selected harmonic order -----------------------------------
  if (row$K == 1L) {
    pilot <- SCP::estCircadianParam(input$expr, input$times, period = row$period,
                                    min_rhythm_pval = row$alpha_pilot,
                                    top_k = row$top_k, paired_sigma = TRUE,
                                    verbose = FALSE)
  } else {
    pilot <- SCP::estCircadianParam2H(input$expr, input$times, period = row$period,
                                      min_rhythm_pval = row$alpha_pilot,
                                      top_k = row$top_k, paired_sigma = TRUE,
                                      verbose = FALSE)
  }

  # 3. Attach source and gene information --------------------------------
  pilot$times <- input$times
  pilot$tested_genes <- rownames(input$expr)
  pilot$time_basis <- "native"
  symbols <- pilot$rhythm_fit$gene
  if (!is.null(input$symbol_map)) {
    matched <- match(symbols, names(input$symbol_map))
    known <- !is.na(matched) & !is.na(input$symbol_map[matched])
    symbols[known] <- unname(input$symbol_map[matched[known]])
  }
  pilot$rhythm_fit$symbol <- symbols
  pilot$source_provenance <- list(accession = row$accession,
    input_sha256 = input$source_hashes, recipe = row$builder,
    period = row$period, seed = row$seed, harmonic_order = row$K,
    alpha_pilot = row$alpha_pilot, top_k = row$top_k,
    paired_sigma = TRUE, preprocessing = input$notes)

  # 4. Save candidates and compare with the unchanged reference -----------
  output_file <- file.path(native_stage, "pilots", row$file)
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(pilot, output_file, compress = "xz")
  reference_file <- file.path("inst/extdata/pilots", row$file)
  equal <- NA
  if (file.exists(reference_file)) {
    reference <- readRDS(reference_file)
    fields <- intersect(c("gene", "pvalue", "A", "phi", "sigma", "A2", "phi2"),
                        names(pilot$rhythm_fit))
    equal <- isTRUE(all.equal(pilot$rhythm_fit[fields], reference$rhythm_fit[fields],
                              check.attributes = FALSE)) &&
      isTRUE(all.equal(pilot$prop_rhythmic, reference$prop_rhythmic))
  }
  audit[[length(audit) + 1L]] <- data.frame(file = row$file,
    output_sha256 = digest::digest(file = output_file, algo = "sha256"),
    numeric_reference_match = equal, n = ncol(input$expr), ngenes = nrow(input$expr))
}
if (!length(audit)) stop("No prepared inputs were found. Run the relevant preparation script first.")
write.csv(do.call(rbind, audit), file.path(native_stage, "fit_comparison.csv"),
          row.names = FALSE)
