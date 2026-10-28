# Exercise preparation, fitting, and stage-only output with synthetic data.

# 1. Create a disposable stage and a two-model source index ---------------
fixture_stage <- Sys.getenv("SCP_NATIVE_SMOKE_STAGE",
  file.path("output", "pilot_rebuild",
            paste0("native_fixture_", format(Sys.time(), "%Y%m%d_%H%M%S"))))
dir.create(fixture_stage, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(SCP_NATIVE_STAGE = fixture_stage)
fixture_index <- data.frame(species = "mouse", dataset = "SYNTHETIC",
  accession = "SYNTHETIC", base_file = "mouse/SYNTHETIC_All.rds",
  file = c("mouse/SYNTHETIC_All.rds", "mouse/SYNTHETIC_All_2H.rds"),
  K = 1:2, period = 24, alpha_pilot = 0.01, top_k = c(300L, 500L),
  seed = 2025L, builder = "examples/pilots/native/04_fit_pilots.R")
fixture_index_file <- file.path(fixture_stage, "source_index.csv")
write.csv(fixture_index, fixture_index_file, row.names = FALSE)
Sys.setenv(SCP_NATIVE_INDEX = fixture_index_file, SCP_NATIVE_DATASETS = "")
source("examples/pilots/native/00_setup.R")
source("examples/pilots/native/geo_helpers.R")

# 2. Supply a known rhythmic waveform and check time-label parsing --------
set.seed(17)
times <- rep(seq(0, 22, by = 2), 2)
expr <- matrix(rnorm(120 * length(times), sd = 0.1), 120)
expr <- expr + rep(seq(0.5, 1.5, length.out = 120), length(times)) *
  rep(cos(2 * pi * times / 24) + 0.3 * cos(4 * pi * times / 24), each = 120)
rownames(expr) <- paste0("gene", seq_len(nrow(expr)))
source_file <- file.path(fixture_stage, "synthetic_expression.rds")
saveRDS(expr, source_file)
native_save_input(expr, times, fixture_index$base_file[1],
                  c(synthetic = source_file), "Synthetic harmonic fixture.")
parsed <- native_time_labels(c("ZT06", "CT18", "12am", "12pm", "1:30pm", "sample123"))
stopifnot(isTRUE(all.equal(parsed[1:5], c(6, 18, 0, 12, 13.5))), is.na(parsed[6]))

# 3. Run the same fitting script used for prepared real inputs ------------
source("examples/pilots/native/04_fit_pilots.R")
one <- readRDS(file.path(fixture_stage, "pilots", fixture_index$file[1]))
two <- readRDS(file.path(fixture_stage, "pilots", fixture_index$file[2]))
stopifnot(one$ngenes == 120L, two$ngenes == 120L,
          identical(one$times, times), isTRUE(two$paired_2h),
          identical(one$tested_genes, rownames(expr)),
          file.exists(file.path(fixture_stage, "fit_comparison.csv")))
# 4. Check the saved harmonic coefficients against an independent QR fit --
design <- cbind(1, cos(2 * pi * times / 24), sin(2 * pi * times / 24),
                cos(4 * pi * times / 24), sin(4 * pi * times / 24))
coef <- stats::lm.fit(design, t(expr))$coefficients
expected_A1 <- sqrt(coef[2, ]^2 + coef[3, ]^2)
expected_A2 <- sqrt(coef[4, ]^2 + coef[5, ]^2)
matched <- match(two$rhythm_fit$gene, rownames(expr))
stopifnot(!anyNA(matched),
          max(abs(two$rhythm_fit$A - expected_A1[matched])) < 1e-8,
          max(abs(two$rhythm_fit$A2 - expected_A2[matched])) < 1e-8)
cat("Native builder fixture passed. Stage:", fixture_stage, "\n")
