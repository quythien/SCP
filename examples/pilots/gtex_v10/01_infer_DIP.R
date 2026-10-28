# Infer donor input time using the CHIRAL package and its published workflow.
# Sample names, fitted tissue times, and donor times stay in the private stage.

# 1. Load the adjusted strict matrices and published preparation functions --
source("examples/pilots/gtex_v10/00_settings.R")
strict <- read_expression(strict_path, "CPM.all.norm")
frames <- lapply(strict, function(item) {
  frame <- as.data.frame(tissue_counts(item))
  colnames(frame) <- gsub("-", ".", colnames(frame), fixed = TRUE)
  frame
})
E <- upstream$CPM_to_E(frames)
rm(strict, frames)
gc()
tissue_names <- names(E)
stopifnot(length(tissue_names) > 1L)
checkpoint_dir <- file.path(stage, "private", "published_chiral")
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
input_sha <- file_hash(strict_path)
source_sha <- file_hash(recipe_file)
run_identity <- digest::digest(list(recipe, input_sha, source_sha), algo = "sha256")
workers <- max(1L, as.integer(Sys.getenv("SCP_GTEX_INFERENCE_WORKERS", "4")))
json_write(c(recipe, list(strict_input_sha256 = input_sha,
  author_DIP_inf_sha256 = source_sha, run_identity = run_identity,
  RNG = "L'Ecuyer-CMRG streams in fixed tissue order", workers = workers)),
  file.path(stage, "qc", "inference_recipe.json"))

# 2. Allocate reproducible streams for one published fit per tissue ---------
RNGkind("L'Ecuyer-CMRG")
set.seed(recipe$inference_seed)
streams <- vector("list", length(E))
streams[[1L]] <- .Random.seed
if (length(E) > 1L) for (i in 2:length(E))
  streams[[i]] <- parallel::nextRNGStream(streams[[i - 1L]])
run_tissue <- function(i) {
  tissue <- tissue_names[i]
  checkpoint <- file.path(checkpoint_dir, paste0(tissue, ".rds"))
  if (file.exists(checkpoint)) {
    old <- readRDS(checkpoint)
    stopifnot(identical(old$run_identity, run_identity))
    return(old)
  }
  assign(".Random.seed", streams[[i]], envir = globalenv())
  began <- proc.time()[["elapsed"]]
  out <- CHIRAL::CHIRAL(E[[i]]$E, iterations = recipe$max_iterations,
    clockgenes = clock_genes, standardize = TRUE, GTEx_names = TRUE,
    TSM = TRUE, q = 0.1, update.q = FALSE, pbar = FALSE)
  stopifnot(!is.null(out), length(out$phi) == ncol(E[[i]]$E),
            all(is.finite(out$phi)))
  out$tissue <- tissue
  converged <- all(c("geni", "clock") %in% names(out))
  fit_seconds <- proc.time()[["elapsed"]] - began
  fitted <- upstream$Fit_OUT(setNames(list(out), tissue), N.cores = 1L)
  oriented <- upstream$Set_OUT(fitted)[[tissue]]
  available <- intersect(clock_genes, toupper(gene_symbols(rownames(E[[i]]$E))))
  result <- list(out = oriented, native_converged = converged, iteration = out$iter,
    available_clock_genes = available, missing_clock_genes = setdiff(clock_genes, available),
    chiral_seconds = fit_seconds, elapsed_seconds = proc.time()[["elapsed"]] - began,
    run_identity = run_identity, tissue = tissue, rng_stream = streams[[i]])
  saveRDS(result, checkpoint, compress = "gzip")
  cat(tissue, "native convergence", converged, "iterations", out$iter,
      "CHIRAL seconds", round(fit_seconds, 1), "total seconds",
      round(result$elapsed_seconds, 1), "\n")
  flush.console()
  result
}
results <- parallel::mclapply(seq_along(E), run_tissue, mc.cores = workers,
  mc.preschedule = FALSE, mc.set.seed = FALSE)
stopifnot(!any(vapply(results, inherits, logical(1), "try-error")))
OUT <- setNames(lapply(results, function(x) x$out), tissue_names)

# 3. Apply the author's donor consensus with the original defaults ---------
people <- upstream$Unique_donors(OUT)
phi_matrix <- upstream$Create_phi_matrix(OUT, people)
phi <- upstream$Phi.from.phi_mat(phi_matrix, ct = 1.95)
stopifnot(all(is.finite(phi)))
donor_names <- paste0("GTEX-", names(phi))
donors <- data.frame(donor = donor_names, DIP = unname(phi * 12 / pi),
  DIP_radians = unname(phi), n_tissues = rowSums(is.finite(phi_matrix[names(phi), , drop = FALSE])),
  status = "available", stringsAsFactors = FALSE)
stopifnot(!anyDuplicated(donors$donor), all(donors$DIP >= 0 & donors$DIP < 24))
main_tips <- phi_matrix
rownames(main_tips) <- paste0("GTEX-", rownames(main_tips))
qc <- do.call(rbind, lapply(results, function(x) data.frame(tissue = x$tissue,
  n_samples = length(x$out$phi), n_genes = nrow(x$out$E),
  n_clock_genes = length(x$available_clock_genes),
  available_clock_genes = paste(x$available_clock_genes, collapse = ","),
  missing_clock_genes = paste(x$missing_clock_genes, collapse = ","),
  native_converged = x$native_converged, iterations = x$iteration,
  chiral_seconds = x$chiral_seconds, elapsed_seconds = x$elapsed_seconds,
  orientation_flipped = x$out$has.been.flipped, stringsAsFactors = FALSE)))
utils::write.csv(qc, file.path(stage, "qc", "tissue_inference.csv"), row.names = FALSE)
saveRDS(list(donors = donors, main_tips = main_tips, recipe = recipe,
  tissue_qc = qc, run_identity = run_identity),
  file.path(stage, "private", "v10_DIP_inference.rds"), compress = "xz")
json_write(list(n_source_donors = nrow(main_tips), n_with_DIP = nrow(donors),
  n_tissues = length(OUT), native_converged_tissues = sum(qc$native_converged),
  iteration_limit_tissues = qc$tissue[!qc$native_converged],
  consensus_min_tissues = 1L, consensus_ct = 1.95,
  consensus_function = "Phi.from.phi_mat", filtering = "Published procedure without extra stability filters"),
  file.path(stage, "qc", "donor_coverage.json"))
cat("Published CHIRAL inference complete:", nrow(donors), "donors from", length(OUT), "tissues.\n")
