pilot_fixture <- function() {
  directory <- tempfile("scp-pilot-fixture-")
  dir.create(directory)
  distribution <- structure(
    list(grid = seq(0, 24, length.out = 512), weights = rep(1, 512),
         period = 24, jitter_sd = 0.1, n_observed = 24L),
    class = "SCPTimeDistribution")
  fit <- data.frame(gene = c("g1", "g2", "g3"),
                    pvalue = c(0.0001, 0.001, 0.05),
                    A = c(1, 0.8, 0.2), phi = c(3, 9, 15),
                    sigma = c(0.5, 0.6, 1), mesor = c(1, 2, 3))
  pilot <- structure(
    list(ngenes = 10L, lBaselineExpr = rep(2, 10), lOD = rep(0, 10),
         amplitude = c(1, 0.8), phase = c(3, 9),
         sigma_rhythmic = c(0.5, 0.6), prop_rhythmic = 0.2,
         paired_sigma = TRUE, period = 24, alpha_pilot = 0.01,
         pilot_top_k = 300L, pilot_cap = 0.2, rhythm_fit = fit,
         time_distribution = distribution, n_pilot = 24L),
    class = "CircadianBioOptions")
  manifest <- data.frame(
    species = c("human", "human", "human", "mouse"),
    dataset = c("GTEx", "GTEx", "GTEx", "Synthetic"),
    tissue = c("Liver", "Liver", "Liver", "Liver"),
    condition = "All", n = 24L, ngenes = 10L,
    design = c("passive", "passive", "passive", "active"),
    sample_type = "in_vivo", assay = "RNA-seq", accession = "fixture",
    citation = "Synthetic test data", license = "CC0",
    file = c("liver_dip.rds", "liver_tod.rds", "liver_v8.rds", "mouse.rds"),
    status = "ingested", gtex_release = c("v10", "v10", "v8", ""),
    time_basis = c("DIP", "TOD", "DIP", "native"),
    time_source = c("CHIRAL", "season-adjusted TOD", "v8", "native"),
    is_default = c(TRUE, FALSE, FALSE, TRUE), stringsAsFactors = FALSE)
  for (i in seq_len(nrow(manifest))) {
    p <- pilot
    p$gtex_release <- manifest$gtex_release[i]
    p$time_basis <- manifest$time_basis[i]
    p$source_provenance <- list(release = p$gtex_release,
                                time_basis = p$time_basis)
    p$fixture_marker <- manifest$file[i]
    saveRDS(p, file.path(directory, manifest$file[i]))
    p$paired_2h <- TRUE
    p$amplitude2 <- c(0.2, 0.3)
    p$phase2 <- c(1, 4)
    p$rhythm_fit$A2 <- c(0.2, 0.3, 0.1)
    p$rhythm_fit$phi2 <- c(1, 4, 7)
    p$fixture_marker <- paste0(manifest$file[i], "_2H")
    saveRDS(p, file.path(directory, sub("[.]rds$", "_2H.rds", manifest$file[i])))
  }
  utils::write.csv(manifest, file.path(directory, "manifest.csv"), row.names = FALSE)
  list(directory = directory, manifest = manifest, pilot = pilot)
}

fixture_database <- function(fixture, env = parent.frame()) {
  testthat::local_mocked_bindings(
    .scp_pilots_dir = function() fixture$directory,
    .package = "SCP", .env = env)
}
