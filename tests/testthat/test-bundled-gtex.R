test_that("the active GTEx bundle contains 46 v10 TOD pilots and both harmonic orders", {
  manifest <- scp_pilots(time_basis = "all")
  gtex <- manifest[manifest$dataset == "GTEx", , drop = FALSE]
  expect_equal(nrow(gtex), 46L)
  expect_true(all(gtex$gtex_release == "v10"))
  expect_true(all(gtex$time_basis == "TOD"))
  expect_equal(length(unique(gtex$tissue)), 46L)
  key <- paste(gtex$species, gtex$dataset, gtex$tissue, gtex$condition, sep = "/")
  expect_false(anyDuplicated(key) > 0L)

  pilot_dir <- SCP:::.scp_pilots_dir()
  tod_files <- list.files(file.path(pilot_dir, "human", "gtex_v10", "tod"),
                          pattern = "[.]rds$", recursive = TRUE)
  expected_files <- c(gtex$file, sub("[.]rds$", "_2H.rds", gtex$file))
  expect_length(tod_files, 92L)
  expect_setequal(tod_files, basename(expected_files))
  expect_length(list.files(file.path(pilot_dir, "human", "gtex_v10", "dip"),
                             pattern = "[.]rds$", recursive = TRUE), 0L)
  expect_error(scp_load_pilot("human", "GTEx", "Adrenal", "All", time_basis = "DIP"),
                 "DIP pilots are not bundled")

  for (i in seq_len(nrow(gtex))) {
    row <- gtex[i, ]
    for (harmonic in 1:2) {
      pilot <- scp_load_pilot(row$species, row$dataset, row$tissue, row$condition,
                               K = harmonic)
      expect_s3_class(pilot, "CircadianBioOptions")
      expect_identical(pilot$gtex_release, "v10")
      expect_identical(pilot$time_basis, "TOD")
      expect_s3_class(pilot$time_distribution, "SCPTimeDistribution")
      expect_true(is.null(pilot$times))
      expect_true(is.null(pilot$cts))
      expect_true(all(is.finite(pilot$amplitude)))
      expect_true(all(is.finite(pilot$sigma_rhythmic)))
      expect_equal(length(pilot$amplitude), length(pilot$sigma_rhythmic))
      full <- pilot$explorer_fit
      expect_s3_class(full, "data.frame")
      expect_equal(nrow(full), row$ngenes)
      expect_equal(full$q_K1, stats::p.adjust(full$p_K1, "BH", n = row$ngenes),
                     tolerance = 1e-12)
      if (harmonic == 2L) {
        expect_equal(full$q_joint, stats::p.adjust(full$p_joint, "BH", n = row$ngenes),
                       tolerance = 1e-12)
      }
    }
  }
})


test_that("GTEx sampling histograms contain only complete aggregate hourly counts", {
  path <- file.path(SCP:::.scp_pilots_dir(), "gtex_v10_TOD_sampling_histograms.csv")
  expect_true(file.exists(path))
  if (file.exists(path)) {
    histogram <- utils::read.csv(path, stringsAsFactors = FALSE)
    expect_identical(names(histogram), c("dataset", "tissue", "condition",
      "gtex_release", "time_basis", "bin_start", "bin_end", "count", "n_samples"))
    expect_equal(nrow(histogram), 46L * 24L)
    expect_false(anyNA(histogram))
    expect_true(all(histogram$dataset == "GTEx"))
    expect_true(all(histogram$gtex_release == "v10"))
    expect_true(all(histogram$time_basis == "TOD"))
    expect_true(all(histogram$count >= 0 & histogram$count == floor(histogram$count)))
    expect_equal(histogram$bin_end, histogram$bin_start + 1)

    manifest <- scp_pilots()
    manifest <- manifest[manifest$dataset == "GTEx", , drop = FALSE]
    columns <- c("dataset", "tissue", "condition", "gtex_release", "time_basis")
    key <- function(x) do.call(paste, c(x[columns], sep = "/"))
    groups <- split(histogram, key(histogram))
    expect_length(groups, 46L)
    expect_setequal(names(groups), key(manifest))
    expect_true(all(vapply(groups, function(x) nrow(x) == 24L, logical(1))))
    expect_true(all(vapply(groups, function(x)
      identical(as.numeric(sort(x$bin_start)), as.numeric(0:23)), logical(1))))
    totals <- vapply(groups, function(x) sum(x$count), numeric(1))
    expect_equal(unname(totals[key(manifest)]), as.numeric(manifest$n))
    expect_true(all(vapply(groups, function(x) length(unique(x$n_samples)) == 1L,
      logical(1))))
    samples <- vapply(groups, function(x) x$n_samples[1], numeric(1))
    expect_equal(unname(samples[key(manifest)]), as.numeric(manifest$n))
  }
})


test_that("the normalized Liver pilot reproduces the v10 reference discoveries", {
  first <- scp_load_pilot("human", "GTEx", "Liver", "All", K = 1)
  second <- scp_load_pilot("human", "GTEx", "Liver", "All", K = 2)
  expect_equal(first$n_samples, 262L)
  expect_equal(first$ngenes, 8917L)
  a <- first$explorer_fit
  b <- second$explorer_fit
  index <- match(a$gene, b$gene)
  expect_false(anyNA(index))
  expect_equal(sum(a$qvalue < 0.05), 210L)
  expect_equal(sum(b$q_joint < 0.05), 725L)
  expect_equal(sum(a$qvalue < 0.05 & b$q_joint[index] < 0.05), 207L)
})
