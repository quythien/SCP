test_that("every active GTEx time basis loads for both harmonic orders", {
  manifest <- scp_pilots(time_basis = "all")
  gtex <- manifest[manifest$dataset == "GTEx", , drop = FALSE]
  expect_gt(nrow(gtex), 0L)
  expect_true(all(gtex$gtex_release == "v10"))
  expect_true(all(gtex$time_basis %in% c("DIP", "TOD")))
  key <- paste(gtex$species, gtex$dataset, gtex$tissue, gtex$condition, sep = "/")
  expect_true(all(vapply(split(gtex$time_basis, key),
                          function(x) setequal(x, c("DIP", "TOD")), logical(1))))
  for (i in seq_len(nrow(gtex))) {
    row <- gtex[i, ]
    for (harmonic in 1:2) {
      pilot <- scp_load_pilot(row$species, row$dataset, row$tissue, row$condition,
                               K = harmonic, time_basis = row$time_basis)
      expect_s3_class(pilot, "CircadianBioOptions")
      expect_identical(pilot$gtex_release, "v10")
      expect_identical(pilot$time_basis, row$time_basis)
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
