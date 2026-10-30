test_that("the installed synthetic example supports both arms of the tutorials", {
  path <- system.file("extdata", "example_caudate_liver.rds", package = "SCP")
  expect_true(nzchar(path), info = "The synthetic Caudate/Liver fixture must be installed.")
  expect_true(file.exists(path))
  if (nzchar(path) && file.exists(path)) {
    example <- readRDS(path)
    expect_true(isTRUE(example$metadata$synthetic))
    expect_false(example$metadata$observed_individual_times_included)
    expect_false(example$metadata$donor_identifiers_included)
    expect_false(example$metadata$observed_expression_included)
    expect_identical(example$metadata$generator$seed, 20261029L)
    schedule <- rep(seq(0, 20, by = 4), each = 8L)
    for (arm in c("caudate", "liver")) {
      value <- example[[arm]]
      expect_identical(dim(value$expr), c(1200L, 48L))
      expect_true(all(is.finite(value$expr)))
      expect_identical(value$times, schedule)
      expect_true(isTRUE(value$synthetic))
      expect_identical(value$time_units, "hours")
      expect_true(all(grepl("^SyntheticGene[0-9]+$", rownames(value$expr))))
      expect_true(all(grepl("^Synthetic(Caudate|Liver)_[0-9]+$", colnames(value$expr))))
    }
    expect_identical(rownames(example$caudate$expr), rownames(example$liver$expr))
    expect_equal(nrow(example$ground_truth), 1200L)
  }
})
