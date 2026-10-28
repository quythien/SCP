test_that("GTEx defaults to v10 DIP and explicit TOD remains available", {
  fixture <- pilot_fixture()
  fixture_database(fixture)
  default <- scp_pilots()
  gtex <- default[default$dataset == "GTEx", , drop = FALSE]
  expect_equal(nrow(gtex), 1L)
  expect_identical(gtex$gtex_release, "v10")
  expect_identical(gtex$time_basis, "DIP")
  expect_identical(scp_load_pilot("human", "GTEx", "Liver")$fixture_marker,
                   "liver_dip.rds")
  expect_identical(scp_load_pilot("human", "GTEx", "Liver", time_basis = "TOD")$fixture_marker,
                   "liver_tod.rds")
  expect_identical(scp_load_pilot("human", "GTEx", "Liver", K = 2,
                                 time_basis = "TOD")$fixture_marker,
                   "liver_tod.rds_2H")
})

test_that("listing both time bases never exposes archived GTEx", {
  fixture <- pilot_fixture()
  fixture_database(fixture)
  manifest <- scp_pilots(status = "all", time_basis = "all")
  expect_setequal(manifest$time_basis[manifest$dataset == "GTEx"], c("DIP", "TOD"))
  expect_true(all(manifest$gtex_release[manifest$dataset == "GTEx"] == "v10"))
  expect_equal(nrow(scp_pilot_search("Liver", species = "human", time_basis = "all")), 2L)
  expect_identical(scp_pilot_search("Liver", species = "human")$time_basis, "DIP")
  expect_error(scp_load_pilot("human", "GTEx", "Liver", time_basis = "all"))
})

test_that("native datasets preserve their default behavior", {
  fixture <- pilot_fixture()
  fixture_database(fixture)
  expect_identical(scp_load_pilot("mouse", "Synthetic", "Liver")$fixture_marker,
                   "mouse.rds")
  manifest <- scp_pilots(time_basis = "TOD")
  expect_true(any(manifest$dataset == "Synthetic"))
  expect_error(scp_load_pilot("mouse", "Synthetic", "Liver", time_basis = "DIP"))
})

test_that("mismatched GTEx object metadata fails closed", {
  fixture <- pilot_fixture()
  fixture_database(fixture)
  p <- readRDS(file.path(fixture$directory, "liver_dip.rds"))
  p$gtex_release <- "v8"
  p$source_provenance$release <- "v8"
  saveRDS(p, file.path(fixture$directory, "liver_dip.rds"))
  expect_error(scp_load_pilot("human", "GTEx", "Liver"))
})

test_that("invalid distribution metadata cannot enter a GTEx pilot", {
  fixture <- pilot_fixture()
  fixture_database(fixture)
  p <- readRDS(file.path(fixture$directory, "liver_dip.rds"))
  p$time_distribution$weights[1] <- -1
  saveRDS(p, file.path(fixture$directory, "liver_dip.rds"))
  expect_error(scp_load_pilot("human", "GTEx", "Liver"))
})
