test_that("GTEx defaults to current v10 seasonal TOD", {
  fixture <- pilot_fixture()
  fixture_database(fixture)
  manifest <- scp_pilots()
  gtex <- manifest[manifest$dataset == "GTEx", , drop = FALSE]
  expect_equal(nrow(gtex), 1L)
  expect_identical(gtex$gtex_release, "v10")
  expect_identical(gtex$time_basis, "TOD")
  expect_identical(scp_load_pilot("human", "GTEx", "Liver")$fixture_marker,
                   "liver_tod.rds")
  expect_identical(scp_load_pilot("human", "GTEx", "Liver", time_basis = "TOD")$fixture_marker,
                   "liver_tod.rds")
  expect_identical(scp_load_pilot("human", "GTEx", "Liver", K = 2)$fixture_marker,
                   "liver_tod.rds_2H")
})

test_that("all listings and searches expose only current GTEx TOD", {
  fixture <- pilot_fixture()
  archived <- fixture$manifest[fixture$manifest$gtex_release == "v8", , drop = FALSE]
  archived$time_basis <- "TOD"
  archived$file <- "liver_v8_tod.rds"
  fixture$manifest <- rbind(fixture$manifest, archived)
  utils::write.csv(fixture$manifest, file.path(fixture$directory, "manifest.csv"), row.names = FALSE)
  old_pilot <- readRDS(file.path(fixture$directory, "liver_v8.rds"))
  old_pilot$time_basis <- old_pilot$source_provenance$time_basis <- "TOD"
  saveRDS(old_pilot, file.path(fixture$directory, archived$file))
  fixture_database(fixture)
  for (basis in list(NULL, "TOD", "all")) {
    manifest <- scp_pilots(status = "all", time_basis = basis)
    gtex <- manifest[manifest$dataset == "GTEx", , drop = FALSE]
    expect_equal(nrow(gtex), 1L)
    expect_identical(gtex$gtex_release, "v10")
    expect_identical(gtex$time_basis, "TOD")
    search <- scp_pilot_search("Liver", species = "human", time_basis = basis)
    expect_equal(nrow(search), 1L)
    expect_identical(search$time_basis, "TOD")
  }
  expect_error(scp_load_pilot("human", "GTEx", "Liver", time_basis = "all"),
                 "Use time_basis")
})

test_that("DIP requests fail clearly without selecting another pilot", {
  fixture <- pilot_fixture()
  fixture_database(fixture)
  expect_error(scp_pilots(time_basis = "DIP"), "DIP pilots are not bundled")
  expect_error(scp_pilot_search("Liver", time_basis = "DIP"), "DIP pilots are not bundled")
  expect_error(scp_load_pilot("human", "GTEx", "Liver", time_basis = "DIP"),
                 "DIP pilots are not bundled")
})

test_that("a missing current TOD row cannot fall back to DIP or v8", {
  fixture <- pilot_fixture()
  fixture$manifest <- fixture$manifest[fixture$manifest$time_basis != "TOD", , drop = FALSE]
  utils::write.csv(fixture$manifest, file.path(fixture$directory, "manifest.csv"), row.names = FALSE)
  fixture_database(fixture)
  manifest <- scp_pilots(status = "all", time_basis = "all")
  expect_false(any(manifest$dataset == "GTEx"))
  expect_error(scp_load_pilot("human", "GTEx", "Liver"), "No ingested pilot")
})

test_that("a missing TOD file cannot fall back to an available DIP file", {
  fixture <- pilot_fixture()
  fixture$manifest$file[fixture$manifest$time_basis == "TOD"] <- "missing_tod.rds"
  utils::write.csv(fixture$manifest, file.path(fixture$directory, "manifest.csv"), row.names = FALSE)
  fixture_database(fixture)
  expect_error(scp_load_pilot("human", "GTEx", "Liver"), "file is missing")
})

test_that("native datasets preserve their default behavior", {
  fixture <- pilot_fixture()
  fixture_database(fixture)
  expect_identical(scp_load_pilot("mouse", "Synthetic", "Liver")$fixture_marker,
                   "mouse.rds")
  manifest <- scp_pilots(time_basis = "TOD")
  expect_true(any(manifest$dataset == "Synthetic"))
  expect_error(scp_load_pilot("mouse", "Synthetic", "Liver", time_basis = "TOD"),
                 "Use time_basis = NULL")
})

test_that("mismatched TOD release and time-basis metadata fail closed", {
  fixture <- pilot_fixture()
  fixture_database(fixture)
  path <- file.path(fixture$directory, "liver_tod.rds")
  original <- readRDS(path)
  changes <- list(
    function(p) {
      p$gtex_release <- "v8"
      p
    },
    function(p) {
      p$source_provenance$release <- "v8"
      p
    },
    function(p) {
      p$time_basis <- "DIP"
      p
    },
    function(p) {
      p$source_provenance$time_basis <- "DIP"
      p
    })
  for (change in changes) {
    saveRDS(change(original), path)
    expect_error(scp_load_pilot("human", "GTEx", "Liver"), "v10 TOD metadata")
  }
})

test_that("invalid TOD sampling distributions are rejected", {
  fixture <- pilot_fixture()
  fixture_database(fixture)
  path <- file.path(fixture$directory, "liver_tod.rds")
  p <- readRDS(path)
  p$time_distribution$weights[1] <- -1
  saveRDS(p, path)
  expect_error(scp_load_pilot("human", "GTEx", "Liver"))
  p$time_distribution <- NULL
  saveRDS(p, path)
  expect_error(scp_load_pilot("human", "GTEx", "Liver"), "no aggregate input-time distribution")
})

test_that("invalid time-basis arguments give a clear error", {
  for (basis in list(NA_character_, character(), c("TOD", "all"), 1, "other"))
    expect_error(scp_pilots(time_basis = basis), "time_basis")
})
