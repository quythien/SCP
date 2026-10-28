test_that("aggregate sampling preserves the numeric circular KDE law", {
  times <- c(-2, 0, 1, 4, 6, 9, 13, 18, 21, 25)
  circular <- times %% 24
  kde <- stats::density(c(circular - 24, circular, circular + 24), from = 0, to = 24)
  distribution <- structure(
    list(grid = kde$x, weights = pmax(kde$y, 1e-10), period = 24,
         jitter_sd = 0.1, n_observed = length(times)),
    class = "SCPTimeDistribution")
  set.seed(4103)
  numeric_draw <- SCP:::sampleTimesFromDist(400, times)
  set.seed(4103)
  aggregate_draw <- SCP:::sampleTimesFromDist(400, distribution)
  expect_identical(aggregate_draw, numeric_draw)
  expect_true(all(aggregate_draw >= 0 & aggregate_draw < 24))
})

test_that("time distribution weights determine sampling", {
  distribution <- structure(
    list(grid = c(2, 12), weights = c(1, 0), period = 24,
         jitter_sd = 0, n_observed = 12L), class = "SCPTimeDistribution")
  expect_identical(SCP:::sampleTimesFromDist(30, distribution), rep(2, 30))
})

test_that("direct simulators require realized numeric schedules", {
  fixture <- pilot_fixture()
  distribution <- fixture$pilot$time_distribution
  expect_error(SCP:::simCircadianSingleCohort(fixture$pilot, distribution),
               "sampleTimesFromDist")
  pilot2 <- fixture$pilot
  pilot2$paired_2h <- TRUE
  pilot2$amplitude2 <- c(0.2, 0.3)
  pilot2$phase2 <- c(1, 4)
  expect_error(simCircadianSingleCohort2H(pilot2, distribution),
               "sampleTimesFromDist")
})
