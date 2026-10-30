diff_simulator <- simCircadianDiff
diff_options <- CircadianBioOptions
diff_runner <- runSimsDiff

diff_fixture <- function(...) {
  args <- list(ngenes = 120L, n1 = 24L, n2 = 24L,
    lBaselineExpr = seq(2, 6, length.out = 120), lBaselineExpr2 = rep(30, 120),
    lOD = log(seq(.3, .8, length.out = 120)), lOD2 = log(seq(1, 2, length.out = 120)),
    amplitude = c(.01, .4, 1.5), sigma_rhythmic = c(.1, .3, .9),
    amplitude2 = c(.02, .7, 2), sigma_rhythmic2 = c(.2, .5, 1.1),
    prop_rhythmic = .8, prop_DR = .3, prop_DP = .2, prop_DM = .1,
    phase_diff = c(-4, 5), dp_shift_mode = "uniform", mesor_diff = c(.4, 1.3),
    design = "active", cts = rep(seq(0, 20, by = 4), each = 4), sim.seed = 1701L)
  replacements <- list(...)
  args[names(replacements)] <- replacements
  args
}

test_that("differential labels match the realized endpoint nulls", {
  sim <- do.call(diff_simulator, diff_fixture())
  gt <- sim$ground_truth
  expect_setequal(unique(gt$diff_type), 0:5)
  expect_identical(gt$is_DR, gt$diff_type %in% c(2L, 3L))
  expect_identical(gt$is_DP, gt$diff_type == 4L)
  expect_identical(gt$is_DM, gt$diff_type == 5L)
  expect_equal(gt$mesor2[gt$diff_type != 5], gt$mesor1[gt$diff_type != 5], tolerance = 0)
  snr_difference <- gt$amplitude1^2 / gt$sigma^2 - gt$amplitude2^2 / gt$sigma2^2
  expect_identical(snr_difference[!gt$is_DR], rep(0, sum(!gt$is_DR)))
  both <- gt$is_rhythmic_g1 & gt$is_rhythmic_g2
  expect_identical(gt$phase_diff[both & !gt$is_DP], rep(0, sum(both & !gt$is_DP)))
  expect_true(all(abs(gt$mesor2[gt$is_DM] - gt$mesor1[gt$is_DM]) >= .4))
  expect_true(all(abs(gt$mesor2[gt$is_DM] - gt$mesor1[gt$is_DM]) <= 1.3))
  expect_true(any(gt$phase_diff[gt$is_DP] < 0))
  expect_true(any(gt$phase_diff[gt$is_DP] > 0))
  expect_true(any(abs(gt$phase_diff[gt$is_DP]) < 2))
  expect_equal(sim$effectsize_DR2, gt$amplitude2 / gt$sigma2, tolerance = 0)
})

test_that("both DR arms retain their source index before the existing bounds", {
  args <- diff_fixture(sigma_rhythmic = c(1e-8, .3, .9),
    sigma_rhythmic2 = c(2e-8, .5, 1.1))
  sim <- do.call(diff_simulator, args)
  gt <- sim$ground_truth
  draws <- sim$parameter_draws
  g1 <- which(gt$diff_type == 2L)
  g2 <- which(gt$diff_type == 3L)
  i1 <- draws$pair_index1[g1]
  i2 <- draws$pair_index2[g2]
  expect_false(anyNA(c(i1, i2)))
  expect_equal(draws$raw_amplitude1[g1], args$amplitude[i1], tolerance = 0)
  expect_equal(draws$raw_sigma1[g1], args$sigma_rhythmic[i1], tolerance = 0)
  expect_equal(draws$raw_amplitude2[g2], args$amplitude2[i2], tolerance = 0)
  expect_equal(draws$raw_sigma2[g2], args$sigma_rhythmic2[i2], tolerance = 0)
  expect_equal(gt$amplitude1[g1], pmax(args$amplitude[i1], .05), tolerance = 0)
  expect_equal(gt$sigma[g1], pmax(args$sigma_rhythmic[i1], 1e-6), tolerance = 0)
  expect_equal(gt$amplitude2[g2], pmax(args$amplitude2[i2], .05), tolerance = 0)
  expect_equal(gt$sigma2[g2], pmax(args$sigma_rhythmic2[i2], 1e-6), tolerance = 0)
  expect_equal(unname(sim$simOptions$transform_counts), c(
    sum(draws$raw_amplitude1 < .05, na.rm = TRUE),
    sum(draws$raw_amplitude2 < .05, na.rm = TRUE),
    sum(draws$raw_sigma1 < 1e-6, na.rm = TRUE),
    sum(draws$raw_sigma2 < 1e-6, na.rm = TRUE)))
  expect_true(all(sim$simOptions$transform_counts > 0))
})

test_that("the constructor resamples each group's selected tuples jointly", {
  A1 <- c(.01, .4, 1.5)
  S1 <- c(.1, .3, .9)
  A2 <- c(.02, .7, 2, 3)
  S2 <- c(.2, .5, 1.1, 1.7)
  bio <- diff_options(ngenes = 80, prop_rhythmic = .7, prop_DR = .2,
    prop_DP = .1, prop_DM = .1, lBaselineExpr = 5, lOD = log(.5),
    amplitude = A1, sigma_rhythmic = S1, amplitude2 = A2,
    sigma_rhythmic2 = S2, paired_sigma = TRUE)
  expect_length(bio$amplitude, 56)
  expect_equal(bio$sigma_rhythmic, S1[match(bio$amplitude, A1)], tolerance = 0)
  expect_equal(bio$sigma_rhythmic2, S2[match(bio$amplitude2, A2)], tolerance = 0)
  expect_true(bio$paired_sigma2)
  expect_error(diff_options(ngenes = 80, prop_rhythmic = .7,
    lBaselineExpr = 5, lOD = 0, amplitude = A1, amplitude2 = A2,
    sigma_rhythmic2 = S2[-1], paired_sigma = TRUE), "matched by gene")
  bad <- diff_fixture(sigma_rhythmic2 = c(.1, .2))
  expect_error(do.call(diff_simulator, bad), "matched by gene")
})

test_that("the two-group estimator carries noise from the same selected genes", {
  set.seed(846)
  times <- rep(seq(0, 22, by = 2), each = 4)
  X <- cbind(1, cos(2 * pi * times / 24), sin(2 * pi * times / 24))
  G <- 40L
  A1 <- seq(.8, 2, length.out = G)
  A2 <- A1
  A1[9:16] <- 0
  A2[1:8] <- 0
  A2[9:16] <- seq(1, 2, length.out = 8)
  phase1 <- seq(0, 20, length.out = G)
  phase2 <- phase1
  phase2[17:24] <- (phase1[17:24] + 4) %% 24
  M1 <- seq(3, 6, length.out = G)
  M2 <- M1
  M2[25:32] <- M1[25:32] + 1.2
  S1 <- seq(.12, .3, length.out = G)
  S2 <- seq(.2, .4, length.out = G)
  residual1 <- t(qr.resid(qr(X), matrix(rnorm(length(times) * G), nrow = length(times))))
  residual2 <- t(qr.resid(qr(X), matrix(rnorm(length(times) * G), nrow = length(times))))
  residual1 <- residual1 * (S1 / sqrt(rowSums(residual1^2) / (length(times) - 3)))
  residual2 <- residual2 * (S2 / sqrt(rowSums(residual2^2) / (length(times) - 3)))
  expr1 <- M1 + A1 * cos(outer(phase1, times, function(p, t) 2 * pi * (t - p) / 24)) + residual1
  expr2 <- M2 + A2 * cos(outer(phase2, times, function(p, t) 2 * pi * (t - p) / 24)) + residual2
  rownames(expr1) <- rownames(expr2) <- paste0("Fixture", seq_len(G))
  bio <- estCircadianParamTwoGroup(expr1, expr2, times, times,
    top_k = 32, verbose = FALSE)
  expect_true(bio$paired_sigma2)
  selected1 <- which(A1 > 0)
  selected2 <- which(A2 > 0)
  for (i in seq_along(bio$amplitude))
    expect_true(any(abs(A1[selected1] - bio$amplitude[i]) < 1e-10 &
      abs(S1[selected1] - bio$sigma_rhythmic[i]) < 1e-10))
  for (i in seq_along(bio$amplitude2))
    expect_true(any(abs(A2[selected2] - bio$amplitude2[i]) < 1e-10 &
      abs(S2[selected2] - bio$sigma_rhythmic2[i]) < 1e-10))
})

test_that("errors are independent while seeds and parameter truth replay", {
  args <- diff_fixture(ngenes = 600L, n1 = 48L, n2 = 48L,
    lBaselineExpr = rep(5, 600), lBaselineExpr2 = rep(20, 600),
    lOD = rep(log(.5), 600), lOD2 = rep(log(.9), 600),
    cts = rep(seq(0, 20, by = 4), each = 8))
  sim <- do.call(diff_simulator, args)
  expect_identical(sim, do.call(diff_simulator, args))
  gt <- sim$ground_truth
  mu1 <- gt$mesor1 + gt$amplitude1 * cos(outer(gt$phase1, sim$times1,
    function(p, t) 2 * pi * (t - p) / 24))
  mu2 <- gt$mesor2 + gt$amplitude2 * cos(outer(gt$phase2, sim$times2,
    function(p, t) 2 * pi * (t - p) / 24))
  z1 <- (sim$expr1 - mu1) / gt$sigma
  z2 <- (sim$expr2 - mu2) / gt$sigma2
  expect_lt(abs(cor(as.vector(z1), as.vector(z2))), .04)
  expect_lt(abs(mean(z1)), .03)
  expect_lt(abs(mean(z2)), .03)
  expect_lt(abs(sd(as.vector(z1)) - 1), .03)
  expect_lt(abs(sd(as.vector(z2)) - 1), .03)
  set.seed(args$sim.seed + args$n1 * 7919L + args$n2 * 104729L)
  backend <- environment(diff_simulator)
  compiled <- exists(".CPP_LOADED", envir = backend, inherits = TRUE) &&
    isTRUE(get(".CPP_LOADED", envir = backend, inherits = TRUE)) &&
    exists("sim_cosinor_expr_fast", envir = backend, mode = "function", inherits = TRUE)
  if (compiled) {
    expected1 <- matrix(rnorm(args$ngenes * args$n1), nrow = args$ngenes, byrow = TRUE)
    expected2 <- matrix(rnorm(args$ngenes * args$n2), nrow = args$ngenes, byrow = TRUE)
  } else {
    expected1 <- expected2 <- matrix(NA_real_, args$ngenes, args$n1)
    for (g in seq_len(args$ngenes)) {
      expected1[g, ] <- rnorm(args$n1)
      expected2[g, ] <- rnorm(args$n2)
    }
  }
  expect_equal(unname(z1), expected1, tolerance = 1e-12)
  expect_equal(unname(z2), expected2, tolerance = 1e-12)
  larger <- args
  larger$n1 <- larger$n2 <- 72L
  next_sim <- do.call(diff_simulator, larger)
  expect_identical(sim$ground_truth, next_sim$ground_truth)
  expect_identical(sim$parameter_draws, next_sim$parameter_draws)
})

test_that("phase truth uses the supplied period and rejects degenerate alternatives", {
  sim <- do.call(diff_simulator, diff_fixture(period = 12,
    cts = rep(seq(0, 10, by = 2), each = 4),
    phase_diff = c(8, 8), dp_shift_mode = "fixed"))
  expect_equal(sim$ground_truth$phase_diff[sim$ground_truth$diff_type == 4],
    rep(-4, sum(sim$ground_truth$diff_type == 4)), tolerance = 1e-12)
  expect_error(do.call(diff_simulator, diff_fixture(phase_diff = c(0, 0),
    dp_shift_mode = "fixed")), "do not match")
  expect_error(do.call(diff_simulator, diff_fixture(mesor_diff = c(0, 0))), "do not match")
})

test_that("all-arrhythmic objects with empty selected vectors remain valid", {
  bio <- diff_options(ngenes = 20, prop_rhythmic = 0,
    prop_DR = 0, prop_DP = 0, prop_DM = 0, lBaselineExpr = 5,
    lOD = log(.5), amplitude = c(.2, .8), sigma_rhythmic = c(.1, .4),
    amplitude2 = c(.3, .9), sigma_rhythmic2 = c(.2, .5), paired_sigma = TRUE)
  sim <- diff_simulator(ngenes = bio$ngenes, n1 = 12, n2 = 12,
    lBaselineExpr = bio$lBaselineExpr, lOD = bio$lOD,
    amplitude = bio$amplitude, sigma_rhythmic = bio$sigma_rhythmic,
    amplitude2 = bio$amplitude2, sigma_rhythmic2 = bio$sigma_rhythmic2,
    prop_rhythmic = 0, prop_DR = 0, prop_DP = 0, prop_DM = 0,
    design = "active")
  expect_true(all(sim$ground_truth$diff_type == 0))
  expect_false(any(sim$ground_truth$is_DR | sim$ground_truth$is_DP | sim$ground_truth$is_DM))
  expect_true(all(is.finite(sim$expr1)))
  expect_true(all(is.finite(sim$expr2)))
})

test_that("runner forwarding preserves both paired arms and full truth", {
  bio <- diff_options(ngenes = 40, prop_rhythmic = .8,
    prop_DR = .4, prop_DP = .2, prop_DM = .1,
    lBaselineExpr = 5, lBaselineExpr2 = 30, lOD = log(.5), lOD2 = log(.9),
    amplitude = c(1, 2), sigma_rhythmic = c(.4, .8),
    amplitude2 = c(2, 3), sigma_rhythmic2 = c(.6, .9), paired_sigma = TRUE)
  design <- CircadianDesignOptions(sample_sizes = c(12, 24), nsims = 1,
    design = "active", cts = seq(0, 20, by = 4), test_types = c("DR", "DP", "DM"))
  analysis <- CircadianAnalysisOptions(alpha = .2, parallel.ncores = 1)
  result <- diff_runner(bio, design, analysis, mc.cores = 1, verbose = FALSE)
  gt <- result$ground_truth[[1]]
  draws <- result$parameter_draws[[1]]
  expect_equal(unname(draws$transform_counts), c(
    sum(draws$raw_amplitude1 < .05, na.rm = TRUE),
    sum(draws$raw_amplitude2 < .05, na.rm = TRUE),
    sum(draws$raw_sigma1 < 1e-6, na.rm = TRUE),
    sum(draws$raw_sigma2 < 1e-6, na.rm = TRUE)))
  g1 <- which(gt$diff_type == 2L)
  g2 <- which(gt$diff_type == 3L)
  expect_equal(gt$amplitude1[g1], bio$amplitude[draws$pair_index1[g1]], tolerance = 0)
  expect_equal(gt$sigma[g1], bio$sigma_rhythmic[draws$pair_index1[g1]], tolerance = 0)
  expect_equal(gt$amplitude2[g2], bio$amplitude2[draws$pair_index2[g2]], tolerance = 0)
  expect_equal(gt$sigma2[g2], bio$sigma_rhythmic2[draws$pair_index2[g2]], tolerance = 0)
  for (endpoint in c("DR", "DP", "DM")) for (j in 1:2)
    expect_equal(result[[paste0("fdr_", endpoint)]][, j, 1],
      p.adjust(result[[paste0("pval_", endpoint)]][, j, 1], "BH"), tolerance = 0)
})

test_that("old objects warn about marginal group2 noise and K2 is separate", {
  bio <- diff_options(ngenes = 20, prop_rhythmic = .8,
    prop_DR = .4, prop_DP = .2, prop_DM = .1, lBaselineExpr = 5,
    lOD = log(.5), amplitude = 1, amplitude2 = 2)
  bio$sigma_rhythmic2 <- bio$paired_sigma2 <- NULL
  design <- CircadianDesignOptions(sample_sizes = 12, nsims = 1,
    design = "active", test_types = "DR")
  analysis <- CircadianAnalysisOptions(alpha = .2)
  expect_warning(result <- diff_runner(bio, design, analysis, mc.cores = 1,
    verbose = FALSE), "no sigma_rhythmic2")
  expect_true(all(is.na(result$parameter_draws[[1]]$pair_index2[
    result$ground_truth[[1]]$diff_type == 3L])))
  bio$paired_2h <- TRUE
  expect_error(diff_runner(bio, design, analysis, verbose = FALSE), "two-harmonic")
})

test_that("zero-amplitude selected pairs retain their indices through the runner", {
  bio <- diff_options(ngenes = 40, prop_rhythmic = .8,
    prop_DR = .4, prop_DP = .2, prop_DM = .1, lBaselineExpr = 5,
    lOD = log(.5), amplitude = c(0, 0), sigma_rhythmic = c(.3, .7),
    amplitude2 = c(0, 0), sigma_rhythmic2 = c(.4, .9), paired_sigma = TRUE)
  design <- CircadianDesignOptions(sample_sizes = 12, nsims = 1,
    design = "active", test_types = "DR")
  analysis <- CircadianAnalysisOptions(alpha = .2)
  result <- diff_runner(bio, design, analysis, mc.cores = 1, verbose = FALSE)
  gt <- result$ground_truth[[1]]
  draws <- result$parameter_draws[[1]]
  g1 <- which(gt$diff_type == 2L)
  g2 <- which(gt$diff_type == 3L)
  expect_equal(draws$raw_amplitude1[g1], bio$amplitude[draws$pair_index1[g1]], tolerance = 0)
  expect_equal(draws$raw_sigma1[g1], bio$sigma_rhythmic[draws$pair_index1[g1]], tolerance = 0)
  expect_equal(draws$raw_amplitude2[g2], bio$amplitude2[draws$pair_index2[g2]], tolerance = 0)
  expect_equal(draws$raw_sigma2[g2], bio$sigma_rhythmic2[draws$pair_index2[g2]], tolerance = 0)
  expect_true(all(gt$amplitude1[g1] == .05))
  expect_true(all(gt$amplitude2[g2] == .05))
  expect_equal(unname(draws$transform_counts[c("amplitude1", "amplitude2")]),
    c(sum(gt$is_rhythmic_g1), sum(gt$is_rhythmic_g2)))
  expect_error(diff_options(ngenes = 40, prop_rhythmic = .8,
    lBaselineExpr = 5, lOD = 0, amplitude = c(-.1, 1),
    sigma_rhythmic = c(.3, .7), paired_sigma = TRUE), "finite matched")
  expect_error(do.call(diff_simulator, diff_fixture(amplitude2 = c(-.1, .7, 2))),
    "nonnegative amplitude2")
})

test_that("invalid designs and inference failures cannot become zero power", {
  bio <- diff_options(ngenes = 20, prop_rhythmic = .8,
    prop_DR = .4, prop_DP = .2, prop_DM = .1, lBaselineExpr = 5,
    lOD = log(.5), amplitude = 1)
  analysis <- CircadianAnalysisOptions(alpha = .2)
  bad_design <- CircadianDesignOptions(sample_sizes = 12, nsims = 1,
    design = "active", cts = rep(0, 12), test_types = "DR")
  expect_error(diff_runner(bio, bad_design, analysis, mc.cores = 1,
    verbose = FALSE), "not identifiable")
  design <- CircadianDesignOptions(sample_sizes = 12, nsims = 1,
    design = "active", test_types = "DR")
  failure_environment <- new.env(parent = environment(diff_runner))
  failure_environment$DCP_Rhythmicity <- function(...) stop("fixture detector failure")
  failing_runner <- diff_runner
  environment(failing_runner) <- failure_environment
  expect_error(failing_runner(bio, design, analysis, mc.cores = 1,
    verbose = FALSE), "DCP pipeline failed")
  failure_environment$DCP_Rhythmicity <- function(...) list(
    rhythm.joint = data.frame(TOJR = rep("both", 20)))
  failure_environment$DCP_DiffR2 <- function(...) data.frame(
    gname = paste0("Gene", 1:20), p.R2 = rep(NA_real_, 20))
  expect_error(failing_runner(bio, design, analysis, mc.cores = 1,
    verbose = FALSE), "Invalid DR p-values")
  failure_environment$DCP_DiffR2 <- function(...) data.frame(
    gname = paste0("Gene", 1:20))
  expect_error(failing_runner(bio, design, analysis, mc.cores = 1,
    verbose = FALSE), "Invalid DCP DR output")
  failure_environment$DCP_Rhythmicity <- function(...) list()
  expect_error(failing_runner(bio, design, analysis, mc.cores = 1,
    verbose = FALSE), "Invalid DCP rhythmicity output")
  failure_environment$simCircadianDiff <- function(...) {
    sim <- diff_simulator(...)
    sim$expr1[1, 1] <- NA_real_
    sim
  }
  expect_error(failing_runner(bio, design, analysis, mc.cores = 1,
    verbose = FALSE), "Invalid simulated expression")
})

test_that("singleton gene pools never select outside their assigned states", {
  configurations <- list(
    list(ngenes = 5L, prop_rhythmic = 1, prop_DR = .8, prop_DP = 0, prop_DM = 0),
    list(ngenes = 5L, prop_rhythmic = 1, prop_DR = .8, prop_DP = .2, prop_DM = 0),
    list(ngenes = 10L, prop_rhythmic = .6, prop_DR = .4, prop_DP = .1, prop_DM = .1))
  for (configuration in configurations) for (seed in 1:20) {
    args <- diff_fixture()
    args[names(configuration)] <- configuration
    args$sim.seed <- seed
    sim <- do.call(diff_simulator, args)
    gt <- sim$ground_truth
    expect_equal(sum(gt$diff_type %in% c(2, 3)), round(args$ngenes * args$prop_DR))
    expect_equal(sum(gt$diff_type == 4), round(args$ngenes * args$prop_DP))
    expect_equal(sum(gt$diff_type == 5), round(args$ngenes * args$prop_DM))
    expect_identical(gt$is_DR, gt$diff_type %in% c(2, 3))
    expect_identical(gt$is_DP, gt$diff_type == 4)
    expect_identical(gt$is_DM, gt$diff_type == 5)
  }
})
