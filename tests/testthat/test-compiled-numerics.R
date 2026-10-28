test_that("compiled cosinor p-values agree with independent QR fits", {
  set.seed(1729)
  times <- sort(stats::runif(41, 0, 24))
  design <- cbind(1, cos(2 * pi * times / 24), sin(2 * pi * times / 24))
  expression <- matrix(stats::rnorm(60 * length(times)), nrow = 60)
  expression[1:20, ] <- expression[1:20, ] +
    outer(seq(0.2, 2, length.out = 20), cos(2 * pi * times / 24))
  fit <- stats::lm.fit(design, t(expression))
  rss <- colSums(fit$residuals^2)
  tss <- rowSums((expression - rowMeans(expression))^2)
  reference <- stats::pf(((tss - rss) / 2) / (rss / (length(times) - 3)),
                         2, length(times) - 3, lower.tail = FALSE)
  expect_equal(as.numeric(SCP:::cosinor_pvals_cpp(expression, times)),
                reference, tolerance = 1e-10)
  expect_true(all(is.na(SCP:::cosinor_pvals_cpp(expression[, 1:3], times[1:3]))))
  expect_identical(as.numeric(SCP:::cosinor_pvals_cpp(expression, rep(0, length(times)))),
                     rep(1, nrow(expression)))
})

test_that("compiled expression simulation agrees with the Gaussian model", {
  mesor <- c(-1, 0, 2, 4)
  amplitude <- c(0, 0.2, 1, 2)
  phase <- c(0, 3, 8, 22)
  sigma <- c(0.2, 0.8, 1.2, 0.5)
  times <- seq(0, 22, by = 2)
  set.seed(2718)
  observed <- SCP:::sim_cosinor_expr_cpp(mesor, amplitude, phase, sigma, times,
                                        alpha2 = 0.3, alpha3 = -0.1)
  set.seed(2718)
  reference <- matrix(NA_real_, nrow = length(mesor), ncol = length(times))
  for (g in seq_along(mesor)) {
    angle <- 2 * pi * (times - phase[g]) / 24
    mu <- mesor[g] + amplitude[g] *
      (cos(angle) + 0.3 * cos(2 * angle) - 0.1 * cos(3 * angle))
    reference[g, ] <- stats::rnorm(length(times), mu, sigma[g])
  }
  expect_equal(observed, reference, tolerance = 1e-13)
})
