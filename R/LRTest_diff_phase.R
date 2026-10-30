#' Finite sample/Large sample Likelihood ratio test for differential phase.
#'
#' Test differential phase of circadian curve fitting using likelihood ratio test
#' @title Likelihood ratio test for detecting differential phase.
#' @param tt1 Time vector of condition 1
#' @param yy1 Expression vector of condition 1
#' @param tt2 Time vector of condition 2
#' @param yy2 Expression vector of condition 2
#' @param period Period of the sine curve. Default is 24.
#' @param FN Use finite sample correction (TRUE) or large sample (FALSE)
#'
#' @return A list with:
#'   \item{phase_1}{Phase estimate of the 1st data}
#'   \item{phase_2}{Phase estimate of the 2nd data}
#'   \item{phase_c}{Phase estimate pooling all data together}
#'   \item{l0}{Log likelihood under the null (same phase between the two groups)}
#'   \item{la}{Log likelihood under the alternative (different phase between the two groups)}
#'   \item{stat}{LR statistics}
#'   \item{pvalue}{P-value from the LR test}
#'
#' @noRd
LRTest_diff_phase <- function(tt1, yy1, tt2, yy2, period = 24, FN = TRUE){
  fit1 <- fitSinCurve(tt1, yy1, period = period)
  fit2 <- fitSinCurve(tt2, yy2, period = period)
  n1 <- length(tt1)
  n2 <- length(tt2)
  n <- n1 + n2
  if (length(FN) != 1L || is.na(FN) || !FN %in% c(TRUE, FALSE))
    stop("FN has to be TRUE or FALSE", call. = FALSE)
  if (n1 <= 3L || n2 <= 3L || (FN && n <= 8L) ||
      fit1$rss <= 0 || fit2$rss <= 0)
    stop("The phase test requires residual degrees of freedom and positive residual variance.",
      call. = FALSE)
  w <- 2 * pi / period

  # Center each group once for the one-dimensional common-phase profile.
  groups <- lapply(list(list(t = tt1, y = yy1), list(t = tt2, y = yy2)), function(g) {
    y <- g$y - mean(g$y)
    s <- sin(w * g$t)
    c <- cos(w * g$t)
    s <- s - mean(s)
    c <- c - mean(c)
    list(n = length(y), y = y, s = s, c = c, yy = sum(y^2),
      ss = sum(s^2), cc = sum(c^2), sc = sum(s * c),
      ys = sum(y * s), yc = sum(y * c))
  })
  profile <- function(phi, direct = FALSE) {
    cv <- cos(w * phi)
    sv <- sin(w * phi)
    value <- numeric(length(phi))
    for (g in groups) {
      q <- g$ss * cv^2 + g$cc * sv^2 + 2 * g$sc * cv * sv
      h <- g$ys * cv + g$yc * sv
      amplitude <- pmax(0, h / q)
      rss <- g$yy - pmax(0, h)^2 / q
      recalculate <- if (direct) seq_along(phi) else
        which(rss < 128 * .Machine$double.eps * max(1, g$yy))
      for (j in recalculate)
        rss[j] <- sum((g$y - amplitude[j] * (cv[j] * g$s + sv[j] * g$c))^2)
      if (any(!is.finite(rss)) || any(rss <= 0))
        stop("The phase profile requires positive residual variance.", call. = FALSE)
      value <- value + g$n / 2 * (1 + log(rss / g$n))
    }
    value
  }

  # Search the circle and refine every detected minimum.
  points <- 4096L
  step <- period / points
  grid <- seq(0, period, length.out = points + 1L)[seq_len(points)]
  values <- profile(grid)
  left <- values[c(points, seq_len(points - 1L))]
  right <- values[c(2:points, 1L)]
  minima <- which((values < left & values <= right) | (values <= left & values < right))
  phases <- c(fit1$phase, fit2$phase)
  boundary_phases <- vapply(groups, function(g) atan2(g$yc, g$ys) / w, numeric(1))
  starts <- unique(c(grid[c(minima, which.min(values))], phases,
    (boundary_phases - period / 4) %% period, (boundary_phases + period / 4) %% period))
  refined <- vapply(starts, function(phi)
    stats::optimize(profile, c(phi - step, phi + step), tol = period * 1e-11)$minimum,
    numeric(1))
  candidates <- c(starts, refined)
  final_values <- profile(candidates, direct = TRUE)
  best <- which.min(final_values)
  phase_c <- candidates[best] %% period
  l0 <- -final_values[best]
  la <- -n1 / 2 * (1 + log(fit1$rss / n1)) - n2 / 2 * (1 + log(fit2$rss / n2))

  # A null model cannot improve on the unrestricted least-squares fit.
  tolerance <- 128 * .Machine$double.eps * max(1, abs(l0), abs(la))
  if (l0 > la + tolerance)
    stop("The common-phase likelihood exceeds the unrestricted likelihood.", call. = FALSE)
  LR_stat <- max(0, 2 * (la - l0))
  pvalue <- if (FN) {
    Fstat <- expm1(LR_stat / n) * (n - 8)
    stats::pf(Fstat, df1 = 1, df2 = n - 8, lower.tail = FALSE)
  } else {
    stats::pchisq(LR_stat, df = 1, lower.tail = FALSE)
  }
  list(phase_1 = fit1$phase, phase_2 = fit2$phase, phase_c = phase_c,
    l0 = l0, la = la, stat = LR_stat, pvalue = pvalue)
}
