#' Finite sample/Large sample Likelihood ratio test for differential amplitude.
#'
#' Test differential amplitude of circadian curve fitting using likelihood ratio test
#' @param tt1 Time vector of condition 1
#' @param yy1 Expression vector of condition 1
#' @param tt2 Time vector of condition 2
#' @param yy2 Expression vector of condition 2
#' @param period Period of the sine curve. Default is 24.
#' @param FN Use finite sample correction (TRUE) or large sample (FALSE)
#'
#' @return A list with:
#'   \item{amp_1}{Amplitude estimate of the 1st data}
#'   \item{amp_2}{Amplitude estimate of the 2nd data}
#'   \item{amp_c}{Amplitude estimate pooling all data together}
#'   \item{l0}{Log likelihood under the null (same amplitude between the two groups)}
#'   \item{la}{Log likelihood under the alternative (different amplitude between the two groups)}
#'   \item{stat}{LR statistics}
#'   \item{pvalue}{P-value from the LR test}
#'
#' @noRd
LRTest_diff_amp <- function(tt1, yy1, tt2, yy2, period = 24, FN = TRUE){
  requireNamespace("nloptr", quietly = TRUE)

  n1 <- length(tt1)
  stopifnot(n1 == length(yy1))
  n2 <- length(tt2)
  stopifnot(length(tt2) == length(yy2))

  w <- 2*pi/period

  fit1 <- fitSinCurve(tt1, yy1, period = period)
  fit2 <- fitSinCurve(tt2, yy2, period = period)

  A1 <- fit1$amp
  A2 <- fit2$amp

  phase1 <- fit1$phase
  phase2 <- fit2$phase

  E1 <- A1 * cos(w * phase1)
  F1 <- A1 * sin(w * phase1)

  E2 <- A2 * cos(w * phase2)
  F2 <- A2 * sin(w * phase2)

  basal1 <- fit1$offset
  basal2 <- fit2$offset

  sigma2_1 <- 1/n1 * fit1$rss
  sigma2_2 <- 1/n2 * fit2$rss

  theta1 <- 1/sigma2_1
  theta2 <- 1/sigma2_2

  p1 <- c(E1, F1, basal1, theta1)
  p2 <- c(E2, F2, basal2, theta2)

  x_Ha <- c(p1, p2)

  asin1 <- sin(w * tt1)
  acos1 <- cos(w * tt1)
  asin2 <- sin(w * tt2)
  acos2 <- cos(w * tt2)

  # Negative log-likelihood function
  eval_f_list <- function(x, asin1, acos1, asin2, acos2) {
    p1 <- x[1:4]
    p2 <- x[5:8]

    E1 <- p1[1]
    F1 <- p1[2]
    basel1 <- p1[3]
    theta1 <- p1[4]
    yhat1 <- E1 * asin1 + F1 * acos1 + basel1

    E2 <- p2[1]
    F2 <- p2[2]
    basel2 <- p2[3]
    theta2 <- p2[4]
    yhat2 <- E2 * asin2 + F2 * acos2 + basel2

    ll1_a <- log(theta1)/2
    ll1_b <- (yy1 - yhat1)^2 * theta1 / 2
    ll1 <- ll1_a - ll1_b

    ll2_a <- log(theta2)/2
    ll2_b <- (yy2 - yhat2)^2 * theta2 / 2
    ll2 <- ll2_a - ll2_b

    partial_E1 <- - theta1 * sum((yy1 - yhat1) * asin1)
    partial_F1 <- - theta1 * sum((yy1 - yhat1) * acos1)
    partial_C1 <- - theta1 * sum(yy1 - yhat1)
    partial_theta1 <- sum((yy1 - yhat1)^2)/2 - n1/2/theta1

    partial_E2 <- - theta2 * sum((yy2 - yhat2) * asin2)
    partial_F2 <- - theta2 * sum((yy2 - yhat2) * acos2)
    partial_C2 <- - theta2 * sum(yy2 - yhat2)
    partial_theta2 <- sum((yy2 - yhat2)^2)/2 - n2/2/theta2

    return( list( "objective" = -sum(ll1) - sum(ll2),
                  "gradient"  = c(partial_E1, partial_F1, partial_C1, partial_theta1,
                                  partial_E2, partial_F2, partial_C2, partial_theta2)
    )
  )
  }

  # Equality constraints: A1^2 = A2^2 (equal amplitudes under H0)
  eval_g_eq <- function(x, asin1, acos1, asin2, acos2)
  {
    p1 <- x[1:4]
    p2 <- x[5:8]

    E1 <- p1[1]
    F1 <- p1[2]
    #basel1 <- p1[3]
    theta1 <- p1[4]

    E2 <- p2[1]
    F2 <- p2[2]
    #basel2 <- p2[3]
    theta2 <- p2[4]

    A2_1 <- (E1^2 + F1^2)
    A2_2 <- (E2^2 + F2^2)
    A2_1 - A2_2
  }

  # Jacobian of equality constraints
  eval_g_eq_jac <- function(x, asin1, acos1, asin2, acos2)
  {
    p1 <- x[1:4]
    p2 <- x[5:8]

    E1 <- p1[1]
    F1 <- p1[2]
    theta1 <- p1[4]

    E2 <- p2[1]
    F2 <- p2[2]
    theta2 <- p2[4]

    A2_1 <- (E1^2 + F1^2)
    A2_2 <- (E2^2 + F2^2)

    c(2 * E1, 2 * F1, 0, 0,
      - 2 * E2, - 2 * F2, 0, 0)
  }

  # Lower and upper bounds
  lb <- c(-Inf,-Inf,-Inf,1e-10, -Inf, -Inf,-Inf, 1e-10)
  ub <- c(Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf)

  # Initial values
  local_opts <- list( "algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-15 )
  opts <- list( "algorithm"= "NLOPT_LD_SLSQP",
                "xtol_rel"= 1.0e-15,
                "maxeval"= 160000,
                "local_opts" = local_opts,
                "print_level" = 0
  )

  res <- nloptr::nloptr ( x0 = x_Ha,
                          eval_f = eval_f_list,
                          lb = lb,
                          ub = ub,
                          eval_g_eq = eval_g_eq,
                          eval_jac_g_eq = eval_g_eq_jac,
                          opts = opts,
                          asin1=asin1,
                          acos1=acos1,
                          asin2=asin2,
                          acos2=acos2)

  x_H0 <- res$solution

  l0 <- - eval_f_list(x_H0, asin1, acos1, asin2, acos2)$objective
  la <- - eval_f_list(x_Ha, asin1, acos1, asin2, acos2)$objective

  LR_stat <- -2*(l0-la)
  if (l0 > la + 1e-6)
    warning("LR_stat < 0 (constrained optimizer may have reached a saddle); clamped to 0")
  LR_stat <- max(0, LR_stat)

  if(!FN){
    pvalue <- pchisq(LR_stat, 1, lower.tail = FALSE)
  } else {
    r <- 1
    k <- 8
    n <- n1+n2
    Fstat <- (exp(LR_stat/n) - 1) * (n-k) / r
    pvalue <- pf(Fstat, df1 = r, df2 = n-k, lower.tail = FALSE)
  }

  amp_c <- sqrt(x_H0[1]^2 + x_H0[2]^2)

  res <- list(amp_1=A1, amp_2=A2, amp_c=amp_c,
              l0=l0,
              la=la,
              stat=-2*(l0-la),
              pvalue=pvalue)
  return(res)
}
