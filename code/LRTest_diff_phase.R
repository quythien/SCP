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

  # Adjust phase for optimizer only — use local copies so originals are returned
  phase1_opt <- phase1
  phase2_opt <- phase2
  if(phase2_opt - phase1_opt > period/2){
    phase2_opt <- phase2_opt - period
  } else if(phase1_opt - phase2_opt > period/2){
    phase1_opt <- phase1_opt - period
  }

  basal1 <- fit1$offset
  basal2 <- fit2$offset

  sigma2_1 <- 1/n1 * fit1$rss
  sigma2_2 <- 1/n2 * fit2$rss

  theta1 <- 1/sigma2_1
  theta2 <- 1/sigma2_2

  p1 <- c(A1, phase1_opt, basal1, theta1)
  p2 <- c(A2, phase2_opt, basal2, theta2)
  x_Ha <- c(p1, p2)

  # Negative log-likelihood function
  eval_f_list <- function(x) {
    p1 <- x[1:4]
    p2 <- x[5:8]
    A1 <- p1[1]
    phase1 <- p1[2]
    basel1 <- p1[3]
    theta1 <- p1[4]

    asin1 <- sin(w * (tt1 + phase1) )
    acos1 <- cos(w * (tt1 + phase1) )
    yhat1 <- A1 * asin1 + basel1

    A2 <- p2[1]
    phase2 <- p2[2]
    basel2 <- p2[3]
    theta2 <- p2[4]

    asin2 <- sin(w * (tt2 + phase2) )
    acos2 <- cos(w * (tt2 + phase2) )
    yhat2 <- A2 * asin2 + basel2

    ll1_a <- log(theta1)/2
    ll1_b <- (yy1 - yhat1)^2 * theta1 / 2
    ll1 <- ll1_a - ll1_b

    ll2_a <- log(theta2)/2
    ll2_b <- (yy2 - yhat2)^2 * theta2 / 2
    ll2 <- ll2_a - ll2_b

    partial_A1 <- - theta1 * sum((yy1 - yhat1) * asin1)
    partial_phase1 <- - theta1 * A1 * w * sum((yy1 - yhat1) * acos1)
    partial_C1 <- - theta1 * sum(yy1 - yhat1)
    partial_theta1 <- sum((yy1 - yhat1)^2)/2 - n1/2/theta1

    partial_A2 <- - theta2 * sum((yy2 - yhat2) * asin2)
    partial_phase2 <- - theta2 * A2 * w * sum((yy2 - yhat2) * acos2)
    partial_C2 <- - theta2 * sum(yy2 - yhat2)
    partial_theta2 <- sum((yy2 - yhat2)^2)/2 - n2/2/theta2

    return( list( "objective" = -sum(ll1) - sum(ll2),
    "gradient" = c(partial_A1, partial_phase1, partial_C1, partial_theta1,
    partial_A2, partial_phase2, partial_C2, partial_theta2)
    )
    )
  }

  # Equality constraints: phase1 = phase2
  eval_g_eq <- function(x)
  {
    p1 <- x[1:4]
    p2 <- x[5:8]
    phase1 <- p1[2]
    phase2 <- p2[2]
    phase1 - phase2
  }

  # Jacobian of equality constraints
  eval_g_eq_jac <- function(x)
  {
    c(0, 1, 0, 0,
    0, -1, 0, 0)
  }

  # Lower and upper bounds
  lb <- c(0,-Inf,-Inf,1e-10, 0, -Inf,-Inf, 1e-10)
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
  opts = opts)

  x_H0 <- res$solution
  l0 <- - eval_f_list(x_H0)$objective
  la <- - eval_f_list(x_Ha)$objective

  LR_stat <- -2*(l0-la)
  if (l0 > la + 1e-6)
    warning("LR_stat < 0 (constrained optimizer may have reached a saddle); clamped to 0")
  LR_stat <- max(0, LR_stat)

  if(!FN){
    pvalue <- pchisq(LR_stat,1,lower.tail = FALSE)
  } else if(FN){
    r <- 1
    k <- 8
    n <- n1+n2
    Fstat <- (exp(LR_stat/n) - 1) * (n-k) / r
    pvalue <- pf(Fstat,df1 = r, df2 = n-k, lower.tail = FALSE)
  } else{
    stop("FN has to be TRUE or FALSE")
  }

  phase_c <- x_H0[2]
  phase_c2 <- x_H0[6]

  res <- list(phase_1=phase1, phase_2=phase2, phase_c=phase_c,
  l0=l0,
  la=la,
  stat=-2*(l0-la),
  pvalue=pvalue)

  return(res)
}
