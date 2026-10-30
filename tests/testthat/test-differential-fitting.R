# Opposite residuals at each time make their contribution orthogonal to every
# fitted time function. The sinusoidal coefficients and residual sum of squares
# are therefore known without fitting a second model.
paired_time_curve <- function(times, amplitude, phase, offset, errors, period=24) {
  tt <- rep(times,each=2)
  residual <- as.vector(rbind(errors,-errors))
  yy <- offset+amplitude*sin(2*pi*(tt+phase)/period)+residual
  list(time=tt,value=yy,rss=2*sum(errors^2))
}

test_that("sinusoidal fitting recovers known coefficients on irregular times", {
  period <- 24
  amplitude <- 1.7
  phase <- 22.4
  offset <- -2.3
  x <- paired_time_curve(c(.4,2.1,4.8,8.3,13.7,18.2,23.6),
    amplitude,phase,offset,c(.2,.4,.3,.7,.5,.6,.8),period)
  fit <- SCP:::fitSinCurve(x$time,x$value,period)
  expect_equal(fit$amp,amplitude,tolerance=1e-12)
  expect_equal(fit$phase,phase,tolerance=1e-12)
  expect_equal(fit$offset,offset,tolerance=1e-12)
  expect_equal(fit$A,amplitude*cos(2*pi*phase/period),tolerance=1e-12)
  expect_equal(fit$B,amplitude*sin(2*pi*phase/period),tolerance=1e-12)
  expect_equal(fit$rss,x$rss,tolerance=1e-12)
  expect_equal(fit$peak,7.6,tolerance=1e-12)
})

test_that("least-squares estimates do not depend on starting values", {
  x <- paired_time_curve(c(.4,2.1,4.8,8.3,13.7,18.2,23.6),
    1.7,22.4,-2.3,c(.2,.4,.3,.7,.5,.6,.8))
  reference <- SCP:::fitSinCurve(x$time,x$value)
  starts <- list(list(amp=0,phase=48,offset=1e6),
    list(amp=100,phase=-36,offset=-1e6),list(amp=.001,phase=11,offset=0))
  for(start in starts)
    expect_equal(SCP:::fitSinCurve(x$time,x$value,parStart=start),reference,tolerance=1e-13)
})

test_that("a known common-phase null has zero likelihood ratio", {
  phase <- 23.8
  x <- paired_time_curve(c(.4,2.1,4.8,8.3,13.7,18.2,23.6),
    1.4,phase,2,c(.2,.4,.3,.7,.5,.6,.8))
  y <- paired_time_curve(c(.2,3.3,6.2,10.8,12.1,16.4,20.3,23.1),
    .7,phase,-4,c(.6,.3,.5,.4,.2,.8,.7,.9))
  n1 <- length(x$time)
  n2 <- length(y$time)
  log_likelihood <- -n1/2*(1+log(x$rss/n1))-n2/2*(1+log(y$rss/n2))
  for(finite_sample in c(TRUE,FALSE)) {
    fit <- SCP:::LRTest_diff_phase(x$time,x$value,y$time,y$value,FN=finite_sample)
    expect_equal(fit$la,log_likelihood,tolerance=1e-10)
    expect_equal(fit$l0,log_likelihood,tolerance=1e-10)
    expect_lt(abs(fit$stat),1e-9)
    expect_equal(fit$pvalue,1,tolerance=1e-6)
    expect_equal(((fit$phase_c-phase+12)%%24)-12,0,tolerance=1e-6)
  }
})

test_that("the phase test finds an analytic zero-amplitude boundary optimum", {
  # With quadrature phases and balanced times, the null RSS product is a
  # concave quadratic in sin(phase)^2. Its minimum is at an endpoint.
  # The stronger first group fixes phase zero and the second amplitude is zero.
  times <- seq(0,21,by=3)
  x <- paired_time_curve(times,1.2,0,2,rep(.3,8))
  y <- paired_time_curve(times,.8,6,-1,rep(.3,8))
  n <- length(times)*2
  expected_LR <- n*log1p((n*.8^2/2)/y$rss)
  fit <- SCP:::LRTest_diff_phase(x$time,x$value,y$time,y$value,FN=FALSE)
  expect_equal(fit$stat,expected_LR,tolerance=1e-8)
  expect_equal(fit$pvalue,stats::pchisq(expected_LR,df=1,lower.tail=FALSE),tolerance=1e-10)
  expect_equal(((fit$phase_c+12)%%24)-12,0,tolerance=1e-6)
})

test_that("phase inference is invariant to group order, units and positive scaling", {
  x <- paired_time_curve(c(.4,2.1,4.8,8.3,13.7,18.2,23.6),
    1.4,22.8,2,c(.2,.4,.3,.7,.5,.6,.8))
  y <- paired_time_curve(c(.2,3.3,6.2,10.8,12.1,16.4,20.3,23.1),
    .9,3.2,-4,c(.6,.3,.5,.4,.2,.8,.7,.9))
  fit <- SCP:::LRTest_diff_phase(x$time,x$value,y$time,y$value)
  transformed <- list(
    swapped=SCP:::LRTest_diff_phase(y$time,y$value,x$time,x$value),
    minutes=SCP:::LRTest_diff_phase(60*x$time,x$value,60*y$time,y$value,period=1440),
    scaled=SCP:::LRTest_diff_phase(x$time,7*x$value,y$time,.2*y$value),
    translated=SCP:::LRTest_diff_phase(x$time,x$value+17,y$time,y$value-31))
  for(result in transformed) {
    expect_equal(result$stat,fit$stat,tolerance=1e-8)
    expect_equal(result$pvalue,fit$pvalue,tolerance=1e-8)
  }
  expect_equal(transformed$swapped$phase_1,fit$phase_2,tolerance=1e-10)
  expect_equal(transformed$swapped$phase_2,fit$phase_1,tolerance=1e-10)
  expect_equal(transformed$minutes$phase_1/60,fit$phase_1,tolerance=1e-10)
  expect_equal(transformed$minutes$phase_2/60,fit$phase_2,tolerance=1e-10)
})

test_that("unidentifiable designs and undefined phase likelihoods are rejected", {
  expect_error(SCP:::fitSinCurve(rep(0,10),seq_len(10)),"rank deficient")
  expect_error(SCP:::LRTest_diff_phase(rep(0,10),seq_len(10),0:9,seq_len(10)),
    "rank deficient")
  expect_error(SCP:::LRTest_diff_phase(0:9,rep(0,10),0:9,seq_len(10)),
    "positive residual variance")
  expect_error(SCP:::LRTest_diff_phase(c(0,5,13,19),c(1,4,0,2),
    c(0,5,13,19),c(1,4,0,2)),"residual degrees of freedom")
  expect_error(SCP:::fitSinCurve(c(NA,1:9),seq_len(10)))
  expect_error(SCP:::LRTest_diff_phase(0:9,seq_len(10),0:9,seq_len(10),FN=NA),
    "FN has to be TRUE or FALSE")
})
