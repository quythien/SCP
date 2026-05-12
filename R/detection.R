#' Differential Circadian Analysis for PowerSim
#'
#' This file provides differential circadian rhythm analysis methods
#' Based on the DiffCircadian package implementation
#'
#' Main Functions:
#'   DCP_Rhythmicity() - Single/two-group rhythmicity analysis with TOJR
#'   DCP_DiffR2()      - Differential rhythmicity test (DR)
#'   DCP_DiffPar()     - Differential parameter tests (DP, DA)
#'   DCP_Analyze()     - Complete differential analysis workflow
#'
#' Author: Based on DiffCircadian package
#' Date: 2025-02-12

# ==============================================================================
# DEPENDENCIES
# ==============================================================================

# Required packages - will be loaded when functions are called
# - limma (linear models)
# - nloptr (non-linear optimization)
# - minpack.lm (Levenberg-Marquardt)
# - Matrix (matrix operations)
# - rootSolve (root finding)
# - lubridate (time zones)

# LRTest_diff_phase.R and LRTest_diff_amp.R are loaded as part of the package;
# the runtime source() calls used in the legacy code/detection.R are not needed here.

# ==============================================================================
# SECTION 1: FITTING FUNCTIONS
# ==============================================================================

#' Fit Sinusoidal Curve
#'
#' Fit a sine curve where tt is time, and yy is expression value.
#'
#' @param tt Time vector
#' @param yy Expression vector
#' @param period Period of the sine curve (default 24)
#' @param parStart Initial value for optimization (default list(amp=3,phase=0,offset=0))
#'
#' @return A list with components:
#'   \item{amp}{Amplitude}
#'   \item{phase}{Phase (0 to period)}
#'   \item{offset}{Basal level}
#'   \item{peak}{Peak time}
#'   \item{A}{A coefficient (sin component)}
#'   \item{B}{B coefficient (cos component)}
#'   \item{tss}{Total sum of squares}
#'   \item{rss}{Residual sum of squares}
#'   \item{R2}{R-squared}
#'
#' @examples
#' set.seed(32608)
#' n <- 50
#' tt <- runif(n,0,24)
#' Amp <- 2
#' Phase <- 6
#' Offset <- 3
#' yy <- Amp * sin(2*pi/24 * (tt + Phase)) + Offset + rnorm(n,0,1)
#' fitSinCurve(tt, yy)
#'
#' @export
fitSinCurve <- function(tt, yy, period = 24, parStart = list(amp=3,phase=0,offset=0)){

  getPred <- function(parS, tt) {
    parS$amp * sin(2*pi/period * (tt + parS$phase)) + parS$offset
  }

  residFun <- function(p, yy, tt) yy - getPred(p,tt)

  nls.out <- minpack.lm::nls.lm(par=parStart, fn = residFun, yy = yy, tt = tt)

  apar <- nls.out$par

  amp0 <- apar$amp
  asign <- sign(amp0)
  amp <- amp0 * asign

  phase0 <- apar$phase
  phase <- (phase0 + ifelse(asign==1,0,period/2)) %% period
  offset <- apar$offset

  peak <- (period/2 * sign(amp0) - period/4 - phase) %%period
  if(peak > period/4*3) peak = peak - period

  A <- amp0 * cos(2*pi/period * phase0)
  B <- amp0 * sin(2*pi/period * phase0)

  rss <- sum(nls.out$fvec^2)
  tss <- sum((yy - mean(yy))^2)
  R2 <- 1 - rss/tss

  res <- list(amp=amp, phase=phase, offset=offset, peak=peak, A=A, B=B, tss=tss, rss=rss, R2=R2)
  res
}

#' Single Group Rhythmicity Analysis (Internal)
#'
#' @param x1 Group 1 data (list with data, time, gname components)
#' @param period Period (default 24)
#' @param alpha Significance threshold
#' @param CI Calculate confidence intervals?
#' @param p.adjust.method P-value adjustment method
#'
#' @return List with rhythmicity results
CP_OneGroup = function(x1, period = 24, alpha = 0.05, CI = FALSE, p.adjust.method = "BH"){
  stopifnot("x1$data must be dataframe" = is.data.frame(x1$data))
  stopifnot("Number of samples in data does not match that in time. " = ncol(x1$data)==length(x1$time))
  stopifnot("Please input the gene labels x1$gname. " = !is.null(x1$gname))
  stopifnot("Number of gnames does not match number of genes in data. " = nrow(x1$data)==length(x1$gname))
  data = x1$data
  time = x1$time
  gname = x1$gname

  #design matrix
  design.vars = data.frame(inphase = cos(2 * pi * time / period),
                           outphase = sin(2 * pi * time / period))
  design = stats::model.matrix(~inphase+outphase, data = design.vars)
  fit = limma::lmFit(data, design)
  fit = limma::eBayes(fit, trend = TRUE, robust = TRUE)
  top = limma::topTable(fit, coef = 2:3, n = nrow(data), sort.by = "none")
  m.top = limma::topTable(fit, coef = 1, n = nrow(data), sort.by = "none")
  all.top = limma::topTable(fit, coef = 1:3, n = nrow(data), sort.by = "none")

  A.hat = apply(top[, 1:2], 1, function(i){sqrt(i[1]^2 + i[2]^2)})
  phase.hat = apply(top[, 1:2], 1, function(i){get_phase(i[1], i[2])$phase})
  RSS = fit$sigma^2
  TSS = apply(data, 1, function(i){sum((i-mean(i))^2)})

  R2 = 1-RSS*(length(time)-3)/TSS

  if(CI){
    XX_inv <- tryCatch(solve(t(design) %*% design),
                       error = function(e) NULL)
    if (is.null(XX_inv)) {
      warning("Design matrix is singular; CIs cannot be computed.")
      CI <- FALSE
    }
  }
  if(CI){
    CI.m.hat.radius = sapply(seq_along(fit$sigma), function(i){
      calculate_CI.M(XX_inv, A.t = matrix(c(1, 0, 0), nrow = 1),
                     r.full = 3, ncol(data), alpha, fit$sigma[i])
    })
    se.hat.A.phase = t(sapply(seq_along(fit$sigma), function(i){
      calculate_CI_A.phase.Taylor(XX_inv,
                                  A.t = rbind(c(0, 1, 0),
                                              c(0, 0, 1)),
                                  phase.hat[i], A.hat[i], fit$sigma[i])
    }))
    se.A.hat = unlist(se.hat.A.phase[,"se.A.hat"]); se.phase.hat = unlist(se.hat.A.phase[, "se.phase.hat"])
    se.phase.hat = ifelse(se.phase.hat>(pi/2), pi/2, se.phase.hat)
    se.qt = stats::qt(1-alpha/2, ncol(data)-3)

    rhythm = data.frame(gname = gname,
                        M = m.top$logFC,
                        A = A.hat,
                        phase = phase.hat,
                        peak = (period-period*phase.hat/(2*pi))%%period,
                        M.ll = m.top$logFC-CI.m.hat.radius, M.ul = m.top$logFC+CI.m.hat.radius,
                        A.ll = A.hat-se.qt*se.A.hat, A.ul = A.hat+se.qt*se.A.hat,
                        phase.ll = phase.hat - se.qt*se.phase.hat, phase.ul = phase.hat + se.qt*se.phase.hat,
                        pvalue = top$P.Value,
                        qvalue = top$adj.P.Val,
                        sigma = fit$sigma, R2 = R2)

  }else{
    rhythm = data.frame(gname = gname,
                        M = m.top$logFC,
                        A = A.hat,
                        phase = phase.hat,
                        peak = (period-period*phase.hat/(2*pi))%%period,
                        pvalue = top$P.Value,
                        qvalue = top$adj.P.Val,
                        sigma = fit$sigma, R2 = R2)

  }
  x1$rhythm = rhythm
  x1$P = period
  return(x1)
}

# ==============================================================================
# SECTION 2: JOINT RHYTHMICITY (TOJR)
# ==============================================================================

#' Rhythmicity Analysis with Cosinor Model
#'
#' This function either takes single-group data and performs only rhythmicity analysis,
#' or takes two-group data and also categorize the genes into types of joint rhythmicity (TOJR).
#'
#' @param x1 Group 1 data (list with data, time, gname components)
#' @param x2 Group 2 data (optional, same format as x1)
#' @param method Algorithm for joint rhythmicity categorization. Should be one of
#'   "Sidak_FS" (default, recommended), "Sidak_BS", "VDA", "AWFisher"
#' @param period Oscillation cycle length (default 24 for circadian)
#' @param amp.cutoff Only genes with amplitude > amp.cutoff are considered rhythmic
#' @param alpha Threshold for rhythmicity p-value
#' @param alpha.FDR Threshold for FDR-adjusted rhythmicity p-value
#' @param CI Calculate confidence intervals?
#' @param p.adjust.method P-value adjustment method
#' @param parallel.ncores Number of cores for parallel processing
#'
#' @return List with:
#'   \item{x1, x2}{Input data with added rhythmicity analysis}
#'   \item{gname_overlap}{Overlapping gene names}
#'   \item{rhythm.joint}{Data frame with joint rhythmicity categories}
#'
#' @details The methods "Sidak_FS" and "Sidak_BS" implement selective sequential model selection
#' with Sidak adjusted p-value. "FS" represents forward stop, and "BS" basic stop.
#' The method "Sidak_FS" has better type I error control compared to Venn diagram
#' analysis (VDA) and adaptively weighted Fisher's method (AWFisher).
#'
#' @references
#' Fithian, W., Taylor, J., Tibshirani, R., & Tibshirani, R. (2015). Selective
#' sequential model selection. arXiv preprint arXiv:1512.02565.
#'
#' @examples
#' x = DCP_sim_data(ngene=1000, nsample=30, A1=c(1, 3), A2=c(1, 3),
#' phase1=c(0, pi/4), phase2=c(pi/4, pi/2),
#' M1=c(4, 6), M2=c(4, 6), sigma1=1, sigma2=1)
#' rhythm.res = DCP_Rhythmicity(x1 = x[[1]], x2 = x[[2]])
#'
#' @export
DCP_Rhythmicity = function(x1, x2=NULL, method = "Sidak_FS", period = 24, amp.cutoff = 0,
                           alpha = 0.05, alpha.FDR = 0.05, CI = FALSE,
                           p.adjust.method = "BH", parallel.ncores  = 1){

  x1 = CP_OneGroup(x1, period, alpha, CI, p.adjust.method)

  if(is.null(x2)){
    return(x1)
  }else{
    gname.overlap = intersect(x1$gname, x2$gname)
    stopifnot("There are no overlapping genes between x1$gname and x2$gname. " = length(gname.overlap)>0)

    x2 = CP_OneGroup(x2, period, alpha, CI, p.adjust.method)

    pM = data.frame(pG1 = x1$rhythm[match(gname.overlap, x1$rhythm$gname), "pvalue"],
                    pG2 = x2$rhythm[match(gname.overlap, x2$rhythm$gname), "pvalue"])
    action1 = apply(pM, 1, which.min)
    action2 = ifelse(action1==1, 2, 1)
    action = data.frame(action1, action2)
    pv = data.frame(pS1 = sapply(seq_along(action1), function(i){pM[i, action1[i]]}),
                    pS2 = sapply(seq_along(action2), function(i){pM[i, action2[i]]}))

    x = list(x1 = x1, x2 = x2,
             gname_overlap = gname.overlap,
             rhythm.joint = cbind.data.frame(gname.overlap, action, pM))
    colnames(x$rhythm.joint) = c("gname", "action1", "action2", "pG1", "pG2")

    x$rhythm.joint$TOJR = toTOJR(x, method, amp.cutoff, alpha, adjustP = FALSE,
                                p.adjust.method, parallel.ncores)
    x$rhythm.joint$TOJR.FDR = toTOJR(x, method, amp.cutoff, alpha.FDR, adjustP = TRUE,
                                  p.adjust.method, parallel.ncores)
    return(x)
  }
}

#' Types of Joint Rhythmicity (TOJR)
#'
#' Categorize genes into four types of joint rhythmicity:
#' arrhy, rhyI (rhythmic in group 1 only), rhyII (rhythmic in group 2 only), both
#'
#' @param x Output from DCP_Rhythmicity(x1, x2) with rhythm.joint component
#' @param method Algorithm for categorization: "Sidak_FS", "Sidak_BS", "VDA", "AWFisher"
#' @param amp.cutoff Amplitude cutoff for considering genes rhythmic
#' @param alpha Threshold for rhythmicity p-value
#' @param adjustP Use adjusted p-values?
#' @param p.adjust.method P-value adjustment method
#' @param parallel.ncores Number of cores for parallel processing
#'
#' @return Vector of TOJR categories
#'
#' @examples
#' # Re-calculate TOJR with q-value cutoff 0.1
#' x = DCP_sim_data(ngene=1000, nsample=30, A1=c(1, 3), A2=c(1, 3),
#' phase1=c(0, pi/4), phase2=c(pi/4, pi/2),
#' M1=c(4, 6), M2=c(4, 6), sigma1=1, sigma2=1)
#' rhythm.res = DCP_Rhythmicity(x1 = x[[1]], x2 = x[[2]])
#' TOJR.new = toTOJR(rhythm.res, alpha = 0.1, adjustP = TRUE)
#'
#' @export
toTOJR = function(x, method = "Sidak_FS", amp.cutoff = 0, alpha = 0.05,
                 adjustP = TRUE, p.adjust.method = "BH", parallel.ncores = 1){

  if(is.null(x$rhythm.joint)){
    stopifnot("Please see examples for correct x input" =
                (length(x)==2)&(!is.null(x[[1]]$rhythm))&(!is.null(x[[2]]$rhythm)))
    x1 = x[[1]]
    x2 = x[[2]]
    gname.overlap = intersect(x1$gname, x2$gname)

    pM = data.frame(pG1 = x1$rhythm[match(gname.overlap, x1$rhythm$gname), "pvalue"],
                    pG2 = x2$rhythm[match(gname.overlap, x2$rhythm$gname), "pvalue"])
    action1 = apply(pM, 1, which.min)
    action2 = ifelse(action1==1, 2, 1)
    action = data.frame(action1, action2)
    pv = data.frame(pS1 = sapply(seq_along(action1), function(i){pM[i, action1[i]]}),
                    pS2 = sapply(seq_along(action2), function(i){pM[i, action2[i]]}))

  }else{
    gname.overlap = x$rhythm.joint$gname
    pM = as.data.frame(x$rhythm.joint[, c("pG1", "pG2")])
    action1 = x$rhythm.joint$action1
    action2 = x$rhythm.joint$action2
    action = data.frame(action1, action2)
    pv = data.frame(pS1 = sapply(seq_along(action1), function(i){pM[i, action1[i]]}),
                    pS2 = sapply(seq_along(action2), function(i){pM[i, action2[i]]}))
  }

  if(adjustP){
    qM = data.frame(p1 = stats::p.adjust(pM$pG1, p.adjust.method),
                    p2 = stats::p.adjust(pM$pG2, p.adjust.method))
    q.action1 = apply(qM, 1, which.min)
    q.action2 = ifelse(q.action1==1, 2, 1)
    qv = data.frame(qS1 = sapply(seq_along(q.action1), function(i){qM[i, q.action1[i]]}),
                    qS2 = sapply(seq_along(q.action2), function(i){qM[i, q.action2[i]]}))
    q.action = data.frame(q.action1, q.action2)
    TOJR_adj = unlist(parallel::mclapply(seq_len(nrow(q.action)), function(i){
      SeqModelSel(q.action[i, ], qv[i, ], alpha, method)
    }, mc.cores = parallel.ncores))
  }else{
    TOJR_adj = unlist(parallel::mclapply(seq_len(nrow(action)), function(i){
      SeqModelSel(action[i, ], pv[i, ], alpha, method)
    }, mc.cores = parallel.ncores))
  }

  if(amp.cutoff!=0){
    amp0.genes1 = x$x1$rhythm$gname[x$x1$rhythm$A<amp.cutoff]
    amp0.genes2 = x$x2$rhythm$gname[x$x2$rhythm$A<amp.cutoff]
    amp0.genes = unique(c(amp0.genes1, amp0.genes2))
    amp0.genes.not.arrhy = amp0.genes[amp0.genes%in% gname.overlap[TOJR_adj!="arrhy"]]
    amp0.status = lapply(amp0.genes.not.arrhy, function(a.gene){
      if((a.gene %in% amp0.genes1)&(a.gene %in% amp0.genes2)){
        return(c(FALSE, FALSE))
      }else if((a.gene %in% amp0.genes1)&(!(a.gene %in% amp0.genes2))){
        return(c(FALSE, TRUE))
      }else{
        return(c(TRUE, FALSE))
      }
    })

    TOJR.to.change = TOJR_adj[match(amp0.genes.not.arrhy, gname.overlap)]
    if(length(amp0.status)>0){
      TOJR_adj[match(amp0.genes.not.arrhy, gname.overlap)] = sapply(seq_along(amp0.status), function(a){
        amp.cut(amp0.status[[a]], TOJR.to.change[a])
      })
    }
  }

  return(TOJR_adj)
}

# ==============================================================================
# SECTION 3: MODEL SELECTION FUNCTIONS
# ==============================================================================

#' Sequential Model Selection for Joint Rhythmicity
#'
#' @param action Vector indicating which group to test first
#' @param pv P-values for the tests
#' @param alpha Significance threshold
#' @param method Selection method: "Fisher_BS", "Fisher_FS", "Sidak_BS", "Sidak_FS",
#'   "Nominal_BS", "Nominal_FS", "VDA", "AWFisher"
#'
#' @return Joint rhythmicity category
SeqModelSel = function(action = c(1, 2), pv = c(0.01, 0.02), alpha = 0.05, method = "Fisher_FS"){
  action = as.numeric(action)
  pv = as.numeric(pv)

  stop.vda = function(nsel, action){
    if(nsel == 2){
      a.stop = "both"
    }else if(nsel == 0){
      a.stop = "arrhy"
    }else if(nsel == 1){
      sel = action[1]
      if(sel == 1){
        a.stop = "rhyI"
      }else{
        a.stop = "rhyII"
      }
    }
    return(a.stop)
  }

  if(method == "Fisher_BS"){
    p_combined = fishers.p(pv)
    stop = StopRule(action, pv, alpha, "BS")
    if(p_combined<alpha){
      stop= ifelse(stop$type=="arrhy", stop$one, stop$type)
    }else{
      stop = "arrhy"
    }
  }else if(method == "Fisher_FS"){
    p_combined = fishers.p(pv)
    stop = StopRule(action, pv, alpha, "FS")
    if(p_combined<alpha){
      stop= ifelse(stop$type=="arrhy", stop$one, stop$type)
    }else{
      stop = "arrhy"
    }
  }else if(method == "Sidak_BS"){
    stop = StopRule(action, pv, alpha.adjust(alpha, length(pv), "sidak"), "BS")$type
  }else if(method == "Sidak_FS"){
    stop = StopRule(action, pv, alpha.adjust(alpha, length(pv), "sidak"), "FS")$type
  }else if(method == "Nominal_BS"){
    stop = StopRule(action, pv, alpha, "BS")$type
  }else if(method == "Nominal_FS"){
    stop = StopRule(action, pv, alpha, "FS")$type
  }else if(method == "VDA"){
    is.rhy = pv<alpha
    if(sum(is.rhy)==2){
      stop = stop.vda(2, action)
    }else if(sum(is.rhy)==0){
      stop = stop.vda(0, action)
    }else{
      stop = stop.vda(1, action)
    }
  }else if(method == "AWFisher"){
    p_combined = AWFisher::AWFisher_pvalue(as.numeric(pv))
    aw.wight = p_combined$weights
    p_combined = p_combined$pvalues
    if(p_combined>alpha){
      stop = stop.vda(0, action)
    }else if(sum(aw.wight)==2){
      stop = stop.vda(2, action)
    }else if(aw.wight[1, 1]==1){
      stop = stop.vda(1, action)
    }else{
      stop = stop.vda(1, action)
    }
  }
  return(stop)
}

#' Stop Rule for Sequential Testing
#'
#' @param action Vector indicating which group to test first
#' @param pv P-values
#' @param alpha Significance threshold
#' @param method "FS" (forward stop) or "BS" (basic stop)
#'
#' @return List with type and one
StopRule = function(action, pv, alpha, method){
  if(method == "FS"){
    a.stop = forwardStop(pv, alpha)$stop
  }else if(method == "BS"){
    if(all(pv<alpha)){
      a.stop=length(action)
    }else{
      a.stop = max(0, min(which(pv>alpha))-1)
    }
  }
  sel = action[seq_len(a.stop)]
  if(a.stop==0){
    a.type = "arrhy"
  }else if(length(sel)==1){
    if(sel==1){
      a.type = "rhyI"
    }else{
      a.type = "rhyII"
    }
  }else{
    a.type = "both"
  }
  return(list(type = a.type,
              one = ifelse(action[1]==1, "rhyI", "rhyII")))
}

#' Forward Stopping Rule
#'
#' @param pv P-values
#' @param alpha Threshold
#'
#' @return List with stop and p
forwardStop = function (pv, alpha = 0.1)
{
  if (alpha < 0 || alpha > 1)
    stop("alpha must be in [0,1]")
  if (min(pv, na.rm = TRUE) < 0 || max(pv, na.rm = TRUE) > 1)
    stop("pvalues must be in [0,1]")
  val = -(1/seq_along(pv)) * cumsum(log(1 - pv))
  oo = which(val <= alpha)
  if (length(oo) == 0)
    out = 0
  else out = oo[length(oo)]
  return(list(stop = out,
              p = if (out == 0) NA_real_ else -(1/out) * sum(log(1 - pv[seq_len(out)]))))
}

#' Fisher's Combined P-value
#'
#' @param ps Vector of p-values
#'
#' @return Combined p-value
fishers.p = function(ps){
  fisher.chi = -2*sum(log(ps))
  p.fisher.chi = stats::pchisq(fisher.chi, 2*length(ps), lower.tail = FALSE)
  return(p.fisher.chi)
}

#' Alpha Adjustment for Multiple Testing
#'
#' @param alpha Original alpha
#' @param k Number of tests
#' @param method "bonferroni" or "sidak"
#'
#' @return Adjusted alpha
alpha.adjust = function(alpha, k, method = "bonferroni"){
  if(method == "bonferroni"){
    alpha.new = alpha/k
  }else if(method == "sidak"){
    alpha.new = 1-(1-alpha)^(1/k)
  }
  return(alpha.new)
}

#' Amplitude Cutoff Adjustment
#'
#' @param amp Vector indicating if amplitudes exceed cutoff
#' @param a.TOJR Current TOJR category
#'
#' @return Adjusted TOJR category
amp.cut = function(amp, a.TOJR){
  if(a.TOJR == "arrhy"){
    return("arrhy")
  }else if(a.TOJR == "rhyI"){
    if(amp[1]){
      return("rhyI")
    }else{
      return("arrhy")
    }
  }else if(a.TOJR == "rhyII"){
    if(amp[2]){
      return("rhyII")
    }else{
      return("arrhy")
    }
  }else if(a.TOJR == "both"){
    if(amp[1]&amp[2]){
      return("both")
    }else if(amp[1]){
      return("rhyI")
    }else if(amp[2]){
      return("rhyII")
    }else{
      return("arrhy")
    }
  }
}

# ==============================================================================
# SECTION 4: DIFFERENTIAL R2 TEST (DR)
# ==============================================================================

#' Likelihood Ratio Test for Delta R-squared
#'
#' @param tt1 Time vector for condition 1
#' @param yy1 Expression vector for condition 1
#' @param tt2 Time vector for condition 2
#' @param yy2 Expression vector for condition 2
#' @param period Period (default 24)
#' @param FN Finite sample correction
#'
#' @return P-value for delta R2
#'
#' @examples
#' set.seed(32608)
#' n <- 50
#' tt1 <- runif(n,0,24)
#' Amp1 <- 2
#' Phase1 <- 6
#' Offset1 <- 3
#' yy1 <- Amp1 * sin(2*pi/24 * (tt1 + Phase1)) + Offset1 + rnorm(n,0,1)
#' tt2 <- runif(n,0,24)
#' Amp2 <- 3
#' Phase2 <- 5
#' Offset2 <- 2
#' yy2 <- Amp2 * sin(2*pi/24 * (tt2 + Phase2)) + Offset2 + rnorm(n,0,1)
#' LR_deltaR2(tt1, yy1, tt2, yy2)
#'
#' @export
LR_deltaR2 <- function(tt1, yy1, tt2, yy2, period = 24, FN=TRUE){

  n1 <- length(tt1)
  stopifnot(n1 == length(yy1))
  n2 <- length(tt2)
  stopifnot(length(tt2) == length(yy2))

  w <- 2*pi/period

  fit1 <- fitSinCurve(tt1, yy1, period)
  fit2 <- fitSinCurve(tt2, yy2, period)

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

  asin1 = sin(w * tt1)
  acos1 = cos(w * tt1)
  asin2 = sin(w * tt2)
  acos2 = cos(w * tt2)

  eval_f_list <- function(x,asin1,acos1,asin2,acos2) {
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
    partial_theta1 <-  sum((yy1 - yhat1)^2)/2 - n1/2/theta1

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

  eval_g_eq <- function(x,asin1,acos1,asin2,acos2)
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
    A2_1 * theta1 - A2_2 * theta2
  }

  eval_g_eq_jac <- function(x,asin1,acos1,asin2,acos2)
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
    A2_1 * theta1 - A2_2 * theta2

    c(theta1 * 2 * E1, theta1 * 2 * F1, 0, A2_1,
      - theta2 * 2 * E2, - theta2 * 2 * F2, 0, - A2_2)
  }

  lb <- c(-Inf,-Inf,-Inf,0, -Inf, -Inf,-Inf, 0)
  ub <- c(Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf)

  local_opts <- list( "algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-15 )
  opts <- list( "algorithm"= "NLOPT_LD_SLSQP",
                "xtol_rel"= 1.0e-15,
                "maxeval"= 160000,
                "local_opts" = local_opts,
                "print_level" = 0)

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

  l0 <- - eval_f_list(x_H0,asin1,acos1,asin2,acos2)$objective
  la <- - eval_f_list(x_Ha,asin1,acos1,asin2,acos2)$objective

  LR_stat <- -2*(l0-la)

  dfdiff <- 1
  if(!FN){
    pvalue <- stats::pchisq(LR_stat,dfdiff,lower.tail = FALSE)
  } else if(FN){
    r <- 1
    k <- 6
    n <- n1+n2
    Fstat <- (exp(LR_stat/n) - 1) * (n-k) / r
    pvalue <- stats::pf(Fstat,df1 = r, df2 = n-k, lower.tail = FALSE)
  } else{
    stop("FN has to be TRUE or FALSE")
  }
  pvalue
}

#' Differential Rhythmicity Test (DR)
#'
#' Tests for differential rhythmicity between two groups using R² difference
#'
#' @param x Output from DCP_Rhythmicity(x1, x2)
#' @param method One of "LR" (recommended), "permutation", "bootstrap", "LR_sigma2"
#' @param TOJR TOJR output. If NULL, rhythm.joint object in x will be used
#' @param alpha P-value cutoff
#' @param nSampling Number of samplings for permutation/bootstrap
#' @param Sampling.save Directory to save sampling results
#' @param Sampling.file.label File label for saved results
#' @param p.adjust.method P-value adjustment method
#' @param parallel.ncores Number of cores for parallel processing
#'
#' @return Data frame with differential R² test results
#'
#' @examples
#' x = DCP_sim_data(ngene=1000, nsample=30, A1=c(1, 3), A2=c(1, 3),
#' phase1=c(0, pi/4), phase2=c(pi/4, pi/2),
#' M1=c(4, 6), M2=c(4, 6), sigma1=1, sigma2=1)
#' rhythm.res = DCP_Rhythmicity(x1 = x[[1]], x2 = x[[2]])
#' rhythm.diffR2 = DCP_DiffR2(rhythm.res)
#'
#' @export
DCP_DiffR2 = function(x, method = "LR", TOJR = NULL,  alpha = 0.05,  nSampling=1000,
                     Sampling.save = NULL, Sampling.file.label = "gII_vs_gI",
                     p.adjust.method = "BH", parallel.ncores = 1){
  stopifnot('method should be one of "LR", "permutation", "bootstrap", "LR_sigma2". ' =
                method %in% c("LR", "permutation", "bootstrap", "LR_sigma2"))

  if(is.null(TOJR)){
    overlap.g = x$rhythm.joint$gname[x$rhythm.joint$TOJR!="arrhy"]
  }else{
    stopifnot('The input number of types of joint rhythmicity does not match that of overlapping genes in two groups ' =
                length(x$rhythm.joint$gname)==length(TOJR))
    overlap.g = x$gname_overlap[TOJR != "arrhy"]
  }

  x1 = x$x1
  x2 = x$x2
  x1.overlap = x1$data[match(overlap.g, x1$gname), ]
  x2.overlap = x2$data[match(overlap.g, x2$gname), ]
  t1 = x1$time
  t2 = x2$time
  x1.rhythm = x1$rhythm[match(overlap.g, x1$rhythm$gname),]
  x2.rhythm = x2$rhythm[match(overlap.g, x2$rhythm$gname),]
  stopifnot("x$x1$P is not equal to x$x2$P. " = x$x1$P==x$x2$P)
  period = x$x1$P

  x.list = lapply(seq_along(overlap.g), function(a){
    list(x1.time = t1,
         x2.time = t2,
         y1 = as.numeric(x1.overlap[a, ]),
         y2 = as.numeric(x2.overlap[a, ]))
  })

  if(method == "LR"){
    res.list = parallel::mclapply(seq_along(x.list), function(a){
      one.res = LR_deltaR2(x.list[[a]]$x1.time, x.list[[a]]$y1,
                           x.list[[a]]$x2.time, x.list[[a]]$y2, period,
                           FN = TRUE)
      one.res.tab = data.frame(gname = overlap.g[a],
                               R2.1 = x1.rhythm$R2[a], R2.2 = x2.rhythm$R2[a],
                               delta.R2 =  x2.rhythm$R2[a]-x1.rhythm$R2[a],
                               p.R2 = one.res)
      return(one.res.tab)
    }, mc.cores = parallel.ncores)
    diffR2.tab = do.call(rbind.data.frame, res.list)
  }else if(method == "LR_sigma2"){

    res.list = parallel::mclapply(seq_along(x.list), function(a){
      one.res = LR_diff(x.list[[a]]$x1.time, x.list[[a]]$y1,
                         x.list[[a]]$x2.time, x.list[[a]]$y2, period,
                         FN = TRUE, type = "fit")
      one.res.tab = data.frame(gname = overlap.g[a],
                               R2.1 = one.res[[1]], R2.2 = one.res[[2]],
                               delta.R2 =  one.res[[2]]-one.res[[1]],
                               p.R2 = one.res$pvalue)
      return(one.res.tab)
    }, mc.cores = parallel.ncores)
    diffR2.tab = do.call(rbind.data.frame, res.list)

  }else if(method == "permutation"){
    if(!is.null(Sampling.save)){
      if(!dir.exists(Sampling.save)){
        dir.create(file.path(Sampling.save), recursive = TRUE)
        message(paste0("Directory created. Permutation results will be saved in ", Sampling.save))
      }
    }
    res.tab = diff_rhythmicity_permutation(x1.overlap,x2.overlap,t1,t2,overlap.g, period,
                                           x1.rhythm, x2.rhythm, nSampling,
                                           Sampling.save, Sampling.file.label, parallel.ncores, alpha)
    diffR2.tab = res.tab[, c("gname", "delta.R2", "p.R2")]
  }else if(method == "bootstrap"){
    if(!is.null(Sampling.save)){
      if(!dir.exists(Sampling.save)){
        dir.create(file.path(Sampling.save), recursive = TRUE)
        message(paste0("Directory created. Bootstrap results will be saved in ", Sampling.save))
      }
    }
    res.tab = diff_rhythmicity_bootstrap(x1.overlap,x2.overlap,t1,t2,overlap.g, period,
                                        x1.rhythm, x2.rhythm, nSampling,
                                        Sampling.save, Sampling.file.label, parallel.ncores, alpha)
    diffR2.tab = res.tab[, c("gname", "delta.R2", "p.R2")]
  }

  diffR2.tab$q.R2 = stats::p.adjust(diffR2.tab$p.R2, p.adjust.method)
  return(diffR2.tab)
}

# ==============================================================================
# SECTION 5: DIFFERENTIAL PARAMETER TESTS (DP and DA)
# ==============================================================================

#' Likelihood Ratio Test for Differential Parameters
#'
#' @param time1 Time points for condition 1
#' @param y1 Expression values for condition 1
#' @param time2 Time points for condition 2
#' @param y2 Expression values for condition 2
#' @param period Period
#' @param FN Finite sample correction
#' @param type Parameter type: "amplitude", "phase", "basal"
#'
#' @return List with p-value and parameter estimates
LR_diff = function(time1, y1, time2, y2, period, FN=TRUE, type="amplitude"){

  fit1 = fitSinCurve(time1, y1, period)
  fit2 = fitSinCurve(time2, y2, period)

  if(type == "fit"){
    return(list(R2.1 = fit1$R2, R2.2 = fit2$R2, pvalue = NA))
  }

  A1 = fit1$amp
  A2 = fit2$amp
  phase1 = fit1$phase
  phase2 = fit2$phase
  basal1 = fit1$offset
  basal2 = fit2$offset

  if(type == "amplitude"){
    # Test A1 != A2
    # Use proper likelihood ratio test from LRTest_diff_amp
    result <- LRTest_diff_amp(tt1 = time1, yy1 = y1, tt2 = time2, yy2 = y2,
                               period = period, FN = FN)

    return(list(pvalue = result$pvalue, A1 = result$amp_1, A2 = result$amp_2))

  }else if(type == "phase"){
    # Test phase1 != phase2
    # Use proper likelihood ratio test from LRTest_diff_phase
    result <- LRTest_diff_phase(tt1 = time1, yy1 = y1, tt2 = time2, yy2 = y2,
                                period = period, FN = FN)

    # Compute phase difference for output
    phase_diff <- abs(result$phase_1 - result$phase_2)
    if(phase_diff > period/2) phase_diff <- period - phase_diff

    return(list(pvalue = result$pvalue,
                phase1 = result$phase_1,
                phase2 = result$phase_2,
                delta_phase = phase_diff))

  }else if(type == "basal"){
    # Test basal1 != basal2
    n1 = length(time1)
    n2 = length(time2)
    RSS1 = fit1$rss
    RSS2 = fit2$rss

    sigma2_1 = RSS1/n1
    sigma2_2 = RSS2/n2

    # Two-sample t-test
    se_diff = sqrt(sigma2_1/n1 + sigma2_2/n2)
    delta = basal2 - basal1
    Z = delta / se_diff
    pvalue = 2 * (1 - pnorm(abs(Z)))

    return(list(pvalue = pvalue, M1 = basal1, M2 = basal2,
                 delta_M = delta))
  }
}

#' Differential Parameter Test (DP and DA)
#'
#' Tests for differences in amplitude, phase, and/or mesor between two groups
#'
#' @param x Output from DCP_Rhythmicity(x1, x2)
#' @param Par One of "A", "phase", "M", "A&phase", or "A&phase&M"
#' @param TOJR TOJR output. If NULL, rhythm.joint object in x will be used
#' @param alpha P-value cutoff
#' @param p.adjust.method P-value adjustment method
#' @param parallel.ncores Number of cores for parallel processing
#'
#' @return Data frame with differential parameter test results
#'
#' @examples
#' x = DCP_sim_data(ngene=1000, nsample=30, A1=c(1, 3), A2=c(1, 3),
#' phase1=c(0, pi/4), phase2=c(pi/4, pi/2),
#' M1=c(4, 6), M2=c(4, 6), sigma1=1, sigma2=1)
#' rhythm.res = DCP_Rhythmicity(x1 = x[[1]], x2 = x[[2]])
#' rhythm.diffPar = DCP_DiffPar(rhythm.res, Par = "A&phase")
#'
#' @export
DCP_DiffPar = function(x, Par = c("A"), TOJR=NULL, alpha = 0.05,
                       p.adjust.method = "BH", parallel.ncores = 1){
  stopifnot('Par should be one of "A", "phase", "M", "A&phase", or "A&phase&M". '
            = Par %in% c("A", "phase", "M", "A&phase", "A&phase&M"))
  stopifnot(
    "x should be output of DCP_Rhythmicity() with both x1 and x2 not NULL" =
      all(c("x1", "x2", "gname_overlap", "rhythm.joint")%in%names(x)))

  if(is.null(TOJR)){
    overlap.g = x$rhythm.joint$gname[x$rhythm.joint$TOJR == "both"]
    stopifnot('There is no RhyBoth genes, please use a lower threshold for TOJR. ' =
                length(overlap.g)>0)
  }else{
    stopifnot('The input number of types of joint rhythmicity does not match that of overlapping genes in two groups ' =
                length(x$rhythm.joint$gname)==length(TOJR))
    overlap.g = x$gname_overlap[TOJR == "both"]
  }

  x1 = x$x1
  x2 = x$x2
  x1.overlap = x1$data[match(overlap.g, x1$gname), ]
  x2.overlap = x2$data[match(overlap.g, x2$gname), ]
  t1 = x1$time
  t2 = x2$time
  x1.rhythm = x1$rhythm[match(overlap.g, x1$rhythm$gname), ]
  x2.rhythm = x2$rhythm[match(overlap.g, x2$rhythm$gname), ]
  stopifnot("x$x1$P is not equal to x$x2$P. " = x$x1$P==x$x2$P)
  period = x$x1$P

  x.list = lapply(seq_along(overlap.g), function(a){
    list(x1.time = t1,
         x2.time = t2,
         y1 = as.numeric(x1.overlap[a, ]),
         y2 = as.numeric(x2.overlap[a, ]))
  })

  if(Par == "A"|Par == "phase"|Par == "M"){
    Par2 = switch(Par, "A" = "amplitude", "phase" = "phase", "M" = "basal")

    test_diffPar = parallel::mclapply(seq_along(x.list), function(a){
      out.diffPar = LR_diff(x.list[[a]]$x1.time, x.list[[a]]$y1,
                           x.list[[a]]$x2.time, x.list[[a]]$y2,
                           period , FN = TRUE, type=Par2)
      one.row = data.frame(gname = overlap.g[a],
                           Par1 = ifelse(Par2=="phase", x1.rhythm[, "peak"][a],
                                         x1.rhythm[, Par][a]),
                           Par2 = ifelse(Par2=="phase", x2.rhythm[, "peak"][a],
                                         x2.rhythm[, Par][a]),
                           delta.Par = ifelse(Par2=="phase",
                                              x2.rhythm[, "peak"][a]-x1.rhythm[, "peak"][a],
                                              x2.rhythm[, Par][a]-x1.rhythm[, Par][a]),
                           pvalue = out.diffPar$pvalue
      )
      return(one.row)
    }, mc.cores = parallel.ncores)
    diffPar.tab = do.call(rbind.data.frame, test_diffPar)
    diffPar.tab$qvalue= stats::p.adjust(diffPar.tab$pvalue, p.adjust.method)
    colnames(diffPar.tab)[2:4] = gsub("Par", Par, colnames(diffPar.tab)[2:4])

  }else{
    data =  cbind.data.frame(x1.overlap, x2.overlap)
    time = c(t1, t2)
    group = c(rep(1, ncol(x1.overlap)), rep(2, ncol(x2.overlap)))
    group.1 = group == 1
    group.2 = group == 2

    design.full = data.frame(
      M1 = 1,
      g2 = group.2,
      inphase = cos(2*pi/period*time), outphase = sin(2*pi/period*time),
      g2.inphase = group.2*cos(2*pi/period*time),
      g2.outphase = group.2*sin(2*pi/period*time))
    fit.full = limma::lmFit(data,
                            stats::model.matrix(~g2+inphase+outphase+g2.inphase+g2.outphase,
                                                data = design.full))
    fit.full = limma::eBayes(fit.full)

    if(Par == "A&phase"){
      diff.top = limma::topTable(fit.full, coef = 5:6, n = nrow(data), sort.by = "none")
      test_diffPar = parallel::mclapply(seq_along(x.list), function(a){
        out.diffA = LR_diff(x.list[[a]]$x1.time, x.list[[a]]$y1,
                           x.list[[a]]$x2.time, x.list[[a]]$y2,
                           period , FN = TRUE, type="amplitude")
        out.diffphase= LR_diff(x.list[[a]]$x1.time, x.list[[a]]$y1,
                              x.list[[a]]$x2.time, x.list[[a]]$y2,
                              period , FN = TRUE, type="phase")
        one.row = data.frame(
          A1 = x1.rhythm$A[a],
          A2 = x2.rhythm$A[a],
          delta.A = x2.rhythm$A[a]-x1.rhythm$A[a],
          p.delta.A = out.diffA$pvalue,
          phase1 = x1.rhythm$peak[a],
          phase2 = x2.rhythm$peak[a],
          delta.phase = x2.rhythm$peak[a]-x1.rhythm$peak[a],
          p.delta.phase = out.diffphase$pvalue
        )
        return(one.row)
      }, mc.cores = parallel.ncores)
      diffPar.tab = do.call(rbind.data.frame, test_diffPar)
      all.tab = data.frame(gname = overlap.g,
                           diff.top[, c("P.Value", "adj.P.Val")], diffPar.tab)
      colnames(all.tab)[2:3] = c("p.overall", "q.overall")
      all.tab$q.delta.A = stats::p.adjust(all.tab$p.delta.A, p.adjust.method)
      all.tab$q.delta.phase = stats::p.adjust(all.tab$p.delta.phase, p.adjust.method)
      all.tab$post.hoc.A.By.q = all.tab$q.overall<alpha&all.tab$q.delta.A<alpha.adjust(alpha, 2, method = "sidak")
      all.tab$post.hoc.phase.By.q = all.tab$q.overall<alpha&all.tab$q.delta.phase<alpha.adjust(alpha, 2, method = "sidak")
      all.tab$post.hoc.A.By.p = all.tab$p.overall<alpha&all.tab$p.delta.A<alpha.adjust(alpha, 2, method = "sidak")
      all.tab$post.hoc.phase.By.p = all.tab$p.overall<alpha&all.tab$p.delta.phase<alpha.adjust(alpha, 2, method = "sidak")
      diffPar.tab = all.tab[, c("gname", "p.overall", "q.overall","post.hoc.A.By.p","post.hoc.phase.By.p", "post.hoc.A.By.q", "post.hoc.phase.By.q",
                                "A1", "A2", "delta.A", "p.delta.A", "q.delta.A", "phase1", "phase2", "delta.phase", "p.delta.phase", "q.delta.phase")]
    }else if(Par == "A&phase&M"){
      diff.top = limma::topTable(fit.full, coef = c(2, 5:6), n = nrow(data), sort.by = "none")
      test_diffPar = parallel::mclapply(seq_along(x.list), function(a){
        out.diffA = LR_diff(x.list[[a]]$x1.time, x.list[[a]]$y1,
                           x.list[[a]]$x2.time, x.list[[a]]$y2,
                           period , FN = TRUE, type="amplitude")
        out.diffphase= LR_diff(x.list[[a]]$x1.time, x.list[[a]]$y1,
                              x.list[[a]]$x2.time, x.list[[a]]$y2,
                              period , FN = TRUE, type="phase")
        out.diffM= LR_diff(x.list[[a]]$x1.time, x.list[[a]]$y1,
                          x.list[[a]]$x2.time, x.list[[a]]$y2,
                          period , FN = TRUE, type="basal")
        one.row = data.frame(
          A1 = x1.rhythm$A[a],
          A2 = x2.rhythm$A[a],
          delta.A = x2.rhythm$A[a]-x1.rhythm$A[a],
          p.delta.A = out.diffA$pvalue,
          phase1 = x1.rhythm$peak[a],
          phase2 = x2.rhythm$peak[a],
          delta.phase = x2.rhythm$peak[a]-x1.rhythm$peak[a],
          p.delta.phase = out.diffphase$pvalue,
          M1 = x1.rhythm$M[a],
          M2 = x2.rhythm$M[a],
          delta.M = x2.rhythm$M[a]-x1.rhythm$M[a],
          p.delta.M = out.diffM$pvalue
        )
        return(one.row)
      }, mc.cores = parallel.ncores)
      diffPar.tab = do.call(rbind.data.frame, test_diffPar)
      all.tab = data.frame(gname = overlap.g, diff.top[, c("P.Value", "adj.P.Val")], diffPar.tab)
      colnames(all.tab)[2:3] = c("p.overall", "q.overall")
      all.tab$q.delta.A = stats::p.adjust(all.tab$p.delta.A, p.adjust.method)
      all.tab$q.delta.phase = stats::p.adjust(all.tab$p.delta.phase, p.adjust.method)
      all.tab$q.delta.M = stats::p.adjust(all.tab$p.delta.M, p.adjust.method)
      all.tab$post.hoc.A.By.q = all.tab$q.overall<alpha&all.tab$q.delta.A<alpha.adjust(alpha, 3, method = "sidak")
      all.tab$post.hoc.phase.By.q = all.tab$q.overall<alpha&all.tab$q.delta.phase<alpha.adjust(alpha, 3, method = "sidak")
      all.tab$post.hoc.M.By.q = all.tab$q.overall<alpha&all.tab$q.delta.M<alpha.adjust(alpha, 3, method = "sidak")
      all.tab$post.hoc.A.By.p = all.tab$p.overall<alpha&all.tab$p.delta.A<alpha.adjust(alpha, 3, method = "sidak")
      all.tab$post.hoc.phase.By.p = all.tab$p.overall<alpha&all.tab$p.delta.phase<alpha.adjust(alpha, 3, method = "sidak")
      all.tab$post.hoc.M.By.p = all.tab$p.overall<alpha&all.tab$p.delta.M<alpha.adjust(alpha, 3, method = "sidak")
      diffPar.tab = all.tab[, c("gname", "p.overall", "q.overall","post.hoc.A.By.p","post.hoc.phase.By.p", "post.hoc.M.By.p",
                                "post.hoc.A.By.q", "post.hoc.phase.By.q", "post.hoc.M.By.q",
                                "A1", "A2", "delta.A", "p.delta.A", "q.delta.A", "phase1", "phase2", "delta.phase", "p.delta.phase", "q.delta.phase",
                                "M1", "M2", "delta.M", "p.delta.M", "q.delta.M")]
    }

  }

  colnames(diffPar.tab) = gsub("phase", "peak", colnames(diffPar.tab))
  diffPar.tab$P = period
  return(diffPar.tab)
}

# ==============================================================================
# SECTION 6: MAIN WRAPPER FUNCTIONS
# ==============================================================================

#' Complete Differential Circadian Analysis
#'
#' Performs differential tests (DR, DP) with classification
#'
#' @param expr1 Expression matrix for group 1 (genes × samples)
#' @param expr2 Expression matrix for group 2 (genes × samples)
#' @param times1 Time points for group 1
#' @param times2 Time points for group 2
#' @param alpha Significance threshold (default 0.05)
#' @param gene_names Gene names (optional)
#' @param parallel.ncores Number of cores for parallel processing
#'
#' @return List with:
#'   \item{DR}{Differential rhythmicity results}
#'   \item{DP}{Differential phase results}
#'   \item{classification}{Gene classifications (DR/DP/NS)}
#'
#' @examples
#' x = DCP_sim_data(ngene=1000, nsample=30, A1=c(1, 3), A2=c(1, 3),
#' phase1=c(0, pi/4), phase2=c(pi/4, pi/2),
#' M1=c(4, 6), M2=c(4, 6), sigma1=1, sigma2=1)
#' dcp_result = DCP_Analyze(x[[1]]$data, x[[2]]$data,
#'                        x[[1]]$time, x[[2]]$time)
#'
#' @export
DCP_Analyze <- function(expr1, expr2, times1, times2, alpha = 0.05,
                        gene_names = NULL, parallel.ncores = 1){

  ngenes <- nrow(expr1)

  if (is.null(gene_names)) {
    gene_names <- paste0("Gene", 1:ngenes)
  }

  # Prepare data in the format expected by DCP_Rhythmicity
  x1 <- list(data = as.data.frame(expr1),
             time = times1,
             gname = gene_names)

  x2 <- list(data = as.data.frame(expr2),
             time = times2,
             gname = gene_names)

  # Run rhythmicity analysis
  rhythm.res <- DCP_Rhythmicity(x1, x2,
                                method = "Sidak_FS",
                                period = 24,
                                amp.cutoff = 0,
                                alpha = alpha,
                                CI = FALSE,
                                parallel.ncores = parallel.ncores)

  # Test 1: Differential Rhythmicity (DR)
  dr_results <- DCP_DiffR2(rhythm.res, method = "LR", alpha = alpha,
                            parallel.ncores = parallel.ncores)

  # Test 2: Differential Parameters (DP and DA)
  dp_results <- DCP_DiffPar(rhythm.res, Par = "A&phase",
                            alpha = alpha,
                            parallel.ncores = parallel.ncores)

  # Classification
  # Priority: DR > DP
  classification <- rep("NS", ngenes)

  # Match gene names
  match_idx <- match(gene_names, dr_results$gname)

  # DR: differential rhythmicity (use BH-adjusted q.R2, not raw p.R2)
  dr_idx <- which(dr_results$q.R2 < alpha)
  if(length(dr_idx) > 0){
    classification[dr_idx] <- "DR"
  }

  # For genes not DR, test DP
  # Match by gene name: dp_results rows correspond to "both"-classified genes
  # (a subset of ngenes), so positional indexing into dp_results by 1:ngenes
  # would silently pick wrong rows. Always index by gname.
  non_dr <- setdiff(seq_len(ngenes), dr_idx)
  if (length(non_dr) > 0L && nrow(dp_results) > 0L &&
      "gname" %in% colnames(dp_results) &&
      "post.hoc.phase.By.q" %in% colnames(dp_results)) {
    dp_genes <- dp_results$gname[dp_results$post.hoc.phase.By.q]
    dp_idx   <- which(gene_names %in% dp_genes & seq_len(ngenes) %in% non_dr)
    if (length(dp_idx) > 0L) classification[dp_idx] <- "DP"
  }

  # Compile results
  results <- list(
    rhythm = rhythm.res,
    DR = dr_results,
    DP = dp_results,
    classification = data.frame(
      gname = gene_names,
      classification = classification,
      stringsAsFactors = FALSE
    )
  )

  return(results)
}

# ==============================================================================
# SECTION 7: HELPER FUNCTIONS (CIs, etc.)
# ==============================================================================

#' Get Phase from Beta Coefficients
#'
get_phase = function(b1.x, b2.x){
  ph.x = atan(-b2.x/b1.x)
  if(b2.x>0){
    if(ph.x<0){
      ph.x = ph.x+2*pi
    }else if (ph.x>0){
      ph.x = ph.x+pi
    }
  }else if(b2.x<0){
    if(ph.x<0){
      ph.x = ph.x+pi
    }
  }else{
    # b2.x == 0: sine coefficient is zero; peak at 0 or pi depending on cosine sign
    ph.x = if (b1.x >= 0) 0 else pi
  }
  return(list(phase = ph.x, tan = -b2.x/b1.x))
}

#' Calculate CI for Mesor
#'
calculate_CI.M = function(XX.inv, A.t, r.full = 6, n, alpha = 0.05, sigma.hat){
  CI.m.hat.radius = stats::qt(1-alpha/2, n-r.full)*sigma.hat*sqrt(A.t%*%XX.inv%*%t(A.t))
  return(CI.m.hat.radius)
}

#' Calculate CI for Amplitude and Phase (Taylor method)
#'
calculate_CI_A.phase.Taylor = function(XX.inv, A.t, phase.hat, A.hat, sigma.hat){
  var.new = A.t%*%XX.inv%*%t(A.t)
  var.beta1 = var.new[1, 1]
  var.beta2 = var.new[2, 2]
  var.beta1.beta2 = var.new[1, 2]
  se.A.hat = sigma.hat*sqrt(var.beta1*cos(phase.hat)^2
                            -2*var.beta1.beta2*sin(phase.hat)*cos(phase.hat)
                            +var.beta2*sin(phase.hat)^2)
  se.phase.hat = sigma.hat*sqrt(var.beta1*sin(phase.hat)^2
                                +2*var.beta1.beta2*sin(phase.hat)*cos(phase.hat)
                                +var.beta2*cos(phase.hat)^2)/A.hat
  return(list(se.A.hat = se.A.hat,
              se.phase.hat = se.phase.hat))
}

#' Choose the minimum absolute value from each row
#'
choose.abs.min = function(mat = matrix(stats::rnorm(21), ncol = 3)){
  abs.min.vec = apply(mat, 1, function(a){
    min.abs.idx = which.min(abs(a))
    return(a[min.abs.idx])
  })
  return(abs.min.vec)
}

# ==============================================================================
# SECTION 8: UTILITY FUNCTIONS
# ==============================================================================

#' Simulate Data for Rhythmicity Comparison
#'
#' Simple two-group study simulation for illustration
#'
#' @param ngene Number of genes
#' @param nsample Number of samples per group
#' @param A1 Amplitude range for group 1 (c(min, max))
#' @param A2 Amplitude range for group 2
#' @param phase1 Phase range for group 1
#' @param phase2 Phase range for group 2
#' @param M1 MESOR range for group 1
#' @param M2 MESOR range for group 2
#' @param sigma1 Noise level for group 1
#' @param sigma2 Noise level for group 2
#'
#' @return List of two elements, each is a list with data, time, gname components
#'
#' @examples
#' x = DCP_sim_data(ngene=1000, nsample=30, A1=c(1, 3), A2=c(1, 3),
#' phase1=c(0, pi/4), phase2=c(pi/4, pi/2),
#' M1=c(4, 6), M2=c(4, 6), sigma1=1, sigma2=1)
#'
#' @export
DCP_sim_data = function(ngene=1000, nsample=30, A1=c(1, 3), A2=c(1, 3),
                       phase1=c(0, 2*pi), phase2=c(0, 2*pi),
                       M1=c(4, 6), M2=c(4, 6), sigma1=1, sigma2=1){
  a.A1 = stats::runif(ngene, A1[1], A1[2])
  a.A2 = stats::runif(ngene, A2[1], A2[2])
  a.m1 = stats::runif(ngene, M1[1], M1[2])
  a.m2 = stats::runif(ngene, M2[1], M2[2])
  a.phase1 = stats::runif(ngene, phase1[1], phase1[2])
  a.phase2 = stats::runif(ngene, phase2[1], phase2[2])
  a.sigma1 = sigma1
  a.sigma2 = sigma2
  x1.time = stats::runif(nsample, min = 0, max = 24)
  x2.time = stats::runif(nsample, min = 0, max = 24)

  noise.mat1 = matrix(stats::rnorm(ngene*nsample, 0, a.sigma1), ncol = nsample, nrow = ngene)
  signal.mat1 = t(sapply(seq_len(ngene), function(a){
    a.m1[a]+a.A1[a]*cos(2*pi/24*x1.time+a.phase1[a])
  }))

  noise.mat2 = matrix(stats::rnorm(ngene*nsample, 0, a.sigma2), ncol = nsample, nrow = ngene)
  signal.mat2 = t(sapply(seq_len(ngene), function(a){
    a.m2[a]+a.A2[a]*cos(2*pi/24*x2.time+a.phase2[a])
  }))

  x1.data = as.data.frame(noise.mat1 + signal.mat1)
  x2.data = as.data.frame(noise.mat2 + signal.mat2)

  x1 = list(data = x1.data,
            time = x1.time,
            gname = paste("gene", seq_len(ngene)))

  x2 = list(data = x2.data,
            time = x2.time,
            gname = paste("gene", seq_len(ngene)))
  return(list(x1, x2))
}

#' Summary of Differential Rhythmicity Tests
#'
#' @param result List of p-value vectors
#' @param test Test names
#' @param test Types of results
#' @param val Cutoffs
#' @param out "long" or "wide" format
#'
#' @return Summary table
#'
#' @export
SummarizeDR = function(result, test = "DRF", type = "p-value", val = 0.05, out = "long"){
  if(!is.list(result)){
    result = list(result)
  }
  l = length(result)
  stopifnot("Arguments `result`, `test`, `type`, `val` must have the same length" =
              (length(test)==l&length(type)==l&length(val)==l))
  tab = lapply(seq_len(l), function(a){
    a.result = result[[a]]
    a.test = test[a]
    a.type = type[a]
    a.val = val[a]
    a.tab0 = table(a.result<a.val)
    a.tab = data.frame(nTRUE = unname(ifelse(is.na(a.tab0["TRUE"]), 0, a.tab0["TRUE"])),
                       nFALSE = unname(ifelse(is.na(a.tab0["FALSE"]), 0, a.tab0["FALSE"])),
                       test = a.test,
                       cutoff = paste0(a.type, "<", a.val))
    return(a.tab)
  })
  tab = do.call(rbind.data.frame, tab)
  tab$nTRUE = paste0(tab$nTRUE, "(/", tab$nTRUE+tab$nFALSE, ")")
  if(out=="long"){
    return(tab[, c(1, 3,4)])
  }else if(out =="wide"){
    tab.list = lapply(seq_len(nrow(tab[, c(1, 3,4)])), function(i){
      a.cell = data.frame(x = c(tab[i, "nTRUE"]))
      rownames(a.cell) = NULL
      colnames(a.cell) = paste0(tab[i, "test"], " ", tab[i, "cutoff"])
      return(a.cell)
    })
    tab2 = do.call(cbind.data.frame, tab.list)
    return(tab2)
  }
}


# ==============================================================================
# SECTION 9: BRIDGE FUNCTIONS FOR SIMULATION FRAMEWORK
# ==============================================================================

#' Bridge: Convert simulation output to DCP format
#'
#' @param expr_matrix Gene x sample matrix
#' @param times Time vector
#' @param gene_names Gene name vector
#'
#' @return DCP-formatted list
#' @export
format_for_DCP <- function(expr_matrix, times, gene_names = NULL) {
  if (is.null(gene_names)) {
    gene_names = paste0("Gene", 1:nrow(expr_matrix))
  }
  
  list(
    data = as.data.frame(expr_matrix),
    time = times,
    gname = gene_names
  )
}

#' Run full DCP differential analysis pipeline
#'
#' @description Uses the complete DiffCircadian pipeline with TOJR,
#' hierarchical testing, and proper LR tests. This is the RIGOROUS approach.
#'
#' @param expr1 Expression matrix for group 1 (genes × samples)
#' @param expr2 Expression matrix for group 2 (genes × samples)
#' @param times1 Time points for group 1
#' @param times2 Time points for group 2
#' @param gene_names Gene names (optional)
#' @param period Period (default 24)
#' @param alpha Significance threshold
#' @param test_DR Test differential rhythmicity?
#' @param test_DP Test differential phase?
#'
#' @return List with p-values for requested tests
#' @export
run_DCP_pipeline <- function(expr1, expr2, times1, times2,
                             gene_names = NULL,
                             period = 24,
                             alpha = 0.05,
                             test_DR = TRUE,
                             test_DP = TRUE) {
  
  # Format data
  x1 = format_for_DCP(expr1, times1, gene_names)
  x2 = format_for_DCP(expr2, times2, gene_names)
  
  # Step 1: TOJR classification
  rhythm_res = DCP_Rhythmicity(x1, x2, 
                               method = "Sidak_FS",
                               period = period,
                               alpha = alpha,
                               CI = FALSE,
                               parallel.ncores = 1)
  
  # Initialize results
  ngenes = nrow(expr1)
  results = list(
    p_DR = rep(1, ngenes),
    p_DP = rep(1, ngenes),
    TOJR = rhythm_res$rhythm.joint$TOJR
  )
  
  # Step 2: Test DR (only for genes not in "arrhy")
  if (test_DR) {
    dr_results = tryCatch({
      DCP_DiffR2(rhythm_res, method = "LR", alpha = alpha)
    }, error = function(e) {
      NULL
    })
    
    if (!is.null(dr_results)) {
      match_idx = match(gene_names, dr_results$gname)
      results$p_DR[!is.na(match_idx)] = dr_results$p.R2[match_idx[!is.na(match_idx)]]
    }
  }
  
  # Step 3: Test DP (only for "both" genes)
  if (test_DP) {
    # Check if there are "both" genes
    n_both = sum(rhythm_res$rhythm.joint$TOJR == "both")

    if (n_both > 0) {
      par_results = tryCatch({
        DCP_DiffPar(rhythm_res, Par = "A&phase", alpha = alpha)
      }, error = function(e) {
        NULL
      })

      if (!is.null(par_results)) {
        match_idx = match(gene_names, par_results$gname)
        results$p_DP[!is.na(match_idx)] = par_results$p.delta.peak[match_idx[!is.na(match_idx)]]
      }
    }
  }
  
  return(results)
}

#' Extract p-values from DCP pipeline for power analysis
#'
#' @description Simplified wrapper that returns just p-values
#' for integration with runSimsDiff()
#'
#' @inheritParams run_DCP_pipeline
#'
#' @return List with p_DR, p_DP
#' @export
testAnyDifferential <- function(y1, y2, times1, times2,
                                period = 24,
                                test_types = c("DR", "DP")) {

  # Create single-gene matrices
  expr1 = matrix(y1, nrow = 1)
  expr2 = matrix(y2, nrow = 1)

  # Run pipeline
  results = run_DCP_pipeline(
    expr1, expr2, times1, times2,
    gene_names = "Gene1",
    period = period,
    test_DR = "DR" %in% test_types,
    test_DP = "DP" %in% test_types
  )

  return(list(
    p_DR = results$p_DR[1],
    p_DP = results$p_DP[1]
  ))
}


# ==============================================================================
# SECTION 10: CIRCACOMPARE WRAPPER
# ==============================================================================

#' CircaCompare Differential Analysis Wrapper
#'
#' Per-gene NLS fitting with group-interaction parameters.
#' Model: y ~ (k + k1*group) + (alpha + alpha1*group) * cos(tau*t - (phi + phi1*group))
#' Provides separate Wald t-tests for amplitude diff (alpha1 -> DA),
#' phase diff (phi1 -> DP), and mesor diff (k1 -> DM).
#'
#' @param expr1 Expression matrix for group 1 (genes x samples)
#' @param times1 Time points for group 1
#' @param expr2 Expression matrix for group 2 (genes x samples)
#' @param times2 Time points for group 2
#' @param gene_names Gene name vector (optional)
#' @param period Circadian period (default 24)
#' @param alpha_threshold Significance level for per-group rhythmicity (default 0.05)
#'
#' @return List with pval_DR (NA — no equivalent), pval_DP, pval_DM
#'         Each is a numeric vector of length ngenes; default 1 for failed/untested.
#'
#' @details CircaCompare calls stop() when either group is arrhythmic, so every
#'          gene call is wrapped in tryCatch. NLS failures also yield p=1 (conservative).
#'          CircaCompare has no differential rhythmicity (DR) test — only per-group
#'          rhythmicity p-values. pval_DR is set to NA for all genes.
#' @export
detect_CircaCompare <- function(expr1, times1, expr2, times2,
                                gene_names = NULL, period = 24,
                                alpha_threshold = 0.05) {

  if (!requireNamespace("circacompare", quietly = TRUE)) {
    stop("Package 'circacompare' is required. Install with: install.packages('circacompare')")
  }

  ngenes <- nrow(expr1)
  if (is.null(gene_names)) gene_names <- paste0("Gene", 1:ngenes)

  n1 <- ncol(expr1)
  n2 <- ncol(expr2)

  pval_DP <- rep(1, ngenes)
  pval_DM <- rep(1, ngenes)

  for (g in seq_len(ngenes)) {
    tryCatch({
      df <- data.frame(
        time    = c(times1, times2),
        measure = c(as.numeric(expr1[g, ]), as.numeric(expr2[g, ])),
        group   = factor(c(rep("g1", n1), rep("g2", n2)))
      )

      res <- circacompare::circacompare(
        x            = df,
        col_time     = "time",
        col_group    = "group",
        col_outcome  = "measure",
        period       = period,
        alpha_threshold = alpha_threshold,
        timeout_n    = 500,   # default 10000 is overkill; most genes converge fast
        suppress_all = TRUE
      )

      summ <- res$summary
      # Match by parameter name for robustness (row indices may shift)
      p_mesor_row <- grep("P-value for mesor", summ$parameter, ignore.case = TRUE)
      p_phase_row <- grep("P-value for.*phase", summ$parameter, ignore.case = TRUE)

      if (length(p_mesor_row) > 0) pval_DM[g] <- as.numeric(summ$value[p_mesor_row[1]])
      if (length(p_phase_row) > 0) pval_DP[g] <- as.numeric(summ$value[p_phase_row[1]])

    }, error = function(e) {
      # NLS failure or arrhythmic group -> p=1 (conservative)
    })
  }

  list(
    pval_DR = rep(NA, ngenes),   # CircaCompare has no DR test
    pval_DP = pval_DP,
    pval_DM = pval_DM
  )
}


# ==============================================================================
# SECTION 11: JTK_CYCLE RHYTHMICITY WRAPPER
# ==============================================================================

#' JTK_CYCLE Rhythmicity Detection
#'
#' Non-parametric rhythmicity test based on Kendall's tau-b.
#' Uses MetaCycle::meta2d with cycMethod = "JTK".
#'
#' @param expr Expression matrix (genes x samples)
#' @param times Time points (numeric vector, hours)
#' @param gene_names Gene name vector (optional)
#' @param period Period to test (default 24)
#'
#' @return Numeric vector of adjusted p-values per gene (length ngenes).
#'         Returns 1 for genes where the test fails.
#'
#' @details JTK_CYCLE is designed for evenly-spaced time series. For passive
#'          designs with uneven sampling, MetaCycle handles the conversion
#'          internally, but power may be reduced compared to evenly-spaced data.
#' @export
detect_JTK <- function(expr, times, gene_names = NULL, period = 24) {

  if (!requireNamespace("MetaCycle", quietly = TRUE)) {
    stop("Package 'MetaCycle' is required. Install with: install.packages('MetaCycle')")
  }

  ngenes <- nrow(expr)
  if (is.null(gene_names)) gene_names <- paste0("Gene", 1:ngenes)

  # MetaCycle::meta2d expects data.frame: first column = gene IDs, then expression
  inDF <- data.frame(CycID = gene_names, as.data.frame(expr),
                     check.names = FALSE, stringsAsFactors = FALSE)

  result <- tryCatch({
    MetaCycle::meta2d(
      infile      = "dummy",
      filestyle   = "txt",
      inDF        = inDF,
      timepoints  = times,
      cycMethod   = "JTK",
      outputFile  = FALSE,
      minper      = max(period - 4, 18),
      maxper      = period + 4,
      releaseNote = FALSE
    )
  }, error = function(e) {
    warning(sprintf("MetaCycle/JTK failed: %s", e$message))
    NULL
  })

  if (!is.null(result) && !is.null(result$JTK)) {
    idx <- match(gene_names, result$JTK$CycID)
    pvals <- result$JTK$ADJ.P[idx]
    pvals[is.na(pvals)] <- 1
    return(pvals)
  }

  rep(1, ngenes)
}


# ==============================================================================
# SECTION 12: RAIN RHYTHMICITY WRAPPER
# ==============================================================================

#' RAIN Rhythmicity Detection
#'
#' Non-parametric rhythmicity test allowing asymmetric waveforms.
#' Uses the rain::rain function from Bioconductor.
#'
#' @param expr Expression matrix (genes x samples)
#' @param times Time points (numeric vector, hours)
#' @param gene_names Gene name vector (optional)
#' @param period Period to test (default 24)
#'
#' @return Numeric vector of p-values per gene (length ngenes).
#'         Returns 1 for genes where the test fails.
#'
#' @details RAIN is designed for evenly-spaced time series. For uneven sampling,
#'          data is sorted by time and the median interval is used as deltat.
#'          Power may be reduced compared to evenly-spaced designs.
#' @export
detect_RAIN <- function(expr, times, gene_names = NULL, period = 24) {

  if (!requireNamespace("rain", quietly = TRUE)) {
    stop("Package 'rain' is required. Install with: BiocManager::install('rain')")
  }

  ngenes <- nrow(expr)
  if (is.null(gene_names)) gene_names <- paste0("Gene", 1:ngenes)

  # Sort by time (RAIN expects time-ordered data)
  ord <- order(times)
  expr_sorted <- expr[, ord, drop = FALSE]
  times_sorted <- times[ord]

  # deltat from unique time points; nr.series from replicate count
  unique_times <- sort(unique(times_sorted))
  B_pts  <- length(unique_times)
  deltat <- if (B_pts > 1) median(diff(unique_times)) else period
  counts <- tabulate(match(times_sorted, unique_times))
  nr     <- if (length(unique(counts)) == 1L && counts[1L] > 1L) counts[1L] else 1L

  result <- tryCatch({
    # rain expects: rows = samples (time-ordered), columns = genes
    rain::rain(
      x         = t(expr_sorted),
      deltat    = deltat,
      period    = period,
      nr.series = nr
    )
  }, error = function(e) {
    warning(sprintf("RAIN failed: %s", e$message))
    NULL
  })

  if (!is.null(result)) {
    pvals <- result[, "pVal"]
    pvals[is.na(pvals)] <- 1
    return(pvals)
  }

  rep(1, ngenes)
}


# ==============================================================================
# SECTION 13: UNIFIED SINGLE-COHORT DETECTION WRAPPERS
# All return numeric[G] of raw p-values (caller handles FDR adjustment).
# ==============================================================================

#' DCP single-cohort rhythmicity detection
#' @param expr   Gene x sample expression matrix
#' @param times  Numeric time vector (length = ncol(expr))
#' @param period Period (default 24)
#' @return Numeric vector of p-values, length nrow(expr)
#' @export
detect_DCP <- function(expr, times, period = 24) {
  pvals <- tryCatch(
    fitCosinorAll(expr, times = times, period = period)$pvalue,
    error = function(e) rep(NA_real_, nrow(expr))
  )
  pvals[is.na(pvals)] <- 1
  pvals
}


#' Multi-harmonic rhythmicity detection (adaptive K = floor((B-1)/2))
#'
#' Fits K harmonics where B = number of unique time points.
#' Tests all harmonic terms jointly via an F-test against intercept-only null.
#'
#' @param expr   Gene x sample expression matrix
#' @param times  Numeric time vector (length = ncol(expr))
#' @param period Period (default 24)
#' @return Numeric vector of p-values, length nrow(expr)
#' @export
detect_MH <- function(expr, times, period = 24) {
  B     <- length(unique(times))
  K     <- max(1L, floor((B - 1L) / 2L))
  omega <- 2 * pi / period
  N     <- length(times)

  dm_full <- matrix(1, nrow = N, ncol = 1L + 2L * K)
  for (k in seq_len(K)) {
    dm_full[, 2L * k]      <- cos(k * omega * times)
    dm_full[, 2L * k + 1L] <- sin(k * omega * times)
  }
  dm_null <- matrix(1, nrow = N, ncol = 1L)
  df1 <- 2L * K
  df2 <- N - 1L - 2L * K
  if (df2 < 1L) return(rep(1, nrow(expr)))

  apply(expr, 1L, function(y) {
    tryCatch({
      rss_full <- sum(lm.fit(dm_full, y)$residuals^2)
      rss_null <- sum(lm.fit(dm_null, y)$residuals^2)
      F_stat   <- ((rss_null - rss_full) / df1) / (rss_full / df2)
      pf(F_stat, df1, df2, lower.tail = FALSE)
    }, error = function(e) 1)
  })
}


# ==============================================================================
# SECTION 14: UNIFIED TWO-GROUP DIFFERENTIAL DETECTION WRAPPERS
# All return list(pval_DR, pval_DP, pval_DM).
# NA for test types the method does not support.
#
# Method × test_type support:
#   DCP          DR  DP  --
#   CircaCompare --  DP  DM
#   LimoRhyde    DR  --  --
#   DODR         DR  --  --
# ==============================================================================

#' DCP two-group differential detection
#'
#' Wraps run_DCP_pipeline() (TOJR → DiffR2 → DiffPar).
#' DM is not tested by DCP; pval_DM is returned as NA.
#'
#' @param expr1,expr2   Gene x sample matrices for groups 1 and 2
#' @param times1,times2 Numeric time vectors
#' @param period        Period (default 24)
#' @return list(pval_DR, pval_DP, pval_DM=NA) each numeric[G]
#' @export
detect_DCP_diff <- function(expr1, times1, expr2, times2, period = 24) {
  ngenes     <- nrow(expr1)
  gene_names <- rownames(expr1) %||% paste0("G", seq_len(ngenes))

  res <- tryCatch(
    run_DCP_pipeline(expr1, expr2, times1, times2,
                     gene_names = gene_names, period = period,
                     test_DR = TRUE, test_DP = TRUE),
    error = function(e) NULL
  )

  if (is.null(res)) {
    return(list(pval_DR = rep(1, ngenes), pval_DP = rep(1, ngenes),
                pval_DM = rep(NA_real_, ngenes)))
  }
  list(
    pval_DR = replace(res$p_DR, is.na(res$p_DR), 1),
    pval_DP = replace(res$p_DP, is.na(res$p_DP), 1),
    pval_DM = rep(NA_real_, ngenes)
  )
}


#' LimoRhyde two-group differential rhythmicity
#'
#' limma interaction model: group × cosinor (cos + sin terms).
#' Tests H0: no group × rhythm interaction. DR only.
#'
#' @param expr1,expr2   Gene x sample matrices
#' @param times1,times2 Numeric time vectors
#' @param period        Period (default 24)
#' @return list(pval_DR, pval_DP=NA, pval_DM=NA)
#' @export
detect_LimoRhyde <- function(expr1, times1, expr2, times2, period = 24) {
  if (!requireNamespace("limma", quietly = TRUE))
    stop("Package 'limma' is required for detect_LimoRhyde()")

  ngenes <- nrow(expr1)
  omega  <- 2 * pi / period
  times  <- c(times1, times2)
  g      <- c(rep(0L, ncol(expr1)), rep(1L, ncol(expr2)))
  expr   <- cbind(expr1, expr2)

  cos_t  <- cos(omega * times)
  sin_t  <- sin(omega * times)
  design <- model.matrix(~ g + cos_t + sin_t + g:cos_t + g:sin_t)

  pval_DR <- tryCatch({
    fit      <- limma::lmFit(expr, design)
    fit      <- limma::eBayes(fit)
    int_cols <- grep("g:cos_t|g:sin_t", colnames(design))
    ct       <- limma::topTable(fit, coef = int_cols, number = ngenes,
                                 sort.by = "none", adjust.method = "none")
    p        <- ct$P.Value
    p[is.na(p)] <- 1
    p
  }, error = function(e) rep(1, ngenes))

  list(pval_DR = pval_DR,
       pval_DP = rep(NA_real_, ngenes),
       pval_DM = rep(NA_real_, ngenes))
}


#' DODR two-group differential oscillation detection
#'
#' Wraps DODR::dodr(). DR only.
#'
#' @param expr1,expr2   Gene x sample matrices
#' @param times1,times2 Numeric time vectors
#' @param period        Period (default 24)
#' @return list(pval_DR, pval_DP=NA, pval_DM=NA)
#' @export
detect_DODR <- function(expr1, times1, expr2, times2, period = 24) {
  if (!requireNamespace("DODR", quietly = TRUE))
    stop("Package 'DODR' is required. Install with: install.packages('DODR')")

  ngenes  <- nrow(expr1)
  pval_DR <- tryCatch({
    res <- DODR::dodr(t(expr1), t(expr2),
                      times1 = times1, times2 = times2,
                      period = period, method = "both")
    p   <- res$p.value
    p[is.na(p)] <- 1
    p
  }, error = function(e) rep(1, ngenes))

  list(pval_DR = pval_DR,
       pval_DP = rep(NA_real_, ngenes),
       pval_DM = rep(NA_real_, ngenes))
}


# ==============================================================================
# detect_FMM: K-harmonic LRT for rhythmicity detection
# ==============================================================================
# Rather than fitting the nonlinear FMM model
# y = M + A*cos(beta + 2*atan(omega*tan((t-alpha)/2))) directly, we use its
# closed-form Fourier expansion
#
#   y(t) = M' + A(1-r^2) * sum_{k>=1} r^(k-1) cos(k(t-alpha) + beta),
#   r = (1-omega)/(1+omega) in (0, 1)
#
# Truncating at K harmonics and releasing the geometric-decay constraint
# yields a (2K+1)-parameter LINEAR regression in which all parameters are
# identifiable under both H0 and H1, avoiding the boundary identifiability
# problem of a direct nonlinear FMM fit (where {beta, alpha, omega} become
# unidentified at A=0, the Davies-supremum problem).
#
# At K=1 this reduces to standard cosinor (DCP). K=2 captures >90% of FMM
# variance for omega >= 0.3, K=3 captures >97% for omega >= 0.3.

#' FMM-Harmonic LRT for Rhythmicity Detection
#'
#' Tests for rhythmicity in time-course data using the K-harmonic Fourier
#' decomposition of the FMM (Frequency Modulated Mobius) signal model.
#' Equivalent to a 2K-parameter joint cosinor regression test under
#' \eqn{H_0: a_1 = b_1 = \cdots = a_K = b_K = 0}. With \code{ebayes = TRUE}
#' (default), variance estimates are moderated by \code{limma::eBayes}
#' following the LimoRhyde framework (Singer & Hughey 2019), extended here
#' to user-specified \eqn{K}.
#'
#' @param expr Numeric matrix (genes x samples) of expression values.
#' @param times Numeric vector of sample collection times in hours
#'   (length = ncol(expr)).
#' @param period Period in hours (default 24).
#' @param K Number of harmonics to include (default 2). K=1 reduces to
#'   standard cosinor (DCP). K=2 is the recommended default for non-cosinor
#'   signals; K=3 may improve power for very sharp peaks (omega < 0.15).
#' @param ebayes Logical. If \code{TRUE} (default), variance is moderated
#'   across genes via \code{limma::eBayes}. If \code{FALSE}, the unshrunk
#'   exact F-test is used (useful for null-calibration verification). The
#'   function automatically falls back to the unshrunk path when
#'   \code{nrow(expr) < 50}.
#' @param adjust.method Multiple-testing adjustment for the BH-style
#'   discovery flag (default "BH"). Set to NULL or "none" to skip.
#' @param fdr.threshold FDR cutoff for the \code{discovery} column
#'   (default 0.05). Ignored if \code{adjust.method} is NULL/"none".
#' @param return.coefficients If TRUE, the returned data frame includes
#'   per-gene harmonic coefficients \code{a_1, b_1, ..., a_K, b_K}.
#' @param mc.cores Parallel cores. Currently unused; the linear-regression
#'   path is fully vectorised across genes. Argument retained for API
#'   compatibility with other detectors.
#'
#' @return Data frame with one row per gene:
#'   \describe{
#'     \item{\code{F.stat}}{Numerator-over-denominator F statistic.}
#'     \item{\code{df1}, \code{df2}}{Numerator (\code{2K}) and denominator
#'           (\code{n-2K-1}, or the moderated denominator df when
#'           \code{ebayes = TRUE}) degrees of freedom.}
#'     \item{\code{p.value}}{F-test p-value (moderated when
#'           \code{ebayes = TRUE}).}
#'     \item{\code{p.adjust}}{Adjusted p-value (if \code{adjust.method}
#'           non-null).}
#'     \item{\code{R2}}{Coefficient of determination of the K-harmonic fit.}
#'     \item{\code{discovery}}{Logical: \code{p.adjust < fdr.threshold}.}
#'   }
#'
#' @details
#' The test fits the linear model
#' \deqn{y_t = M + \sum_{k=1}^K [a_k \cos(k\,2\pi t / T) + b_k \sin(k\,2\pi t / T)] + \varepsilon}
#' and tests the joint null \eqn{H_0: a_1 = b_1 = \cdots = a_K = b_K = 0}.
#' With \code{ebayes = TRUE}, fitting and inference proceed via
#' \code{limma::lmFit} on a basis assembled from per-harmonic calls to
#' \code{limorhyde::limorhyde(times, period = T/k)} for \eqn{k = 1, \ldots, K},
#' followed by \code{limma::eBayes} variance shrinkage and
#' \code{limma::topTable} extraction of the joint harmonic F-test. With
#' \code{ebayes = FALSE}, the same model is fit by a vectorised QR
#' decomposition and tested by the exact nested F-statistic with null
#' distribution \eqn{F(2K, n - 2K - 1)} under Gaussianity.
#'
#' This is motivated by the FMM model's closed-form Fourier expansion:
#' the cosine of a Mobius-transformed angle has harmonic amplitudes that
#' decay geometrically with \eqn{r = (1-\omega)/(1+\omega)}. Truncating
#' at K and releasing the decay constraint trades a small bias (the
#' truncation error \eqn{r^{2K}}) for exact identifiability under H0,
#' avoiding the Davies-supremum calibration problem of a direct
#' nonlinear FMM fit.
#'
#' @references
#' Rueda, C., Larriba, Y., Peddada, S. D. (2019). Frequency Modulated Mobius
#'   Model for the Estimation of Rhythmic Signals. Sci Rep 9, 18138.
#'
#' Singer, J. M., Hughey, J. J. (2019). LimoRhyde: a flexible approach for
#'   differential analysis of rhythmic transcriptome data. J Biol Rhythms
#'   34(1), 5-18.
#'
#' Hughes, M. E. et al. (2017). Guidelines for genome-scale analysis of
#'   biological rhythms. J Biol Rhythms 32(5), 380-393.
#'
#' @examples
#' set.seed(1)
#' n  <- 48
#' tt <- seq(0, 24*(1-1/n), length.out = n)
#' rhy <- 2 * cos(2*pi*tt/24) + 0.6 * cos(4*pi*tt/24)   # 1st + 2nd harmonics
#' expr <- rbind(
#'   rhythmic    = rhy + rnorm(n, 0, 0.5),
#'   arrhythmic  = rnorm(n, 0, 0.5)
#' )
#' detect_FMM(expr, tt, K = 2)
#'
#' @export
detect_FMM <- function(expr,
                        times,
                        period = 24,
                        K = 2L,
                        ebayes = TRUE,
                        adjust.method = "BH",
                        fdr.threshold = 0.05,
                        return.coefficients = FALSE,
                        mc.cores = 1L) {
  # Rhythmicity detection by K-harmonic linear regression, adopting the
  # LimoRhyde framework (Singer & Hughey 2019) extended to user-specified
  # K. With ebayes=TRUE (default), variance estimates are moderated by
  # limma::eBayes; this is the contribution we inherit from LimoRhyde. At
  # K=1 this is exactly LimoRhyde's default; at K=2 it is the extension we
  # recommend, motivated by realistic omega distributions in circadian
  # transcriptomes. Set ebayes=FALSE for the unshrunk exact F-test
  # (useful for null calibration verification).
  if (!requireNamespace("limma", quietly = TRUE))
    stop("Package 'limma' required for detect_FMM. ",
         "Install with BiocManager::install('limma').")

  if (!is.matrix(expr)) expr <- as.matrix(expr)
  K <- as.integer(K)
  if (K < 1L) stop("K must be >= 1.")
  n <- ncol(expr)
  G <- nrow(expr)
  if (length(times) != n)
    stop("length(times) must equal ncol(expr).")
  if (n <= 2 * K + 1)
    stop(sprintf("n=%d is too small for K=%d (need n > 2K+1 = %d).",
                 n, K, 2 * K + 1))

  n_unique <- length(unique(times %% period))
  if (n_unique < 2L * K + 1L) {
    warning(sprintf(
      "K=%d harmonic regression requires at least %d distinct timepoints per period (Nyquist condition); only %d unique observed. Test will be conservative; consider K = %d.",
      K, 2L * K + 1L, n_unique, max(1L, (n_unique - 1L) %/% 2L)
    ))
  }

  # Design matrix: intercept + (cos, sin) per harmonic, k = 1..K.
  # The basis is adapted from LimoRhyde's helper limorhyde::limorhyde()
  # (Singer & Hughey 2019), which constructs the k=1 cosinor basis at a
  # given period via `cbind(cos(x/period * 2*pi), sin(x/period * 2*pi))`.
  # We extend this to user-specified K by calling the helper once per
  # harmonic with period = period/k, which yields the (cos, sin) pair at
  # frequency k*omega_0. The result is column-for-column equivalent to
  # what LimoRhyde would produce if called K times.
  if (!requireNamespace("limorhyde", quietly = TRUE))
    stop("Package 'limorhyde' required for the cosinor basis construction. ",
         "Install with install.packages('limorhyde').")

  basis_list <- list()
  for (k in seq_len(K)) {
    bk <- limorhyde::limorhyde(times,
                                colnamePrefix = sprintf("h%d_", k),
                                period   = period / k,
                                sinusoid = TRUE,
                                intercept = FALSE)
    basis_list[[k]] <- bk
  }
  X <- cbind(intercept = 1, do.call(cbind, basis_list))   # n x (2K+1)

  # eBayes variance shrinkage borrows information across genes. With too
  # few genes the prior degenerates; fall back to unshrunk in that case.
  if (isTRUE(ebayes) && G < 50L) {
    warning(sprintf(
      "detect_FMM: ebayes=TRUE requires sufficient genes for stable variance shrinkage; G=%d is below the threshold of 50. Falling back to unshrunk F-test.",
      G))
    ebayes <- FALSE
  }

  fit <- limma::lmFit(expr, design = X)
  harmonic_cols <- setdiff(colnames(X), "intercept")
  coef_idx      <- match(harmonic_cols, colnames(fit$coefficients))

  if (isTRUE(ebayes)) {
    fit_eb  <- limma::eBayes(fit)
    tt      <- limma::topTable(fit_eb, coef = coef_idx,
                                number = Inf, sort.by = "none",
                                adjust.method = "none")
    F.stat  <- tt$F
    p.value <- tt$P.Value
    df1_out <- length(coef_idx)
    # Per-gene moderated denominator df: limma uses df.total = df.residual + df.prior
    # to evaluate the moderated F. Store the per-gene vector so downstream
    # reproduction of p.value via pf(F, df1, df2) gives gene-by-gene equality.
    df2_out <- if (!is.null(fit_eb$df.total)) fit_eb$df.total
               else fit_eb$df.residual + fit_eb$df.prior

    # R^2 via fitted residuals from limma model
    fitted_mat <- fit$coefficients %*% t(X)
    resid_mat  <- expr - fitted_mat
    SSE1       <- rowSums(resid_mat^2)
    SSE0       <- rowSums((expr - rowMeans(expr))^2)
    R2_vec     <- ifelse(SSE0 < 1e-12, NA_real_, 1 - SSE1 / SSE0)
  } else {
    # Unshrunk exact F-test (multi-harmonic cosinor regression without
    # variance shrinkage). Used for null calibration verification.
    Y         <- t(expr)
    qrX       <- qr(X)
    resid_m   <- Y - qr.fitted(qrX, Y)
    SSE1      <- colSums(resid_m^2)
    SSE0      <- colSums((Y - rep(colMeans(Y), each = n))^2)
    SSE0_safe <- pmax(SSE0, 1e-12)
    SSE1_safe <- pmax(SSE1, 1e-12)
    df1_out   <- 2L * K
    df2_out   <- n - 2L * K - 1L
    F.stat    <- ((SSE0_safe - SSE1_safe) / df1_out) / (SSE1_safe / df2_out)
    F.stat[F.stat < 0] <- 0
    p.value   <- pf(F.stat, df1_out, df2_out, lower.tail = FALSE)
    R2_vec    <- ifelse(SSE0 < 1e-12, NA_real_, 1 - SSE1 / SSE0)
  }

  out <- data.frame(F.stat   = F.stat,
                    df1      = df1_out,
                    df2      = df2_out,
                    p.value  = p.value,
                    R2       = R2_vec)

  if (!is.null(adjust.method) && adjust.method != "none") {
    out$p.adjust  <- p.adjust(p.value, method = adjust.method)
    out$discovery <- !is.na(out$p.adjust) & out$p.adjust < fdr.threshold
  }

  if (return.coefficients) {
    coef_df <- as.data.frame(fit$coefficients)
    colnames(coef_df) <- c("M", as.vector(rbind(paste0("a_", seq_len(K)),
                                                paste0("b_", seq_len(K)))))
    out <- cbind(out, coef_df)
  }

  rownames(out) <- rownames(expr)
  attr(out, "ebayes")    <- isTRUE(ebayes)
  attr(out, "framework") <- "limorhyde-extended"
  out
}
