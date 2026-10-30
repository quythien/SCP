# Frozen-universe calibration for the bounded width comparison.

Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
  MKL_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
if (requireNamespace("RhpcBLASctl",quietly=TRUE)) RhpcBLASctl::blas_set_num_threads(1L)
file_hash <- function(path) digest::digest(file=path,algo="sha256")
object_hash <- function(x) digest::digest(x,algo="sha256")

fit_width_pilot <- function(expression,times,genes) {
  stopifnot(ncol(expression)==length(times),nrow(expression)==length(genes),
    all(is.finite(expression)),all(is.finite(times)),length(times)>3L)
  n <- length(times)
  X <- cbind(1,cos(2*pi*times/24),sin(2*pi*times/24))
  fit <- lm.fit(X,t(expression))
  if (fit$rank!=3L) stop("rank_deficient_pilot")
  rss <- colSums(fit$residuals^2)
  tss <- rowSums((expression-rowMeans(expression))^2)
  sigma <- sqrt(rss/(n-3))
  p <- rep(NA_real_,nrow(expression))
  valid <- tss>0 & sigma>0 & is.finite(sigma)
  p[valid] <- pf((pmax(tss[valid]-rss[valid],0)/2)/(rss[valid]/(n-3)),
    2,n-3,lower.tail=FALSE)
  data.frame(gene=genes,M=fit$coefficients[1,],
    A=sqrt(fit$coefficients[2,]^2+fit$coefficients[3,]^2),
    sigma=sigma,phi=(atan2(fit$coefficients[3,],fit$coefficients[2,])*12/pi)%%24,
    pvalue=p,row.names=NULL)
}

calibrate_width_pilot <- function(fit) {
  passing <- which(is.finite(fit$pvalue)&fit$pvalue<.01&fit$A>0&
    is.finite(fit$A)&is.finite(fit$sigma)&fit$sigma>0&is.finite(fit$phi))
  selected <- head(passing[order(fit$pvalue[passing])],300L)
  metric <- data.frame(genes_fitted=nrow(fit),n_passing=length(passing),
    n_selected=length(selected),prop_rhythmic=length(passing)/nrow(fit),
    median_r=if(length(selected)) median(fit$A[selected]/fit$sigma[selected]) else NA_real_,
    amplitude_floor_fraction=if(length(selected)) mean(fit$A[selected]<.05) else NA_real_)
  bio <- if(!length(selected)) NULL else list(ngenes=5000L,prop_rhythmic=metric$prop_rhythmic,
    period=24,lBaselineExpr=fit$M,lOD=log(fit$sigma[is.finite(fit$sigma)&fit$sigma>0]),
    amplitude=fit$A[selected],sigma_rhythmic=fit$sigma[selected],
    phase=fit$phi[selected],paired_sigma=TRUE)
  list(metric=metric,bio=bio,selected=selected,undefined_gene_count=sum(!is.finite(fit$pvalue)))
}

width_donor_seed <- function(case,outer) as.integer(1200000000+(case-1)*100+outer)
width_inner_seed <- function(case,outer,grid,replicate) {
  stopifnot(case>=1,case<=47,outer>=0,outer<=30,grid>=1,grid<=90,
    replicate>=1,replicate<=100)
  as.integer(1800000000+case*1000000+outer*20000+grid*200+replicate)
}

save_new <- function(x,path) {
  stopifnot(!file.exists(path))
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  temporary <- paste0(path,".partial_",Sys.getpid())
  stopifnot(!file.exists(temporary))
  saveRDS(x,temporary,compress="gzip")
  stopifnot(file.rename(temporary,path))
}
