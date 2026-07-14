#' Fig 3 generator: Efron SUBJECT bootstrap, single-cohort power uncertainty.
#' Caudate control (GSE160521, n=59) vs GTEx skeletal muscle (n=748).
#' Reports the point estimate (plug-in) + bootstrap mean + the 2.5-97.5 percentile
#' 95% CI band (matching runBootstrapDesignGrid's percentile CI). B_out=50, N_sim_in=25.
#' Requires controlled-access raw pilot matrices (server-only). Activate the fast
#' C++ path before running (see FIGURES.md).
suppressWarnings(suppressMessages({ library(Rcpp); library(SCP); library(parallel); library(readxl) }))
dyn.load(system.file("libs","SCP.so",package="SCP")); .CPP_LOADED <- TRUE
NGENES<-5000L; NBOOT<-50L; NSIMS<-25L; NSIMS_INNER<-25L; PVAL<-0.01; SEED<-2025L; CORES<-48L
analysis<-CircadianAnalysisOptions(alpha=0.05,p.adjust.method="BH",fdr_thresholds=0.05)
pc<-function(bio,N,cts,ns,cr=1L){d<-CircadianDesignOptions(sample_sizes=N,nsims=ns,design="passive",cts=cts)
  rowMeans(runSimsSingleCohort(bio,d,analysis,method="cosinor",K=1L,verbose=FALSE,mc.cores=cr)$marginal_power,na.rm=TRUE)}
efron<-function(mat,tod,N,label){
  pp<-prepCircadianData(mat,times=tod,input_type="log2");mat<-pp$data;tod<-pp$times
  set.seed(SEED);if(nrow(mat)>NGENES)mat<-mat[sample(nrow(mat),NGENES),,drop=FALSE];n<-ncol(mat)
  bio0<-estCircadianParam(mat,tod,min_rhythm_pval=PVAL,paired_sigma=TRUE,verbose=FALSE);bio0$ngenes<-NGENES
  set.seed(SEED);plugin<-pc(bio0,N,tod,NSIMS,CORES)
  S<-do.call(cbind,mclapply(seq_len(NBOOT),function(b){set.seed(SEED+b);j<-sample(n,n,replace=TRUE)
    r<-tryCatch({bio<-estCircadianParam(mat[,j,drop=FALSE],tod[j],min_rhythm_pval=PVAL,paired_sigma=TRUE,verbose=FALSE)
      bio$ngenes<-NGENES;pc(bio,N,tod[j],NSIMS_INNER,1L)},error=function(e)NULL)
    if(is.null(r))rep(NA,length(N)) else r},mc.cores=CORES))
  list(label=label,n=n,N=N,plugin=plugin,mean=rowMeans(S,na.rm=T),
       lo=apply(S,1,quantile,.025,na.rm=T),hi=apply(S,1,quantile,.975,na.rm=T))}  # 2.5-97.5 percentile CI
# Controlled-access sources (CAMO server): GSE160521 striatal CPM + matched-donor
# clinical file; GTEx CPM.all.norm (dbGaP). Not redistributed; override via env vars.
KD<-Sys.getenv("SCP_GSE160521_DIR", unset="data/gse160521")
clin<-read.csv(file.path(KD,"DS_clinical_1221_rm97_rm231_matchIndex34.csv"),row.names=1);ctl<-clin[clin$Diagnostic.Category=="CONTROL",]
f<-list.files(KD,pattern="^Caudate_CPMfiltered_logCPM.*csv$",full.names=TRUE);ex<-as.matrix(read.csv(f[1],row.names=1,check.names=FALSE))
cols<-intersect(ctl$pair,colnames(ex));cd<-efron(ex[,cols,drop=FALSE],ctl$CorrectedTOD[match(cols,ctl$pair)]%%24,c(20L,40L,60L,80L,100L,120L,160L,220L),"Caudate control (n=59)")
load(Sys.getenv("SCP_GTEX_CPM_NORM", unset="CPM.all.norm.RData"))
df<-CPM.all.norm[["Muscle - Skeletal"]];ids<-as.character(colnames(df));hh<-sapply(strsplit(ids,"\\."),function(z)if(length(z)>=3)z[3] else NA)
hr<-suppressWarnings(as.numeric(substr(hh,1,2))+as.numeric(substr(hh,3,4))/60);ok<-!is.na(hr)
ms<-efron(as.matrix(df[,ok]),hr[ok],c(40L,80L,120L,160L,200L,240L,280L,320L),"GTEx Muscle-Skeletal (n=748)")
cp<-"#0072B2";cb<-"#D55E00"
pan<-function(p,lab,main,leg,xmax){N<-p$N;keep<-N<=xmax;N<-N[keep]
  plot(N,100*p$plugin[keep],type="l",lwd=2.4,col=cp,ylim=c(0,100),xlim=range(N),xlab="Sample size (n)",ylab="Power (%)",main="",xaxt="n")
  axis(1,at=N,cex.axis=0.8,gap.axis=-1);abline(h=80,lty=2,col="grey60")
  title(main=sprintf("%s   %s",lab,main),line=0.5,cex.main=1.02,font.main=2)
  lines(N,100*p$mean[keep],lwd=1.4,col=cb,lty=2)
  arrows(N,100*p$lo[keep],N,100*p$hi[keep],code=3,angle=90,length=0.035,lwd=1.7,col=cb);points(N,100*p$mean[keep],pch=19,col=cb,cex=0.7)
  if(leg)legend("bottomright",c("Projected power","95% bootstrap CI"),col=c(cp,cb),lwd=c(2.4,1.7),lty=c(1,2),pch=c(NA,19),bty="o",box.col="grey65",cex=0.6,inset=0.02,y.intersp=0.9,seg.len=1.4)}
for(dest in c("submission/figures/Fig3_bootstrap_singlecohort.pdf","output/main_figures/Fig3_bootstrap_singlecohort.pdf")){
  dir.create(dirname(dest),recursive=TRUE,showWarnings=FALSE)
  cairo_pdf(dest,width=7.8,height=3.8);par(mfrow=c(1,2),mar=c(4,4.2,2.6,1.2),mgp=c(2.5,0.7,0),oma=c(0,0,2,0))
  pan(cd,"A","Caudate control (n=59)",FALSE,220);pan(ms,"B","GTEx Muscle-Skeletal (n=748)",TRUE,320)
  mtext("Bootstrap Uncertainty in Single-Cohort Power Estimates",outer=TRUE,side=3,line=0.3,font=2,cex=1.05);dev.off()}
