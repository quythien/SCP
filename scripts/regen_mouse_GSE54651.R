old_wd <- setwd("code"); source("setup.R"); setwd(old_wd)

make_row <- function(p, times, species, dataset, tissue, design, condition) {
  if (is.null(p)) return(NULL)
  pvals<-p$raw$pvalue; qvals<-p.adjust(pvals,"BH"); r_all<-p$raw$r; phi_all<-p$raw$phi
  ok<-!is.na(r_all)&!is.na(pvals)&is.finite(r_all)&r_all>0
  idx<-order(pvals[ok])[seq_len(min(300L,sum(ok)))]; rv<-r_all[ok][idx]; rv<-rv[is.finite(rv)]
  rhy_ok<-!is.na(pvals)&pvals<0.05&!is.na(phi_all); phi_rhy<-phi_all[rhy_ok]%%24
  t_mod<-times%%24
  data.frame(species=species,dataset=dataset,tissue=tissue,design=design,condition=condition,
    n=length(times),ngenes=length(pvals),
    tod_min=round(min(t_mod),1),tod_max=round(max(t_mod),1),tod_sd=round(sd(t_mod),2),
    r_median_top300=round(if(length(rv)>=3)median(rv) else NA_real_,3),
    r_q25_top300=round(if(length(rv)>=3)quantile(rv,.25,names=FALSE) else NA_real_,3),
    r_q75_top300=round(if(length(rv)>=3)quantile(rv,.75,names=FALSE) else NA_real_,3),
    phi_median_rhy=round(if(length(phi_rhy)>=3)median(phi_rhy) else NA_real_,2),
    phi_q25_rhy=round(if(length(phi_rhy)>=3)quantile(phi_rhy,.25,names=FALSE) else NA_real_,2),
    phi_q75_rhy=round(if(length(phi_rhy)>=3)quantile(phi_rhy,.75,names=FALSE) else NA_real_,2),
    rhy_FDR20=sum(qvals<.20,na.rm=TRUE),rhy_FDR10=sum(qvals<.10,na.rm=TRUE),
    rhy_FDR05=sum(qvals<.05,na.rm=TRUE),rhy_FDR01=sum(qvals<.01,na.rm=TRUE),
    rhy_p05=sum(pvals<.05,na.rm=TRUE),rhy_p01=sum(pvals<.01,na.rm=TRUE),
    rhy_p001=sum(pvals<.001,na.rm=TRUE),stringsAsFactors=FALSE)
}

d <- readRDS("data/mice_GSE54651_CPM.RData")
cache_dir <- "output/supp_tissue_summary/cache"
tissues <- names(d$count_clean)
rows <- list()
for (tn in tissues) {
  times <- d$tod[[tn]]
  if (length(times) == 0) { cat(sprintf("  [skip] %s: no TOD\n", tn)); next }
  expr <- as.matrix(d$count_clean[[tn]])
  p <- tryCatch(estimate_circadian_params(expr, times, verbose=FALSE),
                error=function(e){message("ERR ",tn,": ",e$message);NULL})
  if (!is.null(p)) saveRDS(p, file.path(cache_dir, sprintf("mouse_GSE54651_%s.rds", tn)))
  r <- make_row(p, times, "Mouse","GSE54651",tn,"active","All")
  if (!is.null(r)) { rows[[tn]] <- r; cat(sprintf("  [ok] %-5s n=%d r_med=%s rhyFDR05=%d\n",
        tn, r$n, r$r_median_top300, r$rhy_FDR05)) }
}
out <- do.call(rbind, rows)
write.csv(out, "/tmp/mouse_rows.csv", row.names=FALSE)
cat("\nwrote", nrow(out), "mouse rows to /tmp/mouse_rows.csv\n")
