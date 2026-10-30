# 1. Read the plug-in curve and all 100 outer curves for both cohorts.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==3L)
values <- list()
for(i in 1:2) {
  arm <- c("Caudate","BA47")[i]
  paths <- file.path(args[i],sprintf("%s_outer_%03d.rds",arm,0:100))
  stopifnot(all(file.exists(paths)))
  draws <- lapply(paths,readRDS)
  N <- draws[[1]]$points$N
  stopifnot(all(vapply(draws,function(x)identical(x$points$N,N),logical(1))))
  power <- vapply(draws[-1],function(x)x$points$power,numeric(length(N)))

  # 2. Report the bootstrap mean and percentile interval at each N.
  values[[i]] <- data.frame(panel=LETTERS[i],N=N,plugin=draws[[1]]$points$power,
    mean=rowMeans(power),lo=apply(power,1,quantile,probs=.025),
    hi=apply(power,1,quantile,probs=.975))
}
dir.create(args[3],recursive=TRUE,showWarnings=FALSE)
write.csv(do.call(rbind,values),file.path(args[3],"panel_values.csv"),row.names=FALSE)
