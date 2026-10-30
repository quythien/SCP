# 1. Read a prepared full-cohort expression matrix and its recorded configuration.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)[1])
caller <- getwd()
setwd(dirname(dirname(normalizePath(script))))
source("code/exact_F.R")
source("Figure3/code/calibrate.R")
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)>=3L, args[1] %in% c("Caudate","BA47"))
arm <- args[1]
cohort_path <- normalizePath(if (grepl("^/",args[2])) args[2] else file.path(caller,args[2]), mustWork=TRUE)
output <- args[3]
configuration <- readRDS("Figure3/data/configuration.rds")[[arm]]
cohort <- readRDS(cohort_path)
stopifnot(object_hash(unname(cohort$expression)) == configuration$expression_sha256,
  object_hash(as.numeric(cohort$TOD)) == configuration$time_sha256,
  object_hash(cohort$genes) == configuration$genes_sha256)
stopifnot(ncol(cohort$expression) == configuration$n)
outer <- if(length(args)>3L) as.integer(args[-(1:3)]) else 0:100
stopifnot(all(outer %in% 0:100))
dir.create(output,recursive=TRUE,showWarnings=FALSE)

# 2. Resample donors and refit every gene, retaining the expression-time pairs.
for (b in outer) {
  path <- file.path(output,sprintf("%s_outer_%03d.rds",arm,b))
  stopifnot(!file.exists(path))
  index <- seq_len(configuration$n)
  if (b>0L) {
    set.seed(as.integer(configuration$donor_seed_base+b),kind="L'Ecuyer-CMRG")
    index <- sample.int(configuration$n,configuration$n,replace=TRUE)
  }
  times <- cohort$TOD[index]
  fit <- fit_width_pilot(cohort$expression[,index,drop=FALSE],times,cohort$genes)
  calibrated <- calibrate_width_pilot(fit)
  sampler <- scp_time_sampler(times)

  # 3. Evaluate each power curve with 100 inner experiments per sample size.
  points <- vector("list",length(configuration$N))
  for (g in seq_along(configuration$N)) {
    N <- configuration$N[g]
    values <- vapply(seq_len(100L),function(i) {
      seed <- as.integer(configuration$inner_seed_base+b*20000+g*200+i)
      result <- simulate_scp_f_replicate(calibrated$bio,NULL,N,seed,time_sampler=sampler)
      c(power=result$power,FDP=result$FDP)
    },numeric(2))
    points[[g]] <- data.frame(N=N,power=mean(values["power",]),
      power_mcse=sd(values["power",])/10,mean_FDP=mean(values["FDP",]),
      FDP_mcse=sd(values["FDP",])/10)
  }
  saveRDS(list(arm=arm,outer=b,points=do.call(rbind,points),
    calibration=calibrated$bio,sampled_unit_index=index),path)
  cat(arm,"outer",b,"complete\n")
}
