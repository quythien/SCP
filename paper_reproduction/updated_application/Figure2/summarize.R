# 1. Collect all recorded differential case summaries.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2L)
input <- normalizePath(args[1])
out <- args[2]
paths <- list.files(file.path(input,"cases"),pattern="^per_replicate[.]csv$",recursive=TRUE,full.names=TRUE)
stopifnot(length(paths)==32400L)
x <- do.call(rbind,lapply(paths,read.csv))
stopifnot(nrow(x)==518400L,!anyDuplicated(x[c("case_id","endpoint","target_population","BH_threshold")]))
keys <- c("design","B","N_per_group","endpoint","target_population","BH_threshold")
group_key <- do.call(paste,c(lapply(x[keys],function(z) ifelse(is.na(z),"NA",z)),sep="|"))
groups <- split(seq_len(nrow(x)),group_key)

# 2. Retain validity, power and FDP for every cell.
rows <- lapply(groups,function(i) {
  d <- x[i,]
  stopifnot(nrow(d)==600L)
  valid <- d$valid & is.finite(d$power) & is.finite(d$FDP)
  data.frame(pair="Caudate_Putamen",d[1,keys],expected=600L,received=nrow(d),
    valid=sum(valid),invalid=sum(!valid),missing=0L,
    mean_power=if(all(valid)) mean(d$power) else NA_real_,
    power_MCSE=if(all(valid)) sd(d$power)/sqrt(nrow(d)) else NA_real_,
    mean_FDP=if(all(valid)) mean(d$FDP) else NA_real_,
    FDP_MCSE=if(all(valid)) sd(d$FDP)/sqrt(nrow(d)) else NA_real_)
})
curves <- do.call(rbind,rows)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
write.csv(curves,file.path(out,"curves_all_BH.csv"),row.names=FALSE)

# 3. Interpolate the first 80 percent crossing and round up, without extrapolation.
keys <- c("design","B","endpoint","target_population","BH_threshold")
group_key <- do.call(paste,c(lapply(curves[keys],function(z)ifelse(is.na(z),"NA",z)),sep="|"))
guides <- lapply(split(seq_len(nrow(curves)),group_key),function(i) {
  d <- curves[i,];d <- d[order(d$N_per_group),]
  j <- which(d$mean_power>=.8)[1]
  crossing <- if(is.na(j)) NA_real_ else if(j==1L)d$N_per_group[1] else
    d$N_per_group[j-1]+(.8-d$mean_power[j-1])/(d$mean_power[j]-d$mean_power[j-1])*
      (d$N_per_group[j]-d$N_per_group[j-1])
  data.frame(pair="Caudate_Putamen",d[1,keys],raw_crossing=crossing,ceil_N80=ceiling(crossing))
})
write.csv(do.call(rbind,guides),file.path(out,"N80_all_BH.csv"),row.names=FALSE)
