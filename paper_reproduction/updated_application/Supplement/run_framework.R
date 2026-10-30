# 1. Load the frozen simulation engine and paired calibration.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)[1])
bundle <- dirname(normalizePath(script))
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2L)
.libPaths(c(normalizePath(args[1]),.libPaths()))
library(SCP)
source(file.path(bundle,"code/simulate_framework.R"))
psi <- readRDS(file.path(bundle,"data/framework_calibration.rds"))
PERIOD <- 24
NGENES_SIM <- 2000L
N <- c(20,40,60,80,100,120,150,200,250,300)
nsims <- as.integer(Sys.getenv("SCP_REPLICATES","100"))
if(nzchar(Sys.getenv("SCP_N"))) N <- as.numeric(strsplit(Sys.getenv("SCP_N"),",")[[1]])
output <- args[2]
dir.create(output,recursive=TRUE,showWarnings=FALSE)

# 2. Run paired harmonic and residual-noise draws using the recorded seed rule.
results <- setNames(vector("list",length(N)),as.character(N))
status <- list()
for(j in seq_along(N)) {
  results[[j]] <- lapply(seq_len(nsims),function(i) {
    seed <- 20260526+(N[j]*1000+i)*7919
    tryCatch(sim_one_rep(N[j],seed,psi,paired_sigma=TRUE),
      error=function(e)structure(conditionMessage(e),class="try-error"))
  })
  valid <- vapply(results[[j]],is.list,logical(1))
  status[[j]] <- data.frame(N=N[j],requested=nsims,valid=sum(valid),invalid=sum(!valid))
  cat("N",N[j],"valid",sum(valid),"invalid",sum(!valid),"\n")
}

# 3. Keep failed cells explicit instead of treating them as zero power.
saveRDS(results,file.path(output,"paired_framework_results.rds"))
write.csv(do.call(rbind,status),file.path(output,"validity.csv"),row.names=FALSE)
