# 1. Read the full-cohort bootstrap summaries.
source("code/onecol_common.R")
data_dir <- Sys.getenv("SCP_BOOTSTRAP_DATA", "Figure3/data")
source("Figure3/code/panels.R")
values <- read.csv(file.path(data_dir,"panel_values.csv"))
CP <- "#0072B2"
CB <- "#D55E00"

# 2. Preserve the original Figure2.3 geometry and full displayed ranges.
mai_A <- c(0.36,0.46,0.24,0.04)
mai_B <- c(0.36,0.08,0.24,0.04)
cell_B <- (COL_W-((mai_A[2]+mai_A[4])-(mai_B[2]+mai_B[4])))/2
path <- file.path(OUT_DIR,"Figure3.pdf")
pdf(path,width=COL_W,height=2.05,pointsize=PT)
layout(matrix(c(1,2,3,3),nrow=2,byrow=TRUE),
  widths=c(COL_W-cell_B,cell_B),heights=c(1.8,0.25))
base_par(mai_A)
draw_boot(values[values$panel=="A",],"A","Caudate (n = 59)")
cx_A <- plot_centre()
base_par(mai_B)
draw_boot(values[values$panel=="B",],"B","Human BA47 (n = 146)",ylab=FALSE)
cx_B <- plot_centre()
par(mai=c(0.02,0.03,0.02,0.03),xpd=NA)
plot.new()
legend(grconvertX(mean(c(cx_A,cx_B)),"ndc","user"),0.5,xjust=0.5,yjust=0.5,
  c("Projected power (plug-in)","Bootstrap mean and 95% interval"),
  col=c(CP,CB),lty=c(1,2),lwd=c(1.6,1.0),pch=c(NA,19),pt.cex=PT_CEX+0.05,
  ncol=2,bty="n",cex=0.8,seg.len=1.6,x.intersp=0.5)
par(xpd=FALSE)
dev.off()

