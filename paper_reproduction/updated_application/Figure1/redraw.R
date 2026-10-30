# 1. Read the panel results and original style.
source("code/onecol_common.R")
data_dir <- Sys.getenv("SCP_FIGURE1_DATA", "Figure1/data")
source("Figure1/code/panels.R")
current <- readRDS(file.path(data_dir,"panel_data.rds"))
left <- current$NAc
right <- current$Caudate
palette <- c("20"="#CC5252", "40"="#C1CC52", "80"="#52ABCC", "120"="#7E52CC")
left$size_cols <- unname(palette[as.character(left$ss_bc)])
right$size_cols <- unname(palette[as.character(right$ss_bc)])
discovery <- readRDS(file.path(data_dir,"discovery_panel_data.rds"))
stopifnot(!anyNA(left$size_cols), identical(left$size_cols, right$size_cols))

# 2. Draw the power, signal-to-noise ratio and discovery panels.
# Use the original one-column geometry, typography, labels and line weights.
YLAB_LINE <- 2.5
LEFT <- .54
pdf(file.path(OUT_DIR, "Figure1.pdf"), width = COL_W, height = 5, pointsize = PT)
layout(matrix(c(1,2,3,3,4,5,6,7,8,8), nrow = 5, byrow = TRUE),
       heights = c(1.24,.22,1.6,1.6,.34))
par(oma = c(0,0,1.3,0))
base_par(c(.34,LEFT,.20,.05))
draw_A(left, "A"); cx_left <- plot_centre()
draw_A(right, "D"); cx_right <- plot_centre()
center <- mean(c(cx_left,cx_right))
legend_strip(left$thresh_labels,left$thresh_cols,pch=19,ncol=4,at=center)
base_par(c(.66,LEFT,.20,.05))
draw_B(left,"B"); draw_B(right,"E")
draw_C(discovery[[1]],"C"); draw_C(discovery[[2]],"F")
legend_strip(paste0("N=",left$ss_bc),left$size_cols,pch=NA,ncol=4,at=center)
mtext("Human NAc",side=3,outer=TRUE,at=cx_left,line=.1,font=2,cex=1)
mtext("Human Caudate",side=3,outer=TRUE,at=cx_right,line=.1,font=2,cex=1)
dev.off()
