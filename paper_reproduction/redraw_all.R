# Redraw the paper figures from the supplied plotting data.

# 1. Check that this is the archive root ----------------------------------
stopifnot(file.exists("code/onecol_common.R"), file.exists("Figure1/data/panel_data.rds"))

# 2. Draw the five main figures -------------------------------------------
source("Figure1/redraw.R")
source("Figure2/redraw.R")
source("Figure3/redraw.R")
source("Figure4/redraw.R")
source("Figure5/redraw.R")

# 3. Draw the supplementary figure and write the tables --------------------
source("Supplement/FigureS1/redraw.R")
source("Supplement/Tables/redraw.R")
cat("Paper figure outputs are in generated/\n")
