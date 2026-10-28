# Write the supplied supplementary tables as LaTeX fragments.

# 1. Read the table data ---------------------------------------------------
table_s1 <- read.csv("Supplement/Tables/data/TableS1.csv", check.names = FALSE)
table_s2 <- read.csv("Supplement/Tables/data/TableS2.csv", check.names = FALSE)
stopifnot(nrow(table_s1) == 4L, nrow(table_s2) == 13L)

# 2. Write the table source ------------------------------------------------
dir.create("generated", showWarnings = FALSE)
stopifnot(file.copy("Supplement/Tables/code/TableS1_source.tex", "generated/TableS1.tex", overwrite = TRUE))
stopifnot(file.copy("Supplement/Tables/code/TableS2_source.tex", "generated/TableS2.tex", overwrite = TRUE))
cat("Saved generated/TableS1.tex and generated/TableS2.tex\n")
