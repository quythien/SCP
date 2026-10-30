# Adapted preparation blocks. Set SCP_GSE71620_ROOT to the input folder.
setwd(Sys.getenv("SCP_GSE71620_ROOT", "."))

# Reading in 
library(oligo)
pheno <- read.table("data/GSE71620_Phenotype_GEO.txt", 
                             header = TRUE, 
                             sep = "\t", 
                             stringsAsFactors = FALSE)
cel_files <- list.files("data/GSE71620_extracted", 
                        pattern = "CEL.gz", 
                        full.names = TRUE)

data <- read.celfiles(cel_files)

# Checking data 
class(data)
dim(data)         # 1178100 features × 420 samples
sampleNames(data) # list of CEL files / sample IDs
featureNames(data)[1:10] # first 10 probe IDs
pData(data)         # shows the phenoData table
varLabels(phenoData(data))  # variable names

# Normalization 
#  Background correctionn
#  Quantile normalization
#  Log2 transformation

eset <- rma(data)        # ExpressionSet object
print(paste("Normalized expression matrix:", nrow(eset), "genes ×", ncol(eset), "samples"))
#[1] "Normalized expression matrix: 33297 genes × 420 samples" 

library(hugene11sttranscriptcluster.db)
probe_ids <- rownames(eset)  # Your probe IDs
gene_symbols <- mapIds(hugene11sttranscriptcluster.db, 
                       keys = probe_ids,
                       column = "SYMBOL", 
                       keytype = "PROBEID")

annotated_genes <- !is.na(gene_symbols)
eset_annotated <- eset[annotated_genes, ]
gene_names_annotated <- gene_symbols[annotated_genes]

print(paste("Original:", nrow(eset), "probe sets"))  # Original: 33297 probe sets
print(paste("Annotated:", nrow(eset_annotated), "genes with symbols")) # 22025 genes with symbols

# View(exprs(eset_annotated))
library(dplyr)
expr_data <- exprs(eset_annotated)
gene_annotation <- data.frame(
  ProbeID = rownames(expr_data),
  GeneSymbol = gene_names_annotated,
  stringsAsFactors = FALSE
)

# Calculate IQR for each probe
gene_annotation$IQR <- apply(expr_data, 1, IQR)

# Keep only the probe with highest IQR per gene
final_data <- gene_annotation %>%
  group_by(GeneSymbol) %>%
  slice_max(IQR, n = 1, with_ties = FALSE) %>%
  ungroup()

final_expr <- expr_data[final_data$ProbeID, ]
rownames(final_expr) <- final_data$GeneSymbol
nrow(final_expr) # Final set: 19998
# paper: 20237

############################################################################################
# Real data: final_expr (20252 x 420)
# pheno (210 x 8) : 210 subjects 
# PMI: mean postmortem interval 
# RIN: RNA integrity number  
# We have both BA11: Brodmann area 11 and BA47: Brodmann area 47 
# TOD: in range of -6 to 18 
# Split the data:
# Sixty-four subjects were removed based on these criteria: death was not witnessed because their time of death (TOD) cannot be precisely determined
pheno <- (pheno[!is.na(pheno$TOD), ]) # 146 patients 
pheno$AgeGroup <- ifelse(pheno$Age < 40, "younger",
                         ifelse(pheno$Age >= 60, "older", NA))
table(pheno$AgeGroup) # 37 older and 31 younger, match paper 

col_ids <- sub(".*BA(11|47)-([0-9]+)\\.CEL\\.gz", "\\2", colnames(final_expr))
keep_cols <- col_ids %in% pheno$ID
expr_filtered <- final_expr[, keep_cols]
# Filtered data that is also split 
BA11_expr <- expr_filtered[, grepl("BA11", colnames(expr_filtered))]
BA47_expr <- expr_filtered[, grepl("BA47", colnames(expr_filtered))]

saveRDS(list(expr = BA11_expr, pheno = pheno), "data/BA11_data.rds")
saveRDS(list(expr = BA47_expr, pheno = pheno), "data/BA47_data.rds")
