# Load necessary libraries
library(RColorBrewer)
library(SCENIC)
library(doParallel)
library(dplyr)

# Load the SCENIC options from the previous step
scenicOptions <- readRDS("int/scenicOptions.Rds")

# Load GRNBoost output and process it for SCENIC
GRNBoost_output <- read.delim("adj.tsv")
colnames(GRNBoost_output) <- c("TF", "Target", "weight")
saveRDS(GRNBoost_output, file="int/1.4_GENIE3_linkList.Rds")

# Load the filtered expression matrix
exprMat_filtered <- readRDS("int/exprMat_filtered.rds")

# Step 1: Run correlation analysis on the filtered expression matrix
runCorrelation(exprMat_filtered, scenicOptions)

# Step 2: Identify co-expression modules
runSCENIC_1_coexNetwork2modules(scenicOptions)

# Load motif annotations for human genes
data(list="motifAnnotations_hgnc_v9", package="RcisTarget")
motifAnnotations_hgnc <- motifAnnotations_hgnc_v9

# Step 3: Create regulons based on co-expression and motif enrichment
runSCENIC_2_createRegulons(scenicOptions)

# Update the number of cores to 1 for the next steps
scenicOptions@settings$nCores <- 1

# Load the original expression matrix and log-transform the data
exprMat <- readRDS(file = "exprMat_D0Fm-P20NR-P20PR-P15.Rds")
exprMat <- as.matrix(exprMat)
exprMat_log <- log2(exprMat + 1)

# Step 4: Score cells based on regulon activity using the log-transformed expression matrix
runSCENIC_3_scoreCells(scenicOptions, exprMat_log)

# Step 5: Binarize the AUCell results to identify active regulons per cell
runSCENIC_4_aucell_binarize(scenicOptions)

# Load cell information
cellInfo <- readRDS("int/cellInfo.Rds")

# Load the AUC matrix containing regulon activity scores
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)), ]

# Optionally load additional information for downstream analysis
regulonSelection <- loadInt(scenicOptions, "aucell_regulonSelection", ifNotExists = "null", verbose = FALSE)
cellInfo <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, "cellInfo"), ifNotExists = "null")
cellInfo <- data.frame(cellInfo)
colVars <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, "colVars"), ifNotExists = "null")

# Load binarized regulon activity (optional step depending on your analysis needs)
binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_full")

# Re-load regulonAUC for further analysis (ensuring no duplicates)
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)), ]

# Calculate the average regulon activity by cell type and scale the results
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$CellType),
                                     function(cells) rowMeans(getAUC(regulonAUC)[, cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = TRUE, scale = TRUE))

# Generate a heatmap of scaled regulon activity by cell type
pdf(file = "SCENIC_regulonActivity_byCellType_Scaled.pdf", width = 12, height = 48)
ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name = "Regulon activity")
dev.off()