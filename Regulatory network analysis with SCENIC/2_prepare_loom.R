# Load necessary libraries
library(RColorBrewer)
library(SCENIC)
library(doParallel)
library(dplyr)
library(SCopeLoomR)

# Set the number of cores for parallel processing
cores <- 8

# Load the expression matrix and cell information
exprMat <- readRDS(file = "exprMat_D0Fm-P20NR-P20PR-P15.Rds")
cellInfo <- readRDS(file = "cellInfo_D0Fm-P20NR-P20PR-P15.Rds")

# Display basic information about the data
cat("Dimensions of the expression matrix:", dim(exprMat), "\n")
cat("First few rows of the cell information:\n")
print(head(cellInfo))
str(cellInfo)

# Convert CellType to a factor and drop any unused levels
cellInfo$CellType <- factor(cellInfo$CellType)
cellInfo$CellType <- droplevels(cellInfo$CellType)
unique_cell_types <- unique(cellInfo$CellType)

# Create a directory for intermediate files
dir.create("int", showWarnings = FALSE)

# Save the processed cell information
saveRDS(cellInfo, file = "int/cellInfo.Rds")

# Generate a color palette for the cell types
num_colors <- length(unique_cell_types)
colors <- brewer.pal(min(num_colors, 8), "Set1") # Adjust color palette based on number of cell types
if (num_colors > 8) {
  colors <- colorRampPalette(brewer.pal(8, "Set1"))(num_colors)
}

# Assign colors to each cell type
colVars <- list(CellType = setNames(colors, unique_cell_types))

# Save the color variables
saveRDS(colVars, file = "int/colVars.Rds")

# Set up SCENIC options
org <- "hgnc" # Specify the organism: "hgnc" for human
dbDir <- "cisTarget_databases/hg38" # Directory containing the cisTarget databases
myDatasetTitle <- "SCENIC" # A descriptive title for your analysis
data(defaultDbNames)
defaultDbNames$hgnc["500bp"] <- "hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather"
defaultDbNames$hgnc["10kb"] <- "hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather"
data(list="motifAnnotations_hgnc_v9", package="RcisTarget")
motifAnnotations_hgnc <- motifAnnotations_hgnc_v9

# Initialize SCENIC options
scenicOptions <- initializeScenic(
  org = org,
  dbDir = dbDir,
  dbs = defaultDbNames[["hgnc"]],
  datasetTitle = myDatasetTitle,
  nCores = cores
)
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$seed <- 123

# Link the processed cell information and color variables to SCENIC options
scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
scenicOptions@inputDatasetInfo$colVars <- "int/colVars.Rds"

# Save SCENIC options for future use
saveRDS(scenicOptions, file = "int/scenicOptions.Rds")

# Filter genes for the co-expression network analysis
exprMat <- as.matrix(exprMat)
genesKept <- geneFiltering(
  exprMat,
  scenicOptions = scenicOptions,
  minCountsPerGene = 3 * 0.01 * ncol(exprMat),
  minSamples = ncol(exprMat) * 0.01
)
exprMat_filtered <- exprMat[genesKept, ]
cat("Dimensions of the filtered expression matrix:", dim(exprMat_filtered), "\n")
saveRDS(exprMat_filtered, file = "int/exprMat_filtered.rds")

# Build the loom file for downstream analysis with GRNBoost
loom <- build_loom("exprMat_filtered_D0Fm-P20NR-P20PR-P15.loom", dgem = exprMat_filtered)
close_loom(loom)

# Next steps: run GRNBoost on the loom file