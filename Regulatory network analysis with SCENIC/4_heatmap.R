# Load necessary libraries
library(RColorBrewer)
library(SCENIC)
library(doParallel)
library(dplyr)
library(AUCell)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(writexl)

# Load SCENIC options and cell information
scenicOptions <- readRDS("int/scenicOptions.Rds")
cellInfo <- readRDS("int/cellInfo.Rds")

# Load the AUC matrix containing regulon activity scores
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
# Ensure no duplicated regulons
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)), ]

# Load regulon selection and binary regulon activity (optional)
regulonSelection <- loadInt(scenicOptions, "aucell_regulonSelection", ifNotExists = "null", verbose = FALSE)
binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_full")

# Convert the AUC matrix to a data frame and scale the data
auc_mat <- getAUC(regulonAUC)
auc_mat_Scaled <- t(scale(t(auc_mat), center = TRUE, scale = TRUE))
df <- as.data.frame(auc_mat_Scaled)

# Load the list of selected regulons
regulon_list <- read.csv("filtered_regulonActivity_filtered.csv")
regulon_list <- regulon_list$Regulon

# Filter the AUC matrix to include only the selected regulons
df <- df[regulon_list, ]

# Define the color palette for the heatmap
my_palette <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(100))
col <- colorRamp2(seq(-2, 2, length = 100), my_palette)

# Define the desired order of CellType
desired_order <- c("D0-fm", "P20-nr", "P20-pr", "P15")

# Convert the CellType column to a factor with the specified levels and reorder rows
cellInfo$CellType <- factor(cellInfo$CellType, levels = desired_order)
cellInfo <- cellInfo[order(cellInfo$CellType), ]

# Verify the reordering of cell types
str(cellInfo)

# Reorder the columns of the data frame to match the cell type order
df <- df[, row.names(cellInfo)]

# Create the heatmap and save it as a PDF
pdf(file = "SCENIC-D0fm-P20NR-P20PR-P15_cell_orderby_celltype.pdf", width = 6, height = 5)

Heatmap(df,
        name = "Regulon activity",
        col = col,
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        show_row_names = FALSE,
        show_column_names = FALSE,
        use_raster = FALSE)

dev.off()