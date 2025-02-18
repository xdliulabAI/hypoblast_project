# Load necessary libraries
library(Seurat)
library(Matrix)
library(tidyverse)
library(openxlsx)

# Load the Seurat object and update it to the latest version
sobj <- readRDS("readySeu_hrpi.rds")
sobj <- UpdateSeuratObject(sobj)

# Subset the Seurat object to include only specific libraries
sobj <- subset(sobj, subset = library %in% c("D0-fm", "P20-pr", "P20-nr"))

# Update the cell type information in the metadata
sobj$celltype <- sobj$library
md <- sobj@meta.data
sobj@meta.data <- md[, c("orig.ident", "celltype")]

# Convert the RNA count data to a matrix format and update row names based on gene IDs
data_matrix <- as.matrix(sobj@assays$RNA@counts)
name2id_df <- data.frame(GeneName = names(sobj@misc$convert$name2id), 
                         GeneID = sobj@misc$convert$name2id)

# Subset the data matrix to match the gene IDs
data_matrix <- data_matrix[name2id_df$GeneID, ]
row.names(data_matrix) <- name2id_df$GeneName

# Create a new Seurat object with the updated matrix and metadata
sobj <- CreateSeuratObject(counts = data_matrix, meta.data = md)

# Load a second Seurat object and subset it based on specific cluster resolution
P15 <- readRDS("sobj_12_14_clean_harmony.rds")
P15 <- subset(P15, subset = RNA_snn_res.0.05 %in% c("0"))

# Update cell type information and clean up unnecessary layers
P15$celltype <- "P15"
md <- P15@meta.data
P15@meta.data <- md[, c("orig.ident", "celltype")]
P15@assays[["RNA"]]@layers[["data"]] <- NULL
P15@assays[["RNA"]]@layers[["scale.data"]] <- NULL
P15@assays[["RNA"]]@layers[["scale.data.1"]] <- NULL
P15@assays[["RNA"]]@layers[["scale.data.2"]] <- NULL

# Merge the two Seurat objects into one
sobj <- merge(sobj, y = c(P15))
md <- sobj@meta.data

# Join RNA layers after merging the objects
sobj[["RNA"]] <- JoinLayers(sobj[["RNA"]])

# Select a random subset of 10,000 cells for downstream analysis
all_barcodes <- colnames(sobj)
set.seed(123) # Setting seed for reproducibility
selected_barcodes <- sample(all_barcodes, 10000)
sobj <- subset(sobj, cells = selected_barcodes)

# Set cell identities to the cell type
Idents(sobj) <- "celltype"

# Extract the count matrix and metadata for saving
mat <- sobj@assays[["RNA"]]@layers[["counts"]]
row.names(mat) <- row.names(sobj)
colnames(mat) <- colnames(sobj)
cellInfo <- md[, c("orig.ident", "celltype")]
colnames(cellInfo)[2] <- "CellType"

# Save the count matrix and cell information as RDS files
saveRDS(mat, file = "exprMat_D0Fm-P20NR-P20PR-P15.Rds")
saveRDS(cellInfo, file = "cellInfo_D0Fm-P20NR-P20PR-P15.Rds")
