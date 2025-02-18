library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(patchwork)
library(tidyverse)
library(loupeR)

# Command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Integration method and identifier
integration_method <- args[1]
identifier <- args[2]

# Print the parameters
cat("Integration Method:", integration_method, "\n")
cat("Identifier:", identifier, "\n")

# Construct the name for the new reduction
reduction_name <- paste('integrated', identifier, sep = '.')

# Print the new reduction name
cat("Reduction Name:", reduction_name, "\n")

# File paths
outRDS <- paste0("../data/sobj_E1E2_", identifier, ".rds")
outPDF_QC <- paste0("../data/sobj_E1E2_", identifier, "_QC.pdf")
outPDF <- paste0("../data/sobj_E1E2_", identifier, ".pdf")
outLoupe <- paste0("../output/loupe_E1E2_", identifier)

# Print file paths
cat("Output RDS File:", outRDS, "\n")
cat("Output QC PDF File:", outPDF_QC, "\n")
cat("Output PDF File:", outPDF, "\n")
cat("Output Loupe File:", outLoupe, "\n")

# Read the Seurat object
sobj <- readRDS("sobj_E1E2_before_integrate.rds")

# Perform integration dynamically
sobj <- IntegrateLayers(
  object = sobj, method = integration_method,
  new.reduction = reduction_name,
  verbose = TRUE)

# Re-join layers after integration
sobj[["RNA"]] <- JoinLayers(sobj[["RNA"]])

# Further analysis steps
sobj <- FindNeighbors(sobj, reduction = reduction_name, dims = 1:30)
sobj <- FindClusters(sobj, resolution = seq(0.1, 4, 0.1))
sobj <- RunUMAP(sobj, dims = 1:30, reduction = reduction_name)

# Save plots
plot <- DimPlot(sobj, reduction = "umap", split.by = "sample")
ggsave(outPDF_QC, plot, width = 7, height = 7, units = "in")

plot <- DimPlot(sobj, reduction = "umap", group.by = "RNA_snn_res.0.5")
ggsave(outPDF, plot, width = 7, height = 7, units = "in")

create_loupe_from_seurat(sobj, output_name = outLoupe, force = TRUE)

# Save the integrated Seurat object
cat("Progress: Starting to save the object...\n")
saveRDS(sobj, file = outRDS)
cat("Progress: Object saved successfully.\n")
