# Load necessary libraries
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(patchwork)
library(tidyverse)
library(openxlsx)

# Create Seurat object from 10X data
sobj <- CreateSeuratObject(counts = Read10X(data.dir = "count/filtered_feature_bc_matrix/"))

# Assign sample identifiers based on barcode suffixes
sobj$sample <- case_when(
  grepl("-1$", colnames(sobj)) ~ "E1",  # Assign 'E1' to barcodes ending in '-1'
  grepl("-2$", colnames(sobj)) ~ "E2",  # Assign 'E2' to barcodes ending in '-2'
  TRUE ~ NA_character_                  # Assign NA for any other case
)
sobj$sample <- factor(sobj$sample, levels = c("E1", "E2"))  # Set factor levels for samples

# Calculate the percentage of mitochondrial gene content per cell
sobj$percent_MT <- PercentageFeatureSet(sobj, pattern = "^MT-")

# Visualize QC metrics using a violin plot
VlnPlot(sobj, features = c("nFeature_RNA", "nCount_RNA", "percent_MT"), ncol = 3, pt.size = 0)

# Calculate upper and lower bounds for filtering based on nFeature_RNA distribution
upper_bound <- 10^(mean(log10(sobj$nFeature_RNA)) + 2 * sd(log10(sobj$nFeature_RNA)))
lower_bound <- 10^(mean(log10(sobj$nFeature_RNA)) - 2 * sd(log10(sobj$nFeature_RNA)))

# Set the threshold for mitochondrial gene percentage
p_MT <- 10

# Calculate and print the percentage of cells filtered out by each criterion
filtered_by_ub <- sum(sobj$nFeature_RNA > upper_bound) / length(sobj$nFeature_RNA) * 100
filtered_by_lb <- sum(sobj$nFeature_RNA < lower_bound) / length(sobj$nFeature_RNA) * 100
filtered_by_MT <- sum(sobj$percent_MT > p_MT) / length(sobj$percent_MT) * 100

print(paste0("Percentage of cells filtered out by upper bound: ", round(filtered_by_ub, 2), "%"))
print(paste0("Percentage of cells filtered out by lower bound: ", round(filtered_by_lb, 2), "%"))
print(paste0("Percentage of cells filtered out by mitochondrial content: ", round(filtered_by_MT, 2), "%"))

# Visualize the distribution of nFeature_RNA with filtering bounds
ggplot(sobj@meta.data, aes(x = nFeature_RNA, color = sample)) +
  geom_density() +
  geom_vline(xintercept = lower_bound, linetype = "dashed") +
  geom_vline(xintercept = upper_bound, linetype = "dashed")

# Filter cells based on nFeature_RNA and mitochondrial content
sobj <- subset(sobj, subset = nFeature_RNA > 2000 & nFeature_RNA < 12000 & percent_MT < p_MT)

# Split the RNA assay by sample for downstream analysis
sobj[["RNA"]] <- split(sobj[["RNA"]], f = sobj$sample)

# Display the updated Seurat object
sobj

# Run the standard Seurat analysis workflow
sobj <- NormalizeData(sobj)              # Normalize the data
sobj <- FindVariableFeatures(sobj)       # Identify highly variable features
sobj <- ScaleData(sobj)                  # Scale the data
sobj <- RunPCA(sobj)                     # Perform Principal Component Analysis (PCA)

# Save the Seurat object for future use
saveRDS(sobj, file = "../data/sobj_E1E2_before_integrate.rds")