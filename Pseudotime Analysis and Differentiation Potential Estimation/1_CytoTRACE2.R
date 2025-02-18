# Load necessary libraries
library(CytoTRACE2)
library(RColorBrewer)
library(tidyverse)

# Load the Seurat object
seurat_obj <- readRDS("sobj_D0_D4_D8_D20_P15_mnn.rds")

# Run CytoTRACE2 analysis on the Seurat object
# - `is_seurat = TRUE` indicates that the input object is a Seurat object
# - `slot_type = "counts"` specifies that the gene expression data is stored in the "counts" slot
# - `batch_size` and `smooth_batch_size` control the size of batches during processing, for large datasets
# - `species = "human"` indicates that the dataset corresponds to human species
cytotrace2_result <- cytotrace2(
  seurat_obj, 
  is_seurat = TRUE, 
  slot_type = "counts", 
  batch_size = 10000, 
  smooth_batch_size = 1000,
  species = "human"
)

# Save the CytoTRACE2 result for future use
saveRDS(cytotrace2_result, file = "sobj_cytotrace2_D0_D4_D8_D20_P15_mnn.rds")
