# Load necessary libraries
library(Seurat)
library(tidyverse)
library(viridis)
library(openxlsx)
library(scCustomize)
library(RColorBrewer)

# Define a reversed color palette using the RdBu color scale
color_palette <- colorRampPalette(colors = RColorBrewer::brewer.pal(11, "RdBu"))(100) %>% rev()

# Load the Seurat object
sobj <- readRDS("sobj_embryo_model.rds")

# Load the gene list from a CSV file and remove duplicates
gene_list <- read.table("genelist.csv")$V1
duplicated_genes <- gene_list[duplicated(gene_list)]
cat("Duplicated genes removed:", duplicated_genes, "\n")
gene_list <- unique(gene_list)  # Keep only unique genes

# Scale the Seurat object data using the gene list
sobj <- ScaleData(sobj, features = gene_list)

# Set the cell type identity and order the levels
Idents(sobj) <- "celltype"
celltype_order <- rev(levels(sobj$celltype))  # Reverse the order of cell types
sobj$celltype <- factor(sobj$celltype, levels = celltype_order)
Idents(sobj) <- "celltype"  # Re-assign the identity to reflect the updated levels

# Generate a dot plot for the scaled data with a custom color palette
p <- DotPlot(
  object = sobj, features = gene_list
) + 
  scale_colour_gradientn(colours = color_palette) +  # Apply the custom color palette
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)  # Rotate x-axis labels for better readability
  )

# Save the dot plot to a PDF file
pdf("../plots/Dotplot_embryo_model.pdf", width = 15, height = 5.5)
print(p)
dev.off()