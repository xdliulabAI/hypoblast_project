# Load necessary libraries
library(Seurat)
library(patchwork)
library(cowplot)
library(viridis)
library(tidyverse)
library(pheatmap)
library(nichenetr)
library(anndata)
library(Matrix)
library(openxlsx)

# Load the Seurat object
sobj <- readRDS("sobj_embryo_model.rds")

# Load and process the DEG data
deg_data <- read.csv("DEG_embryo_model.csv")

# Extract the top 50 genes from each cluster based on average log2 fold change
top50_genes_per_cluster <- deg_data %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 50) %>%
  ungroup()

# Check for discrepancies between the top genes and the Seurat object
missing_genes <- setdiff(top50_genes_per_cluster$gene, row.names(sobj))
if(length(missing_genes) > 0) {
  cat("Warning: The following genes are missing in the Seurat object:\n", missing_genes, "\n")
}

# Scale data using the top 50 genes per cluster
sobj <- ScaleData(sobj, features = top50_genes_per_cluster$gene)

# Calculate average expression across clusters
avgexp <- AverageExpression(sobj, return.seurat = TRUE)

# Define the cluster names in the desired order
cluster_names <- c("Am", "Early.meso", "Epi", "ExEM 1", "ExEM 2", "ExEM 3", "ExEM 4",
                   "HEP", "Intermediate", "Meso 1", "Meso 2", "Meso 3", "Meso 4", 
                   "Meso 5", "PGC", "PS", "SYS", "VE/YS")

# Update the Seurat object with the cluster information
avgexp$orig.ident <- factor(avgexp$orig.ident, levels = cluster_names)
Idents(avgexp) <- "orig.ident"

# Create a heatmap using the top 50 genes
heatmap_obj <- DoHeatmap(avgexp, features = top50_genes_per_cluster$gene) + NoLegend()
heatmap_data <- heatmap_obj$data %>% drop_na()

# Reshape the heatmap data into a matrix format
expression_matrix <- heatmap_data %>%
  dplyr::select(Feature, Cell, Expression) %>%
  tidyr::spread(key = Cell, value = Expression) %>%
  column_to_rownames("Feature") %>%
  as.matrix()

# Reorder rows based on gene list and cell type order
expression_matrix <- expression_matrix[gene_list, ordered_cell_names]

# Define color palette for heatmap
start_color <- "#DDF1C0"  # Light greenish color for low expression
mid_color <- "#F7F7F7"    # Neutral color for medium expression
end_color <- "#DD1C77"    # Deep purple color for high expression

# Create a color gradient for the heatmap
color_palette_below <- colorRampPalette(colors = c(start_color, mid_color))
color_palette_above <- colorRampPalette(colors = c(mid_color, end_color))
color_palette <- c(color_palette_below(125), color_palette_above(125))

# Define breaks for the color gradient
breaks <- c(seq(-2.5, 0, length.out = 126), seq(0, 2.5, length.out = 126)[-1])

# Define colors for each cell type
celltype_names <- c("Epi", "Am", "PGC", "PS", "Early.meso", "Meso 1", "Meso 2", 
                    "Meso 3", "Meso 4", "Meso 5", "ExEM 1", "ExEM 2", "ExEM 3", 
                    "ExEM 4", "Intermediate", "VE/YS", "SYS", "HEP")

colors_celltype <- c("#65CEF6", "#588198", "#6B6A6B", "#D5E7F7", "#DEE9B1", "#63B472",
                     "#2D563D", "#9FD8C8", "#497C76", "#7AAE93", "#F5C785", "#86382A",
                     "#E9AA44", "#CAA57D", "#F8BFAF", "#F3746C", "#BB4A94", "#C45338")
names(colors_celltype) <- celltype_names

annotation_colors <- list(celltype = colors_celltype)

# Update cell type annotations in the Seurat object
avgexp$celltype <- Idents(avgexp)
annotations <- avgexp@meta.data[, c("celltype"), drop = FALSE]

# Reverse the order of the rows in the matrix
expression_matrix <- expression_matrix[rev(row.names(expression_matrix)), ]

# Generate the heatmap
heatmap_result <- pheatmap(expression_matrix,
                           color = color_palette,
                           breaks = breaks,
                           cluster_rows = FALSE,
                           cluster_cols = FALSE,
                           show_rownames = FALSE,
                           show_colnames = TRUE,
                           annotation_col = annotations,
                           annotation_colors = annotation_colors)

# Save the heatmap as a PDF
ggsave("pHeatmap_embryo_model.pdf", plot = heatmap_result, width = 14, height = 14)