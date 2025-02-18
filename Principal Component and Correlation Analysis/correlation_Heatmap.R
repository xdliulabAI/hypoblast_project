# Load required libraries for data manipulation, matrix operations, and visualization
library(tidyverse)
library(rhdf5)
library(Matrix)
library(writexl)
library(Seurat)
library(loupeR)
library(ggplot2)
library(stats)
library(pheatmap)

# Load the precomputed combined matrix
combined_matrix <- readRDS("combined_matrix.rds")

# Compute the correlation matrix for the filtered matrix
correlation_matrix <- cor(combined_matrix)

# Set up the PDF output for the correlation heatmap
pdf("Correlation_Heatmap_filtered_matrix.pdf", width = 10, height = 8)

# Generate the correlation heatmap with annotations
pheatmap(correlation_matrix, 
         clustering_distance_rows = "euclidean",     # Use Euclidean distance for row clustering
         clustering_distance_cols = "euclidean",     # Use Euclidean distance for column clustering
         color = colorRampPalette(c("#003f5c", "#2f4b7c", "#bc5090", "#ef5675", "#ffa600"))(100),  # Custom color palette
         main = "Correlation Heatmap of Filtered Matrix",  # Title of the heatmap
         display_numbers = TRUE,                     # Display the correlation values on the heatmap
         number_color = "black")                     # Set the color of the displayed numbers to black

# Close the PDF device
dev.off()
