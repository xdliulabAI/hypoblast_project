# Load required libraries
library(tidyverse)
library(rhdf5)
library(Matrix)
library(writexl)
library(Seurat)
library(loupeR)
library(ggplot2)
library(stats)
library(pheatmap)
library(ggrepel)

combined_matrix <- readRDS("combined_matrix.rds")

# Perform Principal Component Analysis (PCA) on the filtered matrix
pca <- prcomp(t(combined_matrix), center = TRUE, scale. = TRUE)

# Calculate the percentage of variance explained by each principal component
percent_var <- pca$sdev^2 / sum(pca$sdev^2) * 100

# Prepare PCA results for visualization
pca_data <- data.frame(Sample = rownames(pca$x),
                       PC1 = pca$x[, 1],
                       PC2 = pca$x[, 2])

# Reorder PCA data to ensure "Hypo_Hypo" is the last row (if present)
pca_data <- pca_data %>%
  arrange(Sample == "Hypo_Hypo")

# Extract the dataset name from the sample names
pca_data$Dataset <- sapply(strsplit(pca_data$Sample, "_"), `[`, 1)

# Define an elegant color palette for the PCA plot
elegant_palette <- c("#4c6ef5", "#BD8537", "#5a4d9b", "#2e86de", "#3c6e71", "#2f4b7c", "#2a3f54", "#1b4965")

# Plot PCA with dataset-specific colors and labeled samples
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Dataset, label = Sample)) +
  geom_point(size = 3) +
  geom_text_repel(size = 3.5, max.overlaps = 100, force = 2) +  # Adjust text label placement
  scale_color_manual(values = elegant_palette) +
  theme_minimal() +
  labs(title = "PCA of Filtered Matrix",
       x = paste0("Principal Component 1 (", round(percent_var[1], 2), "% variance)"),
       y = paste0("Principal Component 2 (", round(percent_var[2], 2), "% variance)")) +
  theme(legend.title = element_blank())

# Save the PCA plot with labels as a PDF
pdf("../plots/PCA_combined_matrix_with_variance.pdf", width = 16, height = 12)
print(pca_plot)
dev.off()

# Plot PCA without labels for a cleaner view
pca_plot_no_label <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Dataset)) +
  geom_point(size = 6) +
  scale_color_manual(values = elegant_palette) +
  theme_minimal() +
  labs(title = "PCA of Filtered Matrix",
       x = paste0("Principal Component 1 (", round(percent_var[1], 2), "% variance)"),
       y = paste0("Principal Component 2 (", round(percent_var[2], 2), "% variance)")) +
  theme(legend.title = element_blank())

# Save the PCA plot without labels as a PDF
pdf("../plots/PCA_combined_matrix_no_label_with_variance.pdf", width = 16, height = 12)
print(pca_plot_no_label)
dev.off()
