library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(Matrix)
library(ggplot2)
library(patchwork)
library(dplyr)
library(RColorBrewer)

# Load the Seurat object with CytoTRACE results
sobj <- readRDS("sobj_cytotrace2_D0_D4_D8_D20_P15_mnn.rds")

# Convert Seurat object to a Monocle3 cell_data_set object
cds <- as.cell_data_set(sobj)

# Cluster cells with a low resolution
cds <- cluster_cells(cds, resolution = 1e-3)

# Plot the cells colored by cluster and partition without showing the trajectory graph
p1 <- plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
wrap_plots(p1, p2)

# Learn the trajectory graph
cds <- learn_graph(cds, use_partition = TRUE, verbose = FALSE, close_loop = FALSE, learn_graph_control = list(ncenter = 400))

# Define a custom color palette using ColorBrewer's RdYlBu palette
color <- rev(brewer.pal(11, "RdYlBu"))

# Plot the cells colored by CytoTRACE2_Relative with a gradient color scale
plot <- plot_cells(cds,
                   color_cells_by = "CytoTRACE2_Relative",
                   group_cells_by = "cluster",
                   label_cell_groups = FALSE,
                   label_groups_by_cluster = FALSE,
                   label_leaves = FALSE,
                   label_branch_points = FALSE,
                   label_roots = FALSE
) + 
  scale_color_gradientn(
    colors = color,
    limits = c(0, 1),
    breaks = c(0, 1),
    labels = c("0.0 (More diff.)", "1.0 (Less diff.)")
  ) + coord_fixed()

# Save the plot as a PDF
ggsave(paste0("CytoSCAPE2_Monocle3_D0_D4_D8_D20_P15.pdf"), plot, width = 7, height = 7, units = "in")
