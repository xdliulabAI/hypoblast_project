# Hypoblast Project - Code Repository

This repository contains the code and scripts used in the analysis presented in the manuscript. Below is a brief overview of the content of the repository, along with descriptions of the main analysis steps.

## Directory Descriptions

### Single-cell RNA sequencing data preprocessing
- **1_cellranger.sh**: Shell script for running Cell Ranger on single-cell RNA-seq data.
- **2_generate_sobj_before_integrate.R**: R script for generating Seurat objects before dataset integration.

### Dimensionality reduction clustering and dataset integration
- **run_prod_E1E2_integrate_harmony.sh**: Shell script for running the integration process.
- **prod_integrate_E1E2.R**: R script for integrating datasets E1 and E2 using Harmony.

### Gene expression analysis and visualization
- **Dotplot.R**: R script for generating dot plots of gene expression.
- **Heatmap.R**: R script for generating heatmaps of gene expression.
- **genelist.csv**: CSV file listing the genes used in the analysis.
- **DEG_embryo_model.csv**: CSV file containing differential expression gene data for the embryo model.

### Pseudotime Analysis and Differentiation Potential Estimation
- **1_CytoTRACE2.R**: R script for running CytoTRACE2 for differentiation potential estimation.
- **2_Monocle3.R**: R script for performing pseudotime analysis using Monocle3.

### Regulatory network analysis with SCENIC
- **1_prepare_data.R**: R script for preparing data for SCENIC analysis.
- **2_prepare_loom.R**: R script for preparing loom files for SCENIC.
- **3.0_run_GRNboost.sh**: Shell script for running the GRNBoost network inference.
- **3.1_arboreto_with_multiprocessing.py**: Python script for running Arboreto with multiprocessing.
- **3.2_generate_regulon.R**: R script for generating regulons.
- **4_heatmap.R**: R script for generating heatmaps of regulon activity.
- **filtered_regulonActivity_filtered.csv**: CSV file of filtered regulon activity data.

### Principal Component and Correlation Analysis
- **PCA.R**: R script for performing Principal Component Analysis (PCA).
- **correlation_Heatmap.R**: R script for generating correlation heatmaps.

## Methods
The methods related to the analysis performed in this project can be found in the manuscript Methods section. 

This document provides detailed descriptions of the data preprocessing and analytical approaches used in the study.

Feel free to explore the code and data provided in this repository. If you have any questions or need further information, please refer to the manuscript or contact the corresponding author.
