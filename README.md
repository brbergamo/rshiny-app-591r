# RShiny App for RNA-Seq Data Analysis

This Shiny app provides an interactive interface to explore sample metadata, gene counts, and differential expression results. Users can visualize, filter, and analyze datasets to support downstream analysis.

## How to Run the App
- **Locally in R**: Run the following command in R:
  ```r
  library(shiny)
  runApp("path/to/app")
  ```
  Replace `"path/to/app"` with the actual path to your Shiny app directory.

## Features

### 1. Sample Information Exploration
Upload and inspect sample metadata to understand dataset distribution.

**Inputs:**
- Sample information matrix (CSV format).

**Shiny Functionalities:**
- **Summary Tab:** Overview of dataset structure, including sample size, metadata types, and key statistics.
- **Data Table Tab:** Sortable table displaying raw sample data.
- **Visualization Tab:** Histograms, density plots, and violin plots with options for variable selection and grouping.

### 2. Counts Matrix Exploration
Assess gene count filtering strategies and data structure.

**Inputs:**
- Normalized counts matrix (CSV format).
- Filtering controls:
  - Slider for minimum variance percentile.
  - Slider for minimum non-zero sample count.

**Shiny Functionalities:**
- **Filtering Summary Tab:** Displays total samples, genes, and filtering impact.
- **Diagnostic Plots Tab:** Scatter plots for variance, count distribution, and zero counts.
- **Clustered Heatmap Tab:** Heatmap of filtered counts with optional log transformation.
- **PCA Tab:** Principal component analysis (PCA) scatter plots and variance-explained summary.

### 3. Differential Expression Analysis
Load and explore differential gene expression results.

**Inputs:**
- Differential expression results (CSV format) from DESeq2, limma, or edgeR.

**Shiny Functionalities:**
- **Results Table Tab:** Sortable table with search functionality.
- **Visualization Tab:** Plots for differential expression analysis insights.

## How to Use
1. Upload CSV files for sample metadata, gene counts, and differential expression results.
2. Navigate tabs to explore summaries, visualizations, and filtering options.
3. Adjust parameters using interactive controls.
4. Download or export results for further analysis.

This app simplifies transcriptomic data exploration, helping users assess sample characteristics, filter genes, and analyze differential expression results efficiently.


