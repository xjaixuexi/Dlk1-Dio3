# Single-cell Analysis of Organ Development in Mouse Models

This repository contains the R code used for the data analysis of the paper: 
*Allele-specific CRISPR perturbation of the imprinted Dlk1–Dio3 domain reveals regulation of BMP–NOTCH–VEGF signaling in embryonic organogenesis*

## Workflow Overview
1. **Preprocessing**: Quality control and filtering of 10X Genomics data.
2. **Integration**: Batch effect correction using Harmony across 12 samples.
3. **Annotation**: Cell type identification using marker genes and scCATCH.
4. **Trajectory**: Pseudotime analysis using Monocle 3.
5. **Communication**: Intercellular signaling analysis using CellChat.

## Requirements
- R >= 4.4
- Seurat v5
- Harmony
- Monocle3
- CellChat

## Usage
1. Clone this repository.
2. Place your `filtered_feature_bc_matrix` folders in the `./data/` directory.
3. Run `analysis_main.R`.
