#' @title Single-cell RNA-seq Analysis Pipeline for Multiorgan Development
#' @description Comprehensive workflow including QC, Integration, Annotation, Trajectory, and Cell-Chat.
#' @author [Jie Xing/Dlk1-Dio3]
#' @date 2026

# 1. Environment Setup ----------------------------------------------------
library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(scCATCH)
library(CellChat)
library(monocle3)
library(CytoTRACE2)
library(ggpubr)
library(RColorBrewer)

# Set project directory - Modify this path to your local data folder
PROJECT_DIR <- "./data/feature_matrix"
OUTPUT_DIR <- "./results"
if(!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR)

# 2. Data Loading and QC -------------------------------------------------
# Define sample mapping
sample_dirs <- c(
  'DGHP0219/filtered_feature_bc_matrix', 'DGMP0409/filtered_feature_bc_matrix',
  'DGP2P321/filtered_feature_bc_matrix', 'DGW1P321/filtered_feature_bc_matrix',
  'DGHH0219/filtered_feature_bc_matrix', 'DGMH0409/filtered_feature_bc_matrix',
  'DGP2H321/filtered_feature_bc_matrix', 'DGW1H321/filtered_feature_bc_matrix',
  'DGHL0219/filtered_feature_bc_matrix', 'DGML0409/filtered_feature_bc_matrix',
  'DGP2L321/filtered_feature_bc_matrix', 'DGW1L321/filtered_feature_bc_matrix'
)
sample_names <- c('PHomo', 'PMK', 'PPK', 'PWT','HHomo', 'HMK', 'HPK', 'HWT','LHomo', 'LMK', 'LPK', 'LWT')
names(sample_dirs) <- sample_names

# Batch reading
seurat_list <- lapply(seq_along(sample_dirs), function(i) {
  counts <- Read10X(data.dir = file.path(PROJECT_DIR, sample_dirs[i]))
  obj <- CreateSeuratObject(counts, project = names(sample_dirs)[i], min.cells = 3, min.features = 200)
  return(obj)
})

# Merge datasets
seurat_merged <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = sample_names)

# Metadata Extraction (Vectorized)
seurat_merged$Tissue <- substr(colnames(seurat_merged), 1, 1) # Extract first character
seurat_merged$Type <- gsub(".*_", "", seurat_merged$orig.ident) 
seurat_merged$Condition <- ifelse(grepl("Homo|MK", seurat_merged$orig.ident), "Lethal", "Vital")

# QC Metrics Calculation
seurat_merged[["percent.mt"]] <- PercentageFeatureSet(seurat_merged, pattern = "^mt-")
seurat_merged[["percent.rp"]] <- PercentageFeatureSet(seurat_merged, pattern = "^Rp[sl]")
seurat_merged[["percent.hb"]] <- PercentageFeatureSet(seurat_merged, pattern = "^Hb[^(P)l]")

# Filtering
seurat_qc <- subset(seurat_merged, 
                    subset = nFeature_RNA > 200 & 
                      nFeature_RNA < 7000 & 
                      percent.mt < 30 & 
                      percent.rp < 18)

# 3. Normalization and Harmony Integration -------------------------------
seurat_integrated <- NormalizeData(seurat_qc) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(verbose = FALSE)

# Integration using Harmony
seurat_integrated <- RunHarmony(seurat_integrated, group.by.vars = "orig.ident")

# Dimensionality Reduction and Clustering
seurat_integrated <- RunUMAP(seurat_integrated, reduction = "harmony", dims = 1:45)
seurat_integrated <- FindNeighbors(seurat_integrated, reduction = "harmony", dims = 1:45) %>%
  FindClusters(resolution = 0.5)

# 4. Cell Type Annotation -------------------------------------------------
# Marker Dictionary
markers_ref <- list(
  B_Cell      = c("Cd19", "Cd79b", "Cd79a"),
  Endothelial = c("Cldn5", "Pecam1", "Vegf1"),
  Hepatocyte  = c("Hnf4a", "Alb"),
  Cardiomyocyte = c("Actc1", "Tnnc1", "Cx43"),
  Trophoblast = c("Prap1", "Lcn2", "Krt7"),
  Erythroid   = c("Hbb-bs", "Hba-a2", "Hba-a1")
)

# Manual Cluster Assignment (Based on marker analysis)
new_ident <- c(
  '0' = "Stromal", '1' = "Erythroid", '2' = "Erythroid", '3' = "Endothelial",
  '10' = "Cardiomyocyte", '12' = "Hepatocyte", '21' = "Spongiotrophoblast"
  # ... Add remaining cluster assignments here
)
seurat_integrated$CellType <- recode(as.character(Idents(seurat_integrated)), !!!new_ident)
Idents(seurat_integrated) <- "CellType"

# 5. Cell-Cell Communication Analysis (CellChat) -------------------------
run_cellchat_analysis <- function(seurat_obj, condition_name) {
  data.input <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")
  meta <- seurat_obj@meta.data
  
  chat <- createCellChat(object = data.input, meta = meta, group.by = "CellType")
  chat@DB <- CellChatDB.mouse
  
  chat <- subsetData(chat) %>%
    identifyOverExpressedGenes() %>%
    identifyOverExpressedInteractions() %>%
    computeCommunProb(population.size = TRUE) %>%
    filterCommunication(min.cells = 10) %>%
    computeCommunProbPathway() %>%
    aggregateNet()
  
  return(chat)
}

# Example: Analysis of Condition Group
# chat_lethal <- run_cellchat_analysis(subset(seurat_integrated, Condition == "Lethal"), "Lethal")

# 6. Trajectory Analysis (Monocle 3) --------------------------------------
cds <- new_cell_data_set(
  expression_data = GetAssayData(seurat_integrated, assay = 'RNA', slot = 'counts'),
  cell_metadata = seurat_integrated@meta.data,
  gene_metadata = data.frame(gene_short_name = rownames(seurat_integrated), row.names = rownames(seurat_integrated))
)

cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
cds <- order_cells(cds) # Interactive or root_cells based

# Plotting Trajectory
pdf(file.path(OUTPUT_DIR, "trajectory_plot.pdf"), width = 8, height = 7)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE)
dev.off()

# 7. Statistics and Visualization -----------------------------------------
# Cell Proportion Analysis
calc_proportion <- function(obj, group_col) {
  prop.table(table(Idents(obj), obj@meta.data[[group_col]]), margin = 2) %>%
    as.data.frame() -> df
  colnames(df) <- c("CellType", "Group", "Proportion")
  return(df)
}

prop_data <- calc_proportion(seurat_integrated, "Condition")

# Barplot for proportions
ggplot(prop_data, aes(x = Group, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  theme_classic() +
  scale_fill_manual(values = brewer.pal(12, "Paired")) +
  labs(title = "Cell Type Proportion by Condition", y = "Percentage (%)")

# 8. Save Data ------------------------------------------------------------
saveRDS(seurat_integrated, file = file.path(OUTPUT_DIR, "processed_seurat_object.rds"))