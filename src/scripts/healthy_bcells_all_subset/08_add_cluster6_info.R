library(Seurat)
library(tidyverse)
library(cowplot)
library(here)
library(scAnalysisR)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

sample <- "healthy_bcells_all_subset"

cell_types <- "RNA_celltype"
clusters <- "RNA_cluster"
pval <- 0.05
logfc <- 0.5

HTO <- FALSE
ADT <- FALSE

if(normalization_method == "SCT"){
  SCT <- TRUE
  seurat_assay <- "SCT"
} else {
  SCT <- FALSE
  seurat_assay <- "RNA"
}

# Set directories
base_dir <- here()

base_dir_proj <- file.path(base_dir, "results", sample)
save_dir <- file.path(base_dir_proj, "R_analysis")

# Read in the data
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed.rds"))
clust_6_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_clust6.rds"))

new_metadata <- clust_6_data[[]] %>%
  dplyr::select(dplyr::contains("cluster_6"))

seurat_data <- AddMetaData(seurat_data, metadata = new_metadata)

new_cluster_df <- seurat_data[[]] %>%
  dplyr::select("RNA_cluster", "cluster_6_combined_cluster") %>%
  dplyr::mutate("new_cluster" = ifelse(RNA_cluster == "6",
                                       paste0(RNA_cluster,  ".",
                                              cluster_6_combined_cluster),
                                       RNA_cluster))

seurat_data <- AddMetaData(seurat_data, metadata = new_cluster_df)

saveRDS(seurat_data, file.path(save_dir, "rda_obj", "seurat_processed.rds"))
