library(Seurat)
library(tidyverse)
library(cowplot)
library(openxlsx)
library(here)
library(scAnalysisR)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

sample <- "healthy_bcells_all_subset"

clusters <- "cluster_6_RNA_cluster"
combined_clusters <- "cluster_6_combined_cluster"

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
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_clust6.rds"))


# RNA cluster DE ---------------------------------------------------------------
seurat_data$cluster_6_RNA_cluster <- seurat_data$RNA_cluster

seurat_data$cluster_6_combined_cluster <- seurat_data$combined_cluster

marker_list <- find_write_markers(seurat_object = seurat_data,
                                  meta_col = clusters,
                                  pval = pval,
                                  logfc = logfc,
                                  assay = "RNA",
                                  save_dir = save_dir)

if(ADT){
  marker_list <- find_write_markers(seurat_object = seurat_data,
                                    meta_col = clusters,
                                    pval = pval,
                                    logfc = logfc,
                                    assay = "ADT",
                                    save_dir = save_dir)
}


marker_list <- find_write_markers(seurat_object = seurat_data,
                                  meta_col = combined_clusters,
                                  pval = pval,
                                  logfc = logfc,
                                  assay = "RNA",
                                  save_dir = save_dir)

if(ADT){
  marker_list <- find_write_markers(seurat_object = seurat_data,
                                    meta_col = combined_clusters,
                                    pval = pval,
                                    logfc = logfc,
                                    assay = "ADT",
                                    save_dir = save_dir)
}


saveRDS(seurat_data, file.path(save_dir, "rda_obj", "seurat_clust6.rds"))
