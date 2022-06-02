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

cell_types <- "RNA_celltype"
clusters <- "new_cluster"
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

# RNA cluster DE ---------------------------------------------------------------

seurat_data$celltype_new_cluster <- paste0(seurat_data[[cell_types]][[1]], "_",
                                       seurat_data[[clusters]][[1]])

seurat_data$celltype_new_cluster <- factor(seurat_data$celltype_new_cluster)

saveRDS(seurat_data, file.path(save_dir, "rda_obj", "seurat_processed.rds"))


# Remove cluster 9
seurat_data <- subset(seurat_data, subset = new_cluster != 9)

marker_list <- find_write_markers(seurat_object = seurat_data,
                                  meta_col = "celltype_new_cluster",
                                  pval = pval,
                                  logfc = logfc,
                                  assay = "RNA",
                                  save_dir = save_dir)

if(ADT){
  marker_list <- find_write_markers(seurat_object = seurat_data,
                                    meta_col = "celltype_new_cluster",
                                    pval = pval,
                                    logfc = logfc,
                                    assay = "ADT",
                                    save_dir = save_dir)
}

saveRDS(seurat_data, file.path(save_dir, "rda_obj",
                               "seurat_processed_final.rds"))
