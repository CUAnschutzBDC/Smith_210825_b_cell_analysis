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

# Cell type DE -----------------------------------------------------------------

marker_list <- find_write_markers(seurat_object = seurat_data,
                                  meta_col = cell_types,
                                  pval = pval,
                                  logfc = logfc,
                                  assay = "RNA",
                                  save_dir = save_dir)

if(ADT){
  marker_list <- find_write_markers(seurat_object = seurat_data,
                                    meta_col = cell_types,
                                    pval = pval,
                                    logfc = logfc,
                                    assay = "ADT",
                                    save_dir = save_dir)
}

# RNA cluster DE ---------------------------------------------------------------

seurat_data$celltype_cluster <- paste0(seurat_data[[cell_types]][[1]], "_",
                                       seurat_data[[clusters]][[1]])

seurat_data$celltype_cluster <- factor(seurat_data$celltype_cluster)

marker_list <- find_write_markers(seurat_object = seurat_data,
                                  meta_col = "celltype_cluster",
                                  pval = pval,
                                  logfc = logfc,
                                  assay = "RNA",
                                  save_dir = save_dir)

if(ADT){
  marker_list <- find_write_markers(seurat_object = seurat_data,
                                    meta_col = "celltype_cluster",
                                    pval = pval,
                                    logfc = logfc,
                                    assay = "ADT",
                                    save_dir = save_dir)
}

# DE IgMlo IgD+ CD21+ vs CD27- IgMlo IgD+ CD21- --------------------------------
seurat_data$cd21_level <- "none"

seurat_data$cd21_level[seurat_data$IgM_cutoff &
                         seurat_data$IgD_cutoff & 
                         seurat_data$CD21_high_cutoff] <- "CD21_pos"

seurat_data$cd21_level[seurat_data$IgM_cutoff &
                         seurat_data$IgD_cutoff & 
                         seurat_data$CD21_low_cutoff] <- "CD21_neg"

table(seurat_data$cd21_level)

Idents(seurat_data) <- "cd21_level"
marker_list <- FindMarkers(seurat_data, ident.1 = "CD21_pos",
                           ident.2 = "CD21_neg") 

save_list <- marker_list %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  tibble::rownames_to_column("gene_name")

plotDimRed(seurat_data, col_by = "cd21_level", plot_type = "rna.umap",
           color = c("none" = "light grey",
                     "CD21_neg" = "blue",
                     "CD21_pos" = "red"))

write.csv(marker_list, file.path(save_dir, "files/DE/cd21_de.csv"))

cd21_wb <- openxlsx::createWorkbook()

openxlsx::addWorksheet(wb = cd21_wb, sheetName = "cd21_de")
openxlsx::writeData(wb = cd21_wb, sheet = "cd21_de", x = save_list)

openxlsx::saveWorkbook(wb = cd21_wb, file = file.path(save_dir,
                                                      "files/DE/cd21_de.xlsx"),
                       overwrite = TRUE)

saveRDS(seurat_data, file.path(save_dir, "rda_obj", "seurat_processed.rds"))
