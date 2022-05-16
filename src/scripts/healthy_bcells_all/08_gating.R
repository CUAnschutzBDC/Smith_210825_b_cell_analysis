library(Seurat)
library(tidyverse)
library(cowplot)
library(here)
library(scAnalysisR)
library(MASS)
library(viridis)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

sample <- "healthy_bcells_all"

vars.to.regress <- NULL

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


ADT_data <- GetAssayData(seurat_data, assay = "ADT", slot = "data") %>%
  t %>%
  data.frame

nColors <- ncol(ADT_data)
protein_colors <- grDevices::colorRampPalette(
  RColorBrewer::brewer.pal(9, "Set1"))(nColors)
names(protein_colors) <- colnames(ADT_data)

# Find cutoff percents

# CD27 - from Zach use top 20%
cd27_cutoff <- quantile(ADT_data$CD27.1, probs = .8)

# IgM - from zach bottom 20%
igm_cutoff <- quantile(ADT_data$IgM, prob = 0.2)

# IgD cutoff - a little less than 1, at the plateau
igd_cutoff <- quantile(ADT_data$IgD, prob = 0.5)

# CD21 cutoff CD21 Bottom 10-15% Maybe the bump around 0.5?
cd21_cutoff <- quantile(ADT_data$CD21, prob = 0.1)

ADT_data_long <- ADT_data %>%
  tidyr::pivot_longer(cols = all_of(colnames(.)), names_to = "protein",
                                    values_to = "value")

# Plots of cutoffs
IgM_plot <- ADT_data_long %>%
  dplyr::filter(protein == "IgM") %>%
  ggplot2::ggplot(ggplot2::aes(x = value, fill = protein)) +
    ggplot2::geom_density() +
    ggplot2::scale_fill_manual(values = protein_colors) + 
    ggplot2::geom_vline(xintercept = igm_cutoff)



IgD_plot <- ADT_data_long %>%
  dplyr::filter(protein == "IgD") %>%
  ggplot2::ggplot(ggplot2::aes(x = value, fill = protein)) +
  ggplot2::geom_density() +
  ggplot2::scale_fill_manual(values = protein_colors) + 
  ggplot2::geom_vline(xintercept = igd_cutoff)


CD21_plot <- ADT_data_long %>%
  dplyr::filter(protein == "CD21") %>%
  ggplot2::ggplot(ggplot2::aes(x = value, fill = protein)) +
  ggplot2::geom_density() +
  ggplot2::scale_fill_manual(values = protein_colors) + 
  ggplot2::geom_vline(xintercept = cd21_cutoff)


CD27_plot <- ADT_data_long %>%
  dplyr::filter(protein == "CD27.1") %>%
  ggplot2::ggplot(ggplot2::aes(x = value, fill = protein)) +
  ggplot2::geom_density() +
  ggplot2::scale_fill_manual(values = protein_colors) + 
  ggplot2::geom_vline(xintercept = cd27_cutoff)

# Add in info about if it is/is not passing cutoffs
ADT_data$IgM_cutoff <- ADT_data$IgM < igm_cutoff
ADT_data$IgD_cutoff <- ADT_data$IgD > igd_cutoff
ADT_data$CD21_cutoff <- ADT_data$CD21 < cd21_cutoff
ADT_data$CD27_cutoff <- ADT_data$CD27.1 > cd27_cutoff

ADT_data$pass_filters <- ADT_data$IgM_cutoff + 
                         ADT_data$IgD_cutoff + 
                         ADT_data$CD21_cutoff

new_metadata <- ADT_data %>%
  dplyr::select(c("IgM_cutoff", "IgD_cutoff", "CD21_cutoff", "CD27_cutoff",
                  "pass_filters"))

seurat_data <- AddMetaData(seurat_data, new_metadata)

seurat_data$pass_filters <- factor(seurat_data$pass_filters)

plotDimRed(seurat_data, "IgM_cutoff", plot_type = "rna.umap")
plotDimRed(seurat_data, "IgD_cutoff", plot_type = "rna.umap")
plotDimRed(seurat_data, "CD21_cutoff", plot_type = "rna.umap")
plotDimRed(seurat_data, "pass_filters", plot_type = "rna.umap")

plotDimRed(seurat_data, "celltype_cluster", plot_type = "rna.umap",
           highlight_group = TRUE, group = 3, meta_data_col = "pass_filters")

plotDimRed(seurat_data, "pass_filters", plot_type = "rna.umap",
           highlight_group = TRUE, group = 2, meta_data_col = "pass_filters")

plotDimRed(seurat_data, "IgM", plot_type = "rna.umap", assay = "ADT")

saveRDS(seurat_data, file.path(save_dir, "rda_obj", "seurat_processed.rds"))
