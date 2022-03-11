library(Seurat)
library(tidyverse)
library(cowplot)
library(here)
library(scAnalysisR)
library(MASS)
library(viridis)

# Get density of points in 2 dimensions.
# @param x A numeric vector.
# @param y A numeric vector.
# @param n Create a square n by n grid to compute density.
# @return The density within each square.
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}


# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

sample <- "healthy_bcells_1"

vars.to.regress <- NULL

HTO <- FALSE
ADT <- TRUE
adt_PCA <- FALSE # Set to true to run PCA on ADT data (if there are enough
# proteins ex: > 100). Otherwise, set to FALSE and the
# normalized values will be used.

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

DefaultAssay(seurat_data) <- "ADT"
adt_list <- rownames(seurat_data)

rna_list <- c(
  "CR2", # CD21,
  "CD27",
  "CXCR5",
  "IGHD", # IgD
  "IGHM" # IgM
)


dsb_adt_vals <- t(GetAssayData(seurat_data, slot = "data", assay = "ADT"))

clr_adt_vals <- t(GetAssayData(seurat_data, slot = "data", assay = "CLR_ADT"))

rna_vals <- t(GetAssayData(seurat_data, slot = "data", assay = "RNA")) %>%
  data.frame %>%
  dplyr::select(all_of(rna_list))


dsb_adt_vals <- cbind(dsb_adt_vals, seurat_data[["RNA_cluster"]])

clr_adt_vals <- cbind(clr_adt_vals, seurat_data[["RNA_cluster"]])

rna_vals <- cbind(rna_vals, seurat_data[["RNA_cluster"]])

nColors <- length(unique(seurat_data$RNA_cluster))
cluster_colors <- grDevices::colorRampPalette(
  RColorBrewer::brewer.pal(9, "Set1"))(nColors)
names(cluster_colors) <- unique(seurat_data$RNA_cluster)


make_plots <- function(plotting_df, plotting_list, data_type = "ADT"){
  return_list <- list()
  if(data_type == "ADT"){
    plotting_df$density <- get_density(plotting_df$IgM,
                                       plotting_df$IgD,
                                       n = 100)
    
    # IgM vs IgD
    plot1 <- ggplot2::ggplot(data = plotting_df,
                             mapping = ggplot2::aes(x = IgM, y = IgD,
                                                    color = RNA_cluster)) +
      ggplot2::scale_color_manual(values = cluster_colors) +
      ggplot2::geom_point()
    
    plot2 <- ggplot2::ggplot(data = plotting_df,
                             mapping = ggplot2::aes(x = IgM, y = IgD,
                                                    color = density)) +
      scale_color_viridis() +
      ggplot2::geom_point()
    
    return_list$igm_igg <- plot1
    return_list$igm_igg_density <- plot2
  }
  
  plotting_df_long <- plotting_df %>%
    tidyr::pivot_longer(cols = all_of(plotting_list),
                        names_to = "protein",
                        values_to = "expression")
  
  plot3 <- ggplot2::ggplot(data = plotting_df_long,
                           mapping = ggplot2::aes(expression, color = protein)) + 
    ggplot2::geom_density() + 
    ggplot2::scale_color_brewer(palette = "Set1")
  
  return_list$density <- plot3

  return(return_list)  
}

dsb_plots <- make_plots(dsb_adt_vals, plotting_list = adt_list)

clr_plots <- make_plots(clr_adt_vals, plotting_list = adt_list)

rna_plots <- make_plots(rna_vals, rna_list, data_type = "RNA")


plot_list <- plotDimRed(seurat_data, col_by = rna_list,
                        plot_type = "rna.umap", assay = "RNA")

plot_list_2 <- featDistPlot(seurat_data, geneset = rna_list,
                            sep_by = "RNA_cluster",
                            combine = FALSE, 
                            assay = "RNA")
# 
# igm_density <- density(dsb_adt_vals$IgM)
# 
# density_df <- data.frame(x = igm_density$x, y = igm_density$y)
# 
# optimize(approxfun(density_df$x, density_df$y),
#          interval=c(1,20))$minimum
# 
# ggplot2::ggplot(density_df, aes(x = x, y = y)) +
#   ggplot2::geom_point() +
#   ggplot2::geom_vline(xintercept = 11.36248)
# 
# 
# # Cutoff for IgM is 11.36
# 
# 
# igd_density <- density(dsb_adt_vals$IgD)
# 
# density_df <- data.frame(x = igd_density$x, y = igd_density$y)
# 
# optimize(approxfun(density_df$x, density_df$y),
#          interval=c(0,10))$minimum
# 
# ggplot2::ggplot(density_df, aes(x = x, y = y)) +
#   ggplot2::geom_point() +
#   ggplot2::geom_vline(xintercept = 9.999944)
