library(Seurat)
library(tidyverse)
library(cowplot)
library(here)
library(scAnalysisR)
library(grid)

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
fig_dir <- file.path(base_dir, "results", "220602_figure")

# Read in the data
seurat_data <- readRDS(file.path(save_dir, "rda_obj",
                                 "seurat_processed_final.rds"))

# Figures ----------------------------------------------------------------------

## Umap of all clusters --------------------------------------------------------
cluster_colors_two <- c("0" = "#924bdb",
                        "1" = "#69ba3d",
                        "3" = "#9a43a4",
                        "2" = "#bf9b31",
                        "4" = "#6477ce",
                        "5" = "#d15131",
                        "6.0" = "#4c9e8e",
                        "6.1" = "#cc4570",
                        "6.2" = "#648d4f",
                        "7" = "#985978",
                        "8" = "#a06846")

pdf(file.path(fig_dir, "umap_distinct.pdf"), width = 6, height = 6)

print(plotDimRed(seurat_data, col_by = "new_cluster", color = cluster_colors_two,
                 plot_type = "rna.umap")[[1]])
dev.off()

## Heatmaps --------------------------------------------------------------------

seurat_data$celltype_new_cluster <- factor(seurat_data$celltype_new_cluster,
                                           levels = c("Mem_B_0", "Bulk_B_1",
                                                      "Mem_B_2", "Bulk_B_3",
                                                      "Bulk_B_4", "Mem_B_5",
                                                      "Mem_B_6.0", "Mem_B_6.1",
                                                      "Mem_B_6.2", "Mem_B_7",
                                                      "Bulk_B_8"))


celltype_cluster_colors_two <- cluster_colors_two
all_celltype_cluster <- unique(seurat_data$celltype_new_cluster)
all_celltype_cluster <- all_celltype_cluster[order(match(gsub(".*_", "",
                                                              all_celltype_cluster),
                                                         names(cluster_colors_two)))]
names(celltype_cluster_colors_two) <- all_celltype_cluster

### Based on DE of just the 9 main clusters
# Get marker genes that are saved
marker_genes_rna <- read.csv(file.path(save_dir,
                                       "files/DE/RNA_markers_celltype_cluster.csv"),
                             row.names = 1)

# Pick out the top 10 makers of each cluster
top10_rna <- marker_genes_rna %>% 
  dplyr::group_by(cluster) %>% 
  dplyr::top_n(n = 10, wt = avg_log2FC) %>%
  dplyr::arrange(cluster)

# Make a heatmap with these genes
rna_heatmap <- plot_heatmap(seurat_data, top10_rna$gene, "celltype_new_cluster",
                            average_expression = FALSE,
                            colors = celltype_cluster_colors_two)

# Make a heatmap using the average values from each cluster
rna_heatmap_ave <- plot_heatmap(seurat_data, top10_rna$gene, "celltype_new_cluster",
                                average_expression = TRUE,
                                colors = celltype_cluster_colors_two)


pdf(file.path(fig_dir, "heatmap_sc_9_clusters_distinct.pdf"),
    width = 6, height = 12)

grid.draw(rectGrob(gp=gpar(fill="white", lwd=0)))
# Draw the heatmap
grid.draw(rna_heatmap)

dev.off()

pdf(file.path(fig_dir, "heatmap_ave_9_clusters_disctinct.pdf"),
    width = 6, height = 12)

grid.draw(rectGrob(gp=gpar(fill="white", lwd=0)))
# Draw the heatmap
grid.draw(rna_heatmap_ave)

dev.off()

## Gene heatmaps ---------------------------------------------------------------
gene_list <- c("ITGAX", "TBX21", "CD1C", "FGR", "FCRL5", "CXCR5", "CR2",
               "HCK", "CXCR3", "CXCR4", "CCR4", "PRDM1", "XBP1", "IRF4",
               "CD38", "CD27", "TRAF5", "FCER2", "CD86", "CD80", "CD69",
               "CD72", "CD22", "NR4A1", "IGHD", "IGHM", "IGHG1", "IGHG2",
               "IL21R")


### Clustered ------------------------------------------------------------------

# Make a heatmap with these genes
rna_heatmap <- plot_heatmap(seurat_data, gene_list, "celltype_new_cluster",
                            average_expression = FALSE,
                            colors = celltype_cluster_colors_two,
                            cluster_rows = TRUE)

# Make a heatmap using the average values from each cluster
rna_heatmap_ave <- plot_heatmap(seurat_data, gene_list, "celltype_new_cluster",
                                average_expression = TRUE,
                                colors = celltype_cluster_colors_two,
                                cluster_rows = TRUE)


pdf(file.path(fig_dir, "heatmap_sc_gene_list_clustered_distinct.pdf"),
    width = 6, height = 6)

grid.draw(rectGrob(gp=gpar(fill="white", lwd=0)))
# Draw the heatmap
grid.draw(rna_heatmap)

dev.off()

pdf(file.path(fig_dir, "heatmap_ave_gene_list_clustered_distinct.pdf"),
    width = 6, height = 6)

grid.draw(rectGrob(gp=gpar(fill="white", lwd=0)))
# Draw the heatmap
grid.draw(rna_heatmap_ave)

dev.off()

## Quality plots ---------------------------------------------------------------
cluster_colors_two <- c("0" = "#924bdb",
                        "1" = "#69ba3d",
                        "3" = "#9a43a4",
                        "2" = "#bf9b31",
                        "4" = "#6477ce",
                        "5" = "#d15131",
                        "6.0" = "#4c9e8e",
                        "6.1" = "#cc4570",
                        "6.2" = "#648d4f",
                        "7" = "#985978",
                        "8" = "#a06846",
                        "9" = "#7184b0")

# Read in the data
seurat_data_all <- readRDS(file.path(save_dir, "rda_obj",
                                     "seurat_processed.rds"))

seurat_data_all$new_cluster <- factor(seurat_data_all$new_cluster,
                                           levels = c("0", "1",
                                                      "2", "3",
                                                      "4", "5",
                                                      "6.0", "6.1",
                                                      "6.2", "7",
                                                      "8", "9"))


violin_plot <- featDistPlot(seurat_data_all,
                            geneset = c("nFeature_RNA", "nCount_RNA","percent.mt"),
                            color = cluster_colors_two, sep_by = "new_cluster")

pdf(file.path(fig_dir, "cluster_9_quality.pdf"),
    width = 10, height = 10)

plot(violin_plot)

dev.off()
