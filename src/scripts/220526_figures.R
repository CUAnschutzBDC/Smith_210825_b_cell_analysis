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
fig_dir <- file.path(base_dir, "results", "220526_figure")

# Read in the data
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed.rds"))

# Figures ----------------------------------------------------------------------

## Umap of all clusters --------------------------------------------------------
all_cluster_colors <- c("6.0" = "#924bdb",
                        "0" = "#69ba3d",
                        "6.1" = "#9a43a4",
                        "1" = "#bf9b31",
                        "2" = "#6477ce",
                        "3" = "#d15131",
                        "4" = "#4c9e8e",
                        "5" = "#cc4570",
                        "7" = "#648d4f",
                        "6.2" = "#985978",
                        "8" = "#a06846",
                        "9" = "#7184b0")

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

pdf(file.path(fig_dir, "umap_similar.pdf"), width = 6, height = 6)
print(plotDimRed(seurat_data, col_by = "new_cluster", color = all_cluster_colors,
           plot_type = "rna.umap")[[1]])

dev.off()

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
                                                  "Bulk_B_8", "Bulk_B_9"))


all_celltype_cluster_colors <- all_cluster_colors
all_celltype_cluster <- unique(seurat_data$celltype_new_cluster)
all_celltype_cluster <- all_celltype_cluster[order(match(gsub(".*_", "",
                                                              all_celltype_cluster),
                                                         names(all_cluster_colors)))]
names(all_celltype_cluster_colors) <- all_celltype_cluster


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
rna_heatmap <- plot_heatmap(seurat_data, top10_rna$gene,
                            "celltype_new_cluster",
                            average_expression = FALSE,
                            colors = all_celltype_cluster_colors)

# Make a heatmap using the average values from each cluster
rna_heatmap_ave <- plot_heatmap(seurat_data, top10_rna$gene,
                                "celltype_new_cluster",
                                average_expression = TRUE,
                                colors = all_celltype_cluster_colors)[[4]]


pdf(file.path(fig_dir, "heatmap_sc_9_clusters_similar.pdf"),
    width = 6, height = 12)

grid.draw(rectGrob(gp=gpar(fill="white", lwd=0)))
# Draw the heatmap
grid.draw(rna_heatmap)

dev.off()

pdf(file.path(fig_dir, "heatmap_ave_9_clusters_similar.pdf"),
    width = 6, height = 12)

grid.draw(rectGrob(gp=gpar(fill="white", lwd=0)))
# Draw the heatmap
grid.draw(rna_heatmap_ave)

dev.off()

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

### Based on DE of all 12 clusters
# Get marker genes that are saved
marker_genes_rna <- read.csv(file.path(save_dir,
                                       "files/DE/RNA_markers_celltype_new_cluster.csv"),
                             row.names = 1)

# Pick out the top 10 makers of each cluster
top10_rna <- marker_genes_rna %>% 
  dplyr::group_by(cluster) %>% 
  dplyr::top_n(n = 10, wt = avg_log2FC) %>%
  dplyr::arrange(cluster)

# Make a heatmap with these genes
rna_heatmap <- plot_heatmap(seurat_data, top10_rna$gene, "celltype_new_cluster",
                            average_expression = FALSE,
                            colors = all_celltype_cluster_colors)

# Make a heatmap using the average values from each cluster
rna_heatmap_ave <- plot_heatmap(seurat_data, top10_rna$gene, "celltype_new_cluster",
                                average_expression = TRUE,
                                colors = all_celltype_cluster_colors)


pdf(file.path(fig_dir, "heatmap_sc_all_clusters_similar.pdf"),
    width = 6, height = 12)

grid.draw(rectGrob(gp=gpar(fill="white", lwd=0)))
# Draw the heatmap
grid.draw(rna_heatmap)

dev.off()

pdf(file.path(fig_dir, "heatmap_ave_all_clusters_similar.pdf"),
    width = 6, height = 12)

grid.draw(rectGrob(gp=gpar(fill="white", lwd=0)))
# Draw the heatmap
grid.draw(rna_heatmap_ave)

dev.off()

# Make a heatmap with these genes
rna_heatmap <- plot_heatmap(seurat_data, top10_rna$gene, "celltype_new_cluster",
                            average_expression = FALSE,
                            colors = celltype_cluster_colors_two)

# Make a heatmap using the average values from each cluster
rna_heatmap_ave <- plot_heatmap(seurat_data, top10_rna$gene, "celltype_new_cluster",
                                average_expression = TRUE,
                                colors = celltype_cluster_colors_two)


pdf(file.path(fig_dir, "heatmap_sc_all_clusters_distinct.pdf"),
    width = 6, height = 12)

grid.draw(rectGrob(gp=gpar(fill="white", lwd=0)))
# Draw the heatmap
grid.draw(rna_heatmap)

dev.off()

pdf(file.path(fig_dir, "heatmap_ave_all_clusters_disctinct.pdf"),
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
                            colors = all_celltype_cluster_colors,
                            cluster_rows = TRUE)

# Make a heatmap using the average values from each cluster
rna_heatmap_ave <- plot_heatmap(seurat_data, gene_list, "celltype_new_cluster",
                                average_expression = TRUE,
                                colors = all_celltype_cluster_colors,
                                cluster_rows = TRUE)


pdf(file.path(fig_dir, "heatmap_sc_gene_list_clustered_similar.pdf"),
    width = 6, height = 6)

grid.draw(rectGrob(gp=gpar(fill="white", lwd=0)))
# Draw the heatmap
grid.draw(rna_heatmap)

dev.off()

pdf(file.path(fig_dir, "heatmap_ave_gene_list_clustered_similar.pdf"),
    width = 6, height = 6)

grid.draw(rectGrob(gp=gpar(fill="white", lwd=0)))
# Draw the heatmap
grid.draw(rna_heatmap_ave)

dev.off()

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


### Not clustered --------------------------------------------------------------
# Make a heatmap with these genes
rna_heatmap <- plot_heatmap(seurat_data, gene_list, "celltype_new_cluster",
                            average_expression = FALSE,
                            colors = all_celltype_cluster_colors)

# Make a heatmap using the average values from each cluster
rna_heatmap_ave <- plot_heatmap(seurat_data, gene_list, "celltype_new_cluster",
                                average_expression = TRUE,
                                colors = all_celltype_cluster_colors)


pdf(file.path(fig_dir, "heatmap_sc_gene_list_ordered_similar.pdf"),
    width = 6, height = 6)

grid.draw(rectGrob(gp=gpar(fill="white", lwd=0)))
# Draw the heatmap
grid.draw(rna_heatmap)

dev.off()

pdf(file.path(fig_dir, "heatmap_ave_gene_list_ordered_similar.pdf"),
    width = 6, height = 6)

grid.draw(rectGrob(gp=gpar(fill="white", lwd=0)))
# Draw the heatmap
grid.draw(rna_heatmap_ave)

dev.off()

# Make a heatmap with these genes
rna_heatmap <- plot_heatmap(seurat_data, gene_list, "celltype_new_cluster",
                            average_expression = FALSE,
                            colors = celltype_cluster_colors_two)

# Make a heatmap using the average values from each cluster
rna_heatmap_ave <- plot_heatmap(seurat_data, gene_list, "celltype_new_cluster",
                                average_expression = TRUE,
                                colors = celltype_cluster_colors_two)


pdf(file.path(fig_dir, "heatmap_sc_gene_list_ordered_distinct.pdf"),
    width = 6, height = 6)

grid.draw(rectGrob(gp=gpar(fill="white", lwd=0)))
# Draw the heatmap
grid.draw(rna_heatmap)

dev.off()

pdf(file.path(fig_dir, "heatmap_ave_gene_list_ordered_distinct.pdf"),
    width = 6, height = 6)

grid.draw(rectGrob(gp=gpar(fill="white", lwd=0)))
# Draw the heatmap
grid.draw(rna_heatmap_ave)

dev.off()

## Volcano plots ---------------------------------------------------------------
Idents(seurat_data) <- "cd21_level"
cd21_markers <- FindMarkers(seurat_data, ident.1 = "CD21_pos",
                           ident.2 = "CD21_neg", min.pct = 0,
                           logfc.threshold = 0) 

cd21_markers <- cd21_markers %>%
  dplyr::mutate(log_pval = -log(p_val_adj),
               color = ifelse(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25,
                              "significant", "not_significant"))

color_list <- c("significant" = "#B3669E",
                "not_significant" = "#000000")


pdf(file.path(fig_dir, "cd21_volcano.pdf"),
    width = 6, height = 6)

print(ggplot(cd21_markers, aes(x = avg_log2FC, y = log_pval, color = color)) +
  geom_point() +
  scale_color_manual(values = color_list))

dev.off()

Idents(seurat_data) <- "new_cluster"

cluster_6_markers <- FindMarkers(seurat_data, ident.1 = "6.1",
                                 ident.2 = "6.2", min.pct = 0,
                                 logfc.threshold = 0) 

cluster_6_markers <- cluster_6_markers %>%
  dplyr::mutate(log_pval = -log(p_val_adj),
                color = ifelse(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25,
                               "significant", "not_significant"))

color_list <- c("significant" = "#B3669E",
                "not_significant" = "#000000")

pdf(file.path(fig_dir, "6.1_vs_6.2_volcano.pdf"),
    width = 6, height = 6)

print(ggplot(cluster_6_markers, aes(x = avg_log2FC, y = log_pval, color = color)) +
  geom_point() +
  scale_color_manual(values = color_list))

dev.off()