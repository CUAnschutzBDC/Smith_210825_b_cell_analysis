library(Seurat)
library(tidyverse)
library(cowplot)
library(here)
library(scAnalysisR)
library(viridis)
library(clustifyr)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

sample <- "healthy_bcells_all"

all_ref_dir <-
  "/Users/wellskr/Documents/Analysis/references/single_cell_references"

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


#-------------------------------------------------------------------------------

#########################
# human cell atlas pbmc #
#########################

# Information for cell mapping
ref_dir <- file.path(all_ref_dir, "pbmc/hca_pbmc")

ref_mat <- read.csv(file.path(ref_dir, "clustifyr_reference.csv"),
                    header = TRUE, row.names = 1)


DefaultAssay(seurat_data) <- "RNA"

cluster_res <- name_clusters(seurat_data, ref_mat,
                             save_dir = save_dir,
                             save_name = "pbmc_celltype", ADT = FALSE,
                             assay = "RNA",
                             nfeatures = 2000, clusters = "RNA_cluster",
                             plot_type = "rna.umap")

seurat_data <- cluster_res$object

seurat_res <- cluster_res$RNA

plotDimRed(seurat_data, col_by = "RNA_pbmc_celltype",
           plot_type = "rna.umap")

write.csv(seurat_res, file.path(save_dir,
                                "files/celltype_hca_pbmc.csv"))

# On combined clusters

cluster_res <- name_clusters(seurat_data, ref_mat,
                             save_dir = save_dir,
                             save_name = "pbmc_combined_celltype", ADT = FALSE,
                             assay = "RNA",
                             nfeatures = 2000, clusters = "combined_cluster",
                             plot_type = "wnn.umap")

seurat_data <- cluster_res$object

seurat_res <- cluster_res$RNA

plotDimRed(seurat_data, col_by = "RNA_pbmc_celltype",
           plot_type = "rna.umap")
plotDimRed(seurat_data, col_by = "RNA_pbmc_combined_celltype",
           plot_type = "rna.umap")
plotDimRed(seurat_data, col_by = "RNA_pbmc_celltype",
           plot_type = "wnn.umap")
plotDimRed(seurat_data, col_by = "RNA_pbmc_combined_celltype",
           plot_type = "wnn.umap")

write.csv(seurat_res, file.path(save_dir,
                                "files/celltype_combined_hca_pbmc.csv"))

#-------------------------------------------------------------------------------

#############
# bulk pbmc #
#############

# Information for cell mapping
ref_dir <- file.path(all_ref_dir, "pbmc/bulk_RNAseq")

ref_mat <- read.table(file.path(ref_dir,
                              "gene_id_GSE118165_RNA_gene_abundance.txt"),
                    header = TRUE, row.names = 1,
                    sep = " ")


DefaultAssay(seurat_data) <- "RNA"

cluster_res <- name_clusters(seurat_data, ref_mat,
                             save_dir = save_dir,
                             save_name = "celltype", ADT = FALSE,
                             assay = "RNA",
                             nfeatures = 2000, clusters = "RNA_cluster",
                             plot_type = "rna.umap")

seurat_data <- cluster_res$object

seurat_res <- cluster_res$RNA

plotDimRed(seurat_data, col_by = "RNA_celltype",
           plot_type = "rna.umap")

write.csv(seurat_res, file.path(save_dir,
                                "files/celltype_bulk_pbmc.csv"))

# On combined clusters


cluster_res <- name_clusters(seurat_data, ref_mat,
                             save_dir = save_dir,
                             save_name = "combined_celltype", ADT = FALSE,
                             assay = "RNA",
                             nfeatures = 2000, clusters = "combined_cluster",
                             plot_type = "wnn.umap")

seurat_data <- cluster_res$object

seurat_res <- cluster_res$RNA

plotDimRed(seurat_data, col_by = "RNA_celltype",
           plot_type = "rna.umap")
plotDimRed(seurat_data, col_by = "RNA_combined_celltype",
           plot_type = "rna.umap")
plotDimRed(seurat_data, col_by = "RNA_celltype",
           plot_type = "wnn.umap")
plotDimRed(seurat_data, col_by = "RNA_combined_celltype",
           plot_type = "wnn.umap")


write.csv(seurat_res, file.path(save_dir,
                                "files/celltype_combined_bulk_pbmc.csv"))

saveRDS(seurat_data, file.path(save_dir, "rda_obj", "seurat_processed.rds"))
