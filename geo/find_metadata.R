library(Seurat)
library(tidyverse)
library(scAnalysisR)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

sample <- "healthy_bcells_all"

# Set directories
base_dir <- here()

base_dir_proj <- file.path(base_dir, "results", sample)
save_dir <- file.path(base_dir_proj, "R_analysis")

# Read in the data
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed.rds"))

# Pull out metadata
full_metadata <- seurat_data[[]] %>%
  dplyr::select(orig.ident, nCount_RNA, nFeature_RNA, nCount_ADT,
                nFeature_ADT, nCount_SCT, nFeature_SCT, nCount_CLR_ADT,
                nFeature_CLR_ADT, percent.mt, Doublet_finder,
                RNA_celltype) %>%
  dplyr::rename(all_cells_celltype = RNA_celltype)

sample <- "healthy_bcells_all_subset"

# Set directories
base_dir <- here()

base_dir_proj <- file.path(base_dir, "results", sample)
save_dir <- file.path(base_dir_proj, "R_analysis")

# Read in the data
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed.rds"))


subset_metadata <- seurat_data[[]] %>%
  dplyr::select(RNA_cluster,
                ADT_cluster, combined_cluster,
                RNA_celltype, RNA_combined_celltype,
                IgM_cutoff, IgD_cutoff, CD21_cutoff, CD27_cutoff, pass_filters,
                CD21_low_cutoff, CD21_high_cutoff, cd21_level, new_cluster) %>%
  dplyr::rename(final_cluster = new_cluster)

full_metadata$barcode <- rownames(full_metadata)
subset_metadata$barcode <-rownames(subset_metadata)

subset_metadata$in_final_obj <- "TRUE"

all_metadata <- merge(full_metadata, subset_metadata, by = "barcode",
                      all.x = TRUE, all.y = TRUE)

all_metadata <- all_metadata %>%
  dplyr::mutate(in_final_obj = case_when(in_final_obj == "TRUE" ~
                                      "TRUE",
                                      is.na(in_final_obj) ~ "FALSE"))

umap_data <- Embeddings(seurat_data, reduction = "rna.umap") %>% data.frame()

umap_data$barcode <- rownames(umap_data)

all_metadata <- merge(all_metadata, umap_data, by = "barcode",
                      all.x = TRUE, all.y = TRUE)

rownames(all_metadata) <- all_metadata$barcode

all_metadata$barcode <- NULL


write.csv(all_metadata, here("geo/final_metadata.csv"))
