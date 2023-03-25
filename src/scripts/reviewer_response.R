library(Seurat)
library(tidyverse)
library(cowplot)
library(openxlsx)
library(here)
library(scAnalysisR)
library(clustree)

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

new_names <- c("0" = "Resting_memory",
               "1" = "Naive_1",
               "3" = "Naive_2",
               "2" = "Memory_IgE_IgG",
               "4" = "Naive_3",
               "5" = "Memory_IgA",
               "6.0" = "Early_memory",
               "6.1" = "BND2",
               "6.2" = "DN2",
               "7" = "Activated_memory",
               "8" = "Activated_naive")

seurat_data$RNA_celltype_final <- new_names[seurat_data$new_cluster]

new_colors <- c("Resting_memory" = "#924bdb", # Resting memory
                "Naive_1" = "#69ba3d", # Naive 1
                "Naive_2" = "#9a43a4", # Naive 2
                "Memory_IgE_IgG" = "#bf9b31", # Memory IgE/IgG1
                "Naive_3" = "#6477ce", # Naive 3
                "Memory_IgA" = "#d15131", # Memory IA
                "Early_memory" = "#4c9e8e", # Early Memory
                "BND2" = "#cc4570", #Bnd2
                "DN2" = "#648d4f", # DN2
                "Activated_memory" = "#985978", # Activated memory
                "Activated_naive" = "#a06846") # Activated naive


print(plotDimRed(seurat_data, col_by = "RNA_celltype_final", color = new_colors,
                 plot_type = "rna.umap")[[1]])


seurat_data$original_clusters <- seurat_data$RNA_cluster

ifelse(!dir.exists(here("results/reviewer_figs")),
       dir.create(here("results/reviewer_figs")),
       FALSE)

# Make original plots ----------------------------------------------------------
# PC number choice

RNA_plots <- plot_PCA(HTO = HTO, assay = seurat_assay,
                      sample_object = seurat_data,
                      jackstraw = TRUE)


pdf(here("results/reviewer_figs/pca_testing.pdf"))

print(RNA_plots[[6]])
print(RNA_plots[[7]])

dev.off()

# Resolution choice
seurat_data <- FindClusters(seurat_data, resolution = c(0.2, 0.5, 0.8, 1))

pdf(here("results/reviewer_figs/resolution_testing.pdf"))
clustree(seurat_data, prefix = "RNA_snn_res.")
dev.off()




# Test different numbers of PCs
# Starting pc # was 14, resolution was 0.8

make_heatmap <- function(seurat_data,
                         pc_val = 14, resolution = 0.8){
  name <- paste0("pc_", pc_val, " | res_", resolution)
  
  set.seed(0)
  
  # Remove previous clustering
  remove_cols <- colnames(seurat_data[[]])[grepl("res\\.[0-9]",
                                                 colnames(seurat_data[[]]))]
  
  for (i in remove_cols){
    seurat_data[[i]] <- NULL
  }
  
  
  umap_data <- group_cells(seurat_data, sample, save_dir, nPCs = pc_val,
                           resolution = resolution, assay = seurat_assay, 
                           HTO = HTO)
  
  seurat_data <- umap_data$object
  
  
  cM <- confusionMatrix(i = seurat_data$original_clusters,
                        j = seurat_data$RNA_cluster)
  
  cM <- cM / rowSums(cM)
  
  heatmap <- pheatmap::pheatmap(cM,
                                main = name,
                                cluster_rows = FALSE,
                                cluster_cols = FALSE)  
}


test_pcs <- c(10, 15, 20)

pc_heatmaps <- lapply(test_pcs, function(x){
  heatmap <- make_heatmap(seurat_data, pc_val = x)
  return(heatmap)
})

test_resolution <- c(0.2, 0.5, 1, 1.2)

resolution_heatmaps <- lapply(test_resolution, function(x){
  heatmap <- make_heatmap(seurat_data, resolution = x)
  return(heatmap)
})

pdf(here("results/reviewer_figs/similarity_heatmaps.pdf"))


for(heatmap in pc_heatmaps){
  print(heatmap)
  grid::grid.newpage()
}

for(heatmap in resolution_heatmaps){
  print(heatmap)
  grid::grid.newpage()
}
dev.off()
