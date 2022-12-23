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

# Read in the data
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed.rds"))


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

plotDimRed(seurat_data, col_by = "new_cluster", color = cluster_colors_two,
           plot_type = "rna.umap")


# BND2 cells = 6.1

# Basically IgD-class switched B cells should be CD27+ (both protein and mRNA) 
# and not express IgM mRNA, which is different from our BND2 cells, which 
# should be CD27- and can express some IgM mRNA, just not a lot by surface 
# protein expression.

# Class switched:
  # CD27+ (protein + RNA)
  # IgD +
  # IgM - (RNA)
# BND2 
  # CD27 - (protein and RNA)
  # IgD +
  # IgM lo (RNA)

seurat_subset <- subset(seurat_data, subset = new_cluster == "6.1")

# IgD vs IgM antibody
# IgD antibody vs CD27 RNA
# IgM antibody vs CD27 RNA
# IgD antibody vs IgM RNA
# CD27 RNA vs IgM RNA

DefaultAssay(seurat_subset) <- "RNA"
all_data <- FetchData(seurat_subset, vars = c("IGHM", "IGHD", "CD27"))

DefaultAssay(seurat_subset) <- "ADT"

new_data <- FetchData(seurat_subset, vars = c("IgM", "IgD", 
                                              "CD27.1"))

metadata <- seurat_subset[[]] %>%
  dplyr::select(IgM_cutoff, IgD_cutoff, CD21_cutoff, CD27_cutoff, pass_filters)
                
identical(rownames(all_data), rownames(new_data))

all_data <- do.call(cbind, list(all_data, new_data, metadata))

make_plot <- function(var_one, var_two, color = NULL){
  # add color for double positive, single positive and negative?
  
  if(is.null(color)){
    base_plot <- ggplot2::ggplot(all_data, aes_string(x = var_one, y = var_two)) 
  } else {
    base_plot <- ggplot2::ggplot(all_data, aes_string(x = var_one, y = var_two,
                                                      color = color)) 
  }
  
  base_plot +
    ggplot2::geom_point() 
}


# IgD vs IgM antibody
make_plot("IgD", "IgM")

# IgD antibody vs CD27 RNA
make_plot("IgD", "CD27")

# IgM antibody vs CD27 RNA
make_plot("IgM", "CD27")

# IgD antibody vs IgM RNA
make_plot("IgD", "IGHM")

# CD27 RNA vs IgM RNA
make_plot("IGHM", "CD27")


make_plot("IGHM", "CD27", color = "pass_filters")

make_plot("IGHM", "CD27", color = "IgM_cutoff")

make_plot("IGHM", "CD27", color = "IgD_cutoff")

make_plot("IGHM", "CD27", color = "CD27_cutoff")


table(all_data$pass_filters)


all_data$cd27_expr <- ifelse(all_data$CD27 > 0, TRUE, FALSE)

# Number that are positive for CD27 expression
table(all_data$cd27_expr) %>%
  data.frame %>% 
  mutate(total = sum(Freq)) %>%
  mutate(fraction = Freq / total)

make_plot("IGHM", "CD27", color = "cd27_expr")

