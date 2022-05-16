library(Seurat)
library(tidyverse)
library(cowplot)
library(here)
library(scAnalysisR)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

sample <- "healthy_bcells_all"

base_dir_proj <- file.path(base_dir, "results", sample)
save_dir <- file.path(base_dir_proj, "R_analysis")


tenx_structure <- "count" # Set to multi if you ran cell ranger multi or
# count if you ran cell ranger count

if(normalization_method == "SCT"){
  SCT <- TRUE
  seurat_assay <- "SCT"
} else {
  SCT <- FALSE
  seurat_assay <- "RNA"
}

ADT <- TRUE
HTO <- FALSE

vars.to.regress <- NULL

# Set directories
base_dir <- here()

mt_pattern <- "^MT-" # "^MT-" for human, "^mt-" for mice

# Create seurat object
seurat_object <- create_seurat_object(sample = sample,
                                      count_path = file.path(base_dir,
                                                             "results"),
                                      ADT = ADT, hashtag = HTO,
                                      tenx_structure = tenx_structure
)

# Add mitochondrial percent
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object,
                                                      pattern = mt_pattern)

# Quality plots to determin cutoffs
rna_qual <- VlnPlot(seurat_object,
                    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                    ncol = 3)

rna_qual

if(ADT){
  adt_qual <- VlnPlot(seurat_object, features = c("nCount_ADT", "nFeature_ADT"))
  adt_qual
}


# mark outliers

seurat_object$passing <- "not_passing"

seurat_object$passing[seurat_object$percent.mt < 15 & 
                        seurat_object$nFeature_RNA > 500 &
                        seurat_object$nFeature_RNA < 5000 &
                        seurat_object$nCount_ADT < 2000] <- "passing"


# Default normalization
DefaultAssay(seurat_object) <- "RNA"
seurat_object <- NormalizeData(seurat_object) %>% 
  FindVariableFeatures() %>%
  ScaleData(vars.to.regress = vars.to.regress)

plot(featDistPlot(seurat_object, geneset = c("IRF4", "JCHAIN", "XBP1"),
                  sep_by = "passing"))


data_info <- FetchData(seurat_object, vars = c("IRF4", "JCHAIN", "XBP1",
                                             "percent.mt", "nFeature_RNA",
                                             "nCount_ADT"))

cutoff_info <- c("nFeature_RNA" = 500,
                 "percent.mt" = 10,
                 "nCount_ADT" = 2000)

make_plots <- function(data_info, cutoff_info, gene){
  plot_list <- lapply(names(cutoff_info), function(variable){
    ggplot(data_info, aes_(x = as.name(gene), y = as.name(variable))) +
      geom_point() +
      ggplot2::geom_hline(yintercept = cutoff_info[[variable]])
  })
  return(plot_list)
}

all_plots <- make_plots(data_info, cutoff_info, "IRF4")

jchain_plots <- make_plots(data_info, cutoff_info, "JCHAIN")

xbp1_plots <- make_plots(data_info, cutoff_info, "XBP1")

# PCA and UMAP
set.seed(0)
DefaultAssay(seurat_object) <- seurat_assay
seurat_object <- PCA_dimRed(seurat_object, assay = seurat_assay)
npcs <- 20

umap_data <- group_cells(seurat_object, assay = seurat_assay, nPCs = npcs)

seurat_object <- umap_data[[1]]

plotDimRed(seurat_object, col_by = "IRF4", plot_type = "rna.umap")
plotDimRed(seurat_object, col_by = "XBP1", plot_type = "rna.umap")
plotDimRed(seurat_object, col_by = "JCHAIN", plot_type = "rna.umap")



# Check doublet removal
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_adt.rds"))

plot(featDistPlot(seurat_data, geneset = c("IRF4", "JCHAIN", "XBP1"),
                  sep_by = "Doublet_finder"))
