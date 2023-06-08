# Author: Matthew Aaron Loberg
# Date: 23-0608
# Script: Fig_3F_Visium_Seurat_Clustering_Analysis.R

#### A note on this script ####
# This script reads in data for a representative anaplastic thyroid carcinoma Visium Spatial Transcriptomics Sample
# Then, the script uses Seurat to perform quality control, dimensionality reduction, and clustering analysis
# Finally, this script finds the marker genes of individual clusters and prints out a heatmap with the top 10 markers of each cluster

#### Load required packages ####
library(Seurat)
library(hdf5r) # required to read in data file
library(patchwork)
library(tidyverse)

#### Setting Directory ####
# Set the data directory to the location of the 10X Genomics SpaceRanger output files for the visium sample of interest
# For Seurat, this should be one folder which contains at minimum the following two subfolders:
# 1. A subfolder with spatial output files ("spatial")
# This spatial subfolder should contain the following files: "aligned_fiducials.jpg", "detected_tissue_image.jpg", "scalefactors_json.json", "spatial_enrichment.csv", "tissue_hires_image.png", tissue_lowres_image.png", "tissue_positions_list.csv"
# 2. A subfolder ("filtered_feature_bc_matrix") containing "barcodes.tsv.gz", "features.tsv.gz", and "matrix.mtz.gz"
Dir <- "Data_in_Use/ATC_Visium_Sample_Location" # Set the directory to load from


#### Reading in ATC_Visium_Obj ####

# Load in ATC_Visium_Obj data
ATC_Visium_Obj <- Load10X_Spatial(data.dir = Dir, slice = "slice1") # Load ATC_Visium_Obj
# Set orig.ident
ATC_Visium_Obj$orig.ident <- "ATC_Visium_Obj" # Set the ident. This way, if ever merged with another object, this object maps back to its identity as ATC_Visium_Obj
# Cleaning up
rm(data_dir)

#### Visualize QC Data (raw count data as violin + spatial) ####
# Visualize raw count data as a violin plot. Note, 'raster = FALSE' for image quality
QC_Violin <- VlnPlot(ATC_Visium_Obj, features = "nCount_Spatial", raster = FALSE) + NoLegend()

# Visiualize spatial location of raw count data
QC_Spatial <- SpatialFeaturePlot(ATC_Visium_Obj, features = "nCount_Spatial") + theme(legend.position = "right")

# Format and save raw reads violin plot
QC_Violin <- QC_Violin + 
  theme(axis.text = element_text(face = "bold", size = 15))
ggsave("outputs/ATC_Visium_Object/QC_Plots/ATC_Visium_Object_Seurat_Raw_Violin_QC.png",
       QC_Violin,
       width = 4, height = 5, dpi = 600)

# Format and save raw reads spatial plot
QC_Spatial <- QC_Spatial + theme(
  legend.text = element_text(face = "bold", size = 15),
  legend.title = element_text(face = "bold", size = 15))
ggsave("outputs/ATC_Visium_Object/QC_Plots/ATC_Visium_Object_Seurat_Raw_Spatial_QC.png",
       QC_Spatial,
       width = 7, height = 5, dpi = 600)

# Cleaning up
rm(QC_Violin, QC_Spatial)

#### Data Transformation ####
# I will perform data transformation with SCTransform
# See SCTransform vignette here: https://satijalab.org/seurat/articles/sctransform_vignette.html
# See SCTransform paper here: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1874-1
# There are many discussions on whether return.only.var.genes should be set to TRUE/FALSE.
# In the vignette from satijalab, they recommend FALSE for best performance. 
ATC_Visium_Obj <- SCTransform(ATC_Visium_Obj, assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)

#### Dimensionality Reduction analysis + Clustering ####
## Principal component analysis + choosing appropriate number of components
# Use an elbow plot to determine the dimensionality of the data. Note: must do this AFTER RunPCA
# First, run PCA on ATC_Visium_Obj w/ the RunPCA command
ATC_Visium_Obj <- RunPCA(ATC_Visium_Obj, verbose = FALSE)
Elbow <- ElbowPlot(ATC_Visium_Obj) # Make Elbow plot
Elbow <- Elbow + theme(
  axis.text = element_text(face = "bold", size = 25),
  axis.title = element_text(face = "bold", size = 30)
)
Elbow # Print elbow plot
ggsave("outputs/ATC_Visium_Object/QC_Plots/ATC_Visium_Object_Seurat_PCA_Dimensionality_Elbow.png",
       Elbow, 
       width = 5, height = 5, dpi = 600)
# From here I need to choose the appropriate # of dimensions 
# I'm going to go with 10 based on the elbow plot
# For dimensionality reduction, I will use the number of dimensions identified by the elblow plot (10)

## Find Neighbors
# Dims set based on PCA elbow plot
ATC_Visium_Obj <- FindNeighbors(ATC_Visium_Obj, dims = 1:10)

## Find Clusters
# Res set to 0.2 based on tissue architecture
ATC_Visium_Obj <- FindClusters(ATC_Visium_Obj, verbose = FALSE, resolution = 0.2)

## Run UMAP
# Dims set based on PCA elbow plot
ATC_Visium_Obj <- RunUMAP(ATC_Visium_Obj, dims = 1:10)

# Make a umap plot
UMAP_Plot <- DimPlot(ATC_Visium_Obj, reduction = "umap", group.by = c("ident"))
UMAP_Plot # Print UMAP_Plot

# Format UMAP plot
UMAP_Plot <- UMAP_Plot +
  theme(
    axis.title = element_text(face = "bold", size = 25),
    axis.text = element_text(face = "bold", size = 20),
    legend.text = element_text(face = "bold", size = 20)
  )

# Save UMAP plot
ggsave("outputs/ATC_Visium_Object/Dimensionality_Reduction/ATC_Visium_Object_Seurat_UMAP_2D_Dim10_Res0.2.png", 
       UMAP_Plot, 
       width = 7, height = 5, dpi = 600)

## Map clusters onto spatial plot of the tissue section
# Adjust pt.size.factor as needed to make the size of dots desired (here maximizing so that take up full spatial capture area)
SpatialDim <- SpatialDimPlot(ATC_Visium_Obj, label = TRUE, label.size = 15, pt.size.factor = 3) #& NoLegend()
# Save spatial plot of clusters
ggsave("outputs/ATC_Visium_Object/Dimensionality_Reduction/ATC_Visium_Object_Seurat_SpatialClusters_Dim10_Res0.2.png",
       SpatialDim,
       width = 10, height = 8, dpi = 600)

#### Find cluster markers + Make heatmap of top 10 markers ####
# Find all markers
ATC_Visium_Obj.Markers <- FindAllMarkers(ATC_Visium_Obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Sort out top 10 markers
ATC_Visium_Obj.Markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

# Make heatmap of top 10 markers of each cluster
Heatmap <- DoHeatmap(ATC_Visium_Obj, features = top10$gene) + NoLegend()
ggsave("outputs/ATC_Visium_Object/Dimensionality_Reduction/ATC_Visium_Object_Seurat_ClusterMarkers_Top10_Dim10_Res0.2.png",
       Heatmap,
       width = 10, height = 10, dpi = 600)

#### cleaning up ####
rm(list = ls())

