# Author: Matthew Aaron Loberg
# Date: August 25th, 2022
# Purpose: New Visium sequencing data just obtained from Vantage
# Here, I will perform dimensionality reduction analysis of Sample8
# I will do this by reading in the SCTransformed Sample8 data generated from the following script: 
# Sample8_Processing_SpaceRanger_Output.R

#### Chapter 1: Load required packages ####
# Load required packages
library(Seurat)
library(hdf5r) # required to read in data file
library(ggplot2)
library(dplyr)
library(patchwork)
library(tidyverse)

#### Chapter 2: Load SCTransformed data and perform dimensionality reduction ####
# Load in SCTransformed Sample8 object 
Sample8 <- readRDS(file = "Data_in_Use/Sample8/Sample8_SCTransformed_All_Genes.rds")

## Principal component analysis + choosing appropriate number of components
# Use an elbow plot to determine the dimensionality of the data. Note: must do this AFTER RunPCA
# First, run PCA on Sample8 w/ the RunPCA command
Sample8 <- RunPCA(Sample8, verbose = FALSE)
Elbow <- ElbowPlot(Sample8) # Make Elbow plot
Elbow <- Elbow + theme(
  axis.text = element_text(face = "bold", size = 25),
  axis.title = element_text(face = "bold", size = 30)
)
Elbow # Print elbow plot
ggsave("outputs/Sample8_QC/Sample8_Dimensionality_Reduction/Sample8_SCTransformed_Sample8_Elbow_Plot.png",
       Elbow, 
       width = 5, height = 5, dpi = 600)
# From here I need to choose the appropriate # of dimensions 
# I'm going to go with 10 based on the elbow plot

#### Chapter 3: Clustering ####
# For dimensionality reduction, I will use the number of dimensions identified by the elblow plot (10)
Sample8 <- FindNeighbors(Sample8, dims = 1:10)
# Resolution of find clusters should be determined based on the tissue morphology, per Luciane
Sample8 <- FindClusters(Sample8, verbose = FALSE, resolution = 0.3) # note, changing to a resolution of 0.3 to decrease the # of clusters 
Sample8 <- RunUMAP(Sample8, dims = 1:10)
UMAP_Plot <- DimPlot(Sample8, reduction = "umap", group.by = c("ident"))
UMAP_Plot # Print UMAP_Plot
UMAP_Plot <- UMAP_Plot +
  theme(
    axis.title = element_text(face = "bold", size = 25),
    axis.text = element_text(face = "bold", size = 20),
    legend.text = element_text(face = "bold", size = 20)
  )
ggsave("outputs/Sample8_QC/Sample8_Dimensionality_Reduction/Sample8_UMAP_ClusterResolution_0.3.png", 
       UMAP_Plot, 
       width = 7, height = 5, dpi = 600)


#### Chapter 4: Save RDS for Future Use ####
# Save the new version of Sample8 that has PCA, UMAP, and clustering data
saveRDS(Sample8, file = "Data_in_Use/Sample8/Sample8_Clustered_0.3Res.rds")

#### Chapter 5: Cleaning up ####
# remove all data from the global environment
rm(list = ls())
