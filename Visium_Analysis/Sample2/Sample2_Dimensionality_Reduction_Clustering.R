# Author: Matthew Aaron Loberg
# Date: August 25th, 2022
# Purpose: New Visium sequencing data just obtained from Vantage
# Here, I will perform dimensionality reduction analysis of Sample2
# I will do this by reading in the SCTransformed Sample2 data generated from the following script: 
# Sample2_Processing_SpaceRanger_Output.R

#### Chapter 1: Load required packages ####
# Load required packages
library(Seurat)
library(hdf5r) # required to read in data file
library(ggplot2)
library(dplyr)
library(patchwork)
library(tidyverse)

#### Chapter 2: Load SCTransformed data and perform dimensionality reduction ####
# Load in SCTransformed Sample2 object 
Sample2 <- readRDS(file = "Data_in_Use/Sample2/Sample2_SCTransformed_All_Genes.rds")

## Principal component analysis + choosing appropriate number of components
# Use an elbow plot to determine the dimensionality of the data. Note: must do this AFTER RunPCA
# First, run PCA on Sample2 w/ the RunPCA command
Sample2 <- RunPCA(Sample2, verbose = FALSE)
Elbow <- ElbowPlot(Sample2) # Make Elbow plot
Elbow <- Elbow + theme(
  axis.text = element_text(face = "bold", size = 25),
  axis.title = element_text(face = "bold", size = 30)
)
Elbow # Print elbow plot
ggsave("outputs/Sample2_QC/Sample2_Dimensionality_Reduction/Sample2_SCTransformed_Sample2_Elbow_Plot.png",
       Elbow, 
       width = 5, height = 5, dpi = 600)
# From here I need to choose the appropriate # of dimensions 
# I'm going to go with 10 based on the elbow plot

#### Chapter 3: Clustering ####
# For dimensionality reduction, I will use the number of dimensions identified by the elblow plot (10)
Sample2 <- FindNeighbors(Sample2, dims = 1:10)
# Resolution of find clusters should be determined based on the tissue morphology, per Luciane
Sample2 <- FindClusters(Sample2, verbose = FALSE, resolution = 0.3) # note, changing to a resolution of 0.3 to decrease the # of clusters 
Sample2 <- RunUMAP(Sample2, dims = 1:10)
UMAP_Plot <- DimPlot(Sample2, reduction = "umap", group.by = c("ident"))
UMAP_Plot # Print UMAP_Plot
UMAP_Plot <- UMAP_Plot +
  theme(
    axis.title = element_text(face = "bold", size = 25),
    axis.text = element_text(face = "bold", size = 20),
    legend.text = element_text(face = "bold", size = 20)
  )
ggsave("outputs/Sample2_QC/Sample2_Dimensionality_Reduction/Sample2_UMAP_ClusterResolution_0.3.png", 
       UMAP_Plot, 
       width = 7, height = 5, dpi = 600)


#### Chapter 4: Save RDS for Future Use ####
# Save the new version of Sample2 that has PCA, UMAP, and clustering data
saveRDS(Sample2, file = "Data_in_Use/Sample2/Sample2_Clustered_0.3Res.rds")

#### Chapter 5: Cleaning up ####
# remove all data from the global environment
rm(list = ls())
