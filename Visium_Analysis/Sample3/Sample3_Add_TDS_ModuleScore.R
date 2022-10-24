# Author: Matthew Aaron Loberg
# Date: 22-0913
# Goal: read in pre-processed Sample3 and add TDS as a module score
# Visiualize TDS spatially

# Load packages:
library(tidyverse)
library(Seurat)

# Load processed Sample3
Sample3 <- readRDS(file = "Data_in_Use/Sample3/Sample3_Clustered_0.3Res.rds")

### Add the TDS Score as a module score
# Create a gene list for TDS (16 genes total)
TDS_Genes <- c("DIO1", "DIO2", "DUOX1", "DUOX2", "FOXE1", 
               "GLIS3", "NKX2-1", "PAX8", "SLC26A4", "SLC5A5", 
               "SLC5A8", "TG", "THRA", "THRB", "TPO", "TSHR")

# Turn into a list
TDS_List <- list(TDS_Genes)
rm(TDS_Genes) # clean up TDS_Genes

# Add TDS score to Sample3                                
Sample3 <- AddModuleScore(object = Sample3,
                          features = TDS_List,
                          assay = "SCT",
                          name = 'TDS')

# Print TDS Spatial Feature Plot
TDS_Spatial <- SpatialFeaturePlot(Sample3, features = c("TDS1"))
ggsave("outputs/Sample3_ModuleScores/22-0913_Sample3_TDS_ModuleScore.png",
       TDS_Spatial,
       width = 4, height = 5, dpi = 600)
