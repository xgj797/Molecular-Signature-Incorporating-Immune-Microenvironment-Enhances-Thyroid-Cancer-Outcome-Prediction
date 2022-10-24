# Author: Matthew Aaron Loberg
# Date: 22-0913
# Goal: read in pre-processed Sample8 and add TDS as a module score
# Visiualize TDS spatially

# Load packages:
library(tidyverse)
library(Seurat)

# Load processed Sample8
Sample8 <- readRDS(file = "Data_in_Use/Sample8/Sample8_Clustered_0.3Res.rds")

### Add the TDS Score as a module score
# Create a gene list for TDS (16 genes total)
TDS_Genes <- c("DIO1", "DIO2", "DUOX1", "DUOX2", "FOXE1", 
               "GLIS3", "NKX2-1", "PAX8", "SLC26A4", "SLC5A5", 
               "SLC5A8", "TG", "THRA", "THRB", "TPO", "TSHR")

# Turn into a list
TDS_List <- list(TDS_Genes)
rm(TDS_Genes) # clean up TDS_Genes

# Add TDS score to Sample8                                
Sample8 <- AddModuleScore(object = Sample8,
                          features = TDS_List,
                          assay = "SCT",
                          name = 'TDS')

# Print TDS Spatial Feature Plot
TDS_Spatial <- SpatialFeaturePlot(Sample8, features = c("TDS1"))
ggsave("outputs/Sample8_ModuleScores/22-0913_Sample8_TDS_ModuleScore.png",
       TDS_Spatial,
       width = 4, height = 5, dpi = 600)
