### Sample1 Spatial Feature plots: Information
# This file reads in SCTransformed Sample1 spatial transcriptomics data and outputs spatial feature plots for FAP and ACTA2

### Chapters: 1-3
# 1. Load packages
# 2. Read in processed Sample1
# 3. Make spatial feature plots of interest (FAP, ACTA2)

### Chapter 1: Load Packages
library(Seurat)
library(patchwork)
library(tidyverse)

### Chapter 2: Read in processed Sample1
Sample1 <- readRDS(file = "data_in_Use/Sample1/Sample1_SCTransformed_All_Genes.rds")

### Chapter 3: Make spatial feature plots of interest (FAP, ACTA2)
## FAP
Sample1_FAP <- SpatialFeaturePlot(Sample1, features = "FAP", pt.size.factor = 3)
ggsave("outputs/Sample1_SpatialFeaturePlots/Sample1_FAP.png",
       Sample1_FAP,
       width = 4, height = 5, dpi = 600)
rm(Sample1_FAP)

## ACTA2
Sample1_ACTA2 <- SpatialFeaturePlot(Sample1, features = "ACTA2", pt.size.factor = 3)
ggsave("outputs/Sample1_SpatialFeaturePlots/Sample1_ACTA2.png",
       Sample1_ACTA2,
       width = 4, height = 5, dpi = 600)
rm(Sample1_ACTA2)

