# Author: Matthew Aaron Loberg
# Date: 22-0913
# Goal: read in pre-processed Sample6 and add fibroblast module scores (normal, myCAF, iCAF)
# Visiualize fibroblast scores (normal, myCAF, iCAF) spatially

# Load packages:
library(tidyverse)
library(Seurat)

# Load processed Sample6
Sample6 <- readRDS(file = "Data_in_Use/Sample6/Sample6_Clustered_0.3Res.rds")

# Read in fibroblast gene sets
Kieffer_Combined <- 
  read.table("Data_in_Use/21-1028_Kieffer_CAF_Subsets_Combined.txt", sep = "\t")

# Subset out normal fibroblasts
Kieffer_Normal_Fibroblast <- list(Kieffer_Combined[2:nrow(Kieffer_Combined), 1])

# Subset out iCAF Combined list
Kieffer_iCAF_Combined <- list(Kieffer_Combined[2:nrow(Kieffer_Combined), 2])
Kieffer_iCAF_Combined <- list(Kieffer_iCAF_Combined[[1]][1:66]) # Restrict iCAF combined to the appropriate length

# Subset out myCAF combined list
Kieffer_myCAF_Combined <- list(Kieffer_Combined[2:nrow(Kieffer_Combined), 3])
Kieffer_myCAF_Combined <- list(Kieffer_myCAF_Combined[[1]][1:75]) # Restrict myCAF combined to the appropriate length

### Add in normal fibroblasts
Sample6 <- AddModuleScore(object = Sample6,
                          features = Kieffer_Normal_Fibroblast,
                          name = 'Kieffer_Normal_Fibroblast')
NormalFibroblastSpatial <- SpatialFeaturePlot(Sample6, features = c("Kieffer_Normal_Fibroblast1"))
ggsave("outputs/Sample6_Normal_Fibroblast_ModuleScore.png",
       NormalFibroblastSpatial,
       width = 4, height = 5, dpi = 600)

### Add in myCAFs
Sample6 <- AddModuleScore(object = Sample6,
                          features = Kieffer_myCAF_Combined,
                          name = 'Kieffer_myCAF_Combined')

# Print out spatial plot for myCAF
MyCAFFibroblastSpatial <- SpatialFeaturePlot(Sample6, features = c("Kieffer_myCAF_Combined1"))
ggsave("outputs/Sample6_myCAF_ModuleScore.png",
       MyCAFFibroblastSpatial,
       width = 4, height = 5, dpi = 600)


### Add in iCAFs
Sample6 <- AddModuleScore(object = Sample6, 
                          features = Kieffer_iCAF_Combined,
                          name = 'Kieffer_iCAF_Combined')


# Print out spatial plot for iCAF
iCAFFibroblastSpatial <- SpatialFeaturePlot(Sample6, features = c("Kieffer_iCAF_Combined1"))
ggsave("outputs/Sample6_iCAF_ModuleScore.png",
       iCAFFibroblastSpatial,
       width = 4, height = 5, dpi = 600)
