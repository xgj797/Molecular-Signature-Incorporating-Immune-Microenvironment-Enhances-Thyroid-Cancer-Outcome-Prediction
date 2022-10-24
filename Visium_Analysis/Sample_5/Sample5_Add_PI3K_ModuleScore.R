# Author: Matthew Aaron Loberg
# Date: 22-0913
# Goal: read in pre-processed Sample5 and add PI3K-AKT-mTOR as a module score
# Visiualize PI3K-AKT-mTOR spatially

# Load packages:
library(tidyverse)
library(Seurat)

# Load processed Sample5
Sample5 <- readRDS(file = "Data_in_Use/Sample5/Sample5_Clustered_0.3Res.rds")

### Applying PI3K Score Gene Set 
PI3K_Score_Genes <- 
  read.table("Data_in_Use/22-0314_HALLMARK_PI3K_AKT_MTOR_SIGNALING.txt", sep = "\t")
PI3K_Score_Genes <- as_tibble(PI3K_Score_Genes)
colnames(PI3K_Score_Genes) <- PI3K_Score_Genes[1,] # Set col names = to first row
PI3K_Score_Genes <- PI3K_Score_Genes[3:nrow(PI3K_Score_Genes),] # Start list with first column as first gene (remove first two rows)
# Need to make into a list for Visium analysis
PI3K_List <- list(PI3K_Score_Genes$HALLMARK_PI3K_AKT_MTOR_SIGNALING)

# Add PI3K score to thyroid merge
Sample5 <- AddModuleScore(object = Sample5,
                          features = PI3K_List,
                          name = 'PI3K_Score')

# Print PI3K Spatial Feature Plot
PI3K_Spatial <- SpatialFeaturePlot(Sample5, features = c("PI3K_Score1"))
ggsave("outputs/Sample5_ModuleScores/Sample5_PI3K_ModuleScore.png",
       PI3K_Spatial,
       width = 4, height = 5, dpi = 600)
