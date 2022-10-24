# Author: Matthew Aaron Loberg
# Date: 22-0913
# Goal: read in pre-processed Sample3 and add BRAF as a module score
# Visiualize BRAF module score spatially

# Load packages:
library(tidyverse)
library(Seurat)

# Load processed Sample3
Sample3 <- readRDS(file = "Data_in_Use/Sample3/Sample3_Clustered_0.3Res.rds")

### Add the BRS Score
## For this I will focus only on genes that are "UP" in BRAF lesions but NOT genes that are "UP" in RAS lesions
## "BRAF Score"
BRS_Genes <- as_tibble(read_csv(file = "Data_in_Use/BRS_Gene_List_71_Annotated.csv"))
BRAF_Genes <- BRS_Genes %>% subset(Annotation == "BRAF")
BRAF_Genes_List <- list(BRAF_Genes$BRS_Gene_list)


# Add BRAF Gene Score to Sample3
Sample3 <- AddModuleScore(object = Sample3,
                          features = BRAF_Genes_List,
                          name = "BRAF_Score")

# Print BRAF Spatial Feature Plot
BRAF_Spatial <- SpatialFeaturePlot(Sample3, features = c("BRAF_Score1"))
ggsave("outputs/Sample3_ModuleScores/Sample3_BRAF_ModuleScore.png",
       BRAF_Spatial,
       width = 4, height = 5, dpi = 600)
