# Matthew Loberg
# GO Analysis

## 22-1019 Update
# Making dot plots for the 549 gene MAP and for the 673 RAS-like upregulated genes

# Main Figure 2G (MAP dot plots) + Supplemental Figure 3H (RAS dot plots)
# Gene ontology analysis performed on geneontology.org
# Results are plotted here as dot plots

## Note: if anywhere there is ATS, it is the same as MAP
# ATS = aggressive tumor score
# MAP = molecular aggression and prediction score
# With resubmission of our manuscript, we changed the name from ATS -> MAP, but some of our code still refers to this score as ATS

# Load Packages
library(ggplot2)
library(tidyverse)
library(stringr)
library(forcats)
library(RColorBrewer)


#### Extracellular Select #####
# load data + format
GO_Biological_Extracellular <- read_csv(file = "data_in_use/22-1018_GO_Biological_Select_Extracellular.csv")
GO_Biological_Extracellular <- GO_Biological_Extracellular %>% dplyr::rename(Fold_Enrichment = 'upload_1 (fold Enrichment)')
GO_Biological_Extracellular$Fold_Enrichment <- log(GO_Biological_Extracellular$Fold_Enrichment, 2)
GO_Biological_Extracellular <- GO_Biological_Extracellular %>% dplyr::rename(FDR = 'upload_1 (FDR)')
GO_Biological_Extracellular <- GO_Biological_Extracellular %>% dplyr::rename(Biological_Process = 'GO biological process complete')
GO_Biological_Extracellular <- GO_Biological_Extracellular %>% dplyr::rename(REFLIST_GENES = 'Homo sapiens - REFLIST (20589)')
GO_Biological_Extracellular <- GO_Biological_Extracellular %>% dplyr::rename(DETECTED_GENES = 'upload_1 (547)')
GO_Biological_Extracellular$GeneRatio <- GO_Biological_Extracellular$DETECTED_GENES/GO_Biological_Extracellular$REFLIST_GENES
GO_Biological_Extracellular$Biological_Process <- str_remove_all(GO_Biological_Extracellular$Biological_Process, "[()01234567891GO:]")

# Add an Order variable
GO_Biological_Extracellular$Order <- 0
for(i in 1:nrow(GO_Biological_Extracellular)){
  GO_Biological_Extracellular$Order[i] <- nrow(GO_Biological_Extracellular) + 1 - i
}

# Plot the GO Biological data for extracellular matrix
DotPlot_GO_Biological_Extracellular <- ggplot(GO_Biological_Extracellular, aes(Fold_Enrichment, fct_reorder(Biological_Process, Order))) + 
  geom_point(aes(size = GeneRatio, fill = FDR), shape = 21, stroke = 1) +
  theme_bw(base_size = 10) +
  scale_size_continuous(range = c(2, 10), limits = c(0,0.5), name = "Gene Ratio") +
  scale_fill_distiller(limits=c(0,.05), type="div", palette= "RdYlBu", direction = 1) +
  scale_x_continuous(name="log2(Fold Enrichment)", seq(0, 5, by = 1), NULL, expand = c(0,0), limits=c(0, 5), NULL) +
  theme(aspect.ratio = 1.5, 
        #legend.position = "none",
        axis.title.y = element_blank()) +
  ggtitle("") + 
  theme(plot.title = element_text(face = "bold", size = 12),
        axis.text = element_text(face = "bold", size = 11),
        axis.text.y = element_blank(),
        title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold", size = 12))
# geom_vline(xintercept = 0, size = 0.25)
ggsave("outputs/22-1019_GO_Analysis_Enrichment_Plots/22-1019_GO_Biological_ATS_549_Extracellular.png", height = 5, width = 3, DotPlot_GO_Biological_Extracellular, dpi = 600)


#### Keratin select ####
# load data + format
GO_Biological_Keratin <- read_csv(file = "data_in_use/22-1018_549_ATS_GO_Biological_Process_Keratin_5.csv")
GO_Biological_Keratin <- GO_Biological_Keratin %>% dplyr::rename(Fold_Enrichment = 'upload_1 (fold Enrichment)')
GO_Biological_Keratin$Fold_Enrichment <- log(GO_Biological_Keratin$Fold_Enrichment, 2)
GO_Biological_Keratin <- GO_Biological_Keratin %>% dplyr::rename(FDR = 'upload_1 (FDR)')
GO_Biological_Keratin <- GO_Biological_Keratin %>% dplyr::rename(Biological_Process = 'GO biological process complete')
GO_Biological_Keratin <- GO_Biological_Keratin %>% dplyr::rename(REFLIST_GENES = 'Homo sapiens - REFLIST (20589)')
GO_Biological_Keratin <- GO_Biological_Keratin %>% dplyr::rename(DETECTED_GENES = 'upload_1 (547)')
GO_Biological_Keratin$GeneRatio <- GO_Biological_Keratin$DETECTED_GENES/GO_Biological_Keratin$REFLIST_GENES
GO_Biological_Keratin$Biological_Process <- str_remove_all(GO_Biological_Keratin$Biological_Process, "[()01234567891GO:]")

# Add an Order variable
GO_Biological_Keratin$Order <- 0
for(i in 1:nrow(GO_Biological_Keratin)){
  GO_Biological_Keratin$Order[i] <- nrow(GO_Biological_Keratin) + 1 - i
}

# Plot the GO Biological data for Keratin matrix
DotPlot_GO_Biological_Keratin <- ggplot(GO_Biological_Keratin, aes(Fold_Enrichment, fct_reorder(Biological_Process, Order))) + 
  geom_point(aes(size = GeneRatio, fill = FDR), shape = 21, stroke = 1) +
  theme_bw(base_size = 10) +
  scale_size_continuous(range = c(2, 10), limits = c(0,0.5), name = "Gene Ratio") +
  scale_fill_distiller(limits=c(0,.05), type="div", palette= "RdYlBu", direction = 1) +
  scale_x_continuous(name="log2(Fold Enrichment)", seq(0, 5, by = 1), NULL, expand = c(0,0), limits=c(0, 5), NULL) +
  theme(aspect.ratio = 1.5, 
        #legend.position = "none",
        axis.title.y = element_blank()) +
  ggtitle("") + 
  theme(plot.title = element_text(face = "bold", size = 12),
        axis.text = element_text(face = "bold", size = 11),
        axis.text.y = element_blank(),
        title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold", size = 12))
# geom_vline(xintercept = 0, size = 0.25)
ggsave("outputs/22-1019_GO_Analysis_Enrichment_Plots/22-1019_GO_Biological_ATS_549_Keratin_5.png", height = 5, width = 3, DotPlot_GO_Biological_Keratin, dpi = 600)

#### Immune ####
# load data + format
GO_Biological_Immune <- read_csv(file = "data_in_use/22-1018_549_ATS_GO_Biological_Process_Immune_5.csv")
GO_Biological_Immune <- GO_Biological_Immune %>% dplyr::rename(Fold_Enrichment = 'upload_1 (fold Enrichment)')
GO_Biological_Immune$Fold_Enrichment <- log(GO_Biological_Immune$Fold_Enrichment, 2)
GO_Biological_Immune <- GO_Biological_Immune %>% dplyr::rename(FDR = 'upload_1 (FDR)')
GO_Biological_Immune <- GO_Biological_Immune %>% dplyr::rename(Biological_Process = 'GO biological process complete')
GO_Biological_Immune <- GO_Biological_Immune %>% dplyr::rename(REFLIST_GENES = 'Homo sapiens - REFLIST (20589)')
GO_Biological_Immune <- GO_Biological_Immune %>% dplyr::rename(DETECTED_GENES = 'upload_1 (549)')
GO_Biological_Immune$GeneRatio <- GO_Biological_Immune$DETECTED_GENES/GO_Biological_Immune$REFLIST_GENES
GO_Biological_Immune$Biological_Process <- str_remove_all(GO_Biological_Immune$Biological_Process, "[()01234567891GO:]")

# Add an Order variable
GO_Biological_Immune$Order <- 0
for(i in 1:nrow(GO_Biological_Immune)){
  GO_Biological_Immune$Order[i] <- nrow(GO_Biological_Immune) + 1 - i
}

# Plot the GO Biological data for Immune matrix
DotPlot_GO_Biological_Immune <- ggplot(GO_Biological_Immune, aes(Fold_Enrichment, fct_reorder(Biological_Process, Order))) + 
  geom_point(aes(size = GeneRatio, fill = FDR), shape = 21, stroke = 1) +
  theme_bw(base_size = 10) +
  scale_size_continuous(range = c(2, 10), limits = c(0,0.5), name = "Gene Ratio") +
  scale_fill_distiller(limits=c(0,.05), type="div", palette= "RdYlBu", direction = 1) +
  scale_x_continuous(name="Log2(Fold Enrichment)", seq(0, 5, by = 1), NULL, expand = c(0,0), limits=c(0, 5), NULL) +
  theme(aspect.ratio = 1.5, 
        #legend.position = "none",
        axis.title.y = element_blank()) +
  ggtitle("") + 
  theme(plot.title = element_text(face = "bold", size = 12),
        axis.text = element_text(face = "bold", size = 11),
        axis.text.y = element_blank(),
        title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold", size = 12))
# geom_vline(xintercept = 0, size = 0.25)
ggsave("outputs/22-1019_GO_Analysis_Enrichment_Plots/22-1019_GO_Biological_ATS_549_Immune_5.png", height = 5, width = 3, DotPlot_GO_Biological_Immune, dpi = 600)

#### Cell Cycle ####
# load data + format
GO_Biological_Cycle <- read_csv(file = "data_in_use/22-1018_549_ATS_GO_Biological_Process_Cell_Cycle_5.csv")
GO_Biological_Cycle <- GO_Biological_Cycle %>% dplyr::rename(Fold_Enrichment = 'upload_1 (fold Enrichment)')
GO_Biological_Cycle$Fold_Enrichment <- log(GO_Biological_Cycle$Fold_Enrichment, 2)
GO_Biological_Cycle <- GO_Biological_Cycle %>% dplyr::rename(FDR = 'upload_1 (FDR)')
GO_Biological_Cycle <- GO_Biological_Cycle %>% dplyr::rename(Biological_Process = 'GO biological process complete')
GO_Biological_Cycle <- GO_Biological_Cycle %>% dplyr::rename(REFLIST_GENES = 'Homo sapiens - REFLIST (20589)')
GO_Biological_Cycle <- GO_Biological_Cycle %>% dplyr::rename(DETECTED_GENES = 'upload_1 (549)')
GO_Biological_Cycle$GeneRatio <- GO_Biological_Cycle$DETECTED_GENES/GO_Biological_Cycle$REFLIST_GENES
GO_Biological_Cycle$Biological_Process <- str_remove_all(GO_Biological_Cycle$Biological_Process, "[()01234567891GO:]")

# Add an Order variable
GO_Biological_Cycle$Order <- 0
for(i in 1:nrow(GO_Biological_Cycle)){
  GO_Biological_Cycle$Order[i] <- nrow(GO_Biological_Cycle) + 1 - i
}

# Plot the GO Biological data for Cycle matrix
DotPlot_GO_Biological_Cycle <- ggplot(GO_Biological_Cycle, aes(Fold_Enrichment, fct_reorder(Biological_Process, Order))) + 
  geom_point(aes(size = GeneRatio, fill = FDR), shape = 21, stroke = 1) +
  theme_bw(base_size = 10) +
  scale_size_continuous(range = c(2, 10), limits = c(0,0.5), name = "Gene Ratio") +
  scale_fill_distiller(limits=c(0,.05), type="div", palette= "RdYlBu", direction = 1) +
  scale_x_continuous(name="Log2(Fold Enrichment)", seq(0, 5, by = 1), NULL, expand = c(0,0), limits=c(0, 5), NULL) +
  theme(aspect.ratio = 1.5, 
        #legend.position = "none",
        axis.title.y = element_blank()) +
  ggtitle("") + 
  theme(plot.title = element_text(face = "bold", size = 12),
        axis.text = element_text(face = "bold", size = 11),
        axis.text.y = element_blank(),
        title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold", size = 12))
# geom_vline(xintercept = 0, size = 0.25)
ggsave("outputs/22-1019_GO_Analysis_Enrichment_Plots/22-1019_GO_Biological_ATS_549_Cycle_5.png", height = 5, width = 3, DotPlot_GO_Biological_Cycle, dpi = 600)

### Adding in RAS-upregulated thyroid ###
# Note: this is supplemental Figure 3H
# This data is from upregulated RAS-gene list
# load data + format
GO_Biological_Thyroid <- read_csv(file = "data_in_use/22-1019_RAS_LIKE_NOT_AGGRESSIVE_GO_Thyroid.csv")
GO_Biological_Thyroid <- GO_Biological_Thyroid %>% dplyr::rename(Fold_Enrichment = 'upload_1 (fold Enrichment)')
GO_Biological_Thyroid$Fold_Enrichment <- log(GO_Biological_Thyroid$Fold_Enrichment, 2)
GO_Biological_Thyroid <- GO_Biological_Thyroid %>% dplyr::rename(FDR = 'upload_1 (FDR)')
GO_Biological_Thyroid <- GO_Biological_Thyroid %>% dplyr::rename(Biological_Process = 'GO biological process complete')
GO_Biological_Thyroid <- GO_Biological_Thyroid %>% dplyr::rename(REFLIST_GENES = 'Homo sapiens - REFLIST (20589)')
GO_Biological_Thyroid <- GO_Biological_Thyroid %>% dplyr::rename(DETECTED_GENES = 'upload_1 (673)')
GO_Biological_Thyroid$GeneRatio <- GO_Biological_Thyroid$DETECTED_GENES/GO_Biological_Thyroid$REFLIST_GENES
GO_Biological_Thyroid$Biological_Process <- str_remove_all(GO_Biological_Thyroid$Biological_Process, "[()01234567891GO:]")

# Add an Order variable
GO_Biological_Thyroid$Order <- 0
for(i in 1:nrow(GO_Biological_Thyroid)){
  GO_Biological_Thyroid$Order[i] <- nrow(GO_Biological_Thyroid) + 1 - i
}

# Plot the GO Biological data for Thyroid matrix
DotPlot_GO_Biological_Thyroid <- ggplot(GO_Biological_Thyroid, aes(Fold_Enrichment, fct_reorder(Biological_Process, Order))) + 
  geom_point(aes(size = GeneRatio, fill = FDR), shape = 21, stroke = 1) +
  theme_bw(base_size = 10) +
  scale_size_continuous(range = c(2, 10), limits = c(0,0.7), name = "Gene Ratio") +
  scale_fill_distiller(limits=c(0,.05), type="div", palette= "RdYlBu", direction = 1) +
  scale_x_continuous(name="Log2(Fold Enrichment)", seq(0, 5, by = 1), NULL, expand = c(0,0), limits=c(0, 5), NULL) +
  theme(aspect.ratio = 1.5, 
        #legend.position = "none",
        axis.title.y = element_blank()) +
  ggtitle("") + 
  theme(plot.title = element_text(face = "bold", size = 12),
        axis.text = element_text(face = "bold", size = 11),
        axis.text.y = element_blank(),
        title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold", size = 12))
# geom_vline(xintercept = 0, size = 0.25)
ggsave("outputs/22-1019_GO_Analysis_Enrichment_Plots/22-1019_GO_Biological_RAS_Thyroid_4.png", height = 5, width = 3, DotPlot_GO_Biological_Thyroid, dpi = 600)
