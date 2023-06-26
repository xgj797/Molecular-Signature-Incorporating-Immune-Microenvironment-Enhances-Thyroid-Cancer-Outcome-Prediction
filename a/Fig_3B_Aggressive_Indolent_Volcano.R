# Matthew Loberg
# DESEQ Results Analysis - Volcano plots of Aggressive/Indolent DESeq2 analysis

### Load packages
library(tidyverse)

### Load data sets
All_DESEQ_Aggressive <- readRDS(file = "data_in_use/22-0908_DESEQ2_All_Primary_LocalDisease_MG_HT_Excluded_resOrdered_Aggressive_Indolent.rds")

### Volcano Plot - load packages
library(dplyr)
library(ggplot2)
library(ggrepel)
library(EnhancedVolcano)

## All plot aggressive
keyvals <- ifelse(
  All_DESEQ_Aggressive$log2FoldChange < -1 & All_DESEQ_Aggressive$padj < 0.05, 'cornflowerblue',
  ifelse(All_DESEQ_Aggressive$log2FoldChange > 1 & All_DESEQ_Aggressive$padj < 0.05, 'firebrick',
         'grey'))

keyvals[is.na(keyvals)] <- 'grey'

names(keyvals)[keyvals == 'firebrick'] <- 'Up'
names(keyvals)[keyvals == 'cornflowerblue'] <- 'Down'
names(keyvals)[keyvals == 'grey'] <- 'NS'

AggressiveVolcano <- EnhancedVolcano(All_DESEQ_Aggressive,
                                      lab = rownames(All_DESEQ_Aggressive),
                                      legendPosition = "none",
                                      x = 'log2FoldChange',
                                      y = 'padj',
                                      pCutoff = 0.05,
                                      FCcutoff = 1.0, 
                                      pointSize = 1.5,
                                      colCustom = keyvals,
                                      title = "Aggressive/Indolent\nDifferential Gene Expression Analysis",
                                      subtitle = NULL,
                                      titleLabSize = 30,
                                      axisLabSize = 25,
                                      # subtitleLabSize = 24,
                                      # labsize = 10,
                                      selectLab = c("FAP",
                                                    # "CXCL12",
                                                    "COL9A3",
                                                    "COL5A2",
                                                    "COL6A3",
                                                    "CXCL5",
                                                    "COL7A1",
                                                    "FGF5",
                                                    "TERT",
                                                    "POSTN",
                                                    "TG",
                                                    "TSHR",
                                                    "TPO",
                                                    # "FN1",
                                                    # "ASPN", # Hurley lab gene of interest
                                                    # "HK3",
                                                    "VCAN",
                                                    "COL11A1",
                                                    "ACTA2",
                                                    "WNT2",
                                                    "LRRC15",
                                                    # "WNT5A",
                                                    "PDPN", 
                                                    "COL1A1",
                                                    "COL1A2",
                                                    "TMSB4X",
                                                    "SLC5A5"),
                                      boxedLabels = TRUE, 
                                      drawConnectors = TRUE, 
                                      colConnectors = 'black',
                                      widthConnectors = 1.0,
                                      labSize = 7)
#lengthConnectors = 1)


AggressiveVolcano <- 
  AggressiveVolcano + theme(axis.title = element_text(face = "bold", size = 40),
                             axis.text = element_text(face = "bold", size = 30))
# ggsave the Aggressive volcano plot
ggsave("outputs/22-0913_DESEQ_Results_Analysis/22-0913_Aggressive_Volcano.png", width = 10, height = 10, AggressiveVolcano, dpi = 600)
