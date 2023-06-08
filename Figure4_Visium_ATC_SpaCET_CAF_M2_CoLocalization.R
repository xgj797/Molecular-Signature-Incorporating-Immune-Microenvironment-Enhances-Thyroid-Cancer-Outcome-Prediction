# Author: Matthew Aaron Loberg
# Date: 23-0608
# Script: 23-0329_Thy8_SpaCET_PANCAN.R

# Please look at the following script: "Figure3_Visium_ATC_SpaCET_PANCAN_Deconvolution.R" prior to this script

#### A note on this script: ####
# This script takes the output from "Figure3_Visium_ATC_SpaCET_PANCAN_Deconvolution.R" and looks at co-localization of M2 macrophages and CAFs
# This script was run for all 8 ATC Visium Spatial Transcriptomics in our manuscript
# First, colocalization is calculated for all SpaCET populations
# Second, an LRInteraction score is calculated
# Third, the colocalization of CAFs and M2 macrophages is plotted

#### Load packages ####
library(tidyverse) # using tidyverse for ggplot and to work with tibbles
library(SpaCET) # Analysis package
library(patchwork) # A few of the ggplots require patchwork for adding plots together

#### Read in saved visium object with SpaCET deconvolution already performed ####
# Deconvolution performed in the following script: "Figure3_Visium_ATC_SpaCET_PANCAN_Deconvolution.R"
# Read in ATC_SpaCET_Object
ATC_SpaCET_Object <- readRDS(file = "Data_in_Use/ATC_Visium_Object_SpaCET_PANCAN_Deconvolved.RDS")

#### Look at co-localization of CAFs and M2 macrophages ####
# Calculate the cell-cell colocalization of all 35 deconvoluted cell types
ATC_SpaCET_Object <- SpaCET.CCI.colocalization(ATC_SpaCET_Object)

# Calculate the ligand-receptor network score across ST spots (required for making the co-localization plot shown if Figure 4D)
ATC_SpaCET_Object <- SpaCET.CCI.LRNetworkScore(ATC_SpaCET_Object, coreNo = 1)

# Ligand-Receptor analysis for a co-localized cell-type pair (CAF and Macrophage M2)
ATC_SpaCET_Object <- SpaCET.CCI.cellTypePair(ATC_SpaCET_Object, cellTypePair=c("CAF","Macrophage M2"))

# Visualize the interaction analysis of a co-localized cell-type pair (CAFs and Macrophage M2).
ATC_SpaCET_CAF_M2_LR_Colocalization_Plot <- SpaCET.visualize.cellTypePair(ATC_SpaCET_Object, cellTypePair=c("CAF","Macrophage M2"))
ggsave(file = "outputs/ATC_Visium_Object/CoLocalization/ATC_Visium_Object_SpaCET_CAF_M2_Localization_Plot.png",
       ATC_SpaCET_CAF_M2_LR_Colocalization_Plot,
       width = 20, height = 7, dpi = 600)
rm(ATC_SpaCET_CAF_M2_LR_Colocalization_Plot)

#### Cleaning up ####
rm(list = ls())

