# Author: Matthew Aaron Loberg
# Date: 23-0608
# Script: Fig_4E_Visium_ATC_SpaCET_Lymphoid_Deconvolution.R

# Please look at the following script: "Figure3_Visium_ATC_SpaCET_PANCAN_Deconvolution.R" prior to this script

#### A note on this script: ####
# This script takes the output from "Figure3_Visium_ATC_SpaCET_PANCAN_Deconvolution.R" and calculates the capture area frequency of Lymphoid deconvolution
# Here, lymphoid deconvolution is defined as the sum of Plasma cells, B cells, CD4 T cells, CD8 T cells, and NK cells
# Lymphoid deconvolution results are then plotted spatially
# Finally, an average capture area %lymphoid is calculated
# This script was run for all 8 ATC Visium Spatial Transcriptomics in our manuscript

#### Load packages ####
library(tidyverse) # using tidyverse for ggplot and to work with tibbles
library(SpaCET) # Analysis package
library(patchwork) # A few of the ggplots require patchwork for adding plots together

#### Read in saved visium object with SpaCET deconvolution already performed ####
# Deconvolution performed in the following script: "Figure3_Visium_ATC_SpaCET_PANCAN_Deconvolution.R"
# Read in ATC_SpaCET_Object
ATC_SpaCET_Object <- readRDS(file = "Data_in_Use/ATC_Visium_Object_SpaCET_PANCAN_Deconvolved.RDS")

#### Extract cell type proportion matrix from ATC_SpaCET_Object ####
# Find the dimensions of the proportion matrix
dim(ATC_SpaCET_Object@results$deconvolution$propMat)

# Print out the rownames of the dimension matrix
rownames(ATC_SpaCET_Object@results$deconvolution$propMat)
# Above indicates that there are 35 rows (corresponding to cell types)
# Row 1 = Malignant
# Row 2 = CAF
# etc. 

# Extract proportion matrix to a new variable
ATC_SpaCET_Object_Prop_Mat <- ATC_SpaCET_Object@results$deconvolution$propMat

#### Add a lymphocyte column to the proportion matrix ####
## Lymphotyes will be the following: 
# Plasma (row 4)
# B cell (row 5)
# T CD4 (row 6)
# T CD8 (row 7)
# NK (row 8)

# First, make a lymphoid matrix that is a combination of all of the other cell types
Lymphoid_Matrix <- as.matrix(ATC_SpaCET_Object_Prop_Mat[4,] + # Plasma
                             ATC_SpaCET_Object_Prop_Mat[5,] + # B Cell
                             ATC_SpaCET_Object_Prop_Mat[6,] + # T CD4
                             ATC_SpaCET_Object_Prop_Mat[7,] + # T CD8
                             ATC_SpaCET_Object_Prop_Mat[8,]) # NK

# Need to flip the new Lymphoid Matrix
Lymphoid_Matrix <- t(Lymphoid_Matrix)

# Set the rowname of Lymphoid_Matrix to "Lymphoid"
rownames(Lymphoid_Matrix) <- c("Lymphoid")
rownames(Lymphoid_Matrix) # View the rownames of lymphoid matrix

# View the first 5 columns of Lymphoid_Matrix
Lymphoid_Matrix[,1:5] # Just making sure it is oriented as desired

# Now bind the Lymphoid_Matrix to the ATC_SpaCET_Object_Prop_Mat
ATC_SpaCET_Object_Prop_Mat <- rbind(ATC_SpaCET_Object_Prop_Mat, Lymphoid_Matrix)
# Check that lymphoid was added as the 36th row name
rownames(ATC_SpaCET_Object_Prop_Mat)
# Success - row 36 is lymphoid

# Clean up the Lymphoid_Matrix
rm(Lymphoid_Matrix)

# Set the SpaCET_obj to have the new Prop_Mat with lymphoid as the 36th row
# This will allow us to plot lymphocytes
ATC_SpaCET_Object@results$deconvolution$propMat <- ATC_SpaCET_Object_Prop_Mat

#### Save a spatial Lymphoid Deconvolution Plot ####
# Plot spatial distribution of Lymphoid deconvolution
ATC_SpaCET_Object_Lymphoid_Fraction_Plot <- SpaCET.visualize.spatialFeature(ATC_SpaCET_Object,
                                                                            spatialType = "CellFraction",
                                                                            spatialFeatures = c("Lymphoid"))
ggsave(file = "outputs/ATC_Visium_Object/Deconvolution/ATC_Visium_Object_SpaCET_Lymphoid_Deconvolution_Plot.png",
       ATC_SpaCET_Object_Lymphoid_Fraction_Plot,
       width = 9, height = 5, dpi = 600)
rm(ATC_SpaCET_Object_Lymphoid_Fraction_Plot)



#### Calculate Lymphoid Population Frequency as an Average of Capture Area Frequencies ####

# Make a variable for outputting population frequencies that is a tibble w/ 1 row and 36 columns
ATC_Frequency_Output <- as_tibble(matrix(data = 0, nrow = 1, ncol = nrow(ATC_SpaCET_Object_Prop_Mat)))

# Make the colnames of the frequency output tibble the rownames (population names) of the deconvolution matrix
colnames(ATC_Frequency_Output) <- rownames(ATC_SpaCET_Object_Prop_Mat)

# Calculate the average frequency of each deconvoluted cell population across the sample
# 1. Sum the total deconvoluted amount of all of the spatial capture arease
for(i in 1:nrow(ATC_SpaCET_Object_Prop_Mat)){ # i = rows of PropMat (1 - 36) - 36th row = lymphocyte group I just added
  for(n in 1:ncol(ATC_SpaCET_Object_Prop_Mat)){ # n = cols of PropMat (the spatial capture areas to be summed)
    ATC_Frequency_Output[1,i] <- ATC_Frequency_Output[1,i] + ATC_SpaCET_Object_Prop_Mat[i,n]
  }
}
# 2.Divide Frequency output by the number of spatial capture areas to get a mean frequency
# 3. Multiply by 100 to make it a percentage
ATC_Frequency_Output <- ATC_Frequency_Output/ncol(ATC_SpaCET_Object_Prop_Mat)*100

#### Save Data & Clean Up ####
# Save the frequency tibble for combining with other data/merging with other samples
saveRDS(ATC_Frequency_Output, file = "Data_in_Use/ATC_Visium_Object_SpaCET_Deconvolved_Frequency_Averages_Lymphoid_Included.RDS")

# Cleaning Up
rm(list = ls())



