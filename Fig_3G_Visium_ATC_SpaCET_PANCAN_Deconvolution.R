# Author: Matthew Aaron Loberg
# Date: 23-0608
# Script: Fig_3G_Visium_ATC_SpaCET_PANCAN_Deconvolution.R

#### A note on this script: ####
# This script contains information for running SpaCET analysis on one anaplastic thyroid carcinoma (ATC) Visium Spatial Transcriptomics sample using output data from the 10X Genomics SpaceRanger pipeline
# This script was run for all 8 ATC Visium Spatial Transcriptomics in our manuscript
# While there is a thyroid cancer (THCA) option for copy number alterations in the detection of malignant cell frequency, THCA copy number alterations are based on differentiated thyroid carcinoma and do not resemble what is seen in ATC 
# Thus, we used the PANCAN option to better portray what would be expected to be seen in ATC 

#### A note on SpaCET: ####
# Here are available tutorials on using SpaCET:
# Paper here: https://www.nature.com/articles/s41467-023-36062-6
# Github here: https://github.com/data2intelligence/SpaCET
# Tutorial here: https://data2intelligence.github.io/SpaCET/articles/visium_BC.html (Cell type deconvolution and interaction analysis)
# Tutorial here: https://data2intelligence.github.io/SpaCET/articles/oldST_PDAC.html (Deconvolution with a matched scRNA-seq data set)
# SpaCET deconvolutes the following MAJOR cell fractions (these should add up to one total):
# Malignant, CAF, Endothelial, Plasma, B cell, T CD4, T CD8, NK, cDC, pDC, Macrophage, Mast, Neutrophil, Unidentifiable
# SpaCET deconvolutes the following MINOR cell fractions:
# B cell naive, B cell non-switched memory, B cell switched memory, B cell exhausted, T CD4 naive, Th1, Th2, Th17, Tfh, Treg, TCD8 naive, TCD8 central memory, T CD8 effector memory, T CD8 effector, T CD8 exhausted, cDC1 CLEC9A, cDC2 CD1C, cDC3 LAMP3, Macrophage M1, Macrophage M2, Macrophage other

#### Load packages ####
library(tidyverse) # contains ggplot, other packages
library(SpaCET) # Analysis package for deconvolution
library(patchwork) # A few of the ggplots require patchwork for adding plots together
library(beepr) # I will be adding beepr to let me know when parts of the code finish running

#### Setting Directory ####
# Set the data directory to the location of the 10X Genomics SpaceRanger output files
# For SpaCET, this should be one folder which contains at minimum the following two subfolders:
# 1. A subfolder with spatial output files ("spatial")
# This spatial subfolder should contain the following files: "aligned_fiducials.jpg", "detected_tissue_image.jpg", "scalefactors_json.json", "spatial_enrichment.csv", "tissue_hires_image.png", tissue_lowres_image.png", "tissue_positions_list.csv"
# 2. A subfolder ("filtered_feature_bc_matrix") containing "barcodes.tsv.gz", "features.tsv.gz", and "matrix.mtz.gz"
Dir <- "Data_in_Use/ATC_Visium_Sample_Location"

#### Read in Data ####
# Use the Visium specific SpaCET function to read in an ATC Visium object at the indicated directory
ATC_Visium_Object <- create.SpaCET.object.10X(visiumPath = Dir)
# Clean up directory path variable
rm(Dir)

#### QC Metrics ####
# Calculate QC metrics
ATC_Visium_Object <- SpaCET.quality.control(ATC_Visium_Object) 

# Plot the QC Metrics
# Note: the way that this command is written, patchwork must be loaded as "library(patchwork)" -> added to load packages above
# This is not indicated in the SpaCET tutorial; however, I found that it only works if patchwork is loaded
# See stack overflow post here: https://stackoverflow.com/questions/62599918/ggplot2-cant-add-to-a-ggplot-object
ATC_Visium_Object_QC_Metrics <- SpaCET.visualize.spatialFeature(
  ATC_Visium_Object,
  spatialType = "QualityControl",
  spatialFeatures = c("UMI", "Gene")
)
ggsave(filename = "outputs/ATC_Visium_Object/QC_Plots/ATC_Visium_Object_QC_Metrics.png",
       ATC_Visium_Object_QC_Metrics,
       width = 20, height = 5, dpi = 600)
# Clean up
rm(ATC_Visium_Object_QC_Metrics)

#### Deconvolve Spatial Transcriptomics Data with SpaCET ####

### Background:
## See info on deconvolution process here: https://data2intelligence.github.io/SpaCET/articles/visium_BC.html
## In brief, SpaCET deconvolution process works as follows:
# 1. SpaCET first estimates malignant cell fraction by a gene pattern dictionary of copy number alterations (CNA) and expression changes (based on those copy number alterations) in common malignancies
# 2. Based on a SpaCET in-house hierarchical cell lineage derived from single-cell RNA-seq data sets from diverse cancer types, immune and stromal cell fractions are determined
# For this to work, users need to specify the cancer type of this tumor Spatial Transcriptomics data set by using cancerType parameter
# The cancerType parameter selects cancer type-specific CNA or expression signatures to infer malignant cell fraction. 
# One of the available cancer types is "THCA" (thyroid cancer) -> this does not work well for ATC as it is likely based on differentiated thyroid tumors
# For cancer types NOT included in the SpaCET dictionarly, there is a pan-cancer expression signature
# The pan-cancer signature was created by averaging all cancer type-specific expression signatures
# Here, I will use PANCAN for ATC

# deconvolve ATC_Visium_Object using "PANCAN" cancer type
ATC_Visium_Object <- SpaCET.deconvolution(ATC_Visium_Object, cancerType = "PANCAN", coreNo = 1)

# Save the deconvoluted object for future scripts: 
saveRDS(ATC_Visium_Object, file = "Data_in_Use/ATC_Visium_Object_SpaCET_PANCAN_Deconvolved.RDS")

## Notes regarding results: 
# Values represent the fraction of cell types
# The fraction sum of all cell types in a spot (column) may be beyond 1
# This is because SpaCET outputs the cell fractions for both major lineages and their corresponding sublineages (e.g., T CD4 as well as Th1, Th2)

#### Visualize the cell type proportion for CAFs and Macrophages ####

# Plot spatial distribution of CAFs
ATC_Visium_Object_CAF_Fraction_Plot <- SpaCET.visualize.spatialFeature(ATC_Visium_Object,
                                                                       spatialType = "CellFraction",
                                                                       spatialFeatures = c("CAF"))
ggsave(file = "outputs/ATC_Visium_Object/Deconvolution/ATC_Visium_Object_SpaCET_CAF_Deconvolution_Plot.png",
       ATC_Visium_Object_CAF_Fraction_Plot,
       width = 9, height = 5, dpi = 600)
rm(ATC_Visium_Object_CAF_Fraction_Plot)

# Plot spatial distribution of Macrophages
ATC_Visium_Object_Macrophage_Fraction_Plot <- SpaCET.visualize.spatialFeature(ATC_Visium_Object,
                                                                              spatialType = "CellFraction",
                                                                              spatialFeatures = c("Macrophage"))
ggsave(file = "outputs/ATC_Visium_Object/Deconvolution/ATC_Visium_Object_SpaCET_Macrophage_Deconvolution_Plot.png",
       ATC_Visium_Object_Macrophage_Fraction_Plot,
       width = 9, height = 5, dpi = 600)
rm(ATC_Visium_Object_Macrophage_Fraction_Plot)

#### Cleaning up ####
rm(ATC_Visium_Object)
