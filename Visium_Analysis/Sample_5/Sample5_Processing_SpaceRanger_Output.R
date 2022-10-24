### Sample5 Processing: Information
# This script reads in Sample5 data, visualizes the raw counts, and performs SCTransform.
# The SCTransformed R object is saved for further downstream analyses. 
# Note: This sequencing performed at Vanderbilt University Medical Center (VUMC). Pilot run at JHU.

### Chapters: 1-3
# 1. Load packages
# 2. Load sample and visualize raw counts
# 3. SCTransform

#### Chapter 1: Loading Packages ####
# Load required packages
library(Seurat)
library(hdf5r) # required to read in data file
library(ggplot2)
library(dplyr)
library(patchwork)
library(tidyverse)

### Chapter 2: Reading in Sample5 and looking at raw count data by violin and SpatialFeaturePlot ####

# Load in Sample5 Data
data_dir <- 'data_in_Use/Sample5' # Set directory to load from
Sample5 <- Load10X_Spatial(data.dir = data_dir, slice = "slice1") # Load Sample 1
# Cleaning up
rm(data_dir)

# I was having a problem with spatial feature plot coordinates being "characters" instead of "integers" (note: with data sequenced at VUMC, not JHU)
# When I ran spatial feature plot I was getting the following error: 
# "Error in FUN(left, right) : non-numeric argument to binary operator
# According to stack overflow, the following code should fix the issue 
# See line: https://stackoverflow.com/questions/73131436/error-in-funleft-right-non-numeric-argument-to-binary-operator-when-runni
Sample5@images[["slice1"]]@coordinates[["tissue"]] <- as.integer(Sample5@images[["slice1"]]@coordinates[["tissue"]])
Sample5@images[["slice1"]]@coordinates[["row"]] <- as.integer(Sample5@images[["slice1"]]@coordinates[["row"]])
Sample5@images[["slice1"]]@coordinates[["col"]] <- as.integer(Sample5@images[["slice1"]]@coordinates[["col"]])
Sample5@images[["slice1"]]@coordinates[["imagerow"]] <- as.integer(Sample5@images[["slice1"]]@coordinates[["imagerow"]])
Sample5@images[["slice1"]]@coordinates[["imagecol"]] <- as.integer(Sample5@images[["slice1"]]@coordinates[["imagecol"]])

# Visualize raw count data as a violin plot. Note, 'raster = FALSE' for image quality
plot1 <- VlnPlot(Sample5, features = "nCount_Spatial", raster = FALSE) + NoLegend()

# Visiualize spatial location of raw count data
plot2 <- SpatialFeaturePlot(Sample5, features = "nCount_Spatial") + theme(legend.position = "right")

# Format and save plot1 and plot2 with raw count data
plot1 <- plot1 + 
  theme(
    axis.text = element_text(face = "bold", size = 15)
  )
ggsave("outputs/Sample5/Raw_Counts_Violin.png",
       plot1,
       width = 4, height = 5, dpi = 600)

plot2 <- plot2 + theme(
  legend.text = element_text(face = "bold", size = 15),
  legend.title = element_text(face = "bold", size = 15)
)
ggsave("outputs/Sample5/Raw_Counts_Spatial.png",
       plot2,
       width = 7, height = 5, dpi = 600)

# Cleaning up
rm(plot1, plot2, test)

# Save raw Sample5 as an RDS
saveRDS(Sample5, "data_in_Use/Sample5/Sample5_Raw.rds")

#### Chapter 3: Data Transformation ####
# I will perform data transformation with SCTransform
# See SCTransform vignette here: https://satijalab.org/seurat/articles/sctransform_vignette.html
# See SCTransform paper here: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1874-1
# There are many discussions on whether return.only.var.genes should be set to TRUE/FALSE.
# In the vignette from satijalab, they recommend FALSE for best performance. 
Sample5 <- SCTransform(Sample5, assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)

# Save SCTransformed Sample5 as an RDS
saveRDS(Sample5, "data_in_Use/Sample5/Sample5_SCTransformed_All_Genes.rds")

# Cleaning up
rm(Sample5)

# I will start subsequent analysis by loading the SCTransformed/processed version of Sample5
