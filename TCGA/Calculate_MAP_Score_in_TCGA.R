# Author: Matthew Aaron Loberg
# Date: 23-0608 (date that script was cleaned up for the GitHub for our resubmission to Cell Genomics)
# Script: Calculate_MAP_Score_in_TCGA.R

#### Notes on this script ####
# This script will take the 549 genes that are in the molecular aggression and prediction (MAP) score and use them to calculate a MAP score for the samples in the TCGA papillary thyroid carcinoma (PTC) cohort
# The TCGA PTC cohort acts as a validation cohort for us throughout the paper
# Calculating 

#### Load Required Packages ####
library(tidyverse)

