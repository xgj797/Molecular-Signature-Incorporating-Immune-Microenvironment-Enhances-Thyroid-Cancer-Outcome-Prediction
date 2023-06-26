# Author: Matthew Aaron Loberg
# Date: 23-0608 (date that script was cleaned up for the GitHub for our resubmission to Cell Genomics)
# Script: Calculate_MAP_Score_in_TCGA_Validation_Cohort.R

#### Notes on this script ####
# This script will take the 549 genes that are in the molecular aggression and prediction (MAP) score and use them to calculate a MAP score for the samples in the TCGA papillary thyroid carcinoma (PTC) cohort
# The TCGA PTC cohort acts as a validation cohort for us throughout the paper
# Due to differences in gene coverage, the TCGA MAP Score will be based on 520 of the 549 genes that were identified in our cohort

#### Load Required Packages ####
library(tidyverse)

#### Load TCGA Patient and Sample Data ####

# Load patient data file -> make into tibble for easier cleaning
TCGA_Patient_Data <- as_tibble(read.table(file = "data_in_use/22-0107_TCGA_data_clinical_patient.txt", sep = '\t', header = TRUE))

# Load sample data file -> make into tibble for easier cleaning
TCGA_Sample_Data <- as_tibble(read.table(file = "data_in_use/22-0107_TCGA_data_clinical_sample.txt", sep = '\t', header = TRUE)) # Load sample data file

#### Merge TCGA Patient and Sample Data ####

# PATIEND_ID has an "01" appended to the end of each sample. Here, I will remove it from each sample. 
TCGA_Sample_Data$PATIENT_ID <- substr(TCGA_Sample_Data$PATIENT_ID, 1, nchar(TCGA_Sample_Data$PATIENT_ID)-3)

# Combine sample data and patient data
Patient_Sample_Combined <- as_tibble(merge(TCGA_Sample_Data, TCGA_Patient_Data, by.x = "PATIENT_ID", by.y = "PATIENT_ID", all = T))

# Subset to only samples with RNA Sequencing data (496 -> 482 samples)
Patient_Sample_Combined <- Patient_Sample_Combined %>% subset(RNASEQ_DATA == "Yes")

# Clean up environment - only leave the data that I need (Patient_Sample_Combined)
remove(list = c("TCGA_Patient_Data", "TCGA_Sample_Data"))

#### Load TCGA sequencing RSEM Expression Data ####

# Read in TCGA RNA Data (Z-Score format) and clean the data to the desired format
TCGA_RNA_Z_Scores <- read.table(file = "data_in_use/data_RNA_Seq_v2_mRNA_median_Zscores.txt", sep = '\t', header = FALSE) # Read in TCGA RNA data in Z-Score format
Gene_Symbols <- TCGA_RNA_Z_Scores$V1 # Create variable for storing Gene Symbols -> will become column names
TCGA_RNA_Z_Scores <- TCGA_RNA_Z_Scores[,3:ncol(TCGA_RNA_Z_Scores)] # Delete first two columns, don't need these
TCGA_Transposed <- t(TCGA_RNA_Z_Scores) # Transpose the data so that each column pertains to a particular gene
colnames(TCGA_Transposed) <- Gene_Symbols # Make Gene Symbols the column names
TCGA_RNA_Tibble <- as_tibble(TCGA_Transposed) # Convert type to tibble
TCGA_RNA_Tibble <- dplyr::rename(TCGA_RNA_Tibble, SAMPLE_ID = Hugo_Symbol) # dplyr::rename "Hugo_Symbol" to "SAMPLE_ID" to match Patient_Sample_Combined for future merge

# Cleaning up: remove the data I no longer need
remove(list = c("TCGA_RNA_Z_Scores", "TCGA_Transposed", "Gene_Symbols"))

#### Format MAP Score gene list (549) for TCGA Analysis ####

# Read in 549 gene BRAF-Poor outcome list (MAP Score genes) generated on 22-0617
MAPScore_Genes <- read_csv(file = "data_in_use/22-0927_All_BRAF_All_Aggressive_Overlap_Genes.csv")

# I need to make sure that this gene set is compatible with the TCGA data
# A number of the genes in the 549 gene MAP Score are not present in the TCGA Sequencing. I will remove those here:
# There are 18 genes not present in TCGA that are in our cohort. This makes the TCGA MAP score gene list 531 (of the 549 from our cohort)
MAPScore_Genes_Filtered <- MAPScore_Genes
MAPScore_Genes_Filtered <- MAPScore_Genes_Filtered %>% filter(value != "LRRC38" &
                                                              value != "ARL14EPL" &
                                                              value != "NT5DC4" &
                                                              value != "IZUMO3" &
                                                              value != "STRIT1" &
                                                              value != "KLF18" &
                                                              value != "SMIM31" &
                                                              value != "TAF11L11" &
                                                              value != "OOSP1" &
                                                              value != "GOLGA6L22" &
                                                              value != "KRTAP2-3" &
                                                              value != "PRSS56" &
                                                              value != "ACOD1" &
                                                              value != "OR5BS1P" &
                                                              value != "LACTBL1" &
                                                              value != "LY6L" &
                                                              value != "MDFIC2" &
                                                              value != "RFPL4AL1")

# A number of the genes use different gene symbols in the TCGA Dataset
# I will change the names of the TCGA dataset here to match the MAP Score Gene List

# No H1-5 -> found old symbol name, HIST1H1B
TCGA_RNA_Tibble <- TCGA_RNA_Tibble %>% dplyr::rename('H1-5' = HIST1H1B)
# No H3C2 -> found old symbol name, HIST1H3B
TCGA_RNA_Tibble <- TCGA_RNA_Tibble %>% dplyr::rename(H3C2 = HIST1H3B)
# No H2BC17 -> found old symbol name, HIST1H2BO
TCGA_RNA_Tibble <- TCGA_RNA_Tibble %>% dplyr::rename(H2BC17 = HIST1H2BO)
# No PLPP4 -> found old symbol name, PPAPDC1A
TCGA_RNA_Tibble <- TCGA_RNA_Tibble %>% dplyr::rename(PLPP4 = PPAPDC1A)
# No H2AC16 -> found old symbol name, HIST1H2AL
TCGA_RNA_Tibble <- TCGA_RNA_Tibble %>% dplyr::rename(H2AC16 = HIST1H2AL)
# No ACP7 -> found old symbol name, PAPL
TCGA_RNA_Tibble <- TCGA_RNA_Tibble %>% dplyr::rename(ACP7 = PAPL)
# No LORICRIN -> found old symbol name, LOR
TCGA_RNA_Tibble <- TCGA_RNA_Tibble %>% dplyr::rename(LORICRIN = LOR)
# No OPRPN -> found old symbol name, PROL1
TCGA_RNA_Tibble <- TCGA_RNA_Tibble %>% dplyr::rename(OPRPN = PROL1)
# No CCDC197 -> found old symbol name, LINC00521
TCGA_RNA_Tibble <- TCGA_RNA_Tibble %>% dplyr::rename(CCDC197 = LINC00521)
# No TAFA4 -> found old symbol name, FAM19A4
TCGA_RNA_Tibble <- TCGA_RNA_Tibble %>% dplyr::rename(TAFA4 = FAM19A4)
# No OR8U3 -> found old symbol name, OR5R1
TCGA_RNA_Tibble <- TCGA_RNA_Tibble %>% dplyr::rename(OR8U3 = OR5R1)
# No H2AC4 -> found old symbol name, HIST1H2AB
TCGA_RNA_Tibble <- TCGA_RNA_Tibble %>% dplyr::rename(H2AC4 = HIST1H2AB)
# No H3C12 -> found old symbol name, HIST1H3J
TCGA_RNA_Tibble <- TCGA_RNA_Tibble %>% dplyr::rename(H3C12 = HIST1H3J)
# No H3C3 -> found old symbol name, HIST1H3C
TCGA_RNA_Tibble <- TCGA_RNA_Tibble %>% dplyr::rename(H3C3 = HIST1H3C)
# No ADGRF4 -> found old symbol name, GPR115
TCGA_RNA_Tibble <- TCGA_RNA_Tibble %>% dplyr::rename(ADGRF4 = GPR115)
# No CALHM6 -> found old symbol name, FAM26F
TCGA_RNA_Tibble <- TCGA_RNA_Tibble %>% dplyr::rename(CALHM6 = FAM26F)
# No RTL3 -> found old symbol name, ZCCHC5
TCGA_RNA_Tibble <- TCGA_RNA_Tibble %>% dplyr::rename(RTL3 = ZCCHC5)

# Restrict TCGA data by the filtered MAP Score Genes and make into a numeric
TCGA_RNA_Tibble_MAP <- TCGA_RNA_Tibble[c(MAPScore_Genes_Filtered$value)] # Restrict genes
MAPScore_Gene_Names <- colnames(TCGA_RNA_Tibble_MAP) # save col names for future use
TCGA_RNA_Matrix_MAP <- as.matrix(TCGA_RNA_Tibble_MAP) # Make into a matrix
TCGA_RNA_Matrix_MAP <- matrix(as.numeric(TCGA_RNA_Matrix_MAP), ncol = ncol(TCGA_RNA_Matrix_MAP)) # Make the matrix a numeric
TCGA_RNA_Tibble_MAP <- as_tibble(TCGA_RNA_Matrix_MAP) # Make back into a tibble (now it will have numeric columns instead of chr)
colnames(TCGA_RNA_Tibble_MAP) <- MAPScore_Gene_Names # Add back column names
TCGA_RNA_Tibble_MAP <- TCGA_RNA_Tibble_MAP %>% add_column(SAMPLE_ID = TCGA_RNA_Tibble$SAMPLE_ID, .before = "ADAMTS14") # ADAMTS14 is the first gene in PoorOutcome. Adding Sample ID as a column before that

# cleaning up
rm(MAPScore_Genes, MAPScore_Genes_Filtered, TCGA_RNA_Matrix_MAP, TCGA_RNA_Tibble)

# Before calculating MAP Score in TCGA, I will remove all of the columns with NAs for sequencing data:
# This is 11 total and brings the TCGA MAP Gene list from 531 to 520 
TCGA_RNA_Tibble_MAP <- TCGA_RNA_Tibble_MAP %>% subset(select = -c(ANXA8L1,
                                                                  PRAMEF14,
                                                                  PRAMEF20,
                                                                  LYPD8,
                                                                  OR5K4,
                                                                  TXNDC8,
                                                                  XAGE5,
                                                                  PRAMEF18,
                                                                  OR8U3,
                                                                  OR6C1,
                                                                  DEFB121))

#### Calculate MAP score in TCGA cohort ####
# To do this I will do the sum of the Z-Scores divided by the number of genes
TCGA_RNA_Tibble_MAP <- TCGA_RNA_Tibble_MAP %>% rowwise() %>% mutate(MAP_Score = mean(c_across(where(is.numeric))))
head(TCGA_RNA_Tibble_MAP$MAP_Score)

# Create a new tibble with just the MAP Score and the sample ID
TCGA_Tibble_MAP <- TCGA_RNA_Tibble_MAP[c("SAMPLE_ID", "MAP_Score")]

# Save this new TCGA MAP score
write_csv(TCGA_Tibble_MAP, file = "data_in_use/22-0930_TCGA_MAP_Scores.csv")

# Cleaning up
rm(TCGA_Tibble_MAP, 
   TCGA_RNA_Tibble_MAP,
   Patient_Sample_Combined,
   MAPScore_Gene_Names)

