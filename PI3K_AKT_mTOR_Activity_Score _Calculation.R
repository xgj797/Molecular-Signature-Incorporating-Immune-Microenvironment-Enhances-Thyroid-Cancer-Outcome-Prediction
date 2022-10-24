### PI3K Score: Information:
# This script reads in DESeq2 normalized expression data for our thyroid cohort and a 105 gene list for PI3K-AKT-mTOR signaling activity and uses the gene list to create a PI3K score for each sample (HT, MNG excluded).
# The gene list is from HALLMARK_PI3K_AKT_MTOR_SIGNALING

#### Information extended: Method of quantification ####
# 1. Log-transform the RNA-seq expression value of a given gene in sample i (Li = log transformed expression value)
# 2. Compute the 'pooled' median and standard deviation of the gene (across all samples)
# 4. Compute the Z-score for sample i by subtracting the pooled median from Li and dividing by the pooled standard deviation.
# 5. Once the Z-scores for all genes are computed, the PIK3-AKT-mTOR Activity score for sample i is the sum of the Z-scores for the genes in the gene set

#### Load packages ####
library(tidyverse)

###### Step 1: Read in & Subset Clinical Data to exclude MNG and HT (will not calculate PI3K activity score on these samples######
# Read in clinical data
ClinicalData <- read_csv(file = "data_in_use/Clindata_Complete.csv")
ClinicalData <- ClinicalData[c("RNA.ID", "Diagnosis")] # Restrict clinical data to RNA Sequencing ID and Diagnosis
ClinicalData <- ClinicalData %>% subset(Diagnosis != "MNG" & Diagnosis != "HT") # Remove MNG, HT

###### Step 2: Read in Gene List ######
PI3K_AKT_MTOR <- read.table(file = "data_in_use/22-0314_HALLMARK_PI3K_AKT_MTOR_SIGNALING.txt", sep = "\t")
PI3K_AKT_MTOR <- as_tibble(PI3K_AKT_MTOR)
colnames(PI3K_AKT_MTOR) <- PI3K_AKT_MTOR[1,] # Set col names = to first row
PI3K_AKT_MTOR <- PI3K_AKT_MTOR[3:nrow(PI3K_AKT_MTOR),] # Start list with first column as first gene (remove first two rows)

###### Step 4: Read in RNA Sequencing Data and Restrict to RNA.ID, PI3K Genes, and Samples of interest ######
# Read in RNA sequencing data
Thyroid_Cohort_RNA <- read.table(file = "data_in_use/DESeq2_Normalized_RNA_Counts.txt", sep = '\t', header = TRUE)
Thyroid_Cohort_RNA <- as_tibble(Thyroid_Cohort_RNA) # Change type to tibble

## Restrict to samples and genes of interest
# Make cohort A that includes the 48 genes that go down with MEK inhibition
Thyroid_Cohort_RNA_Restricted <- Thyroid_Cohort_RNA[c("RNA.ID", PI3K_AKT_MTOR$HALLMARK_PI3K_AKT_MTOR_SIGNALING)] # Restrict to WNT Genes + RNA IDs
Thyroid_Cohort_RNA_Restricted <- Thyroid_Cohort_RNA_Restricted %>% filter(RNA.ID %in% c(ClinicalData$RNA.ID)) # Restrict to Clinical Data Samples
Thyroid_Cohort_RNA_Gene_Set <- Thyroid_Cohort_RNA_Restricted[c(PI3K_AKT_MTOR$HALLMARK_PI3K_AKT_MTOR_SIGNALING)] # Just genes, no IDs (optimal for math)

###### Step 5: Log2 transform the data ######
Thyroid_Cohort_RNA_Log2 <- Thyroid_Cohort_RNA_Gene_Set + 1 # Add 1 to avoid log2(0)
Thyroid_Cohort_RNA_Log2 <- log(Thyroid_Cohort_RNA_Log2, 2) # log2 transform the data

###### Step 6: Calculate Pooled Medians and Standard Deviations ######
# Cohort A medians and STDevs
Cohort_Column_Medians <- apply(Thyroid_Cohort_RNA_Log2, 2, median) # Takes the median for each gene across the cohort and makes it into a vector for all medians
Cohort_STDev <- apply(Thyroid_Cohort_RNA_Log2, 2, sd) # Takes the STDev for each gene across the cohort and makes it into a vector of all STDev values

###### Step 7: Use Pooled Medians and Standard Deviations to calculate Z-Scores ######
# Cohort Z-Scores
Thyroid_Cohort_Z_Scores <- as.matrix(Thyroid_Cohort_RNA_Log2) # Make a Z-scores matrix to put calculated values in
for(y in 1:ncol(Thyroid_Cohort_RNA_Log2)){
  for (x in 1:nrow(Thyroid_Cohort_RNA_Log2)){
    Thyroid_Cohort_Z_Scores[x,y] <- (Thyroid_Cohort_RNA_Log2[x,y]-Cohort_Column_Medians[y])/Cohort_STDev[y] # (GeneX[i] - GeneX[mean])/GeneX[STDev])
  }
}

###### Step 8: sum the Z-scores for each sample #######
Cohort_Sums <- apply(Thyroid_Cohort_Z_Scores, 1, sum) # Sum Z-scores in Cohort A genes
Cohort_Sums_Average <- Cohort_Sums/105 # Divide by # of genes in gene set -> average Z-score

###### Step 9: Combine PI3K Score with RNA.ID and Output as CSV to data_in_use ######
PI3K_AKT_MTOR_Output <- data.frame(Thyroid_Cohort_RNA_Restricted$RNA.ID, Cohort_Sums) # Combine RNA ID and Cohort Sums (PI3K Scores)
PI3K_AKT_MTOR_Output <- PI3K_AKT_MTOR_Output %>% dplyr::rename(RNA.ID = Thyroid_Cohort_RNA_Restricted.RNA.ID)
write.csv(PI3K_AKT_MTOR_Output, "data_in_use/22-0314_HALLMARK_PI3K_AKT_MTOR_SIGNALING_Score_by_TCGA_Method.csv", row.names = FALSE) # Output as CSV

