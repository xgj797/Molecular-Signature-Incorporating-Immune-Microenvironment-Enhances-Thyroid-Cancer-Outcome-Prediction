# Goal: Calculate an ERK Score for our thyroid cohort similar to the 2014 TCGA Paper ERK score calculation

# Reference: https://www.cell.com/cms/10.1016/j.cell.2014.09.050/attachment/5fa1f8ce-934c-4869-b3b7-b3f74526e829/mmc1.pdf
# ^ This reference is the 2014 TCGA Paper supplemental methods
# See page 22 - 24 for BRAF, TDS, and ERK score methods (and summary below)

#### TCGA Paper Method Notes:####
# Derived from 52 genes previously shown to be responsive to MEK inhibition (Pratilas et al., 2009)
# The 52 genes were differentially expressed in BRAFV600E mutant melanomas before AND after MEK inhibition (PD0325901)
# 48 downregulated with MEK inhibition (Set A)
# 4 up-regulated with MEK inhibition (Set B)
# ERK Activity signature was operationally defined by: 
# Quantification of the composite expression of this 52-gene set in each sample. 

#### Pratilas et al. method of quantification ####
# Z-score of the 52 genes in each sample 
# expression of each gene is first mean-centered and then divided by the corresponding standard deviation

#### TCGA Method of quantification ####
# 1. Log-transform the RNA-seq expression value of a given gene in sample i (Li = log transformed expression value)
# 3. Compute the 'pooled' mean and standard deviation
# 4. Compute the Z-score for sample i by subtracting the pooled mean from Li and dividing by the pooled standard deviation.
# 5. Once the Z-scores for all genes are computed, the ERK-activity score for sample i is the sum of the Z-scores for the genes in set A minus the Z-scores for the genes in set B

#### Original publication reference: ####
# https://www.pnas.org/content/106/11/4519.long
# https://www.pnas.org/content/pnas/suppl/2009/02/27/0900780106.DCSupplemental/0900780106SI.pdf (gene list is on page 10/11 of this PDF)
# Note that list is length 60 but that is because there are duplicate probe sets

#### My methods/modifications ####
# 9 Steps performed below

# Load packages
library(tidyverse)

#### Step 1: Read in & Subset Clinical Data to exclude MNG, HT (no need for ERK Score on these samples, would skew data) ####
# Read in clinical data
ClinicalData <- read_csv(file = "data_in_use/Clindata_Complete.csv")
ClinicalData <- ClinicalData[c("RNA.ID", "Diagnosis")] # Restrict clinical data to RNA Sequencing ID and Diagnosis
ClinicalData <- ClinicalData %>% subset(Diagnosis != "MNG" & Diagnosis != "HT") # Remove MNG, HT

###### Step 2: Read in Gene List and update Gene to current HGNC symbols to be compatible with our data ######
ERK_Genes <- read_csv(file = "data_in_use/22-0117_ERK_Gene_List.csv")

# Rename genes with new HGNC symbols
for(i in 1:nrow(ERK_Genes)){
  if(ERK_Genes$ERK_Genes[i] == "PYCRL"){ # Reference: https://www.genecards.org/cgi-bin/carddisp.pl?gene=PYCR3
    ERK_Genes$ERK_Genes[i] <- "PYCR3"
  }
  if(ERK_Genes$ERK_Genes[i] == "IL8"){ # https://www.genecards.org/cgi-bin/carddisp.pl?gene=CXCL8
    ERK_Genes$ERK_Genes[i] <- "CXCL8"
  }
  if(ERK_Genes$ERK_Genes[i] == "LNK"){ # Reference: https://www.genecards.org/cgi-bin/carddisp.pl?gene=SH2B3
    ERK_Genes$ERK_Genes[i] <- "SH2B3"
  }
  if(ERK_Genes$ERK_Genes[i] == "CD3EAP"){ # Reference: https://www.genecards.org/cgi-bin/carddisp.pl?gene=POLR1G
    ERK_Genes$ERK_Genes[i] <- "POLR1G"
  }
  if(ERK_Genes$ERK_Genes[i] == "HSPC111"){ # Reference: https://www.genecards.org/cgi-bin/carddisp.pl?gene=NOP16
    ERK_Genes$ERK_Genes[i] <- "NOP16"
  }
  if(ERK_Genes$ERK_Genes[i] == "BXDC2"){ # Reference: https://www.genecards.org/cgi-bin/carddisp.pl?gene=BRIX1
    ERK_Genes$ERK_Genes[i] <- "BRIX1"
  }
  if(ERK_Genes$ERK_Genes[i] == "FLJ10534"){ # Reference: https://www.genecards.org/cgi-bin/carddisp.pl?gene=TSR1
    ERK_Genes$ERK_Genes[i] <- "TSR1"
  }
  if(ERK_Genes$ERK_Genes[i] == "ALF"){ # Reference: https://www.genecards.org/cgi-bin/carddisp.pl?gene=GTF2A1L
    ERK_Genes$ERK_Genes[i] <- "GTF2A1L"
  }
}


###### Step 3: Subset Gene List to Set A (ERK upregulated/MEK inhibitor downregulated) and Set B (ERK downregulated/MEK inhibitor upregulated) ######
ERK_Genes_A <- ERK_Genes %>% subset (ERK_Upregulated == 1)
ERK_Genes_B <- ERK_Genes %>% subset (ERK_Upregulated == 0)

###### Step 4: Read in RNA Sequencing Data and Restrict to RNA.ID, ERK Genes, and Samples of interest ######
# Read in RNA sequencing data
Thyroid_Cohort_RNA <- read.table(file = "data_in_use/DESeq2_Normalized_RNA_Counts.txt", sep = '\t', header = TRUE)
Thyroid_Cohort_RNA <- as_tibble(Thyroid_Cohort_RNA) # Change type to tibble

## Restrict to samples and genes of interest
# Make cohort A that includes the 48 genes that go down with MEK inhibition
Thyroid_Cohort_RNA_A <- Thyroid_Cohort_RNA[c("RNA.ID", ERK_Genes_A$ERK_Genes)] # Restrict to ERK Genes + RNA IDs
Thyroid_Cohort_RNA_A <- Thyroid_Cohort_RNA_A %>% filter(RNA.ID %in% c(ClinicalData$RNA.ID)) # Restrict to Clinical Data Samples
Thyroid_Cohort_RNA_A_Gene_Set <- Thyroid_Cohort_RNA_A[c(ERK_Genes_A$ERK_Genes)] # Just genes, no IDs (optimal for math)

# Make cohort B that includes the 4 genes that go up with MEK inhibition 
Thyroid_Cohort_RNA_B <- Thyroid_Cohort_RNA[c("RNA.ID", ERK_Genes_B$ERK_Genes)] # Restrict to ERK genes + RNA IDs
Thyroid_Cohort_RNA_B <- Thyroid_Cohort_RNA_B %>% filter(RNA.ID %in% c(ClinicalData$RNA.ID)) # Restrict to Clinical Data Samples
Thyroid_Cohort_RNA_B_Gene_Set <- Thyroid_Cohort_RNA_B[c(ERK_Genes_B$ERK_Genes)] # Just genes, no IDs (optimal for math)

###### Step 5: Log2 transform the data ######
Thyroid_Cohort_RNA_A_Log2 <- Thyroid_Cohort_RNA_A_Gene_Set + 1 # Add 1 to avoid log2(0)
Thyroid_Cohort_RNA_B_Log2 <- Thyroid_Cohort_RNA_B_Gene_Set + 1 # Add 1 to avoid log2(0)
Thyroid_Cohort_RNA_A_Log2 <- log(Thyroid_Cohort_RNA_A_Log2, 2) # log2 transform the data
Thyroid_Cohort_RNA_B_Log2 <- log(Thyroid_Cohort_RNA_B_Log2, 2) # log2 transform the data

###### Step 6: Calculate Pooled Means and Standard Deviations ######
# Cohort A means and STDevs
Cohort_A_Column_Means <- colMeans(Thyroid_Cohort_RNA_A_Log2) # Takes the mean for each gene across the cohort and makes it into a vector of all mean values
Cohort_A_STDev <- apply(Thyroid_Cohort_RNA_A_Log2, 2, sd) # Takes the STDev for each gene across the cohort and makes it into a vector of all STDev values

# Cohort B means and STDevs
Cohort_B_Column_Means <- colMeans(Thyroid_Cohort_RNA_B_Log2) # Takes the mean for each gene across the cohort and makes it into a vector of all mean values
Cohort_B_STDev <- apply(Thyroid_Cohort_RNA_B_Log2, 2, sd) # Takes the STDev for each gene across the cohort and makes it into a vector of all STDev values

###### Step 7: Use Pooled Means and Standard Deviations to calculate Z-Scores ######
# Cohort A Z-Scores
Thyroid_Cohort_A_Z_Scores <- as.matrix(Thyroid_Cohort_RNA_A_Log2) # Make a Z-scores matrix to put calculated values in
for(y in 1:ncol(Thyroid_Cohort_RNA_A_Log2)){
  for (x in 1:nrow(Thyroid_Cohort_RNA_A_Log2)){
    Thyroid_Cohort_A_Z_Scores[x,y] <- (Thyroid_Cohort_RNA_A_Log2[x,y]-Cohort_A_Column_Means[y])/Cohort_A_STDev[y] # (GeneX[i] - GeneX[mean])/GeneX[STDev])
  }
}

# Cohort B Z-Scores
Thyroid_Cohort_B_Z_Scores <- as.matrix(Thyroid_Cohort_RNA_B_Log2) # Makes a Z-score matrix to put calculated values in
for(y in 1:ncol(Thyroid_Cohort_RNA_B_Log2)){
  for (x in 1:nrow(Thyroid_Cohort_RNA_B_Log2)){
    Thyroid_Cohort_B_Z_Scores[x,y] <- (Thyroid_Cohort_RNA_B_Log2[x,y]-Cohort_B_Column_Means[y])/Cohort_B_STDev[y] # (GeneX[i] - GeneX[mean])/GeneX[STDev])
  }
}

###### Step 8: sum the Z-scores for each sample #######
Cohort_A_Sums <- apply(Thyroid_Cohort_A_Z_Scores, 1, sum) # Sum Z-scores in Cohort A genes
Cohort_B_Sums <- apply(Thyroid_Cohort_B_Z_Scores, 1, sum) # Sum Z-scores in Cohort B genes
Cohort_Sums <- Cohort_A_Sums - Cohort_B_Sums # Subtract the Cohort B sums from the Cohort A sums to get ERK Scores

###### Step 9: Combine ERK Score with RNA.ID and Output as CSV to data_in_use ######
ERK_Score_Output <- data.frame(Thyroid_Cohort_RNA_A$RNA.ID, Cohort_Sums) # Combine RNA ID and Cohort Sums (ERK Scores)
Test <- data.frame(Thyroid_Cohort_RNA_A$RNA.ID, Cohort_A_Sums, Cohort_B_Sums, Cohort_Sums) # Just looking at math  to make sure vector subtraction worked
write.csv(ERK_Score_Output, "data_in_use/22-0131_ERK_Score_by_TCGA_Method.csv", row.names = FALSE) # Output as CSV


