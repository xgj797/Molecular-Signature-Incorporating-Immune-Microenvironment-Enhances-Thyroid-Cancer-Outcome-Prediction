# DESEQ for all samples (non-metastatic) based on BRAF-Like (BRS < 0) vs RAS-Like (BRS > 0) determination 

### Load packages
library(data.table)
library(DESeq2)
library(FactoMineR)
library(calibrate)
library(ggfortify)
library(viridis)
library(paletteer)
library(scico)
library(tidyverse)

# Load protein coding count table
count.df.master=fread("data_in_use/21-0823_weiss_human_wnt.proteincoding.count", header = TRUE, sep = "\t")

# Change rownames to first column and convert to matrix
x=count.df.master$Feature_gene_name

count.df.master2=count.df.master[,8:ncol(count.df.master)]
count.df.master2=as.matrix(count.df.master2)
row.names(count.df.master2) = x

# Read in clinical Data: 
# Note that this most recent file does NOT include the most recent BRS, will have to load that as well
ClinicalData <- read_csv("data_in_use/VUMC.cohort.GX_9-8-22_V2.csv")

# Subset Clinical Data -> new file with all samples excluding MNG, HT, and HTPTC
ClinicalData_Restricted <- ClinicalData %>% subset(Diagnosis != "MNG" & 
                                                   Diagnosis != "HT" & 
                                                   Diagnosis != "HTPTC" &
                                                   Diagnosis != "normal" &
                                                   !(is.na(BRS)))

# Add an X to the RNA.ID column to be compatible with the gene expression data
ClinicalData_Restricted$RNA.ID <- sub('', 'X', ClinicalData_Restricted$RNA.ID)

ClinicalData_Restricted$BRS_Status <- "BRAF_Like"
for(i in 1:nrow(ClinicalData_Restricted)){
  if(ClinicalData_Restricted$BRS[i] > 0){
    ClinicalData_Restricted$BRS_Status[i] <- "RAS_Like"
  }
}


# Convert to data frame and add row names that match with count.df.master2
ClinicalData_Restricted <- as.data.frame(ClinicalData_Restricted)
rownames(ClinicalData_Restricted) <- ClinicalData_Restricted$RNA.ID

# Make count.df.master2 and ClinicalData_Restricted contain same samples
count.df <- count.df.master2[,colnames(count.df.master2) %in% rownames(ClinicalData_Restricted)] # Restrict count.df to samples in ClinicalData
Meta <- ClinicalData_Restricted[rownames(ClinicalData_Restricted) %in% colnames(count.df),] # Restrict meta to ClinicalData_Local_ATCs in count.df

# Test order
test <- colnames(count.df) == rownames(Meta)
test

# Out of order, need to re-order
count.df <- count.df[, rownames(Meta)]

# Test order take 2
test <- colnames(count.df) == rownames(Meta)
test

# Create DESeq2 object with count.df, Meta, as a function of BRS_Status
dds <- DESeqDataSetFromMatrix(countData = count.df, colData = Meta, design = ~BRS_Status)
dds <- DESeq(dds)
res <- results(dds)
summary(res)
resOrdered <- res[order(res$padj), ]
res_df <- as.data.frame(resOrdered)

saveRDS(res_df, file = "data_in_use/22-0908_DESEQ2_All_Primary_LocalDisease_MG_HT_Excluded_resOrdered_BRAF_RAS.rds")
rm(list = ls())

