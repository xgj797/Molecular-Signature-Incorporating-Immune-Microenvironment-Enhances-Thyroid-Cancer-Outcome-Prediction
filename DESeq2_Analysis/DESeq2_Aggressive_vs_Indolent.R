# DESEQ for all samples (non-metastatic) based on aggressive/indolent determination 

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

# Read in clinical Data
ClinicalData <- read_csv("data_in_use/VUMC.cohort.GX_9-8-22_V2.csv")

# Subset Clinical Data -> new file with just samples of interest
ClinicalData_Local <- ClinicalData %>% subset(Diagnosis.with.HTPTC != "MNG" & 
                                              Diagnosis.with.HTPTC != "HT" &
                                              Diagnosis.with.HTPTC != "HTPTC" & # removing HTPTCs due to backgrount immune signal
                                              (Location.type == "NonMet")) 
# Note: numbers at this step break down to the same as prior (193), good to proceed

# Subset clinical data local to just aggressive vs indolent
ClinicalData_Local <- ClinicalData_Local %>% subset(Aggression.category.for.Matt.Deseq == 1 |
                                                    Aggression.category.for.Matt.Deseq == 0)
# This leads to 187 samples, which is the same as prior

# Giving names "A" = poor outcome and "B" = not poor outcome
# This will ensure that group "A" is the positive log-fold change (poor outcome)
# And group "B" is the negative log-fold change (poor outcome)
ClinicalData_Local$DESEQ_Group <- "B"
for(i in 1:nrow(ClinicalData_Local)){
  if(ClinicalData_Local$Aggression.category.for.Matt.Deseq[i] == 1){
    ClinicalData_Local$DESEQ_Group[i] <- "A"
  }
}

# Convert to data frame and add row names that match with count.df.master2
ClinicalData_Local <- as.data.frame(ClinicalData_Local)
rownames(ClinicalData_Local) <- ClinicalData_Local$RNA.ID

# Make count.df.master2 and ClinicalData_Local contain the same samples 
# need to rename the colnames of count.df.master2 to remove the x in front of the RNA ID
colnames(count.df.master2) <- sub('.', '', colnames(count.df.master2))
count.df <- count.df.master2[,colnames(count.df.master2) %in% rownames(ClinicalData_Local)] # Restrict count.df to samples in ClinicalData_Local
Meta <- ClinicalData_Local[rownames(ClinicalData_Local) %in% colnames(count.df),] # Restrict meta to ClinicalData_Local_ATCs in count.df

# Test order
test <- colnames(count.df) == rownames(Meta)
test

# Out of order, need to re-order
count.df <- count.df[, rownames(Meta)]

# Test order take 2
test <- colnames(count.df) == rownames(Meta)
test
# Note that Meta has the same number of observations (178) as prior 
# This is exciting!

# Create DESeq2 object with count.df, Meta, as a function of BRS_Status
dds <- DESeqDataSetFromMatrix(countData = count.df, colData = Meta, design = ~DESEQ_Group)
dds <- DESeq(dds)
res <- results(dds)
# If you print results(dds), it shows Group B vs A, indicating that group B is the positive log2FoldChange
summary(res)
resOrdered <- res[order(res$padj), ]
res_df <- as.data.frame(resOrdered)

# Note: Negative log-fold change = aggressive (group A)
# I now see this consistantly where the negative log2fold change is the group that is alphabetically first
# Not sure why this happened, but I'm going to fix it as follows: 
res_df$log2FoldChange <- res_df$log2FoldChange*-1

# Now I have them ordered so that positive = up in poor outcome cases

# Will save data frame as RDS that can be easily loaded into other R scripts for future use
saveRDS(res_df, file = "data_in_use/22-0908_DESEQ2_All_Primary_LocalDisease_MG_HT_Excluded_resOrdered_Aggressive_Indolent.rds")
rm(list = ls())
# example loading
All_DESEQ_Aggressive <- readRDS(file = "data_in_use/22-0908_DESEQ2_All_Primary_LocalDisease_MG_HT_Excluded_resOrdered_Aggressive_Indolent.rds")
# Example cleaning
rm(list = ls())
