
### This script: Normalization of count data with DESeq2

# For DESeq installation, see the following: 
# https://bioconductor.org/packages/release/bioc/html/DESeq2.html

### Load necessary packages:
library(data.table)
library(tidyverse)
library(DESeq2)

# Load protein coding count table
count.df.master=fread("data_in_use/21-0823_weiss_human_wnt.proteincoding.count", header = TRUE, sep = "\t")

# Change rownames to first column and convert to matrix
x=count.df.master$Feature_gene_name
count.df.master2=count.df.master[,8:ncol(count.df.master)]
count.df.master2=as.matrix(count.df.master2)
row.names(count.df.master2) = x

# Read in metadatafile
meta <- read.csv("data_in_use/Clindata_Complete.csv", row.names = 1, as.is = TRUE)

# Restric to only samples in meta
count.df <- count.df.master2[,colnames(count.df.master2) %in% rownames(meta)]

# Restric meta to only samples in sequencing data
meta2 <- meta[rownames(meta) %in% colnames(count.df),]

# Check meta is in correct order
test = colnames(count.df) == row.names(meta2)
test # returned TRUE

# Creat a column "Prim"
meta2$Prim <- "Primary"
for(i in 1:nrow(meta2)){
  if (meta2$Primary[i] == 0){
    meta2$Prim[i] <- "Met"
  }
}

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = count.df, colData = meta2, design = ~Prim)

dds <- estimateSizeFactors(dds)
sizeFactors(dds)

# Create normalized countrs matrix
normalized_counts <- counts(dds, normalized = TRUE)

# Transpose normalized matrix: 
normalized_counts_transposed <- t(normalized_counts)

# Convert normalized counts to tibble
normalized_counts_tibble <- as_tibble(normalized_counts_transposed, rownames = "RNA.ID")

# Save normalized as text file
write.table(normalized_counts_tibble,"data_in_use/DESeq2_Normalized_RNA_Counts.txt",sep="\t",row.names=FALSE)

# Make a variance stabilising transformation 
# E.g., for heat maps, Z-scoring, or plotting differential gene expression analysis

VSD <- varianceStabilizingTransformation(dds, blind=TRUE)

# Make the VSD into a readable matrix
VSD_Output <- assay(VSD)

# From the internet regarding calculating Z-scores from VSDs: 
# We apply scale which is the actual Z-scoring function
# In R data.frames and matrices are basically lists of vectors (every column is a vector) but we want rowwise (per gene) Z-scores.
# Scale by default (like any operation on matrices/dfs) operates column-wise though, so we first t (transpose) the matrix, then scale it, and then 
# t transpose it back so we again have a column = sample and row = gene matrix (or data.frame). That's it.


# Save Z scores as a matrix
# Note, I'm not doing the second transpose because I want my genes as columns
Z <- scale(t(VSD_Output))

# Save Z scores as a tibble with RNA.ID labeling the first column
Z_Scores_Output <- as_tibble(Z, rownames = "RNA.ID")

# Save Z-scores as text file
write.table(Z_Scores_Output,"data_in_use/DESeq2_Z-Scores.txt",sep="\t",row.names=FALSE)
