### Script Information: PCA of bulk RNA sequencing from all Samples excluding MNG, HT, PTCs with hashimoto's thyroiditis background
# Color by Diagnosis, BRS, PI3K Score, and TDS

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
# I don't know why there are 300 and not 299. One RNA.ID
x=count.df.master$Feature_gene_name

count.df.master2=count.df.master[,8:ncol(count.df.master)]
count.df.master2=as.matrix(count.df.master2)
row.names(count.df.master2) = x
# Note on this data: I don't know why there are 300 and not 299. 
# One RNA.ID is not present in our clinical data file: 6129CP23. 
# I don't know why this is, but I have emailed George today (22-0323) to ask

# Read in clinical Data: 
ClinicalData <- read_csv("data_in_use/VUMC.cohort.GX_4-26-22.csv")

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

# Make count.df.master2 and ClinicalData_Restricted contain the same samples 
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

# Can write results to CSV as follows: 
# write.csv(resOrdered, "resOrdered.csv)

# Make VSD from dds
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)

# Can write VSD to CSV as follows: 
# write.csv(assay(vsd), "x_vsd.csv")

# Perform PCA with our own new VSD file (instead of Tiger's)
pcDat <- pcDat <- prcomp(t(assay(vsd)))

sampleinfo <- Meta

# Plot PCA
autoplot(pcDat) + scale_y_reverse()

# Plot the PCA with the following annotations: Diagnosis, BRS, TDS, PI3K Score

# Diagnosis Simplified -> new annotation
sampleinfo$Diagnosis_Simplified <- sampleinfo$Diagnosis
for(i in 1:nrow(sampleinfo)){
  if(sampleinfo$Diagnosis_Simplified[i] == "FA" | sampleinfo$Diagnosis_Simplified[i] == "HA" | sampleinfo$Diagnosis_Simplified[i] == "FTC" | sampleinfo$Diagnosis_Simplified[i] == "HC" | sampleinfo$Diagnosis_Simplified[i] == "NIFTP"){
    sampleinfo$Diagnosis_Simplified[i] <- "Follicular"
  }
}


####


# Diagnosis Simplified (All Follicular grouped)
# Triangles on PDTC/ATC, outline on all, size 5 test
plot <- autoplot(pcDat,
                 data = sampleinfo,
                 fill = "Diagnosis_Simplified",
                 colour = "black",
                 shape = "Diagnosis_Simplified",
                 size = 5)
plot <- plot + scale_fill_manual(values = c(
  "royalblue1", # Follicular color
  "lightskyblue1", # FVPTC color
  "red", # PTC color
  "purple", # PDC color
  "magenta"), # ATC color - note that this is a new color
  "Histotype",
  breaks = c("Follicular", "FVPTC", "PTC", "PDTC", "ATC")) + 
  
  scale_shape_manual(values = c(21, 21, 21, 24, 24), "Histotype", 
                     breaks = c("Follicular", "FVPTC", "PTC", "PDTC", "ATC")) +
  
  scale_y_reverse() +
  # plot + labs(col = "Histotype") 
  theme_classic() +
  ggtitle("Diagnosis") +
  labs(x = "PC1", y = "PC2") +
  theme(#legend.box.background = element_rect(size = 2),
    legend.position = c(.16,.84),
    title = element_text(face = "bold", size = 30),
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(face = "bold", size = 30),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    legend.background = element_blank(),
    legend.text = element_text(face = "bold", size = 25),
    legend.title = element_blank()) 
ggsave(file = "outputs/22-0926_Fig2A_DESeq_All_Samples_PCA/22-0926_All_Samples_Diagnosis_Simplified_Triangles_Outline_Size5.png",
       width = 7,
       height = 6.5,
       plot, dpi = 600)

####

sampleinfo$BRS <- as.numeric(sampleinfo$BRS)
### BRAF-RAS Score w/ outline
plot <- autoplot(pcDat,
                 data = sampleinfo,
                 fill = "BRS",
                 shape = 21,
                 size = 5)
plot <- plot + scale_y_reverse()  + 
  ggtitle("BRAF-RAS Score") +
  labs(x = "PC1", y = "PC2") +
  scale_fill_gradient2(low = "red", mid = "grey", high = "blue", midpoint = 0, "BRAF-RAS \nScore", 
                       breaks = c(-.5,.5),
                       labels = c("BRAF", "RAS")) + 
  theme_classic() + 
  theme(#legend.box.background = element_rect(size = 2),
    legend.position = c(.12,.85),
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    title = element_text(face = "bold", size = 30), 
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(face = "bold", size = 30),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.text = element_text(face = "bold", size = 25),
    legend.title = element_blank())
#plot # Print plot
ggsave(file = "outputs/22-0926_Fig2A_DESeq_All_Samples_PCA/22-0926_All_Samples_BRAF-RAS_Score_Outline_size5.png",
       width = 7,
       height = 6.5,
       plot, dpi = 600)

### PI3K-AKT-MTOR Score outlined
sampleinfo$PI3K_AKT_MTOR <- as.numeric(sampleinfo$PI3K_AKT_MTOR)
plot <- autoplot(pcDat,
                 data = sampleinfo,
                 fill = "PI3K_AKT_MTOR",
                 shape = 21,
                 size = 5)
plot <- plot + scale_y_reverse()  + 
  scale_fill_gradient2(low = "green", mid = "grey", high = "red", midpoint = 0, "PI3K\nAKT\nMTOR\nScore",
                       breaks = c(-60, 60),
                       labels = c("Low", "High")) + 
  theme_classic() + 
  ggtitle("PI3K Activity Score") + 
  labs(x = "PC1", y = "PC2") +
  theme(#legend.box.background = element_rect(size = 2),
    legend.position = c(.12,.85),
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    title = element_text(face = "bold", size = 30), 
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(face = "bold", size = 30),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.text = element_text(face = "bold", size = 25),
    legend.title = element_blank())
#plot
ggsave(file = "outputs/22-0926_Fig2A_DESeq_All_Samples_PCA/22-0926_All_Samples_PI3K-AKT-MTOR_Score_Outline_size5.png",
       width = 7,
       height = 6.5,
       plot, dpi = 600)



### TDS Score w/ outline
sampleinfo$TDS <- as.numeric(sampleinfo$TDS)
plot <- autoplot(pcDat,
                 data = sampleinfo,
                 fill = "TDS",
                 shape = 21,
                 size = 5)
plot <- plot + scale_y_reverse()  + 
  scale_fill_gradient2(low = "yellow", mid = "orange", high = "brown", midpoint = 0, "TDS\nScore",
                       breaks = c(-.8,.8),
                       labels = c("Low", "High")) + 
  theme_classic() + 
  ggtitle("Thyroid Differentiation Score")+
  labs(x = "PC1", y = "PC2") +
  theme(#legend.box.background = element_rect(size = 2),
    legend.position = c(.12,.85),
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    title = element_text(face = "bold", size = 28), 
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(face = "bold", size = 30),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.text = element_text(face = "bold", size = 25),
    legend.title = element_blank())
#plot
ggsave(file = "outputs/22-0926_Fig2A_DESeq_All_Samples_PCA/22-0926_All_Samples_TDS_Score_Outline_size5.png",
       width = 7,
       height = 6.5,
       plot, dpi = 600)
