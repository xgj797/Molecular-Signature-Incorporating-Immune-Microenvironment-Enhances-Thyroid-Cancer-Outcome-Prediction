######### TDS_gene_heatmap_with_mets.R #########

###### Table of Contents ######

### I. Load Packages
### II. Differential Gene Expression Analysis and VSD Creation
### III. Heatmap of genes of differentiation score


###### I. Load Packages ######

# Load in packages
library(data.table)
library(FactoMineR)
library(calibrate)
library(ggfortify)
library(viridis)
library(paletteer)
library(scico)
library(ConsensusClusterPlus)
library(DESeq2)
library(EnhancedVolcano)
library(ComplexHeatmap)
library(circlize)



###### II. Differential Gene Expression Analysis and VSD Creation ######

# Set directory
setwd("_")

# Load protein coding count table
count.df.master=fread("weiss_human_wnt.proteincoding.count_diag_order_final_withmets_WITHOUT_HTPTCs_9-30-22.txt", header=TRUE, sep="\t")

# Change rownames to first column and convert to matrix
x=count.df.master$Feature_gene_name
# Subset, leaving out first 7 columns (non-sample columns)
count.df.master2=count.df.master[,8:ncol(count.df.master)]
# Convert to matrix
count.df.master2=as.matrix(count.df.master2)
# Add row names
row.names(count.df.master2)=x

# Read in metadatafile, edited order, note HTPTC samples NOT dropped
meta <- read.csv("meta_diag_and_CAF-EPIC_sort_thyroid_AND_met_remove_outlier_and_MNG_and_HT_merge_E_and_I_FVPTCs_WITHOUT_HTPTCs_9-30-22.csv", header=TRUE, row.names=1, as.is=TRUE)

# Restrict to only samples in meta
count.df <- count.df.master2[, colnames(count.df.master2) %in% rownames(meta)]

# Check meta is in correct order
test=colnames(count.df)==row.names(meta)
test

# Try varying "design", can try binary variables in meta such as sex, poor outcome.
# Can even add multiple ones together.
# The last variable of "design" is the one you're testing for, the ones prior are just being controlled for.
# For example, controlling for sex before testing for poor outcome, so you don't get sex related genes.
# Just the first variable in "design" has "~" preceding it
# Could control for thyroid and non-thyroid (met)

# We will control for sex for now
#dds <- DESeqDataSetFromMatrix(countData = count.df, colData = meta, design = ~Sex)

#Try testing:
dds <- DESeqDataSetFromMatrix(countData = count.df, colData = meta, design = ~Sex + Aggressive.disease)
#dds <- DESeqDataSetFromMatrix(countData = count.df, colData = meta, design = ~Sex + Diagnosis_Malig)
#dds <- DESeqDataSetFromMatrix(countData = count.df, colData = meta, design = ~Sex + Diagnosis_ATC)
#dds <- DESeqDataSetFromMatrix(countData = count.df, colData = meta, design = ~Sex + Diagnosis_transformed)

# Conduct differential gene expression analysis and create results dataframe
dds <- DESeq(dds)
res <- results(dds)
summary(res)
resOrdered <- res[order(res$padj),]
res_df <- as.data.frame(resOrdered)

#write.csv(resOrdered, "resOrdered_sex.csv")
write.csv(resOrdered, "resOrdered_sex_aggressivedisease.csv")
#write.csv(resOrdered, "resOrdered_sex_malig.csv")
#write.csv(resOrdered, "resOrdered_sex_ATC.csv")
#write.csv(resOrdered, "resOrdered_sex_transformed.csv")

# Make VSD from dds
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
#write.csv(assay(vsd), "sex.vsd.csv")
write.csv(assay(vsd), "sex_aggressivedisease.vsd.csv")
#write.csv(assay(vsd), "sex_malig.vsd.csv")
#write.csv(assay(vsd), "sex_ATC.vsd.csv")
#write.csv(assay(vsd), "sex_transformed.vsd.csv")


# Read in protein coding vsd normalized matrix
#thyroid=fread("sex.vsd.csv", header = TRUE,  sep = ",")
thyroid=fread("sex_aggressivedisease.vsd.csv", header = TRUE,  sep = ",")

# Set first column as rownames
x=thyroid$V1
thyroid2=thyroid[,2:ncol(thyroid)]
thyroid2=as.matrix(thyroid2)
row.names(thyroid2)=x


###### III. Heatmap of genes of differentiation score ######

#load ComplexHeatmap if not already loaded
#library(ComplexHeatmap)

#Data is "thyroid2" matrix
#Let's trying subsetting this matrix data to select genes of interest 

#first, define the subset
TDSgene <- c("DIO1", "DIO2", "DUOX1", "DUOX2", "FOXE1", "GLIS3", "NKX2-1", "PAX8",
             "SLC26A4", "SLC5A5", "SLC5A8", "TG", "THRA", "THRB", "TPO", "TSHR")
#then, extract these rows from the data
thyroid2TDSgene = thyroid2[rownames(thyroid2) %in% TDSgene, ]

#scale the data to Z-scores (by row)
mat <- t(scale(t(thyroid2TDSgene)))

# set heatmap colors
mypalette <- colorRamp2(c(-4, 0, 4), c("turquoise3", "white", "violetred3"))

# Top annotations
HA_Top <- HeatmapAnnotation(Diagnosis = meta$Diagnosis, # Note: Can add Diagnosis as a label instead of an annotation...gives the ability to write in the box
                            Location = meta$Location.type,
                            BRAF.mutation = meta$BRAF.mutation,
                            RAS.mutation = meta$RAS.mutation,
                            col = list(Diagnosis = c("ATC" = "magenta",  # These are the colors for adding Diagnosis as an annotation. Remove if adding as a label.
                                                     "PDTC" = "purple", 
                                                     "PTC or IFVPTC" = "red",
                                                     "FTC" = "blue",
                                                     "OTC" = "dodgerblue2",
                                                     "NIFTP or EFVPTC" = "lightskyblue",
                                                     "OA" = "cadetblue1",
                                                     "FA" = "lightcyan"),
                                       Location = c("Local disease"="Black",
                                                     "Regional LN met"="Green",
                                                     "Distant met"="Red"),
                                       BRAF.mutation = c("Yes" = "Red",
                                                           "No" = "Black",
                                                           "NA" = "Grey"),
                                       RAS.mutation = c("Yes" = "Blue",
                                                        "No" = "Black",
                                                        "NA"= "Grey")),
                            #gp = gpar(col = "black"), #draws black borders around annotation cells, if desired
                            annotation_name_gp = gpar(col = "black", fontsize = 10),
                            annotation_name_side = "right",
                            annotation_legend_param = list(Diagnosis=list(at=c("ATC","PDTC","PTC or IFVPTC","FTC","OTC","NIFTP or EFVPTC","OA","FA")),
                                                           Location=list(at=c("Local disease","Regional LN met", "Distant met")),
                                                           BRAF.mutation=list(at=c("Yes","No","NA")),
                                                           RAS.mutation=list(at=c("Yes","No","NA")))) #see https://jokergoo.github.io/ComplexHeatmap-reference/book/

# define heatmap (no column clustering and dendrogram)
thyroid2TDSgeneheatmap<-function(mat, meta) {
  Heatmap(mat,
          name = "mat", #name that appears above heatmap color scale legend
          col=mypalette, #set heatmap color to settings previously defined by mypalette
          column_order = 1:254, #prints out heatmap columns in the same order as meta
          row_names_gp = gpar(fontsize = 12), #row (gene) name font size
          column_names_gp = gpar(fontsize = 2), #column name (RNA.ID) font size
          top_annotation = HA_Top # Set the top annotation equal to HA_Top
  )
}

# make heatmap
png(filename="thyroid_TDS_genes_heatmap.png", #heatmap output filename
    res = 300,
    units = "in",
    pointsize = 4,
    width = 8,
    height = 6)
thyroid2TDSgeneheatmap(mat, meta)
dev.off()
