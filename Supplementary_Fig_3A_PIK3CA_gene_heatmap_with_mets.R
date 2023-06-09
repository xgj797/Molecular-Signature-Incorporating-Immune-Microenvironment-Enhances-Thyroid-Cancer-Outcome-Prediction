######### PIK3CA_gene_heatmap_with_mets.R #########

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
count.df.master=fread("weiss_human_wnt.proteincoding.count_diag_order_final_withmets_and_HTPTC_9-30-22.txt", header=TRUE, sep="\t")

# Change rownames to first column and convert to matrix
x=count.df.master$Feature_gene_name
# Subset, leaving out first 7 columns (non-sample columns)
count.df.master2=count.df.master[,8:ncol(count.df.master)]
# Convert to matrix
count.df.master2=as.matrix(count.df.master2)
# Add row names
row.names(count.df.master2)=x

# Read in metadatafile, edited order, note HTPTC samples NOT dropped
meta <- read.csv("meta_diag_and_CAF-EPIC_sort_thyroid_AND_met_remove_outlier_and_MNG_and_HT_merge_E_and_I_FVPTCs_include_HTPTCs_9-30-22.csv", header=TRUE, row.names=1, as.is=TRUE)

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
PIK3CAgene <- c("ACACA",
              "ACTR2",
              "ACTR3",
              "ADCY2",
              "AKT1",
              "AKT1S1",
              "AP2M1",
              "ARF1",
              "ARHGDIA",
              "ARPC3",
              "ATF1",
              "CAB39",
              "CAB39L",
              "CALR",
              "CAMK4",
              "CDK1",
              "CDK2",
              "CDK4",
              "CDKN1A",
              "CDKN1B",
              "CFL1",
              "CLTC",
              "CSNK2B",
              "CXCR4",
              "DAPP1",
              "DDIT3",
              "DUSP3",
              "E2F1",
              "ECSIT",
              "EGFR",
              "EIF4E",
              "FASLG",
              "FGF17",
              "FGF22",
              "FGF6",
              "GNA14",
              "GNGT1",
              "GRB2",
              "GRK2",
              "GSK3B",
              "HRAS",
              "HSP90B1",
              "IL2RG",
              "IL4",
              "IRAK4",
              "ITPR2",
              "LCK",
              "MAP2K3",
              "MAP2K6",
              "MAP3K7",
              "MAPK1",
              "MAPK10",
              "MAPK8",
              "MAPK9",
              "MAPKAP1",
              "MKNK1",
              "MKNK2",
              "MYD88",
              "NCK1",
              "NFKBIB",
              "NGF",
              "NOD1",
              "PAK4",
              "PDK1",
              "PFN1",
              "PIK3R3",
              "PIKFYVE",
              "PIN1",
              "PITX2",
              "PLA2G12A",
              "PLCB1",
              "PLCG1",
              "PPP1CA",
              "PPP2R1B",
              "PRKAA2",
              "PRKAG1",
              "PRKAR2A",
              "PRKCB",
              "PTEN",
              "PTPN11",
              "RAC1",
              "RAF1",
              "RALB",
              "RIPK1",
              "RIT1",
              "RPS6KA1",
              "RPS6KA3",
              "RPTOR",
              "SFN",
              "SLA",
              "SLC2A1",
              "SMAD2",
              "SQSTM1",
              "STAT2",
              "TBK1",
              "THEM4",
              "TIAM1",
              "TNFRSF1A",
              "TRAF2",
              "TRIB3",
              "TSC2",
              "UBE2D3",
              "UBE2N",
              "VAV3",
              "YWHAB")
#then, extract these rows from the data
thyroid2PIK3CAgene = thyroid2[rownames(thyroid2) %in% PIK3CAgene, ]

#scale the data to Z-scores (by row)
mat <- t(scale(t(thyroid2PIK3CAgene)))

# set heatmap colors
mypalette <- colorRamp2(c(-4, 0, 4), c("green", "white", "red"))

# Top annotations
HA_Top <- HeatmapAnnotation(Diagnosis = meta$Diagnosis, # Note: Can add Diagnosis as a label instead of an annotation...gives the ability to write in the box
                            Location = meta$Location.type,
                            BRAF.mutation = meta$BRAF.mutation,
                            RAS.mutation = meta$RAS.mutation,
                            col = list(Diagnosis = c("ATC" = "magenta",  # These are the colors for adding Diagnosis as an annotation. Remove if adding as a label.
                                                     "PDTC" = "purple", 
                                                     "PTC or IFVPTC" = "red",
                                                     "FTC" = "blue",
                                                     "HC" = "dodgerblue2",
                                                     "NIFTP or EFVPTC" = "lightskyblue",
                                                     "HA" = "cadetblue1",
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
                            annotation_legend_param = list(Diagnosis=list(at=c("ATC","PDTC","PTC or IFVPTC","FTC","HC","NIFTP or EFVPTC","HA","FA")),
                                                           Location=list(at=c("Local disease","Regional LN met", "Dist metastasis")),
                                                           BRAF.mutation=list(at=c("Yes","No","NA")),
                                                           RAS.mutation=list(at=c("Yes","No","NA")))) #see https://jokergoo.github.io/ComplexHeatmap-reference/book/

# define heatmap (no column clustering and dendrogram)
thyroid2PIK3CAgeneheatmap<-function(mat, meta) {
  Heatmap(mat,
          name = "mat", #name that appears above heatmap color scale legend
          col=mypalette, #set heatmap color to settings previously defined by mypalette
          column_order = 1:265, #prints out heatmap columns in the same order as meta
          row_names_gp = gpar(fontsize = 4), #row (gene) name font size
          column_names_gp = gpar(fontsize = 2), #column name (RNA.ID) font size
          top_annotation = HA_Top # Set the top annotation equal to HA_Top
  )
}

# make heatmap
png(filename="thyroid_PIK3CA_genes_heatmap.png", #heatmap output filename
    res = 300,
    units = "in",
    pointsize = 4,
    width = 8,
    height = 6)
thyroid2PIK3CAgeneheatmap(mat, meta)
dev.off()
