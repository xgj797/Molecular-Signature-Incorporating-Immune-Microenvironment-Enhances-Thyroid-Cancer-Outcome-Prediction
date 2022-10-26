# Script: Oncoplot.R

#### Guide for ComplexHeatmap ####
# https://jokergoo.github.io/ComplexHeatmap-reference/book/
# See chapter 7 for oncoPrint instructions

#### Table of Contents ####
# 1. Load packages
# 2. Read maf file
# 3. Change clinical data of interest to numerics -> this allows for them to be plotted as numerics
# 4. Subset maf to samples of interest
# 5. Use maf tools to print a matrix to feed into ComplexHeatmap's oncoprint function
# 6. Read in Oncoplot Matrix generated with maftools
# 7. Format and subset Oncoplot Matrix and Clinical Data
# 8. Add in fusion data to Oncoplot Data
# 9. Make Diagnosis into a numerical variable for sorting
# 10. Formatting mutation grid of ComplexHeatmap oncoPrint plot
# 11. Formatting ComplexHeatmap annotations
# 12. Make the ComplexHeatmap with oncoPrint

#ComplexHeatmap install may require XQuartz installation

#### 1. Load Packages ####
library(tidyverse)
library(maftools)
library(ComplexHeatmap)
library(circlize) # R package for making gradient heatmap annotations

# Set your working directory
setwd("_")

#### 2. Read in Data ####
# read maf file
maf <- read.maf(maf = "listgenes_vaclass_diagnosis_2_24_22_update.maf")

# read clin dataframe
Clinical_Data <- read.csv("Diagnosis_guide_9-28-22.csv", header = TRUE)
Clinical_Data <- as_tibble(Clinical_Data)

# remove normals -> no resection data (e.g., sequencing, etc.)
Clinical_Data <- Clinical_Data %>% subset(Diagnosis != "normal")

#### 3. Change clinical data to numerics ####
### Some of the variables need to be numerics to generate gradients in the oncoplot ###
# Change age at resection to numeric
Age.resection.integer <- as.numeric(Clinical_Data$Age.resection)
Clinical_Data$Age.resection <- Age.resection.integer

# Change ERK Score variable to numeric
ERK_Score_numeric <- as.numeric(Clinical_Data$ERK)
Clinical_Data$ERK <- ERK_Score_numeric

# Change BRS variable as numeric
BRS_numeric <- as.numeric(Clinical_Data$BRS)
Clinical_Data$BRS <- BRS_numeric

# Change TDS variable to numeric
TDS_numeric <- as.numeric(Clinical_Data$TDS)
Clinical_Data$TDS <- TDS_numeric

# Change PI3K_AKT_MTOR variable to numeric
PI3K_AKT_MTOR_numeric <- as.numeric(Clinical_Data$PI3K_AKT_MTOR)
Clinical_Data$PI3K_AKT_MTOR <- PI3K_AKT_MTOR_numeric

# Needs to be in vector format for oncoplot
Clinical_Data <- as.data.frame(Clinical_Data)

#### 4. Subset maf to samples of interest ####

#Diagnoses of interest
Clinical_Data <- Clinical_Data %>% subset(Diagnosis == "ATC" | Diagnosis == "PDTC" | Diagnosis == "PTC and IFVPTC" | Diagnosis == "NIFTP and EFVPTC" | Diagnosis == "FTC" | Diagnosis == "HC")

maf_thyroidorigin = subsetMaf(maf = maf, tsb = c(Clinical_Data$Tumor_Sample_Barcode), mafObj = TRUE)

#### 5. Use MAF tools to print a Matrix to feed into ComplexHeatmap's oncoPrint function ####
oncoplot(maf = maf_thyroidorigin, top = 20, annotationDat = Clinical_Data, # top = set the number of genes to look at, annotationDat = data for annotations
         clinicalFeatures = c("Diagnosis", # Annotation features to show at bottom
                              "Location.type", 
                              "Sex",
                              "Age.resection",
                              "Aggressive.disease"),
         anno_height = 2,
         removeNonMutated = FALSE, # Shows "Dark Matter"
         sortByAnnotation = TRUE, # Must sort by annotation of Diagnosis ... because diagnosis is listed first it is automatically used
         annotationOrder = c("ATC","PDTC","PTC and IFVPTC", "NIFTP and EFVPTC", "FTC","HC"), # Here is the order that I want the diagnoses in 
         numericAnnoCol = NULL, # This is for changing the color of annotation gradients...I couldn't get it to work in maf tools though
         writeMatrix = TRUE) # This is the most important step -> I will use writeMatrix as input for the complex heatmap package

#### 6. Read in Oncoplot Matrix generated with maftools ####
Oncoplot_Data <- as.matrix(read.table(file = "onco_matrix.txt", sep = '\t', header = TRUE))

#### 7. Format and subset Oncoplot Matrix and Clinical Data ####
Sample_Names <- colnames(Oncoplot_Data) # Create a variable "Sample_Names" that contains a list of the samples in the OncoPlot Matrix
Clinical_Data_Oncoplot_Restricted <- Clinical_Data %>% subset(Tumor_Sample_Barcode %in% c(Sample_Names)) # restrict clinical data to the samples in the OncoPlot matrix
rownames(Clinical_Data_Oncoplot_Restricted) <- Clinical_Data_Oncoplot_Restricted$Tumor_Sample_Barcode # Assign row names as the tumor sample barcode
Clinical_Data_Oncoplot_Restricted <- Clinical_Data_Oncoplot_Restricted[colnames(Oncoplot_Data),] # this reorders the clinical data so that it is in the same order as the OncoPlot data...they must match for OncoPrint to work
Clinical_Data_Oncoplot_Restricted # print Clinical Data for the restricted oncoplot -> just for myself to visualize it
test <- rownames(Clinical_Data_Oncoplot_Restricted) == colnames(Oncoplot_Data) # must be same order for annotations to work
test # print out test -> should be all "true"

# Replace "0" with "" in Oncoplot data -> this is the format that is read by ComplexHeatmap
for(y in 1:nrow(Oncoplot_Data)){
  for(x in 1:ncol(Oncoplot_Data)){
    if(Oncoplot_Data[y,x] == "0"){
      Oncoplot_Data[y,x] <- ""
    }
  }
}

#### 8. Add in fusion data to Oncoplot_Data ####
# Note: can do similar for CNV data if desired
# MET fusions
Oncoplot_Data_Fusions_Included <- Oncoplot_Data # make a new matrix of Oncoplot_Data with fusions_included designated
MET_fusion <- c(Clinical_Data_Oncoplot_Restricted$MET_fusion) # Pull MET fusion data from Clinical Data
names(MET_fusion) <- names(Oncoplot_Data_Fusions_Included) # assign sample names to MET fusion data
Oncoplot_Data_Fusions_Included <- rbind(Oncoplot_Data_Fusions_Included, MET_fusion) # bind MET fusion data onto OncoPlot Data

# BRAF fusions - the rest of the fusions will be a repeat of what was done for MET
BRAF_fusion <- c(Clinical_Data_Oncoplot_Restricted$BRAF_fusion)
names(BRAF_fusion) <- names(Oncoplot_Data_Fusions_Included)
Oncoplot_Data_Fusions_Included <- rbind(Oncoplot_Data_Fusions_Included, BRAF_fusion)

# PPARG fusions
PPARG_fusion <- c(Clinical_Data_Oncoplot_Restricted$PPARG_fusion)
names(PPARG_fusion) <- names(Oncoplot_Data_Fusions_Included)
Oncoplot_Data_Fusions_Included <- rbind(Oncoplot_Data_Fusions_Included, PPARG_fusion)

# RET fusions
RET_fusion <- c(Clinical_Data_Oncoplot_Restricted$RET_fusion)
names(RET_fusion) <- names(Oncoplot_Data_Fusions_Included)
Oncoplot_Data_Fusions_Included <- rbind(Oncoplot_Data_Fusions_Included, RET_fusion)

# NTRK3 fusions
NTRK3_fusion <- c(Clinical_Data_Oncoplot_Restricted$NTRK3_fusion)
names(NTRK3_fusion) <- names(Oncoplot_Data_Fusions_Included)
Oncoplot_Data_Fusions_Included <- rbind(Oncoplot_Data_Fusions_Included, NTRK3_fusion)

# TERT fusions
TERT_fusion <- c(Clinical_Data_Oncoplot_Restricted$TERT_fusion)
names(TERT_fusion) <- names(Oncoplot_Data_Fusions_Included)
Oncoplot_Data_Fusions_Included <- rbind(Oncoplot_Data_Fusions_Included, TERT_fusion)

# NRG1 fusions
NRG1_fusion <- c(Clinical_Data_Oncoplot_Restricted$NRG1_fusion)
names(NRG1_fusion) <- names(Oncoplot_Data_Fusions_Included)
Oncoplot_Data_Fusions_Included <- rbind(Oncoplot_Data_Fusions_Included, NRG1_fusion)

# Change "Gene_fusion" to "Gene fusion" for better oncoplot labels, be sure to update if different fusions shown
Fusion_rownames = as_tibble(rownames(Oncoplot_Data_Fusions_Included))
for(i in 1:nrow(Fusion_rownames)){
  if (Fusion_rownames$value[i] == "MET_fusion"){
    Fusion_rownames$value[i] <- "MET fusion" 
  }
  else if (Fusion_rownames$value[i] == "BRAF_fusion"){
    Fusion_rownames$value[i] <- "BRAF fusion" 
  }
  else if (Fusion_rownames$value[i] == "PPARG_fusion"){
    Fusion_rownames$value[i] <- "PPARG fusion" 
  }
  else if (Fusion_rownames$value[i] == "RET_fusion"){
    Fusion_rownames$value[i] <- "RET fusion" 
  }
  else if (Fusion_rownames$value[i] == "NTRK3_fusion"){
    Fusion_rownames$value[i] <- "NTRK3 fusion" 
  }
  else if (Fusion_rownames$value[i] == "TERT_fusion"){
    Fusion_rownames$value[i] <- "TERT fusion" 
  }
  else if (Fusion_rownames$value[i] == "NRG1_fusion"){
    Fusion_rownames$value[i] <- "NRG1 fusion" 
  }
}

rownames(Oncoplot_Data_Fusions_Included) <- Fusion_rownames$value

# Replace "No" with "" and "Yes" with "Fusion"
for(y in 1:nrow(Oncoplot_Data_Fusions_Included)){
  for(x in 1:ncol(Oncoplot_Data_Fusions_Included)){
    if(Oncoplot_Data_Fusions_Included[y,x] == "No"){
      Oncoplot_Data_Fusions_Included[y,x] <- ""
    }
    else if(Oncoplot_Data_Fusions_Included[y,x] == "Yes"){
      Oncoplot_Data_Fusions_Included[y,x] <- "Fusion"
    }
  }
}

#### 9. Make Diagnosis into a numerical variable for sorting ####
# For sorting purposes
# The diagnosis that you want on the left should be 1
Clinical_Data_Oncoplot_Restricted$Diagnosis_Sorting <- 0
for(i in 1:nrow(Clinical_Data_Oncoplot_Restricted)){
  # Diagnosis sorting
  if(Clinical_Data_Oncoplot_Restricted$Diagnosis[i] == "ATC"){
    Clinical_Data_Oncoplot_Restricted$Diagnosis_Sorting[i] <- 1
  }
  else if(Clinical_Data_Oncoplot_Restricted$Diagnosis[i] == "PDTC"){
    Clinical_Data_Oncoplot_Restricted$Diagnosis_Sorting[i] <- 2
  }
  else if(Clinical_Data_Oncoplot_Restricted$Diagnosis[i] == "PTC and IFVPTC"){
    Clinical_Data_Oncoplot_Restricted$Diagnosis_Sorting[i] <- 3
  }
  else if(Clinical_Data_Oncoplot_Restricted$Diagnosis[i] == "NIFTP and EFVPTC"){
    Clinical_Data_Oncoplot_Restricted$Diagnosis_Sorting[i] <- 4
  }
  else if(Clinical_Data_Oncoplot_Restricted$Diagnosis[i] == "FTC"){
    Clinical_Data_Oncoplot_Restricted$Diagnosis_Sorting[i] <- 5
  }
  else if(Clinical_Data_Oncoplot_Restricted$Diagnosis[i] == "HC"){
    Clinical_Data_Oncoplot_Restricted$Diagnosis_Sorting[i] <- 6
  }
}

# Testing rowname congruence
test <- rownames(Clinical_Data_Oncoplot_Restricted) == colnames(Oncoplot_Data) # must be same order for annotations to work
test

#### 10. Formatting mutation grid of ComplexHeatmap oncoPrint plot ####

### Choosing colors for mutations ###
col = c("In_Frame_Del" = "gold1",             # Change IN_Frame_Del color
        "Nonsense_Mutation" = "red",           # Change Nonsense_Mutation color
        "Missense_Mutation" = "springgreen4",    # Change Missense_Mutation color
        "Frame_Shift_Ins" = "yellow",            # Change Frame_Shift_Del color
        "Frame_Shift_Del" = "blue",            # Change Frame_Shift_Del color
        "Splice_Site" = "darkorange",              # Change Splice_Site color
        "Multi_Hit" = "black",                 # Change Multi_Hit color
        "Fusion" = "purple4")                  # Change Fusion color

### Formatting mutation box sizes and background with alter_fun
# Can modify these to change the box styles as desired
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "pt"), h-unit(.5, "pt"), 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  # big blue
  In_Frame_Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "pt"), h-unit(0.5, "pt"), 
              gp = gpar(fill = col["In_Frame_Del"], col = NA))
  },
  # big red
  Nonsense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "pt"), h-unit(0.5, "pt"), 
              gp = gpar(fill = col["Nonsense_Mutation"], col = NA))
  },
  # big green
  Missense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "pt"), h-unit(0.5, "pt"),
              gp = gpar(fill = col["Missense_Mutation"], col = NA))
  },
  # big yellow
  Frame_Shift_Ins = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "pt"), h-unit(0.5, "pt"), 
              gp = gpar(fill = col["Frame_Shift_Ins"], col = NA))
  },
  # big red
  Frame_Shift_Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "pt"), h-unit(0.5, "pt"), 
              gp = gpar(fill = col["Frame_Shift_Del"], col = NA))
  },
  # big red
  Splice_Site = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "pt"), h-unit(0.5, "pt"), 
              gp = gpar(fill = col["Splice_Site"], col = NA))
  },
  Multi_Hit = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "pt"), h-unit(0.5, "pt"), 
              gp = gpar(fill = col["Multi_Hit"], col = NA))
  },
  # Fusions as a triangle
  Fusion = function(x, y, w, h){
    grid.polygon(
      unit.c(x + 0.25*w, x - 0.25*w, x - 0.25*w),
      unit.c(y         , y - 0.40*h, y + 0.40*h), 
      gp = gpar(fill = col["Fusion"], col = NA)
    )
  }
)



# Test alter_fun -> prints out the style of mutation boxes that have been generated
test_alter_fun(alter_fun)

column_title = "Oncoplot" # Change to include desired title name 
heatmap_legend_param = list(title = "Alterations", at = c("In_Frame_Del", "Nonsense_Mutation", "Missense_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del", "Splice_Site", "Multi_Hit", "Fusion"),
                            labels = c("In Frame Deletion", "Nonsense Mutation", "Missense Mutation", "Frame Shift Insertion", "Frame Shift Deletion", "Splice Site or Promoter", "Multi Hit", "Fusions"))

#### 11. Formatting ComplexHeatmap Annotations ####

### Annotation color gradients
BRS_Color <- colorRamp2(c(-0.8, 0, .8), c("red", "mediumorchid3", "blue")) # old color combo that I liked
TDS_Color <- colorRamp2(c(-1, 0, 1), c("yellow", "orange", "brown"))
ERK_Color <- colorRamp2(c(-50, 0, 50), c("aquamarine", "lavender", "#c51b8a"))
Age_Color <- colorRamp2(c(18, 54, 90), c("cadetblue1", "cornflowerblue", "blue4"))
PI3K_AKT_MTOR_Color <- colorRamp2(c(-100, 0, 100), c("green", "white", "red"))

### Make the annotations

# Top annotations
HA_Top <- HeatmapAnnotation(column_barplot = anno_oncoprint_barplot(border = TRUE,
                                                                    show_fraction = FALSE),
                            Diagnosis = Clinical_Data_Oncoplot_Restricted$Diagnosis, # Note: Can add Diagnosis as a label instead of an annotation...gives the ability to write in the box
                            'Location' = Clinical_Data_Oncoplot_Restricted$Location.type,
                            Sex = Clinical_Data_Oncoplot_Restricted$Sex,
                            Age = Clinical_Data_Oncoplot_Restricted$Age.resection,
                            'Aggressive disease' = Clinical_Data_Oncoplot_Restricted$Aggressive.disease,
                            TDS = Clinical_Data_Oncoplot_Restricted$TDS,
                            ERK = Clinical_Data_Oncoplot_Restricted$ERK,
                            BRS = Clinical_Data_Oncoplot_Restricted$BRS,
                            'PI3K-AKT-MTOR'  = Clinical_Data_Oncoplot_Restricted$PI3K_AKT_MTOR,
                            col = list(Diagnosis = c("ATC" = "magenta",  # These are the colors for adding Diagnosis as an annotation. Remove if adding as a label.
                                                     "PDTC" = "purple", 
                                                     "PTC and IFVPTC" = "red",
                                                     "NIFTP and EFVPTC" = "lightskyblue",
                                                     "FTC" = "blue",
                                                     "HC" = "dodgerblue2"),
                                       'Location' = c("Non-metastatic" = "Black",
                                                      "Local metastasis" = "Green",
                                                      "Distant metastasis" = "Red"),
                                       Sex = c("M" = "mediumpurple1", 
                                               "F" = "palevioletred1"),
                                       Age = Age_Color,
                                       'Aggressive disease' = c("Aggressive" = "Red",
                                                      "Indolent" = "Blue",
                                                      "NA" = "Grey"),
                                       BRS = BRS_Color,
                                       TDS = TDS_Color,
                                       ERK = ERK_Color,
                                       'PI3K-AKT-MTOR' = PI3K_AKT_MTOR_Color),
                            gp = gpar(col = "black"),
                            annotation_name_gp = gpar(col = "black", fontsize = 16, fontface = "bold"),
                            annotation_name_side = "left",
                            annotation_legend_param = list(Diagnosis=list(at=c("ATC","PDTC","PTC and IFVPTC","NIFTP and EFVPTC","FTC","HC")),
                                                          'Location'=list(at=c("Non-metastatic","Local metastasis", "Distant metastasis")),
                                                          'Aggressive disease'=list(at=c("Aggressive","Indolent")))) #recommend by Matt 4-5-22, see https://jokergoo.github.io/ComplexHeatmap-reference/book/

# Left annotations
HA_Left = rowAnnotation(foo = anno_block(gp = gpar(fill = c("mediumspringgreen", "mediumslateblue"), alpha = 0.5),
                                                   labels = c("Mutations", "Fusions"),
                                                   labels_gp = gpar(col = "black", fontsize = 20, fontface = "bold")))

#### 12. Make the ComplexHeatmap plot with oncoPrint ####

png(filename="ThyroidOrigin_Spaced_by_Diagnosis.png",
    res = 300,
    units = "in",
    pointsize = 20,
    width = 18,
    height = 12
)
ht <- oncoPrint(Oncoplot_Data_Fusions_Included,
                alter_fun = alter_fun, col = col,
                column_split = factor(Clinical_Data_Oncoplot_Restricted$Diagnosis_Sorting), # Creates column splits -> remove to remove column splits
                top_annotation = HA_Top, # Set the top annotation equal to HA_Top
                pct_side = "right", # Set percents to appear on the right of the plot
                pct_gp = gpar(col = "black", fontsize = 16, fontface = "bold"), # change 
                row_names_side = "left", # set row names (gene names) to appear on the left side of the plot
                row_names_gp = gpar(col = "black", fontsize = 16, fontface = "bold"), # change gene name font sizes/color/style (bold)
                column_title = column_title,
                column_title_gp = gpar(col = "black", fontsize = 16, fontface = "bold"),
                heatmap_legend_param = heatmap_legend_param, 
                row_split = factor(rep(c("Mutation", "Fusion"), times = c(20,7)), levels = c("Mutation", "Fusion")),
                column_order = 1:239, #holds column order in the order fed in by the oncoplot matrix, will need adjusting depending on cohort subsetting used (see environment, Clinical_Data_Oncoplot_Restricted, number of obs.)
                left_annotation = HA_Left,
                row_title = NULL) # Prevents double print of "mutation, fusion"
draw(ht, merge_legend = TRUE) # Merges legends from mutations and annotations into one column
dev.off()
