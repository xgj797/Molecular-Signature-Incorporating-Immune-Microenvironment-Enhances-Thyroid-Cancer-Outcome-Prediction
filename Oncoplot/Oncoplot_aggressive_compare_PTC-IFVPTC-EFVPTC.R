# Script: Oncoplot_aggressive_compare_PTC-IFVPTC-EFVPTC.R

#### Guide for ComplexHeatmap ####
# https://jokergoo.github.io/ComplexHeatmap-reference/book/
# See chapter 7 for oncoPrint instructions

#### Table of Contents ####
# 1. Load packages
# 2. Read maf file
# 3. Subset maf to samples of interest
# 4. Use maf tools to print a matrix to feed into ComplexHeatmap's oncoprint function
# 5. Read in Oncoplot Matrix generated with maftools
# 6. Format and subset Oncoplot Matrix and Clinical Data
# 7. Add in fusion data to Oncoplot Data
# 8. Make Diagnosis into a numerical variable for sorting
# 9. Formatting mutation grid of ComplexHeatmap oncoPrint plot
# 10. Formatting ComplexHeatmap annotations
# 11. Make the ComplexHeatmap with oncoPrint

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
Clinical_Data <- read.csv("Diagnosis_guide_9-28-22_for_supfig.csv", header = TRUE)
Clinical_Data <- as_tibble(Clinical_Data)

# remove normals -> no resection data (e.g., sequencing, etc.)
Clinical_Data <- Clinical_Data %>% subset(Diagnosis != "normal")

# subset to either aggressive or indolent
#Clinical_Data <- Clinical_Data %>% subset(Aggressive.disease == "Indolent")
Clinical_Data <- Clinical_Data %>% subset(Aggressive.disease == "Aggressive")

# Needs to be in vector format for oncoplot
Clinical_Data <- as.data.frame(Clinical_Data)

#### 3. Subset maf to samples of interest ####
#Non-metastatic only, comment out if want to include all samples
Clinical_Data <- Clinical_Data %>% subset(Location.type == "Primary" | Location.type == "Localdisease")

#Diagnoses of interest
Clinical_Data <- Clinical_Data %>% subset(Diagnosis == "PTC" | Diagnosis == "IFVPTC"| Diagnosis == "EFVPTC")

maf_thyroidorigin = subsetMaf(maf = maf, tsb = c(Clinical_Data$Tumor_Sample_Barcode), mafObj = TRUE)

#### 4. Use MAF tools to print a Matrix to feed into ComplexHeatmap's oncoPrint function ####
oncoplot(maf = maf_thyroidorigin, top = 20, annotationDat = Clinical_Data, # top = set the number of genes to look at, annotationDat = data for annotations
         clinicalFeatures = c("Diagnosis", # Annotation features to show at bottom
                              "Sex",
                              "Age.resection"),
         anno_height = 2,
         removeNonMutated = FALSE, # Shows "Dark Matter"
         sortByAnnotation = TRUE, # Must sort by annotation of Diagnosis ... because Diagnosis is listed first it is automatically used
         annotationOrder = c("PTC","IFVPTC","EFVPTC"), # Here is the order that I want the diagnoses in (save for ref: "FC","FCOT","NIFTP")
         numericAnnoCol = NULL, # This is for changing the color of annotation gradients...I couldn't get it to work in maf tools though
         writeMatrix = TRUE) # This is the most important step -> I will use writeMatrix as input for the complex heatmap package

#### 5. Read in Oncoplot Matrix generated with maftools ####
Oncoplot_Data <- as.matrix(read.table(file = "onco_matrix.txt", sep = '\t', header = TRUE))

#### 6. Format and subset Oncoplot Matrix and Clinical Data ####
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

#### 7. Add in fusion data to Oncoplot_Data ####
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

# Testing rowname congruence
test <- rownames(Clinical_Data_Oncoplot_Restricted) == colnames(Oncoplot_Data) # must be same order for annotations to work
test

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

#### 8. Make Diagnosis into a numerical variable for sorting ####
# For sorting purposes
# The Diagnosis that you want on the left should be 1
Clinical_Data_Oncoplot_Restricted$Diagnosis_Sorting <- 0
for(i in 1:nrow(Clinical_Data_Oncoplot_Restricted)){
   if(Clinical_Data_Oncoplot_Restricted$Diagnosis[i] == "PTC"){
    Clinical_Data_Oncoplot_Restricted$Diagnosis_Sorting[i] <- 1
   }
  if(Clinical_Data_Oncoplot_Restricted$Diagnosis[i] == "IFVPTC"){
    Clinical_Data_Oncoplot_Restricted$Diagnosis_Sorting[i] <- 2
  }
  else if(Clinical_Data_Oncoplot_Restricted$Diagnosis[i] == "EFVPTC"){
    Clinical_Data_Oncoplot_Restricted$Diagnosis_Sorting[i] <- 3
  }
}

# Testing rowname congruence
test <- rownames(Clinical_Data_Oncoplot_Restricted) == colnames(Oncoplot_Data) # must be same order for annotations to work
test

#### 9. Formatting mutation grid of ComplexHeatmap oncoPrint plot ####

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

#### 10. Formatting ComplexHeatmap Annotations ####

### Annotation color gradients
BRS_Color <- colorRamp2(c(-0.8, 0, .8), c("red", "mediumorchid3", "blue")) # old color combo that I liked
TDS_Color <- colorRamp2(c(-1, 0, 1), c("yellow", "orange", "brown"))
ERK_Color <- colorRamp2(c(-50, 0, 50), c("aquamarine", "lavender", "#c51b8a"))
Age_Color <- colorRamp2(c(18, 54, 90), c("cadetblue1", "cornflowerblue", "blue4"))
PI3K_AKT_MTOR_Color <- colorRamp2(c(-100, 0, 100), c("green", "white", "red"))
WNT_Canon_Color <- colorRamp2(c(-100, 0, 100), c("mediumspringgreen", "white", "mediumslateblue"))
WNT_NonCanon_Color <- colorRamp2(c(-50, 0, 50), c("mediumspringgreen", "white", "mediumslateblue"))
Cancer.associated.fibroblast_EPIC_Color <- colorRamp2(c(0, 0.5, 1), c("white", "red", "blue"))

### Make the annotations

# Top annotations
HA_Top <- HeatmapAnnotation(column_barplot = anno_oncoprint_barplot(border = TRUE,
                                                                    show_fraction = FALSE),
                            Diagnosis = Clinical_Data_Oncoplot_Restricted$Diagnosis, # Note: Can add Diagnosis as a label instead of an annotation...gives the ability to write in the box
                            #'Location type' = Clinical_Data_Oncoplot_Restricted$Location.type,
                            #Sex = Clinical_Data_Oncoplot_Restricted$Sex,
                            #Age = Clinical_Data_Oncoplot_Restricted$Age.resection,
                            #TDS = Clinical_Data_Oncoplot_Restricted$TDS,
                            #ERK = Clinical_Data_Oncoplot_Restricted$ERK,
                            #BRS = Clinical_Data_Oncoplot_Restricted$BRS,
                            #'PI3K-AKT-MTOR'  = Clinical_Data_Oncoplot_Restricted$PI3K_AKT_MTOR,
                            #'WNT Canonical' = Clinical_Data_Oncoplot_Restricted$WNT_Canon,
                            #'WNT NonCanonical' = Clinical_Data_Oncoplot_Restricted$WNT_NonCanon,
                            #'CAF EPIC' = Clinical_Data_Oncoplot_Restricted$Cancer.associated.fibroblast_EPIC,
                            col = list(Diagnosis = c("PTC" = "red",
                                                     "IFVPTC" = "red4",
                                                     "EFVPTC" = "lightskyblue")),
                                       #'Location type' = c("Primary" = "Black",
                                        #                 "Localdisease" = "Blue")),
                                       #Sex = c("M" = "mediumpurple1", 
                                        #       "F" = "palevioletred1"),
                                       #Age = Age_Color,
                                       #BRS = BRS_Color,
                                       #TDS = TDS_Color,
                                       #ERK = ERK_Color,
                                       #'PI3K-AKT-MTOR' = PI3K_AKT_MTOR_Color,
                                       #'WNT Canonical' = WNT_Canon_Color,
                                       #'WNT NonCanonical' = WNT_NonCanon_Color,
                                       #'CAF EPIC' = Cancer.associated.fibroblast_EPIC_Color
                            gp = gpar(col = "black"),
                            annotation_name_gp = gpar(col = "black", fontsize = 11, fontface = "bold"),
                            annotation_name_side = "left") #recommend by Matt 3-14-22

# Left annotations
HA_Left = rowAnnotation(foo = anno_block(gp = gpar(fill = c("mediumspringgreen", "mediumslateblue"), alpha = 0.5),
                                                   labels = c("Mutations", "Fusions"),
                                                   labels_gp = gpar(col = "black", fontsize = 12, fontface = "bold")))

#### 11. Make the ComplexHeatmap plot with oncoPrint ####

#png(filename="PTC-FVPTC_Indolent.png",
png(filename="PTC-FVPTC_Aggressive.png",
    res = 300,
    units = "in",
    pointsize = 18,
    width = 9,
    height = 6
)
ht <- oncoPrint(Oncoplot_Data_Fusions_Included,
                alter_fun = alter_fun, col = col,
                column_split = factor(Clinical_Data_Oncoplot_Restricted$Diagnosis_Sorting), # Creates column splits -> remove to remove column splits
                top_annotation = HA_Top, # Set the top annotation equal to HA_Top
                pct_side = "right", # Set percents to appear on the right of the plot
                pct_gp = gpar(col = "black", fontsize = 11, fontface = "bold"), # change 
                row_names_side = "left", # set row names (gene names) to appear on the left side of the plot
                row_names_gp = gpar(col = "black", fontsize = 11, fontface = "bold"), # change gene name font sizes/color/style (bold)
                column_title = column_title,
                column_title_gp = gpar(col = "black", fontsize = 12, fontface = "bold"),
                heatmap_legend_param = heatmap_legend_param, 
                row_split = factor(rep(c("Mutation", "Fusion"), times = c(20,7)), levels = c("Mutation", "Fusion")),
#                column_order = 1:58, # holds column order in the order fed in by the oncoplot matrix, will need adjusting depending on cohort subsetting used (see environment, Clinical_Data_Oncoplot_Restricted, number of obs.)
                column_order = 1:21, # holds column order in the order fed in by the oncoplot matrix, will need adjusting depending on cohort subsetting used (see environment, Clinical_Data_Oncoplot_Restricted, number of obs.)
                left_annotation = HA_Left,
                row_title = NULL) # Prevents double print of "mutation, fusion"
draw(ht, merge_legend = TRUE) # Merges legends from mutations and annotations into one column
dev.off()
