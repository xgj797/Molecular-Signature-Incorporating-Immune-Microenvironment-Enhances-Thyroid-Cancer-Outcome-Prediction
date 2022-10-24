### Figure 4C Violin Plots: Information
# This script reads in normalized expression data for our cohort of thyroid lesions, filters to non-metastatic lesions, and prints out violin plots for the general CAF Markers FAP and ACTA2

### Chapters: 
# 1. Load necessary packages
# 2. Read in data file
# 3. Prepare file for making violin plots of general CAF markers (FAP, ACTA2)
# 4. Restrict to non-metastatic lesions of interest
# 5. Make general CAF plots (FAP, ACTA2)

### Chapter 1. Load necessary packages
library(tidyverse)
library(RColorBrewer)

### Chapter 2. Read in data file
CleanedMergedData <- readRDS(file = "data_in_use/CleanedMergedData_DESeq2NormalizedReads.rds")

### Chapter 3. Prepare file for making violin plots of general CAF markers (FAP, ACTA2)
# Establish a column for labeling as Follicular, Papillary, or Transformed (will use to color plots)
CleanedMergedData$Category <- "Benign"

# Follicular lesions (these are grouped due to their RAS-like status by BRAF-RAS score)
for(i in 1:nrow(CleanedMergedData)){
  if(CleanedMergedData$Diagnosis.with.iFVPTC.and.eFVPTC[i] == "FA" |      # FA = Follicular Adenoma
     CleanedMergedData$Diagnosis.with.iFVPTC.and.eFVPTC[i] == "HA" |      # HA = Hurthle Cell Adenoma
     CleanedMergedData$Diagnosis.with.iFVPTC.and.eFVPTC[i] == "NIFTP" |   # NIFTP = Noninvasive follicular thyroid neoplasm with papillary-like nuclear features (NIFTP)
     CleanedMergedData$Diagnosis.with.iFVPTC.and.eFVPTC[i] == "FTC" |     # FTC = Follicular Thyroid Carcinoma
     CleanedMergedData$Diagnosis.with.iFVPTC.and.eFVPTC[i] == "HC" |      # HC = Hurthle Cell Carcinoma 
     CleanedMergedData$Diagnosis.with.iFVPTC.and.eFVPTC[i] == "EFVPTC"){  # EFVPTC = Encapsulated Follicular Variant of Papillary Thyroid Carcinoma
        CleanedMergedData$Category[i] <- "Follicular" # Set all of the above variants as "Follicular" for the purpose of coloring plots
  }
}

# Papillary lesions (These are groupeddue to their BRAF-like status by BRAF-RAS score)
for(i in 1:nrow(CleanedMergedData)){
  if(CleanedMergedData$Diagnosis.with.iFVPTC.and.eFVPTC[i] == "IFVPTC" | # IFVPTC = Infiltrative Follicular Variant of Papillary Thyroid Carcinoma
     CleanedMergedData$Diagnosis.with.iFVPTC.and.eFVPTC[i] == "PTC"){    # PTC = Papillary Thyroid Carcinoma
        CleanedMergedData$Category[i] <- "Papillary" # Set all of the above variants as "Follicular" for the purpose of coloring plots
  }
}

# Transformed lesion
for(i in 1:nrow(CleanedMergedData)){
  if(CleanedMergedData$Diagnosis.with.iFVPTC.and.eFVPTC[i] == "PDTC" | # PDTC = Poorly Differentiated Thyroid Carcinoma
     CleanedMergedData$Diagnosis.with.iFVPTC.and.eFVPTC[i] == "ATC"){  # ATC = Anaplastic Thyroid Carcinoma
        CleanedMergedData$Category[i] <- "Transformed" # Set all of the above as "Transformed" for the purpose of coloring plots
  }
}

# Create a new variable to group diagnoses for plotting (Diagnosis_Simplified)
# Grouping FA, HA, FTC, HC, NIFTP, and EFVPTC under "Follicular"
CleanedMergedData$Diagnosis_Simplified <- CleanedMergedData$Diagnosis.with.iFVPTC.and.eFVPTC
for(i in 1:nrow(CleanedMergedData)){
  if(CleanedMergedData$Diagnosis_Simplified[i] == "FA" | 
     CleanedMergedData$Diagnosis_Simplified[i] == "HA" |
     CleanedMergedData$Diagnosis_Simplified[i] == "FTC" |
     CleanedMergedData$Diagnosis_Simplified[i] == "HC" |
     CleanedMergedData$Diagnosis_Simplified[i] == "NIFTP" |
     CleanedMergedData$Diagnosis_Simplified[i] == "EFVPTC"){
        CleanedMergedData$Diagnosis_Simplified[i] <- "Follicular"
  }
}

# Additional grouping of diagnoses for plotting useing the "Diagnosis_Simplified" variable created above
# Grouping PTC with IFVPTC (PTC/IFVPTC)
for(i in 1:nrow(CleanedMergedData)){
  if(CleanedMergedData$Diagnosis_Simplified[i] == "PTC" | CleanedMergedData$Diagnosis_Simplified[i] == "IFVPTC"){
    CleanedMergedData$Diagnosis_Simplified[i] <- "PTC/\nIFVPTC" # \n creates a space between PTC and IFVPTC on the plot
  }
}

# Create a new variable to add a BRAF and RAS component to the simplified diagnosis
CleanedMergedData$Diagnosis_Simplified_BRS <- CleanedMergedData$Diagnosis_Simplified
CleanedMergedData$Category_BRS <- CleanedMergedData$Category
for(i in 1:nrow(CleanedMergedData)){
  if(CleanedMergedData$Diagnosis_Simplified_BRS[i] == "ATC"){
    if(CleanedMergedData$BRS[i] > 0){ # ATC samples with BRAF-RAS Score (BRS) > 0 are RAS-like
      CleanedMergedData$Diagnosis_Simplified_BRS[i] <- "  ATC \n RAS" # \n creates a space between ATC and RAS for the purpose of plotting
    }
    else if(CleanedMergedData$BRS[i] < 0){ # ATC samples with BRAF-RAS Score (BRS) < 0 are BRAF-like
      CleanedMergedData$Diagnosis_Simplified_BRS[i] <- " ATC \n BRAF" # \n creates a space between ATC and BRAF for the purpose of plotting
    }
  }
}

### Chapter 4. Restrict to non-metastatic lesions of interest
# Restrict to non-metastatic lesions
CleanedMergedData_NonMet <- CleanedMergedData %>% subset(Location.type == "Non-Metastic")

# Restrict to histotypes of interest
CleanedMergedData_NonMet_Restricted <- CleanedMergedData_NonMet %>% subset(Diagnosis_Simplified != "HT" & # Remove # Not plotting Hashimoto's thyroiditis (HT)
                                                                           Diagnosis_Simplified != "MNG" &  # Not plotting multinodular goiter (MNG)
                                                                           Diagnosis.with.HTPTC != "HTPTC") # Excluding PTC samples with background Hashimoto's thyroiditis due to confounding expression of iCAF markers

### Chapter 5. Make general CAF plots (FAP, ACTA2)
## FAP
# Log2 transform the FAP count data for plotting
CleanedMergedData_NonMet_Restricted$FAP <- log(CleanedMergedData_NonMet_Restricted$FAP+1, 2)

# Create a variable Plot_FAP that is a ggplot of the FAP data by diagnosis
Plot_FAP <- ggplot(CleanedMergedData_NonMet_Restricted, aes(Diagnosis_Simplified_BRS, FAP)) +
  geom_violin(outlier.size = -1, 
              aes(fill = Category_BRS),
              alpha = 0.4, 
              show.legend = FALSE, 
              scale = "width") + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5, 
              alpha = 0.7,
              show.legend = FALSE) +
  scale_fill_manual(values = c("blue", "red", "purple")) + # Blue = Follicular; Red = Papillary; Puruple = Transformed
  labs (x = "Diagnosis", y = "Log2\n(Normalized FAP Counts)") + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 13),
    axis.text.x = element_text(face = "bold", size = 10.5),
    axis.text.y = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="Diagnosis", limits = c("Follicular", "PTC/\nIFVPTC", "PDTC", "  ATC \n RAS", " ATC \n BRAF")) +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10, 12),
                     limits = c(0, 12))
                     
# Save the Plot_FAP as a png                       
ggsave("outputs/Transformation_Series_FAP_BRS_Violin.png", 
       width = 4, height = 3, 
       Plot_FAP, dpi = 600)

# Non-parametric statistics for FAP by diagnosis w/ bonferroni correction on the pairwise.wilcox.test
kruskal.test(FAP ~ Diagnosis_Simplified_BRS, data = CleanedMergedData_NonMet_Restricted) # 
pairwise.wilcox.test(CleanedMergedData_NonMet_Restricted$FAP, CleanedMergedData_NonMet_Restricted$Diagnosis_Simplified_BRS, 
                     p.adjust.method = "bonferroni")
                     

## ACTA2
# Log2 transform the ACTA2 count data for plotting
CleanedMergedData_NonMet_Restricted$ACTA2 <- log(CleanedMergedData_NonMet_Restricted$ACTA2+1, 2)

# Create a variable Plot_ACTA2 that is a ggplot of the ACTA2 data by diagnosis
Plot_ACTA2 <- ggplot(CleanedMergedData_NonMet_Restricted, aes(Diagnosis_Simplified_BRS, ACTA2)) +
  geom_violin(outlier.size = -1, 
              aes(fill = Category_BRS),
              alpha = 0.4, 
              show.legend = FALSE, 
              scale = "width") + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5, 
              alpha = 0.7,
              show.legend = FALSE) +
  scale_fill_manual(values = c("blue", "red", "purple")) + # Blue = Follicular; Red = Papillary; Puruple = Transformed
  labs (x = "Diagnosis", y = "Log2\n(Normalized ACTA2 Counts)") + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 13),
    axis.text.x = element_text(face = "bold", size = 10.5),
    axis.text.y = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="Diagnosis", limits = c("Follicular", "PTC/\nIFVPTC", "PDTC", "  ATC \n RAS", " ATC \n BRAF")) +
  scale_y_continuous(breaks = c(7, 8, 9, 10, 11, 12),
                     limits = c(6.9, 12)) 
                     
# Save the Plot_ACTA2 as a png                       
ggsave("outputs/Transformation_Series_ACTA2_BRS_Violin.png", 
       width = 4, height = 3, 
       Plot_ACTA2, dpi = 600)

# Non-parametric statistics for ACTA2 by diagnosis w/ bonferroni correction on the pairwise.wilcox.test
kruskal.test(ACTA2 ~ Diagnosis_Simplified_BRS, data = CleanedMergedData_NonMet_Restricted) # 
pairwise.wilcox.test(CleanedMergedData_NonMet_Restricted$ACTA2, CleanedMergedData_NonMet_Restricted$Diagnosis_Simplified_BRS, 
                     p.adjust.method = "bonferroni")

