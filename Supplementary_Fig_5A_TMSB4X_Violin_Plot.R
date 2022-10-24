### Supplementary Fig 5A Violin Plot for TMSB4X: Information
# This script reads in normalized expression data for our cohort of thyroid lesions, filters to non-metastatic lesions, and prints out a violin plot for TMSB4X
# TMSB4X has previously been shown as a predictor of aggressive disease in PTCs (Nature Communications PTC scRNA-Seq paper)

### Chapters: 
# 1. Load necessary packages
# 2. Read in data file
# 3. Prepare file for making violin plot of TMSB4X
# 4. Restrict to non-metastatic lesions of interest
# 5. Make TMSB4X Violin Plot

### Chapter 1. Load necessary packages
library(tidyverse)
library(RColorBrewer)

### Chapter 2. Read in data file
CleanedMergedData <- readRDS(file = "data_in_use/CleanedMergedData_DESeq2NormalizedReads.rds")

### Chapter 3. Prepare file for making violin plot of TMSB4X
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
# Grouping FA with HA (FA/HA) and FTC with HC (FTC/HC)
CleanedMergedData$Diagnosis_Simplified <- CleanedMergedData$Diagnosis.with.iFVPTC.and.eFVPTC
for(i in 1:nrow(CleanedMergedData)){
  if(CleanedMergedData$Diagnosis_Simplified[i] == "FA" | CleanedMergedData$Diagnosis_Simplified[i] == "HA"){
    CleanedMergedData$Diagnosis_Simplified[i] <- "FA/HA"
  }
  else if(CleanedMergedData$Diagnosis_Simplified[i] == "FTC" | CleanedMergedData$Diagnosis_Simplified[i] == "HC"){
    CleanedMergedData$Diagnosis_Simplified[i] <- "FTC/\nHC" # \n creates a space between FTC and HC on the plot
  }
}

# Additional grouping of diagnoses for plotting useing the "Diagnosis_Simplified" variable created above
# Grouping EFVPTC with NIFTP (EFVPTC/NIFTP) and PTC with IFVPTC (PTC/IFVPTC)
for(i in 1:nrow(CleanedMergedData)){
  if(CleanedMergedData$Diagnosis_Simplified[i] == "NIFTP" | CleanedMergedData$Diagnosis_Simplified[i] == "EFVPTC"){
    CleanedMergedData$Diagnosis_Simplified[i] <- "EFVPTC/\nNIFTP" # \n creates a space between EFVPTC and NIFTP on the plot
  }
  else if(CleanedMergedData$Diagnosis_Simplified[i] == "PTC" | CleanedMergedData$Diagnosis_Simplified[i] == "IFVPTC"){
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

### Chapter 5. Make TMSB4X violin plot
## TMSB4X
# Log2 transform the TMSB4X count data for plotting
CleanedMergedData_NonMet_Restricted$TMSB4X <- log(CleanedMergedData_NonMet_Restricted$TMSB4X+1, 2)

# Create a variable Plot_TMSB4X that is a ggplot of the TMSB4X data 
Plot_TMSB4X <- ggplot(CleanedMergedData_NonMet_Restricted, aes(Diagnosis_Simplified_BRS, TMSB4X)) +
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
  labs (x = "Diagnosis", y = "Log2\n(Normalized TMSB4X Counts") + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 13),
    axis.text.x = element_text(face = "bold", size = 10.5),
    axis.text.y = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="Diagnosis", limits = c("FA/HA", "FTC/\nHC", "EFVPTC/\nNIFTP", "PTC/\nIFVPTC", "PDTC", "  ATC \n RAS", " ATC \n BRAF")) +
  scale_y_continuous(breaks = c(9, 11, 13, 15),
                     limits = c(9, 16)) 
                     
# Save the Plot_TMSB4X as a png                    
ggsave("outputs/Transformation_Series_TMSB4X_BRS_Violin.png", 
       width = 5, height = 3, 
       Plot_TMSB4X, dpi = 600)

# Non-parametric statistics for TMSB4X w/ bonferroni correction on the pairwise.wilcox.test
kruskal.test(TMSB4X ~ Diagnosis_Simplified_BRS, data = CleanedMergedData_NonMet_Restricted) # shouldn't be necessary for only two groups
pairwise.wilcox.test(CleanedMergedData_NonMet_Restricted$TMSB4X, CleanedMergedData_NonMet_Restricted$Diagnosis_Simplified_BRS, 
                     p.adjust.method = "bonferroni")
