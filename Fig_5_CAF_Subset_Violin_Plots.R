### Figure 5 Violin Plots: Information
# This script reads in normalized expression data for our cohort of thyroid lesions, filters to non-metastatic lesions, and prints out violin plots for iCAF and myCAF markers
# myCAF Markers: POSTN, COL1A1, WNT2, LRRC15
# iCAF Markers: PDGFRA, CXCL12

### Chapters: 
# 1. Load necessary packages
# 2. Read in data file
# 3. Prepare file for making violin plots of myCAF and iCAF markers
# 4. Restrict to non-metastatic lesions of interest
# 5. Make myCAF Plots
# 6. Make iCAF Plots

### Chapter 1. Load necessary packages
library(tidyverse)
library(RColorBrewer)

### Chapter 2. Read in data file
CleanedMergedData <- readRDS(file = "data_in_use/CleanedMergedData_DESeq2NormalizedReads.rds")

### Chapter 3. Prepare file for making violin plots of myCAF and iCAF Markers
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

### Chapter 5. Make myCAF Plots
## POSTN
# Log2 transform the POSTN count data for plotting
CleanedMergedData_NonMet_Restricted$POSTN <- log(CleanedMergedData_NonMet_Restricted$POSTN+1, 2)

# Create a variable Plot_POSTN that is a ggplot of the POSTN data 
Plot_POSTN <- ggplot(CleanedMergedData_NonMet_Restricted, aes(Diagnosis_Simplified_BRS, POSTN)) +
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
  labs (x = "Diagnosis", y = "Log2\n(Normalized POSTN Counts") + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 13),
    axis.text.x = element_text(face = "bold", size = 10.5),
    axis.text.y = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="Diagnosis", limits = c("FA/HA", "FTC/\nHC", "EFVPTC/\nNIFTP", "PTC/\nIFVPTC", "PDTC", "  ATC \n RAS", " ATC \n BRAF")) +
  scale_y_continuous(breaks = c(5, 7, 9, 11, 13, 15),
                     limits = c(4.5,15.5)) 
                     
# Save the Plot_POSTN as a png                    
ggsave("outputs/Transformation_Series_POSTN_BRS_Violin.png", 
       width = 5, height = 3, 
       Plot_POSTN, dpi = 600)

# Non-parametric statistics for POSTN w/ bonferroni correction on the pairwise.wilcox.test
kruskal.test(POSTN ~ Diagnosis_Simplified_BRS, data = CleanedMergedData_NonMet_Restricted) # shouldn't be necessary for only two groups
pairwise.wilcox.test(CleanedMergedData_NonMet_Restricted$POSTN, CleanedMergedData_NonMet_Restricted$Diagnosis_Simplified_BRS, 
                     p.adjust.method = "bonferroni")
                     

## WNT2
# Log2 transform the WNT2 count data for plotting
CleanedMergedData_NonMet_Restricted$WNT2 <- log(CleanedMergedData_NonMet_Restricted$WNT2+1, 2)

# Create a variable Plot_WNT2 that is a ggplot of the WNT2 data 
Plot_WNT2 <- ggplot(CleanedMergedData_NonMet_Restricted, aes(Diagnosis_Simplified_BRS, WNT2)) +
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
  labs (x = "Diagnosis", y = "Log2\n(Normalized WNT2 Counts") + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 13),
    axis.text.x = element_text(face = "bold", size = 10.5),
    axis.text.y = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="Diagnosis", limits = c("FA/HA", "FTC/\nHC", "EFVPTC/\nNIFTP", "PTC/\nIFVPTC", "PDTC", "  ATC \n RAS", " ATC \n BRAF")) +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8),
                     limits = c(0, 8.5)) 
                     
# Save the Plot_WNT2 as a png                    
ggsave("outputs/Transformation_Series_WNT2_BRS_Violin.png", 
       width = 5, height = 3, 
       Plot_WNT2, dpi = 600)

# Non-parametric statistics for WNT2 w/ bonferroni correction on the pairwise.wilcox.test
kruskal.test(WNT2 ~ Diagnosis_Simplified_BRS, data = CleanedMergedData_NonMet_Restricted) # shouldn't be necessary for only two groups
pairwise.wilcox.test(CleanedMergedData_NonMet_Restricted$WNT2, CleanedMergedData_NonMet_Restricted$Diagnosis_Simplified_BRS, 
                     p.adjust.method = "bonferroni")
                     

## COL1A1
# Log2 transform the COL1A1 count data for plotting
CleanedMergedData_NonMet_Restricted$COL1A1 <- log(CleanedMergedData_NonMet_Restricted$COL1A1+1, 2)

# Create a variable Plot_COL1A1 that is a ggplot of the COL1A1 data 
Plot_COL1A1 <- ggplot(CleanedMergedData_NonMet_Restricted, aes(Diagnosis_Simplified_BRS, COL1A1)) +
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
  labs (x = "Diagnosis", y = "Log2\n(Normalized COL1A1 Counts") + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 13),
    axis.text.x = element_text(face = "bold", size = 10.5),
    axis.text.y = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="Diagnosis", limits = c("FA/HA", "FTC/\nHC", "EFVPTC/\nNIFTP", "PTC/\nIFVPTC", "PDTC", "  ATC \n RAS", " ATC \n BRAF")) +
  scale_y_continuous(breaks = c(8, 10, 12, 14, 16, 18),
                     limits = c(8, 18)) 

# Save the Plot_COL1A1 as a png                    
ggsave("outputs/Transformation_Series_COL1A1_BRS_Violin.png", 
       width = 5, height = 3, 
       Plot_COL1A1, dpi = 600)

# Non-parametric statistics for COL1A1 w/ bonferroni correction on the pairwise.wilcox.test
kruskal.test(COL1A1 ~ Diagnosis_Simplified_BRS, data = CleanedMergedData_NonMet_Restricted) # shouldn't be necessary for only two groups
pairwise.wilcox.test(CleanedMergedData_NonMet_Restricted$COL1A1, CleanedMergedData_NonMet_Restricted$Diagnosis_Simplified_BRS, 
                     p.adjust.method = "bonferroni")
                     
                     
  ## LRRC15
# Log2 transform the LRRC15 count data for plotting
CleanedMergedData_NonMet_Restricted$LRRC15 <- log(CleanedMergedData_NonMet_Restricted$LRRC15+1, 2)

# Create a variable Plot_LRRC15 that is a ggplot of the LRRC15 data 
Plot_LRRC15 <- ggplot(CleanedMergedData_NonMet_Restricted, aes(Diagnosis_Simplified_BRS, LRRC15)) +
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
  labs (x = "Diagnosis", y = "Log2\n(Normalized LRRC15 Counts") + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 13),
    axis.text.x = element_text(face = "bold", size = 10.5),
    axis.text.y = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="Diagnosis", limits = c("FA/HA", "FTC/\nHC", "EFVPTC/\nNIFTP", "PTC/\nIFVPTC", "PDTC", "  ATC \n RAS", " ATC \n BRAF")) +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10),
                     limits = c(0, 11.2)) 

# Save the Plot_LRRC15 as a png                    
ggsave("outputs/Transformation_Series_LRRC15_BRS_Violin.png", 
       width = 5, height = 3, 
       Plot_LRRC15, dpi = 600)

# Non-parametric statistics for LRRC15 w/ bonferroni correction on the pairwise.wilcox.test
kruskal.test(LRRC15 ~ Diagnosis_Simplified_BRS, data = CleanedMergedData_NonMet_Restricted) # shouldn't be necessary for only two groups
pairwise.wilcox.test(CleanedMergedData_NonMet_Restricted$LRRC15, CleanedMergedData_NonMet_Restricted$Diagnosis_Simplified_BRS, 
                     p.adjust.method = "bonferroni")                   


### Chapter 6. Make iCAF Plots
## PDGFRA
# Log2 transform the PDGFRA count data for plotting
CleanedMergedData_NonMet_Restricted$PDGFRA <- log(CleanedMergedData_NonMet_Restricted$PDGFRA+1, 2)

# Create a variable Plot_PDGFRA that is a ggplot of the PDGFRA data 
Plot_PDGFRA <- ggplot(CleanedMergedData_NonMet_Restricted, aes(Diagnosis_Simplified_BRS, PDGFRA)) +
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
  labs (x = "Diagnosis", y = "Log2\n(Normalized PDGFRA Counts") + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 13),
    axis.text.x = element_text(face = "bold", size = 10.5),
    axis.text.y = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="Diagnosis", limits = c("FA/HA", "FTC/\nHC", "EFVPTC/\nNIFTP", "PTC/\nIFVPTC", "PDTC", "  ATC \n RAS", " ATC \n BRAF")) +
  scale_y_continuous(breaks = c(0, 3, 6, 9, 12),
                     limits = c(0, 12)) 

# Save the Plot_PDGFRA as a png                    
ggsave("outputs/Transformation_Series_PDGFRA_BRS_Violin.png", 
       width = 5, height = 3, 
       Plot_PDGFRA, dpi = 600)

# Non-parametric statistics for PDGFRA w/ bonferroni correction on the pairwise.wilcox.test
kruskal.test(PDGFRA ~ Diagnosis_Simplified_BRS, data = CleanedMergedData_NonMet_Restricted) # shouldn't be necessary for only two groups
pairwise.wilcox.test(CleanedMergedData_NonMet_Restricted$PDGFRA, CleanedMergedData_NonMet_Restricted$Diagnosis_Simplified_BRS, 
                     p.adjust.method = "bonferroni")


## CXCL12
# Log2 transform the CXCL12 count data for plotting
CleanedMergedData_NonMet_Restricted$CXCL12 <- log(CleanedMergedData_NonMet_Restricted$CXCL12+1, 2)

# Create a variable Plot_CXCL12 that is a ggplot of the CXCL12 data 
Plot_CXCL12 <- ggplot(CleanedMergedData_NonMet_Restricted, aes(Diagnosis_Simplified_BRS, CXCL12)) +
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
  labs (x = "Diagnosis", y = "Log2\n(Normalized CXCL12 Counts") + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 13),
    axis.text.x = element_text(face = "bold", size = 10.5),
    axis.text.y = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="Diagnosis", limits = c("FA/HA", "FTC/\nHC", "EFVPTC/\nNIFTP", "PTC/\nIFVPTC", "PDTC", "  ATC \n RAS", " ATC \n BRAF")) +
  scale_y_continuous(breaks = c(2, 5, 8, 11, 14),
                     limits = c(1.8, 14)) 

# Save the Plot_CXCL12 as a png                    
ggsave("outputs/Transformation_Series_CXCL12_BRS_Violin.png", 
       width = 5, height = 3, 
       Plot_CXCL12, dpi = 600)

# Non-parametric statistics for CXCL12 w/ bonferroni correction on the pairwise.wilcox.test
kruskal.test(CXCL12 ~ Diagnosis_Simplified_BRS, data = CleanedMergedData_NonMet_Restricted) # shouldn't be necessary for only two groups
pairwise.wilcox.test(CleanedMergedData_NonMet_Restricted$CXCL12, CleanedMergedData_NonMet_Restricted$Diagnosis_Simplified_BRS, 
                     p.adjust.method = "bonferroni")

                                                                           
                                                                           
