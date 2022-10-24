
### Figure 6 CAF High vs CAF Low Analysis: Information
# This script reads in normalized expression data for our cohort of thyroid lesions, filters to CAF high vs CAF low ATCs by CAF EPIC (TIMER 2.0) and plots the following markers:
# TIMER immune markers: Neutrophils, Macrophages
# CIBERSORT absolute value markers: M1 Macrophages, M2 Macrophages
# M2 Normalized Marker Expression: MRC1, CD163
# TIDE: TIDE Score, Exclusion Score, Dysfunction Score

### Chapters: 
# 1. Load necessary packages
# 2. Read in data file
# 3. Rename variables of interest and set to numerics
# 4. Log2 transform variables of interest
# 5. Establish ATC-High and ATC-Low CAF definitions (EPIC)
# 6. Make CAF-High vs CAF-Low plots of interest for Figure 6 (Neutrophil, Macrophage, M1, M2, MRC1, CD163, TIDE, Exclusion, Dysfunction)

### Chapter 1: Load necessary packages
library(tidyverse)

### Chapter 2: Read in data file
CleanedMergedData <- readRDS(file = "data_in_use/CleanedMergedData_DESeq2NormalizedReads.rds")

### Chapter 3: Rename variables of interest and set to numerics
# Rename CAF EPIC in data file for ease of use
CleanedMergedData <- CleanedMergedData %>% dplyr::rename("CAF_EPIC" = "Cancer associated fibroblast_EPIC")
# Set CAF EPIC as a numeric (data read in has it as a string)
CleanedMergedData$CAF_EPIC <- as.numeric(CleanedMergedData$CAF_EPIC)

# Rename other variables of interest (e.g., macrophage variables)
CleanedMergedData <- CleanedMergedData %>% dplyr::rename("M1" = "Macrophage M1_CIBERSORT-ABS")
CleanedMergedData <- CleanedMergedData %>% dplyr::rename("M2" = "Macrophage M2_CIBERSORT-ABS")

# Set variables of interest to numerics
CleanedMergedData$Macrophage_TIMER <- as.numeric(CleanedMergedData$Macrophage_TIMER)
CleanedMergedData$M1 <- as.numeric(CleanedMergedData$M1)
CleanedMergedData$M2 <- as.numeric(CleanedMergedData$M2)
CleanedMergedData$Neutrophil_TIMER <- as.numeric(CleanedMergedData$Neutrophil_TIMER)
CleanedMergedData$TIDE <- as.numeric(CleanedMergedData$TIDE)
CleanedMergedData$Exclusion <- as.numeric(CleanedMergedData$Exclusion)
CleanedMergedData$Dysfunction <- as.numeric(CleanedMergedData$Dysfunction)

### Chapter 4: Log2 transform variables of interest (MRC1 and CD163 expression data
# log2 transform variables of interest
# MRC1 (M2 Marker)
CleanedMergedData$MRC1 <- log(CleanedMergedData$MRC1+1, 2) # +1 to prevent log2(0)
# CD163 (M2 Marker)
CleanedMergedData$CD163 <- log(CleanedMergedData$CD163+1, 2) # +1 to prevent log2(0)

### Chapter 5: Establish ATC-High and ATC-Low CAF definitions (EPIC)

# Subset data to only ATC samples
ClinicalData_ATC <- CleanedMergedData %>% subset(Diagnosis == "ATC" & !is.na(CAF_EPIC)) 

# Create CAF quantiles from the ATC data for dividing into 50th percentile
ATC_Quantiles <- quantile(ClinicalData_ATC$CAF_EPIC)

# Group ATCs into > 50th percentile and < 50th percentile for CAFs
ClinicalData_ATC$CAF_Level_Binary <- "ATC"
for(i in 1:nrow(ClinicalData_ATC)){
  if(ClinicalData_ATC$CAF_EPIC[i] >= ATC_Quantiles[3]){
    ClinicalData_ATC$CAF_Level_Binary[i] <- "CAF-High"
  }
  else if(ClinicalData_ATC$CAF_EPIC[i] < ATC_Quantiles[3]){
    ClinicalData_ATC$CAF_Level_Binary[i] <- "CAF-Low"
  }
}

### Chapter 6: Make CAF-High vs CAF-Low plots of interest for Figure 6: 

## Plot Neutrophil_TIMER by CAF level for ATCs
plot <- ggplot(ClinicalData_ATC, aes(CAF_Level_Binary, Neutrophil_TIMER)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = CAF_Level_Binary, alpha = CAF_Level_Binary),
               show.legend = FALSE) + 
  scale_alpha_discrete(range = c(0.7, 0.7)) +
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 4, 
              alpha = 0.9,
              show.legend = FALSE) +
  labs (x = "ATC CAF Level", y = "Neutrophils TIMER") + 
  scale_fill_manual(values = c("green", "seagreen")) + # green = CAF-High; seagreen = CAF-Low
  theme_classic() + 
  ggtitle("Neutrophils TIMER") +
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 35), 
    axis.title.x = element_text(face = "bold", size = 35), 
    axis.title.y = element_blank(),
    axis.text = element_text(face = "bold", size = 30)) +
  scale_x_discrete(name ="ATC", limits = c("CAF-Low", "CAF-High")) # +
# geom_hline(yintercept = 0.0, linetype = 2, colour = "red")

# Save the Neutrophil_TIMER plot to an outputs folder
ggsave("outputs/Neutrophil_TIMER_by_CAFs_ATCs.png",
       width = 5.5,
       height = 5,
       plot, dpi = 600)
       
# Stats for Neutrophil_TIMER by CAF-Level: 
wilcox.test(data = ClinicalData_ATC, Neutrophil_TIMER ~ CAF_Level_Binary) 


## Plot Macrophage_TIMER by CAF level for ATCs
plot <- ggplot(ClinicalData_ATC, aes(CAF_Level_Binary, Macrophage_TIMER)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = CAF_Level_Binary, alpha = CAF_Level_Binary),
               show.legend = FALSE) + 
  scale_alpha_discrete(range = c(0.7, 0.7)) +
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 4, 
              alpha = 0.9,
              show.legend = FALSE) +
  labs (x = "ATC CAF Level", y = "Macrophage TIMER") + 
  scale_fill_manual(values = c("green", "seagreen")) + # green = CAF-High; seagreen = CAF-Low
  ggtitle("Macrophages TIMER") +
  theme_classic() + 
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 35), 
    axis.title.x = element_text(face = "bold", size = 35), 
    axis.title.y = element_blank(),
    axis.text = element_text(face = "bold", size = 30)) +
  scale_x_discrete(name ="ATC", limits = c("CAF-Low", "CAF-High")) #+
# geom_hline(yintercept = 0.0, linetype = 2, colour = "red")

# Save the TIMER Macrophage plot to an outputs folder
ggsave("outputs/Macrophage_TIMER_by_CAFs_ATCs.png",
       width = 5.5,
       height = 5,
       plot, dpi = 600)
       
# Stats for Macrophage_TIMER by CAF-Level: 
wilcox.test(data = ClinicalData_ATC, Macrophage_TIMER ~ CAF_Level_Binary) 


## Plot M1 Cibersort-Abs by CAF level for ATCs
plot <- ggplot(ClinicalData_ATC, aes(CAF_Level_Binary, M1)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = CAF_Level_Binary, alpha = CAF_Level_Binary),
               show.legend = FALSE) + 
  scale_alpha_discrete(range = c(0.7, 0.7)) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 4, 
              alpha = 0.9,
              show.legend = FALSE) +
  labs (x = "ATC CAF Level", y = "M1 Macrophages") + 
  scale_fill_manual(values = c("green", "seagreen")) + # green = CAF-High; seagreen = CAF-Low
  ggtitle("M1 CIBERSORT-Abs") +
  theme_classic() + 
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 34), 
    axis.title.x = element_text(face = "bold", size = 35), 
    axis.title.y = element_blank(),
    axis.text = element_text(face = "bold", size = 30)) +
  scale_x_discrete(name ="ATC", limits = c("CAF-Low", "CAF-High")) # +
# geom_hline(yintercept = 0.0, linetype = 2, colour = "red")

# Save the M1 Cibersort Abs plot to an outputs folder:
ggsave("outputs/M1_CibersortAbs_by_CAFs_ATCs.png",
       width = 5.5,
       height = 5,
       plot, dpi = 600)
       
# Stats for the M1 Cibersort Abs plot by CAF-Level: 
wilcox.test(data = ClinicalData_ATC, M1 ~ CAF_Level_Binary) 


## Plot M2 Cibersort-Abs by CAF level for ATCs
plot <- ggplot(ClinicalData_ATC, aes(CAF_Level_Binary, M2)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = CAF_Level_Binary, alpha = CAF_Level_Binary),
               show.legend = FALSE) + 
  scale_alpha_discrete(range = c(0.7, 0.7)) +
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 4, 
              alpha = 0.9,
              show.legend = FALSE) +
  labs (x = "ATC CAF Level", y = "M2 Macrophages") + 
  scale_fill_manual(values = c("green", "seagreen")) + # green = CAF-High; seagreen = CAF-Low
  ggtitle("M2 CIBERSORT-Abs") +
  theme_classic() + 
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 35), 
    axis.title.x = element_text(face = "bold", size = 35), 
    axis.title.y = element_blank(),
    axis.text = element_text(face = "bold", size = 30)) +
  scale_x_discrete(name ="ATC", limits = c("CAF-Low", "CAF-High")) # +
# geom_hline(yintercept = 0.0, linetype = 2, colour = "red")

# Save the M2 Cibersort Abs plot to an outputs folder:
ggsave("outputs/M2_CibersortAbs_by_CAFs_ATCs.png",
       width = 5.5,
       height = 5,
       plot, dpi = 600)
       
# Stats for the M2 Cibersort Abs plot by CAF-Level: 
wilcox.test(data = ClinicalData_ATC, M2 ~ CAF_Level_Binary) 


## Plot MRC1 Expression by CAF level for ATCs
plot <- ggplot(ClinicalData_ATC, aes(CAF_Level_Binary, MRC1)) +
  geom_violin(outlier.size = -1, 
               aes(fill = CAF_Level_Binary, alpha = CAF_Level_Binary),
               scale = "width",
               show.legend = FALSE) + 
  scale_alpha_discrete(range = c(0.7, 0.7)) +
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 4, 
              alpha = 0.9,
              show.legend = FALSE) +
  labs (x = "ATC CAF Level", y = "Dysfunction Score") + 
  scale_fill_manual(values = c("green", "seagreen")) + # green = CAF-High; seagreen = CAF-Low
  ggtitle("MRC1") +
  scale_y_continuous(breaks = c(8, 9, 10, 11, 12, 13, 14), limits = c(7.8, 14)) +
  theme_classic() + 
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 35), 
    axis.title.x = element_text(face = "bold", size = 35), 
    axis.title.y = element_blank(),
    axis.text = element_text(face = "bold", size = 30)) +
  scale_x_discrete(name ="ATC", limits = c("CAF-Low", "CAF-High"))
  #geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)
  
# Save the MRC1 by CAF-Level plot to an outputs folder:  
ggsave("outputs/MRC1_by_CAFs_ATCs.png",
       width = 5.5,
       height = 5,
       plot, dpi = 600)
       
# Stats for the MRC1 by CAF-Level plot:
wilcox.test(data = ClinicalData_ATC, MRC1 ~ CAF_Level_Binary) # This returned 7.069e-05 (equivalent w/ above)


## Plot CD163 Expression by CAF level for ATCs
plot <- ggplot(ClinicalData_ATC, aes(CAF_Level_Binary, CD163)) +
  geom_violin(outlier.size = -1, 
              aes(fill = CAF_Level_Binary, alpha = CAF_Level_Binary),
              scale = "width",
              show.legend = FALSE) + 
  scale_alpha_discrete(range = c(0.7, 0.7)) +
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 4, 
              alpha = 0.9,
              show.legend = FALSE) +
  labs (x = "ATC CAF Level", y = "Dysfunction Score") + 
  scale_fill_manual(values = c("green", "seagreen")) + # green = CAF-High; seagreen = CAF-Low
  ggtitle("CD163") +
  scale_y_continuous(breaks = c(9, 10, 11, 12, 13, 14), limits = c(8.8, 14.2)) +
  theme_classic() + 
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 35), 
    axis.title.x = element_text(face = "bold", size = 35), 
    axis.title.y = element_blank(),
    axis.text = element_text(face = "bold", size = 30)) +
  scale_x_discrete(name ="ATC", limits = c("CAF-Low", "CAF-High"))
#geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)

# Save the CD163 by CAF-Level plot to an outputs folder:  
ggsave("outputs/CD163_by_CAFs_ATCs.png",
       width = 5.5,
       height = 5,
       plot, dpi = 600)
       
# Stats for the CD163 by CAF-Level plot:
wilcox.test(data = ClinicalData_ATC, CD163 ~ CAF_Level_Binary)


## Plot TIDE Score by CAF-Level for ATCs
plot <- ggplot(ClinicalData_ATC, aes(CAF_Level_Binary, TIDE)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = CAF_Level_Binary, alpha = CAF_Level_Binary),
               show.legend = FALSE) + 
  scale_alpha_discrete(range = c(0.7, 0.7)) +
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 4, 
              alpha = 0.9,
              show.legend = FALSE) +
  labs (x = "ATC CAF Level", y = "TIDE Score") + 
  scale_fill_manual(values = c("green", "seagreen")) + # green = CAF-High; seagreen = CAF-Low
  ggtitle("TIDE Score") +
  theme_classic() + 
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 35), 
    axis.title.x = element_text(face = "bold", size = 35), 
    axis.title.y = element_blank(),
    axis.text = element_text(face = "bold", size = 30)) +
  scale_x_discrete(name ="ATC", limits = c("CAF-Low", "CAF-High")) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)
  
# Save the TIDE Score by CAF-Level plot to an outputs folder:   
ggsave("outputs/TIDE_Score_by_CAFs_ATCs.png",
       width = 5.5,
       height = 5,
       plot, dpi = 600)
       
# Stats for the TIDE Score by CAF-Level plot:
wilcox.test(data = ClinicalData_ATC, TIDE ~ CAF_Level_Binary) 


## Plot Exclusion Score  by CAF level for ATCs
plot <- ggplot(ClinicalData_ATC, aes(CAF_Level_Binary, Exclusion)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = CAF_Level_Binary, alpha = CAF_Level_Binary),
               show.legend = FALSE) + 
  scale_alpha_discrete(range = c(0.7, 0.7)) +
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 4, 
              alpha = 0.9,
              show.legend = FALSE) +
  labs (x = "ATC CAF Level", y = "Exclusion Score") + 
  scale_fill_manual(values = c("green", "seagreen")) + # green = CAF-High; seagreen = CAF-Low
  ggtitle("Exclusion Score") +
  theme_classic() + 
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 35), 
    axis.title.x = element_text(face = "bold", size = 35), 
    axis.title.y = element_blank(),
    axis.text = element_text(face = "bold", size = 30)) +
  scale_x_discrete(name ="ATC", limits = c("CAF-Low", "CAF-High")) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)
  
# Save the TIDE Exclusion Score by CAF-Level plot to an outputs folder:     
ggsave("outputs/Exclusion_Score_by_CAFs_ATCs.png",
       width = 5.5,
       height = 5,
       plot, dpi = 600)
       
# Stats for the TIDE Exlusion Score by CAF-Level plot:
wilcox.test(data = ClinicalData_ATC, Exclusion ~ CAF_Level_Binary) # This returned 5.086e-6 (equivalent w/ above)


## Plot Dysfunction Score by CAF level for ATCs
plot <- ggplot(ClinicalData_ATC, aes(CAF_Level_Binary, Dysfunction)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = CAF_Level_Binary, alpha = CAF_Level_Binary),
               show.legend = FALSE) + 
  scale_alpha_discrete(range = c(0.7, 0.7)) +
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 4, 
              alpha = 0.9,
              show.legend = FALSE) +
  labs (x = "ATC CAF Level", y = "Dysfunction Score") + 
  scale_fill_manual(values = c("green", "seagreen")) + # green = CAF-High; seagreen = CAF-Low
  ggtitle("Dysfunction Score") +
  theme_classic() + 
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 35), 
    axis.title.x = element_text(face = "bold", size = 35), 
    axis.title.y = element_blank(),
    axis.text = element_text(face = "bold", size = 30)) +
  scale_x_discrete(name ="ATC", limits = c("CAF-Low", "CAF-High")) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)
  
# Save the TIDE Dysfunction Score by CAF-Level plot to an outputs folder:   
ggsave("outputs/Dysfunction_Score_by_CAFs_ATCs.png",
       width = 5.5,
       height = 5,
       plot, dpi = 600)
       
# Stats for the TIDE Dysfunction Score by CAF-Level plot:
wilcox.test(data = ClinicalData_ATC, Dysfunction ~ CAF_Level_Binary) 

