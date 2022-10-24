### Chapters:
# 1. Load necessary packages
# 2. Load TCGA PTC Data (TIMER, Patient Data, Sample Data)
# 3. Clean and merge loaded TCGA data
# 4. Establish CAF-High and CAF-Low Groups for TCGA PTC Data
# 5. Looking at Specific TIMER Populations by CAF Level
# 6. Repeat above with PTCs from our cohort

# Note: was originally supplementary 12, so some of the outputs say Sup12 instead of Sup14

### Chapter 1: Load necessary packages ###
library(tidyverse)

### Chapter 2: Load TCGA PTC Data (TIMER, Patient Data, Sample Data) ###
# Load TIMER 2.0 TCGA Data
TCGA_Timer <- read_csv("data_in_use/infiltration_estimation_for_tcga.csv")

# Load TCGA Patient Data
TCGA_Patient_Data <- read.table(file = "data_in_use/TCGA_data_clinical_patient.txt", sep = '\t', header = TRUE) # Load patient data file
TCGA_Patient_Data <- as_tibble(TCGA_Patient_Data) # Make patient data into tibble for easier data cleaning

# Load TCGA Sample Data
TCGA_Sample_Data <- read.table(file = "data_in_use/TCGA_data_clinical_sample.txt", sep = '\t', header = TRUE) # Load sample data file
TCGA_Sample_Data <- as_tibble(TCGA_Sample_Data) # Make sample data into tibble for easier data cleaning


### Chapter 3: Clean and merge loaded TCGA data ###

# PATIEND_ID has an "01" appended to the end of each sample. Here, I will remove it from each sample. 
TCGA_Sample_Data$PATIENT_ID <- substr(TCGA_Sample_Data$PATIENT_ID, 1, nchar(TCGA_Sample_Data$PATIENT_ID)-3)
# TCGA_Sample_Data_Restricted$PATIENT_ID <- substr(TCGA_Sample_Data_Restricted$PATIENT_ID, 1, nchar(TCGA_Sample_Data_Restricted$PATIENT_ID)-3) # Old line of code from when I had a restricted file
# Combine sample data and patient data
Patient_Sample_Combined <- as_tibble(merge(TCGA_Sample_Data, TCGA_Patient_Data, by.x = "PATIENT_ID", by.y = "PATIENT_ID", all = T))

# Restrict TCGA_Timer to the patients in the Patient_Sample_Combined
TCGA_Timer <- TCGA_Timer %>% subset(cell_type %in% Patient_Sample_Combined$SAMPLE_ID)
TCGA_Timer <- TCGA_Timer %>% dplyr::rename("SAMPLE_ID" = "cell_type")

# Merge TCGA Timer data
TCGA_All_Data <- merge(TCGA_Timer, Patient_Sample_Combined)


### Chapter 4: Establish CAF-High and CAF-Low Groups for TCGA PTC Data ###
# Rename CAF prediction algorithm (EPIC)
TCGA_All_Data <- TCGA_All_Data %>% dplyr::rename("CAF_EPIC" = "Cancer associated fibroblast_EPIC")

# Log2 Transform CAF EPIC data
TCGA_All_Data$CAF_EPIC <- log(TCGA_All_Data$CAF_EPIC*100+1, 2)

# Establish CAF-Low and CAF-High definitions 
TCGA_Quantiles <- quantile(TCGA_All_Data$CAF_EPIC)

TCGA_All_Data$CAF_Level_Binary <- "Unknown"
for(i in 1:nrow(TCGA_All_Data)){
  if(TCGA_All_Data$CAF_EPIC[i] >= TCGA_Quantiles[3]){
    TCGA_All_Data$CAF_Level_Binary[i] <- "CAF-High"
  }
  else if(TCGA_All_Data$CAF_EPIC[i] < TCGA_Quantiles[3]){
    TCGA_All_Data$CAF_Level_Binary[i] <- "CAF-Low"
  }
}

### Chapter 5: Looking at Specific TIMER Populations by CAF Level ###
# TCGA Timer Macrophage by CAF_Level
plot <- ggplot(TCGA_All_Data, aes(CAF_Level_Binary, Macrophage_TIMER)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = CAF_Level_Binary),
               alpha = 0.4, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5, 
              alpha = 0.9,
              show.legend = FALSE) +
  labs (x = "TCGA PTCs", y = "TIMER Macrophages") + 
  # scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12), limits = c(2,12)) +
  scale_fill_manual(values = c("green", "seagreen")) + 
  theme_classic() + 
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 35), 
    axis.title.y = element_blank(),
    axis.text = element_text(face = "bold", size = 30)) +
  scale_x_discrete(name ="TCGA PTCs", limits = c("CAF-Low", "CAF-High")) +
  scale_y_continuous(limits = c(0, 0.45), breaks = c(0, 0.1, 0.2, 0.3, 0.4))
# geom_hline(yintercept = 0.0, linetype = 2, colour = "red")
ggsave("outputs/Sup12_Plots/TCGA_TIMER_Macrophage_by_CAF_Level.png",
       width = 7,
       height = 5,
       plot, dpi = 600)     
wilcox.test(data = TCGA_All_Data, Macrophage_TIMER ~ CAF_Level_Binary) # This output: P = 0.53

# TCGA Timer Neutrophil by CAF_Level
plot <- ggplot(TCGA_All_Data, aes(CAF_Level_Binary, Neutrophil_TIMER)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = CAF_Level_Binary),
               alpha = 0.4, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5, 
              alpha = 0.9,
              show.legend = FALSE) +
  labs (x = "TCGA PTCs", y = "TIMER Neutrophils") + 
  # scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12), limits = c(2,12)) +
  scale_fill_manual(values = c("green", "seagreen")) + 
  theme_classic() + 
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 35), 
    axis.title.y = element_blank(),
    axis.text = element_text(face = "bold", size = 30)) +
  scale_x_discrete(name ="TCGA PTCs", limits = c("CAF-Low", "CAF-High")) +
  scale_y_continuous(limits = c(0, 0.4), breaks = c(0, 0.1, 0.2, 0.3, 0.4))
# geom_hline(yintercept = 0.0, linetype = 2, colour = "red")
ggsave("outputs/Sup12_Plots/TIMER_Neutrophil_by_CAF_Level.png",
       width = 7,
       height = 5,
       plot, dpi = 600)
wilcox.test(data = TCGA_All_Data, Neutrophil_TIMER ~ CAF_Level_Binary) # This output: P < 2.2e-16

TCGA_All_Data <- TCGA_All_Data %>% dplyr::rename("M1" = "Macrophage M1_CIBERSORT-ABS")
TCGA_All_Data <- TCGA_All_Data %>% dplyr::rename("M2" = "Macrophage M2_CIBERSORT-ABS")

# TCGA Timer M1 Macrophage by CAF_Level
plot <- ggplot(TCGA_All_Data, aes(CAF_Level_Binary, M1)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = CAF_Level_Binary),
               alpha = 0.4, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5, 
              alpha = 0.9,
              show.legend = FALSE) +
  labs (x = "TCGA PTCs", y = "M1 Macrophages") + 
  # scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12), limits = c(2,12)) +
  scale_fill_manual(values = c("green", "seagreen")) + 
  theme_classic() + 
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 35), 
    axis.title.y = element_blank(),
    axis.text = element_text(face = "bold", size = 30)) +
  scale_x_discrete(name ="TCGA PTCs", limits = c("CAF-Low", "CAF-High")) +
  scale_y_continuous(limits = c(0, 0.15), breaks = c(0, 0.03, 0.06, 0.09, 0.12, 0.15))
# geom_hline(yintercept = 0.0, linetype = 2, colour = "red")
ggsave("outputs/Sup12_Plots/TCGA_CIBERSORT_M1_by_CAF_Level.png",
       width = 7,
       height = 5,
       plot, dpi = 600)
wilcox.test(data = TCGA_All_Data, M1 ~ CAF_Level_Binary) # This output: P = 3.074e-09

# M2 Macrophages
plot <- ggplot(TCGA_All_Data, aes(CAF_Level_Binary, M2)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = CAF_Level_Binary),
               alpha = 0.4, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5, 
              alpha = 0.9,
              show.legend = FALSE) +
  labs (x = "TCGA PTCs", y = "M2 Macrophages") + 
  # scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12), limits = c(2,12)) +
  scale_fill_manual(values = c("green", "seagreen")) + 
  theme_classic() + 
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 35), 
    axis.title.y = element_blank(),
    axis.text = element_text(face = "bold", size = 30)) +
  scale_x_discrete(name ="TCGA PTCs", limits = c("CAF-Low", "CAF-High")) +
  scale_y_continuous(limits = c(0, 0.45), breaks = c(0, 0.1, 0.2, 0.3, 0.4))
# geom_hline(yintercept = 0.0, linetype = 2, colour = "red")
ggsave("outputs/Sup12_Plots/TCGA_CIBERSORT_M2_by_CAF_Level.png",
       width = 7,
       height = 5,
       plot, dpi = 600)
wilcox.test(data = TCGA_All_Data, M2 ~ CAF_Level_Binary) # This output: < 2.2e-16

#### TIDE Integration ####
# Read in the TCGA TIDE File - Generated on 22-0601 on the TIDE website (http://tide.dfci.harvard.edu/)
TCGA_TIDE_Data <- read_csv(file = "data_in_use/22-0601_TCGA_Tide_Results.csv")

# Rename Patient -> SAMPLE_ID
TCGA_TIDE_Data <- TCGA_TIDE_Data %>% dplyr::rename("SAMPLE_ID" = "Patient")

# Switch . for - in TCGA_TIDE_Data

TCGA_TIDE_Data$SAMPLE_ID <- lapply(TCGA_TIDE_Data$SAMPLE_ID, gsub, pattern = ".", replacement = "-", fixed = TRUE)

TCGA_TIDE_Merge <- merge(TCGA_All_Data, TCGA_TIDE_Data)

### Exlusion, Dysfunction, TIDE Score plotting ###
# Now that TIDE data is integrated, I can plot these values for the TCGA

# TIDE Score
plot <- ggplot(TCGA_TIDE_Merge, aes(CAF_Level_Binary, TIDE)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = CAF_Level_Binary),
               alpha = 0.4, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5, 
              alpha = 0.9,
              show.legend = FALSE) +
  labs (x = "TCGA PTCs", y = "TIDE Score") + 
  # scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12), limits = c(2,12)) +
  scale_fill_manual(values = c("green", "seagreen")) + 
  theme_classic() + 
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 35), 
    axis.title.y = element_blank(),
    axis.text = element_text(face = "bold", size = 30)) +
  scale_x_discrete(name ="TCGA PTCs", limits = c("CAF-Low", "CAF-High")) +
  scale_y_continuous(limits = c(-2, 2.1), breaks = c(-2, -1, 0, 1, 2)) + 
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)
ggsave("outputs/Sup12_Plots/TCGA_TIDE_Score_by_CAF_Level.png",
       width = 7,
       height = 6,
       plot, dpi = 600)
wilcox.test(data = TCGA_TIDE_Merge, TIDE ~ CAF_Level_Binary) # This output: 2.1745e-15

# Dysfunction Score
plot <- ggplot(TCGA_TIDE_Merge, aes(CAF_Level_Binary, Dysfunction)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = CAF_Level_Binary),
               alpha = 0.4, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5, 
              alpha = 0.9,
              show.legend = FALSE) +
  labs (x = "TCGA PTCs", y = "Dysfunction Score") + 
  # scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12), limits = c(2,12)) +
  scale_fill_manual(values = c("green", "seagreen")) + 
  theme_classic() + 
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 35), 
    axis.title.y = element_blank(),
    axis.text = element_text(face = "bold", size = 30)) +
  scale_x_discrete(name ="TCGA PTCs", limits = c("CAF-Low", "CAF-High")) +
  scale_y_continuous(limits = c(-2.1, 2), breaks = c(-2, -1, 0, 1, 2)) + 
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)
ggsave("outputs/Sup12_Plots/TCGA_Dysfunction_Score_by_CAF_Level.png",
       width = 7,
       height = 6,
       plot, dpi = 600)
wilcox.test(data = TCGA_TIDE_Merge, Dysfunction ~ CAF_Level_Binary) # This output: <2.2e-16

# Exclusion Score
plot <- ggplot(TCGA_TIDE_Merge, aes(CAF_Level_Binary, Exclusion)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = CAF_Level_Binary),
               alpha = 0.4, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5, 
              alpha = 0.9,
              show.legend = FALSE) +
  labs (x = "TCGA PTCs", y = "Exclusion Score") + 
  # scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12), limits = c(2,12)) +
  scale_fill_manual(values = c("green", "seagreen")) + 
  theme_classic() + 
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 35), 
    axis.title.y = element_blank(),
    axis.text = element_text(face = "bold", size = 30)) +
  scale_x_discrete(name ="TCGA PTCs", limits = c("CAF-Low", "CAF-High")) +
  scale_y_continuous(limits = c(-2.15, 2.1), breaks = c(-2, -1, 0, 1, 2)) + 
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)
ggsave("outputs/Sup12_Plots/TCGA_Exclusion_Score_by_CAF_Level.png",
       width = 7,
       height = 6,
       plot, dpi = 600)
wilcox.test(data = TCGA_TIDE_Merge, Exclusion ~ CAF_Level_Binary) # This output: 1.299e-05

# Read in TCGA RSEM Z-Score RNA-sequencing data to add plots of specific markers of interest
TCGA_RNA_Z_Scores <- read.table(file = "data_in_use/data_RNA_Seq_v2_mRNA_median_Zscores.txt", sep = '\t', header = FALSE) # Read in TCGA RNA data in Z-Score format
Gene_Symbols <- TCGA_RNA_Z_Scores$V1 # Create variable for storing Gene Symbols -> will become column names
TCGA_RNA_Z_Scores <- TCGA_RNA_Z_Scores[,3:ncol(TCGA_RNA_Z_Scores)] # Delete first two columns, don't need these
TCGA_Transposed <- t(TCGA_RNA_Z_Scores) # Transpose the data so that each column pertains to a particular gene
colnames(TCGA_Transposed) <- Gene_Symbols # Make Gene Symbols the column names
TCGA_RNA_Tibble <- as_tibble(TCGA_Transposed) # Convert type to tibble
TCGA_RNA_Tibble <- dplyr::rename(TCGA_RNA_Tibble, SAMPLE_ID = Hugo_Symbol) # Rename "Hugo_Symbol" to "SAMPLE_ID" to match Patient_Sample_Combined for future merge
# TCGA_RNA_Tibble[,2:ncol(TCGA_RNA_Tibble)] <- as.numeric(unlist(TCGA_RNA_Tibble[,2:ncol(TCGA_RNA_Tibble)])) 
# Note the above command was very resource intensive...maybe it would be better to do this on the subset list of data of interest...also didn't work as I intended it to

# Cleaning up: remove the data I no longer need
remove(list = c("TCGA_RNA_Z_Scores", "TCGA_Transposed", "Gene_Symbols"))

# Merge data sets
TCGA_M2_Markers <- merge(TCGA_All_Data, TCGA_RNA_Tibble)

# Look at M2 Markers by CAF-Level
# MRC1 by CAF level
# First, change MRC1 to a numeric
TCGA_M2_Markers$MRC1 <- as.numeric(TCGA_M2_Markers$MRC1)
# Now, plot MRC1 by CAF level
# Now plot MRC1 by CAF level as a violin plot: 
plot <- ggplot(TCGA_M2_Markers, aes(CAF_Level_Binary, MRC1)) +
  geom_violin(outlier.size = -1, 
              aes(fill = CAF_Level_Binary),
              alpha = 0.4, 
              show.legend = FALSE,
              scale = "width") + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5, 
              alpha = 0.9,
              show.legend = FALSE) +
  labs (x = "TCGA PTCs", y = "MRC1 Z-Score") + 
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8), limits = c(-1,8.1)) +
  scale_fill_manual(values = c("green", "seagreen")) + 
  theme_classic() + 
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 35), 
    axis.title.y = element_blank(),
    axis.text = element_text(face = "bold", size = 30)) +
  scale_x_discrete(name ="TCGA PTCs", limits = c("CAF-Low", "CAF-High")) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)
ggsave("outputs/Sup12_Plots/TCGA_MRC1_Z_Score_by_CAF_Level_violin.png",
       width = 7,
       height = 6,
       plot, dpi = 600)

# Statistics for MRC1 by CAF level
wilcox.test(data = TCGA_M2_Markers, MRC1 ~ CAF_Level_Binary) # This output: < 2.2e-16

# CD163 by CAF level
# First, change CD163 to a numeric
TCGA_M2_Markers$CD163 <- as.numeric(TCGA_M2_Markers$CD163)
# Now, plot CD163 by CAF level
# CD163 by CAF level: violin plot
plot <- ggplot(TCGA_M2_Markers, aes(CAF_Level_Binary, CD163)) +
  geom_violin(outlier.size = -1, 
              aes(fill = CAF_Level_Binary),
              alpha = 0.4, 
              show.legend = FALSE,
              scale = "width") + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5, 
              alpha = 0.9,
              show.legend = FALSE) +
  labs (x = "TCGA PTCs", y = "CD163 Z-Score") + 
  scale_y_continuous(breaks = c(0, 2, 4, 6), limits = c(-1, 6.1)) +
  scale_fill_manual(values = c("green", "seagreen")) + 
  theme_classic() + 
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 35), 
    axis.title.y = element_blank(),
    axis.text = element_text(face = "bold", size = 30)) +
  scale_x_discrete(name ="TCGA PTCs", limits = c("CAF-Low", "CAF-High")) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)
ggsave("outputs/Sup12_Plots/TCGA_CD163_Z_Score_by_CAF_Level_Violin.png",
       width = 7,
       height = 6,
       plot, dpi = 600)

# CD163 by CAFs Statistics
wilcox.test(data = TCGA_M2_Markers, CD163 ~ CAF_Level_Binary) # This output: = 1.512e-14

remove(list = ls())

### Chapter 6: Now repeat with PTCs from our cohort ####
# Read in Cleaned Merged Data
CleanedMergedData <- readRDS(file = "data_in_use/CleanedMergedData_DESeq2NormalizedReads_Diagnosis_Simplified.rds")

## Establish PTC-High and PTC-Low CAF definitions (EPIC)
# Rename CAF EPIC
CleanedMergedData <- CleanedMergedData %>% dplyr::rename("CAF_EPIC" = "Cancer associated fibroblast_EPIC")
# Set CAF EPIC as a numeric
CleanedMergedData$CAF_EPIC <- as.numeric(CleanedMergedData$CAF_EPIC)

# Rename other variables of interest (e.g., macrophage variables)
CleanedMergedData <- CleanedMergedData %>% dplyr::rename("M1" = "Macrophage M1_CIBERSORT-ABS")
CleanedMergedData <- CleanedMergedData %>% dplyr::rename("M2" = "Macrophage M2_CIBERSORT-ABS")

# Set variables of interest to numerics
CleanedMergedData$BRS <- as.numeric(CleanedMergedData$BRS)
CleanedMergedData$Macrophage_TIMER <- as.numeric(CleanedMergedData$Macrophage_TIMER)
CleanedMergedData$M1 <- as.numeric(CleanedMergedData$M1)
CleanedMergedData$M2 <- as.numeric(CleanedMergedData$M2)
CleanedMergedData$Neutrophil_TIMER <- as.numeric(CleanedMergedData$Neutrophil_TIMER)
CleanedMergedData$TIDE <- as.numeric(CleanedMergedData$TIDE)
CleanedMergedData$Exclusion <- as.numeric(CleanedMergedData$Exclusion)
CleanedMergedData$Dysfunction <- as.numeric(CleanedMergedData$Dysfunction)

# log transform variables of interest
# CXCL5 (Ligand of interest)
CleanedMergedData$CXCL5 <- log(CleanedMergedData$CXCL5+1, 2)
# MRC1 (M2 Marker)
CleanedMergedData$MRC1 <- log(CleanedMergedData$MRC1+1, 2)
# CD163 (M2 Marker)
CleanedMergedData$CD163 <- log(CleanedMergedData$CD163+1, 2)
# CD80 (M1 Marker)
CleanedMergedData$CD80 <- log(CleanedMergedData$CD80+1, 2)
# CD86 (M1 Marker)
CleanedMergedData$CD86 <- log(CleanedMergedData$CD86+1, 2)
# FCGR1A (CD64, M1 Marker)
CleanedMergedData$FCGR1A <- log(CleanedMergedData$FCGR1A+1, 2)
# FCGR1B (CD64, M1 Marker)
CleanedMergedData$FCGR1B <- log(CleanedMergedData$FCGR1B+1, 2)
# FCGR2A (CD32, M1 Marker)
CleanedMergedData$FCGR2A <- log(CleanedMergedData$FCGR2A+1, 2)
# FCGR2B (CD32, M1 Marker)
CleanedMergedData$FCGR2B <- log(CleanedMergedData$FCGR2B+1, 2)
# CD40 (M1 Marker)
CleanedMergedData$CD40 <- log(CleanedMergedData$CD40+1, 2)
# NOS2 (M1 Marker)
CleanedMergedData$NOS2 <- log(CleanedMergedData$NOS2+1, 2)
# SOCS1 (M1 Marker)
CleanedMergedData$SOCS1 <- log(CleanedMergedData$SOCS1+1, 2)
# MARCO (M1 Marker)
CleanedMergedData$MARCO <- log(CleanedMergedData$MARCO+1, 2)
# CXCL9 (M1 Marker)
CleanedMergedData$CXCL9 <- log(CleanedMergedData$CXCL9+1, 2)
# CXCL10 (M1 Marker)
CleanedMergedData$CXCL10 <- log(CleanedMergedData$CXCL10+1, 2)
# CXCL11 (M1 Marker)
CleanedMergedData$CXCL11 <- log(CleanedMergedData$CXCL11+1, 2)
# IL1B (M1 Marker)
CleanedMergedData$IL1B <- log(CleanedMergedData$IL1B+1, 2)


# Subset data to create PTC data sets
ClinicalData_PTC_Local <- CleanedMergedData %>% subset(Diagnosis_Simplified == "PTC/\nIFVPTC" & !is.na(CAF_EPIC)) # Note that local here actually encompasses all samples
PTC_Quantiles <- quantile(ClinicalData_PTC_Local$CAF_EPIC)

# Group IFVPTCs/PTCs into > 50th percentile and < 50th percentile for CAFs
ClinicalData_PTC_Local$CAF_Level_Binary <- "PTC"
for(i in 1:nrow(ClinicalData_PTC_Local)){
  if(ClinicalData_PTC_Local$CAF_EPIC[i] >= PTC_Quantiles[3]){
    ClinicalData_PTC_Local$CAF_Level_Binary[i] <- "CAF-High"
  }
  else if(ClinicalData_PTC_Local$CAF_EPIC[i] < PTC_Quantiles[3]){
    ClinicalData_PTC_Local$CAF_Level_Binary[i] <- "CAF-Low"
  }
}

# Make the plots of interest for Figure 6: 
# Plot Neutrophil_TIMER by CAF level for PTCs
plot <- ggplot(ClinicalData_PTC_Local, aes(CAF_Level_Binary, Neutrophil_TIMER)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = CAF_Level_Binary),
               alpha = 0.2,
               show.legend = FALSE) + 
  scale_alpha_discrete(range = c(0.2, 0.7)) +
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 4, 
              alpha = 0.9,
              show.legend = FALSE) +
  labs (x = "PTC CAF Level", y = "Neutrophils TIMER") + 
  scale_fill_manual(values = c("green", "seagreen")) + 
  theme_classic() + 
  ggtitle("Neutrophils TIMER") +
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_blank(), 
    axis.title.x = element_text(face = "bold", size = 35), 
    axis.title.y = element_blank(),
    axis.text = element_text(face = "bold", size = 30)) +
  scale_x_discrete(name ="PTC/IFVPTC", limits = c("CAF-Low", "CAF-High")) +
  scale_y_continuous(limits = c(0, 0.25), breaks = c(0, 0.05, 0.1, 0.15, 0.2, 0.25))
# geom_hline(yintercept = 0.0, linetype = 2, colour = "red")
ggsave("outputs/22-0929_Sup12_Plots_All_Sample/22-0929_Neutrophil_TIMER_by_CAFs_PTCs.png",
       width = 7,
       height = 5,
       plot, dpi = 600)
# Stats: 
kruskal.test(Neutrophil_TIMER ~ CAF_Level_Binary, data = ClinicalData_PTC_Local) # shouldn't be necessary for only two groups
pairwise.wilcox.test(ClinicalData_PTC_Local$Neutrophil_TIMER, ClinicalData_PTC_Local$CAF_Level_Binary, # this returned 0.48
                     p.adjust.method = "BH")
wilcox.test(data = ClinicalData_PTC_Local, Neutrophil_TIMER ~ CAF_Level_Binary) # This returned 0.4809

# Plot Macrophage_TIMER by CAF level for PTCs
plot <- ggplot(ClinicalData_PTC_Local, aes(CAF_Level_Binary, Macrophage_TIMER)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = CAF_Level_Binary),
               alpha = 0.2,
               show.legend = FALSE) + 
  scale_alpha_discrete(range = c(0.2, 0.7)) +
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 4, 
              alpha = 0.9,
              show.legend = FALSE) +
  labs (x = "PTC CAF Level", y = "Macrophage TIMER") + 
  scale_fill_manual(values = c("green", "seagreen")) + 
  ggtitle("Macrophages TIMER") +
  theme_classic() + 
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_blank(), 
    axis.title.x = element_text(face = "bold", size = 35), 
    axis.title.y = element_blank(),
    axis.text = element_text(face = "bold", size = 30)) +
  scale_x_discrete(name ="PTC/IFVPTC", limits = c("CAF-Low", "CAF-High")) +
  scale_y_continuous(limits = c(0, 0.2), breaks = c(0, 0.05, 0.1, 0.15, 0.2))
# geom_hline(yintercept = 0.0, linetype = 2, colour = "red")
ggsave("outputs/22-0929_Sup12_Plots_All_Sample/22-0929_Macrophage_TIMER_by_CAFs_PTCs.png",
       width = 7,
       height = 5,
       plot, dpi = 600)
# Stats: 
kruskal.test(Macrophage_TIMER ~ CAF_Level_Binary, data = ClinicalData_PTC_Local) # shouldn't be necessary for only two groups
pairwise.wilcox.test(ClinicalData_PTC_Local$Macrophage_TIMER, ClinicalData_PTC_Local$CAF_Level_Binary, # this returned 0.34
                     p.adjust.method = "BH")
wilcox.test(data = ClinicalData_PTC_Local, Macrophage_TIMER ~ CAF_Level_Binary) # This returned 0.3393 (equivalent w/ above)

# Plot M1 Cibersort-Abs by CAF level for PTCs
plot <- ggplot(ClinicalData_PTC_Local, aes(CAF_Level_Binary, M1)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = CAF_Level_Binary),
               alpha = 0.2,
               show.legend = FALSE) + 
  scale_alpha_discrete(range = c(0.2, 0.7)) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 4, 
              alpha = 0.9,
              show.legend = FALSE) +
  labs (x = "PTC CAF Level", y = "M1 Macrophages") + 
  scale_fill_manual(values = c("green", "seagreen")) + 
  ggtitle("M1 CIBERSORT-Abs") +
  theme_classic() + 
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_blank(), 
    axis.title.x = element_text(face = "bold", size = 35), 
    axis.title.y = element_blank(),
    axis.text = element_text(face = "bold", size = 30)) +
  scale_x_discrete(name ="PTC/IFVPTC", limits = c("CAF-Low", "CAF-High")) +
  scale_y_continuous(limits = c(0, 0.07), breaks = c(0, 0.02, 0.04, 0.06))
# geom_hline(yintercept = 0.0, linetype = 2, colour = "red")
ggsave("outputs/22-0929_Sup12_Plots_All_Sample/22-0929_M1_CibersortAbs_by_CAFs_PTCs.png",
       width = 7,
       height = 5,
       plot, dpi = 600)
# Stats: 
kruskal.test(M1 ~ CAF_Level_Binary, data = ClinicalData_PTC_Local) # shouldn't be necessary for only two groups
pairwise.wilcox.test(ClinicalData_PTC_Local$M1, ClinicalData_PTC_Local$CAF_Level_Binary, # this returned 0.88 (< 0.0001)
                     p.adjust.method = "BH")
wilcox.test(data = ClinicalData_PTC_Local, M1 ~ CAF_Level_Binary) # This returned 0.8771 (equivalent w/ above)

# Plot M2 Cibersort-Abs by CAF level for PTCs
plot <- ggplot(ClinicalData_PTC_Local, aes(CAF_Level_Binary, M2)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = CAF_Level_Binary),
               alpha = 0.2,
               show.legend = FALSE) + 
  scale_alpha_discrete(range = c(0.2, 0.7)) +
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 4, 
              alpha = 0.9,
              show.legend = FALSE) +
  labs (x = "PTC CAF Level", y = "M2 Macrophages") + 
  scale_fill_manual(values = c("green", "seagreen")) + 
  ggtitle("M2 CIBERSORT-Abs") +
  theme_classic() + 
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_blank(), 
    axis.title.x = element_text(face = "bold", size = 35), 
    axis.title.y = element_blank(),
    axis.text = element_text(face = "bold", size = 30)) +
  scale_x_discrete(name ="PTC/IFVPTC", limits = c("CAF-Low", "CAF-High")) +
  scale_y_continuous(limits = c(0, 0.15), breaks = c(0, 0.03, 0.06, 0.09, 0.12, 0.15))
# geom_hline(yintercept = 0.0, linetype = 2, colour = "red")
ggsave("outputs/22-0929_Sup12_Plots_All_Sample/22-0929_M2_CibersortAbs_by_CAFs_PTCs.png",
       width = 7,
       height = 5,
       plot, dpi = 600)
# Stats: 
kruskal.test(M2 ~ CAF_Level_Binary, data = ClinicalData_PTC_Local) # shouldn't be necessary for only two groups
pairwise.wilcox.test(ClinicalData_PTC_Local$M2, ClinicalData_PTC_Local$CAF_Level_Binary, # this returned 0.38
                     p.adjust.method = "BH")
wilcox.test(data = ClinicalData_PTC_Local, M2 ~ CAF_Level_Binary) # This returned 0.3776 (equivalent w/ above)

# Plot TIDE Score by CAF level for PTCs
plot <- ggplot(ClinicalData_PTC_Local, aes(CAF_Level_Binary, TIDE)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = CAF_Level_Binary),
               alpha = 0.2,
               show.legend = FALSE) + 
  scale_alpha_discrete(range = c(0.2, 0.7)) +
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 4, 
              alpha = 0.9,
              show.legend = FALSE) +
  labs (x = "PTC CAF Level", y = "TIDE Score") + 
  scale_fill_manual(values = c("green", "seagreen")) + 
  ggtitle("TIDE Score") +
  theme_classic() + 
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_blank(), 
    axis.title.x = element_text(face = "bold", size = 35), 
    axis.title.y = element_blank(),
    axis.text = element_text(face = "bold", size = 30)) +
  scale_x_discrete(name ="PTC/IFVPTC", limits = c("CAF-Low", "CAF-High")) +
  scale_y_continuous(limits = c(-2,2), breaks = c(-2, -1, 0, 1, 2)) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)
ggsave("outputs/22-0929_Sup12_Plots_All_Sample/22-0929_TIDE_Score_by_CAFs_PTCs.png",
       width = 7,
       height = 6,
       plot, dpi = 600)
# Stats: 
kruskal.test(TIDE ~ CAF_Level_Binary, data = ClinicalData_PTC_Local) # shouldn't be necessary for only two groups
pairwise.wilcox.test(ClinicalData_PTC_Local$TIDE, ClinicalData_PTC_Local$CAF_Level_Binary, # this returned 0.072
                     p.adjust.method = "BH")
wilcox.test(data = ClinicalData_PTC_Local, TIDE ~ CAF_Level_Binary) # This returned 0.07202 (equivalent w/ above)

# Plot Exclusion Score  by CAF level for PTCs
plot <- ggplot(ClinicalData_PTC_Local, aes(CAF_Level_Binary, Exclusion)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = CAF_Level_Binary),
               alpha = 0.2,
               show.legend = FALSE) + 
  scale_alpha_discrete(range = c(0.2, 0.7)) +
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 4, 
              alpha = 0.9,
              show.legend = FALSE) +
  labs (x = "PTC CAF Level", y = "Exclusion Score") + 
  scale_fill_manual(values = c("green", "seagreen")) + 
  ggtitle("Exclusion Score") +
  theme_classic() + 
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_blank(), 
    axis.title.x = element_text(face = "bold", size = 35), 
    axis.title.y = element_blank(),
    axis.text = element_text(face = "bold", size = 30)) +
  scale_x_discrete(name ="PTC/IFVPTC", limits = c("CAF-Low", "CAF-High")) +
  scale_y_continuous(limits = c(-2, 2), breaks = c(-2, -1, 0, 1, 2)) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)
ggsave("outputs/22-0929_Sup12_Plots_All_Sample/22-0929_Exclusion_Score_by_CAFs_PTCs.png",
       width = 7,
       height = 6,
       plot, dpi = 600)
# Stats: 
kruskal.test(Exclusion ~ CAF_Level_Binary, data = ClinicalData_PTC_Local) # shouldn't be necessary for only two groups
pairwise.wilcox.test(ClinicalData_PTC_Local$Exclusion, ClinicalData_PTC_Local$CAF_Level_Binary, # this returned 0.068
                     p.adjust.method = "BH")
wilcox.test(data = ClinicalData_PTC_Local, Exclusion ~ CAF_Level_Binary) # This returned 0.06792 (equivalent w/ above)

# Plot Dysfunction Score by CAF level for PTCs
plot <- ggplot(ClinicalData_PTC_Local, aes(CAF_Level_Binary, Dysfunction)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = CAF_Level_Binary),
               alpha = 0.2, 
               show.legend = FALSE) + 
  scale_alpha_discrete(range = c(0.2, 0.7)) +
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 4, 
              alpha = 0.9,
              show.legend = FALSE) +
  labs (x = "PTC CAF Level", y = "Dysfunction Score") + 
  scale_fill_manual(values = c("green", "seagreen")) + 
  ggtitle("Dysfunction Score") +
  theme_classic() + 
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_blank(), 
    axis.title.x = element_text(face = "bold", size = 35), 
    axis.title.y = element_blank(),
    axis.text = element_text(face = "bold", size = 30)) +
  scale_x_discrete(name ="PTC/IFVPTC", limits = c("CAF-Low", "CAF-High")) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2) +
  scale_y_continuous(limits = c(-2, 2), breaks = c(-2, -1, 0, 1, 2))
ggsave("outputs/22-0929_Sup12_Plots_All_Sample/22-0929_Dysfunction_Score_by_CAFs_PTCs.png",
       width = 7,
       height = 6,
       plot, dpi = 600)
# Stats: 
kruskal.test(Dysfunction ~ CAF_Level_Binary, data = ClinicalData_PTC_Local) # shouldn't be necessary for only two groups
pairwise.wilcox.test(ClinicalData_PTC_Local$Dysfunction, ClinicalData_PTC_Local$CAF_Level_Binary, # this returned 0.029
                     p.adjust.method = "BH")
wilcox.test(data = ClinicalData_PTC_Local, Dysfunction ~ CAF_Level_Binary) # This returned 0.02865 (equivalent w/ above)

# Plot MRC1 Expression by CAF level for PTCs
plot <- ggplot(ClinicalData_PTC_Local, aes(CAF_Level_Binary, MRC1)) +
  geom_violin(outlier.size = -1, 
              aes(fill = CAF_Level_Binary),
              alpha = 0.2, 
              scale = "width",
              show.legend = FALSE) + 
  scale_alpha_discrete(range = c(0.2, 0.7)) +
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 4, 
              alpha = 0.9,
              show.legend = FALSE) +
  labs (x = "PTC CAF Level", y = "Dysfunction Score") + 
  scale_fill_manual(values = c("green", "seagreen")) + 
  ggtitle("MRC1") +
  scale_y_continuous(breaks = c(5, 7, 9, 11), limits = c(4.9, 11.9)) +
  theme_classic() + 
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_blank(), 
    axis.title.x = element_text(face = "bold", size = 35), 
    axis.title.y = element_blank(),
    axis.text = element_text(face = "bold", size = 30)) +
  scale_x_discrete(name ="PTC/IFVPTC", limits = c("CAF-Low", "CAF-High"))
#geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)
ggsave("outputs/22-0929_Sup12_Plots_All_Sample/22-0929_MRC1_by_CAFs_PTCs.png",
       width = 7,
       height = 6,
       plot, dpi = 600)
# Stats: 
kruskal.test(MRC1 ~ CAF_Level_Binary, data = ClinicalData_PTC_Local) # shouldn't be necessary for only two groups
pairwise.wilcox.test(ClinicalData_PTC_Local$MRC1, ClinicalData_PTC_Local$CAF_Level_Binary, # this returned 6.9e-05
                     p.adjust.method = "BH")
wilcox.test(data = ClinicalData_PTC_Local, MRC1 ~ CAF_Level_Binary) # This returned 6.943e-05 (equivalent w/ above)

# Plot CD163 Expression by CAF level for PTCs
plot <- ggplot(ClinicalData_PTC_Local, aes(CAF_Level_Binary, CD163)) +
  geom_violin(outlier.size = -1, 
              aes(fill = CAF_Level_Binary),
              alpha = 0.2, 
              scale = "width",
              show.legend = FALSE) + 
  scale_alpha_discrete(range = c(0.2, 0.7)) +
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 4, 
              alpha = 0.9,
              show.legend = FALSE) +
  labs (x = "PTC CAF Level", y = "Dysfunction Score") + 
  scale_fill_manual(values = c("green", "seagreen")) + 
  ggtitle("CD163") +
  scale_y_continuous(breaks = c(4, 6, 8, 10, 12), limits = c(3.3, 12.6)) +
  theme_classic() + 
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_blank(), 
    axis.title.x = element_text(face = "bold", size = 35), 
    axis.title.y = element_blank(),
    axis.text = element_text(face = "bold", size = 30)) +
  scale_x_discrete(name ="PTC/IFVPTC", limits = c("CAF-Low", "CAF-High"))
#geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)
ggsave("outputs/Sup12_Plots/CD163_by_CAFs_PTCs.png",
       width = 7,
       height = 6,
       plot, dpi = 600)
# Stats: 
wilcox.test(data = ClinicalData_PTC_Local, CD163 ~ CAF_Level_Binary) # This returned 0.2871 
