# Author: Matthew Aaron Loberg
# 22-0531 TCGA Data Exploration

# Looking at CAF Levels and CAF Markers by mutation status, BRAF-RAS Status, and histology


# Load required packages
library(tidyverse)

#### SETUP ####

# Read in TIMER 2.0 TCGA Data
TCGA_Timer <- read_csv("data_in_use/infiltration_estimation_for_tcga.csv")

# Read in list of BRAF Mutant Patients/Samples
BRAF_Patient_ID <- read.table(file = "data_in_use/ListOfBRAFMutantPatients_with_RNA-Seq.tsv", sep = '\t', header = TRUE)

# Read in list of RAS mutant Patients/Samples
RAS_Patient_ID <- read.table(file = "data_in_use/ListOfRASMutantPatients_with_RNA-Seq.tsv", sep = '\t', header = TRUE)


# The following is from 22-0426 Script: 
# Load Data
TCGA_Patient_Data <- read.table(file = "data_in_use/22-0107_TCGA_data_clinical_patient.txt", sep = '\t', header = TRUE) # Load patient data file
TCGA_Patient_Data <- as_tibble(TCGA_Patient_Data) # Make patient data into tibble for easier data cleaning

TCGA_Sample_Data <- read.table(file = "data_in_use/22-0107_TCGA_data_clinical_sample.txt", sep = '\t', header = TRUE) # Load sample data file
TCGA_Sample_Data <- as_tibble(TCGA_Sample_Data) # Make sample data into tibble for easier data cleaning

# Restrict data to the columns needed
# This probably isn't entirely necessary, but I'm doing it to familiarize myself with the variables and think about the things that would be worthwhile to look at
# For patient data, I want the following columns: Patient ID, Tumor Status, Overall survival (OS) status, OS months, disease free status (DFS), DFS months
# TCGA_Patient_Data_Restricted <- TCGA_Patient_Data[c("PATIENT_ID", "TUMOR_STATUS", "OS_STATUS", "OS_MONTHS", "DFS_STATUS", "DFS_MONTHS")]
# 22-0405 update: Not sure yet what columns I want. I need to think of which columns might be predictors of aggressive disease.
# For now working with all columns/non-restricted file

# For sample data, I want the following columns:
# 22-0405 update: Not sure yet what columns I want. I need to think of which columns might be predictors of aggressive disease.
# TCGA_Sample_Data_Restricted <- TCGA_Sample_Data[c("PATIENT_ID", "SAMPLE_ID", "BRAFV600E_RAS", "EXTRATHYROIDAL_EXTENSION", "DIFFERENTIATION_SCORE", "ERK_SCORE", "DISEASE_STAGE", "BRAFV600E_RAF_SCORE", "RNASEQ_DATA", "SAMPLE_TYPE")]
# For now working with all columns/non-restricted file

# PATIEND_ID has an "01" appended to the end of each sample. Here, I will remove it from each sample. 
TCGA_Sample_Data$PATIENT_ID <- substr(TCGA_Sample_Data$PATIENT_ID, 1, nchar(TCGA_Sample_Data$PATIENT_ID)-3)
# TCGA_Sample_Data_Restricted$PATIENT_ID <- substr(TCGA_Sample_Data_Restricted$PATIENT_ID, 1, nchar(TCGA_Sample_Data_Restricted$PATIENT_ID)-3) # Old line of code from when I had a restricted file
# Combine sample data and patient data
Patient_Sample_Combined <- as_tibble(merge(TCGA_Sample_Data, TCGA_Patient_Data, by.x = "PATIENT_ID", by.y = "PATIENT_ID", all = T))

# Restrict TCGA_Timer to the patients in the Patient_Sample_Combined
TCGA_Timer <- TCGA_Timer %>% subset(cell_type %in% Patient_Sample_Combined$SAMPLE_ID)
TCGA_Timer <- TCGA_Timer %>% dplyr::rename("SAMPLE_ID" = "cell_type")

TCGA_All_Data <- merge(TCGA_Timer, Patient_Sample_Combined)
TCGA_All_Data$BRAF_RAS_Mutant_Status <- "No BRAF/RAS\nMutation"
for(i in 1:nrow(TCGA_All_Data)){
  if(TCGA_All_Data$SAMPLE_ID[i] %in% RAS_Patient_ID$Sample.ID){
    TCGA_All_Data$BRAF_RAS_Mutant_Status[i] <- "RAS\nMutant"
  }
  else if(TCGA_All_Data$SAMPLE_ID[i] %in% BRAF_Patient_ID$Sample.ID){
    TCGA_All_Data$BRAF_RAS_Mutant_Status[i] <- "BRAF\nMutant"
  }
}

TCGA_All_Data <- TCGA_All_Data %>% dplyr::rename("CAF_EPIC" = "Cancer associated fibroblast_EPIC")
TCGA_All_Data <- TCGA_All_Data %>% dplyr::rename("CAF_MCPCOUNTER" = "Cancer associated fibroblast_MCPCOUNTER")

TCGA_All_Data$CAF_EPIC <- log(TCGA_All_Data$CAF_EPIC*100+1, 2)
TCGA_All_Data$CAF_MCPCOUNTER <- log(TCGA_All_Data$CAF_MCPCOUNTER+1, 2)

#### Plots by Mutation Status ####

# CAFs by BRAF-RAS Mutation Status, EPIC
plot <- ggplot(TCGA_All_Data, aes(BRAF_RAS_Mutant_Status, CAF_EPIC)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = BRAF_RAS_Mutant_Status),
               alpha = 0.4, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5, 
              alpha = 0.9,
              show.legend = FALSE) +
  labs (x = "TCGA PTCs", y = "log2(CAF EPIC)") + 
  scale_fill_manual(values = c("red", "grey", "blue")) + 
  theme_classic() + 
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 30), 
    axis.title.y = element_text(face = "bold", size = 30),
    axis.text = element_text(face = "bold", size = 23)) +
  scale_x_discrete(name ="TCGA PTCs", limits = c("RAS\nMutant", "No BRAF/RAS\nMutation", "BRAF\nMutant")) #+
  # geom_hline(yintercept = 0.0, linetype = 2, colour = "red")
ggsave("outputs/22-0531_TCGA_TIMER_Exploration/22-0531_CAF_EPIC_by_Mutation_Status.png",
       width = 7,
       height = 5,
       plot, dpi = 600)
# CAFs by BRAF-RAS Mutation Status, EPIC: violin
plot <- ggplot(TCGA_All_Data, aes(BRAF_RAS_Mutant_Status, CAF_EPIC)) +
  geom_violin(outlier.size = -1, 
               aes(fill = BRAF_RAS_Mutant_Status),
               alpha = 0.2, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5, 
              alpha = 0.9,
              show.legend = FALSE) +
  labs (x = "TCGA PTCs", y = "CAF EPIC") + 
  scale_fill_manual(values = c("red", "grey", "blue")) + 
  theme_classic() + 
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 30), 
    axis.title.y = element_text(face = "bold", size = 30),
    axis.text = element_text(face = "bold", size = 23)) +
  scale_x_discrete(name ="TCGA PTCs", limits = c("RAS\nMutant", "No BRAF/RAS\nMutation", "BRAF\nMutant")) #+
# geom_hline(yintercept = 0.0, linetype = 2, colour = "red")
ggsave("outputs/22-0531_TCGA_TIMER_Exploration/22-0531_CAF_EPIC_by_Mutation_Status_Violin.png",
       width = 7,
       height = 5,
       plot, dpi = 600)
# Stats: 
kruskal.test(CAF_EPIC ~ BRAF_RAS_Mutant_Status, data = TCGA_All_Data) # 
pairwise.wilcox.test(TCGA_All_Data$CAF_EPIC, TCGA_All_Data$BRAF_RAS_Mutant_Status, # 
                     p.adjust.method = "BH")
# BRAF vs RAS: < 2e-16; BRAF vs none: 1.9e-08; RAS vs none: 1.1e-08
# wilcox.test(data = TCGA_All_Data, CAF_EPIC ~ BRAF_RAS_Mutant_Status) # This test not valid w/ > 2 groups

# CAFs by BRAF-RAS Mutation Status, MCPCOUNTER
plot <- ggplot(TCGA_All_Data, aes(BRAF_RAS_Mutant_Status, CAF_MCPCOUNTER)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = BRAF_RAS_Mutant_Status),
               alpha = 0.4, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5, 
              alpha = 0.9,
              show.legend = FALSE) +
  labs (x = "TCGA PTCs", y = "log2\n(CAF MCPCOUNTER)") + 
  scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12), limits = c(2,12)) +
  scale_fill_manual(values = c("red", "grey", "blue")) + 
  theme_classic() + 
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 30), 
    axis.title.y = element_text(face = "bold", size = 27),
    axis.text = element_text(face = "bold", size = 23)) +
  scale_x_discrete(name ="TCGA PTCs", limits = c("RAS\nMutant", "No BRAF/RAS\nMutation", "BRAF\nMutant")) #+
# geom_hline(yintercept = 0.0, linetype = 2, colour = "red")
ggsave("outputs/22-0531_TCGA_TIMER_Exploration/22-0531_CAF_MCPCOUNTER_by_Mutation_Status.png",
       width = 7,
       height = 5,
       plot, dpi = 600)
# CAFs by BRAF-RAS Mutation Status, MCPCOUNTER: violin plot
plot <- ggplot(TCGA_All_Data, aes(BRAF_RAS_Mutant_Status, CAF_MCPCOUNTER)) +
  geom_violin(outlier.size = -1, 
               aes(fill = BRAF_RAS_Mutant_Status),
               alpha = 0.2, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5, 
              alpha = 0.9,
              show.legend = FALSE) +
  labs (x = "TCGA PTCs", y = "CAF MCPCOUNTER") + 
  scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12), limits = c(2,12)) +
  scale_fill_manual(values = c("red", "grey", "blue")) + 
  theme_classic() + 
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 30), 
    axis.title.y = element_text(face = "bold", size = 29.5),
    axis.text = element_text(face = "bold", size = 23)) +
  scale_x_discrete(name ="TCGA PTCs", limits = c("RAS\nMutant", "No BRAF/RAS\nMutation", "BRAF\nMutant")) #+
# geom_hline(yintercept = 0.0, linetype = 2, colour = "red")
ggsave("outputs/22-0531_TCGA_TIMER_Exploration/22-0531_CAF_MCPCOUNTER_by_Mutation_Status_Violin.png",
       width = 7,
       height = 5,
       plot, dpi = 600)
# Stats: 
kruskal.test(CAF_MCPCOUNTER ~ BRAF_RAS_Mutant_Status, data = TCGA_All_Data) # shouldn't be necessary for only two groups
pairwise.wilcox.test(TCGA_All_Data$CAF_MCPCOUNTER, TCGA_All_Data$BRAF_RAS_Mutant_Status, # 
                     p.adjust.method = "BH")
# BRAF vs RAS: < 1.1e-09; BRAF vs none: 0.0086; RAS vs none: 4.2e-07
# wilcox.test(data = TCGA_All_Data, CAF_MCPCOUNTER ~ BRAF_RAS_Mutant_Status) # Not relevant w/ > 2 groups

#### End: Plots by BRAF-RAS Mutation Status ####

#### Start: Plots by BRAF-RAS Status #####

for(i in 1:nrow(TCGA_All_Data)){
  if(TCGA_All_Data$BRAFV600E_RAS[i] == "Braf-like"){
    TCGA_All_Data$BRAFV600E_RAS[i] <- "BRAF-Like"
  }
  else if(TCGA_All_Data$BRAFV600E_RAS[i] == "Ras-like"){
    TCGA_All_Data$BRAFV600E_RAS[i] <- "RAS-Like"
  }
}
TCGA_BRS_Data <- TCGA_All_Data %>% subset(BRAFV600E_RAS != "")
# CAFs by BRAF-RAS Score Status, EPIC
plot <- ggplot(TCGA_BRS_Data, aes(BRAFV600E_RAS, CAF_EPIC)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = BRAFV600E_RAS),
               alpha = 0.4, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5, 
              alpha = 0.9,
              show.legend = FALSE) +
  labs (x = "TCGA PTCs", y = "log2(CAF EPIC)") + 
  # scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12), limits = c(2,12)) + # did this for MCPCOUNTER axis
  scale_fill_manual(values = c("red", "blue")) + 
  theme_classic() + 
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 30), 
    axis.title.y = element_text(face = "bold", size = 30),
    axis.text = element_text(face = "bold", size = 25)) +
  scale_x_discrete(name ="TCGA PTCs", limits = c("RAS-Like", "BRAF-Like")) #+
# geom_hline(yintercept = 0.0, linetype = 2, colour = "red")
ggsave("outputs/22-0531_TCGA_TIMER_Exploration/22-0531_CAF_EPIC_by_BRAF_RAS_Status.png",
       width = 7,
       height = 5,
       plot, dpi = 600)
wilcox.test(data = TCGA_BRS_Data, CAF_EPIC ~ BRAFV600E_RAS) # This output: 2.2e-16

# CAFs by BRAF-RAS Score Status, MCPCOUNTER
plot <- ggplot(TCGA_BRS_Data, aes(BRAFV600E_RAS, CAF_MCPCOUNTER)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = BRAFV600E_RAS),
               alpha = 0.4, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5, 
              alpha = 0.9,
              show.legend = FALSE) +
  labs (x = "TCGA PTCs", y = "log2\n(CAF MCPCOUNTER)") + 
  scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12), limits = c(2,12)) +
  scale_fill_manual(values = c("red", "blue")) + 
  theme_classic() + 
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 30), 
    axis.title.y = element_text(face = "bold", size = 29),
    axis.text = element_text(face = "bold", size = 25)) +
  scale_x_discrete(name ="TCGA PTCs", limits = c("RAS-Like", "BRAF-Like")) #+
# geom_hline(yintercept = 0.0, linetype = 2, colour = "red")
ggsave("outputs/22-0531_TCGA_TIMER_Exploration/22-0531_CAF_MCPCOUNTER_by_BRAF_RAS_Status.png",
       width = 7,
       height = 5,
       plot, dpi = 600)
wilcox.test(data = TCGA_BRS_Data, CAF_MCPCOUNTER ~ BRAFV600E_RAS) # This output: 2.734e-13

#### Histology Setup ####
# I pulled the following block of code from 22-0426_TCGA_Survival_526_BRAF_PoorOutcome.R lines 565-597
# It takes the TCGA data and makes histologic groups for plotting 
TCGA_All_Data$HISTOLOGICAL_TYPE_COMPLEX <- TCGA_All_Data$HISTOLOGICAL_TYPE
for(i in 1:nrow(TCGA_All_Data)){
  if(TCGA_All_Data$HISTOLOGICAL_TYPE[i] == ""){
    TCGA_All_Data$HISTOLOGICAL_TYPE_COMPLEX[i] <- "Unknown\nNot Specified"
  }
  else if(TCGA_All_Data$HISTOLOGICAL_TYPE_OTHER[i] == "Diffuse Sclerosing Variant" |
          TCGA_All_Data$HISTOLOGICAL_TYPE_OTHER[i] == "Diffuse sclerosing variant" |
          TCGA_All_Data$HISTOLOGICAL_TYPE_OTHER[i] == "Papillary carcinoma, diffuse sclerosing variant"){
    TCGA_All_Data$HISTOLOGICAL_TYPE_COMPLEX[i] <- "Diffuse\nSclerosing"
  }
  else if(TCGA_All_Data$HISTOLOGICAL_TYPE_OTHER[i] == "cribriform morular" |
          TCGA_All_Data$HISTOLOGICAL_TYPE_OTHER[i] == "cribiform morular"){
    TCGA_All_Data$HISTOLOGICAL_TYPE_COMPLEX[i] <- "Cribriform\nMorular"
  }
  else if(TCGA_All_Data$HISTOLOGICAL_TYPE_OTHER[i] == "encapsulated follicular"){
    TCGA_All_Data$HISTOLOGICAL_TYPE_COMPLEX[i] <- "Follicular" 
  }
  else if(TCGA_All_Data$HISTOLOGICAL_TYPE[i] == "Tall Cell"){
    TCGA_All_Data$HISTOLOGICAL_TYPE_COMPLEX[i] <- "Tall\nCell"
  }
}

# Subset the data to just the following groups: 
# Follicular 
# Diffuse\nSclerosing
# Tall Cell
# Classical
Histology_Plot_Cohort <- TCGA_All_Data %>% subset(HISTOLOGICAL_TYPE_COMPLEX == "Classical" |
                                                  HISTOLOGICAL_TYPE_COMPLEX == "Diffuse\nSclerosing" |
                                                  HISTOLOGICAL_TYPE_COMPLEX == "Follicular" | 
                                                  HISTOLOGICAL_TYPE_COMPLEX == "Tall\nCell")
# End of line 597 from 22-0426

# This is also from the same R file but lines 616-623
### Make a new plot combining diffuse sclerosing and tall cell
Histology_Plot_Cohort$HISTOLOGICAL_TYPE_SIMPLIFIED <- Histology_Plot_Cohort$HISTOLOGICAL_TYPE_COMPLEX
for(i in 1:nrow(Histology_Plot_Cohort)){
  if(Histology_Plot_Cohort$HISTOLOGICAL_TYPE_COMPLEX[i] == "Diffuse\nSclerosing" |
     Histology_Plot_Cohort$HISTOLOGICAL_TYPE_COMPLEX[i] == "Tall\nCell"){
    Histology_Plot_Cohort$HISTOLOGICAL_TYPE_SIMPLIFIED[i] <- "Diffuse\nSclerosing/\nTall Cell"
  }
}

#### Plots by Histology ####
# Histology simplified plot: EPIC CAFs
plot <- ggplot(Histology_Plot_Cohort, aes(HISTOLOGICAL_TYPE_SIMPLIFIED, CAF_EPIC)) + 
  geom_boxplot(outlier.size = -1,
               aes(fill = HISTOLOGICAL_TYPE_SIMPLIFIED),
               alpha = 0.4,
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5,
              alpha = c(0.8),
              show.legend = FALSE) +
  # scale_fill_manual(values = c("#F8766D", "#C77CFF", "#00BFC4")) + 
  scale_x_discrete(name = "Histology", limits = c("Follicular", "Classical", "Diffuse\nSclerosing/\nTall Cell")) + 
  scale_fill_manual(values = c("red", "purple", "blue")) +
  labs(y = "log2(CAF EPIC)") + 
  theme_classic() + 
  # geom_hline(yintercept = 0.0, linetype = 2, colour = "red") + 
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 30), 
    axis.title.y = element_text(face = "bold", size = 30),
    axis.text = element_text(face = "bold", size = 25))
plot
ggsave("outputs/22-0531_TCGA_TIMER_Exploration/22-0531_CAF_EPIC_by_Histology.png",
       width = 7,
       height = 5,
       plot, dpi = 600)
kruskal.test(CAF_EPIC ~ HISTOLOGICAL_TYPE_SIMPLIFIED, data = Histology_Plot_Cohort) # 
pairwise.wilcox.test(Histology_Plot_Cohort$CAF_EPIC, Histology_Plot_Cohort$HISTOLOGICAL_TYPE_SIMPLIFIED, # 
                     p.adjust.method = "BH")

# Histology simplified plot: MCPCOUNTER CAFs
plot <- ggplot(Histology_Plot_Cohort, aes(HISTOLOGICAL_TYPE_SIMPLIFIED, CAF_MCPCOUNTER)) + 
  geom_boxplot(outlier.size = -1,
               aes(fill = HISTOLOGICAL_TYPE_SIMPLIFIED),
               alpha = 0.4,
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5,
              alpha = c(0.8),
              show.legend = FALSE) +
  # scale_fill_manual(values = c("#F8766D", "#C77CFF", "#00BFC4")) + 
  scale_x_discrete(name = "Histology", limits = c("Follicular", "Classical", "Diffuse\nSclerosing/\nTall Cell")) + 
  scale_fill_manual(values = c("red", "purple", "blue")) +
  labs(y = "log2\n(CAF MCPCOUNTER)") + 
  theme_classic() + 
  scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12), limits = c(2,12)) +
  # geom_hline(yintercept = 0.0, linetype = 2, colour = "red") + 
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 30), 
    axis.title.y = element_text(face = "bold", size = 24),
    axis.text = element_text(face = "bold", size = 25))
plot
ggsave("outputs/22-0531_TCGA_TIMER_Exploration/22-0531_CAF_MCPCOUNTER_by_Histology.png",
       width = 7,
       height = 5,
       plot, dpi = 600)
kruskal.test(CAF_MCPCOUNTER ~ HISTOLOGICAL_TYPE_SIMPLIFIED, data = Histology_Plot_Cohort) # 
pairwise.wilcox.test(Histology_Plot_Cohort$CAF_MCPCOUNTER, Histology_Plot_Cohort$HISTOLOGICAL_TYPE_SIMPLIFIED, # 
                     p.adjust.method = "BH")


#### CAF-High vs CAF-Low Defining + Plotting ####
TCGA_Quantiles <- quantile(TCGA_All_Data$CAF_EPIC)

TCGA_All_Data$CAF_Level_Binary <- "Unknown"
for(i in 1:nrow(TCGA_All_Data)){
  if(TCGA_All_Data$CAF_EPIC[i] >= TCGA_Quantiles[3]){
    TCGA_All_Data$CAF_Level_Binary[i] <- "High-CAF"
  }
  else if(TCGA_All_Data$CAF_EPIC[i] < TCGA_Quantiles[3]){
    TCGA_All_Data$CAF_Level_Binary[i] <- "Low-CAF"
  }
}

# Plot for showing 50th percentile of CAFs
TCGA_All_Data$Unifier <- "Unified"
plot <- ggplot(TCGA_All_Data, aes(Unifier, CAF_EPIC)) +
  annotate("rect", xmin = -Inf, xmax = Inf, 
                             ymin = -Inf, ymax = TCGA_Quantiles[3],
            alpha = 0.4,
            fill = "blue") +
  annotate("rect", xmin = -Inf, xmax = Inf,
                   ymin = TCGA_Quantiles[3], ymax = Inf,
           alpha = 0.4,
           fill = "red") +
  geom_boxplot(outlier.size = -1, 
               aes(),
               alpha = 0.2, 
               show.legend = FALSE,
               lwd = 2) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 3, 
              alpha = 0.9,
              show.legend = FALSE) +
  labs (x = "TCGA PTCs", y = "log2(CAF EPIC)") + 
  # scale_fill_manual(values = c("red")) + 
  theme_classic() + 
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 35), 
    axis.title.y = element_text(face = "bold", size = 35),
    axis.text.y = element_text(face = "bold", size = 25),
    axis.text.x = element_blank()) +
  # scale_x_discrete(name ="PTC", limits = c("PTC")) +
  geom_hline(yintercept = TCGA_Quantiles[3], linetype = 2, colour = "red", size = 3)
ggsave("outputs/22-0531_TCGA_TIMER_Exploration/22-0531_TCGA_CAFs_50th_Percentile_Line.png",
       width = 7,
       height = 5,
       plot, dpi = 600)

# Plot for showing 50th percentile of CAFs w/ alternate alpha
TCGA_All_Data$Unifier <- "Unified"
plot <- ggplot(TCGA_All_Data, aes(Unifier, CAF_EPIC)) +
  annotate("rect", xmin = -Inf, xmax = Inf, 
           ymin = -Inf, ymax = TCGA_Quantiles[3],
           alpha = 0.2,
           fill = "blue") +
  annotate("rect", xmin = -Inf, xmax = Inf,
           ymin = TCGA_Quantiles[3], ymax = Inf,
           alpha = 0.2,
           fill = "red") +
  geom_boxplot(outlier.size = -1, 
               aes(),
               alpha = 0.2, 
               show.legend = FALSE,
               lwd = 2) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 3, 
              alpha = 0.9,
              show.legend = FALSE) +
  labs (x = "TCGA PTCs", y = "log2(CAF EPIC)") + 
  # scale_fill_manual(values = c("red")) + 
  theme_classic() + 
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 35), 
    axis.title.y = element_text(face = "bold", size = 35),
    axis.text.y = element_text(face = "bold", size = 25),
    axis.text.x = element_blank()) +
  # scale_x_discrete(name ="PTC", limits = c("PTC")) +
  geom_hline(yintercept = TCGA_Quantiles[3], linetype = 2, colour = "red", size = 3)
ggsave("outputs/22-0531_TCGA_TIMER_Exploration/22-0531_TCGA_CAFs_50th_alpha_2.png",
       width = 7,
       height = 5,
       plot, dpi = 600)

# 50th percentile plot w/ alternate color scheme
plot <- ggplot(TCGA_All_Data, aes(Unifier, CAF_EPIC)) +
  annotate("rect", xmin = -Inf, xmax = Inf, 
           ymin = -Inf, ymax = TCGA_Quantiles[3],
           alpha = 0.7,
           fill = "darkorchid4") +
  annotate("rect", xmin = -Inf, xmax = Inf,
           ymin = TCGA_Quantiles[3], ymax = Inf,
           alpha = 0.4,
           fill = "magenta") +
  geom_boxplot(outlier.size = -1, 
               aes(),
               alpha = 0.4, 
               show.legend = FALSE,
               lwd = 2) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 3, 
              alpha = 0.9,
              show.legend = FALSE) +
  labs (x = "TCGA PTCs", y = "log2(CAF EPIC)") + 
  # scale_fill_manual(values = c("red")) + 
  theme_classic() + 
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 35), 
    axis.title.y = element_text(face = "bold", size = 35),
    axis.text.y = element_text(face = "bold", size = 25),
    axis.text.x = element_blank()) +
  # scale_x_discrete(name ="PTC", limits = c("PTC")) +
  geom_hline(yintercept = TCGA_Quantiles[3], linetype = 2, colour = "red", size = 3)
ggsave("outputs/22-0531_TCGA_TIMER_Exploration/22-0531_TCGA_CAFs_50th_Percentile_Line_Orange_Cyan_Alternate.png",
       width = 7,
       height = 5,
       plot, dpi = 600)

#### Looking at Specific TIMER Populations by CAF Level ####
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
  scale_fill_manual(values = c("red", "blue")) + 
  theme_classic() + 
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 30), 
    axis.title.y = element_text(face = "bold", size = 29),
    axis.text = element_text(face = "bold", size = 25)) +
  scale_x_discrete(name ="TCGA PTCs", limits = c("Low-CAF", "High-CAF")) #+
# geom_hline(yintercept = 0.0, linetype = 2, colour = "red")
ggsave("outputs/22-0531_TCGA_TIMER_Exploration/22-0531_TCGA_TIMER_Macrophage_by_CAF_Level.png",
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
  scale_fill_manual(values = c("red", "blue")) + 
  theme_classic() + 
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 30), 
    axis.title.y = element_text(face = "bold", size = 29),
    axis.text = element_text(face = "bold", size = 25)) +
  scale_x_discrete(name ="TCGA PTCs", limits = c("Low-CAF", "High-CAF")) #+
# geom_hline(yintercept = 0.0, linetype = 2, colour = "red")
ggsave("outputs/22-0531_TCGA_TIMER_Exploration/22-0531_TCGA_TIMER_Neutrophil_by_CAF_Level.png",
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
  scale_fill_manual(values = c("red", "blue")) + 
  theme_classic() + 
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 30), 
    axis.title.y = element_text(face = "bold", size = 29),
    axis.text = element_text(face = "bold", size = 25)) +
  scale_x_discrete(name ="TCGA PTCs", limits = c("Low-CAF", "High-CAF")) #+
# geom_hline(yintercept = 0.0, linetype = 2, colour = "red")
ggsave("outputs/22-0531_TCGA_TIMER_Exploration/22-0531_TCGA_CIBERSORT_M1_by_CAF_Level.png",
       width = 7,
       height = 5,
       plot, dpi = 600)
wilcox.test(data = TCGA_All_Data, M1 ~ CAF_Level_Binary) # This output: P = 3.74e-09

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
  scale_fill_manual(values = c("red", "blue")) + 
  theme_classic() + 
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 30), 
    axis.title.y = element_text(face = "bold", size = 29),
    axis.text = element_text(face = "bold", size = 25)) +
  scale_x_discrete(name ="TCGA PTCs", limits = c("Low-CAF", "High-CAF")) #+
# geom_hline(yintercept = 0.0, linetype = 2, colour = "red")
ggsave("outputs/22-0531_TCGA_TIMER_Exploration/22-0531_TCGA_CIBERSORT_M2_by_CAF_Level.png",
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
  scale_fill_manual(values = c("red", "blue")) + 
  theme_classic() + 
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 30), 
    axis.title.y = element_text(face = "bold", size = 29),
    axis.text = element_text(face = "bold", size = 25)) +
  scale_x_discrete(name ="TCGA PTCs", limits = c("Low-CAF", "High-CAF")) +
geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)
ggsave("outputs/22-0531_TCGA_TIMER_Exploration/22-0531_TCGA_TIDE_Score_by_CAF_Level.png",
       width = 7,
       height = 5,
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
  scale_fill_manual(values = c("red", "blue")) + 
  theme_classic() + 
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 30), 
    axis.title.y = element_text(face = "bold", size = 29),
    axis.text = element_text(face = "bold", size = 25)) +
  scale_x_discrete(name ="TCGA PTCs", limits = c("Low-CAF", "High-CAF")) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)
ggsave("outputs/22-0531_TCGA_TIMER_Exploration/22-0531_TCGA_Dysfunction_Score_by_CAF_Level.png",
       width = 7,
       height = 5,
       plot, dpi = 600)
wilcox.test(data = TCGA_TIDE_Merge, Dysfunction ~ CAF_Level_Binary) # This output: 2.2e-16

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
  scale_fill_manual(values = c("red", "blue")) + 
  theme_classic() + 
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 30), 
    axis.title.y = element_text(face = "bold", size = 29),
    axis.text = element_text(face = "bold", size = 25)) +
  scale_x_discrete(name ="TCGA PTCs", limits = c("Low-CAF", "High-CAF")) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)
ggsave("outputs/22-0531_TCGA_TIMER_Exploration/22-0531_TCGA_Exclusion_Score_by_CAF_Level.png",
       width = 7,
       height = 5,
       plot, dpi = 600)
wilcox.test(data = TCGA_TIDE_Merge, Exclusion ~ CAF_Level_Binary) # This output: 1.299e-05


#### Adding in RNA-Seq data to look specifically at M2 Markers ####
# Borrowing code for reading in sequencing data from 22-0426_TCGA_Survival_526_BRAF_PoorOutcome.R
# Specifically, lines 76-88

## old code
# Read in TCGA RNA Data (Z-Score format) and clean the data to the desired format
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

## New code
# Merge data sets
TCGA_M2_Markers <- merge(TCGA_All_Data, TCGA_RNA_Tibble)

# Look at M2 Markers by CAF-Level
# MRC1 by CAF level
# First, change MRC1 to a numeric
TCGA_M2_Markers$MRC1 <- as.numeric(TCGA_M2_Markers$MRC1)
# Now, plot MRC1 by CAF level
plot <- ggplot(TCGA_M2_Markers, aes(CAF_Level_Binary, MRC1)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = CAF_Level_Binary),
               alpha = 0.4, 
               show.legend = FALSE,
               lwd = 2) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5, 
              alpha = 0.9,
              show.legend = FALSE) +
  labs (x = "TCGA PTCs", y = "MRC1 Z-Score") + 
  # scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12), limits = c(2,12)) +
  scale_fill_manual(values = c("red", "blue")) + 
  theme_classic() + 
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 30), 
    axis.title.y = element_text(face = "bold", size = 29),
    axis.text = element_text(face = "bold", size = 25)) +
  scale_x_discrete(name ="TCGA PTCs", limits = c("Low-CAF", "High-CAF")) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)
ggsave("outputs/22-0531_TCGA_TIMER_Exploration/22-0531_TCGA_MRC1_Z_Score_by_CAF_Level.png",
       width = 7,
       height = 5,
       plot, dpi = 600)

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
  # scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12), limits = c(2,12)) +
  scale_fill_manual(values = c("red", "blue")) + 
  theme_classic() + 
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 30), 
    axis.title.y = element_text(face = "bold", size = 29),
    axis.text = element_text(face = "bold", size = 25)) +
  scale_x_discrete(name ="TCGA PTCs", limits = c("Low-CAF", "High-CAF")) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)
ggsave("outputs/22-0531_TCGA_TIMER_Exploration/22-0531_TCGA_MRC1_Z_Score_by_CAF_Level_violin.png",
       width = 7,
       height = 5,
       plot, dpi = 600)

# Statistics for MRC1 by CAF level
wilcox.test(data = TCGA_M2_Markers, MRC1 ~ CAF_Level_Binary) # This output: < 2.2e-16

# CD163 by CAF level
# First, change CD163 to a numeric
TCGA_M2_Markers$CD163 <- as.numeric(TCGA_M2_Markers$CD163)
# Now, plot CD163 by CAF level
plot <- ggplot(TCGA_M2_Markers, aes(CAF_Level_Binary, CD163)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = CAF_Level_Binary),
               alpha = 0.4, 
               show.legend = FALSE,
               lwd = 2) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5, 
              alpha = 0.9,
              show.legend = FALSE) +
  labs (x = "TCGA PTCs", y = "CD163 Z-Score") + 
  # scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12), limits = c(2,12)) +
  scale_fill_manual(values = c("red", "blue")) + 
  theme_classic() + 
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 30), 
    axis.title.y = element_text(face = "bold", size = 29),
    axis.text = element_text(face = "bold", size = 25)) +
  scale_x_discrete(name ="TCGA PTCs", limits = c("Low-CAF", "High-CAF")) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)
ggsave("outputs/22-0531_TCGA_TIMER_Exploration/22-0531_TCGA_CD163_Z_Score_by_CAF_Level.png",
       width = 7,
       height = 5,
       plot, dpi = 600)

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
  # scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12), limits = c(2,12)) +
  scale_fill_manual(values = c("red", "blue")) + 
  theme_classic() + 
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 30), 
    axis.title.y = element_text(face = "bold", size = 29),
    axis.text = element_text(face = "bold", size = 25)) +
  scale_x_discrete(name ="TCGA PTCs", limits = c("Low-CAF", "High-CAF")) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)
ggsave("outputs/22-0531_TCGA_TIMER_Exploration/22-0531_TCGA_CD163_Z_Score_by_CAF_Level_Violin.png",
       width = 7,
       height = 5,
       plot, dpi = 600)

# CD163 by CAFs Statistics
wilcox.test(data = TCGA_M2_Markers, CD163 ~ CAF_Level_Binary) # This output: = 1.512e-14

### End of M2 data: Remove data frames/tibbles added for this section ###
remove(TCGA_M2_Markers, TCGA_RNA_Tibble)

#### Use BRAF-RAS Status and Histological Status to look at FAP/ACTA-2 Level ####
# For this, load in RNA data again 
## old code
# Read in TCGA RNA Data (Z-Score format) and clean the data to the desired format
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

# Merge
TCGA_BRS <- merge(TCGA_All_Data, TCGA_RNA_Tibble)
TCGA_BRS <- TCGA_BRS %>% subset(BRAFV600E_RAS != "")

# Plot FAP, ACTA2 by BRAF-RAS Status
# FAP by BRS Status
# First, change FAP to a numeric
TCGA_BRS$FAP <- as.numeric(TCGA_BRS$FAP)
# Now, plot FAP by BRS Status
plot <- ggplot(TCGA_BRS, aes(BRAFV600E_RAS, FAP)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = BRAFV600E_RAS),
               alpha = 0.4, 
               show.legend = FALSE,
               lwd = 2) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 2, 
              alpha = 0.9,
              show.legend = FALSE) +
  labs (x = "TCGA PTCs", y = "FAP Z-Score") + 
  # scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12), limits = c(2,12)) +
  scale_fill_manual(values = c("red", "blue")) + 
  theme_classic() + 
  ylim(-1,5) +
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 30), 
    axis.title.y = element_text(face = "bold", size = 29),
    axis.text = element_text(face = "bold", size = 25)) +
  scale_x_discrete(name ="TCGA PTCs", limits = c("RAS-Like", "BRAF-Like")) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)

# Now, plot FAP by BRS Status as a violin plot
plot <- ggplot(TCGA_BRS, aes(BRAFV600E_RAS, FAP)) +
  geom_violin(outlier.size = -1, 
               aes(fill = BRAFV600E_RAS),
               alpha = 0.4, 
               show.legend = FALSE,
               lwd = 1,
               scale = "width") + 
  #geom_jitter(aes(),
              #position = position_jitter(width = 0.1, height = 0),
              #size = 1, 
              #alpha = 0.3,
              #show.legend = FALSE) +
  labs (x = "TCGA PTCs", y = "FAP Z-Score") + 
  # scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12), limits = c(2,12)) +
  scale_fill_manual(values = c("red", "blue")) + 
  theme_classic() + 
  ylim(-1,5) +
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 30), 
    axis.title.y = element_text(face = "bold", size = 29),
    axis.text = element_text(face = "bold", size = 25)) +
  scale_x_discrete(name ="TCGA PTCs", limits = c("RAS-Like", "BRAF-Like")) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)
ggsave("outputs/22-0531_TCGA_TIMER_Exploration/22-0531_TCGA_FAP_Z_Score_by_BRS_Status_Violin.png",
       width = 7,
       height = 5,
       plot, dpi = 600)

# FAP by BRS Status Statistics
wilcox.test(data = TCGA_BRS, FAP ~ BRAFV600E_RAS) # This output: < 2.2e-16

# ACTA2 by BRS Status
# First, change ACTA2 to a numeric
TCGA_BRS$ACTA2 <- as.numeric(TCGA_BRS$ACTA2)
# Now, plot ACTA2 by BRS Status
plot <- ggplot(TCGA_BRS, aes(BRAFV600E_RAS, ACTA2)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = BRAFV600E_RAS),
               alpha = 0.4, 
               show.legend = FALSE,
               lwd = 2) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 2, 
              alpha = 0.9,
              show.legend = FALSE) +
  labs (x = "TCGA PTCs", y = "ACTA2 Z-Score") + 
  # scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12), limits = c(2,12)) +
  scale_fill_manual(values = c("red", "blue")) + 
  theme_classic() + 
  ylim(-2,5) +
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 30), 
    axis.title.y = element_text(face = "bold", size = 29),
    axis.text = element_text(face = "bold", size = 25)) +
  scale_x_discrete(name ="TCGA PTCs", limits = c("RAS-Like", "BRAF-Like")) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)
ggsave("outputs/22-0531_TCGA_TIMER_Exploration/22-0531_TCGA_ACTA2_Z_Score_by_BRS_Status.png",
       width = 7,
       height = 5,
       plot, dpi = 600)

# Now, plot ACTA2 by BRS Status as a violin plot
plot <- ggplot(TCGA_BRS, aes(BRAFV600E_RAS, ACTA2)) +
  geom_violin(outlier.size = -1, 
               aes(fill = BRAFV600E_RAS),
               alpha = 0.4, 
               show.legend = FALSE,
               lwd = 1,
               scale = "width") + 
  # geom_jitter(aes(),
              # position = position_jitter(width = 0.1, height = 0),
              # size = 1, 
              # alpha = 0.3,
              # show.legend = FALSE) +
  labs (x = "TCGA PTCs", y = "ACTA2 Z-Score") + 
  # scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12), limits = c(2,12)) +
  scale_fill_manual(values = c("red", "blue")) + 
  theme_classic() + 
  ylim(-2,5) +
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 30), 
    axis.title.y = element_text(face = "bold", size = 29),
    axis.text = element_text(face = "bold", size = 25)) +
  scale_x_discrete(name ="TCGA PTCs", limits = c("RAS-Like", "BRAF-Like")) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)
ggsave("outputs/22-0531_TCGA_TIMER_Exploration/22-0531_TCGA_ACTA2_Z_Score_by_BRS_Status_Violin.png",
       width = 7,
       height = 5,
       plot, dpi = 600)

# ACTA2 by BRS Status Statistics
wilcox.test(data = TCGA_BRS, ACTA2 ~ BRAFV600E_RAS) # This output: = 0.5844

## myCAF and iCAF markers by BRAF-RAS Status
# myCAF: POSTN, WNT2, COL1A1, LRRC15
# iCAF: PDGFRA, CXCL12
# myCAF First

# POSTN by BRS Status
# First, change POSTN to a numeric
TCGA_BRS$POSTN <- as.numeric(TCGA_BRS$POSTN)
# Now, plot POSTN by BRS Status
plot <- ggplot(TCGA_BRS, aes(BRAFV600E_RAS, POSTN)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = BRAFV600E_RAS),
               alpha = 0.4, 
               show.legend = FALSE,
               lwd = 1) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 2, 
              alpha = 0.9,
              show.legend = FALSE) +
  labs (x = "TCGA PTCs", y = "POSTN Z-Score") + 
  # scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12), limits = c(2,12)) +
  scale_fill_manual(values = c("red", "blue")) + 
  theme_classic() + 
  ylim(-1,5) +
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 30), 
    axis.title.y = element_text(face = "bold", size = 29),
    axis.text = element_text(face = "bold", size = 25)) +
  scale_x_discrete(name ="TCGA PTCs", limits = c("RAS-Like", "BRAF-Like")) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)
ggsave("outputs/22-0531_TCGA_TIMER_Exploration/22-0531_TCGA_POSTN_Z_Score_by_BRS_Status.png",
       width = 7,
       height = 5,
       plot, dpi = 600)

# Now, plot POSTN by BRS Status as a violin plot
plot <- ggplot(TCGA_BRS, aes(BRAFV600E_RAS, POSTN)) +
  geom_violin(outlier.size = -1, 
              aes(fill = BRAFV600E_RAS),
              alpha = 0.4, 
              show.legend = FALSE,
              lwd = 1,
              scale = "width") + 
  # geom_jitter(aes(),
  # position = position_jitter(width = 0.1, height = 0),
  # size = 1, 
  # alpha = 0.3,
  # show.legend = FALSE) +
  labs (x = "TCGA PTCs", y = "POSTN Z-Score") + 
  # scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12), limits = c(2,12)) +
  scale_fill_manual(values = c("red", "blue")) + 
  theme_classic() + 
  ylim(-1,5) +
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 30), 
    axis.title.y = element_text(face = "bold", size = 29),
    axis.text = element_text(face = "bold", size = 25)) +
  scale_x_discrete(name ="TCGA PTCs", limits = c("RAS-Like", "BRAF-Like")) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)
ggsave("outputs/22-0531_TCGA_TIMER_Exploration/22-0531_TCGA_POSTN_Z_Score_by_BRS_Status_Violin.png",
       width = 7,
       height = 5,
       plot, dpi = 600)

# POSTN by BRS Status Statistics
wilcox.test(data = TCGA_BRS, POSTN ~ BRAFV600E_RAS) # This output: P = 7.171e-9

# WNT2 by BRS Status
# First, change WNT2 to a numeric
TCGA_BRS$WNT2 <- as.numeric(TCGA_BRS$WNT2)
# Now, plot WNT2 by BRS Status
plot <- ggplot(TCGA_BRS, aes(BRAFV600E_RAS, WNT2)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = BRAFV600E_RAS),
               alpha = 0.4, 
               show.legend = FALSE,
               lwd = 1) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 2, 
              alpha = 0.9,
              show.legend = FALSE) +
  labs (x = "TCGA PTCs", y = "WNT2 Z-Score") + 
  # scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12), limits = c(2,12)) +
  scale_fill_manual(values = c("red", "blue")) + 
  theme_classic() + 
  ylim(-1,5) +
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 30), 
    axis.title.y = element_text(face = "bold", size = 29),
    axis.text = element_text(face = "bold", size = 25)) +
  scale_x_discrete(name ="TCGA PTCs", limits = c("RAS-Like", "BRAF-Like")) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)
ggsave("outputs/22-0531_TCGA_TIMER_Exploration/22-0531_TCGA_WNT2_Z_Score_by_BRS_Status.png",
       width = 7,
       height = 5,
       plot, dpi = 600)

# Now, plot WNT2 by BRS Status as a violin plot
plot <- ggplot(TCGA_BRS, aes(BRAFV600E_RAS, WNT2)) +
  geom_violin(outlier.size = -1, 
              aes(fill = BRAFV600E_RAS),
              alpha = 0.4, 
              show.legend = FALSE,
              lwd = 1,
              scale = "width") + 
  # geom_jitter(aes(),
  # position = position_jitter(width = 0.1, height = 0),
  # size = 1, 
  # alpha = 0.3,
  # show.legend = FALSE) +
  labs (x = "TCGA PTCs", y = "WNT2 Z-Score") + 
  # scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12), limits = c(2,12)) +
  scale_fill_manual(values = c("red", "blue")) + 
  theme_classic() + 
  ylim(-1,5) +
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 30), 
    axis.title.y = element_text(face = "bold", size = 29),
    axis.text = element_text(face = "bold", size = 25)) +
  scale_x_discrete(name ="TCGA PTCs", limits = c("RAS-Like", "BRAF-Like")) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)
ggsave("outputs/22-0531_TCGA_TIMER_Exploration/22-0531_TCGA_WNT2_Z_Score_by_BRS_Status_Violin.png",
       width = 7,
       height = 5,
       plot, dpi = 600)

# WNT2 by BRS Status Statistics
wilcox.test(data = TCGA_BRS, WNT2 ~ BRAFV600E_RAS) # This output: P <2.2e-16

# COL1A1 by BRS Status
# First, change COL1A1 to a numeric
TCGA_BRS$COL1A1 <- as.numeric(TCGA_BRS$COL1A1)
# Now, plot COL1A1 by BRS Status
plot <- ggplot(TCGA_BRS, aes(BRAFV600E_RAS, COL1A1)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = BRAFV600E_RAS),
               alpha = 0.4, 
               show.legend = FALSE,
               lwd = 1) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 2, 
              alpha = 0.9,
              show.legend = FALSE) +
  labs (x = "TCGA PTCs", y = "COL1A1 Z-Score") + 
  # scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12), limits = c(2,12)) +
  scale_fill_manual(values = c("red", "blue")) + 
  theme_classic() + 
  ylim(-1,5) +
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 30), 
    axis.title.y = element_text(face = "bold", size = 29),
    axis.text = element_text(face = "bold", size = 25)) +
  scale_x_discrete(name ="TCGA PTCs", limits = c("RAS-Like", "BRAF-Like")) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)
ggsave("outputs/22-0531_TCGA_TIMER_Exploration/22-0531_TCGA_COL1A1_Z_Score_by_BRS_Status.png",
       width = 7,
       height = 5,
       plot, dpi = 600)

# Now, plot COL1A1 by BRS Status as a violin plot
plot <- ggplot(TCGA_BRS, aes(BRAFV600E_RAS, COL1A1)) +
  geom_violin(outlier.size = -1, 
              aes(fill = BRAFV600E_RAS),
              alpha = 0.4, 
              show.legend = FALSE,
              lwd = 1,
              scale = "width") + 
  # geom_jitter(aes(),
  # position = position_jitter(width = 0.1, height = 0),
  # size = 1, 
  # alpha = 0.3,
  # show.legend = FALSE) +
  labs (x = "TCGA PTCs", y = "COL1A1 Z-Score") + 
  # scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12), limits = c(2,12)) +
  scale_fill_manual(values = c("red", "blue")) + 
  theme_classic() + 
  ylim(-1,5) +
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 30), 
    axis.title.y = element_text(face = "bold", size = 29),
    axis.text = element_text(face = "bold", size = 25)) +
  scale_x_discrete(name ="TCGA PTCs", limits = c("RAS-Like", "BRAF-Like")) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)
ggsave("outputs/22-0531_TCGA_TIMER_Exploration/22-0531_TCGA_COL1A1_Z_Score_by_BRS_Status_Violin.png",
       width = 7,
       height = 5,
       plot, dpi = 600)

# COL1A1 by BRS Status Statistics
wilcox.test(data = TCGA_BRS, COL1A1 ~ BRAFV600E_RAS) # This output: P <2.2e-16

# LRRC15 by BRS Status
# First, change LRRC15 to a numeric
TCGA_BRS$LRRC15 <- as.numeric(TCGA_BRS$LRRC15)
# Now, plot LRRC15 by BRS Status
plot <- ggplot(TCGA_BRS, aes(BRAFV600E_RAS, LRRC15)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = BRAFV600E_RAS),
               alpha = 0.4, 
               show.legend = FALSE,
               lwd = 1) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 2, 
              alpha = 0.9,
              show.legend = FALSE) +
  labs (x = "TCGA PTCs", y = "LRRC15 Z-Score") + 
  # scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12), limits = c(2,12)) +
  scale_fill_manual(values = c("red", "blue")) + 
  theme_classic() + 
  ylim(-1,5) +
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 30), 
    axis.title.y = element_text(face = "bold", size = 29),
    axis.text = element_text(face = "bold", size = 25)) +
  scale_x_discrete(name ="TCGA PTCs", limits = c("RAS-Like", "BRAF-Like")) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)
ggsave("outputs/22-0531_TCGA_TIMER_Exploration/22-0531_TCGA_LRRC15_Z_Score_by_BRS_Status.png",
       width = 7,
       height = 5,
       plot, dpi = 600)

# Now, plot LRRC15 by BRS Status as a violin plot
plot <- ggplot(TCGA_BRS, aes(BRAFV600E_RAS, LRRC15)) +
  geom_violin(outlier.size = -1, 
              aes(fill = BRAFV600E_RAS),
              alpha = 0.4, 
              show.legend = FALSE,
              lwd = 1,
              scale = "width") + 
  # geom_jitter(aes(),
  # position = position_jitter(width = 0.1, height = 0),
  # size = 1, 
  # alpha = 0.3,
  # show.legend = FALSE) +
  labs (x = "TCGA PTCs", y = "LRRC15 Z-Score") + 
  # scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12), limits = c(2,12)) +
  scale_fill_manual(values = c("red", "blue")) + 
  theme_classic() + 
  ylim(-1,5) +
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 30), 
    axis.title.y = element_text(face = "bold", size = 29),
    axis.text = element_text(face = "bold", size = 25)) +
  scale_x_discrete(name ="TCGA PTCs", limits = c("RAS-Like", "BRAF-Like")) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)
ggsave("outputs/22-0531_TCGA_TIMER_Exploration/22-0531_TCGA_LRRC15_Z_Score_by_BRS_Status_Violin.png",
       width = 7,
       height = 5,
       plot, dpi = 600)

# LRRC15 by BRS Status Statistics
wilcox.test(data = TCGA_BRS, LRRC15 ~ BRAFV600E_RAS) # This output: P = 3.145e-14

## Next iCAF Markers
# PDGFRA
# PDGFRA by BRS Status
# First, change PDGFRA to a numeric
TCGA_BRS$PDGFRA <- as.numeric(TCGA_BRS$PDGFRA)
# Now, plot PDGFRA by BRS Status
plot <- ggplot(TCGA_BRS, aes(BRAFV600E_RAS, PDGFRA)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = BRAFV600E_RAS),
               alpha = 0.4, 
               show.legend = FALSE,
               lwd = 1) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 2, 
              alpha = 0.9,
              show.legend = FALSE) +
  labs (x = "TCGA PTCs", y = "PDGFRA Z-Score") + 
  # scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12), limits = c(2,12)) +
  scale_fill_manual(values = c("red", "blue")) + 
  theme_classic() + 
  ylim(-1,5) +
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 30), 
    axis.title.y = element_text(face = "bold", size = 29),
    axis.text = element_text(face = "bold", size = 25)) +
  scale_x_discrete(name ="TCGA PTCs", limits = c("RAS-Like", "BRAF-Like")) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)
ggsave("outputs/22-0531_TCGA_TIMER_Exploration/22-0531_TCGA_PDGFRA_Z_Score_by_BRS_Status.png",
       width = 7,
       height = 5,
       plot, dpi = 600)

# Now, plot PDGFRA by BRS Status as a violin plot
plot <- ggplot(TCGA_BRS, aes(BRAFV600E_RAS, PDGFRA)) +
  geom_violin(outlier.size = -1, 
              aes(fill = BRAFV600E_RAS),
              alpha = 0.4, 
              show.legend = FALSE,
              lwd = 1,
              scale = "width") + 
  # geom_jitter(aes(),
  # position = position_jitter(width = 0.1, height = 0),
  # size = 1, 
  # alpha = 0.3,
  # show.legend = FALSE) +
  labs (x = "TCGA PTCs", y = "PDGFRA Z-Score") + 
  # scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12), limits = c(2,12)) +
  scale_fill_manual(values = c("red", "blue")) + 
  theme_classic() + 
  ylim(-1,5) +
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 30), 
    axis.title.y = element_text(face = "bold", size = 29),
    axis.text = element_text(face = "bold", size = 25)) +
  scale_x_discrete(name ="TCGA PTCs", limits = c("RAS-Like", "BRAF-Like")) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)
ggsave("outputs/22-0531_TCGA_TIMER_Exploration/22-0531_TCGA_PDGFRA_Z_Score_by_BRS_Status_Violin.png",
       width = 7,
       height = 5,
       plot, dpi = 600)
# PDGFRA by BRS Status Statistics
wilcox.test(data = TCGA_BRS, PDGFRA ~ BRAFV600E_RAS) # This output: P < 2.2e-16

# CXCL12
# CXCL12
# CXCL12 by BRS Status
# First, change CXCL12 to a numeric
TCGA_BRS$CXCL12 <- as.numeric(TCGA_BRS$CXCL12)
# Now, plot CXCL12 by BRS Status
plot <- ggplot(TCGA_BRS, aes(BRAFV600E_RAS, CXCL12)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = BRAFV600E_RAS),
               alpha = 0.4, 
               show.legend = FALSE,
               lwd = 1) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 2, 
              alpha = 0.9,
              show.legend = FALSE) +
  labs (x = "TCGA PTCs", y = "CXCL12 Z-Score") + 
  # scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12), limits = c(2,12)) +
  scale_fill_manual(values = c("red", "blue")) + 
  theme_classic() + 
  ylim(-2,5) +
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 30), 
    axis.title.y = element_text(face = "bold", size = 29),
    axis.text = element_text(face = "bold", size = 25)) +
  scale_x_discrete(name ="TCGA PTCs", limits = c("RAS-Like", "BRAF-Like")) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)
ggsave("outputs/22-0531_TCGA_TIMER_Exploration/22-0531_TCGA_CXCL12_Z_Score_by_BRS_Status.png",
       width = 7,
       height = 5,
       plot, dpi = 600)

# Now, plot CXCL12 by BRS Status as a violin plot
plot <- ggplot(TCGA_BRS, aes(BRAFV600E_RAS, CXCL12)) +
  geom_violin(outlier.size = -1, 
              aes(fill = BRAFV600E_RAS),
              alpha = 0.4, 
              show.legend = FALSE,
              lwd = 1,
              scale = "width") + 
  # geom_jitter(aes(),
  # position = position_jitter(width = 0.1, height = 0),
  # size = 1, 
  # alpha = 0.3,
  # show.legend = FALSE) +
  labs (x = "TCGA PTCs", y = "CXCL12 Z-Score") + 
  # scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12), limits = c(2,12)) +
  scale_fill_manual(values = c("red", "blue")) + 
  theme_classic() + 
  ylim(-2,5) +
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 30), 
    axis.title.y = element_text(face = "bold", size = 29),
    axis.text = element_text(face = "bold", size = 25)) +
  scale_x_discrete(name ="TCGA PTCs", limits = c("RAS-Like", "BRAF-Like")) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)
ggsave("outputs/22-0531_TCGA_TIMER_Exploration/22-0531_TCGA_CXCL12_Z_Score_by_BRS_Status_Violin.png",
       width = 7,
       height = 5,
       plot, dpi = 600)
# CXCL12 by BRS Status Statistics
wilcox.test(data = TCGA_BRS, CXCL12 ~ BRAFV600E_RAS) # This output: P = 8.849e-5

### Add in Histology
# Subset the data to just the following groups: 
# Follicular 
# Diffuse\nSclerosing
# Tall Cell
# Classical
TCGA_Histology <- TCGA_All_Data %>% subset(HISTOLOGICAL_TYPE_COMPLEX == "Classical" |
                                           HISTOLOGICAL_TYPE_COMPLEX == "Diffuse\nSclerosing" |
                                           HISTOLOGICAL_TYPE_COMPLEX == "Follicular" | 
                                           HISTOLOGICAL_TYPE_COMPLEX == "Tall\nCell")
# End of line 597 from 22-0426

# This is also from the same R file but lines 616-623
### Make a new plot combining diffuse sclerosing and tall cell
TCGA_Histology$HISTOLOGICAL_TYPE_SIMPLIFIED <- TCGA_Histology$HISTOLOGICAL_TYPE_COMPLEX
for(i in 1:nrow(Histology_Plot_Cohort)){
  if(TCGA_Histology$HISTOLOGICAL_TYPE_COMPLEX[i] == "Diffuse\nSclerosing" |
     TCGA_Histology$HISTOLOGICAL_TYPE_COMPLEX[i] == "Tall\nCell"){
     TCGA_Histology$HISTOLOGICAL_TYPE_SIMPLIFIED[i] <- "Diffuse\nSclerosing/\nTall Cell"
  }
}

# For this, load in RNA data again 
## old code
# Read in TCGA RNA Data (Z-Score format) and clean the data to the desired format
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

TCGA_Histology <- merge(TCGA_Histology, TCGA_RNA_Tibble)

# Change FAP in TCGA_Histology to a numeric
TCGA_Histology$FAP <- as.numeric(TCGA_Histology$FAP)
# Now, plot FAP by HISTOLOGY
plot <- ggplot(TCGA_Histology, aes(HISTOLOGICAL_TYPE_SIMPLIFIED, FAP)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = HISTOLOGICAL_TYPE_SIMPLIFIED),
               alpha = 0.4, 
               show.legend = FALSE,
               lwd = 2) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 2, 
              alpha = 0.9,
              show.legend = FALSE) +
  labs (x = "TCGA PTCs", y = "FAP Z-Score") + 
  # scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12), limits = c(2,12)) +
  scale_fill_manual(values = c("red", "purple", "blue")) + 
  theme_classic() + 
  ggtitle("FAP Expression by Histology") +
  ylim(-1,5) +
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 30), 
    axis.title.x = element_text(face = "bold", size = 30), 
    axis.title.y = element_text(face = "bold", size = 29),
    axis.text = element_text(face = "bold", size = 25)) +
  scale_x_discrete(name ="TCGA PTCs", limits = c("Follicular", "Classical", "Diffuse\nSclerosing/\nTall Cell")) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)
ggsave("outputs/22-0531_TCGA_TIMER_Exploration/22-0531_TCGA_FAP_Z_Score_by_Histology.png",
       width = 7,
       height = 5,
       plot, dpi = 600)
# Now, plot FAP by HISTOLOGY as a violin plot
plot <- ggplot(TCGA_Histology, aes(HISTOLOGICAL_TYPE_SIMPLIFIED, FAP)) +
  geom_violin(outlier.size = -1, 
               aes(fill = HISTOLOGICAL_TYPE_SIMPLIFIED),
               alpha = 0.4, 
               show.legend = FALSE,
               lwd = 1, 
               scale = "width") + 
  #geom_jitter(aes(),
              #position = position_jitter(width = 0.1, height = 0),
              #size = 2, 
              #alpha = 0.9,
              #show.legend = FALSE) +
  labs (x = "TCGA PTCs", y = "FAP Z-Score") + 
  # scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12), limits = c(2,12)) +
  scale_fill_manual(values = c("red", "purple", "blue")) + 
  theme_classic() + 
  ggtitle("FAP Expression by Histology") +
  ylim(-1,5) +
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 30), 
    axis.title.x = element_text(face = "bold", size = 30), 
    axis.title.y = element_text(face = "bold", size = 29),
    axis.text = element_text(face = "bold", size = 25)) +
  scale_x_discrete(name ="TCGA PTCs", limits = c("Follicular", "Classical", "Diffuse\nSclerosing/\nTall Cell")) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)
ggsave("outputs/22-0531_TCGA_TIMER_Exploration/22-0531_TCGA_FAP_Z_Score_by_Histology_violin.png",
       width = 7,
       height = 5,
       plot, dpi = 600)
# Stats for FAP by histology
kruskal.test(FAP ~ HISTOLOGICAL_TYPE_SIMPLIFIED, data = TCGA_Histology) # 
pairwise.wilcox.test(TCGA_Histology$FAP, TCGA_Histology$HISTOLOGICAL_TYPE_SIMPLIFIED, # 
                     p.adjust.method = "BH")

# Change ACTA2 in TCGA_Histology to a numeric
TCGA_Histology$ACTA2 <- as.numeric(TCGA_Histology$ACTA2)
# Now, plot ACTA2 by HISTOLOGY
plot <- ggplot(TCGA_Histology, aes(HISTOLOGICAL_TYPE_SIMPLIFIED, ACTA2)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = HISTOLOGICAL_TYPE_SIMPLIFIED),
               alpha = 0.4, 
               show.legend = FALSE,
               lwd = 2) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 2, 
              alpha = 0.9,
              show.legend = FALSE) +
  labs (x = "TCGA PTCs", y = "ACTA2 Z-Score") + 
  # scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12), limits = c(2,12)) +
  scale_fill_manual(values = c("red", "purple", "blue")) + 
  theme_classic() + 
  ggtitle("ACTA2 Expression by Histology") +
  ylim(-2,5) +
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 29), 
    axis.title.x = element_text(face = "bold", size = 30), 
    axis.title.y = element_text(face = "bold", size = 29),
    axis.text = element_text(face = "bold", size = 25)) +
  scale_x_discrete(name ="TCGA PTCs", limits = c("Follicular", "Classical", "Diffuse\nSclerosing/\nTall Cell")) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)
ggsave("outputs/22-0531_TCGA_TIMER_Exploration/22-0531_TCGA_ACTA2_Z_Score_by_Histology.png",
       width = 7,
       height = 5,
       plot, dpi = 600)

# Now, plot ACTA2 by HISTOLOGY as a violin plot
plot <- ggplot(TCGA_Histology, aes(HISTOLOGICAL_TYPE_SIMPLIFIED, ACTA2)) +
  geom_violin(outlier.size = -1, 
              aes(fill = HISTOLOGICAL_TYPE_SIMPLIFIED),
              alpha = 0.4, 
              show.legend = FALSE,
              lwd = 1,
              scale = "width") + 
  #geom_jitter(aes(),
  #position = position_jitter(width = 0.1, height = 0),
  #size = 2, 
  #alpha = 0.9,
  #show.legend = FALSE) +
  labs (x = "TCGA PTCs", y = "ACTA2 Z-Score") + 
  # scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12), limits = c(2,12)) +
  scale_fill_manual(values = c("red", "purple", "blue")) + 
  theme_classic() + 
  ggtitle("ACTA2 Expression by Histology") +
  ylim(-2,5) +
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 29), 
    axis.title.x = element_text(face = "bold", size = 30), 
    axis.title.y = element_text(face = "bold", size = 29),
    axis.text = element_text(face = "bold", size = 25)) +
  scale_x_discrete(name ="TCGA PTCs", limits = c("Follicular", "Classical", "Diffuse\nSclerosing/\nTall Cell")) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)
ggsave("outputs/22-0531_TCGA_TIMER_Exploration/22-0531_TCGA_ACTA2_Z_Score_by_Histology_violin.png",
       width = 7,
       height = 5,
       plot, dpi = 600)

# Stats for ACTA2 by histology
kruskal.test(ACTA2 ~ HISTOLOGICAL_TYPE_SIMPLIFIED, data = TCGA_Histology) # 
pairwise.wilcox.test(TCGA_Histology$ACTA2, TCGA_Histology$HISTOLOGICAL_TYPE_SIMPLIFIED, # 
                     p.adjust.method = "BH")

#### myCAF and iCAF markers by histology ####
## myCAF: POSTN, WNT2, COL1A1, LRRC15
## iCAF: PDGFRA, CXCL12
# Change ACTA2 in TCGA_Histology to a numeric
TCGA_Histology$POSTN <- as.numeric(TCGA_Histology$POSTN)
# Now, plot POSTN by HISTOLOGY
plot <- ggplot(TCGA_Histology, aes(HISTOLOGICAL_TYPE_SIMPLIFIED, POSTN)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = HISTOLOGICAL_TYPE_SIMPLIFIED),
               alpha = 0.4, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 2, 
              alpha = 0.9,
              show.legend = FALSE) +
  labs (x = "TCGA PTCs", y = "POSTN Z-Score") + 
  # scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12), limits = c(2,12)) +
  scale_fill_manual(values = c("red", "purple", "blue")) + 
  theme_classic() + 
  ggtitle("POSTN Expression by Histology") +
  ylim(-1,5) +
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 29), 
    axis.title.x = element_text(face = "bold", size = 30), 
    axis.title.y = element_text(face = "bold", size = 29),
    axis.text = element_text(face = "bold", size = 25)) +
  scale_x_discrete(name ="TCGA PTCs", limits = c("Follicular", "Classical", "Diffuse\nSclerosing/\nTall Cell")) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)
ggsave("outputs/22-0531_TCGA_TIMER_Exploration/22-0531_TCGA_POSTN_Z_Score_by_Histology.png",
       width = 7,
       height = 5,
       plot, dpi = 600)

# Now, plot POSTN by HISTOLOGY as a violin plot
plot <- ggplot(TCGA_Histology, aes(HISTOLOGICAL_TYPE_SIMPLIFIED, POSTN)) +
  geom_violin(outlier.size = -1, 
              aes(fill = HISTOLOGICAL_TYPE_SIMPLIFIED),
              alpha = 0.4, 
              show.legend = FALSE,
              lwd = 1, 
              scale = "width") + 
  #geom_jitter(aes(),
  #position = position_jitter(width = 0.1, height = 0),
  #size = 2, 
  #alpha = 0.9,
  #show.legend = FALSE) +
  labs (x = "TCGA PTCs", y = "POSTN Z-Score") + 
  # scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12), limits = c(2,12)) +
  scale_fill_manual(values = c("red", "purple", "blue")) + 
  theme_classic() + 
  ggtitle("POSTN Expression by Histology") +
  ylim(-1,5) +
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 29), 
    axis.title.x = element_text(face = "bold", size = 30), 
    axis.title.y = element_text(face = "bold", size = 29),
    axis.text = element_text(face = "bold", size = 25)) +
  scale_x_discrete(name ="TCGA PTCs", limits = c("Follicular", "Classical", "Diffuse\nSclerosing/\nTall Cell")) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)
ggsave("outputs/22-0531_TCGA_TIMER_Exploration/22-0531_TCGA_POSTN_Z_Score_by_Histology_violin.png",
       width = 7,
       height = 5,
       plot, dpi = 600)

# Stats for POSTN by histology
kruskal.test(POSTN ~ HISTOLOGICAL_TYPE_SIMPLIFIED, data = TCGA_Histology) # p-value = 5.259e-13
pairwise.wilcox.test(TCGA_Histology$POSTN, TCGA_Histology$HISTOLOGICAL_TYPE_SIMPLIFIED, # F vs C: 2.3e-6; D vs C: 1.8e-7
                     p.adjust.method = "BH")

# Change WNT2 in TCGA_Histology to a numeric
TCGA_Histology$WNT2 <- as.numeric(TCGA_Histology$WNT2)
# Now, plot WNT2 by HISTOLOGY
plot <- ggplot(TCGA_Histology, aes(HISTOLOGICAL_TYPE_SIMPLIFIED, WNT2)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = HISTOLOGICAL_TYPE_SIMPLIFIED),
               alpha = 0.4, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 2, 
              alpha = 0.9,
              show.legend = FALSE) +
  labs (x = "TCGA PTCs", y = "WNT2 Z-Score") + 
  # scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12), limits = c(2,12)) +
  scale_fill_manual(values = c("red", "purple", "blue")) + 
  theme_classic() + 
  ggtitle("WNT2 Expression by Histology") +
  ylim(-1,5) +
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 29), 
    axis.title.x = element_text(face = "bold", size = 30), 
    axis.title.y = element_text(face = "bold", size = 29),
    axis.text = element_text(face = "bold", size = 25)) +
  scale_x_discrete(name ="TCGA PTCs", limits = c("Follicular", "Classical", "Diffuse\nSclerosing/\nTall Cell")) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)
ggsave("outputs/22-0531_TCGA_TIMER_Exploration/22-0531_TCGA_WNT2_Z_Score_by_Histology.png",
       width = 7,
       height = 5,
       plot, dpi = 600)

# Now, plot WNT2 by HISTOLOGY as a violin plot
plot <- ggplot(TCGA_Histology, aes(HISTOLOGICAL_TYPE_SIMPLIFIED, WNT2)) +
  geom_violin(outlier.size = -1, 
              aes(fill = HISTOLOGICAL_TYPE_SIMPLIFIED),
              alpha = 0.4, 
              show.legend = FALSE,
              lwd = 1,
              scale = "width") + 
  #geom_jitter(aes(),
  #position = position_jitter(width = 0.1, height = 0),
  #size = 2, 
  #alpha = 0.9,
  #show.legend = FALSE) +
  labs (x = "TCGA PTCs", y = "WNT2 Z-Score") + 
  # scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12), limits = c(2,12)) +
  scale_fill_manual(values = c("red", "purple", "blue")) + 
  theme_classic() + 
  ggtitle("WNT2 Expression by Histology") +
  ylim(-1,5) +
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 29), 
    axis.title.x = element_text(face = "bold", size = 30), 
    axis.title.y = element_text(face = "bold", size = 29),
    axis.text = element_text(face = "bold", size = 25)) +
  scale_x_discrete(name ="TCGA PTCs", limits = c("Follicular", "Classical", "Diffuse\nSclerosing/\nTall Cell")) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)
ggsave("outputs/22-0531_TCGA_TIMER_Exploration/22-0531_TCGA_WNT2_Z_Score_by_Histology_violin.png",
       width = 7,
       height = 5,
       plot, dpi = 600)

# Stats for WNT2 by histology
kruskal.test(WNT2 ~ HISTOLOGICAL_TYPE_SIMPLIFIED, data = TCGA_Histology) # P < 2.2e-16
pairwise.wilcox.test(TCGA_Histology$WNT2, TCGA_Histology$HISTOLOGICAL_TYPE_SIMPLIFIED, # F vs C: 6.8e-11; D vs C: 1.7e-6
                     p.adjust.method = "BH")

# Change COL1A1 in TCGA_Histology to a numeric
TCGA_Histology$COL1A1 <- as.numeric(TCGA_Histology$COL1A1)
# Now, plot COL1A1 by HISTOLOGY
plot <- ggplot(TCGA_Histology, aes(HISTOLOGICAL_TYPE_SIMPLIFIED, COL1A1)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = HISTOLOGICAL_TYPE_SIMPLIFIED),
               alpha = 0.4, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 2, 
              alpha = 0.9,
              show.legend = FALSE) +
  labs (x = "TCGA PTCs", y = "COL1A1 Z-Score") + 
  # scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12), limits = c(2,12)) +
  scale_fill_manual(values = c("red", "purple", "blue")) + 
  theme_classic() + 
  ggtitle("COL1A1 Expression by Histology") +
  ylim(-1,5) +
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 29), 
    axis.title.x = element_text(face = "bold", size = 30), 
    axis.title.y = element_text(face = "bold", size = 29),
    axis.text = element_text(face = "bold", size = 25)) +
  scale_x_discrete(name ="TCGA PTCs", limits = c("Follicular", "Classical", "Diffuse\nSclerosing/\nTall Cell")) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)
ggsave("outputs/22-0531_TCGA_TIMER_Exploration/22-0531_TCGA_COL1A1_Z_Score_by_Histology.png",
       width = 7,
       height = 5,
       plot, dpi = 600)

# Now, plot COL1A1 by HISTOLOGY as a violin plot
plot <- ggplot(TCGA_Histology, aes(HISTOLOGICAL_TYPE_SIMPLIFIED, COL1A1)) +
  geom_violin(outlier.size = -1, 
              aes(fill = HISTOLOGICAL_TYPE_SIMPLIFIED),
              alpha = 0.4, 
              show.legend = FALSE,
              lwd = 1, 
              scale = "width") + 
  #geom_jitter(aes(),
  #position = position_jitter(width = 0.1, height = 0),
  #size = 2, 
  #alpha = 0.9,
  #show.legend = FALSE) +
  labs (x = "TCGA PTCs", y = "COL1A1 Z-Score") + 
  # scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12), limits = c(2,12)) +
  scale_fill_manual(values = c("red", "purple", "blue")) + 
  theme_classic() + 
  ggtitle("COL1A1 Expression by Histology") +
  ylim(-1,5) +
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 29), 
    axis.title.x = element_text(face = "bold", size = 30), 
    axis.title.y = element_text(face = "bold", size = 29),
    axis.text = element_text(face = "bold", size = 25)) +
  scale_x_discrete(name ="TCGA PTCs", limits = c("Follicular", "Classical", "Diffuse\nSclerosing/\nTall Cell")) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)
ggsave("outputs/22-0531_TCGA_TIMER_Exploration/22-0531_TCGA_COL1A1_Z_Score_by_Histology_violin.png",
       width = 7,
       height = 5,
       plot, dpi = 600)

# Stats for COL1A1 by histology
kruskal.test(COL1A1 ~ HISTOLOGICAL_TYPE_SIMPLIFIED, data = TCGA_Histology) # P < 2.2e-16
pairwise.wilcox.test(TCGA_Histology$COL1A1, TCGA_Histology$HISTOLOGICAL_TYPE_SIMPLIFIED, # F Vs C: < 2e-16; C vs D: 7.7e-7
                     p.adjust.method = "BH")

# Change LRRC15 in TCGA_Histology to a numeric
TCGA_Histology$LRRC15 <- as.numeric(TCGA_Histology$LRRC15)
# Now, plot LRRC15 by HISTOLOGY
plot <- ggplot(TCGA_Histology, aes(HISTOLOGICAL_TYPE_SIMPLIFIED, LRRC15)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = HISTOLOGICAL_TYPE_SIMPLIFIED),
               alpha = 0.4, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 2, 
              alpha = 0.9,
              show.legend = FALSE) +
  labs (x = "TCGA PTCs", y = "LRRC15 Z-Score") + 
  # scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12), limits = c(2,12)) +
  scale_fill_manual(values = c("red", "purple", "blue")) + 
  theme_classic() + 
  ggtitle("LRRC15 Expression by Histology") +
  ylim(-1,5) +
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 29), 
    axis.title.x = element_text(face = "bold", size = 30), 
    axis.title.y = element_text(face = "bold", size = 29),
    axis.text = element_text(face = "bold", size = 25)) +
  scale_x_discrete(name ="TCGA PTCs", limits = c("Follicular", "Classical", "Diffuse\nSclerosing/\nTall Cell")) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)
ggsave("outputs/22-0531_TCGA_TIMER_Exploration/22-0531_TCGA_LRRC15_Z_Score_by_Histology.png",
       width = 7,
       height = 5,
       plot, dpi = 600)

# Now, plot LRRC15 by HISTOLOGY as a violin plot
plot <- ggplot(TCGA_Histology, aes(HISTOLOGICAL_TYPE_SIMPLIFIED, LRRC15)) +
  geom_violin(outlier.size = -1, 
              aes(fill = HISTOLOGICAL_TYPE_SIMPLIFIED),
              alpha = 0.4, 
              show.legend = FALSE,
              lwd = 1, 
              scale = "width") + 
  #geom_jitter(aes(),
  #position = position_jitter(width = 0.1, height = 0),
  #size = 2, 
  #alpha = 0.9,
  #show.legend = FALSE) +
  labs (x = "TCGA PTCs", y = "LRRC15 Z-Score") + 
  # scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12), limits = c(2,12)) +
  scale_fill_manual(values = c("red", "purple", "blue")) + 
  theme_classic() + 
  ggtitle("LRRC15 Expression by Histology") +
  ylim(-1,5) +
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 29), 
    axis.title.x = element_text(face = "bold", size = 30), 
    axis.title.y = element_text(face = "bold", size = 29),
    axis.text = element_text(face = "bold", size = 25)) +
  scale_x_discrete(name ="TCGA PTCs", limits = c("Follicular", "Classical", "Diffuse\nSclerosing/\nTall Cell")) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)
ggsave("outputs/22-0531_TCGA_TIMER_Exploration/22-0531_TCGA_LRRC15_Z_Score_by_Histology_violin.png",
       width = 7,
       height = 5,
       plot, dpi = 600)

# Stats for LRRC15 by histology
kruskal.test(LRRC15 ~ HISTOLOGICAL_TYPE_SIMPLIFIED, data = TCGA_Histology) # P = 1.515e-13
pairwise.wilcox.test(TCGA_Histology$LRRC15, TCGA_Histology$HISTOLOGICAL_TYPE_SIMPLIFIED, # F vs C: 1.3e-7; D vs C: 4.5e-7
                     p.adjust.method = "BH")

## iCAF Markers: PDGFRA, CXCL12
# Change PDGFRA in TCGA_Histology to a numeric
TCGA_Histology$PDGFRA <- as.numeric(TCGA_Histology$PDGFRA)
# Now, plot PDGFRA by HISTOLOGY
plot <- ggplot(TCGA_Histology, aes(HISTOLOGICAL_TYPE_SIMPLIFIED, PDGFRA)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = HISTOLOGICAL_TYPE_SIMPLIFIED),
               alpha = 0.4, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 2, 
              alpha = 0.9,
              show.legend = FALSE) +
  labs (x = "TCGA PTCs", y = "PDGFRA Z-Score") + 
  # scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12), limits = c(2,12)) +
  scale_fill_manual(values = c("red", "purple", "blue")) + 
  theme_classic() + 
  ggtitle("PDGFRA Expression by Histology") +
  ylim(-1,5) +
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 28), 
    axis.title.x = element_text(face = "bold", size = 30), 
    axis.title.y = element_text(face = "bold", size = 29),
    axis.text = element_text(face = "bold", size = 25)) +
  scale_x_discrete(name ="TCGA PTCs", limits = c("Follicular", "Classical", "Diffuse\nSclerosing/\nTall Cell")) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)
ggsave("outputs/22-0531_TCGA_TIMER_Exploration/22-0531_TCGA_PDGFRA_Z_Score_by_Histology.png",
       width = 7,
       height = 5,
       plot, dpi = 600)

# Now, plot PDGFRA by HISTOLOGY as a violin plot
plot <- ggplot(TCGA_Histology, aes(HISTOLOGICAL_TYPE_SIMPLIFIED, PDGFRA)) +
  geom_violin(outlier.size = -1, 
              aes(fill = HISTOLOGICAL_TYPE_SIMPLIFIED),
              alpha = 0.4, 
              show.legend = FALSE,
              lwd = 1, 
              scale = "width") + 
  #geom_jitter(aes(),
  #position = position_jitter(width = 0.1, height = 0),
  #size = 2, 
  #alpha = 0.9,
  #show.legend = FALSE) +
  labs (x = "TCGA PTCs", y = "PDGFRA Z-Score") + 
  # scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12), limits = c(2,12)) +
  scale_fill_manual(values = c("red", "purple", "blue")) + 
  theme_classic() + 
  ggtitle("PDGFRA Expression by Histology") +
  ylim(-1,5) +
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 29), 
    axis.title.x = element_text(face = "bold", size = 30), 
    axis.title.y = element_text(face = "bold", size = 29),
    axis.text = element_text(face = "bold", size = 25)) +
  scale_x_discrete(name ="TCGA PTCs", limits = c("Follicular", "Classical", "Diffuse\nSclerosing/\nTall Cell")) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)
ggsave("outputs/22-0531_TCGA_TIMER_Exploration/22-0531_TCGA_PDGFRA_Z_Score_by_Histology_violin.png",
       width = 7,
       height = 5,
       plot, dpi = 600)

# Stats for PDGFRA by histology
kruskal.test(PDGFRA ~ HISTOLOGICAL_TYPE_SIMPLIFIED, data = TCGA_Histology) # P = 1.954e-12
pairwise.wilcox.test(TCGA_Histology$PDGFRA, TCGA_Histology$HISTOLOGICAL_TYPE_SIMPLIFIED, # F vs C: 2.0e-10; D vs C: 0.017
                     p.adjust.method = "BH")

# Change CXCL12 in TCGA_Histology to a numeric
TCGA_Histology$CXCL12 <- as.numeric(TCGA_Histology$CXCL12)
# Now, plot CXCL12 by HISTOLOGY
plot <- ggplot(TCGA_Histology, aes(HISTOLOGICAL_TYPE_SIMPLIFIED, CXCL12)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = HISTOLOGICAL_TYPE_SIMPLIFIED),
               alpha = 0.4, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 2, 
              alpha = 0.9,
              show.legend = FALSE) +
  labs (x = "TCGA PTCs", y = "CXCL12 Z-Score") + 
  # scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12), limits = c(2,12)) +
  scale_fill_manual(values = c("red", "purple", "blue")) + 
  theme_classic() + 
  ggtitle("CXCL12 Expression by Histology") +
  ylim(-2,5) +
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 28), 
    axis.title.x = element_text(face = "bold", size = 30), 
    axis.title.y = element_text(face = "bold", size = 29),
    axis.text = element_text(face = "bold", size = 25)) +
  scale_x_discrete(name ="TCGA PTCs", limits = c("Follicular", "Classical", "Diffuse\nSclerosing/\nTall Cell")) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)
ggsave("outputs/22-0531_TCGA_TIMER_Exploration/22-0531_TCGA_CXCL12_Z_Score_by_Histology.png",
       width = 7,
       height = 5,
       plot, dpi = 600)

# Now, plot CXCL12 by HISTOLOGY as a violin plot
plot <- ggplot(TCGA_Histology, aes(HISTOLOGICAL_TYPE_SIMPLIFIED, CXCL12)) +
  geom_violin(outlier.size = -1, 
              aes(fill = HISTOLOGICAL_TYPE_SIMPLIFIED),
              alpha = 0.4, 
              show.legend = FALSE,
              lwd = 1, 
              scale = "width") + 
  #geom_jitter(aes(),
  #position = position_jitter(width = 0.1, height = 0),
  #size = 2, 
  #alpha = 0.9,
  #show.legend = FALSE) +
  labs (x = "TCGA PTCs", y = "CXCL12 Z-Score") + 
  # scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12), limits = c(2,12)) +
  scale_fill_manual(values = c("red", "purple", "blue")) + 
  theme_classic() + 
  ggtitle("CXCL12 Expression by Histology") +
  ylim(-2,5) +
  theme(
    panel.border = element_rect(colour = "black", size = 4, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 28), 
    axis.title.x = element_text(face = "bold", size = 30), 
    axis.title.y = element_text(face = "bold", size = 29),
    axis.text = element_text(face = "bold", size = 25)) +
  scale_x_discrete(name ="TCGA PTCs", limits = c("Follicular", "Classical", "Diffuse\nSclerosing/\nTall Cell")) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)
ggsave("outputs/22-0531_TCGA_TIMER_Exploration/22-0531_TCGA_CXCL12_Z_Score_by_Histology_violin.png",
       width = 7,
       height = 5,
       plot, dpi = 600)

# Stats for CXCL12 by histology
kruskal.test(CXCL12 ~ HISTOLOGICAL_TYPE_SIMPLIFIED, data = TCGA_Histology) # P = 0.009115
pairwise.wilcox.test(TCGA_Histology$CXCL12, TCGA_Histology$HISTOLOGICAL_TYPE_SIMPLIFIED, # F vs C: 0.012; D vs C: 0.603
                     p.adjust.method = "BH")

##### Now adding survival curves for High-CAF vs Low-CAF #####
library(survival)
library(survminer)
library(cowplot) # Using cowplot due to ggsave issues with survival curves

# DFS_STATUS needs to be numeric, so separate column
TCGA_All_Data <-separate(data = TCGA_All_Data, col = DFS_STATUS, into = c("DFS_Status_Numeric", "DFS_Status_Char"), sep = ":")
TCGA_All_Data$DFS_Status_Numeric <- as.numeric(TCGA_All_Data$DFS_Status_Numeric) # Change to numeric

# OS_STATUS needs to be numeric, so separate column
TCGA_All_Data <-separate(data = TCGA_All_Data, col = OS_STATUS, into = c("OS_Status_Numeric", "OS_Status_Char"), sep = ":")
TCGA_All_Data$OS_Status_Numeric <- as.numeric(TCGA_All_Data$OS_Status_Numeric) # Change to numeric

# Looking at DFS and OS
# Make a new tibble for each of the above restricted to samples with DFS month data
DFS_CAFs <- TCGA_All_Data %>% subset(!(is.na(DFS_MONTHS)))

# Make a new tibble for each of the above restricted to samples with OS month data
OS_CAFs <- TCGA_All_Data %>% subset(!(is.na(OS_MONTHS)))

#### EPIC CAF High vs EPIC CAF Low Disease Free Survival Plot ####
# Summary statistics on the DFS Samples CAF high vs CAF Low
summary(DFS_CAFs$CAF_EPIC)
DFS_CAF_Quantile <- quantile(DFS_CAFs$CAF_EPIC)
DFS_CAF_Quantile # Print DFS PoorOutcome Quantiles
# Sort samples into greater than or less than 50th percentile based on PoorOutcome Score
# Rather than the 50th percentile, I will be using > 0 and < 0
DFS_CAFs <- mutate(DFS_CAFs, Cat=ifelse(DFS_CAFs$CAF_EPIC > DFS_CAF_Quantile[3], "> 50th", "< 50th"))

#Survival curve calculation
fit<-survfit(Surv(DFS_MONTHS, DFS_Status_Numeric) ~ Cat, data=DFS_CAFs)
print(fit)

# Most basic plot
ggsurvplot(fit, data=DFS_CAFs)

# More involved, customized plot
DFS_Plot <- ggsurvplot(fit,
                       data=DFS_CAFs,
                       pval = TRUE, pval.method = TRUE, conf.int = FALSE,
                       risk.table = TRUE, 
                       risk.table.col = "strata", 
                       legend = c(.80,.2),
                       pval.coord = c(.15,.15),
                       pval.size = 8,
                       break.time.by = 24,
                       pval.method.coord = c(.15,.06),
                       font.x = c(23, "bold", "black"),
                       font.y = c(23, "bold", "black"),
                       font.title = c(25, "bold", "black"),
                       # linetype = "strata", # Change line type by groups
                       # surv.median.line = "hv", # Specify median survival
                       ggtheme = theme_classic(),
                       font.tickslab=c(15, "bold", "black"),
                       font.legend = c(20, "bold", "black"),
                       risk.table.fontsize = 8,
                       risk.table.y.text = FALSE,
                       risk.table.legend = FALSE,
                       legend.title = "CAF Score",
                       title = "CAF EPIC DFS",
                       legend.labs = c("Low", "High"),
                       tables.height = .33,
                       palette = c("#2E9FDF", "red", "#86AA00"))+
  labs(y = "Disease Free Survival", x = "Months")  
DFS_Plot$plot <- DFS_Plot$plot + theme(panel.border = element_rect(colour = "black", size = 4, fill = NA))
DFS_Plot$table <- DFS_Plot$table + 
  theme(panel.border = element_rect(colour = "black", size = 4, fill = NA),
        plot.title = element_text(size = 23, color = "black", face = "bold"),
        axis.title = element_text(size = 23, color = "black", face = "bold"),
        axis.text = element_text(size = 15, color = "black", face = "bold"),
        legend.position = "none")+
  labs(y = NULL) 
DFS_Plot

# Saving the plots with ggsave
p1 = DFS_Plot$plot
p2 = DFS_Plot$table
plotp = cowplot::plot_grid(p1,p2,align = "v",axis = "b", ncol =1, rel_heights = c(3,1))
ggsave("outputs/22-0531_TCGA_TIMER_Exploration/22-0531_TCGA_CAF_DFS.png",
       plot = plotp, 
       width = 7, height = 7,
       dpi = 600)

#### EPIC CAF High vs EPIC CAF Low Overall Survival Plot ####
# Summary statistics on the OS Samples CAF high vs CAF Low
summary(OS_CAFs$CAF_EPIC)
OS_CAF_Quantile <- quantile(OS_CAFs$CAF_EPIC)
OS_CAF_Quantile # Print OS PoorOutcome Quantiles
# Sort samples into greater than or less than 50th percentile based on PoorOutcome Score
# Rather than the 50th percentile, I will be using > 0 and < 0
OS_CAFs <- mutate(OS_CAFs, Cat=ifelse(OS_CAFs$CAF_EPIC > OS_CAF_Quantile[3], "> 50th", "< 50th"))

#Survival curve calculation
fit<-survfit(Surv(OS_MONTHS, OS_Status_Numeric) ~ Cat, data=OS_CAFs)
print(fit)

# Most basic plot
ggsurvplot(fit, data=OS_CAFs)

# More involved, customized plot
OS_Plot <- ggsurvplot(fit,
                       data=OS_CAFs,
                       pval = TRUE, pval.method = TRUE, conf.int = FALSE,
                       risk.table = TRUE, 
                       risk.table.col = "strata", 
                       legend = c(.80,.50),
                       pval.coord = c(.15,.15),
                       pval.size = 8,
                       break.time.by = 24,
                       pval.method.coord = c(.15,.06),
                       font.x = c(23, "bold", "black"),
                       font.y = c(23, "bold", "black"),
                       font.title = c(25, "bold", "black"),
                       # linetype = "strata", # Change line type by groups
                       # surv.median.line = "hv", # Specify median survival
                       ggtheme = theme_classic(),
                       font.tickslab=c(15, "bold", "black"),
                       font.legend = c(20, "bold", "black"),
                       risk.table.fontsize = 8,
                       risk.table.y.text = FALSE,
                       risk.table.legend = FALSE,
                       legend.title = "CAF EPIC Level",
                       title = "CAF EPIC OS",
                       legend.labs = c("Low", "High"),
                       tables.height = .33,
                       palette = c("#2E9FDF", "red", "#86AA00"))+
  labs(y = "Disease Free Survival", x = "Months")  
OS_Plot$plot <- OS_Plot$plot + theme(panel.border = element_rect(colour = "black", size = 4, fill = NA))
OS_Plot$table <- OS_Plot$table + 
  theme(panel.border = element_rect(colour = "black", size = 4, fill = NA),
        plot.title = element_text(size = 23, color = "black", face = "bold"),
        axis.title = element_text(size = 23, color = "black", face = "bold"),
        axis.text = element_text(size = 15, color = "black", face = "bold"),
        legend.position = "none")+
  labs(y = NULL) 
OS_Plot

# Saving the plots with ggsave
p1 = OS_Plot$plot
p2 = OS_Plot$table
plotp = cowplot::plot_grid(p1,p2,align = "v",axis = "b", ncol =1, rel_heights = c(3,1))
ggsave("outputs/22-0531_TCGA_TIMER_Exploration/22-0531_TCGA_CAF_OS.png",
       plot = plotp, 
       width = 7, height = 7,
       dpi = 600)

##### Repeat the above survival plots with MCPCOUNTER instead of EPIC #####
### MCPCOUNTER CAF High vs MCPCOUNTER CAF Low Disease Free Survival Plot
# Summary statistics on the DFS Samples CAF high vs CAF Low
summary(DFS_CAFs$CAF_MCPCOUNTER)
DFS_CAF_Quantile <- quantile(DFS_CAFs$CAF_MCPCOUNTER)
DFS_CAF_Quantile # Print DFS PoorOutcome Quantiles
# Sort samples into greater than or less than 50th percentile based on PoorOutcome Score
# Rather than the 50th percentile, I will be using > 0 and < 0
DFS_CAFs <- mutate(DFS_CAFs, Cat=ifelse(DFS_CAFs$CAF_MCPCOUNTER > DFS_CAF_Quantile[3], "> 50th", "< 50th"))

#Survival curve calculation
fit<-survfit(Surv(DFS_MONTHS, DFS_Status_Numeric) ~ Cat, data=DFS_CAFs)
print(fit)

# Most basic plot
ggsurvplot(fit, data=DFS_CAFs)

# More involved, customized plot
DFS_Plot <- ggsurvplot(fit,
                       data=DFS_CAFs,
                       pval = TRUE, pval.method = TRUE, conf.int = FALSE,
                       risk.table = TRUE, 
                       risk.table.col = "strata", 
                       legend = c(.8,.20),
                       pval.coord = c(.15,.15),
                       pval.size = 8,
                       break.time.by = 24,
                       pval.method.coord = c(.15,.06),
                       font.x = c(23, "bold", "black"),
                       font.y = c(23, "bold", "black"),
                       font.title = c(25, "bold", "black"),
                       # linetype = "strata", # Change line type by groups
                       # surv.median.line = "hv", # Specify median survival
                       ggtheme = theme_classic(),
                       font.tickslab=c(15, "bold", "black"),
                       font.legend = c(20, "bold", "black"),
                       risk.table.fontsize = 8,
                       risk.table.y.text = FALSE,
                       risk.table.legend = FALSE,
                       legend.title = "CAF Score",
                       title = "CAF MCPCOUNTER DFS",
                       legend.labs = c("Low", "High"),
                       tables.height = .33,
                       palette = c("#2E9FDF", "red", "#86AA00"))+
  labs(y = "Disease Free Survival", x = "Months")  
DFS_Plot$plot <- DFS_Plot$plot + theme(panel.border = element_rect(colour = "black", size = 4, fill = NA))
DFS_Plot$table <- DFS_Plot$table + 
  theme(panel.border = element_rect(colour = "black", size = 4, fill = NA),
        plot.title = element_text(size = 23, color = "black", face = "bold"),
        axis.title = element_text(size = 23, color = "black", face = "bold"),
        axis.text = element_text(size = 15, color = "black", face = "bold"),
        legend.position = "none")+
  labs(y = NULL) 
DFS_Plot

# Saving the plots with ggsave
p1 = DFS_Plot$plot
p2 = DFS_Plot$table
plotp = cowplot::plot_grid(p1,p2,align = "v",axis = "b", ncol =1, rel_heights = c(3,1))
ggsave("outputs/22-0531_TCGA_TIMER_Exploration/22-0531_TCGA_CAF_DFS_MCPCOUNTER.png",
       plot = plotp, 
       width = 7, height = 7,
       dpi = 600)

### MCPCOUNTER CAF High vs MCPCOUNTER CAF Low Overall Survival Plot
# Summary statistics on the OS Samples CAF high vs CAF Low
summary(OS_CAFs$CAF_MCPCOUNTER)
OS_CAF_Quantile <- quantile(OS_CAFs$CAF_MCPCOUNTER)
OS_CAF_Quantile # Print OS PoorOutcome Quantiles
# Sort samples into greater than or less than 50th percentile based on PoorOutcome Score
# Rather than the 50th percentile, I will be using > 0 and < 0
OS_CAFs <- mutate(OS_CAFs, Cat=ifelse(OS_CAFs$CAF_MCPCOUNTER > OS_CAF_Quantile[3], "> 50th", "< 50th"))

#Survival curve calculation
fit<-survfit(Surv(OS_MONTHS, OS_Status_Numeric) ~ Cat, data=OS_CAFs)
print(fit)

# Most basic plot
ggsurvplot(fit, data=OS_CAFs)

# More involved, customized plot
OS_Plot <- ggsurvplot(fit,
                      data=OS_CAFs,
                      pval = TRUE, pval.method = TRUE, conf.int = FALSE,
                      risk.table = TRUE, 
                      risk.table.col = "strata", 
                      legend = c(.80,.50),
                      pval.coord = c(.15,.45),
                      pval.size = 8,
                      break.time.by = 24,
                      pval.method.coord = c(.15,.60),
                      font.x = c(23, "bold", "black"),
                      font.y = c(23, "bold", "black"),
                      font.title = c(25, "bold", "black"),
                      # linetype = "strata", # Change line type by groups
                      # surv.median.line = "hv", # Specify median survival
                      ggtheme = theme_classic(),
                      font.tickslab=c(15, "bold", "black"),
                      font.legend = c(20, "bold", "black"),
                      risk.table.fontsize = 8,
                      risk.table.y.text = FALSE,
                      risk.table.legend = FALSE,
                      legend.title = "CAF MCPCOUNTER Level",
                      title = "CAF MCPCOUNTER OS",
                      legend.labs = c("Low CAF", "High CAF"),
                      tables.height = .33,
                      palette = c("#2E9FDF", "red", "#86AA00"))+
  labs(y = "Disease Free Survival", x = "Months")  
OS_Plot$plot <- OS_Plot$plot + theme(panel.border = element_rect(colour = "black", size = 4, fill = NA))
OS_Plot$table <- OS_Plot$table + 
  theme(panel.border = element_rect(colour = "black", size = 4, fill = NA),
        plot.title = element_text(size = 23, color = "black", face = "bold"),
        axis.title = element_text(size = 23, color = "black", face = "bold"),
        axis.text = element_text(size = 15, color = "black", face = "bold"),
        legend.position = "none")+
  labs(y = NULL) 
OS_Plot

# Saving the plots with ggsave
p1 = OS_Plot$plot
p2 = OS_Plot$table
plotp = cowplot::plot_grid(p1,p2,align = "v",axis = "b", ncol =1, rel_heights = c(3,1))
ggsave("outputs/22-0531_TCGA_TIMER_Exploration/22-0531_TCGA_CAF_OS_MCPCOUNTER.png",
       plot = plotp, 
       width = 7, height = 7,
       dpi = 600)
