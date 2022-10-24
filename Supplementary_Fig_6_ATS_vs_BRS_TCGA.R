# Author: Matthew Loberg
# Date: 22-0617
# Script: Supplemental Figure 6 (TCGA Violins)

## Description
# Script for making the plots in supplemental Figure 6 - TCGA Violins comparing ATS and BRS

## 22-0804 Update
# using bonferroni's correction instead of benjamini-hochberg procedure

## 22-0930 Update
# Using new overlap genes for the aggressive tumor score (ATS) definition
# Should be an input of 549 genes

# Load required packages
library(tidyverse)

#### Load TCGA Patient and Sample Data ####

TCGA_Patient_Data <- read.table(file = "data_in_use/22-0107_TCGA_data_clinical_patient.txt", sep = '\t', header = TRUE) # Load patient data file
TCGA_Patient_Data <- as_tibble(TCGA_Patient_Data) # Make patient data into tibble for easier data cleaning

TCGA_Sample_Data <- read.table(file = "data_in_use/22-0107_TCGA_data_clinical_sample.txt", sep = '\t', header = TRUE) # Load sample data file
TCGA_Sample_Data <- as_tibble(TCGA_Sample_Data) # Make sample data into tibble for easier data cleaning

#### Merge TCGA Patient and Sample Data ####

# PATIEND_ID has an "01" appended to the end of each sample. Here, I will remove it from each sample. 
TCGA_Sample_Data$PATIENT_ID <- substr(TCGA_Sample_Data$PATIENT_ID, 1, nchar(TCGA_Sample_Data$PATIENT_ID)-3)
# TCGA_Sample_Data_Restricted$PATIENT_ID <- substr(TCGA_Sample_Data_Restricted$PATIENT_ID, 1, nchar(TCGA_Sample_Data_Restricted$PATIENT_ID)-3) # Old line of code from when I had a restricted file
# Combine sample data and patient data
Patient_Sample_Combined <- as_tibble(merge(TCGA_Sample_Data, TCGA_Patient_Data, by.x = "PATIENT_ID", by.y = "PATIENT_ID", all = T))

Patient_Sample_Combined <- Patient_Sample_Combined %>% subset(RNASEQ_DATA == "Yes")

# Clean up environment - only leave the data that I need (Patient_Sample_Combined)
remove(list = c("TCGA_Patient_Data", "TCGA_Sample_Data"))

#### Load TCGA sequencing RSEM Expression Data ####

# Read in TCGA RNA Data (Z-Score format) and clean the data to the desired format
TCGA_RNA_Z_Scores <- read.table(file = "data_in_use/data_RNA_Seq_v2_mRNA_median_Zscores.txt", sep = '\t', header = FALSE) # Read in TCGA RNA data in Z-Score format
Gene_Symbols <- TCGA_RNA_Z_Scores$V1 # Create variable for storing Gene Symbols -> will become column names
TCGA_RNA_Z_Scores <- TCGA_RNA_Z_Scores[,3:ncol(TCGA_RNA_Z_Scores)] # Delete first two columns, don't need these
TCGA_Transposed <- t(TCGA_RNA_Z_Scores) # Transpose the data so that each column pertains to a particular gene
colnames(TCGA_Transposed) <- Gene_Symbols # Make Gene Symbols the column names
TCGA_RNA_Tibble <- as_tibble(TCGA_Transposed) # Convert type to tibble
TCGA_RNA_Tibble <- dplyr::rename(TCGA_RNA_Tibble, SAMPLE_ID = Hugo_Symbol) # dplyr::rename "Hugo_Symbol" to "SAMPLE_ID" to match Patient_Sample_Combined for future merge

# Cleaning up: remove the data I no longer need
remove(list = c("TCGA_RNA_Z_Scores", "TCGA_Transposed", "Gene_Symbols"))

#### Format poor outcome gene list (526) for TCGA Analysis ####
# Read in 526 gene BRAF-Poor outcome list generated on 22-0617
PoorOutcome_Genes <- read_csv(file = "data_in_use/22-0927_All_BRAF_All_Aggressive_Overlap_Genes.csv")

# I need to make sure that this gene set is compatible with the TCGA data
# I will start by testing the PoorOutcome gene set to see if all of the genes are present
Test <- TCGA_RNA_Tibble[c(PoorOutcome_Genes$value)]
PoorOutcome_Genes_Filtered <- PoorOutcome_Genes
Test <- TCGA_RNA_Tibble[c(PoorOutcome_Genes_Filtered$value)]

# No H1-5 -> found old symbol name, HIST1H1B
TCGA_RNA_Tibble <- TCGA_RNA_Tibble %>% dplyr::rename('H1-5' = HIST1H1B)
# No H3C2 -> found old symbol name, HIST1H3B
TCGA_RNA_Tibble <- TCGA_RNA_Tibble %>% dplyr::rename(H3C2 = HIST1H3B)
# No H2BC17 -> found old symbol name, HIST1H2BO
TCGA_RNA_Tibble <- TCGA_RNA_Tibble %>% dplyr::rename(H2BC17 = HIST1H2BO)
# No PLPP4 -> found old symbol name, PPAPDC1A
TCGA_RNA_Tibble <- TCGA_RNA_Tibble %>% dplyr::rename(PLPP4 = PPAPDC1A)
# No H2AC16 -> found old symbol name, HIST1H2AL
TCGA_RNA_Tibble <- TCGA_RNA_Tibble %>% dplyr::rename(H2AC16 = HIST1H2AL)
# No ACP7 -> found old symbol name, PAPL
TCGA_RNA_Tibble <- TCGA_RNA_Tibble %>% dplyr::rename(ACP7 = PAPL)
# No DISP3 -> found old symbol name, PTCHD2
TCGA_RNA_Tibble <- TCGA_RNA_Tibble %>% dplyr::rename(DISP3 = PTCHD2)
# No SERTM2 - Entrez 401613 - not in data
PoorOutcome_Genes_Filtered <- PoorOutcome_Genes_Filtered %>% filter(value != "SERTM2")
# LRRC38 - Entrez 126755 - not in data
PoorOutcome_Genes_Filtered <- PoorOutcome_Genes_Filtered %>% filter(value != "LRRC38")
# ARL14EPL - Entrez 644100 - not in data
PoorOutcome_Genes_Filtered <- PoorOutcome_Genes_Filtered %>% filter(value != "ARL14EPL")
# No NT5DC4 - Entrez 284958 - not in data
PoorOutcome_Genes_Filtered <- PoorOutcome_Genes_Filtered %>% filter(value != "NT5DC4")
# No RP11-314A5.7 - could not find any info about this gene?? Removing for now
PoorOutcome_Genes_Filtered <- PoorOutcome_Genes_Filtered %>% filter(value != "RP11-314A5.7")
# No ADGRG7 -> found old symbol name, GPR128
TCGA_RNA_Tibble <- TCGA_RNA_Tibble %>% dplyr::rename(ADGRG7 = GPR128)
# No LORICRIN -> found old symbol name, LOR
TCGA_RNA_Tibble <- TCGA_RNA_Tibble %>% dplyr::rename(LORICRIN = LOR)
# No IZUMO3 -> Entrez 100129669 - not in data
PoorOutcome_Genes_Filtered <- PoorOutcome_Genes_Filtered %>% filter(value != "IZUMO3")
# No OPRPN -> found old symbol name, PROL1
TCGA_RNA_Tibble <- TCGA_RNA_Tibble %>% dplyr::rename(OPRPN = PROL1)
# No STRIT1 -> Entrez 100507537 - not in data
PoorOutcome_Genes_Filtered <- PoorOutcome_Genes_Filtered %>% filter(value != "STRIT1")
# No GPR33 -> Entrez 2856 -> not in data
PoorOutcome_Genes_Filtered <- PoorOutcome_Genes_Filtered %>% filter(value != "GPR33")
# No OR14L1P -> Entrez 127617 -> not in data
PoorOutcome_Genes_Filtered <- PoorOutcome_Genes_Filtered %>% filter(value != "OR14L1P")
# No KLF18 -> Entrez 105378952 -> not in data
PoorOutcome_Genes_Filtered <- PoorOutcome_Genes_Filtered %>% filter(value != "KLF18")
# No CCDC197 -> found old symbol name, LINC00521
TCGA_RNA_Tibble <- TCGA_RNA_Tibble %>% dplyr::rename(CCDC197 = LINC00521)
# No SMIM31 -> Entrez 100505989 -> not in data
PoorOutcome_Genes_Filtered <- PoorOutcome_Genes_Filtered %>% filter(value != "SMIM31")
# No TAF11L11 -> Entrez 112488746 -> not in data
PoorOutcome_Genes_Filtered <- PoorOutcome_Genes_Filtered %>% filter(value != "TAF11L11")
# No OOSP1 -> Entrez 255649 -> not in data
PoorOutcome_Genes_Filtered <- PoorOutcome_Genes_Filtered %>% filter(value != "OOSP1")
# No OR1R1P -> Entrez 9596 -> not in data
PoorOutcome_Genes_Filtered <- PoorOutcome_Genes_Filtered %>% filter(value != "OR1R1P")
# No TAFA4 -> found old symbol name, FAM19A4
TCGA_RNA_Tibble <- TCGA_RNA_Tibble %>% dplyr::rename(TAFA4 = FAM19A4)
# No GOLGA6L22 -> Entrez 440243 -> not in data
PoorOutcome_Genes_Filtered <- PoorOutcome_Genes_Filtered %>% filter(value != "GOLGA6L22")
# No OR4C5 -> Entrez 79346 -> not in data
PoorOutcome_Genes_Filtered <- PoorOutcome_Genes_Filtered %>% filter(value != "OR4C5")
# No KRTAP2-3 -> Entrez 730755 -> not in data
PoorOutcome_Genes_Filtered <- PoorOutcome_Genes_Filtered %>% filter(value != "KRTAP2-3")
# No PNMA6F -> Entrez 105373377 -> not in data
PoorOutcome_Genes_Filtered <- PoorOutcome_Genes_Filtered %>% filter(value != "PNMA6F")
# No PRSS56 -> Entrez 646960 -> not in data
PoorOutcome_Genes_Filtered <- PoorOutcome_Genes_Filtered %>% filter(value != "PRSS56")
# No ACOD1 -> Entrez 730249 -> not in data
PoorOutcome_Genes_Filtered <- PoorOutcome_Genes_Filtered %>% filter(value != "ACOD1")
# No PLPPR5 -> found alias symbol name, LPPR5
TCGA_RNA_Tibble <- TCGA_RNA_Tibble %>% dplyr::rename(PLPPR5 = LPPR5)
# No OR8U3 -> found old symbol name, OR5R1
TCGA_RNA_Tibble <- TCGA_RNA_Tibble %>% dplyr::rename(OR8U3 = OR5R1)
# No OR5BS1P -> Entrez 390313 -> not in data
PoorOutcome_Genes_Filtered <- PoorOutcome_Genes_Filtered %>% filter(value != "OR5BS1P")
# No RP11-388F6.5 -> could not find Entrez ID
PoorOutcome_Genes_Filtered <- PoorOutcome_Genes_Filtered %>% filter(value != "RP11-388F6.5")
# No LACTBL1 -> Entrez 646262 -> not in data
PoorOutcome_Genes_Filtered <- PoorOutcome_Genes_Filtered %>% filter(value != "LACTBL1")
# No LY6L -> Entrez 101928108 -> not in data
PoorOutcome_Genes_Filtered <- PoorOutcome_Genes_Filtered %>% filter(value != "LY6L")
# No H2AC4 -> found old symbol name, HIST1H2AB
TCGA_RNA_Tibble <- TCGA_RNA_Tibble %>% dplyr::rename(H2AC4 = HIST1H2AB)
# No H3C12 -> found old symbol name, HIST1H3J
TCGA_RNA_Tibble <- TCGA_RNA_Tibble %>% dplyr::rename(H3C12 = HIST1H3J)
# No H3C3 -> found old symbol name, HIST1H3C
TCGA_RNA_Tibble <- TCGA_RNA_Tibble %>% dplyr::rename(H3C3 = HIST1H3C)
# No ADGRF4 -> found old symbol name, GPR115
TCGA_RNA_Tibble <- TCGA_RNA_Tibble %>% dplyr::rename(ADGRF4 = GPR115)
# No CALHM6 -> found old symbol name, FAM26F
TCGA_RNA_Tibble <- TCGA_RNA_Tibble %>% dplyr::rename(CALHM6 = FAM26F)
# No MDFIC2 -> Entrez 101928108 -> not in data (need to double check this)
PoorOutcome_Genes_Filtered <- PoorOutcome_Genes_Filtered %>% filter(value != "MDFIC2")
# No RTL3 -> found old symbol name, ZCCHC5
TCGA_RNA_Tibble <- TCGA_RNA_Tibble %>% dplyr::rename(RTL3 = ZCCHC5)
# No RFPL4AL1 -> Entrez 101928108 -> not in data (need to double check this)
PoorOutcome_Genes_Filtered <- PoorOutcome_Genes_Filtered %>% filter(value != "RFPL4AL1")

Test <- TCGA_RNA_Tibble[c(PoorOutcome_Genes_Filtered$value)]
rm(Test)

# Restrict TCGA data by PoorOutcome and make into a numeric
TCGA_RNA_Tibble_PoorOutcome <- TCGA_RNA_Tibble[c(PoorOutcome_Genes_Filtered$value)] # Restrict genes
PoorOutcome_Gene_Names <- colnames(TCGA_RNA_Tibble_PoorOutcome) # save col names for future use
TCGA_RNA_Matrix_PoorOutcome <- as.matrix(TCGA_RNA_Tibble_PoorOutcome) # Make into a matrix
TCGA_RNA_Matrix_PoorOutcome <- matrix(as.numeric(TCGA_RNA_Matrix_PoorOutcome), ncol = ncol(TCGA_RNA_Matrix_PoorOutcome)) # Make the matrix a numeric
TCGA_RNA_Tibble_PoorOutcome <- as_tibble(TCGA_RNA_Matrix_PoorOutcome) # Make back into a tibble (now it will have numeric columns instead of chr)
colnames(TCGA_RNA_Tibble_PoorOutcome) <- PoorOutcome_Gene_Names # Add back column names
TCGA_RNA_Tibble_PoorOutcome <- TCGA_RNA_Tibble_PoorOutcome %>% add_column(SAMPLE_ID = TCGA_RNA_Tibble$SAMPLE_ID, .before = "ADAMTS14") # KCNN4 is the first gene in PoorOutcome. Adding Sample ID as a column before that

# cleaning up
rm(PoorOutcome_Genes, PoorOutcome_Genes_Filtered, TCGA_RNA_Matrix_PoorOutcome, TCGA_RNA_Tibble)

#### Calculate BRAF-Poor Outcome score in TCGA cohort ####
# So, to do this I will do the sum of the Z-Scores divided by the number of genes
TCGA_RNA_Tibble_PoorOutcome <- TCGA_RNA_Tibble_PoorOutcome %>% rowwise() %>% mutate(PoorOutcome_Score = mean(c_across(where(is.numeric))))
TCGA_RNA_Tibble_PoorOutcome$PoorOutcome_Score
# For some reason this is causing "NAs"
# I think this is probably because some value within the data has NA values
# Testing for columns with NAs:
total <- 0
for(i in 1:ncol(TCGA_RNA_Tibble_PoorOutcome)){
  if(is.na(TCGA_RNA_Tibble_PoorOutcome[1,i])){
    print(i)
    total <- total+1
  }
}
total
# I found that the following columns have NAs: 
# ANXA8L1
# OR8D1
# PRAMEF14
# OR5AS1
# PRAMEF20
# LYPD8
# OR5K4
# SPPL2C
# DAOA
# TXNDC8
# XAGE5
# PRAMEF18
# OR8U3
# OR6C1
# RTP5

# new aggrissive: DEFB121

# I will remove all of the columns with NAs: 
TCGA_RNA_Tibble_PoorOutcome <- TCGA_RNA_Tibble_PoorOutcome %>% subset(select = -c(ANXA8L1,
                                                                                  # OR8D1, from old 1695 set
                                                                                  PRAMEF14,
                                                                                  # OR5AS1, from old 1695 set
                                                                                  PRAMEF20,
                                                                                  LYPD8,
                                                                                  OR5K4,
                                                                                  #SPPL2C, from old 526 set
                                                                                  #DAOA, from old 526 set
                                                                                  TXNDC8,
                                                                                  XAGE5,
                                                                                  PRAMEF18,
                                                                                  OR8U3,
                                                                                  OR6C1,
                                                                                  DEFB121,
                                                                                  # RTP5, from old 1695 set
                                                                                  PoorOutcome_Score))

# Recalculate Poor Ourcome Score after removing NAs
TCGA_RNA_Tibble_PoorOutcome <- TCGA_RNA_Tibble_PoorOutcome %>% rowwise() %>% mutate(PoorOutcome_Score = mean(c_across(where(is.numeric))))
TCGA_RNA_Tibble_PoorOutcome$PoorOutcome_Score

# Create a new tibble with just the PoorOutcome Score and the sample ID
TCGA_Tibble_PoorOutcome <- TCGA_RNA_Tibble_PoorOutcome[c("SAMPLE_ID", "PoorOutcome_Score")]

#### Merge Poor Outcome Score with patient and sample data ####
# Merge this tibble with the clinical patient data
TCGA_Tibble_PoorOutcome <- TCGA_Tibble_PoorOutcome %>% merge(Patient_Sample_Combined)

#### Disease Stage Plots: Supplemental Figure 5a and Supplemental Figure 5b ####
## Disease stage in depth 

# Will start by subsetting just to lesions with disease stage data
Disease_Stage_Tibble <- TCGA_Tibble_PoorOutcome %>% subset(DISEASE_STAGE != "" & DISEASE_STAGE != "[Not Available]")
# dplyr::rename Stage IVA and Stage IVC as just stage IV
for(i in 1:nrow(Disease_Stage_Tibble)){
  if(Disease_Stage_Tibble$DISEASE_STAGE[i] == "Stage IVA" |
     Disease_Stage_Tibble$DISEASE_STAGE[i] == "Stage IVC"){
    Disease_Stage_Tibble$DISEASE_STAGE[i] <- "Stage IV"
  }
}

# Make the plot with restricted groups
# Poor outcome score by disease stage: violin plot
plot <- ggplot(Disease_Stage_Tibble, aes(DISEASE_STAGE, PoorOutcome_Score)) + 
  geom_violin(outlier.size = -1,
              aes(fill = DISEASE_STAGE),
              alpha = 0.2,
              show.legend = FALSE,
              scale = "width") + 
  geom_jitter(aes(),
  position = position_jitter(width = 0.1, height = 0),
  size = 1.5,
  alpha = 0.7,
  show.legend = FALSE) + 
  scale_fill_manual(values = c(
    "black", # Stage I Color
    "#00BFC4", # Stage II Color
    "red", # Stage III Color
    "purple")) + # Stage IV Color
  # scale_colour_manual(values = c("black", "red")) + 
  scale_x_discrete(name = "Disease Stage", limits = c("Stage I", "Stage II", "Stage III", "Stage IV")) + 
  labs(y = "Aggressive Tumor\nScore") + 
  ggtitle("Aggressive Tumor Score\nby Disease Stage") +
  scale_y_continuous(breaks = c(0, 0.5, 1, 1.5),
                     limits = c(-0.4, 1.5)) +
  theme_classic() + 
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)  + 
  theme(    
    panel.border = element_rect(colour = "black", size = 3, fill = NA),
    plot.title = element_text(size = 30, color = "black", face = "bold", hjust = 0.5, vjust = 1.8),
    axis.title = element_text(size = 30, color = "black", face = "bold"),
    axis.text = element_text(size = 25, color = "black", face = "bold"))
plot
ggsave("outputs/22-0930_Supplemental_Fig_5_TCGA/22-0930_Sup5_TCGA_DiseaseStage_PoorOutcome_Score_jitter.png",
       plot, 
       width = 8, height = 7, dpi = 600)
# Run stats 
res.aov <- aov(PoorOutcome_Score ~ DISEASE_STAGE, data = Disease_Stage_Tibble)
summary(res.aov)
TukeyHSD(res.aov)

kruskal.test(PoorOutcome_Score ~ DISEASE_STAGE, data = Disease_Stage_Tibble) # shouldn't be necessary for only two groups
pairwise.wilcox.test(Disease_Stage_Tibble$PoorOutcome_Score, Disease_Stage_Tibble$DISEASE_STAGE, # this returned 0.2475
                     p.adjust.method = "BH")
pairwise.wilcox.test(Disease_Stage_Tibble$PoorOutcome_Score, Disease_Stage_Tibble$DISEASE_STAGE, # this returned 0.2475
                     p.adjust.method = "bonferroni")

# Make the plot with restricted groups for BRS (Instead of BRAF_PoorOutcome_Score)
# BRS by Disease Stage: violin
plot <- ggplot(Disease_Stage_Tibble, aes(DISEASE_STAGE, BRAFV600E_RAF_SCORE)) + 
  geom_violin(outlier.size = -1,
              aes(fill = DISEASE_STAGE),
              alpha = 0.2,
              show.legend = FALSE,
              scale = "width") + 
  geom_jitter(aes(),
  position = position_jitter(width = 0.1, height = 0),
  size = 1.5,
  alpha = 0.7,
  show.legend = FALSE) + 
  scale_fill_manual(values = c(
    "black", # Stage I Color
    "#00BFC4", # Stage II Color
    "red", # Stage III Color
    "purple")) + # Stage IV Color
  # scale_colour_manual(values = c("black", "red")) + 
  scale_x_discrete(name = "Disease Stage", limits = c("Stage I", "Stage II", "Stage III", "Stage IV")) + 
  labs(y = "BRAF-RAS Score") + 
  ggtitle("BRAF-RAS Score\nby Disease Stage") +
  theme_classic() + 
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)  + 
  theme(
    panel.border = element_rect(colour = "black", size = 3, fill = NA),
    plot.title = element_text(size = 30, color = "black", face = "bold", hjust = 0.5, vjust = 1.8),
    axis.title = element_text(size = 30, color = "black", face = "bold"),
    axis.text = element_text(size = 25, color = "black", face = "bold"))
plot
ggsave("outputs/22-0930_Supplemental_Fig_5_TCGA/22-0930_Sup5_TCGA_DiseaseStage_BRS_jitter.png",
       plot, 
       width = 8, height = 7, dpi = 600)

# Stats for BRS by Disease Stage
kruskal.test(BRAFV600E_RAF_SCORE ~ DISEASE_STAGE, data = Disease_Stage_Tibble) # shouldn't be necessary for only two groups
pairwise.wilcox.test(Disease_Stage_Tibble$BRAFV600E_RAF_SCORE, Disease_Stage_Tibble$DISEASE_STAGE, # this returned 0.2475
                     p.adjust.method = "BH")
pairwise.wilcox.test(Disease_Stage_Tibble$BRAFV600E_RAF_SCORE, Disease_Stage_Tibble$DISEASE_STAGE, # this returned 0.2475
                     p.adjust.method = "bonferroni")


#### Risk Group Plots: Supplemental Figure 5a and Supplemental Figure 5b ####
# Subset out samples without a risk group 
TCGA_Tibble_PoorOutcome_RISK <- TCGA_Tibble_PoorOutcome %>% subset(RISK_GROUP != "")

# Plot poor outcome score by risk group: violin
plot <- ggplot(TCGA_Tibble_PoorOutcome_RISK, aes(RISK_GROUP, PoorOutcome_Score)) + 
  geom_violin(outlier.size = -1,
              aes(fill = RISK_GROUP),
              alpha = 0.2,
              show.legend = FALSE,
              scale = "width") + 
  geom_jitter(aes(),
  position = position_jitter(width = 0.1, height = 0),
  size = 1.5,
  alpha = 0.7,
  show.legend = FALSE) + 
  scale_fill_manual(values = c("purple", # High Risk group color
                               "red", # Intermediate risk group color
                               "#00BFC4")) + # Low risk group color 
  # scale_colour_manual(values = c("black", "red")) + 
  scale_x_discrete(name = "Risk Group", limits = c("Low", "Intermediate", "High")) + 
  labs(y = "Aggressive Tumor\nScore") + 
  ggtitle("Aggressive Tumor Score\nby Risk Group") +
  scale_y_continuous(breaks = c(0, 0.5, 1, 1.5),
                     limits = c(-0.4, 1.5)) +
  theme_classic() + 
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)  + 
  theme(
    panel.border = element_rect(colour = "black", size = 3, fill = NA),
    plot.title = element_text(size = 30, color = "black", face = "bold", hjust = 0.5, vjust = 1.8),
    axis.title = element_text(size = 30, color = "black", face = "bold"),
    axis.text = element_text(size = 25, color = "black", face = "bold"))
plot
ggsave("outputs/22-0930_Supplemental_Fig_5_TCGA/22-0930_Sup5_TCGA_RiskGroup_Poor_Outcome_jitter.png",
       plot, 
       width = 8, height = 7, dpi = 600)

# Risk Group Stats
res.aov <- aov(PoorOutcome_Score ~ RISK_GROUP, data = TCGA_Tibble_PoorOutcome_RISK)
summary(res.aov)
TukeyHSD(res.aov)

# Stats for Poor Outcome by RISK Group
kruskal.test(PoorOutcome_Score ~ RISK_GROUP, data = TCGA_Tibble_PoorOutcome_RISK) # shouldn't be necessary for only two groups
pairwise.wilcox.test(TCGA_Tibble_PoorOutcome_RISK$PoorOutcome_Score, TCGA_Tibble_PoorOutcome_RISK$RISK_GROUP, # this returned 0.2475
                     p.adjust.method = "BH")
pairwise.wilcox.test(TCGA_Tibble_PoorOutcome_RISK$PoorOutcome_Score, TCGA_Tibble_PoorOutcome_RISK$RISK_GROUP, # this returned 0.2475
                     p.adjust.method = "bonferroni")

# BRS by RISK Group: Violin
plot <- ggplot(TCGA_Tibble_PoorOutcome_RISK, aes(RISK_GROUP, BRAFV600E_RAF_SCORE)) + 
  geom_violin(outlier.size = -1,
              aes(fill = RISK_GROUP),
              alpha = 0.2,
              show.legend = FALSE,
              scale = "width") + 
  geom_jitter(aes(),
  position = position_jitter(width = 0.1, height = 0),
  size = 1.5,
  alpha = 0.7,
  show.legend = FALSE) + 
  scale_fill_manual(values = c("purple", # High Risk group color
                               "red", # Intermediate risk group color
                               "#00BFC4")) + # Low risk group color 
  # scale_colour_manual(values = c("black", "red")) + 
  scale_x_discrete(name = "Risk Group", limits = c("Low", "Intermediate", "High")) + 
  labs(y = "BRAF-RAS Score") + 
  ggtitle("BRAF-RAS Score\nby Risk Group") +
  theme_classic() + 
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)  + 
  theme(
    panel.border = element_rect(colour = "black", size = 3, fill = NA),
    plot.title = element_text(size = 30, color = "black", face = "bold", hjust = 0.5, vjust = 1.8),
    axis.title = element_text(size = 30, color = "black", face = "bold"),
    axis.text = element_text(size = 25, color = "black", face = "bold"))
plot
ggsave("outputs/22-0930_Supplemental_Fig_5_TCGA/22-0930_Sup5_TCGA_RiskGroup_BRS_jitter.png",
       plot, 
       width = 8, height = 7, dpi = 600)

# Stats for BRS by RISK Group
kruskal.test(BRAFV600E_RAF_SCORE ~ RISK_GROUP, data = TCGA_Tibble_PoorOutcome_RISK) # shouldn't be necessary for only two groups
pairwise.wilcox.test(TCGA_Tibble_PoorOutcome_RISK$BRAFV600E_RAF_SCORE, TCGA_Tibble_PoorOutcome_RISK$RISK_GROUP, # this returned 0.2475
                     p.adjust.method = "BH")
pairwise.wilcox.test(TCGA_Tibble_PoorOutcome_RISK$BRAFV600E_RAF_SCORE, TCGA_Tibble_PoorOutcome_RISK$RISK_GROUP, # this returned 0.2475
                     p.adjust.method = "bonferroni")

#### Extrathyroidal Extension Plots: Supplemental Figure 5a and Supplemental Figure 5b ####
# Simplifying the TCGA extrathyroidal extension category
TCGA_Tibble_PoorOutcome$Extrathyroidal_Modified <- TCGA_Tibble_PoorOutcome$EXTRATHYROIDAL_EXTENSION
for(i in 1:nrow(TCGA_Tibble_PoorOutcome)){
  if(TCGA_Tibble_PoorOutcome$Extrathyroidal_Modified[i] == "Moderate/Advanced (T4a)"){
    TCGA_Tibble_PoorOutcome$Extrathyroidal_Modified[i] <- "T4"
  }
  else if(TCGA_Tibble_PoorOutcome$Extrathyroidal_Modified[i] == "Very Advanced (T4b)"){
    TCGA_Tibble_PoorOutcome$Extrathyroidal_Modified[i] <- "T4"
  }
  else if(TCGA_Tibble_PoorOutcome$Extrathyroidal_Modified[i] == "Minimal (T3)"){
    TCGA_Tibble_PoorOutcome$Extrathyroidal_Modified[i] <- "T3"
  }
}

# Extrathyroidal Modified Extension Stats
Extrathyroidal_Stats <- TCGA_Tibble_PoorOutcome %>% subset(Extrathyroidal_Modified == "None" |
                                                             Extrathyroidal_Modified == "T3" |
                                                             Extrathyroidal_Modified == "T4")


# Now I will plot the new extrathyroidal modified category by 
# Poor outcome score by extrathyroidal extension: Violin
plot <- ggplot(Extrathyroidal_Stats, aes(Extrathyroidal_Modified, PoorOutcome_Score)) + 
  geom_violin(outlier.size = -1,
              aes(fill = Extrathyroidal_Modified),
              alpha = 0.2,
              show.legend = FALSE,
              scale = "width") + 
  geom_jitter(aes(),
  position = position_jitter(width = 0.1, height = 0),
  size = 1.5,
  alpha = 0.7,
  show.legend = FALSE) + 
  scale_fill_manual(values = c("#00BFC4", "red", "purple")) + 
  # scale_colour_manual(values = c("black", "red")) + 
  scale_x_discrete(name = "Extrathyroidal Extension", limits = c("None", "T3", "T4")) + 
  labs(y = "Aggressive Tumor\nScore") + 
  ggtitle("Aggressive Tumor Score\nby Extrathyroidal Extension") +
  theme_classic() + 
  scale_y_continuous(breaks = c(0, 0.5, 1, 1.5),
                     limits = c(-0.4, 1.5)) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)  + 
  theme(
    panel.border = element_rect(colour = "black", size = 3, fill = NA),
    plot.title = element_text(size = 30, color = "black", face = "bold", hjust = 0.5, vjust = 1.8),
    axis.title = element_text(size = 30, color = "black", face = "bold"),
    axis.text = element_text(size = 25, color = "black", face = "bold"))
plot
ggsave("outputs/22-0930_Supplemental_Fig_5_TCGA/22-0930_Sup5_TCGA_ExtrathyroidalExtension_Poor_Outcome_jitter.png",
       plot, 
       width = 8, height = 7, dpi = 600)

# Stats for Poor Outcome by Extrathyroidal Extension
kruskal.test(PoorOutcome_Score ~ Extrathyroidal_Modified, data = Extrathyroidal_Stats) # shouldn't be necessary for only two groups
pairwise.wilcox.test(Extrathyroidal_Stats$PoorOutcome_Score, Extrathyroidal_Stats$Extrathyroidal_Modified, # this returned 0.2475
                     p.adjust.method = "BH")
pairwise.wilcox.test(Extrathyroidal_Stats$PoorOutcome_Score, Extrathyroidal_Stats$Extrathyroidal_Modified, # this returned 0.2475
                     p.adjust.method = "bonferroni")

# BRS Score by extrathyroidal extension: Violin
plot <- ggplot(Extrathyroidal_Stats, aes(Extrathyroidal_Modified, BRAFV600E_RAF_SCORE)) + 
  geom_violin(outlier.size = -1,
              aes(fill = Extrathyroidal_Modified),
              alpha = 0.2,
              show.legend = FALSE,
              scale = "width") + 
  geom_jitter(aes(),
  position = position_jitter(width = 0.1, height = 0),
  size = 1.5,
  alpha = 0.7,
  show.legend = FALSE) + 
  scale_fill_manual(values = c("#00BFC4", "red", "purple")) + 
  # scale_colour_manual(values = c("black", "red")) + 
  scale_x_discrete(name = "Extrathyroidal Extension", limits = c("None", "T3", "T4")) + 
  labs(y = "BRAF-RAS Score") + 
  ggtitle("BRAF-RAS Score\nby Extrathyroidal Extension") +
  theme_classic() + 
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)  + 
  theme(    
    panel.border = element_rect(colour = "black", size = 3, fill = NA),
    plot.title = element_text(size = 30, color = "black", face = "bold", hjust = 0.5, vjust = 1.8),
    axis.title = element_text(size = 30, color = "black", face = "bold"),
    axis.text = element_text(size = 25, color = "black", face = "bold"))
plot
ggsave("outputs/22-0930_Supplemental_Fig_5_TCGA/22-0930_Sup5_TCGA_ExtrathyroidalExtension_BRS_jitter.png",
       plot, 
       width = 8, height = 7, dpi = 600)

# Stats for BRS by Extrathyroidal Extension
kruskal.test(BRAFV600E_RAF_SCORE ~ Extrathyroidal_Modified, data = Extrathyroidal_Stats) # shouldn't be necessary for only two groups
pairwise.wilcox.test(Extrathyroidal_Stats$BRAFV600E_RAF_SCORE, Extrathyroidal_Stats$Extrathyroidal_Modified, # this returned 0.2475
                     p.adjust.method = "BH")
pairwise.wilcox.test(Extrathyroidal_Stats$BRAFV600E_RAF_SCORE, Extrathyroidal_Stats$Extrathyroidal_Modified, # this returned 0.2475
                     p.adjust.method = "bonferroni")

#### Histology Type Plots: Supplemental Figure 5c and Supplemental Figure 5d ####
# Here I will look at histological type and how it compares to both the PoorOutcome score and how it compares to BRS
TCGA_Tibble_PoorOutcome$HISTOLOGICAL_TYPE_COMPLEX <- TCGA_Tibble_PoorOutcome$HISTOLOGICAL_TYPE
for(i in 1:nrow(TCGA_Tibble_PoorOutcome)){
  if(TCGA_Tibble_PoorOutcome$HISTOLOGICAL_TYPE[i] == ""){
    TCGA_Tibble_PoorOutcome$HISTOLOGICAL_TYPE_COMPLEX[i] <- "Unknown\nNot Specified"
  }
  else if(TCGA_Tibble_PoorOutcome$HISTOLOGICAL_TYPE_OTHER[i] == "Diffuse Sclerosing Variant" |
          TCGA_Tibble_PoorOutcome$HISTOLOGICAL_TYPE_OTHER[i] == "Diffuse sclerosing variant" |
          TCGA_Tibble_PoorOutcome$HISTOLOGICAL_TYPE_OTHER[i] == "Papillary carcinoma, diffuse sclerosing variant"){
    TCGA_Tibble_PoorOutcome$HISTOLOGICAL_TYPE_COMPLEX[i] <- "Diffuse\nSclerosing"
  }
  else if(TCGA_Tibble_PoorOutcome$HISTOLOGICAL_TYPE_OTHER[i] == "cribriform morular" |
          TCGA_Tibble_PoorOutcome$HISTOLOGICAL_TYPE_OTHER[i] == "cribiform morular"){
    TCGA_Tibble_PoorOutcome$HISTOLOGICAL_TYPE_COMPLEX[i] <- "Cribriform\nMorular"
  }
  else if(TCGA_Tibble_PoorOutcome$HISTOLOGICAL_TYPE_OTHER[i] == "encapsulated follicular"){
    TCGA_Tibble_PoorOutcome$HISTOLOGICAL_TYPE_COMPLEX[i] <- "Follicular" 
  }
  else if(TCGA_Tibble_PoorOutcome$HISTOLOGICAL_TYPE[i] == "Tall Cell"){
    TCGA_Tibble_PoorOutcome$HISTOLOGICAL_TYPE_COMPLEX[i] <- "Tall\nCell"
  }
}

# Subset the data to just the following groups: 
# Follicular 
# Diffuse\nSclerosing
# Tall Cell
# Classical
Histology_Plot_Cohort <- TCGA_Tibble_PoorOutcome %>% subset(HISTOLOGICAL_TYPE_COMPLEX == "Classical" |
                                                              HISTOLOGICAL_TYPE_COMPLEX == "Diffuse\nSclerosing" |
                                                              HISTOLOGICAL_TYPE_COMPLEX == "Follicular" | 
                                                              HISTOLOGICAL_TYPE_COMPLEX == "Tall\nCell")

### Make a new plot combining diffuse sclerosing and tall cell
Histology_Plot_Cohort$HISTOLOGICAL_TYPE_SIMPLIFIED <- Histology_Plot_Cohort$HISTOLOGICAL_TYPE_COMPLEX
for(i in 1:nrow(Histology_Plot_Cohort)){
  if(Histology_Plot_Cohort$HISTOLOGICAL_TYPE_COMPLEX[i] == "Diffuse\nSclerosing" |
     Histology_Plot_Cohort$HISTOLOGICAL_TYPE_COMPLEX[i] == "Tall\nCell"){
    Histology_Plot_Cohort$HISTOLOGICAL_TYPE_SIMPLIFIED[i] <- "Diffuse\nSclerosing\n+\nTall Cell"
  }
}

# Plot new histology definitions
# Poor outcome score by histology simplified: Violin
plot <- ggplot(Histology_Plot_Cohort, aes(HISTOLOGICAL_TYPE_SIMPLIFIED, PoorOutcome_Score)) + 
  geom_violin(outlier.size = -1,
              aes(fill = HISTOLOGICAL_TYPE_SIMPLIFIED),
              alpha = 0.4,
              show.legend = FALSE,
              scale = "width") + 
  geom_jitter(aes(),
  position = position_jitter(width = 0.1, height = 0),
  size = 1.5,
  alpha = c(0.8),
  show.legend = FALSE) +
  scale_fill_manual(values = c("red", "purple", "blue")) + 
  scale_x_discrete(name = "Histology", limits = c("Follicular", "Classical", "Diffuse\nSclerosing\n+\nTall Cell")) + 
  labs(y = "Aggressive Tumor\nScore") + 
  ggtitle("Aggressive Tumor Score\nby Histology") +
  scale_y_continuous(breaks = c(0, 0.5, 1, 1.5),
                     limits = c(-.4, 1.5)) +
  theme_classic() + 
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2) + 
  theme(
    panel.border = element_rect(colour = "black", size = 3, fill = NA),    
    plot.title = element_text(size = 30, color = "black", face = "bold"),
    axis.title = element_text(size = 30, color = "black", face = "bold"),
    axis.text = element_text(size = 25, color = "black", face = "bold"))
plot
ggsave("outputs/22-0930_Supplemental_Fig_5_TCGA/22-0930_Sup5_TCGA_Histology_Poor_Outcome_jitter.png",
       plot, 
       width = 8, height = 7, dpi = 600)

# Histology Plot Stats
# ANOVA old
# res.aov <- aov(PoorOutcome_Score ~ HISTOLOGICAL_TYPE_SIMPLIFIED, data = Histology_Plot_Cohort)
# summary(res.aov)
# TukeyHSD(res.aov)

# Stats for Poor Outcome Score by Histology
kruskal.test(PoorOutcome_Score ~ HISTOLOGICAL_TYPE_SIMPLIFIED, data = Histology_Plot_Cohort) # shouldn't be necessary for only two groups
pairwise.wilcox.test(Histology_Plot_Cohort$PoorOutcome_Score, Histology_Plot_Cohort$HISTOLOGICAL_TYPE_SIMPLIFIED, # this returned 0.2475
                     p.adjust.method = "BH")
pairwise.wilcox.test(Histology_Plot_Cohort$PoorOutcome_Score, Histology_Plot_Cohort$HISTOLOGICAL_TYPE_SIMPLIFIED, # this returned 0.2475
                     p.adjust.method = "bonferroni")

# Repeat the above histology plot with BRAF-RAS Score
# BRS Score by Histotype: Violin
plot <- ggplot(Histology_Plot_Cohort, aes(HISTOLOGICAL_TYPE_SIMPLIFIED, BRAFV600E_RAF_SCORE)) + 
  geom_violin(outlier.size = -1,
              aes(fill = HISTOLOGICAL_TYPE_SIMPLIFIED),
              alpha = 0.4,
              show.legend = FALSE,
              scale = "width") + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5,
              alpha = c(0.8),
              show.legend = FALSE) +
  scale_fill_manual(values = c("red", "purple", "blue")) + 
  scale_x_discrete(name = "Histology", limits = c("Follicular", "Classical", "Diffuse\nSclerosing\n+\nTall Cell")) + 
  labs(y = "BRAF-RAS Score") + 
  ggtitle("BRAF-RAS Score\nby Histology") +
  theme_classic() + 
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2) + 
  theme(
    panel.border = element_rect(colour = "black", size = 3, fill = NA),
    plot.title = element_text(size = 30, color = "black", face = "bold"),
    axis.title = element_text(size = 30, color = "black", face = "bold"),
    axis.text = element_text(size = 25, color = "black", face = "bold"))
plot
ggsave("outputs/22-0930_Supplemental_Fig_5_TCGA/22-0930_Sup5_TCGA_Histology_BRS_jitter.png",
       plot, 
       width = 8, height = 7, dpi = 600)

# Stats for BRS by Histology
kruskal.test(BRAFV600E_RAF_SCORE ~ HISTOLOGICAL_TYPE_SIMPLIFIED, data = Histology_Plot_Cohort) # shouldn't be necessary for only two groups
pairwise.wilcox.test(Histology_Plot_Cohort$BRAFV600E_RAF_SCORE, Histology_Plot_Cohort$HISTOLOGICAL_TYPE_SIMPLIFIED, # this returned 0.2475
                     p.adjust.method = "BH")
pairwise.wilcox.test(Histology_Plot_Cohort$BRAFV600E_RAF_SCORE, Histology_Plot_Cohort$HISTOLOGICAL_TYPE_SIMPLIFIED, # this returned 0.2475
                     p.adjust.method = "bonferroni")
