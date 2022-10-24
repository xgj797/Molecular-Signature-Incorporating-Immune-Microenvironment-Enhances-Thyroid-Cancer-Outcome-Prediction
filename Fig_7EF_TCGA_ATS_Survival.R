# Author: Matthew A Loberg
# Date: 22-0913
# Description: TCGA Survival Curves by ATS

# 22-0913 update
# using new ATS score

# Load required packages
library(ggplot2)
library(survival)
library(survminer)
library(tidyverse)
library(cowplot) # Using cowplot due to ggsave issues with survival curves

### Load TCGA Patient and Sample Data ###
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

# DFS_STATUS needs to be numeric, so separate column
Patient_Sample_Combined <-separate(data = Patient_Sample_Combined, col = DFS_STATUS, into = c("DFS_Status_Numeric", "DFS_Status_Char"), sep = ":")
Patient_Sample_Combined$DFS_Status_Numeric <- as.numeric(Patient_Sample_Combined$DFS_Status_Numeric) # Change to numeric

# OS_STATUS needs to be numeric, so separate column
Patient_Sample_Combined <-separate(data = Patient_Sample_Combined, col = OS_STATUS, into = c("OS_Status_Numeric", "OS_Status_Char"), sep = ":")
Patient_Sample_Combined$OS_Status_Numeric <- as.numeric(Patient_Sample_Combined$OS_Status_Numeric) # Change to numeric

#### Read in TCGA RNA-Sequencing Data in Z-Score format ####
# Read in TCGA RNA Data (Z-Score format) and clean the data to the desired format
TCGA_RNA_Z_Scores <- read.table(file = "data_in_use/data_RNA_Seq_v2_mRNA_median_Zscores.txt", sep = '\t', header = FALSE) # Read in TCGA RNA data in Z-Score format
Gene_Symbols <- TCGA_RNA_Z_Scores$V1 # Create variable for storing Gene Symbols -> will become column names
TCGA_RNA_Z_Scores <- TCGA_RNA_Z_Scores[,3:ncol(TCGA_RNA_Z_Scores)] # Delete first two columns, don't need these
TCGA_Transposed <- t(TCGA_RNA_Z_Scores) # Transpose the data so that each column pertains to a particular gene
colnames(TCGA_Transposed) <- Gene_Symbols # Make Gene Symbols the column names
TCGA_RNA_Tibble <- as_tibble(TCGA_Transposed) # Convert type to tibble
TCGA_RNA_Tibble <- rename(TCGA_RNA_Tibble, SAMPLE_ID = Hugo_Symbol) # Rename "Hugo_Symbol" to "SAMPLE_ID" to match Patient_Sample_Combined for future merge

# Cleaning up: remove the data I no longer need
remove(list = c("TCGA_RNA_Z_Scores", "TCGA_Transposed", "Gene_Symbols"))

#### Calculating BRAF-PoorOutcome Score ####
### Read in Gene Set of interest ###
# 1695 genes generated on 22-0331
# Today updating to 526 genes generated on 22-0426 (old as of 22-0913)
# Today (22-0913) updating to new ATS score genes
PoorOutcome_Genes <- read_csv(file = "data_in_use/22-0426_All_BRAF_All_PoorOutcome_Overlap_Genes.csv")
# New gene set
PoorOutcome_Genes <- read_csv(file = "data_in_use/22-0908_All_BRAF_All_Aggressive_Overlap_Genes.csv")

# I need to make sure that this gene set is compatible with the TCGA data
# I will start by testing the PoorOutcome gene set to see if all of the genes are present
Test <- TCGA_RNA_Tibble[c(PoorOutcome_Genes$value)]
PoorOutcome_Genes_Filtered <- PoorOutcome_Genes
Test <- TCGA_RNA_Tibble[c(PoorOutcome_Genes_Filtered$value)]

# No H1-5 -> found old symbol name, HIST1H1B
TCGA_RNA_Tibble <- TCGA_RNA_Tibble %>% rename('H1-5' = HIST1H1B)
# No H3C2 -> found old symbol name, HIST1H3B
TCGA_RNA_Tibble <- TCGA_RNA_Tibble %>% rename(H3C2 = HIST1H3B)
# No H2BC17 -> found old symbol name, HIST1H2BO
TCGA_RNA_Tibble <- TCGA_RNA_Tibble %>% rename(H2BC17 = HIST1H2BO)
# No PLPP4 -> found old symbol name, PPAPDC1A
TCGA_RNA_Tibble <- TCGA_RNA_Tibble %>% rename(PLPP4 = PPAPDC1A)
# No H2AC16 -> found old symbol name, HIST1H2AL
TCGA_RNA_Tibble <- TCGA_RNA_Tibble %>% rename(H2AC16 = HIST1H2AL)
# No ACP7 -> found old symbol name, PAPL
TCGA_RNA_Tibble <- TCGA_RNA_Tibble %>% rename(ACP7 = PAPL)
# No DISP3 -> found old symbol name, PTCHD2
TCGA_RNA_Tibble <- TCGA_RNA_Tibble %>% rename(DISP3 = PTCHD2)
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
TCGA_RNA_Tibble <- TCGA_RNA_Tibble %>% rename(ADGRG7 = GPR128)
# No LORICRIN -> found old symbol name, LOR
TCGA_RNA_Tibble <- TCGA_RNA_Tibble %>% rename(LORICRIN = LOR)
# No IZUMO3 -> Entrez 100129669 - not in data
PoorOutcome_Genes_Filtered <- PoorOutcome_Genes_Filtered %>% filter(value != "IZUMO3")
# No OPRPN -> found old symbol name, PROL1
TCGA_RNA_Tibble <- TCGA_RNA_Tibble %>% rename(OPRPN = PROL1)
# No STRIT1 -> Entrez 100507537 - not in data
PoorOutcome_Genes_Filtered <- PoorOutcome_Genes_Filtered %>% filter(value != "STRIT1")
# No GPR33 -> Entrez 2856 -> not in data
PoorOutcome_Genes_Filtered <- PoorOutcome_Genes_Filtered %>% filter(value != "GPR33")
# No OR14L1P -> Entrez 127617 -> not in data
PoorOutcome_Genes_Filtered <- PoorOutcome_Genes_Filtered %>% filter(value != "OR14L1P")
# No KLF18 -> Entrez 105378952 -> not in data
PoorOutcome_Genes_Filtered <- PoorOutcome_Genes_Filtered %>% filter(value != "KLF18")
# No CCDC197 -> found old symbol name, LINC00521
TCGA_RNA_Tibble <- TCGA_RNA_Tibble %>% rename(CCDC197 = LINC00521)
# No SMIM31 -> Entrez 100505989 -> not in data
PoorOutcome_Genes_Filtered <- PoorOutcome_Genes_Filtered %>% filter(value != "SMIM31")
# No TAF11L11 -> Entrez 112488746 -> not in data
PoorOutcome_Genes_Filtered <- PoorOutcome_Genes_Filtered %>% filter(value != "TAF11L11")
# No OOSP1 -> Entrez 255649 -> not in data
PoorOutcome_Genes_Filtered <- PoorOutcome_Genes_Filtered %>% filter(value != "OOSP1")
# No OR1R1P -> Entrez 9596 -> not in data
PoorOutcome_Genes_Filtered <- PoorOutcome_Genes_Filtered %>% filter(value != "OR1R1P")
# No TAFA4 -> found old symbol name, FAM19A4
TCGA_RNA_Tibble <- TCGA_RNA_Tibble %>% rename(TAFA4 = FAM19A4)
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
TCGA_RNA_Tibble <- TCGA_RNA_Tibble %>% rename(PLPPR5 = LPPR5)
# No OR8U3 -> found old symbol name, OR5R1
TCGA_RNA_Tibble <- TCGA_RNA_Tibble %>% rename(OR8U3 = OR5R1)
# No OR5BS1P -> Entrez 390313 -> not in data
PoorOutcome_Genes_Filtered <- PoorOutcome_Genes_Filtered %>% filter(value != "OR5BS1P")
# No RP11-388F6.5 -> could not find Entrez ID
PoorOutcome_Genes_Filtered <- PoorOutcome_Genes_Filtered %>% filter(value != "RP11-388F6.5")
# No LACTBL1 -> Entrez 646262 -> not in data
PoorOutcome_Genes_Filtered <- PoorOutcome_Genes_Filtered %>% filter(value != "LACTBL1")
# No LY6L -> Entrez 101928108 -> not in data
PoorOutcome_Genes_Filtered <- PoorOutcome_Genes_Filtered %>% filter(value != "LY6L")
# No H2AC4 -> found old symbol name, HIST1H2AB
TCGA_RNA_Tibble <- TCGA_RNA_Tibble %>% rename(H2AC4 = HIST1H2AB)
# No H3C12 -> found old symbol name, HIST1H3J
TCGA_RNA_Tibble <- TCGA_RNA_Tibble %>% rename(H3C12 = HIST1H3J)
# No H3C3 -> found old symbol name, HIST1H3C
TCGA_RNA_Tibble <- TCGA_RNA_Tibble %>% rename(H3C3 = HIST1H3C)
# No ADGRF4 -> found old symbol name, GPR115
TCGA_RNA_Tibble <- TCGA_RNA_Tibble %>% rename(ADGRF4 = GPR115)
# No CALHM6 -> found old symbol name, FAM26F
TCGA_RNA_Tibble <- TCGA_RNA_Tibble %>% rename(CALHM6 = FAM26F)
# No MDFIC2 -> Entrez 101928108 -> not in data (need to double check this)
PoorOutcome_Genes_Filtered <- PoorOutcome_Genes_Filtered %>% filter(value != "MDFIC2")
# No RTL3 -> found old symbol name, ZCCHC5
TCGA_RNA_Tibble <- TCGA_RNA_Tibble %>% rename(RTL3 = ZCCHC5)
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

# From Maggie Axelrod: "I think this code is just for a single gene, but for gene scores, I tend to do sum of z-scores divided by number of genes - but there are many ways of doing this."
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
# OR6C4
# OR7G2
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
                                                                                  # OR6C4, from old 1695 set
                                                                                  # OR7G2, from old 1695 set
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

# Merge this tibble with the clinical patient data
TCGA_Tibble_PoorOutcome <- TCGA_Tibble_PoorOutcome %>% merge(Patient_Sample_Combined)




#### DFS and OS Curves ####
# Looking at DFS and OS
# Make a new tibble for each of the above restricted to samples with DFS month data
DFS_PoorOutcome <- TCGA_Tibble_PoorOutcome %>% subset(!(is.na(DFS_MONTHS)))
DFS_BRAF_RAS_Score <- DFS_PoorOutcome %>% subset(BRAFV600E_RAS != "")

# Make a new tibble for each of the above restricted to samples with OS month data
OS_PoorOutcome <- TCGA_Tibble_PoorOutcome %>% subset(!(is.na(OS_MONTHS)))
OS_BRAF_RAS_Score <- OS_PoorOutcome %>% subset(BRAFV600E_RAS != "")

### Poor Outcome Disease Free Survival Plot
# Summary statistics on the DFS Samples PoorOutcomes Score
summary(DFS_PoorOutcome$PoorOutcome_Score)
DFS_PoorOutcome_Quantile <- quantile(DFS_PoorOutcome$PoorOutcome_Score)
DFS_PoorOutcome_Quantile # Print DFS PoorOutcome Quantiles
# Sort samples into greater than or less than 50th percentile based on PoorOutcome Score
# Rather than the 50th percentile, I will be using > 0 and < 0
DFS_PoorOutcome <- mutate(DFS_PoorOutcome, Cat=ifelse(DFS_PoorOutcome$PoorOutcome_Score > 0, "Positive", "Negative"))

#Survival curve calculation
fit<-survfit(Surv(DFS_MONTHS, DFS_Status_Numeric) ~ Cat, data=DFS_PoorOutcome)
print(fit)

# Most basic plot
ggsurvplot(fit, data=DFS_PoorOutcome)

# More involved, customized plot
DFS_Plot <- ggsurvplot(fit,
                       data=DFS_PoorOutcome,
                       pval = TRUE, pval.method = TRUE, conf.int = FALSE,
                       risk.table = TRUE, 
                       risk.table.col = "strata", 
                       legend = c(.80,.20),
                       pval.coord = c(.15,.18),
                       pval.size = 8,
                       break.time.by = 24,
                       pval.method.coord = c(.15,.09),
                       font.x = c(25, "bold", "black"),
                       font.y = c(25, "bold", "black"),
                       font.title = c(27, "bold", "black"),
                       # linetype = "strata", # Change line type by groups
                       # surv.median.line = "hv", # Specify median survival
                       ggtheme = theme_classic(),
                       font.tickslab=c(20, "bold", "black"),
                       font.legend = c(20, "bold", "black"),
                       risk.table.fontsize = 8,
                       risk.table.y.text = FALSE,
                       risk.table.legend = FALSE,
                       legend.title = "Mean Z-Score",
                       title = "Aggressive Thyroid Score DFS",
                       legend.labs = c("Z-Score < 0", "Z-Score > 0"),
                       tables.height = .33,
                       palette = c("mediumpurple1", # Z-Score < 0
                                   "orchid1", # Z-Score > 0
                                   "#86AA00"))+
  labs(y = "Disease Free Survival", x = "Months")  
DFS_Plot$plot <- DFS_Plot$plot + theme(panel.border = element_rect(colour = "black", size = 2, fill = NA))
DFS_Plot$table <- DFS_Plot$table + 
  theme(panel.border = element_rect(colour = "black", size = 2, fill = NA),
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
ggsave("outputs/22-0913_TCGA_Aggressive_Score_Survival/22-0913_TCGA_Aggressive_DFS_526.png", 
       plot = plotp, 
       width = 7, height = 7,
       dpi = 600)

### Poor Outcome Overall Survival Plot
# Summary statistics on the DFS Samples PoorOutcomes Score
summary(DFS_PoorOutcome$PoorOutcome_Score)
DFS_PoorOutcome_Quantile <- quantile(DFS_PoorOutcome$PoorOutcome_Score)
DFS_PoorOutcome_Quantile # Print DFS PoorOutcome Quantiles
# Sort samples into greater than or less than 50th percentile based on PoorOutcome Score
# Rather than the 50th percentile, I will be using > 0 and < 0
OS_PoorOutcome <- mutate(OS_PoorOutcome, Cat=ifelse(OS_PoorOutcome$PoorOutcome_Score > 0, "Positive", "Negative"))

#Survival curve calculation
fit<-survfit(Surv(OS_MONTHS, OS_Status_Numeric) ~ Cat, data=OS_PoorOutcome)
print(fit)

# Most basic plot
ggsurvplot(fit, data=OS_PoorOutcome)

# More involved, customized plot
OS_Plot <- ggsurvplot(fit,
                      data=OS_PoorOutcome,
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
                      ggtheme = theme_bw(),
                      font.tickslab=c(15, "bold", "black"),
                      font.legend = c(20, "bold", "black"),
                      risk.table.fontsize = 8,
                      risk.table.y.text = FALSE,
                      risk.table.legend = FALSE,
                      legend.title = "Mean Z-Score",
                      title = "Aggressive Thyroid Score OS",
                      legend.labs = c("Z-Score < 0", "Z-Score > 0"),
                      tables.height = .33,
                      palette = c("mediumpurple1", # Z-Score < 0
                                  "orchid1", # Z-Score > 0
                                  "#86AA00"))+
  labs(y = "Overall Survival", x = "Months")  
OS_Plot$table <- OS_Plot$table + 
  theme(plot.title = element_text(size = 23, color = "black", face = "bold"),
        axis.title = element_text(size = 23, color = "black", face = "bold"),
        axis.text = element_text(size = 15, color = "black", face = "bold"),
        legend.position = "none")+
  labs(y = NULL) 
OS_Plot

# Saving the plots with ggsave
p1 = OS_Plot$plot
p2 = OS_Plot$table
plotp = cowplot::plot_grid(p1,p2,align = "v",axis = "b", ncol =1, rel_heights = c(3,1))
ggsave("outputs/22-0913_TCGA_ATS_Score_Survival/22-0913_TCGA_ATS_OS.png", 
       plot = plotp, 
       width = 7, height = 7,
       dpi = 600)

### BRAF-RAS Disease Free Survival Plot
#Disease Free Survival curve calculation
fit<-survfit(Surv(DFS_MONTHS, DFS_Status_Numeric) ~ BRAFV600E_RAS, data=DFS_BRAF_RAS_Score)
print(fit)

# Most basic plot
ggsurvplot(fit, data=DFS_BRAF_RAS_Score)

# More involved, customized plot
DFS_Plot <- ggsurvplot(fit,
                       data=DFS_BRAF_RAS_Score,
                       pval = TRUE, pval.method = TRUE, conf.int = FALSE,
                       risk.table = TRUE, 
                       risk.table.col = "strata", 
                       legend = c(.74,.20),
                       pval.coord = c(.15,.18),
                       pval.size = 8,
                       break.time.by = 24,
                       pval.method.coord = c(.15,.09),
                       font.x = c(25, "bold", "black"),
                       font.y = c(25, "bold", "black"),
                       font.title = c(27, "bold", "black"),
                       # linetype = "strata", # Change line type by groups
                       # surv.median.line = "hv", # Specify median survival
                       ggtheme = theme_classic(),
                       font.tickslab=c(20, "bold", "black"),
                       font.legend = c(20, "bold", "black"),
                       risk.table.fontsize = 8,
                       risk.table.y.text = FALSE,
                       risk.table.legend = FALSE,
                       legend.title = "BRAF-RAS Category",
                       legend.labs = c("BRAF-like", "RAS-like"),
                       title = "BRAF-RAS Score DFS",
                       tables.height = .33,
                       palette = c("red", "#2E9FDF","#86AA00"))+
  labs(y = "Disease Free Survival", x = "Months")  
DFS_Plot$plot <- DFS_Plot$plot + theme(panel.border = element_rect(colour = "black", size = 2, fill = NA))
DFS_Plot$table <- DFS_Plot$table + 
  theme(panel.border = element_rect(colour = "black", size = 2, fill = NA),
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
ggsave("outputs/22-0617_TCGA_PoorOutcome_Score_Survival/22-0617_TCGA_BRAF-RAS_DFS.png", 
       plot = plotp, 
       width = 7, height = 7,
       dpi = 600)





### BRAF-RAS Overall Survival Plot
#Survival curve calculation
fit<-survfit(Surv(OS_MONTHS, OS_Status_Numeric) ~ BRAFV600E_RAS, data=OS_BRAF_RAS_Score)
print(fit)

# Most basic plot
ggsurvplot(fit, data=OS_BRAF_RAS_Score)

# More involved, customized plot
OS_Plot <- ggsurvplot(fit,
                      data=OS_BRAF_RAS_Score,
                      pval = TRUE, pval.method = TRUE, conf.int = FALSE,
                      risk.table = TRUE, 
                      risk.table.col = "strata", 
                      legend = c(.75,.50),
                      pval.coord = c(.15,.45),
                      pval.size = 8,
                      break.time.by = 24,
                      pval.method.coord = c(.15,.60),
                      font.x = c(23, "bold", "black"),
                      font.y = c(23, "bold", "black"),
                      font.title = c(25, "bold", "black"),
                      # linetype = "strata", # Change line type by groups
                      # surv.median.line = "hv", # Specify median survival
                      ggtheme = theme_bw(),
                      font.tickslab=c(15, "bold", "black"),
                      font.legend = c(20, "bold", "black"),
                      risk.table.fontsize = 8,
                      risk.table.y.text = FALSE,
                      risk.table.legend = FALSE,
                      legend.title = "BRAF-RAS Category",
                      legend.labs = c("BRAF-like", "RAS-like"),
                      title = "BRAF-RAS Score DFS",
                      tables.height = .33,
                      palette = c("red", "#2E9FDF","#86AA00"))+
  labs(y = "Overall Survival", x = "Months")  
OS_Plot$table <- DFS_Plot$table + 
  theme(plot.title = element_text(size = 23, color = "black", face = "bold"),
        axis.title = element_text(size = 23, color = "black", face = "bold"),
        axis.text = element_text(size = 15, color = "black", face = "bold"),
        legend.position = "none")+
  labs(y = NULL) 
OS_Plot

# Saving the plots with ggsave
p1 = OS_Plot$plot
p2 = OS_Plot$table
plotp = cowplot::plot_grid(p1,p2,align = "v",axis = "b", ncol =1, rel_heights = c(3,1))
ggsave("outputs/22-0617_TCGA_PoorOutcome_Score_Survival/22-0617_TCGA_BRAF-RAS_OS.png", 
       plot = plotp, 
       width = 7, height = 7,
       dpi = 600)

