# Matthew A. Loberg
# May 3, 2022
# Purpose: 
# Exploration of the 526 gene BRAF-Poor Outcome score that I generated on 22-0426

# 22-0913_Update
# Using new aggressive tumor score that I calculated on 22-0908 (22-0908_Aggressive_Tumor_Score)

# 22-0927 Update
# Making a new color for aggressive lesions - a brighter pink coler (per Vivian/George)

# Load required packages
library(tidyverse)

# Read in data
ClinicalData <- read_csv(file = "data_in_use/VUMC.cohort.GX_9-8-22_V2.csv") # Latest version of the clinical data file

# Restrict to only sample types that were used in score calculation. Remove the following: 
ClinicalData <- ClinicalData %>% subset(Diagnosis != "normal" &
                                        Diagnosis != "FA" &
                                        Diagnosis != "OA" &
                                        Diagnosis != "MNG" &
                                        Diagnosis != "HT")

MAP <- read_csv(file = "data_in_use/22-0908_BRAF_Aggressive_Overlap_Genes_Score_NoFAOA_526.csv")

ClinicalData <- merge(ClinicalData, MAP)

# Merge PTCs and IFVPTCs, NIFTPs and EFVPTCs; also merge FA/OA and FTC/OTC
ClinicalData$Diagnosis_Merged <- ClinicalData$Diagnosis.with.iFVPTC.and.eFVPTC
for(i in 1:nrow(ClinicalData)){
  if(ClinicalData$Diagnosis_Merged[i] == "PTC" | ClinicalData$Diagnosis_Merged[i] == "IFVPTC"){
    ClinicalData$Diagnosis_Merged[i] <- "IFVPTC+\nPTC"
  }
  else if(ClinicalData$Diagnosis_Merged[i] == "NIFTP" | ClinicalData$Diagnosis_Merged[i] == "EFVPTC"){
    ClinicalData$Diagnosis_Merged[i] <- "NIFTP+\nEFVPTC"
  }
}


# Make plot for Fig3D with poor outcome labeled and restricted to primary and local disease
ClinicalData_Local <- ClinicalData %>% subset(Location.type == "Primary" |
                                              Location.type == "Localdisease")
plot <- ggplot(ClinicalData_Local, aes(Diagnosis_Merged, BRAFPoorOutcome_526)) +
  geom_boxplot(outlier.size = -1, 
               aes(),
               alpha = 0.2, 
               show.legend = FALSE) + 
  geom_jitter(aes(color = Aggression.category.for.Matt.Deseq),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5, 
              alpha = 0.9,
              show.legend = TRUE) +
  labs (x = "Diagnosis", y = "MAP Score") + 
  scale_color_manual(values = c("black", "magenta", "grey"), name = "Aggressive/Indolent", labels = c("Indolent", "Aggressive", "Unknown")) +
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 25), 
    axis.title.y = element_text(face = "bold", size = 22),
    axis.text = element_text(face = "bold", size = 16),
    legend.text = element_text(face = "bold", size = 15),
    legend.title = element_text(face = "bold", size = 20),
    legend.position = c(.25,.82),
    legend.box.background = element_rect(size = 2)) +
  scale_x_discrete(name ="Diagnosis", limits = c("FTC", "OTC", "NIFTP+\nEFVPTC", "IFVPTC+\nPTC", "PDTC", "ATC")) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red")
ggsave("outputs/22-0927_Fig3D_Aggressive_Tumor_Score_Exploration/22-0927_MAP_Local_by_Histotype_Aggressive_Coloring.png",
       width = 7,
       height = 5,
       plot, dpi = 600)

# Make plot w/ larger sized dots
plot <- ggplot(ClinicalData_Local, aes(Diagnosis_Merged, BRAFPoorOutcome_526)) +
  geom_boxplot(outlier.size = -1, 
               aes(),
               alpha = 0.2, 
               show.legend = FALSE) + 
  geom_jitter(aes(color = Aggression.category.for.Matt.Deseq),
              position = position_jitter(width = 0.1, height = 0),
              size = 2, 
              alpha = 0.9,
              show.legend = TRUE) +
  labs (x = "Diagnosis", y = "MAP Score") + 
  scale_color_manual(values = c("black", "magenta", "grey"), name = "Aggressive/Indolent", labels = c("Indolent", "Aggressive", "Unknown")) +
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 25), 
    axis.title.y = element_text(face = "bold", size = 22),
    axis.text = element_text(face = "bold", size = 16),
    legend.text = element_text(face = "bold", size = 15),
    legend.title = element_text(face = "bold", size = 20),
    legend.position = c(.25,.82),
    legend.box.background = element_rect(size = 2)) +
  scale_x_discrete(name ="Diagnosis", limits = c("FTC", "HC", "NIFTP+\nEFVPTC", "IFVPTC+\nPTC", "PDTC", "ATC")) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red")
ggsave("outputs/22-0927_Fig3D_Aggressive_Tumor_Score_Exploration/22-0927_MAP_Local_by_Histotype_Aggressive_Coloring_size2.png",
       width = 7,
       height = 5,
       plot, dpi = 600)
