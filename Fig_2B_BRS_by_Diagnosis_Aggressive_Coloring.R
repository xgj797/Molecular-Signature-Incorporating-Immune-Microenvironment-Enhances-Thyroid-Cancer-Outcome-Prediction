### Figure 2B: BRS by aggressive yes/no 

### Load packages
library(tidyverse)

# Read in clinical data and restrict to Diagnosis and location
# Read in clinical data
ClinicalData <- read_csv(file = "data_in_use/VUMC.cohort.GX_9-8-22.csv")
ClinicalData <- ClinicalData %>% subset(Diagnosis != "normal")
ClinicalData <- ClinicalData[c("IP", 
                               "Location.type", 
                               "Diagnosis", 
                               "Poor.outcome",
                               "Aggression.category.for.Matt.Deseq",
                               "Diagnosis.with.HTPTC",
                               "Diagnosis.with.iFVPTC.and.eFVPTC",
                               "RNA.ID", 
                               "BRS", "TDS", "ERK", "PI3K_AKT_MTOR", "WNT_Canon", "WNT_NonCanon")] # Restrict columns
ClinicalData$BRS <- as.numeric(ClinicalData$BRS)

ClinicalData <- ClinicalData %>% dplyr::rename(Sequencing.ID = RNA.ID) # dplyr::rename RNA.ID column to Sequencing.ID

# Function for Plotting BRS by Histotypes (no MNG, HT) with a red line at 0.0 (denoting switch from BRAF-like to RAS-like)
Plot_BRS_Histotypes_PoorOutcome <- function(DataSet, XVar, YVar, YAxisTitle){
  plot <- ggplot(DataSet, aes(XVar, YVar)) +
    geom_boxplot(outlier.size = -1, 
                 aes(),
                 alpha = 0.2, 
                 show.legend = FALSE) + 
    geom_jitter(aes(color = Aggression.category.for.Matt.Deseq),
                position = position_jitter(width = 0.1, height = 0),
                size = 2, 
                alpha = 0.9,
                show.legend = TRUE) +
    scale_color_manual(values = c("black", "magenta", "grey"), name = "Aggressive/Indolent", labels = c("Indolent", "Aggressive", "Unknown")) +
    ylim(-.8,.8) +
    labs (x = "Diagnosis", y = YAxisTitle, color = "Aggression.category.for.Matt.Deseq") + 
    theme_classic() + 
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5), 
      axis.title = element_text(face = "bold", size = 30),
      axis.text = element_text(face = "bold", size = 20),
      legend.position = c(.15,.2),
      legend.box.background = element_rect(size = 2),
      legend.text = element_text(face = "bold", size = 15),
      legend.title = element_text(face = "bold", size = 15)) +
    scale_x_discrete(name ="Diagnosis", limits = c("FA\n+\nOA","FTC\n+\nOTC", "EFVPTC\n+\nNIFTP", "PTC\n+\nIFVPTC", "PDTC", "ATC")) +
    geom_hline(yintercept = 0.0, linetype = 2, colour = "red")
  #geom_hline(yintercept = 0.5, linetype = 2, colour = "red") +# Adding line at 0.5 to see changes in score
  #geom_hline(yintercept = -0.5, linetype = 2, colour = "red")# Adding line at 0.5 to see changes in score
  return(plot)
}

### Data formatting for plots
ClinicalData$Diagnosis_Simplified <- ClinicalData$Diagnosis.with.iFVPTC.and.eFVPTC
for(i in 1:nrow(ClinicalData)){
  if(ClinicalData$Diagnosis[i] == "OA" | ClinicalData$Diagnosis[i] == "FA"){
    ClinicalData$Diagnosis_Simplified[i] <- "FA\n+\nOA"
  }
  else if(ClinicalData$Diagnosis[i] == "OTC" | ClinicalData$Diagnosis[i] == "FTC"){
    ClinicalData$Diagnosis_Simplified[i] <- "FTC\n+\nOTC"
  }
  else if(ClinicalData$Diagnosis.with.iFVPTC.and.eFVPTC[i] == "PTC" | ClinicalData$Diagnosis.with.iFVPTC.and.eFVPTC[i] == "IFVPTC"){
    ClinicalData$Diagnosis_Simplified[i] <- "PTC\n+\nIFVPTC"
  }
  else if(ClinicalData$Diagnosis.with.iFVPTC.and.eFVPTC[i] == "NIFTP" | ClinicalData$Diagnosis.with.iFVPTC.and.eFVPTC[i] == "EFVPTC"){
    ClinicalData$Diagnosis_Simplified[i] <- "EFVPTC\n+\nNIFTP"
  }
}

ClinicalData <- ClinicalData %>% subset(!is.na(BRS)) # remove samples without BRS (no RNA data)
ClinicalData_Dx_Restricted <- ClinicalData %>% subset(Diagnosis != "MNG" & Diagnosis != "HT") # Restrict diagnoses
ClinicalData_NonMet <- ClinicalData_Dx_Restricted %>% subset(Location.type == "NonMet")

# Plot Internal BRS for all histotypes (NonMet only)
# Colored by Aggressive/Indolent

Internal_BRS_by_Histotype_NonMet <- Plot_BRS_Histotypes_PoorOutcome(ClinicalData_NonMet,
                                                                   ClinicalData_NonMet$Diagnosis_Simplified,
                                                                   ClinicalData_NonMet$BRS,
                                                                   "BRAF-RAS Score")
ggsave("outputs/22-0927_BRS_Fig2B/22-0927_BRS_Fig2B_NonMet.png",
       width = 9, height = 5, 
       Internal_BRS_by_Histotype_NonMet, dpi = 600)
