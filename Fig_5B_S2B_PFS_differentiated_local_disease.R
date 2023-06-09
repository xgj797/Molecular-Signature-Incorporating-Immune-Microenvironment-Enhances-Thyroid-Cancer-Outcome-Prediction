# PFS_all_welldiff.R

# Load required packages
library(ggplot2)
library(survival)
library(survminer)
library(tidyverse)
library(cowplot)

# Set your working directory
setwd("_")

# Load Data
ClinicalData <- read_csv("VUMC.cohort.GX_9-15-22_for-Fei.csv") # All Data
ClinicalData <- as_tibble(ClinicalData)
# File from Fei

# Restrict data to the rows needed
# Start by restricting to local disease
ClinicalData_localdisease <- ClinicalData %>% subset(Location.type == "Primary" | Location.type == "Localdisease")
ClinicalData_localdisease_Malignant <- ClinicalData_localdisease %>% subset(is.na(Progressive.disease) |
                                                                  Progressive.disease == "0" |
                                                                  Progressive.disease == "1") #This excludes non-malignant samples which only have non-0 or -1 values in this column

# If restricting to well-diff, exclude ATC and PDTC 
ClinicalData_welldiff_Malignant <- ClinicalData_localdisease_Malignant %>% subset(Diagnosis != "ATC" & Diagnosis != "PDTC")

# Calculating length of PROGRESSION FREE survival from initial therapy complete date
ClinicalData_welldiff_Malignant$Recur.or.progress.date <- as.Date(ClinicalData_welldiff_Malignant$Recur.or.progress.date, "%m/%d/%y")
ClinicalData_welldiff_Malignant$Initial.therapy.complete.date <- as.Date(ClinicalData_welldiff_Malignant$Initial.therapy.complete.date, "%m/%d/%y")
ClinicalData_welldiff_Malignant$Date.last.followup <- as.Date(ClinicalData_welldiff_Malignant$Date.last.followup, "%m/%d/%y")
ClinicalData_welldiff_Malignant$ProgressionFree_Length <- NA
ClinicalData_welldiff_Malignant$ProgressionStatus <- 0
for(i in 1:nrow(ClinicalData_welldiff_Malignant)){
  if(!is.na(ClinicalData_welldiff_Malignant$Recur.or.progress.date[i])){
    ClinicalData_welldiff_Malignant$ProgressionFree_Length[i] <- ClinicalData_welldiff_Malignant$Recur.or.progress.date[i] - 
      ClinicalData_welldiff_Malignant$Initial.therapy.complete.date[i]
    ClinicalData_welldiff_Malignant$ProgressionStatus[i] <- 1
  }
  else{
    ClinicalData_welldiff_Malignant$ProgressionFree_Length[i] <- ClinicalData_welldiff_Malignant$Date.last.followup[i] - 
      ClinicalData_welldiff_Malignant$Initial.therapy.complete.date[i]
  }
}

ClinicalData_welldiff_Malignant$ProgressionFree_Length <- ClinicalData_welldiff_Malignant$ProgressionFree_Length/(365/12)

# Subset out samples if NA for progressionfree_length
ClinicalData_welldiff_Malignant <- ClinicalData_welldiff_Malignant  %>% subset(!is.na(ClinicalData_welldiff_Malignant$ProgressionFree_Length))

# Look at PFS in MAP score (high vs low)
PFS_MAPscore_welldiff_Malignant <- ClinicalData_welldiff_Malignant %>% subset(!is.na(`MAP.score`))

# Look at PFS in BRS
PFS_BRS_All_welldiff_Malignant <- ClinicalData_welldiff_Malignant %>% subset(!is.na(BRS))

# PFS Plots
PFS_MAPscore_welldiff_Malignant <- mutate(PFS_MAPscore_welldiff_Malignant, Cat=ifelse(PFS_MAPscore_welldiff_Malignant$`MAP.score` >= 0, "High", "Low"))
PFS_BRS_All_welldiff_Malignant <- mutate(PFS_BRS_All_welldiff_Malignant, Cat=ifelse(PFS_BRS_All_welldiff_Malignant$BRS > 0, "RAS-Like", "BRAF-Like"))


### MAP score high/low
#Survival curve calculation
fit<-survfit(Surv(ProgressionFree_Length, ProgressionStatus) ~ Cat, data= PFS_MAPscore_welldiff_Malignant)
print(fit)

# Most basic plot
ggsurvplot(fit, data=PFS_MAPscore_welldiff_Malignant)

# More involved, customized plot
PFS_Plot <- ggsurvplot(fit,
                       data=PFS_MAPscore_welldiff_Malignant,
                       pval = TRUE, pval.method = TRUE, conf.int = FALSE,
                       risk.table = TRUE, 
                       risk.table.col = "strata", 
                       legend = c(.75,.20), # Move the legand
                       pval.coord = c(.10,.20), # Move the p-value
                       pval.size = 8,
                       pval.method.coord = c(.10,.10),  # Move the p-value method
                       font.x = c(23, "bold", "black"),
                       font.y = c(23, "bold", "black"),
                       font.title = c(25, "bold", "black"),
                       # linetype = "strata", # Change line type by groups
                       # surv.median.line = "hv", # Specify median survival
                       ggtheme = theme_classic(),
                       font.tickslab=c(17, "bold", "black"),
                       font.legend = c(20, "bold", "black"),
                       risk.table.fontsize = 7,
                       risk.table.y.text = FALSE,
                       risk.table.legend = FALSE,
                       legend.title = "MAP score",
                       title = "MAP score PFS",
                       legend.labs = c("Positive", "Negative"),
                       tables.height = .33,
                       break.time.by = 24,
                       palette = c("#df2ead", "#3e00aa", "#86AA00"))+ # Can switch plot color palette order
  labs(y = "Progression Free Survival", x = "Months")  
PFS_Plot$plot <- PFS_Plot$plot + theme(panel.border = element_rect(colour = "black", size = 2, fill = NA))
PFS_Plot$table <- PFS_Plot$table + 
  theme(plot.title = element_text(size = 23, color = "black", face = "bold"),
        axis.title = element_text(size = 23, color = "black", face = "bold"),
        axis.text = element_text(size = 20, color = "black", face = "bold"),
        legend.position = "none")+
  theme(plot.margin = margin(0.1,1,0.1,0.1, "cm"))+ #gives space for last x-axis tick label, avoiding cutoff
  labs(y = NULL) 
PFS_Plot

# Saving the plots with ggsave
p1 = PFS_Plot$plot
p2 = PFS_Plot$table
plotp = cowplot::plot_grid(p1,p2,align = "v",axis = "b", ncol =1, rel_heights = c(3,1))
ggsave("outputs/welldiff/PFS/22-0909_All_Well_Diff_Malignant_Aggression_PFS.png", 
       plot = plotp, 
       width = 7, height = 7,
       dpi = 600)


### BRS PFS
#Survival curve calculation
fit<-survfit(Surv(ProgressionFree_Length, ProgressionStatus) ~ Cat, data= PFS_BRS_All_welldiff_Malignant)
print(fit)

# Most basic plot
ggsurvplot(fit, data=PFS_BRS_All_welldiff_Malignant)

# More involved, customized plot
PFS_Plot <- ggsurvplot(fit,
                       data=PFS_BRS_All_welldiff_Malignant,
                       pval = TRUE, pval.method = TRUE, conf.int = FALSE,
                       risk.table = TRUE, 
                       risk.table.col = "strata", 
                       legend = c(.70,.20), # Move the legend
                       pval.coord = c(.15,.20), # Move the p-value
                       pval.size = 8,
                       pval.method.coord = c(.15,.10), # Move the p-value method
                       font.x = c(23, "bold", "black"),
                       font.y = c(23, "bold", "black"),
                       font.title = c(25, "bold", "black"),
                       # linetype = "strata", # Change line type by groups
                       # surv.median.line = "hv", # Specify median survival
                       ggtheme = theme_classic(),
                       font.tickslab=c(17, "bold", "black"),
                       font.legend = c(20, "bold", "black"),
                       risk.table.fontsize = 7,
                       risk.table.y.text = FALSE,
                       risk.table.legend = FALSE,
                       legend.title = "BRAF or RAS-like",
                       title = "BRAF-RAS Score PFS",
                       legend.labs = c("BRAF-Like", "RAS-Like"),
                       tables.height = .33,
                       break.time.by = 24,
                       palette = c("red", "#2E9FDF", "#86AA00"))+ # Can switch plot color palette order
  labs(y = "Progression Free Survival", x = "Months")  
PFS_Plot$plot <- PFS_Plot$plot + theme(panel.border = element_rect(colour = "black", size = 2, fill = NA))
PFS_Plot$table <- PFS_Plot$table + 
  theme(plot.title = element_text(size = 23, color = "black", face = "bold"),
        axis.title = element_text(size = 23, color = "black", face = "bold"),
        axis.text = element_text(size = 20, color = "black", face = "bold"),
        legend.position = "none")+
  theme(plot.margin = margin(0.1,1,0.1,0.1, "cm"))+ #gives space for last x-axis tick label, avoiding cutoff
  labs(y = NULL) 
PFS_Plot

# Saving the plots with ggsave
p1 = PFS_Plot$plot
p2 = PFS_Plot$table
plotp = cowplot::plot_grid(p1,p2,align = "v",axis = "b", ncol =1, rel_heights = c(3,1))
ggsave("outputs/welldiff/PFS/22-0909_All_Well_Diff_Malignant_BRS_PFS.png", 
       plot = plotp, 
       width = 7, height = 7,
       dpi = 600)
