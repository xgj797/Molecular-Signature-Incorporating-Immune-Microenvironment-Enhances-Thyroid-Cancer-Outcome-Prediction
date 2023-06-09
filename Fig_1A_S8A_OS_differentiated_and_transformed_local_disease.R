# OS_all_local.R

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

# Calculating length of OVERALL survival from initial therapy complete date
ClinicalData_localdisease_Malignant$Date.of.death <- as.Date(ClinicalData_localdisease_Malignant$Date.of.death, "%m/%d/%y")
ClinicalData_localdisease_Malignant$Initial.therapy.complete.date <- as.Date(ClinicalData_localdisease_Malignant$Initial.therapy.complete.date, "%m/%d/%y")
ClinicalData_localdisease_Malignant$Date.last.followup <- as.Date(ClinicalData_localdisease_Malignant$Date.last.followup, "%m/%d/%y")
ClinicalData_localdisease_Malignant$overallsurvival_Length <- NA
ClinicalData_localdisease_Malignant$ProgressionStatus <- 0
for(i in 1:nrow(ClinicalData_localdisease_Malignant)){
  if(!is.na(ClinicalData_localdisease_Malignant$Date.of.death[i])){
    ClinicalData_localdisease_Malignant$overallsurvival_Length[i] <- ClinicalData_localdisease_Malignant$Date.of.death[i] - 
      ClinicalData_localdisease_Malignant$Initial.therapy.complete.date[i]
    ClinicalData_localdisease_Malignant$ProgressionStatus[i] <- 1
  }
  else{
    ClinicalData_localdisease_Malignant$overallsurvival_Length[i] <- ClinicalData_localdisease_Malignant$Date.last.followup[i] - 
      ClinicalData_localdisease_Malignant$Initial.therapy.complete.date[i]
  }
}

ClinicalData_localdisease_Malignant$overallsurvival_Length <- ClinicalData_localdisease_Malignant$overallsurvival_Length/(365/12)

# Subset out samples if NA for overallsurvival_length
ClinicalData_localdisease_Malignant <- ClinicalData_localdisease_Malignant  %>% subset(!is.na(ClinicalData_localdisease_Malignant$overallsurvival_Length))

# Look at OS in known mutations
OS_TERT.mut_All_localdiseaseMalignant <- ClinicalData_localdisease_Malignant %>% subset(!is.na(TERT.mutation))
OS_TP53.mut_All_localdiseaseMalignant <- ClinicalData_localdisease_Malignant %>% subset(!is.na(TP53.mutation))
OS_PIK3CA.mut_All_localdiseaseMalignant <- ClinicalData_localdisease_Malignant %>% subset(!is.na(PIK3CA.mutation))

# Look at OS in MAP score (high vs low)
OS_MAPscore_All_localdiseaseMalignant <- ClinicalData_localdisease_Malignant %>% subset(!is.na(`MAP.score`))

# OS Plots
OS_TERT.mut_All_localdiseaseMalignant <- mutate(OS_TERT.mut_All_localdiseaseMalignant, Cat=ifelse(OS_TERT.mut_All_localdiseaseMalignant$TERT.mutation == 1, "Mutant", "Not Mutant"))
OS_TP53.mut_All_localdiseaseMalignant <- mutate(OS_TP53.mut_All_localdiseaseMalignant, Cat=ifelse(OS_TP53.mut_All_localdiseaseMalignant$TP53.mutation == 1, "Mutant", "Not Mutant"))
OS_PIK3CA.mut_All_localdiseaseMalignant <- mutate(OS_PIK3CA.mut_All_localdiseaseMalignant, Cat=ifelse(OS_PIK3CA.mut_All_localdiseaseMalignant$PIK3CA.mutation == 1, "Mutant", "Not Mutant"))
OS_MAPscore_All_localdiseaseMalignant <- mutate(OS_MAPscore_All_localdiseaseMalignant, Cat=ifelse(OS_MAPscore_All_localdiseaseMalignant$`MAP.score` >= 0, "High", "Low"))


### TERTp
#Survival curve calculation
fit<-survfit(Surv(overallsurvival_Length, ProgressionStatus) ~ Cat, data= OS_TERT.mut_All_localdiseaseMalignant)
print(fit)

# Most basic plot
ggsurvplot(fit, data=OS_TERT.mut_All_localdiseaseMalignant)

# More involved, customized plot
OS_Plot <- ggsurvplot(fit,
                       data=OS_TERT.mut_All_localdiseaseMalignant,
                       pval = TRUE, pval.method = TRUE, conf.int = FALSE,
                       risk.table = TRUE, 
                       risk.table.col = "strata", 
                       legend = c(.80,.20), # Move the legend
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
                       legend.title = "Mutation",
                       title = "TERT promoter OS",
                       legend.labs = c("YES", "NO"),
                       tables.height = .33,
                       break.time.by = 24,
                       palette = c("red", "#2E9FDF", "#86AA00"))+ # Can switch plot color palette order
  labs(y = "Progression Free Survival", x = "Months")  
OS_Plot$plot <- OS_Plot$plot + theme(panel.border = element_rect(colour = "black", size = 2, fill = NA))
OS_Plot$table <- OS_Plot$table + 
  theme(plot.title = element_text(size = 23, color = "black", face = "bold"),
        axis.title = element_text(size = 23, color = "black", face = "bold"),
        axis.text = element_text(size = 20, color = "black", face = "bold"),
        legend.position = "none")+
  theme(plot.margin = margin(0.1,1,0.1,0.1, "cm"))+ #gives space for last x-axis tick label, avoiding cutoff
  labs(y = NULL) 
OS_Plot

# Saving the plots with ggsave
p1 = OS_Plot$plot
p2 = OS_Plot$table
plotp = cowplot::plot_grid(p1,p2,align = "v",axis = "b", ncol =1, rel_heights = c(3,1))
ggsave("outputs/all_localdisease/OS/22-0909_All_Malignant_TERTp_OS.png", 
       plot = plotp, 
       width = 7, height = 7,
       dpi = 600)


### TP53
#Survival curve calculation
fit<-survfit(Surv(overallsurvival_Length, ProgressionStatus) ~ Cat, data= OS_TP53.mut_All_localdiseaseMalignant)
print(fit)

# Most basic plot
ggsurvplot(fit, data=OS_TP53.mut_All_localdiseaseMalignant)

# More involved, customized plot
OS_Plot <- ggsurvplot(fit,
                       data=OS_TP53.mut_All_localdiseaseMalignant,
                       pval = TRUE, pval.method = TRUE, conf.int = FALSE,
                       risk.table = TRUE, 
                       risk.table.col = "strata", 
                       legend = c(.85,.20), # Move the legend
                       pval.coord = c(55,.40), # Move the p-value
                       pval.size = 8,
                       pval.method.coord = c(55,.30), # Move the p-value method
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
                       legend.title = "Mutation",
                       title = "TP53 OS",
                       legend.labs = c("YES", "NO"),
                       tables.height = .33,
                       break.time.by = 24,
                       palette = c("red", "#2E9FDF", "#86AA00"))+ # Can switch plot color palette order
  labs(y = "Progression Free Survival", x = "Months")  
OS_Plot$plot <- OS_Plot$plot + theme(panel.border = element_rect(colour = "black", size = 2, fill = NA))
OS_Plot$table <- OS_Plot$table + 
  theme(plot.title = element_text(size = 23, color = "black", face = "bold"),
        axis.title = element_text(size = 23, color = "black", face = "bold"),
        axis.text = element_text(size = 20, color = "black", face = "bold"),
        legend.position = "none")+
  theme(plot.margin = margin(0.1,1,0.1,0.1, "cm"))+ #gives space for last x-axis tick label, avoiding cutoff
  labs(y = NULL) 
OS_Plot

# Saving the plots with ggsave
p1 = OS_Plot$plot
p2 = OS_Plot$table
plotp = cowplot::plot_grid(p1,p2,align = "v",axis = "b", ncol =1, rel_heights = c(3,1))
ggsave("outputs/all_localdisease/OS/22-0909_All_Malignant_TP53_OS.png", 
       plot = plotp, 
       width = 7, height = 7,
       dpi = 600)


### PIK3CA
#Survival curve calculation
fit<-survfit(Surv(overallsurvival_Length, ProgressionStatus) ~ Cat, data= OS_PIK3CA.mut_All_localdiseaseMalignant)
print(fit)

# Most basic plot
ggsurvplot(fit, data=OS_PIK3CA.mut_All_localdiseaseMalignant)

# More involved, customized plot
OS_Plot <- ggsurvplot(fit,
                       data=OS_PIK3CA.mut_All_localdiseaseMalignant,
                       pval = TRUE, pval.method = TRUE, conf.int = FALSE,
                       risk.table = TRUE, 
                       risk.table.col = "strata", 
                       legend = c(.85,.20), # Move the legand
                       pval.coord = c(60,.40), # Move the p-value
                       pval.size = 8,
                       pval.method.coord = c(60,.30),  # Move the p-value method
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
                       legend.title = "Mutation",
                       title = "PIK3CA OS",
                       legend.labs = c("YES", "NO"),
                       tables.height = .33,
                       break.time.by = 24,
                       palette = c("red", "#2E9FDF", "#86AA00"))+ # Can switch plot color palette order
  labs(y = "Progression Free Survival", x = "Months")  
OS_Plot$plot <- OS_Plot$plot + theme(panel.border = element_rect(colour = "black", size = 2, fill = NA))
OS_Plot$table <- OS_Plot$table + 
  theme(plot.title = element_text(size = 23, color = "black", face = "bold"),
        axis.title = element_text(size = 23, color = "black", face = "bold"),
        axis.text = element_text(size = 20, color = "black", face = "bold"),
        legend.position = "none")+
  theme(plot.margin = margin(0.1,1,0.1,0.1, "cm"))+ #gives space for last x-axis tick label, avoiding cutoff
  labs(y = NULL) 
OS_Plot

# Saving the plots with ggsave
p1 = OS_Plot$plot
p2 = OS_Plot$table
plotp = cowplot::plot_grid(p1,p2,align = "v",axis = "b", ncol =1, rel_heights = c(3,1))
ggsave("outputs/all_localdisease/OS/22-0909_All_Malignant_PIK3CA_OS.png", 
       plot = plotp, 
       width = 7, height = 7,
       dpi = 600)


### MAP score high/low
#Survival curve calculation
fit<-survfit(Surv(overallsurvival_Length, ProgressionStatus) ~ Cat, data= OS_MAPscore_All_localdiseaseMalignant)
print(fit)

# Most basic plot
ggsurvplot(fit, data=OS_MAPscore_All_localdiseaseMalignant)

# More involved, customized plot
OS_Plot <- ggsurvplot(fit,
                       data=OS_MAPscore_All_localdiseaseMalignant,
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
                       title = "MAP score OS",
                       legend.labs = c("Positive", "Negative"),
                       tables.height = .33,
                       break.time.by = 24,
                       palette = c("#df2ead", "#3e00aa", "#86AA00"))+ # Can switch plot color palette order
  labs(y = "Progression Free Survival", x = "Months")  
OS_Plot$plot <- OS_Plot$plot + theme(panel.border = element_rect(colour = "black", size = 2, fill = NA))
OS_Plot$table <- OS_Plot$table + 
  theme(plot.title = element_text(size = 23, color = "black", face = "bold"),
        axis.title = element_text(size = 23, color = "black", face = "bold"),
        axis.text = element_text(size = 20, color = "black", face = "bold"),
        legend.position = "none")+
  theme(plot.margin = margin(0.1,1,0.1,0.1, "cm"))+ #gives space for last x-axis tick label, avoiding cutoff
  labs(y = NULL) 
OS_Plot

# Saving the plots with ggsave
p1 = OS_Plot$plot
p2 = OS_Plot$table
plotp = cowplot::plot_grid(p1,p2,align = "v",axis = "b", ncol =1, rel_heights = c(3,1))
ggsave("outputs/all_localdisease/OS/22-0909_All_Malignant_Aggression_OS.png", 
       plot = plotp, 
       width = 7, height = 7,
       dpi = 600)
