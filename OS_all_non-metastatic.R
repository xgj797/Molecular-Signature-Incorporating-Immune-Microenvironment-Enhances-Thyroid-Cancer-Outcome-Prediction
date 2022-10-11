# OS_all_non-metastatic.R
#using 50th percentile CAF high and low based on OVERALL cohort

# Load required packages
library(ggplot2)
library(survival)
library(survminer)
library(tidyverse)
library(cowplot)

# Set your working directory
setwd("/Users/georgexu/Dropbox/Vanderbilt\ Grad\ School/Lab/Projects/Tumor\ aggression\ project\ with\ Fei/PFS_Code_9-9-22")


# Load Data
ClinicalData <- read_csv("VUMC.cohort.GX_9-15-22_for-Fei.csv") # All Data
ClinicalData <- as_tibble(ClinicalData)
# File from Fei

# Restrict data to the rows needed
# Start by restricting to non-metastatic
ClinicalData_Nonmet <- ClinicalData %>% subset(Location.type == "Primary" | Location.type == "Localdisease")
ClinicalData_Nonmet_Malignant <- ClinicalData_Nonmet %>% subset(is.na(Progressive.disease) |
                                                              Progressive.disease == "0" |
                                                              Progressive.disease == "1") #This excludes non-malignant samples which only have non-0 or -1 values in this column

## Excluding 170 and 173 because sample taken long after initial diagnosis
ClinicalData_Nonmet_Malignant <- ClinicalData_Nonmet_Malignant %>% subset(IP != "IP170" & IP != "IP173")

## Excluding IP 272, and IP 325 -> Local Disease that we dropped
ClinicalData_Nonmet_Malignant <- ClinicalData_Nonmet_Malignant %>% subset(IP != "IP215" & IP != "IP272" & IP != "IP325")

# Calculating length of OVERALL survival from initial therapy complete date
ClinicalData_Nonmet_Malignant$Date.of.death <- as.Date(ClinicalData_Nonmet_Malignant$Date.of.death, "%m/%d/%y")
ClinicalData_Nonmet_Malignant$Initial.therapy.complete.date <- as.Date(ClinicalData_Nonmet_Malignant$Initial.therapy.complete.date, "%m/%d/%y")
ClinicalData_Nonmet_Malignant$Date.last.followup <- as.Date(ClinicalData_Nonmet_Malignant$Date.last.followup, "%m/%d/%y")
ClinicalData_Nonmet_Malignant$ProgressionFree_Length <- NA
ClinicalData_Nonmet_Malignant$ProgressionStatus <- 0
for(i in 1:nrow(ClinicalData_Nonmet_Malignant)){
  if(!is.na(ClinicalData_Nonmet_Malignant$Date.of.death[i])){
    ClinicalData_Nonmet_Malignant$ProgressionFree_Length[i] <- ClinicalData_Nonmet_Malignant$Date.of.death[i] - 
      ClinicalData_Nonmet_Malignant$Initial.therapy.complete.date[i]
    ClinicalData_Nonmet_Malignant$ProgressionStatus[i] <- 1
  }
  else{
    ClinicalData_Nonmet_Malignant$ProgressionFree_Length[i] <- ClinicalData_Nonmet_Malignant$Date.last.followup[i] - 
      ClinicalData_Nonmet_Malignant$Initial.therapy.complete.date[i]
  }
}

ClinicalData_Nonmet_Malignant$ProgressionFree_Length <- ClinicalData_Nonmet_Malignant$ProgressionFree_Length/(365/12)


# Subset out samples if NA for progressionfree_length
ClinicalData_Nonmet_Malignant <- ClinicalData_Nonmet_Malignant  %>% subset(!is.na(ClinicalData_Nonmet_Malignant$ProgressionFree_Length))


# Look at OS in known mutations
OS_TERT.TP53.PIK3CA.mut_All_NonmetMalignant <- ClinicalData_Nonmet_Malignant %>% subset(!is.na(TERT.TP53.PIK3CA.mut))
OS_TERT.mut_All_NonmetMalignant <- ClinicalData_Nonmet_Malignant %>% subset(!is.na(TERT.mutation))
OS_TP53.mut_All_NonmetMalignant <- ClinicalData_Nonmet_Malignant %>% subset(!is.na(TP53.mutation))
OS_PIK3CA.mut_All_NonmetMalignant <- ClinicalData_Nonmet_Malignant %>% subset(!is.na(PIK3CA.mutation))

# Look at OS in BRS
OS_BRS_All_NonmetMalignant <- ClinicalData_Nonmet_Malignant %>% subset(!is.na(BRS))

# Look at OS in CAFs (high vs low)
OS_CAF_EPIC_All_NonmetMalignant <- ClinicalData_Nonmet_Malignant %>% subset(!is.na(`Cancer associated fibroblast_EPIC_high`))
OS_CAF_MCPCOUNTER_All_NonmetMalignant <- ClinicalData_Nonmet_Malignant %>% subset(!is.na(`Cancer associated fibroblast_MCPCOUNTER_high`))

# Look at PFS in Aggression score (high vs low)
OS_Aggression_All_NonmetMalignant <- ClinicalData_Nonmet_Malignant %>% subset(!is.na(`ATS`))

# Look at OS in Jennifer groups 
OS_Jennifer_All_NonmetMalignant <- ClinicalData_Nonmet_Malignant %>% subset(!is.na(Aggression.category.detailed))

# OS Plots
OS_BRS_All_NonmetMalignant <- mutate(OS_BRS_All_NonmetMalignant, Cat=ifelse(OS_BRS_All_NonmetMalignant$BRS > 0, "RAS-Like", "BRAF-Like"))
OS_TERT.TP53.PIK3CA.mut_All_NonmetMalignant <- mutate(OS_TERT.TP53.PIK3CA.mut_All_NonmetMalignant, Cat=ifelse(OS_TERT.TP53.PIK3CA.mut_All_NonmetMalignant$TERT.TP53.PIK3CA.mut == 1, "Mutant", "Not Mutant"))
OS_TERT.mut_All_NonmetMalignant <- mutate(OS_TERT.mut_All_NonmetMalignant, Cat=ifelse(OS_TERT.mut_All_NonmetMalignant$TERT.mutation == 1, "Mutant", "Not Mutant"))
OS_TP53.mut_All_NonmetMalignant <- mutate(OS_TP53.mut_All_NonmetMalignant, Cat=ifelse(OS_TP53.mut_All_NonmetMalignant$TP53.mutation == 1, "Mutant", "Not Mutant"))
OS_PIK3CA.mut_All_NonmetMalignant <- mutate(OS_PIK3CA.mut_All_NonmetMalignant, Cat=ifelse(OS_PIK3CA.mut_All_NonmetMalignant$PIK3CA.mutation == 1, "Mutant", "Not Mutant"))
OS_CAF_EPIC_All_NonmetMalignant <- mutate(OS_CAF_EPIC_All_NonmetMalignant, Cat=ifelse(OS_CAF_EPIC_All_NonmetMalignant$`Cancer associated fibroblast_EPIC_high` == 1, "High", "Low"))
OS_CAF_MCPCOUNTER_All_NonmetMalignant <- mutate(OS_CAF_MCPCOUNTER_All_NonmetMalignant, Cat=ifelse(OS_CAF_MCPCOUNTER_All_NonmetMalignant$`Cancer associated fibroblast_MCPCOUNTER_high` == 1, "High", "Low"))
OS_Aggression_All_NonmetMalignant <- mutate(OS_Aggression_All_NonmetMalignant, Cat=ifelse(OS_Aggression_All_NonmetMalignant$`ATS` >= 0, "High", "Low"))


### TERTp/TP53/PIK3CA
#Survival curve calculation
fit<-survfit(Surv(ProgressionFree_Length, ProgressionStatus) ~ Cat, data= OS_TERT.TP53.PIK3CA.mut_All_NonmetMalignant)
print(fit)

# Most basic plot
ggsurvplot(fit, data=OS_TERT.TP53.PIK3CA.mut_All_NonmetMalignant)

# More involved, customized plot
OS_Plot <- ggsurvplot(fit,
           data=OS_TERT.TP53.PIK3CA.mut_All_NonmetMalignant,
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
           title = "TP53/TERTp/PI3KCA OS",
           legend.labs = c("YES", "NO"),
           tables.height = .33,
           break.time.by = 24,
           palette = c("red", "#2E9FDF", "#86AA00"))+ # Can switch plot color palette order
    labs(y = "Overall Survival", x = "Months")  
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
ggsave("outputs/all_Nonmet/OS/22-0909_All_Malignant_TERTp_TP53_PIK3CA_OS.png", 
       plot = plotp, 
       width = 7, height = 7,
       dpi = 600)



### BRS OS
#Survival curve calculation
fit<-survfit(Surv(ProgressionFree_Length, ProgressionStatus) ~ Cat, data= OS_BRS_All_NonmetMalignant)
print(fit)

# Most basic plot
ggsurvplot(fit, data=OS_BRS_All_NonmetMalignant)

# More involved, customized plot
OS_Plot <- ggsurvplot(fit,
                       data=OS_BRS_All_NonmetMalignant,
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
                       title = "BRAF-RAS Score OS",
                       legend.labs = c("BRAF-Like", "RAS-Like"),
                       tables.height = .33,
                       break.time.by = 24,
                       palette = c("red", "#2E9FDF", "#86AA00"))+ # Can switch plot color palette order
  labs(y = "Overall Survival", x = "Months")  
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
ggsave("outputs/all_Nonmet/OS/22-0909_All_Malignant_BRS_OS.png", 
       plot = plotp, 
       width = 7, height = 7,
       dpi = 600)


### TERTp
#Survival curve calculation
fit<-survfit(Surv(ProgressionFree_Length, ProgressionStatus) ~ Cat, data= OS_TERT.mut_All_NonmetMalignant)
print(fit)

# Most basic plot
ggsurvplot(fit, data=OS_TERT.mut_All_NonmetMalignant)

# More involved, customized plot
OS_Plot <- ggsurvplot(fit,
                       data=OS_TERT.mut_All_NonmetMalignant,
                       pval = TRUE, pval.method = TRUE, conf.int = FALSE,
                       risk.table = TRUE, 
                       risk.table.col = "strata", 
                       legend = c(.50,.20), # Move the legend
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
  labs(y = "Overall Survival", x = "Months")  
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
ggsave("outputs/all_Nonmet/OS/22-0909_All_Malignant_TERTp_OS.png", 
       plot = plotp, 
       width = 7, height = 7,
       dpi = 600)


### TP53
#Survival curve calculation
fit<-survfit(Surv(ProgressionFree_Length, ProgressionStatus) ~ Cat, data= OS_TP53.mut_All_NonmetMalignant)
print(fit)

# Most basic plot
ggsurvplot(fit, data=OS_TP53.mut_All_NonmetMalignant)

# More involved, customized plot
OS_Plot <- ggsurvplot(fit,
                       data=OS_TP53.mut_All_NonmetMalignant,
                       pval = TRUE, pval.method = TRUE, conf.int = FALSE,
                       risk.table = TRUE, 
                       risk.table.col = "strata", 
                       legend = c(.50,.20), # Move the legend
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
                       title = "TP53 OS",
                       legend.labs = c("YES", "NO"),
                       tables.height = .33,
                       break.time.by = 24,
                       palette = c("red", "#2E9FDF", "#86AA00"))+ # Can switch plot color palette order
  labs(y = "Overall Survival", x = "Months")  
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
ggsave("outputs/all_Nonmet/OS/22-0909_All_Malignant_TP53_OS.png", 
       plot = plotp, 
       width = 7, height = 7,
       dpi = 600)


### PIK3CA
#Survival curve calculation
fit<-survfit(Surv(ProgressionFree_Length, ProgressionStatus) ~ Cat, data= OS_PIK3CA.mut_All_NonmetMalignant)
print(fit)

# Most basic plot
ggsurvplot(fit, data=OS_PIK3CA.mut_All_NonmetMalignant)

# More involved, customized plot
OS_Plot <- ggsurvplot(fit,
                       data=OS_PIK3CA.mut_All_NonmetMalignant,
                       pval = TRUE, pval.method = TRUE, conf.int = FALSE,
                       risk.table = TRUE, 
                       risk.table.col = "strata", 
                       legend = c(.50,.20), # Move the legend
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
                       title = "PIK3CA OS",
                       legend.labs = c("YES", "NO"),
                       tables.height = .33,
                       break.time.by = 24,
                       palette = c("red", "#2E9FDF", "#86AA00"))+ # Can switch plot color palette order
  labs(y = "Overall Survival", x = "Months")  
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
ggsave("outputs/all_Nonmet/OS/22-0909_All_Malignant_PIK3CA_OS.png", 
       plot = plotp, 
       width = 7, height = 7,
       dpi = 600)


### CAF EPIC high/low
#Survival curve calculation
fit<-survfit(Surv(ProgressionFree_Length, ProgressionStatus) ~ Cat, data= OS_CAF_EPIC_All_NonmetMalignant)
print(fit)

# Most basic plot
ggsurvplot(fit, data=OS_CAF_EPIC_All_NonmetMalignant)

# More involved, customized plot
OS_Plot <- ggsurvplot(fit,
                       data=OS_CAF_EPIC_All_NonmetMalignant,
                       pval = TRUE, pval.method = TRUE, conf.int = FALSE,
                       risk.table = TRUE, 
                       risk.table.col = "strata", 
                       legend = c(.85,.20), # Move the legand
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
                       legend.title = "CAF score",
                       title = "Cancer-Associated Fibroblast\nEPIC OS",
                       legend.labs = c("High", "Low"),
                       tables.height = .33,
                       break.time.by = 24,
                       palette = c("#02f02e", "#48827a", "#86AA00"))+ # Can switch plot color palette order
  labs(y = "Overall Survival", x = "Months")  
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
ggsave("outputs/all_Nonmet/OS/22-0909_All_Malignant_CAF_EPIC_OS.png", 
       plot = plotp, 
       width = 7, height = 7,
       dpi = 600)


### CAF MCPCOUNTER high/low
#Survival curve calculation
fit<-survfit(Surv(ProgressionFree_Length, ProgressionStatus) ~ Cat, data= OS_CAF_MCPCOUNTER_All_NonmetMalignant)
print(fit)

# Most basic plot
ggsurvplot(fit, data=OS_CAF_MCPCOUNTER_All_NonmetMalignant)

# More involved, customized plot
OS_Plot <- ggsurvplot(fit,
                       data=OS_CAF_MCPCOUNTER_All_NonmetMalignant,
                       pval = TRUE, pval.method = TRUE, conf.int = FALSE,
                       risk.table = TRUE, 
                       risk.table.col = "strata", 
                       legend = c(.85,.20), # Move the legand
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
                       legend.title = "CAF score",
                       title = "Cancer-Associated Fibroblast\nMCPCOUNTER OS",
                       legend.labs = c("High", "Low"),
                       tables.height = .33,
                       break.time.by = 24,
                       palette = c("#02f02e", "#48827a", "#86AA00"))+ # Can switch plot color palette order
  labs(y = "Overall Survival", x = "Months")  
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
ggsave("outputs/all_Nonmet/OS/22-0909_All_Malignant_CAF_MCPCOUNTER_OS.png", 
       plot = plotp, 
       width = 7, height = 7,
       dpi = 600)


### Aggression score high/low
#Survival curve calculation
fit<-survfit(Surv(ProgressionFree_Length, ProgressionStatus) ~ Cat, data= OS_Aggression_All_NonmetMalignant)
print(fit)

# Most basic plot
ggsurvplot(fit, data=OS_Aggression_All_NonmetMalignant)

# More involved, customized plot
OS_Plot <- ggsurvplot(fit,
                       data=OS_Aggression_All_NonmetMalignant,
                       pval = TRUE, pval.method = TRUE, conf.int = FALSE,
                       risk.table = TRUE, 
                       risk.table.col = "strata", 
                       legend = c(.55,.20), # Move the legand
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
                       legend.title = "ATS",
                       title = "Aggression Score OS",
                       legend.labs = c("Positive", "Negative"),
                       tables.height = .33,
                       break.time.by = 24,
                       palette = c("#df2ead", "#3e00aa", "#86AA00"))+ # Can switch plot color palette order
  labs(y = "Overall Survival", x = "Months")  
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
ggsave("outputs/all_Nonmet/OS/22-0909_All_Malignant_Aggression_OS.png", 
       plot = plotp, 
       width = 7, height = 7,
       dpi = 600)


### Jennifer categories
#Survival curve calculation
fit<- survfit(Surv(ProgressionFree_Length, ProgressionStatus) ~ Aggression.category.detailed, data = OS_Jennifer_All_NonmetMalignant)
print(fit)

# Most basic plot
ggsurvplot(fit, data=OS_Jennifer_All_NonmetMalignant)

# More involved, customized plot
OS_Plot <- ggsurvplot(fit,
                       data=OS_Jennifer_All_NonmetMalignant,
                       #pval = TRUE, pval.method = TRUE, conf.int = FALSE,
                       risk.table = TRUE, 
                       risk.table.col = "strata", 
                       legend = c(.3,.2), # Move the legand
                       #pval.coord = c(.10,.20), # Move the p-value
                       #pval.size = 8,
                       #pval.method.coord = c(.10,.10),  # Move the p-value method
                       font.x = c(23, "bold", "black"),
                       font.y = c(23, "bold", "black"),
                       font.title = c(25, "bold", "black"),
                       # linetype = "strata", # Change line type by groups
                       # surv.median.line = "hv", # Specify median survival
                       ggtheme = theme_classic(),
                       font.tickslab=c(17, "bold", "black"),
                       font.legend = c(7, "bold", "black"),
                       risk.table.fontsize = 5,
                       risk.table.y.text = FALSE,
                       risk.table.legend = FALSE,
                       legend.title = "Jennifer categories",
                       title = "Jennifer categories OS",
                       #legend.labs = c("High", "Low"),
                       tables.height = .33,
                       break.time.by = 24,
                       palette = c("blue", "pink", "red", "orange", "black"))+ # Can switch plot color palette order
  labs(y = "Overall Survival", x = "Months")  
OS_Plot$plot <- OS_Plot$plot + theme(panel.border = element_rect(colour = "black", size = 2, fill = NA))
OS_Plot$table <- OS_Plot$table + 
  theme(plot.title = element_text(size = 13, color = "black", face = "bold"),
        axis.title = element_text(size = 13, color = "black", face = "bold"),
        axis.text = element_text(size = 10, color = "black", face = "bold"),
        legend.position = "none")+
  theme(plot.margin = margin(0.1,1,0.1,0.1, "cm"))+ #gives space for last x-axis tick label, avoiding cutoff
  labs(y = NULL) 
OS_Plot

# Saving the plots with ggsave
p1 = OS_Plot$plot
p2 = OS_Plot$table
plotp = cowplot::plot_grid(p1,p2,align = "v",axis = "b", ncol =1, rel_heights = c(3,1))
ggsave("outputs/all_Nonmet/OS/22-0909_All_Malignant_Jennifer_categories_OS.png", 
       plot = plotp, 
       width = 7, height = 7,
       dpi = 600)



