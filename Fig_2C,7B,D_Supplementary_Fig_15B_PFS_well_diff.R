# PFS_well_diff.R
#using 50th percentile CAF high and low based on OVERALL cohort

# Load required packages
library(ggplot2)
library(survival)
library(survminer)
library(tidyverse)
library(cowplot)

# Set your working directory
#setwd("_")


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

# If restricting to well-diff, exclude ATC and PDTC 
ClinicalData_Nonmet_Malignant_Welldiff <- ClinicalData_Nonmet_Malignant %>% subset(Diagnosis != "ATC" & Diagnosis != "PDTC")

# Calculating length of PROGRESSION FREE survival from initial therapy complete date
ClinicalData_Nonmet_Malignant_Welldiff$Recur.or.progress.date <- as.Date(ClinicalData_Nonmet_Malignant_Welldiff$Recur.or.progress.date, "%m/%d/%y")
ClinicalData_Nonmet_Malignant_Welldiff$Initial.therapy.complete.date <- as.Date(ClinicalData_Nonmet_Malignant_Welldiff$Initial.therapy.complete.date, "%m/%d/%y")
ClinicalData_Nonmet_Malignant_Welldiff$Date.last.followup <- as.Date(ClinicalData_Nonmet_Malignant_Welldiff$Date.last.followup, "%m/%d/%y")
ClinicalData_Nonmet_Malignant_Welldiff$ProgressionFree_Length <- NA
ClinicalData_Nonmet_Malignant_Welldiff$ProgressionStatus <- 0
for(i in 1:nrow(ClinicalData_Nonmet_Malignant_Welldiff)){
  if(!is.na(ClinicalData_Nonmet_Malignant_Welldiff$Recur.or.progress.date[i])){
    ClinicalData_Nonmet_Malignant_Welldiff$ProgressionFree_Length[i] <- ClinicalData_Nonmet_Malignant_Welldiff$Recur.or.progress.date[i] - 
                                                              ClinicalData_Nonmet_Malignant_Welldiff$Initial.therapy.complete.date[i]
    ClinicalData_Nonmet_Malignant_Welldiff$ProgressionStatus[i] <- 1
  }
  else{
    ClinicalData_Nonmet_Malignant_Welldiff$ProgressionFree_Length[i] <- ClinicalData_Nonmet_Malignant_Welldiff$Date.last.followup[i] - 
    ClinicalData_Nonmet_Malignant_Welldiff$Initial.therapy.complete.date[i]
  }
}

ClinicalData_Nonmet_Malignant_Welldiff$ProgressionFree_Length <- ClinicalData_Nonmet_Malignant_Welldiff$ProgressionFree_Length/(365/12)


# Subset out samples if NA for progressionfree_length
ClinicalData_Nonmet_Malignant_Welldiff <- ClinicalData_Nonmet_Malignant_Welldiff  %>% subset(!is.na(ClinicalData_Nonmet_Malignant_Welldiff$ProgressionFree_Length))


# Look at PFS in known mutations
PFS_TERT.TP53.PIK3CA.mut_Well_DiffMalignant <- ClinicalData_Nonmet_Malignant_Welldiff %>% subset(!is.na(TERT.TP53.PIK3CA.mut))
PFS_TERT.mut_Well_DiffMalignant <- ClinicalData_Nonmet_Malignant_Welldiff %>% subset(!is.na(TERT.mutation))
PFS_TP53.mut_Well_DiffMalignant <- ClinicalData_Nonmet_Malignant_Welldiff %>% subset(!is.na(TP53.mutation))
PFS_PIK3CA.mut_Well_DiffMalignant <- ClinicalData_Nonmet_Malignant_Welldiff %>% subset(!is.na(PIK3CA.mutation))

# Look at PFS in BRS
PFS_BRS_Well_DiffMalignant <- ClinicalData_Nonmet_Malignant_Welldiff %>% subset(!is.na(BRS))

# Look at PFS in CAFs (high vs low)
PFS_CAF_EPIC_Well_DiffMalignant <- ClinicalData_Nonmet_Malignant_Welldiff %>% subset(!is.na(`Cancer associated fibroblast_EPIC_high`))
PFS_CAF_MCPCOUNTER_Well_DiffMalignant <- ClinicalData_Nonmet_Malignant_Welldiff %>% subset(!is.na(`Cancer associated fibroblast_MCPCOUNTER_high`))

# Look at PFS in Aggression score (high vs low)
PFS_Aggression_Well_DiffMalignant <- ClinicalData_Nonmet_Malignant_Welldiff %>% subset(!is.na(`ATS`))

# Look at PFS in Jennifer groups 
PFS_Jennifer_Well_DiffMalignant <- ClinicalData_Nonmet_Malignant_Welldiff %>% subset(!is.na(Aggression.category.detailed))

# PFS Plots
PFS_BRS_Well_DiffMalignant <- mutate(PFS_BRS_Well_DiffMalignant, Cat=ifelse(PFS_BRS_Well_DiffMalignant$BRS > 0, "RAS-Like", "BRAF-Like"))
PFS_TERT.TP53.PIK3CA.mut_Well_DiffMalignant <- mutate(PFS_TERT.TP53.PIK3CA.mut_Well_DiffMalignant, Cat=ifelse(PFS_TERT.TP53.PIK3CA.mut_Well_DiffMalignant$TERT.TP53.PIK3CA.mut == 1, "Mutant", "Not Mutant"))
PFS_TERT.mut_Well_DiffMalignant <- mutate(PFS_TERT.mut_Well_DiffMalignant, Cat=ifelse(PFS_TERT.mut_Well_DiffMalignant$TERT.mutation == 1, "Mutant", "Not Mutant"))
PFS_TP53.mut_Well_DiffMalignant <- mutate(PFS_TP53.mut_Well_DiffMalignant, Cat=ifelse(PFS_TP53.mut_Well_DiffMalignant$TP53.mutation == 1, "Mutant", "Not Mutant"))
PFS_PIK3CA.mut_Well_DiffMalignant <- mutate(PFS_PIK3CA.mut_Well_DiffMalignant, Cat=ifelse(PFS_PIK3CA.mut_Well_DiffMalignant$PIK3CA.mutation == 1, "Mutant", "Not Mutant"))
PFS_CAF_EPIC_Well_DiffMalignant <- mutate(PFS_CAF_EPIC_Well_DiffMalignant, Cat=ifelse(PFS_CAF_EPIC_Well_DiffMalignant$`Cancer associated fibroblast_EPIC_high` == 1, "High", "Low"))
PFS_CAF_MCPCOUNTER_Well_DiffMalignant <- mutate(PFS_CAF_MCPCOUNTER_Well_DiffMalignant, Cat=ifelse(PFS_CAF_MCPCOUNTER_Well_DiffMalignant$`Cancer associated fibroblast_MCPCOUNTER_high` == 1, "High", "Low"))
PFS_Aggression_Well_DiffMalignant <- mutate(PFS_Aggression_Well_DiffMalignant, Cat=ifelse(PFS_Aggression_Well_DiffMalignant$`ATS` >= 0, "High", "Low"))


### TERTp/TP53/PIK3CA
#Survival curve calculation
fit<-survfit(Surv(ProgressionFree_Length, ProgressionStatus) ~ Cat, data= PFS_TERT.TP53.PIK3CA.mut_Well_DiffMalignant)
print(fit)

# Most basic plot
ggsurvplot(fit, data=PFS_TERT.TP53.PIK3CA.mut_Well_DiffMalignant)

# More involved, customized plot
PFS_Plot <- ggsurvplot(fit,
           data=PFS_TERT.TP53.PIK3CA.mut_Well_DiffMalignant,
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
           title = "TP53/TERTp/PI3KCA PFS",
           legend.labs = c("YES", "NO"),
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
ggsave("outputs/well_diff/PFS/22-0909_Well_Diff_TERTp_TP53_PIK3CA_PFS.png", 
       plot = plotp, 
       width = 7, height = 7,
       dpi = 600)



### BRS PFS
#Survival curve calculation
fit<-survfit(Surv(ProgressionFree_Length, ProgressionStatus) ~ Cat, data= PFS_BRS_Well_DiffMalignant)
print(fit)

# Most basic plot
ggsurvplot(fit, data=PFS_BRS_Well_DiffMalignant)

# More involved, customized plot
PFS_Plot <- ggsurvplot(fit,
                       data=PFS_BRS_Well_DiffMalignant,
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
ggsave("outputs/well_diff/PFS/22-0909_Well_Diff_BRS_PFS.png", 
       plot = plotp, 
       width = 7, height = 7,
       dpi = 600)


### TERTp
#Survival curve calculation
fit<-survfit(Surv(ProgressionFree_Length, ProgressionStatus) ~ Cat, data= PFS_TERT.mut_Well_DiffMalignant)
print(fit)

# Most basic plot
ggsurvplot(fit, data=PFS_TERT.mut_Well_DiffMalignant)

# More involved, customized plot
PFS_Plot <- ggsurvplot(fit,
                       data=PFS_TERT.mut_Well_DiffMalignant,
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
                       title = "TERT promoter PFS",
                       legend.labs = c("YES", "NO"),
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
ggsave("outputs/well_diff/PFS/22-0909_Well_Diff_TERTp_PFS.png", 
       plot = plotp, 
       width = 7, height = 7,
       dpi = 600)


### TP53
#Survival curve calculation
fit<-survfit(Surv(ProgressionFree_Length, ProgressionStatus) ~ Cat, data= PFS_TP53.mut_Well_DiffMalignant)
print(fit)

# Most basic plot
ggsurvplot(fit, data=PFS_TP53.mut_Well_DiffMalignant)

# More involved, customized plot
PFS_Plot <- ggsurvplot(fit,
                       data=PFS_TP53.mut_Well_DiffMalignant,
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
                       title = "TP53 PFS",
                       legend.labs = c("YES", "NO"),
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
ggsave("outputs/well_diff/PFS/22-0909_Well_Diff_TP53_PFS.png", 
       plot = plotp, 
       width = 7, height = 7,
       dpi = 600)


### PIK3CA
#Survival curve calculation
fit<-survfit(Surv(ProgressionFree_Length, ProgressionStatus) ~ Cat, data= PFS_PIK3CA.mut_Well_DiffMalignant)
print(fit)

# Most basic plot
ggsurvplot(fit, data=PFS_PIK3CA.mut_Well_DiffMalignant)

# More involved, customized plot
PFS_Plot <- ggsurvplot(fit,
                       data=PFS_PIK3CA.mut_Well_DiffMalignant,
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
                       title = "PIK3CA PFS",
                       legend.labs = c("YES", "NO"),
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
ggsave("outputs/well_diff/PFS/22-0909_Well_Diff_PIK3CA_PFS.png", 
       plot = plotp, 
       width = 7, height = 7,
       dpi = 600)


### CAF EPIC high/low
#Survival curve calculation
fit<-survfit(Surv(ProgressionFree_Length, ProgressionStatus) ~ Cat, data= PFS_CAF_EPIC_Well_DiffMalignant)
print(fit)

# Most basic plot
ggsurvplot(fit, data=PFS_CAF_EPIC_Well_DiffMalignant)

# More involved, customized plot
PFS_Plot <- ggsurvplot(fit,
                       data=PFS_CAF_EPIC_Well_DiffMalignant,
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
                       title = "Cancer-Associated Fibroblast\nEPIC PFS",
                       legend.labs = c("High", "Low"),
                       tables.height = .33,
                       break.time.by = 24,
                       palette = c("#02f02e", "#48827a", "#86AA00"))+ # Can switch plot color palette order
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
ggsave("outputs/well_diff/PFS/22-0909_Well_Diff_CAF_EPIC_PFS.png", 
       plot = plotp, 
       width = 7, height = 7,
       dpi = 600)


### CAF MCPCOUNTER high/low
#Survival curve calculation
fit<-survfit(Surv(ProgressionFree_Length, ProgressionStatus) ~ Cat, data= PFS_CAF_MCPCOUNTER_Well_DiffMalignant)
print(fit)

# Most basic plot
ggsurvplot(fit, data=PFS_CAF_MCPCOUNTER_Well_DiffMalignant)

# More involved, customized plot
PFS_Plot <- ggsurvplot(fit,
                       data=PFS_CAF_MCPCOUNTER_Well_DiffMalignant,
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
                       title = "Cancer-Associated Fibroblast\nMCPCOUNTER PFS",
                       legend.labs = c("High", "Low"),
                       tables.height = .33,
                       break.time.by = 24,
                       palette = c("#02f02e", "#48827a", "#86AA00"))+ # Can switch plot color palette order
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
ggsave("outputs/well_diff/PFS/22-0909_Well_Diff_CAF_MCPCOUNTER_PFS.png", 
       plot = plotp, 
       width = 7, height = 7,
       dpi = 600)


### Aggression score high/low
#Survival curve calculation
fit<-survfit(Surv(ProgressionFree_Length, ProgressionStatus) ~ Cat, data= PFS_Aggression_Well_DiffMalignant)
print(fit)

# Most basic plot
ggsurvplot(fit, data=PFS_Aggression_Well_DiffMalignant)

# More involved, customized plot
PFS_Plot <- ggsurvplot(fit,
                       data=PFS_Aggression_Well_DiffMalignant,
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
                       legend.title = "ATS",
                       title = "Aggression Score PFS",
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
ggsave("outputs/well_diff/PFS/22-0909_Well_Diff_Aggression_PFS.png", 
       plot = plotp, 
       width = 7, height = 7,
       dpi = 600)
