# Author: Matthew Aaron Loberg
# Script: Score Exploration by Histotype (BRS, TDS, PI3K, ERK)
# Original script date: 22-0509

## 22-1004 Update
# Removing HTPTC from TDS
# Removing HTPTC from TDS ERK correlation
# Updating to latest clinical file (9-8-22)

# Load required packages
library(tidyverse)

# Read in data
ClinicalData <- read_csv(file = "data_in_use/VUMC.cohort.GX_9-8-22_V2.csv") # Latest version of the clinical data file
ClinicalData <- ClinicalData %>% subset(Diagnosis != "normal")
ClinicalData$BRS <- as.numeric(ClinicalData$BRS)

# Merge samples to make a new diagnosis column 
# Merge PTCs and IFVPTCs, NIFTPs and EFVPTCs
ClinicalData$Diagnosis_Merged <- ClinicalData$Diagnosis.with.iFVPTC.and.eFVPTC
for(i in 1:nrow(ClinicalData)){
  if(ClinicalData$Diagnosis_Merged[i] == "PTC" | ClinicalData$Diagnosis_Merged[i] == "IFVPTC"){
    ClinicalData$Diagnosis_Merged[i] <- "IFVPTC+\nPTC"
  }
  else if(ClinicalData$Diagnosis_Merged[i] == "NIFTP" | ClinicalData$Diagnosis_Merged[i] == "EFVPTC"){
    ClinicalData$Diagnosis_Merged[i] <- "NIFTP+\nEFVPTC"
  }
  else if(ClinicalData$Diagnosis_Merged[i] == "FA" | ClinicalData$Diagnosis_Merged[i] == "HA"){
    ClinicalData$Diagnosis_Merged[i] <- "FA/HA"
  }
  else if(ClinicalData$Diagnosis_Merged[i] == "FTC" | ClinicalData$Diagnosis_Merged[i] == "HC"){
    ClinicalData$Diagnosis_Merged[i] <- "FTC/HC"
  }
}

# BRS Score Plot
ClinicalData_BRS <- ClinicalData %>% subset(!is.na(ClinicalData$BRS) & Diagnosis_Merged != "MNG" & Diagnosis_Merged != "HT")
ClinicalData_BRS_Local <- ClinicalData_BRS %>% subset(Location.type == "Primary" | Location.type == "Localdisease")
# All BRS
plot <- ggplot(ClinicalData_BRS, aes(Diagnosis_Merged, BRS)) +
  geom_boxplot(outlier.size = -1, 
               aes(),
               alpha = 0.2, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5, 
              alpha = 0.9,
              show.legend = TRUE) +
  labs (x = "Diagnosis", y = "BRAF-RAS Score") + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 24), 
    axis.title.y = element_text(face = "bold", size = 24),
    axis.text = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="Diagnosis", limits = c("FA/HA", "FTC/HC", "NIFTP+\nEFVPTC", "IFVPTC+\nPTC", "PDTC", "ATC")) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red")
ggsave("outputs/22-1004_Plotting_Scores_BRS_TDS_ERK_PI3K/22-1004_BRS_All_Samples.png",
       width = 7,
       height = 5,
       plot, dpi = 600)
# Local only BRS
plot <- ggplot(ClinicalData_BRS_Local, aes(Diagnosis_Merged, BRS)) +
  geom_boxplot(outlier.size = -1, 
               aes(),
               alpha = 0.2, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5, 
              alpha = 0.9,
              show.legend = TRUE) +
  labs (x = "Diagnosis", y = "BRAF-RAS Score") + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 24), 
    axis.title.y = element_text(face = "bold", size = 24),
    axis.text = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="Diagnosis", limits = c("FA/HA", "FTC/HC", "NIFTP+\nEFVPTC", "IFVPTC+\nPTC", "PDTC", "ATC")) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red")
ggsave("outputs/22-1004_Plotting_Scores_BRS_TDS_ERK_PI3K/22-1004_BRS_Local_Samples.png",
       width = 7,
       height = 5,
       plot, dpi = 600)

# TDS Score Plot
# Note: Including MNG for TDS -> though it seems that TDS wasn't calculated for MNG? (so removed by the is.na command)
# Note: adjusting this from 22-0509 -> removing HTPTCs
ClinicalData$TDS <- as.numeric(ClinicalData$TDS)
ClinicalData_TDS <- ClinicalData %>% subset(!is.na(ClinicalData$TDS) & Diagnosis_Merged != "HT" & Diagnosis.with.HTPTC != "HTPTC")
ClinicalData_TDS_Local <- ClinicalData_TDS %>% subset(Location.type == "Primary" | Location.type == "Localdisease")
# All TDS
plot <- ggplot(ClinicalData_TDS, aes(Diagnosis_Merged, TDS)) +
  geom_boxplot(outlier.size = -1, 
               aes(),
               alpha = 0.2, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5, 
              alpha = 0.9,
              show.legend = TRUE) +
  labs (x = "Diagnosis", y = "Thyroid Differentation Score") + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 24), 
    axis.title.y = element_text(face = "bold", size = 22),
    axis.text = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="Diagnosis", limits = c("FA/HA", "FTC/HC", "NIFTP+\nEFVPTC", "IFVPTC+\nPTC", "PDTC", "ATC")) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red")
ggsave("outputs/22-1004_Plotting_Scores_BRS_TDS_ERK_PI3K/22-1004_TDS_All_Samples_NO_HTPTC.png",
       width = 7,
       height = 5,
       plot, dpi = 600)

# Local only TDS
plot <- ggplot(ClinicalData_TDS_Local, aes(Diagnosis_Merged, TDS)) +
  geom_boxplot(outlier.size = -1, 
               aes(),
               alpha = 0.2, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5, 
              alpha = 0.9,
              show.legend = TRUE) +
  labs (x = "Diagnosis", y = "Thyroid Differentation Score") + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 24), 
    axis.title.y = element_text(face = "bold", size = 22),
    axis.text = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="Diagnosis", limits = c("FA/HA", "FTC/HC", "NIFTP+\nEFVPTC", "IFVPTC+\nPTC", "PDTC", "ATC")) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red")
ggsave("outputs/22-1004_Plotting_Scores_BRS_TDS_ERK_PI3K/22-1004_TDS_Local_Samples_NO_HTPTC.png",
       width = 7,
       height = 5,
       plot, dpi = 600)

# ERK Score Plot
# Note: Including MNG for ERK -> though it seems that ERK wasn't calculated for MNG?
ClinicalData$ERK <- as.numeric(ClinicalData$ERK)
ClinicalData_ERK <- ClinicalData %>% subset(!is.na(ClinicalData$ERK) & Diagnosis_Merged != "HT")
ClinicalData_ERK_Local <- ClinicalData_ERK %>% subset(Location.type == "Primary" | Location.type == "Localdisease")
# All ERK
plot <- ggplot(ClinicalData_ERK, aes(Diagnosis_Merged, ERK)) +
  geom_boxplot(outlier.size = -1, 
               aes(),
               alpha = 0.2, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5, 
              alpha = 0.9,
              show.legend = TRUE) +
  labs (x = "Diagnosis", y = "ERK Activity Score") + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 24), 
    axis.title.y = element_text(face = "bold", size = 24),
    axis.text = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="Diagnosis", limits = c("FA/HA", "FTC/HC", "NIFTP+\nEFVPTC", "IFVPTC+\nPTC", "PDTC", "ATC")) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red")
ggsave("outputs/22-1004_Plotting_Scores_BRS_TDS_ERK_PI3K/22-1004_ERK_All_Samples.png",
       width = 7,
       height = 5,
       plot, dpi = 600)

# Local only ERK
plot <- ggplot(ClinicalData_ERK_Local, aes(Diagnosis_Merged, ERK)) +
  geom_boxplot(outlier.size = -1, 
               aes(),
               alpha = 0.2, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5, 
              alpha = 0.9,
              show.legend = TRUE) +
  labs (x = "Diagnosis", y = "ERK Activity Score") + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 24), 
    axis.title.y = element_text(face = "bold", size = 24),
    axis.text = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="Diagnosis", limits = c("FA/HA", "FTC/HC", "NIFTP+\nEFVPTC", "IFVPTC+\nPTC", "PDTC", "ATC")) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red")
ggsave("outputs/22-1004_Plotting_Scores_BRS_TDS_ERK_PI3K/22-1004_ERK_Local_Samples.png",
       width = 7,
       height = 5,
       plot, dpi = 600)

## PI3K-AKT-MTOR Plot
ClinicalData$PI3K_AKT_MTOR <- as.numeric(ClinicalData$PI3K_AKT_MTOR)
ClinicalData_PI3K <- ClinicalData %>% subset(!is.na(ClinicalData$PI3K_AKT_MTOR) & Diagnosis_Merged != "HT")
ClinicalData_PI3K_Local <- ClinicalData_PI3K %>% subset(Location.type == "Primary" | Location.type == "Localdisease")
# All PI3K
plot <- ggplot(ClinicalData_PI3K, aes(Diagnosis_Merged, PI3K_AKT_MTOR)) +
  geom_boxplot(outlier.size = -1, 
               aes(),
               alpha = 0.2, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5, 
              alpha = 0.9,
              show.legend = TRUE) +
  labs (x = "Diagnosis", y = "PI3K-AKT-mTOR Activity Score") + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 24), 
    axis.title.y = element_text(face = "bold", size = 20),
    axis.text = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="Diagnosis", limits = c("FA/HA", "FTC/HC", "NIFTP+\nEFVPTC", "IFVPTC+\nPTC", "PDTC", "ATC")) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red")
ggsave("outputs/22-1004_Plotting_Scores_BRS_TDS_ERK_PI3K/22-1004_PI3K_AKT_mTOR_All_Samples.png",
       width = 7,
       height = 5,
       plot, dpi = 600)

# Local only PI3K
plot <- ggplot(ClinicalData_PI3K_Local, aes(Diagnosis_Merged, PI3K_AKT_MTOR)) +
  geom_boxplot(outlier.size = -1, 
               aes(),
               alpha = 0.2, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5, 
              alpha = 0.9,
              show.legend = TRUE) +
  labs (x = "Diagnosis", y = "PI3K-AKT-mTOR Activity Score") + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 24), 
    axis.title.y = element_text(face = "bold", size = 20),
    axis.text = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="Diagnosis", limits = c("FA/HA", "FTC/HC", "NIFTP+\nEFVPTC", "IFVPTC+\nPTC", "PDTC", "ATC")) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red")
ggsave("outputs/22-1004_Plotting_Scores_BRS_TDS_ERK_PI3K/22-1004_PI3K_Local_Samples.png",
       width = 7,
       height = 5,
       plot, dpi = 600)

## BRS TDS Association
# Linear Model
# Change to numerics
ClinicalData_BRS$TDS <- as.numeric(ClinicalData_BRS$TDS)
ClinicalData_BRS$ERK <- as.numeric(ClinicalData_BRS$ERK)
ClinicalData_BRS$PI3K_AKT_MTOR <- as.numeric(ClinicalData_BRS$PI3K_AKT_MTOR)
# Maker BRS/TDS linear model
All.lm.BRS_TDS <- lm(TDS ~ BRS, ClinicalData_BRS)
# Make BRS/TDS plot
Plot <- ggplot(ClinicalData_BRS, aes(BRS, TDS))+
  geom_jitter(alpha = 0.6, size = 2) + 
  theme_classic() + 
  xlab("BRAF-RAS Score") + 
  ylab("Thyroid Differentiation Score") + 
  theme(panel.border = element_rect(colour = "black", size = 4, fill = NA),
        axis.title.x = element_text(face = "bold", size = 20),
        axis.title.y = element_text(face = "bold", size = 18),
        axis.text = element_text(face = "bold", size = 15)) + 
  geom_abline(slope = coef(All.lm.BRS_TDS)[["BRS"]],
              intercept = coef(All.lm.BRS_TDS)[["(Intercept)"]],
              lwd = 2) 
ggsave("outputs/22-1004_Plotting_Scores_BRS_TDS_ERK_PI3K/22-1004_BRS_TDS_Correlation.png",
       width = 5,
       height = 4,
       Plot, dpi = 600)

# Make a new BRS TDS associaiton plot with HTPTC excluded
ClinicalData_BRS_TDS <- ClinicalData_BRS %>% subset(Diagnosis.with.HTPTC != "HTPTC")
# Maker BRS/TDS linear model
All.lm.BRS_TDS <- lm(TDS ~ BRS, ClinicalData_BRS_TDS)
# Make BRS/TDS plot
Plot <- ggplot(ClinicalData_BRS_TDS, aes(BRS, TDS))+
  geom_jitter(alpha = 0.6, size = 2) + 
  theme_classic() + 
  xlab("BRAF-RAS Score") + 
  ylab("Thyroid Differentiation Score") + 
  theme(panel.border = element_rect(colour = "black", size = 4, fill = NA),
        axis.title.x = element_text(face = "bold", size = 20),
        axis.title.y = element_text(face = "bold", size = 18),
        axis.text = element_text(face = "bold", size = 15)) + 
  geom_abline(slope = coef(All.lm.BRS_TDS)[["BRS"]],
              intercept = coef(All.lm.BRS_TDS)[["(Intercept)"]],
              lwd = 2) 
ggsave("outputs/22-1004_Plotting_Scores_BRS_TDS_ERK_PI3K/22-1004_BRS_TDS_Correlation_NO_HTPTC.png",
       width = 5,
       height = 4,
       Plot, dpi = 600)


# Maker BRS/ERK linear model
All.lm.BRS_ERK <- lm(ERK ~ BRS, ClinicalData_BRS)
# Make BRS/ERK Plot
Plot <- ggplot(ClinicalData_BRS, aes(BRS, ERK))+
  geom_jitter(alpha = 0.6, size = 2) + 
  theme_classic() + 
  xlab("BRAF-RAS Score") + 
  ylab("ERK Activity Score") + 
  theme(panel.border = element_rect(colour = "black", size = 4, fill = NA),
        axis.title.x = element_text(face = "bold", size = 20),
        axis.title.y = element_text(face = "bold", size = 20),
        axis.text = element_text(face = "bold", size = 15)) + 
  geom_abline(slope = coef(All.lm.BRS_ERK)[["BRS"]],
              intercept = coef(All.lm.BRS_ERK)[["(Intercept)"]],
              lwd = 2) 
ggsave("outputs/22-1004_Plotting_Scores_BRS_TDS_ERK_PI3K/22-1004_BRS_ERK_Correlation.png",
       width = 5,
       height = 4,
       Plot, dpi = 600)
# Make BRS/PI3K_AKT_mTOR Linear Model
All.lm.BRS_PI3K <- lm(PI3K_AKT_MTOR ~ BRS, ClinicalData_BRS)
# Make BRS/PI3K_AKT_mTOR Plot
Plot <- ggplot(ClinicalData_BRS, aes(BRS, PI3K_AKT_MTOR))+
  geom_jitter(alpha = 0.6, size = 2) + 
  theme_classic() + 
  xlab("BRAF-RAS Score") + 
  ylab("PI3K-AKT-mTOR Activity Score") + 
  theme(panel.border = element_rect(colour = "black", size = 4, fill = NA),
        axis.title.x = element_text(face = "bold", size = 20),
        axis.title.y = element_text(face = "bold", size = 16),
        axis.text = element_text(face = "bold", size = 15)) + 
  geom_abline(slope = coef(All.lm.BRS_PI3K)[["BRS"]],
              intercept = coef(All.lm.BRS_PI3K)[["(Intercept)"]],
              lwd = 2) 
ggsave("outputs/22-1004_Plotting_Scores_BRS_TDS_ERK_PI3K/22-1004_BRS_PI3K_Correlation.png",
       width = 5,
       height = 4,
       Plot, dpi = 600)

# Note: can retrieve R-squared and P-values from linear models with the following command: 
# summary(All.lm.BRS_TDS)
