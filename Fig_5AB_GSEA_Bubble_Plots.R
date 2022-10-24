### Figure 5A GSEA Bubble Plots: Info
# Data from GSEA run on Broad institute program


### Load Packages
library(ggplot2)
library(dplyr)
library(stringr)
library(forcats)
library(RColorBrewer)

### Load data
dot_data_Local <- read.csv("Data/22-1006_Local_Results.csv")

### Subset data
dot_data_Local_normal <- dot_data_Local[c(1,7,10),] # Isolate normal fibroblast rows (note, removing PDTC)
dot_data_Local_myCAF <- dot_data_Local[c(2,8,11),] # Isolate myCAF rows (note, removing PDTC)
dot_data_Local_iCAF <- dot_data_Local[c(3,9,12),] # Isolate iCAF rows (note, removing PDTC)
  
### Show just normal fibroblasts
  dot_data_Local_normal$GeneSet <- c("ATC", "PTC/IFVPTC", "FTC/HC")
  dot_plot <- ggplot(dot_data_Local_normal, aes(NES, fct_reorder(GeneSet, Order))) + 
    geom_point(aes(size = GeneRatio, fill = FDR_qvalue), shape = 21, stroke = 1) +
    theme_bw(base_size = 10) +
    scale_fill_distiller(limits=c(0,1), type="div", palette= "RdYlBu", direction = 1, "FDR q Value") +
    scale_size_continuous(range = c(0, 20), limits = c(0,1), "Gene Ratio (tags)") +
    scale_x_continuous(name="NES", seq(-2, 2.0, by = 1), NULL, expand = c(0,0), limits=c(-2.3, 2.3), NULL) +
    theme(aspect.ratio = 1,
          axis.title = element_text(face = "bold", size = 13),
          axis.text = element_text(face = "bold", size = 12),
          title = element_text(face = "bold", size = 11)) +
    labs(y = "") +
    ggtitle("Normal Fibroblast GSEA\nRelative to MNG") + 
    geom_vline(xintercept = 0, size = 0.25)
  ggsave("ggplot_output/22-1007_ggplot/22-1007_GSEA_Normal_Fibroblast_Local_Thyroid_Malignancies.png", 
         width = 5, height = 5,
         dot_plot, dpi = 600)
  
  ### Show just myofibroblasts
  dot_data_Local_myCAF$GeneSet <- c("ATC", "PTC/IFVPTC", "FTC/HC")
  dot_plot <- ggplot(dot_data_Local_myCAF, aes(NES, fct_reorder(GeneSet, Order))) + 
    geom_point(aes(size = GeneRatio, fill = FDR_qvalue), shape = 21, stroke = 1) +
    theme_bw(base_size = 10) +
    scale_size_continuous(range = c(0, 20), limits = c(0,1), "Gene Ratio (tags)") +
    scale_fill_distiller(limits=c(0,1), type="div", palette= "RdYlBu", direction = 1, "FDR q Value") +
    scale_x_continuous(name="NES", seq(-2, 2.0, by = 1), NULL, expand = c(0,0), limits=c(-2.3, 2.3), NULL) +
    theme(aspect.ratio = 1,
          axis.title = element_text(face = "bold", size = 13),
          axis.text = element_text(face = "bold", size = 12),
          title = element_text(face = "bold", size = 11)) +
    ggtitle("myCAF GSEA\nRelative to MNG") + 
    geom_vline(xintercept = 0, size = 0.25) + 
    labs(y = "")
  ggsave("ggplot_output/22-1007_ggplot/22-1007_GSEA_myCAF_Local_Thyroid_Malignancies.png", 
         width = 5, height = 5,
         dot_plot, dpi = 600)
  
  ### Show just iCAFs
  dot_data_Local_iCAF$GeneSet <- c("ATC", "PTC/IFVPTC", "FTC/HC")
  dot_plot <- ggplot(dot_data_Local_iCAF, aes(NES, fct_reorder(GeneSet, Order))) + 
    geom_point(aes(size = GeneRatio, fill = FDR_qvalue), shape = 21, stroke = 1) +
    theme_bw(base_size = 10) +
    scale_size_continuous(range = c(0, 20), limits = c(0,1), "Gene Ratio (tags)") +
    scale_fill_distiller(limits=c(0,1), type="div", palette= "RdYlBu", direction = 1, "FDR q Value") +
    scale_x_continuous(name="NES", seq(-2, 2.0, by = 1), NULL, expand = c(0,0), limits=c(-2.3, 2.3), NULL) +
    theme(aspect.ratio = 1,
          axis.title = element_text(face = "bold", size = 13),
          axis.text = element_text(face = "bold", size = 12),
          title = element_text(face = "bold", size = 11)) +
    ggtitle("iCAF GSEA\nRelative to MNG") + 
    geom_vline(xintercept = 0, size = 0.25) + 
    labs(y = "")
  ggsave("ggplot_output/22-1007_ggplot/22-1007_GSEA_iCAF_Local_Thyroid_Malignancies.png", 
         width = 5, height = 5,
         dot_plot, dpi = 600)



#### Figure 5B: ATC-BRAF vs ATC-RAS 
ATC_dot_data <- read.csv(file = "Data/22-0324_Local_ATC_BRAF_vs_RAS.csv")

# Plot ATC dot data
dot_plot <- ggplot(ATC_dot_data, aes(NES, fct_reorder(GeneSet, Order))) + 
  geom_point(aes(size = GeneRatio, fill = FDR_qvalue), shape = 21, stroke = 1) +
  theme_bw(base_size = 10) +
  scale_size_continuous(range = c(0, 15), limits = c(0,1), "Gene Ratio (tags)") +
  scale_fill_distiller(limits=c(0,1), type="div", palette= "RdYlBu", direction = 1, "FDR q Value") +
  scale_x_continuous(name="NES", seq(-3, 3, by = 1), NULL, expand = c(0,0), limits=c(-3.7, 3.7), NULL) +
  theme(aspect.ratio = 1,
        axis.title = element_text(face = "bold", size = 10),
        axis.text = element_text(face = "bold", size = 10),
        title = element_text(face = "bold", size = 9)) +
  ggtitle("RAS-Like \t\t\t\tBRAF-Like\nLocal ATCs \tLocal ATCs") + 
  geom_vline(xintercept = 0, size = 0.25) + 
  labs(y = "")
ggsave("ggplot_output/22-0324_ggplot/22-0324_GSEA_ATC_BRAF-Like_RAS-Like_Local_Thyroid_Malignancies.png", 
       width = 4, height = 4,
       dot_plot, dpi = 600)
  
