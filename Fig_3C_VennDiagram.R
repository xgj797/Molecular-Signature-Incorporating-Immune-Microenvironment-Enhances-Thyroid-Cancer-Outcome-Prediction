### Venn Diagram: Information
# Load in the aggressive vs indolent and the BRAF vs RAS overlap lists and make venn diagrams of the overla

### load packages
library(tidyverse)
library(ggvenn)

### Load data sets
All_DESEQ_BRAF <- readRDS(file = "data_in_use/22-0908_DESEQ2_All_Primary_LocalDisease_MG_HT_Excluded_resOrdered_BRAF_RAS.rds")
All_DESEQ_Aggressive <- readRDS(file = "data_in_use/22-0908_DESEQ2_All_Primary_LocalDisease_MG_HT_Excluded_resOrdered_Aggressive_Indolent.rds")

### subset out to BRAF-associated genes with an adjusted p-value < 0.05 and log2FoldChange > 1
All_DESEQ_BRAF_Upregulated <- All_DESEQ_BRAF %>% subset(log2FoldChange > 1) 
All_DESEQ_BRAF_Upregulated <- All_DESEQ_BRAF_Upregulated %>% subset(padj < 0.05)

# ##Subset out the aggressive associated-genes with an adjusted p-value < 0.05 and log2FoldChange > 2
All_DESEQ_Aggressive_Upregulated <- All_DESEQ_Aggressive %>% subset(log2FoldChange > 2)
All_DESEQ_Aggressive_Upregulated <- All_DESEQ_Aggressive_Upregulated %>% subset(padj < 0.05)

### Create a list of the genes of the overlap between BRAF and Aggressive
Aggressive_BRAF_Overlap <- All_DESEQ_BRAF_Upregulated
Aggressive_BRAF_Overlap$Genes <- rownames(Aggressive_BRAF_Overlap)
Aggressive_BRAF_Overlap <- Aggressive_BRAF_Overlap %>% subset(Genes %in% rownames(All_DESEQ_Aggressive_Upregulated))

Aggressive_BRAF_OverlapList <- as_tibble(Aggressive_BRAF_Overlap$Genes)
Aggressive_BRAF_OverlapList # Print out to test

write_csv(Aggressive_BRAF_OverlapList, file = "data_in_use/22-0927_All_BRAF_All_Aggressive_Overlap_Genes.csv")

# Comparing BRAF all vs Aggressive all
Overlap_BRAF_Aggressive_List <- list('BRAF-Like\nUpregulated' = rownames(All_DESEQ_BRAF_Upregulated), 'Aggressive\nUpregulated' = rownames(All_DESEQ_Aggressive_Upregulated))
#Make the Venn Diagram
BRAF_AggressiveVenn <- ggvenn(Overlap_BRAF_Aggressive_List, c("BRAF-Like\nUpregulated", "Aggressive\nUpregulated"), fill_color = c("red", "grey"), text_size = 10, set_name_size = 12)
# Change fonts
BRAF_AggressiveVenn <- BRAF_AggressiveVenn + theme(text = element_text(face = "bold"), 
                                                     axis.title = element_text(face = "bold"))
# ggsave the venn diagram
ggsave("outputs/22-0927_DESEQ_Results_Analysis_Venn/22-0927_BRAF_Aggressive_Venn_549_Gene_Set.png", BRAF_AggressiveVenn, dpi = 600) 

# Rerunning w/ RAS-like overlap 
# Subset out the poor-outcome associated-genes with an adjusted p-value < 0.05 & log2FoldChange >1
All_DESEQ_Aggressive_Upregulated <- All_DESEQ_Aggressive %>% subset(log2FoldChange > 1)
All_DESEQ_Aggressive_Upregulated <- All_DESEQ_Aggressive_Upregulated %>% subset(padj < 0.05)
# SUbset out to RAS-associated genes with an adjusted p-value < 0.05, log-fold change > 1
All_DESEQ_RAS_Upregulated <- All_DESEQ_BRAF %>% subset(log2FoldChange < -1 & padj < 0.05)
# Comparing RAS all vs Aggressive Outcome all
Overlap_RAS_Aggressive <- list('RAS-Like\nUpregulated' = rownames(All_DESEQ_RAS_Upregulated), 'Aggressive\nUpregulated' = rownames(All_DESEQ_Aggressive_Upregulated))
# Make the Venn Diagram
RAS_AggressiveVenn <- ggvenn(Overlap_RAS_Aggressive, c("RAS-Like\nUpregulated", "Aggressive\nUpregulated"), fill_color = c("blue", "grey"), text_size = 10, set_name_size = 12)
# ggsave the venn diagram
ggsave("outputs/22-0927_DESEQ_Results_Analysis_Venn/22-0927_RAS_Aggressive_Venn.png", RAS_AggressiveVenn, dpi = 600) 
