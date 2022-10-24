#Immune data in ComplexHeatmap - TCGA data

library(ComplexHeatmap)


setwd("_")
data = read.csv("TCGA_data.csv", header = TRUE)
annot = read.csv("TCGA_annotation_data.csv", header = TRUE)
#don't forget to remove first column's header before loading into R, make sure no extra comma at end of first row

data <- scale(data)

#convert to dataframe
data.df <- as.data.frame(data)

#invert BRS data so that color scale will be reversed, and BRAF-like samples will appear blue
data.df$BRS <- -(data.df$BRS)

#convert from datafram to matrix needed for heatmap3
data.ma <- as.matrix(data.df)


HA_Left <- rowAnnotation(HISTOLOGICAL_TYPE_COMPLEX = annot$HISTOLOGICAL_TYPE_COMPLEX, # Note: Can add Diagnosis as a label instead of an annotation...gives the ability to write in the box
                            'BRAF_Mutation' = annot$BRAF_Mutation,
                            'RAS_Mutation' = annot$RAS_Mutation,
                            col = list(HISTOLOGICAL_TYPE_COMPLEX = c("Classical" = "Red",  # These are the colors for adding Diagnosis as an annotation. Remove if adding as a label.
                                                     "Columnar Cell" = "Orange", 
                                                     "Cribriform Morular" = "Yellow",
                                                     "Follicular" = "lightskyblue",
                                                     "Diffuse Sclerosing" = "magenta",
                                                     "Tall Cell" = "Purple",
                                                     "Mixed" = "Grey",
                                                     "Unknown" = "Black"),
                                       'BRAF_Mutation' = c("YES" = "Red",
                                                       "NO" = "Black"),
                                       'RAS_Mutation' = c("YES" = "lightskyblue",
                                                       "NO" = "Black")),
                            annotation_legend_param = list(HISTOLOGICAL_TYPE_COMPLEX=list(at=c("Classical","Columnar Cell","Cribriform Morular","Follicular","Diffuse Sclerosing","Tall Cell","Mixed","Unknown")),
                                                           'BRAF_Mutation'=list(at=c("YES", "NO")),
                                                           'RAS_Mutation'=list(at=c("YES", "NO")))) #recommend by Matt 4-5-22, see https://jokergoo.github.io/ComplexHeatmap-reference/book/


#heatmap, sample cluster only
Heatmap(data.ma, name = "Scaled score", width = unit(12, "cm"), height = unit(12, "cm"), #adjust size of heatmap
        show_row_names = FALSE, #set false to hid row (sample) names
        column_order = 1:33, #order columns as they are in the annot input
        column_names_gp = gpar(fontsize = 8), #adjust size of column labels
        #column_names_rot = 45, #adjust size and angle of column labels
        left_annotation = HA_Left) #put diagnosis annotation on the left


