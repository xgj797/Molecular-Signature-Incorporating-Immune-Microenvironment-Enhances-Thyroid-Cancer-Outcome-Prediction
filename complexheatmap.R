#Immune data in ComplexHeatmap.R

library(ComplexHeatmap)


setwd("/Users/georgexu/Dropbox/Vanderbilt\ Grad\ School/Lab/Projects/Data/Immune\ Clustering/9-30-22\ -\ immune\ heatmap\ with\ mets\ added\ -\ HTPTC\ included")
data = read.csv("6_TIMER_and_2_CAFs_with_cibersort_and_TIDE_include_mets_and_HTPTC.csv", header = TRUE)
annot = read.csv("annotation_data_9-30-22_include_mets_and_HTPTC.csv", header = TRUE)
#don't forget to remove first column's header before loading into R

data <- scale(data)

#convert from matrix to dataframe
data.df <- as.data.frame(data)

#invert BRS data so that color scale will be reversed, and BRAF-like samples will appear blue
data.df$BRS <- -(data.df$BRS)

#convert from datafram to matrix needed for heatmap3
data.ma <- as.matrix(data.df)

HA_Left <- rowAnnotation(Diagnosis = annot$Diagnosis, # Note: Can add Diagnosis as a label instead of an annotation...gives the ability to write in the box
                            'Location' = annot$Location.type,
                            col = list(Diagnosis = c("ATC" = "magenta",  # These are the colors for adding Diagnosis as an annotation. Remove if adding as a label.
                                                     "PDTC" = "purple", 
                                                     "PTC and IFVPTC" = "red",
                                                     "HT with PTC" = "orange",
                                                     "FTC" = "blue"),
                                       'Location' = c("Non-metastatic" = "Black",
                                                      "Local metastasis" = "Green",
                                                      "Distant metastasis" = "Red")),
                            annotation_legend_param = list(Diagnosis=list(at=c("ATC","PDTC","PTC and IFVPTC","HT with PTC","FTC")),
                                                           'Location'=list(at=c("Non-metastatic","Local metastasis", "Distant metastasis")))) #recommend by Matt 4-5-22, see https://jokergoo.github.io/ComplexHeatmap-reference/book/

#heatmap, sample cluster only
Heatmap(data.ma, name = "Scaled score", width = unit(12, "cm"), height = unit(12, "cm"), #adjust size of heatmap
        show_row_names = FALSE, #set false to hid row (sample) names
        column_order = 1:34, #order columns as they are in the annot input
        column_names_gp = gpar(fontsize = 8), #adjust size of column labels
        left_annotation = HA_Left) #put diagnosis annotation on the left


