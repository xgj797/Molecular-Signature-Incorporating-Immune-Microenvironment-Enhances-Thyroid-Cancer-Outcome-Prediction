#Immune data in ComplexHeatmap.R

library(ComplexHeatmap)
library(colorRamp2)

setwd("_")
data = read.csv("6_TIMER_and_2_CAFs_with_cibersort_and_TIDE_include_mets_exclude_HTPTC_sort_MAP_within_diagnosis.csv", header = TRUE)
annot = read.csv("annotation_data_4-27-23_include_mets_exclude_HTPTC_sort_MAP_within_diagnosis_MAP_categorical.csv", header = TRUE)

data <- scale(data)

#convert from matrix to dataframe
data.df <- as.data.frame(data)

#convert from datafram to matrix needed for heatmap3
data.ma <- as.matrix(data.df)

HA_Left <- rowAnnotation(Diagnosis = annot$Diagnosis, # Note: Can add Diagnosis as a label instead of an annotation...gives the ability to write in the box
                         'Location' = annot$Location.type,
                         'Aggressive.disease' = annot$Aggressive.disease,
                         'MAP.category' = annot$MAP.category,
                         'MAP' = annot$MAP,
                         col = list(Diagnosis = c("ATC" = "magenta",  # These are the colors for adding Diagnosis as an annotation. Remove if adding as a label.
                                                  "PDTC" = "purple", 
                                                  "PTC and IFVPTC" = "red",
                                                  "FTC" = "blue"),
                                    'Location' = c("Local disease" = "Black",
                                                   "Regional LN met" = "Green",
                                                   "Distant met" = "Red"),
                                    'Aggressive.disease' = c("Aggressive" = "mediumorchid",
                                                             "Benign" = "seashell"),
                                    'MAP.category' = c("Negative" = "cyan",
                                                       "Positive" = "magenta"),
                                    'MAP' = colorRamp2(c(-0.5,2.5), c("white", "black"))),
                         #'MAP' = colorRamp2(c(-0.5,0,2.5), c("cyan", "white", "magenta"))), old version
                         annotation_legend_param = list(Diagnosis=list(at=c("ATC","PDTC","PTC and IFVPTC","FTC")),
                                                        'Location'=list(at=c("Local disease","Regional LN met", "Distant met"))))

#heatmap, sample cluster only
Heatmap(data.ma, name = "Scaled score", width = unit(12, "cm"), height = unit(12, "cm"), #adjust size of heatmap
        show_row_names = FALSE, #set false to hid row (sample) names
        row_order = 1:180, #set order manually, sorted on MAP score in input file
        column_order = 1:32, #order columns as they are in the annot input
        column_names_gp = gpar(fontsize = 8), #adjust size of column labels
        #column_names_rot = 45, #adjust size and angle of column labels
        left_annotation = HA_Left) #put diagnosis annotation on the left


