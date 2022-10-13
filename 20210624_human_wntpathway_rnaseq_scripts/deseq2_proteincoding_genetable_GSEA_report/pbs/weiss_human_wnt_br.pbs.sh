
cd '/scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/deseq2_proteincoding_genetable_GSEA_report/result'

set -o pipefail



R -e "library(knitr);rmarkdown::render('weiss_human_wnt.Rmd');"
