
cd '/scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/sequencetask/result'

set -o pipefail


R --vanilla --slave -f /scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/sequencetask/result/weiss_human_wnt_summary.r 
R --vanilla --slave -f /scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/sequencetask/result/weiss_human_wnt_report.r 
