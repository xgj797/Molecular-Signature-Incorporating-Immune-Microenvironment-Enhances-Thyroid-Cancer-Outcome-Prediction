
cd '/scratch/weissvl/shengq2/20210616_human_exomeseq/sequencetask/result'

set -o pipefail


R --vanilla --slave -f /scratch/weissvl/shengq2/20210616_human_exomeseq/sequencetask/result/human_exomeseq_summary.r 
R --vanilla --slave -f /scratch/weissvl/shengq2/20210616_human_exomeseq/sequencetask/result/human_exomeseq_report.r 
