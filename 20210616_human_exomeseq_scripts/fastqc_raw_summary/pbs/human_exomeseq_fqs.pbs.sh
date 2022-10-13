
cd '/scratch/weissvl/shengq2/20210616_human_exomeseq/fastqc_raw_summary/result'

set -o pipefail



#qcimg2pdf.sh -o human_exomeseq

python3 /data/cqs/softwares/ngsperl/lib/QC/fastQCSummary.py -i fileList1.txt -o human_exomeseq.FastQC

R --vanilla -f /data/cqs/softwares/ngsperl/lib/QC/fastQCSummary.r --args human_exomeseq.FastQC human_exomeseq.FastQC.Rmd
