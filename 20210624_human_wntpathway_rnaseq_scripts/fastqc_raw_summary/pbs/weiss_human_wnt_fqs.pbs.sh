
cd '/scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/fastqc_raw_summary/result'

set -o pipefail



#qcimg2pdf.sh -o weiss_human_wnt

python3 /data/cqs/softwares/ngsperl/lib/QC/fastQCSummary.py -i fileList1.txt -o weiss_human_wnt.FastQC

R --vanilla -f /data/cqs/softwares/ngsperl/lib/QC/fastQCSummary.r --args weiss_human_wnt.FastQC weiss_human_wnt.FastQC.Rmd
