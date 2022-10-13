
cd '/scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/genetable/result'

set -o pipefail


cqstools data_table -k 0 -v 6 -e --fillMissingWithZero -o ./weiss_human_wnt.count -l /scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/genetable/pbs/weiss_human_wnt_tb.filelist -m /data/cqs/references/gencode/GRCh38.p13/gencode.v38.annotation.gtf.map 
