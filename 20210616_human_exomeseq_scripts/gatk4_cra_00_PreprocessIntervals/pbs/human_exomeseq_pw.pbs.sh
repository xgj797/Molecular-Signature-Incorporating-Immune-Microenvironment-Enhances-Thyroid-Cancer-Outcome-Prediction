
cd '/scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_00_PreprocessIntervals/result'

set -o pipefail





  
gatk PreprocessIntervals \
    -L /data/cqs/references/exomeseq/IDT/Exome-IDT_V1.hg38.slop50.bed.interval_list \
    -R /data/cqs/references/broad/hg38/v0/Homo_sapiens_assembly38.fasta \
    --bin-length 0 \
    --padding 200 \
    --interval-merging-rule OVERLAPPING_ONLY \
    -O human_exomeseq.preprocessed.interval_list

#human_exomeseq__fileList1.list human_exomeseq.preprocessed.interval_list
  
