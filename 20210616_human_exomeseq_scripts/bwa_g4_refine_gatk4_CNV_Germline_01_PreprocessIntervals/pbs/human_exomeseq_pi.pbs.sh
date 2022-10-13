
cd '/scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_01_PreprocessIntervals/result'

set -o pipefail


  

cd /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_01_PreprocessIntervals/result

gatk --java-options "-Xmx20G" PreprocessIntervals  \
  -L /data/cqs/references/exomeseq/IDT/Exome-IDT_V1.hg38.slop50.bed  \
  --sequence-dictionary /data/cqs/references/broad/hg38/v0/Homo_sapiens_assembly38.dict \
  --reference /data/cqs/references/broad/hg38/v0/Homo_sapiens_assembly38.fasta \
  --interval-merging-rule OVERLAPPING_ONLY  \
  --output human_exomeseq.preprocessed.interval_list
  
rm -rf .conda

