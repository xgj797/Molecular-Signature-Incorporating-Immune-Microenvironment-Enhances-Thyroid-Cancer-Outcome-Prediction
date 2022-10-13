
cd '/scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_08_SizeFactor/result'

set -o pipefail





python3 /data/cqs/softwares/ngsperl/lib/CQS/../GATK4/getBackgroundCount.py   -i /data/cqs/references/exomeseq/IDT/Exome-IDT_V1.hg38.slop50.bed -c /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_07_CombineGCNV/result/human_exomeseq.txt -b human_exomeseq__fileList1.list -o human_exomeseq.txt
