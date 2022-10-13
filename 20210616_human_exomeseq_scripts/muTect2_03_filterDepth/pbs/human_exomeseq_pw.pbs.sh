
cd '/scratch/weissvl/shengq2/20210616_human_exomeseq/muTect2_03_filterDepth/result'

set -o pipefail





python3 /data/cqs/softwares/ngsperl/lib/CQS/../GATK/filterMutect.py   -i /scratch/weissvl/shengq2/20210616_human_exomeseq/muTect2_02_merge/result/human_exomeseq_pass.combined.vcf   -o human_exomeseq.filtered.vcf
