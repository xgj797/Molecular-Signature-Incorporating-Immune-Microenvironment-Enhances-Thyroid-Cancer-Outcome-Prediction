
cd '/scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_1_scatterIntervals/result'

set -o pipefail



python3 /data/cqs/softwares/ngsperl/lib/CQS/../GATK4/scatterIntervals.py  -n 100  -i /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_03_FilterIntervals/result/human_exomeseq.filtered.interval_list   -o human_exomeseq

