
cd '/scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_07_CombineGCNV/result'

set -o pipefail





python3 /data/cqs/softwares/ngsperl/lib/CQS/../GATK4/combineGCNV.py   -b /data/cqs/references/exomeseq/IDT/Exome-IDT_V1.hg38.slop50.bed --annovar_db /data/cqs/references/annovar/humandb/ --annovar_buildver hg38 -i human_exomeseq__fileList1.list -o human_exomeseq.txt
