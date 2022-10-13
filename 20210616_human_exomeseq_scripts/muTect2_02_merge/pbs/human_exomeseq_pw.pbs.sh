
cd '/scratch/weissvl/shengq2/20210616_human_exomeseq/muTect2_02_merge/result'

set -o pipefail





python3 /data/cqs/softwares/ngsperl/lib/CQS/../GATK/mergeMutect.py   -i human_exomeseq__fileList1.list -o human_exomeseq_pass.combined.vcf
