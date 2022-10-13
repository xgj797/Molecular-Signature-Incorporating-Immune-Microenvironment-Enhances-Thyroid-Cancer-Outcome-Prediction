
cd '/scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_03_filterMAF/result'

set -o pipefail





python3 /data/cqs/softwares/ngsperl/lib/CQS/../Annotation/filterVcf.py -p 0.9 -f 0.3  -i /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_02_vqsr/result/human_exomeseq.indels.snp.recal.pass.norm.nospan.vcf.gz   -o human_exomeseq.maf_filtered.vcf.gz
