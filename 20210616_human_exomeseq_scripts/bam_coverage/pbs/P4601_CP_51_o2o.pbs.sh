
cd '/scratch/weissvl/shengq2/20210616_human_exomeseq/bam_coverage/result'

set -o pipefail



    

  
#mosdepth --by /data/cqs/references/exomeseq/IDT/Exome-IDT_V1.hg38.bed P4601_CP_51 /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine/result/P4601_CP_51.rmdup.recal.bam

rm P4601_CP_51.per-base.bed.gz P4601_CP_51.per-base.bed.gz.csi P4601_CP_51.mosdepth.global.dist.txt P4601_CP_51.mosdepth.region.dist.txt P4601_CP_51.regions.bed.gz.csi

zcat P4601_CP_51.regions.bed.gz | awk '{ cursize = $3-$2; total += $5 * cursize; size += cursize } END { print total/size }' > P4601_CP_51.coverage.txt
#zcat P4601_CP_51.regions.bed.gz | awk '{ total += $5; count++ } END { print total/count }' > P4601_CP_51.regions.coverage.mean.txt
  

