
cd '/scratch/weissvl/shengq2/20210616_human_exomeseq/bam_coverage/result'

set -o pipefail



    

  
#mosdepth --by /data/cqs/references/exomeseq/IDT/Exome-IDT_V1.hg38.bed P6121_CP_38 /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine/result/P6121_CP_38.rmdup.recal.bam

rm P6121_CP_38.per-base.bed.gz P6121_CP_38.per-base.bed.gz.csi P6121_CP_38.mosdepth.global.dist.txt P6121_CP_38.mosdepth.region.dist.txt P6121_CP_38.regions.bed.gz.csi

zcat P6121_CP_38.regions.bed.gz | awk '{ cursize = $3-$2; total += $5 * cursize; size += cursize } END { print total/size }' > P6121_CP_38.coverage.txt
#zcat P6121_CP_38.regions.bed.gz | awk '{ total += $5; count++ } END { print total/count }' > P6121_CP_38.regions.coverage.mean.txt
  

