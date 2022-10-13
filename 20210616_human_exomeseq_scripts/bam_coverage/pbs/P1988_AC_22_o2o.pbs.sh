
cd '/scratch/weissvl/shengq2/20210616_human_exomeseq/bam_coverage/result'

set -o pipefail



    

  
#mosdepth --by /data/cqs/references/exomeseq/IDT/Exome-IDT_V1.hg38.bed P1988_AC_22 /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine/result/P1988_AC_22.rmdup.recal.bam

rm P1988_AC_22.per-base.bed.gz P1988_AC_22.per-base.bed.gz.csi P1988_AC_22.mosdepth.global.dist.txt P1988_AC_22.mosdepth.region.dist.txt P1988_AC_22.regions.bed.gz.csi

zcat P1988_AC_22.regions.bed.gz | awk '{ cursize = $3-$2; total += $5 * cursize; size += cursize } END { print total/size }' > P1988_AC_22.coverage.txt
#zcat P1988_AC_22.regions.bed.gz | awk '{ total += $5; count++ } END { print total/count }' > P1988_AC_22.regions.coverage.mean.txt
  

