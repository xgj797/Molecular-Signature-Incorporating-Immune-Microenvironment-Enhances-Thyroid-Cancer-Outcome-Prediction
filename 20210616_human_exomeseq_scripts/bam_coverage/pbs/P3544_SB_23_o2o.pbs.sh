
cd '/scratch/weissvl/shengq2/20210616_human_exomeseq/bam_coverage/result'

set -o pipefail



    

  
#mosdepth --by /data/cqs/references/exomeseq/IDT/Exome-IDT_V1.hg38.bed P3544_SB_23 /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine/result/P3544_SB_23.rmdup.recal.bam

rm P3544_SB_23.per-base.bed.gz P3544_SB_23.per-base.bed.gz.csi P3544_SB_23.mosdepth.global.dist.txt P3544_SB_23.mosdepth.region.dist.txt P3544_SB_23.regions.bed.gz.csi

zcat P3544_SB_23.regions.bed.gz | awk '{ cursize = $3-$2; total += $5 * cursize; size += cursize } END { print total/size }' > P3544_SB_23.coverage.txt
#zcat P3544_SB_23.regions.bed.gz | awk '{ total += $5; count++ } END { print total/count }' > P3544_SB_23.regions.coverage.mean.txt
  

