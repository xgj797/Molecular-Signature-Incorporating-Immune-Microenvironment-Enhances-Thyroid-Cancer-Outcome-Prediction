
cd '/scratch/weissvl/shengq2/20210616_human_exomeseq/bam_coverage/result'

set -o pipefail



    

  
#mosdepth --by /data/cqs/references/exomeseq/IDT/Exome-IDT_V1.hg38.bed P1896_AC_27 /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine/result/P1896_AC_27.rmdup.recal.bam

rm P1896_AC_27.per-base.bed.gz P1896_AC_27.per-base.bed.gz.csi P1896_AC_27.mosdepth.global.dist.txt P1896_AC_27.mosdepth.region.dist.txt P1896_AC_27.regions.bed.gz.csi

zcat P1896_AC_27.regions.bed.gz | awk '{ cursize = $3-$2; total += $5 * cursize; size += cursize } END { print total/size }' > P1896_AC_27.coverage.txt
#zcat P1896_AC_27.regions.bed.gz | awk '{ total += $5; count++ } END { print total/count }' > P1896_AC_27.regions.coverage.mean.txt
  

