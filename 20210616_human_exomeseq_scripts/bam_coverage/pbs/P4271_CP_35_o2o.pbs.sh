
cd '/scratch/weissvl/shengq2/20210616_human_exomeseq/bam_coverage/result'

set -o pipefail



    

  
#mosdepth --by /data/cqs/references/exomeseq/IDT/Exome-IDT_V1.hg38.bed P4271_CP_35 /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine/result/P4271_CP_35.rmdup.recal.bam

rm P4271_CP_35.per-base.bed.gz P4271_CP_35.per-base.bed.gz.csi P4271_CP_35.mosdepth.global.dist.txt P4271_CP_35.mosdepth.region.dist.txt P4271_CP_35.regions.bed.gz.csi

zcat P4271_CP_35.regions.bed.gz | awk '{ cursize = $3-$2; total += $5 * cursize; size += cursize } END { print total/size }' > P4271_CP_35.coverage.txt
#zcat P4271_CP_35.regions.bed.gz | awk '{ total += $5; count++ } END { print total/count }' > P4271_CP_35.regions.coverage.mean.txt
  

