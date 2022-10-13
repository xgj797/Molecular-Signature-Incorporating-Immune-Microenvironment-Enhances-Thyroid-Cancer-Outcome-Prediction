
cd '/scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result'

set -o pipefail



    

  

gatk --java-options "-Xmx40g" CollectReadCounts  \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine/result/P1809_AC_104.rmdup.recal.bam \
    --interval-merging-rule OVERLAPPING_ONLY \
    -O P1809_AC_104.counts.hdf5 \
    -L /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_00_PreprocessIntervals/result/human_exomeseq.preprocessed.interval_list 

