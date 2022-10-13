
cd '/scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_05_CollectAllelicCounts/result'

set -o pipefail



    

  

gatk --java-options "-Xmx40g" CollectAllelicCounts \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine/result/P3544_SB_32.rmdup.recal.bam \
    -R /data/cqs/references/broad/hg38/v0/Homo_sapiens_assembly38.fasta \
    -O P3544_SB_32.allelicCounts.tsv \
    -L /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_00_PreprocessIntervals/result/human_exomeseq.preprocessed.interval_list 

