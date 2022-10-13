
cd '/scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_08_PlotModeledSegments/result'

set -o pipefail



    

  

gatk --java-options "-Xmx40g" PlotModeledSegments \
    --denoised-copy-ratios /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_03_DenoiseReadCounts/result/P3544_SB_01.denoisedCR.tsv \
    --segments /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_06_ModelSegments/result/P3544_SB_01/P3544_SB_01.modelFinal.seg \
    --sequence-dictionary /data/cqs/references/broad/hg38/v0/Homo_sapiens_assembly38.main.dict \
    --minimum-contig-length 46709983 \
    --output . \
    --output-prefix P3544_SB_01 
    

