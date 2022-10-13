
cd '/scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_04_PlotDenoisedCopyRatios/result'

set -o pipefail



    

  

gatk --java-options "-Xmx40g" PlotDenoisedCopyRatios \
    --denoised-copy-ratios /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_03_DenoiseReadCounts/result/P1809_AC_109.denoisedCR.tsv \
    --standardized-copy-ratios /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_03_DenoiseReadCounts/result/P1809_AC_109.standardizedCR.tsv \
    --sequence-dictionary /data/cqs/references/broad/hg38/v0/Homo_sapiens_assembly38.main.dict \
    --minimum-contig-length 46709983 \
    --output . \
    --output-prefix P1809_AC_109 
    

