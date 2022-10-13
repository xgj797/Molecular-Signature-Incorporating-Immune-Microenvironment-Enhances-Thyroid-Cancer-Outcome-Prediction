
cd '/scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_06_ModelSegments/result/P1809_AC_065'

set -o pipefail



    

  

gatk --java-options "-Xmx40g" ModelSegments \
    --denoised-copy-ratios /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_03_DenoiseReadCounts/result/P1809_AC_065.denoisedCR.tsv \
    --output . \
    --output-prefix P1809_AC_065
    

