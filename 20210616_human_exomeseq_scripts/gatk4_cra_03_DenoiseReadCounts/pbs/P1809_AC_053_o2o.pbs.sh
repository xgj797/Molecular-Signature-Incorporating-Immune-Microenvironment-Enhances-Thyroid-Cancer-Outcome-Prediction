
cd '/scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_03_DenoiseReadCounts/result'

set -o pipefail



    

  

gatk --java-options "-Xmx40g" DenoiseReadCounts \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P1809_AC_053.counts.hdf5 \
    --standardized-copy-ratios P1809_AC_053.standardizedCR.tsv \
    --denoised-copy-ratios P1809_AC_053.denoisedCR.tsv \
    --count-panel-of-normals /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_02_CreateReadCountPanelOfNormals/result/human_exomeseq.pon.hdf5 

