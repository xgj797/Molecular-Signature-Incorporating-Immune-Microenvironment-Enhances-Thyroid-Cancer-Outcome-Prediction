
cd '/scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_07_CallCopyRatioSegments/result'

set -o pipefail



    

  

gatk --java-options "-Xmx10g" CallCopyRatioSegments \
    --input /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_06_ModelSegments/result/P1896_AC_29/P1896_AC_29.cr.seg \
    --output P1896_AC_29.called.seg
    

