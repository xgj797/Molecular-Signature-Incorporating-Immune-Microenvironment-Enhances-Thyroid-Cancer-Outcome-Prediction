
cd '/scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_07_CallCopyRatioSegments/result'

set -o pipefail



    

  

gatk --java-options "-Xmx10g" CallCopyRatioSegments \
    --input /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_06_ModelSegments/result/P1809_AC_079/P1809_AC_079.cr.seg \
    --output P1809_AC_079.called.seg
    

