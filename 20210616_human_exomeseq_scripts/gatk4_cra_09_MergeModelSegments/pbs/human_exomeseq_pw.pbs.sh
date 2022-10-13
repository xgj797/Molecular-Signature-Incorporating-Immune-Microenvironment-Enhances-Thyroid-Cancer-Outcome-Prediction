
cd '/scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_09_MergeModelSegments/result'

set -o pipefail





python3 ../GATK4/combineCRA.py  -i human_exomeseq__fileList1.list -o human_exomeseq.modelFinal.seg
    
