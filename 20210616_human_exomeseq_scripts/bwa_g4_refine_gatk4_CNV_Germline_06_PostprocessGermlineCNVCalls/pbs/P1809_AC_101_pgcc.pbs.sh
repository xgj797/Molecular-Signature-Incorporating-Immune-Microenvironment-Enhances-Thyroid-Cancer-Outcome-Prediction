
cd '/scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_06_PostprocessGermlineCNVCalls/result/P1809_AC_101'

set -o pipefail


  

cd /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_06_PostprocessGermlineCNVCalls/result/P1809_AC_101

gatk --java-options "-Xmx40G" PostprocessGermlineCNVCalls  \
   --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_001/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_002/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_003/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_004/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_005/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_006/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_007/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_008/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_009/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_010/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_011/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_012/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_013/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_014/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_015/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_016/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_017/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_018/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_019/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_020/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_021/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_022/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_023/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_024/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_025/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_026/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_027/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_028/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_029/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_030/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_031/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_032/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_033/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_034/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_035/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_036/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_037/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_038/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_039/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_040/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_041/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_042/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_043/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_044/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_045/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_046/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_047/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_048/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_049/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_050/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_051/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_052/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_053/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_054/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_055/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_056/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_057/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_058/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_059/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_060/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_061/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_062/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_063/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_064/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_065/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_066/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_067/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_068/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_069/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_070/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_071/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_072/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_073/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_074/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_075/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_076/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_077/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_078/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_079/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_080/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_081/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_082/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_083/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_084/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_085/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_086/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_087/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_088/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_089/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_090/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_091/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_092/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_093/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_094/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_095/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_096/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_097/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_098/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_099/gcc-calls --calls-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_100/gcc-calls \
   --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_001/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_002/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_003/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_004/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_005/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_006/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_007/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_008/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_009/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_010/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_011/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_012/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_013/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_014/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_015/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_016/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_017/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_018/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_019/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_020/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_021/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_022/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_023/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_024/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_025/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_026/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_027/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_028/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_029/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_030/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_031/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_032/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_033/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_034/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_035/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_036/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_037/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_038/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_039/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_040/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_041/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_042/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_043/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_044/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_045/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_046/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_047/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_048/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_049/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_050/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_051/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_052/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_053/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_054/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_055/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_056/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_057/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_058/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_059/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_060/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_061/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_062/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_063/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_064/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_065/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_066/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_067/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_068/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_069/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_070/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_071/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_072/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_073/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_074/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_075/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_076/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_077/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_078/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_079/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_080/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_081/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_082/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_083/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_084/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_085/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_086/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_087/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_088/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_089/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_090/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_091/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_092/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_093/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_094/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_095/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_096/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_097/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_098/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_099/gcc-model --model-shard-path /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_100/gcc-model \
  --sample-index 76 \
  --allosomal-contig chrX \
  --allosomal-contig chrY  \
  --autosomal-ref-copy-number 2 \
  --contig-ploidy-calls /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_04_DetermineGermlineContigPloidyCohortMode/result/human_exomeseq-calls \
  --output-denoised-copy-ratios P1809_AC_101.denoised_copy_ratios.tsv \
  --output-genotyped-intervals P1809_AC_101.genotyped_intervals.vcf.gz \
  --output-genotyped-segments P1809_AC_101.genotyped_segments.vcf.gz
            
rm -rf .cache .conda .config .theano

