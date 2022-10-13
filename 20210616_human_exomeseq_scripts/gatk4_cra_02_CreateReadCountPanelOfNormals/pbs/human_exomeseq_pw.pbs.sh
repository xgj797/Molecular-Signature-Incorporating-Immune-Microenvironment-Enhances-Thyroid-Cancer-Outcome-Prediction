
cd '/scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_02_CreateReadCountPanelOfNormals/result'

set -o pipefail





  

gatk --java-options "-Xmx40g" CreateReadCountPanelOfNormals \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P1809_AC_031.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P1809_AC_032.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P1809_AC_033.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P1809_AC_034.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P1809_AC_035.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P1809_AC_036.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P1809_AC_037.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P1809_AC_038.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P1809_AC_039.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P1809_AC_040.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P1809_AC_041.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P1809_AC_042.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P1809_AC_043.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P1809_AC_044.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P1809_AC_045.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P1809_AC_046.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P1809_AC_047.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P1809_AC_048.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P1809_AC_049.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P1809_AC_050.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P1809_AC_122.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P1809_AC_123.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P1809_AC_124.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P1809_AC_125.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P1809_AC_126.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P1809_AC_127.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P1809_AC_128.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P1809_AC_129.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P1809_AC_130.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P1809_AC_131.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P1809_AC_132.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P1809_AC_133.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P1809_AC_134.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P1809_AC_135.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P4601_CP_01.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P4601_CP_02.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P4601_CP_03.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P4601_CP_04.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P4601_CP_05.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P4601_CP_06.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P4601_CP_07.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P4601_CP_08.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P4601_CP_09.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P4601_CP_10.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P4601_CP_11.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P4601_CP_12.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P4601_CP_13.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P4601_CP_14.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P4601_CP_15.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P4601_CP_16.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P4601_CP_17.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P4601_CP_19.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P4601_CP_20.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P4601_CP_21.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P4601_CP_22.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P4601_CP_23.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P4601_CP_24.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P4601_CP_25.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P4601_CP_26.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P4601_CP_27.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P4601_CP_28.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P4601_CP_29.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P4601_CP_30.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P4601_CP_31.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P4601_CP_32.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P4601_CP_33.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P4601_CP_34.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P4601_CP_35.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P4601_CP_36.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P4601_CP_37.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P4601_CP_38.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P4601_CP_39.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P4601_CP_40.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P4601_CP_41.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P4601_CP_42.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P4601_CP_43.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P4601_CP_44.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P4601_CP_45.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P4601_CP_46.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P4601_CP_47.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P4601_CP_48.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P4601_CP_49.counts.hdf5 \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/gatk4_cra_01_CollectReadCounts/result/P4601_CP_50.counts.hdf5 \
    --minimum-interval-median-percentile 5.0 \
    -O human_exomeseq.pon.hdf5
        
