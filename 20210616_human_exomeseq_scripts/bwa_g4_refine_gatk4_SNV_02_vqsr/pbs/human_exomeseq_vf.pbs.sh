
cd '/scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_02_vqsr/result'

set -o pipefail


  
if [ ! -s human_exomeseq.merged.vcf.gz ]; then
  echo CombineGVCFs=`date` 
  gatk --java-options "-Xmx40G" \
    CombineGVCFs \
    -R /data/cqs/references/broad/hg38/v0/bwa_index_0.7.17/Homo_sapiens_assembly38.fasta \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_003.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_004.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_005.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_006.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_008.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_010.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_011.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_012.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_013.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_014.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_015.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_016.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_017.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_018.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_019.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_020.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_021.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_022.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_023.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_024.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_025.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_026.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_027.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_028.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_029.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_030.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_031.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_032.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_033.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_034.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_035.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_036.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_037.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_038.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_039.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_040.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_041.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_042.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_043.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_044.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_045.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_046.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_047.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_048.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_049.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_050.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_051.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_052.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_053.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_054.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_055.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_056.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_057.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_058.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_059.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_060.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_061.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_062.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_063.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_064.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_065.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_066.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_067.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_068.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_069.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_070.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_071.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_072.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_073.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_074.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_075.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_076.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_077.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_078.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_079.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_080.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_081.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_082.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_083.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_084.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_085.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_086.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_087.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_088.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_089.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_090.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_091.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_092.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_093.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_094.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_095.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_096.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_097.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_098.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_099.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_100.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_101.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_102.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_103.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_104.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_105.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_106.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_107.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_108.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_109.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_110.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_111.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_112.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_113.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_114.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_115.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_116.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_117.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_118.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_119.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_120.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_121.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_122.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_123.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_124.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_125.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_126.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_127.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_128.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_129.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_130.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_131.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_132.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_133.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_134.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1809_AC_135.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1896_AC_01.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1896_AC_02.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1896_AC_03.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1896_AC_04.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1896_AC_05.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1896_AC_06.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1896_AC_07.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1896_AC_08.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1896_AC_09.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1896_AC_10.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1896_AC_11.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1896_AC_12.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1896_AC_13.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1896_AC_14.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1896_AC_15.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1896_AC_16.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1896_AC_17.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1896_AC_18.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1896_AC_19.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1896_AC_20.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1896_AC_21.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1896_AC_22.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1896_AC_23.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1896_AC_24.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1896_AC_25.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1896_AC_26.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1896_AC_27.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1896_AC_28.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1896_AC_29.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1896_AC_30.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1896_AC_31.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1896_AC_32.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1896_AC_33.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1896_AC_34.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1896_AC_35.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1896_AC_36.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1896_AC_37.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1896_AC_38.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1896_AC_39.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1988_AC_01.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1988_AC_02.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1988_AC_03.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1988_AC_04.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1988_AC_05.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1988_AC_06.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1988_AC_07.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1988_AC_08.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1988_AC_09.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1988_AC_10.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1988_AC_11.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1988_AC_12.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1988_AC_13.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1988_AC_14.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1988_AC_15.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1988_AC_16.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1988_AC_17.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1988_AC_18.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1988_AC_19.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1988_AC_20.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1988_AC_21.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1988_AC_22.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1988_AC_23.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1988_AC_24.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1988_AC_25.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P1988_AC_26.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P3544_SB_01.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P3544_SB_02.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P3544_SB_03.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P3544_SB_04.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P3544_SB_05.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P3544_SB_06.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P3544_SB_07.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P3544_SB_08.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P3544_SB_09.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P3544_SB_10.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P3544_SB_11.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P3544_SB_12.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P3544_SB_13.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P3544_SB_14.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P3544_SB_15.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P3544_SB_16.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P3544_SB_17.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P3544_SB_18.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P3544_SB_19.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P3544_SB_20.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P3544_SB_21.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P3544_SB_22.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P3544_SB_23.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P3544_SB_24.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P3544_SB_25.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P3544_SB_26.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P3544_SB_27.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P3544_SB_28.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P3544_SB_29.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P3544_SB_30.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P3544_SB_31.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P3544_SB_32.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P3544_SB_33.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P3544_SB_34.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P3544_SB_35.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4271_CP_01.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4271_CP_02.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4271_CP_03.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4271_CP_04.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4271_CP_05.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4271_CP_06.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4271_CP_07.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4271_CP_08.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4271_CP_09.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4271_CP_10.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4271_CP_11.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4271_CP_12.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4271_CP_13.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4271_CP_14.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4271_CP_15.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4271_CP_16.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4271_CP_17.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4271_CP_18.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4271_CP_19.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4271_CP_20.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4271_CP_21.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4271_CP_22.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4271_CP_23.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4271_CP_24.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4271_CP_25.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4271_CP_26.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4271_CP_27.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4271_CP_28.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4271_CP_29.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4271_CP_30.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4271_CP_32.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4271_CP_33.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4271_CP_34.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4271_CP_35.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4271_CP_36.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4271_CP_37.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4271_CP_38.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4271_CP_39.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_01.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_02.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_03.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_04.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_05.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_06.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_07.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_08.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_09.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_10.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_11.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_12.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_13.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_14.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_15.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_16.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_17.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_19.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_20.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_21.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_22.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_23.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_24.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_25.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_26.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_27.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_28.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_29.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_30.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_31.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_32.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_33.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_34.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_35.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_36.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_37.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_38.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_39.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_40.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_41.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_42.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_43.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_44.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_45.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_46.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_47.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_48.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_49.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_50.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_51.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_53.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_54.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_55.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_56.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_57.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_58.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_59.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_60.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_61.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_62.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_63.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_64.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_65.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_66.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_67.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_68.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_70.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_71.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_72.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P4601_CP_73.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P6121_CP_21.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P6121_CP_22.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P6121_CP_23.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P6121_CP_24.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P6121_CP_25.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P6121_CP_26.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P6121_CP_27.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P6121_CP_28.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P6121_CP_29.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P6121_CP_30.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P6121_CP_31.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P6121_CP_32.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P6121_CP_33.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P6121_CP_34.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P6121_CP_35.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P6121_CP_36.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P6121_CP_37.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P6121_CP_38.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P6121_CP_39.g.vcf.gz \
      -V /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result/P6121_CP_40.g.vcf.gz \
      -O human_exomeseq.merged.vcf.gz
fi

if [[ -s human_exomeseq.merged.vcf.gz && ! -s human_exomeseq.raw.vcf.gz ]]; then
  echo GenotypeGVCFs=`date` 
  gatk --java-options "-Xmx40G" \
    GenotypeGVCFs \
    -R /data/cqs/references/broad/hg38/v0/bwa_index_0.7.17/Homo_sapiens_assembly38.fasta \
    -D /data/cqs/references/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz \
    -O human_exomeseq.raw.vcf.gz \
    -V human_exomeseq.merged.vcf.gz 
fi

if [[ -s human_exomeseq.raw.vcf.gz && ! -s human_exomeseq.variant_filtered.vcf.gz ]]; then
  echo VariantFiltration=`date` 
  gatk --java-options "-Xmx40G" \
    VariantFiltration \
    --filter-expression "ExcessHet > 54.69" \
    --filter-name ExcessHet \
    -O human_exomeseq.variant_filtered.vcf.gz \
    -V human_exomeseq.raw.vcf.gz
fi

if [[ -s human_exomeseq.variant_filtered.vcf.gz && ! -s human_exomeseq.variant_filtered.sites_only.vcf.gz ]]; then
  echo MakeSitesOnlyVcf=`date` 
  gatk --java-options "-Xmx40G" \
    MakeSitesOnlyVcf \
    --INPUT human_exomeseq.variant_filtered.vcf.gz \
    --OUTPUT human_exomeseq.variant_filtered.sites_only.vcf.gz
fi

if [[ -s human_exomeseq.variant_filtered.sites_only.vcf.gz && ! -s human_exomeseq.indels.recal.vcf.gz ]]; then
  echo IndelVariantRecalibrator=`date`
  gatk --java-options "-Xmx40G" \
    VariantRecalibrator \
    -V human_exomeseq.variant_filtered.sites_only.vcf.gz \
    -O human_exomeseq.indels.recal.vcf.gz \
    --tranches-file human_exomeseq.indels.tranches \
    --trust-all-polymorphic \
    -tranche 100.0 \
    -tranche 99.95 \
    -tranche 99.9 \
    -tranche 99.5 \
    -tranche 99.0 \
    -tranche 97.0 \
    -tranche 96.0 \
    -tranche 95.0 \
    -tranche 94.0 \
    -tranche 93.5 \
    -tranche 93.0 \
    -tranche 92.0 \
    -tranche 91.0 \
    -tranche 90.0 \
    -an QD \
    -an MQRankSum \
    -an ReadPosRankSum \
    -an FS \
    -an MQ \
    -an SOR \
    -an DP \
    --max-gaussians 4 \
    --resource:mills,known=false,training=true,truth=true,prior=12 /data/cqs/references/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
     --resource:axiomPoly,known=false,training=true,truth=false,prior=10 /data/cqs/references/broad/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz \
     --resource:dbsnp,known=true,training=false,truth=false,prior=2 /data/cqs/references/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz \
     -mode INDEL 
fi

if [[ -s human_exomeseq.variant_filtered.sites_only.vcf.gz && ! -s human_exomeseq.snp.recal.vcf.gz ]]; then
  echo SnpVariantRecalibrator=`date`
  gatk --java-options "-Xmx40G" \
    VariantRecalibrator \
    -V human_exomeseq.variant_filtered.sites_only.vcf.gz \
    -O human_exomeseq.snp.recal.vcf.gz \
    --tranches-file human_exomeseq.snp.tranches \
    --trust-all-polymorphic \
    -tranche 100.0 \
    -tranche 99.95 \
    -tranche 99.9 \
    -tranche 99.8 \
    -tranche 99.6 \
    -tranche 99.5 \
    -tranche 99.4 \
    -tranche 99.3 \
    -tranche 99.0 \
    -tranche 98.0 \
    -tranche 97.0 \
    -tranche 90.0 \
    -an QD \
    -an MQRankSum \
    -an ReadPosRankSum \
    -an FS \
    -an MQ \
    -an SOR \
    -an DP \
    --max-gaussians 6 \
    --resource:hapmap,known=false,training=true,truth=true,prior=15 /data/cqs/references/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz \
     --resource:omni,known=false,training=true,truth=true,prior=12 /data/cqs/references/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz \
     --resource:1000G,known=false,training=true,truth=false,prior=10 /data/cqs/references/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
     --resource:dbsnp,known=true,training=false,truth=false,prior=2 /data/cqs/references/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz \
     -mode SNP
fi
  

if [[ ! -s human_exomeseq.indels.recal.tmp.vcf.gz ]]; then
  echo IndelApplyVQSR=`date`
  gatk --java-options "-Xmx40G" \
    ApplyVQSR \
    -O human_exomeseq.indels.recal.tmp.vcf.gz \
    -V human_exomeseq.variant_filtered.vcf.gz \
    --recal-file human_exomeseq.indels.recal.vcf.gz \
    --tranches-file human_exomeseq.indels.tranches \
    --truth-sensitivity-filter-level 99.7 \
    --create-output-variant-index true \
    -mode INDEL
fi

if [[ -s human_exomeseq.snp.recal.vcf.gz && ! -s human_exomeseq.indels.snp.recal.vcf.gz ]]; then
  echo SnpApplyVQSR=`date`
  gatk --java-options "-Xmx40G" \
    ApplyVQSR \
    -O human_exomeseq.indels.snp.recal.vcf.gz \
    -V human_exomeseq.indels.recal.tmp.vcf.gz \
    --recal-file human_exomeseq.snp.recal.vcf.gz \
    --tranches-file human_exomeseq.snp.tranches \
    --truth-sensitivity-filter-level 99.7 \
    --create-output-variant-index true \
    -mode SNP
fi

if [[ -s human_exomeseq.indels.snp.recal.vcf.gz && ! -s human_exomeseq.indels.snp.recal.pass.vcf.gz ]]; then
  echo SelectVariant=`date`
  gatk --java-options "-Xmx40G" \
    SelectVariants \
    -O human_exomeseq.indels.snp.recal.pass.vcf.gz \
    -V human_exomeseq.indels.snp.recal.vcf.gz \
    --exclude-filtered
fi


if [[ -s human_exomeseq.indels.snp.recal.pass.vcf.gz && ! -s human_exomeseq.norm.vcf ]]; then
  echo LeftAlignAndNorm=`date`
  bcftools norm -m- -o human_exomeseq.split.vcf human_exomeseq.indels.snp.recal.pass.vcf.gz 
  bcftools norm -f /data/cqs/references/broad/hg38/v0/bwa_index_0.7.17/Homo_sapiens_assembly38.fasta -o human_exomeseq.norm.vcf human_exomeseq.split.vcf 
fi

if [[ -s human_exomeseq.norm.vcf && ! -s human_exomeseq.indels.snp.recal.pass.norm.nospan.vcf.gz ]]; then
  echo noSpanDeletion=`date`
  python3 /data/cqs/softwares/ngsperl/lib/GATK4/fixLeftTrimDeletion.py -i human_exomeseq.norm.vcf -o human_exomeseq.indels.snp.recal.pass.norm.nospan.vcf
  bgzip human_exomeseq.indels.snp.recal.pass.norm.nospan.vcf
  tabix -p vcf human_exomeseq.indels.snp.recal.pass.norm.nospan.vcf.gz
fi

if [[ -s human_exomeseq.indels.snp.recal.pass.norm.nospan.vcf.gz ]]; then
  rm human_exomeseq.merged.vcf.gz human_exomeseq.merged.vcf.gz.tbi \
    human_exomeseq.raw.vcf.gz human_exomeseq.raw.vcf.gz.tbi \
    human_exomeseq.indels.recal.vcf.gz human_exomeseq.indels.recal.vcf.gz.tbi \
    human_exomeseq.snp.recal.vcf.gz human_exomeseq.snp.recal.vcf.gz.tbi \
    human_exomeseq.indels.snp.recal.vcf.gz human_exomeseq.indels.snp.recal.vcf.gz.tbi \
    human_exomeseq.indels.recal.tmp.vcf.gz human_exomeseq.indels.recal.tmp.vcf.gz.tbi \
    human_exomeseq.indels.recal.vcf.gz human_exomeseq.indels.recal.vcf.gz.tbi human_exomeseq.indels.tranches  \
    human_exomeseq.snp.recal.vcf.gz human_exomeseq.snp.recal.vcf.gz.tbi human_exomeseq.snp.tranches \
    human_exomeseq.split.vcf human_exomeseq.norm.vcf \
    .conda
fi

