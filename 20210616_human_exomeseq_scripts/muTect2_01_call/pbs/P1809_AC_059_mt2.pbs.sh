
cd '/scratch/weissvl/shengq2/20210616_human_exomeseq/muTect2_01_call/result/P1809_AC_059'

set -o pipefail


 

mkdir tmp_P1809_AC_059    

if [ ! -s /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine/result/P1809_AC_059.rmdup.recal.bam.bai ]; then
  samtools index /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine/result/P1809_AC_059.rmdup.recal.bam
fi

echo Mutect2 ...
gatk --java-options "-Djava.io.tmpdir=`pwd`/tmp_P1809_AC_059 -Xms40g" Mutect2 --downsampling-stride 20 --max-reads-per-alignment-start 6 --max-suspicious-reads-per-alignment-start 6 --germline-resource /data/cqs/references/broad/hg38/af-only-gnomad.hg38.vcf.gz  --panel-of-normals /data/h_vangard_1/references/broad/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz  -L /data/cqs/references/exomeseq/IDT/Exome-IDT_V1.hg38.slop50.bed \
  -R /data/cqs/references/broad/hg38/v0/Homo_sapiens_assembly38.fasta \
  -I /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine/result/P1809_AC_059.rmdup.recal.bam -tumor P1809_AC_059 --f1r2-tar-gz f1r2.tar.gz \
  -O P1809_AC_059.unfiltered.vcf.gz

m2_exit_code=$?
if [[ $m2_exit_code -eq 0 ]]; then

  echo LearnReadOrientationModel ...
  gatk --java-options "-Djava.io.tmpdir=`pwd`/tmp_P1809_AC_059 -Xms40g" LearnReadOrientationModel -I f1r2.tar.gz -O read-orientation-model.tar.gz

  echo GetPileupSummaries tumor ...
  gatk --java-options "-Djava.io.tmpdir=`pwd`/tmp_P1809_AC_059 -Xms40g" GetPileupSummaries \
    -R /data/cqs/references/broad/hg38/v0/Homo_sapiens_assembly38.fasta --interval-set-rule INTERSECTION -L /data/cqs/references/exomeseq/IDT/Exome-IDT_V1.hg38.slop50.bed \
    -I /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine/result/P1809_AC_059.rmdup.recal.bam \
    -V /data/h_vangard_1/references/broad/gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz \
    -L /data/h_vangard_1/references/broad/gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz \
    -O P1809_AC_059.tumor_pileups.table 

  echo CalculateContamination ...
  gatk --java-options "-Djava.io.tmpdir=`pwd`/tmp_P1809_AC_059 -Xms40g" CalculateContamination \
    -I P1809_AC_059.tumor_pileups.table \
    -O P1809_AC_059.contamination.table \
    --tumor-segmentation P1809_AC_059.segments.table 

  echo FilterMutectCalls ...
  gatk --java-options "-Djava.io.tmpdir=`pwd`/tmp_P1809_AC_059 -Xms40g" FilterMutectCalls \
    -V P1809_AC_059.unfiltered.vcf.gz \
    -R /data/cqs/references/broad/hg38/v0/Homo_sapiens_assembly38.fasta \
    -O P1809_AC_059.filtered.vcf.gz --contamination-table P1809_AC_059.contamination.table --tumor-segmentation P1809_AC_059.segments.table  --ob-priors read-orientation-model.tar.gz \
    --stats P1809_AC_059.unfiltered.vcf.gz.stats \
    --filtering-stats P1809_AC_059.filtered.vcf.gz.stats

  filter_exit_code=$?
  if [[ $filter_exit_code -ne 0 ]]; then
    rm P1809_AC_059.filtered.vcf.gz.*
    touch P1809_AC_059.filtered.vcf.gz.failed
  else
    rm P1809_AC_059.unfiltered.vcf.gz P1809_AC_059.unfiltered.vcf.gz.tbi
    echo SelectVariants ...
    gatk --java-options "-Djava.io.tmpdir=`pwd`/tmp_P1809_AC_059 -Xms40g" SelectVariants \
      --exclude-filtered \
      -V P1809_AC_059.filtered.vcf.gz \
      -O P1809_AC_059.pass.vcf.gz

    select_exit_code=$?
    if [[ $select_exit_code -ne 0 ]]; then
      rm P1809_AC_059.pass.vcf.gz.*
      touch P1809_AC_059.pass.vcf.gz.failed
    else
      touch P1809_AC_059.pass.vcf.gz.succeed
    fi
  fi
else
  rm P1809_AC_059.unfiltered.vcf.gz.*
  touch P1809_AC_059.unfiltered.vcf.gz.failed
fi

rm -rf tmp_P1809_AC_059 f1r2.tar.gz read-orientation-model.tar.gz

