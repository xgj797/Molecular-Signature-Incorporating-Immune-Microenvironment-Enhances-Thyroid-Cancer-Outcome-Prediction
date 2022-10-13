
res_dir='/scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine/result'
tmp_dir=$(mktemp -d -t ci-$(date +%Y-%m-%d-%H-%M-%S)-XXXXXXXXXX)

tmp_cleaner()
{
rm -rf ${tmp_dir}
exit -1
}
trap 'tmp_cleaner' TERM

echo using tmp_dir=$tmp_dir
cd $tmp_dir


set -o pipefail



echo localize start at `date`

echo /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa/result/P4601_CP_05.sortedByCoord.bam      
if [[ ! -s /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa/result/P4601_CP_05.sortedByCoord.bam ]]; then
  echo file not exists: /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa/result/P4601_CP_05.sortedByCoord.bam
  touch $res_dir/P4601_CP_05.sortedByCoord.bam.not.exist
  rm P4601_CP_05.sortedByCoord.bam
  exit 1
fi

for i in {1..5}; do 
  echo iteration $i ...
  cp -f /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa/result/P4601_CP_05.sortedByCoord.bam P4601_CP_05.sortedByCoord.bam
  diff /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa/result/P4601_CP_05.sortedByCoord.bam P4601_CP_05.sortedByCoord.bam
  status=$?
  if [[ $status -eq 0 ]]; then
    break
  fi
done

if [[ $status -ne 0 ]]; then
  echo file copied is not idential to original file, quit.
  touch $res_dir/P4601_CP_05.sortedByCoord.bam.copy.failed
  rm P4601_CP_05.sortedByCoord.bam
  exit 1
fi

if [[ ! -s /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa/result/P4601_CP_05.sortedByCoord.bam.bai ]]; then
  echo file not exists: /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa/result/P4601_CP_05.sortedByCoord.bam.bai
  touch $res_dir/P4601_CP_05.sortedByCoord.bam.bai.not.exist
  rm P4601_CP_05.sortedByCoord.bam
  exit 1
fi

for i in {1..5}; do 
  echo iteration $i ...
  cp -f /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa/result/P4601_CP_05.sortedByCoord.bam.bai P4601_CP_05.sortedByCoord.bam.bai
  diff /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa/result/P4601_CP_05.sortedByCoord.bam.bai P4601_CP_05.sortedByCoord.bam.bai
  status=$?
  if [[ $status -eq 0 ]]; then
    break
  fi
done

if [[ $status -ne 0 ]]; then
  echo file copied is not idential to original file, quit.
  touch $res_dir/P4601_CP_05.sortedByCoord.bam.bai.copy.failed
  rm P4601_CP_05.sortedByCoord.bam
  exit 1
fi

ls *
echo localize end at `date`



echo MarkDuplicates=`date` 
gatk --java-options "-Xmx40G" \
  MarkDuplicates \
  --INPUT P4601_CP_05.sortedByCoord.bam \
  --OUTPUT P4601_CP_05.rmdup.bam \
  --METRICS_FILE P4601_CP_05.rmdup.bam.metrics \
  --VALIDATION_STRINGENCY SILENT \
  --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
  --REMOVE_DUPLICATES true \
  --CREATE_INDEX true \
  --ASSUME_SORTED true \
  --CREATE_MD5_FILE false

status=$?
if [[ $status -ne 0 ]]; then
  touch P4601_CP_05.MarkDuplicates.failed
  rm -f P4601_CP_05.rmdup.bam P4601_CP_05.rmdup.bam.metrics P4601_CP_05.rmdup.bai
else
  echo BaseRecalibrator=`date` 
  gatk --java-options "-Xmx40G" \
    BaseRecalibrator \
    -R /data/cqs/references/broad/hg38/v0/bwa_index_0.7.17/Homo_sapiens_assembly38.fasta \
    -I P4601_CP_05.rmdup.bam \
    --use-original-qualities  \
     --known-sites /data/cqs/references/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    -O P4601_CP_05.rmdup.recal_data.csv  

  status=$?
  if [[ $status -ne 0 ]]; then
    touch P4601_CP_05.BaseRecalibrator.failed
  else
    echo ApplyBQSR=`date`
    gatk --java-options "-Xmx40G" \
      ApplyBQSR \
      -R /data/cqs/references/broad/hg38/v0/bwa_index_0.7.17/Homo_sapiens_assembly38.fasta \
      -I P4601_CP_05.rmdup.bam \
      -O P4601_CP_05.rmdup.recal.bam \
      -bqsr P4601_CP_05.rmdup.recal_data.csv \
      --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
      --add-output-sam-program-record \
      --create-output-bam-md5 \
      --use-original-qualities

    status=$?
    if [[ $status -ne 0 ]]; then
      touch P4601_CP_05.ApplyBQSR.failed
    else
      touch P4601_CP_05.succeed
      ln P4601_CP_05.rmdup.recal.bai P4601_CP_05.rmdup.recal.bam.bai

      echo flagstat = `date` 
      samtools flagstat P4601_CP_05.rmdup.recal.bam > P4601_CP_05.rmdup.recal.bam.stat 

      echo flagstat = `date` 
      samtools idxstats P4601_CP_05.rmdup.recal.bam > P4601_CP_05.rmdup.recal.bam.chromosome.count 

      rm /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine/result/P4601_CP_05.*.failed
      rm  P4601_CP_05.rmdup.bam P4601_CP_05.rmdup.bai
    fi
  fi
fi


rm P4601_CP_05.sortedByCoord.bam P4601_CP_05.sortedByCoord.bam.bai

if [[ -d $tmp_dir && $tmp_dir != '/' ]]; then
  echo copy result from $tmp_dir to $res_dir
  #if the pbs was generated again during task is running, copy may be unpredictable. 
  #make sure to change to tmp_dir before copy result

  cd $tmp_dir
  cd $tmp_dir
  cd $tmp_dir
  cd $tmp_dir
  cd $tmp_dir
  cd $tmp_dir
  cd $tmp_dir
  cd $tmp_dir
  cd $tmp_dir
  cd $tmp_dir

  cp -p -r * $res_dir
  cd $res_dir
  echo delete tmp folder $tmp_dir
  rm -rf $tmp_dir
  echo move file and clean tmp folder done.
fi
