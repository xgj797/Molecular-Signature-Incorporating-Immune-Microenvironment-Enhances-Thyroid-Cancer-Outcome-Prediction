
res_dir='/scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_bam_validation/result'
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

echo /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa/result/P6121_CP_31.sortedByCoord.bam      
if [[ ! -s /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa/result/P6121_CP_31.sortedByCoord.bam ]]; then
  echo file not exists: /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa/result/P6121_CP_31.sortedByCoord.bam
  touch $res_dir/P6121_CP_31.sortedByCoord.bam.not.exist
  rm P6121_CP_31.sortedByCoord.bam
  exit 1
fi

for i in {1..5}; do 
  echo iteration $i ...
  cp -f /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa/result/P6121_CP_31.sortedByCoord.bam P6121_CP_31.sortedByCoord.bam
  diff /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa/result/P6121_CP_31.sortedByCoord.bam P6121_CP_31.sortedByCoord.bam
  status=$?
  if [[ $status -eq 0 ]]; then
    break
  fi
done

if [[ $status -ne 0 ]]; then
  echo file copied is not idential to original file, quit.
  touch $res_dir/P6121_CP_31.sortedByCoord.bam.copy.failed
  rm P6121_CP_31.sortedByCoord.bam
  exit 1
fi

if [[ ! -s /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa/result/P6121_CP_31.sortedByCoord.bam.bai ]]; then
  echo file not exists: /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa/result/P6121_CP_31.sortedByCoord.bam.bai
  touch $res_dir/P6121_CP_31.sortedByCoord.bam.bai.not.exist
  rm P6121_CP_31.sortedByCoord.bam
  exit 1
fi

for i in {1..5}; do 
  echo iteration $i ...
  cp -f /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa/result/P6121_CP_31.sortedByCoord.bam.bai P6121_CP_31.sortedByCoord.bam.bai
  diff /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa/result/P6121_CP_31.sortedByCoord.bam.bai P6121_CP_31.sortedByCoord.bam.bai
  status=$?
  if [[ $status -eq 0 ]]; then
    break
  fi
done

if [[ $status -ne 0 ]]; then
  echo file copied is not idential to original file, quit.
  touch $res_dir/P6121_CP_31.sortedByCoord.bam.bai.copy.failed
  rm P6121_CP_31.sortedByCoord.bam
  exit 1
fi

ls *
echo localize end at `date`


    

 gatk ValidateSamFile -I P6121_CP_31.sortedByCoord.bam -O P6121_CP_31.txt --IGNORE_WARNINGS --SKIP_MATE_VALIDATION --VALIDATE_INDEX false --INDEX_VALIDATION_STRINGENCY NONE
    
status=$?
if [[ $status -ne 0 ]]; then
  touch P6121_CP_31.failed
fi

  


rm P6121_CP_31.sortedByCoord.bam P6121_CP_31.sortedByCoord.bam.bai

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
