
res_dir='/scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_bam_validation/result'
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

echo /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine/result/P1809_AC_079.rmdup.recal.bam      
if [[ ! -s /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine/result/P1809_AC_079.rmdup.recal.bam ]]; then
  echo file not exists: /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine/result/P1809_AC_079.rmdup.recal.bam
  touch $res_dir/P1809_AC_079.rmdup.recal.bam.not.exist
  rm P1809_AC_079.rmdup.recal.bam
  exit 1
fi

for i in {1..5}; do 
  echo iteration $i ...
  cp -f /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine/result/P1809_AC_079.rmdup.recal.bam P1809_AC_079.rmdup.recal.bam
  diff /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine/result/P1809_AC_079.rmdup.recal.bam P1809_AC_079.rmdup.recal.bam
  status=$?
  if [[ $status -eq 0 ]]; then
    break
  fi
done

if [[ $status -ne 0 ]]; then
  echo file copied is not idential to original file, quit.
  touch $res_dir/P1809_AC_079.rmdup.recal.bam.copy.failed
  rm P1809_AC_079.rmdup.recal.bam
  exit 1
fi

if [[ ! -s /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine/result/P1809_AC_079.rmdup.recal.bam.bai ]]; then
  echo file not exists: /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine/result/P1809_AC_079.rmdup.recal.bam.bai
  touch $res_dir/P1809_AC_079.rmdup.recal.bam.bai.not.exist
  rm P1809_AC_079.rmdup.recal.bam
  exit 1
fi

for i in {1..5}; do 
  echo iteration $i ...
  cp -f /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine/result/P1809_AC_079.rmdup.recal.bam.bai P1809_AC_079.rmdup.recal.bam.bai
  diff /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine/result/P1809_AC_079.rmdup.recal.bam.bai P1809_AC_079.rmdup.recal.bam.bai
  status=$?
  if [[ $status -eq 0 ]]; then
    break
  fi
done

if [[ $status -ne 0 ]]; then
  echo file copied is not idential to original file, quit.
  touch $res_dir/P1809_AC_079.rmdup.recal.bam.bai.copy.failed
  rm P1809_AC_079.rmdup.recal.bam
  exit 1
fi

ls *
echo localize end at `date`


    

  
gatk ValidateSamFile -I P1809_AC_079.rmdup.recal.bam -O P1809_AC_079.txt --IGNORE_WARNINGS --SKIP_MATE_VALIDATION --VALIDATE_INDEX false --INDEX_VALIDATION_STRINGENCY NONE
    
status=$?
if [[ $status -ne 0 ]]; then
  if [[ -e P1809_AC_079.txt ]]; then
    mv P1809_AC_079.txt P1809_AC_079.failed
  else
    touch P1809_AC_079.failed
  fi
fi

  


rm P1809_AC_079.rmdup.recal.bam P1809_AC_079.rmdup.recal.bam.bai

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
