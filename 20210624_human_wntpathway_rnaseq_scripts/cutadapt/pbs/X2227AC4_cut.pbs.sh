
res_dir='/scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/cutadapt/result'
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

echo /data/h_vivian_weiss/2227/2227-AC-4-TGACCAAT-TCTTTCCC_S12_R1_001.fastq.gz      
if [[ ! -s /data/h_vivian_weiss/2227/2227-AC-4-TGACCAAT-TCTTTCCC_S12_R1_001.fastq.gz ]]; then
  echo file not exists: /data/h_vivian_weiss/2227/2227-AC-4-TGACCAAT-TCTTTCCC_S12_R1_001.fastq.gz
  touch $res_dir/2227-AC-4-TGACCAAT-TCTTTCCC_S12_R1_001.fastq.gz.not.exist
  rm 2227-AC-4-TGACCAAT-TCTTTCCC_S12_R1_001.fastq.gz
  exit 1
fi

for i in {1..5}; do 
  echo iteration $i ...
  cp -fL /data/h_vivian_weiss/2227/2227-AC-4-TGACCAAT-TCTTTCCC_S12_R1_001.fastq.gz 2227-AC-4-TGACCAAT-TCTTTCCC_S12_R1_001.fastq.gz
  diff /data/h_vivian_weiss/2227/2227-AC-4-TGACCAAT-TCTTTCCC_S12_R1_001.fastq.gz 2227-AC-4-TGACCAAT-TCTTTCCC_S12_R1_001.fastq.gz
  status=$?
  if [[ $status -eq 0 ]]; then
    break
  fi
done

if [[ $status -ne 0 ]]; then
  echo file copied is not idential to original file, quit.
  touch $res_dir/2227-AC-4-TGACCAAT-TCTTTCCC_S12_R1_001.fastq.gz.copy.failed
  rm 2227-AC-4-TGACCAAT-TCTTTCCC_S12_R1_001.fastq.gz
  exit 1
fi

echo /data/h_vivian_weiss/2227/2227-AC-4-TGACCAAT-TCTTTCCC_S12_R2_001.fastq.gz      
if [[ ! -s /data/h_vivian_weiss/2227/2227-AC-4-TGACCAAT-TCTTTCCC_S12_R2_001.fastq.gz ]]; then
  echo file not exists: /data/h_vivian_weiss/2227/2227-AC-4-TGACCAAT-TCTTTCCC_S12_R2_001.fastq.gz
  touch $res_dir/2227-AC-4-TGACCAAT-TCTTTCCC_S12_R2_001.fastq.gz.not.exist
  rm 2227-AC-4-TGACCAAT-TCTTTCCC_S12_R1_001.fastq.gz 2227-AC-4-TGACCAAT-TCTTTCCC_S12_R2_001.fastq.gz
  exit 1
fi

for i in {1..5}; do 
  echo iteration $i ...
  cp -fL /data/h_vivian_weiss/2227/2227-AC-4-TGACCAAT-TCTTTCCC_S12_R2_001.fastq.gz 2227-AC-4-TGACCAAT-TCTTTCCC_S12_R2_001.fastq.gz
  diff /data/h_vivian_weiss/2227/2227-AC-4-TGACCAAT-TCTTTCCC_S12_R2_001.fastq.gz 2227-AC-4-TGACCAAT-TCTTTCCC_S12_R2_001.fastq.gz
  status=$?
  if [[ $status -eq 0 ]]; then
    break
  fi
done

if [[ $status -ne 0 ]]; then
  echo file copied is not idential to original file, quit.
  touch $res_dir/2227-AC-4-TGACCAAT-TCTTTCCC_S12_R2_001.fastq.gz.copy.failed
  rm 2227-AC-4-TGACCAAT-TCTTTCCC_S12_R1_001.fastq.gz 2227-AC-4-TGACCAAT-TCTTTCCC_S12_R2_001.fastq.gz
  exit 1
fi

ls *
echo localize end at `date`


cutadapt -j 2 -q 20 -a AGATCGGAAGAG -A AGATCGGAAGAG -n 3 --trim-n -a "A{50}" -a "T{50}" -a "G{50}" -a "C{50}" -A "A{50}" -A "T{50}" -A "G{50}" -A "C{50}" -o X2227AC4_clipped.1.fastq.gz -p X2227AC4_clipped.2.fastq.gz  -m 30  --too-short-output=X2227AC4_clipped.1.fastq.short.gz --too-short-paired-output=X2227AC4_clipped.2.fastq.short.gz 2227-AC-4-TGACCAAT-TCTTTCCC_S12_R1_001.fastq.gz 2227-AC-4-TGACCAAT-TCTTTCCC_S12_R2_001.fastq.gz
status=$?
if [[ $status -eq 0 ]]; then
  touch X2227AC4.succeed
  md5sum X2227AC4_clipped.1.fastq.gz > X2227AC4_clipped.1.fastq.gz.md5
  md5sum X2227AC4_clipped.2.fastq.gz > X2227AC4_clipped.2.fastq.gz.md5
else
  rm X2227AC4_clipped.1.fastq.gz X2227AC4_clipped.2.fastq.gz
  touch X2227AC4.failed
fi


rm 2227-AC-4-TGACCAAT-TCTTTCCC_S12_R1_001.fastq.gz 2227-AC-4-TGACCAAT-TCTTTCCC_S12_R2_001.fastq.gz

cutadapt --version 2>&1 | awk '{print "Cutadapt,v"$1}' > X2227AC4.version

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
