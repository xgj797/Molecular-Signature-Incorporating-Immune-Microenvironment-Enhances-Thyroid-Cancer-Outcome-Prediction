
res_dir='/scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/fastqc_raw/result/X2021AC5'
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

echo /data/h_vivian_weiss/2021/2021-AC-5-ACAGTGAT-TCTTTCCC_S124_R1_001.fastq.gz      
if [[ ! -s /data/h_vivian_weiss/2021/2021-AC-5-ACAGTGAT-TCTTTCCC_S124_R1_001.fastq.gz ]]; then
  echo file not exists: /data/h_vivian_weiss/2021/2021-AC-5-ACAGTGAT-TCTTTCCC_S124_R1_001.fastq.gz
  touch $res_dir/2021-AC-5-ACAGTGAT-TCTTTCCC_S124_R1_001.fastq.gz.not.exist
  rm 2021-AC-5-ACAGTGAT-TCTTTCCC_S124_R1_001.fastq.gz
  exit 1
fi

for i in {1..5}; do 
  echo iteration $i ...
  cp -f /data/h_vivian_weiss/2021/2021-AC-5-ACAGTGAT-TCTTTCCC_S124_R1_001.fastq.gz 2021-AC-5-ACAGTGAT-TCTTTCCC_S124_R1_001.fastq.gz
  diff /data/h_vivian_weiss/2021/2021-AC-5-ACAGTGAT-TCTTTCCC_S124_R1_001.fastq.gz 2021-AC-5-ACAGTGAT-TCTTTCCC_S124_R1_001.fastq.gz
  status=$?
  if [[ $status -eq 0 ]]; then
    break
  fi
done

if [[ $status -ne 0 ]]; then
  echo file copied is not idential to original file, quit.
  touch $res_dir/2021-AC-5-ACAGTGAT-TCTTTCCC_S124_R1_001.fastq.gz.copy.failed
  rm 2021-AC-5-ACAGTGAT-TCTTTCCC_S124_R1_001.fastq.gz
  exit 1
fi

echo /data/h_vivian_weiss/2021/2021-AC-5-ACAGTGAT-TCTTTCCC_S124_R2_001.fastq.gz      
if [[ ! -s /data/h_vivian_weiss/2021/2021-AC-5-ACAGTGAT-TCTTTCCC_S124_R2_001.fastq.gz ]]; then
  echo file not exists: /data/h_vivian_weiss/2021/2021-AC-5-ACAGTGAT-TCTTTCCC_S124_R2_001.fastq.gz
  touch $res_dir/2021-AC-5-ACAGTGAT-TCTTTCCC_S124_R2_001.fastq.gz.not.exist
  rm 2021-AC-5-ACAGTGAT-TCTTTCCC_S124_R1_001.fastq.gz 2021-AC-5-ACAGTGAT-TCTTTCCC_S124_R2_001.fastq.gz
  exit 1
fi

for i in {1..5}; do 
  echo iteration $i ...
  cp -f /data/h_vivian_weiss/2021/2021-AC-5-ACAGTGAT-TCTTTCCC_S124_R2_001.fastq.gz 2021-AC-5-ACAGTGAT-TCTTTCCC_S124_R2_001.fastq.gz
  diff /data/h_vivian_weiss/2021/2021-AC-5-ACAGTGAT-TCTTTCCC_S124_R2_001.fastq.gz 2021-AC-5-ACAGTGAT-TCTTTCCC_S124_R2_001.fastq.gz
  status=$?
  if [[ $status -eq 0 ]]; then
    break
  fi
done

if [[ $status -ne 0 ]]; then
  echo file copied is not idential to original file, quit.
  touch $res_dir/2021-AC-5-ACAGTGAT-TCTTTCCC_S124_R2_001.fastq.gz.copy.failed
  rm 2021-AC-5-ACAGTGAT-TCTTTCCC_S124_R1_001.fastq.gz 2021-AC-5-ACAGTGAT-TCTTTCCC_S124_R2_001.fastq.gz
  exit 1
fi

ls *
echo localize end at `date`


rm -f X2021AC5.fastqc.failed

fastqc  --extract -t 2 -o `pwd` "2021-AC-5-ACAGTGAT-TCTTTCCC_S124_R1_001.fastq.gz" "2021-AC-5-ACAGTGAT-TCTTTCCC_S124_R2_001.fastq.gz" 2> >(tee X2021AC5.fastqc.stderr.log >&2)

status=$?
if [[ $status -ne 0 ]]; then
  touch X2021AC5.fastqc.failed
else
  touch X2021AC5.fastqc.succeed
fi

fastqc --version | cut -d ' ' -f2 | awk '{print "FastQC,"$1}' > `pwd`/fastqc.version

rm 2021-AC-5-ACAGTGAT-TCTTTCCC_S124_R1_001.fastq.gz 2021-AC-5-ACAGTGAT-TCTTTCCC_S124_R2_001.fastq.gz

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
