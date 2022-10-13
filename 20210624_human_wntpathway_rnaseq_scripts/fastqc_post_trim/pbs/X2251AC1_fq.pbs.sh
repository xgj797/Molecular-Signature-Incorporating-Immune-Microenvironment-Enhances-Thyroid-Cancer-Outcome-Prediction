
res_dir='/scratch/weissvl/shengq2/20200831_human_wntpathway_rnaseq.v2/fastqc_post_trim/result/X2251AC1'
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

echo /scratch/weissvl/shengq2/20200831_human_wntpathway_rnaseq.v2/cutadapt/result/X2251AC1_clipped.1.fastq.gz      
if [[ ! -s /scratch/weissvl/shengq2/20200831_human_wntpathway_rnaseq.v2/cutadapt/result/X2251AC1_clipped.1.fastq.gz ]]; then
  echo file not exists: /scratch/weissvl/shengq2/20200831_human_wntpathway_rnaseq.v2/cutadapt/result/X2251AC1_clipped.1.fastq.gz
  touch $res_dir/X2251AC1_clipped.1.fastq.gz.not.exist
  rm X2251AC1_clipped.1.fastq.gz
  exit 1
fi

for i in {1..5}; do 
  echo iteration $i ...
  cp -f /scratch/weissvl/shengq2/20200831_human_wntpathway_rnaseq.v2/cutadapt/result/X2251AC1_clipped.1.fastq.gz X2251AC1_clipped.1.fastq.gz
  diff /scratch/weissvl/shengq2/20200831_human_wntpathway_rnaseq.v2/cutadapt/result/X2251AC1_clipped.1.fastq.gz X2251AC1_clipped.1.fastq.gz
  status=$?
  if [[ $status -eq 0 ]]; then
    break
  fi
done

if [[ $status -ne 0 ]]; then
  echo file copied is not idential to original file, quit.
  touch $res_dir/X2251AC1_clipped.1.fastq.gz.copy.failed
  rm X2251AC1_clipped.1.fastq.gz
  exit 1
fi

echo /scratch/weissvl/shengq2/20200831_human_wntpathway_rnaseq.v2/cutadapt/result/X2251AC1_clipped.2.fastq.gz      
if [[ ! -s /scratch/weissvl/shengq2/20200831_human_wntpathway_rnaseq.v2/cutadapt/result/X2251AC1_clipped.2.fastq.gz ]]; then
  echo file not exists: /scratch/weissvl/shengq2/20200831_human_wntpathway_rnaseq.v2/cutadapt/result/X2251AC1_clipped.2.fastq.gz
  touch $res_dir/X2251AC1_clipped.2.fastq.gz.not.exist
  rm X2251AC1_clipped.1.fastq.gz X2251AC1_clipped.2.fastq.gz
  exit 1
fi

for i in {1..5}; do 
  echo iteration $i ...
  cp -f /scratch/weissvl/shengq2/20200831_human_wntpathway_rnaseq.v2/cutadapt/result/X2251AC1_clipped.2.fastq.gz X2251AC1_clipped.2.fastq.gz
  diff /scratch/weissvl/shengq2/20200831_human_wntpathway_rnaseq.v2/cutadapt/result/X2251AC1_clipped.2.fastq.gz X2251AC1_clipped.2.fastq.gz
  status=$?
  if [[ $status -eq 0 ]]; then
    break
  fi
done

if [[ $status -ne 0 ]]; then
  echo file copied is not idential to original file, quit.
  touch $res_dir/X2251AC1_clipped.2.fastq.gz.copy.failed
  rm X2251AC1_clipped.1.fastq.gz X2251AC1_clipped.2.fastq.gz
  exit 1
fi

ls *
echo localize end at `date`


fastqc  --extract -t 2 -o `pwd` "X2251AC1_clipped.1.fastq.gz" "X2251AC1_clipped.2.fastq.gz" 2> >(tee X2251AC1.fastqc.stderr.log >&2)

fastqc --version | cut -d ' ' -f2 | awk '{print "FastQC,"$1}' > `pwd`/fastqc.version

rm X2251AC1_clipped.1.fastq.gz X2251AC1_clipped.2.fastq.gz

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
