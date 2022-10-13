
res_dir='/scratch/weissvl/shengq2/20210616_human_exomeseq/fastqc_post_trim/result/P1809_AC_062'
tmp_dir=$(mktemp -d -t ci-$(date +%Y-%m-%d-%H-%M-%S)-XXXXXXXXXX)

tmp_cleaner()
{
rm -rf ${tmp_dir}
exit -1
}
trap 'tmp_cleaner' TERM

echo using tmp_dir=$tmp_dir
cd $tmp_dir





echo localize start at `date`

echo /scratch/weissvl/shengq2/20210616_human_exomeseq/cutadapt/result/P1809_AC_062_clipped.1.fastq.gz      
if [[ ! -s /scratch/weissvl/shengq2/20210616_human_exomeseq/cutadapt/result/P1809_AC_062_clipped.1.fastq.gz ]]; then
  echo file not exists: /scratch/weissvl/shengq2/20210616_human_exomeseq/cutadapt/result/P1809_AC_062_clipped.1.fastq.gz
  exit 1
fi
cp /scratch/weissvl/shengq2/20210616_human_exomeseq/cutadapt/result/P1809_AC_062_clipped.1.fastq.gz P1809_AC_062_clipped.1.fastq.gz

echo /scratch/weissvl/shengq2/20210616_human_exomeseq/cutadapt/result/P1809_AC_062_clipped.2.fastq.gz      
if [[ ! -s /scratch/weissvl/shengq2/20210616_human_exomeseq/cutadapt/result/P1809_AC_062_clipped.2.fastq.gz ]]; then
  echo file not exists: /scratch/weissvl/shengq2/20210616_human_exomeseq/cutadapt/result/P1809_AC_062_clipped.2.fastq.gz
  exit 1
fi
cp /scratch/weissvl/shengq2/20210616_human_exomeseq/cutadapt/result/P1809_AC_062_clipped.2.fastq.gz P1809_AC_062_clipped.2.fastq.gz

ls *
echo localize end at `date`

fastqc  --extract -t 2 -o `pwd` "P1809_AC_062_clipped.1.fastq.gz" "P1809_AC_062_clipped.2.fastq.gz" 1> >(tee P1809_AC_062.fastqc.stdout.log) 2> >(tee P1809_AC_062.fastqc.stderr.log >&2) 
fastqc --version | cut -d ' ' -f2 | awk '{print "FastQC,"$1}' > `pwd`/fastqc.version

rm P1809_AC_062_clipped.1.fastq.gz P1809_AC_062_clipped.2.fastq.gz

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
