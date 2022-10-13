
cd '/scratch/weissvl/shengq2/20210616_human_exomeseq/fastqc_post_trim/result/P3544_SB_17'

set -o pipefail



rm -f P3544_SB_17.fastqc.failed

fastqc  --extract -t 2 -o `pwd` "/scratch/weissvl/shengq2/20210616_human_exomeseq/cutadapt/result/P3544_SB_17_clipped.1.fastq.gz" "/scratch/weissvl/shengq2/20210616_human_exomeseq/cutadapt/result/P3544_SB_17_clipped.2.fastq.gz" 2> >(tee P3544_SB_17.fastqc.stderr.log >&2)

status=$?
if [[ $status -ne 0 ]]; then
  touch P3544_SB_17.fastqc.failed
else
  touch P3544_SB_17.fastqc.succeed
fi

fastqc --version | cut -d ' ' -f2 | awk '{print "FastQC,"$1}' > `pwd`/fastqc.version
