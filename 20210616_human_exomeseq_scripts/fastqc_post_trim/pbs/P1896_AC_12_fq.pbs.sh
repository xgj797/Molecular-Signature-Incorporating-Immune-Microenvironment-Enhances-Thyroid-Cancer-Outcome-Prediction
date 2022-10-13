
cd '/scratch/weissvl/shengq2/20210616_human_exomeseq/fastqc_post_trim/result/P1896_AC_12'

set -o pipefail



rm -f P1896_AC_12.fastqc.failed

fastqc  --extract -t 2 -o `pwd` "/scratch/weissvl/shengq2/20210616_human_exomeseq/cutadapt/result/P1896_AC_12_clipped.1.fastq.gz" "/scratch/weissvl/shengq2/20210616_human_exomeseq/cutadapt/result/P1896_AC_12_clipped.2.fastq.gz" 2> >(tee P1896_AC_12.fastqc.stderr.log >&2)

status=$?
if [[ $status -ne 0 ]]; then
  touch P1896_AC_12.fastqc.failed
else
  touch P1896_AC_12.fastqc.succeed
fi

fastqc --version | cut -d ' ' -f2 | awk '{print "FastQC,"$1}' > `pwd`/fastqc.version
