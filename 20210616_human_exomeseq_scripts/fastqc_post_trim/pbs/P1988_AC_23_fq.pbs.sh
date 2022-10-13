
cd '/scratch/weissvl/shengq2/20210616_human_exomeseq/fastqc_post_trim/result/P1988_AC_23'

set -o pipefail



rm -f P1988_AC_23.fastqc.failed

fastqc  --extract -t 2 -o `pwd` "/scratch/weissvl/shengq2/20210616_human_exomeseq/cutadapt/result/P1988_AC_23_clipped.1.fastq.gz" "/scratch/weissvl/shengq2/20210616_human_exomeseq/cutadapt/result/P1988_AC_23_clipped.2.fastq.gz" 2> >(tee P1988_AC_23.fastqc.stderr.log >&2)

status=$?
if [[ $status -ne 0 ]]; then
  touch P1988_AC_23.fastqc.failed
else
  touch P1988_AC_23.fastqc.succeed
fi

fastqc --version | cut -d ' ' -f2 | awk '{print "FastQC,"$1}' > `pwd`/fastqc.version
