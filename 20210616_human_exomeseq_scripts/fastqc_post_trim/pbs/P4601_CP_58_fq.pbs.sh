
cd '/scratch/weissvl/shengq2/20210616_human_exomeseq/fastqc_post_trim/result/P4601_CP_58'

set -o pipefail



rm -f P4601_CP_58.fastqc.failed

fastqc  --extract -t 2 -o `pwd` "/scratch/weissvl/shengq2/20210616_human_exomeseq/cutadapt/result/P4601_CP_58_clipped.1.fastq.gz" "/scratch/weissvl/shengq2/20210616_human_exomeseq/cutadapt/result/P4601_CP_58_clipped.2.fastq.gz" 2> >(tee P4601_CP_58.fastqc.stderr.log >&2)

status=$?
if [[ $status -ne 0 ]]; then
  touch P4601_CP_58.fastqc.failed
else
  touch P4601_CP_58.fastqc.succeed
fi

fastqc --version | cut -d ' ' -f2 | awk '{print "FastQC,"$1}' > `pwd`/fastqc.version
