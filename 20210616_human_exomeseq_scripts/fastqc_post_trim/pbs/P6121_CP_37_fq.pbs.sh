
cd '/scratch/weissvl/shengq2/20210616_human_exomeseq/fastqc_post_trim/result/P6121_CP_37'

set -o pipefail



rm -f P6121_CP_37.fastqc.failed

fastqc  --extract -t 2 -o `pwd` "/scratch/weissvl/shengq2/20210616_human_exomeseq/cutadapt/result/P6121_CP_37_clipped.1.fastq.gz" "/scratch/weissvl/shengq2/20210616_human_exomeseq/cutadapt/result/P6121_CP_37_clipped.2.fastq.gz" 2> >(tee P6121_CP_37.fastqc.stderr.log >&2)

status=$?
if [[ $status -ne 0 ]]; then
  touch P6121_CP_37.fastqc.failed
else
  touch P6121_CP_37.fastqc.succeed
fi

fastqc --version | cut -d ' ' -f2 | awk '{print "FastQC,"$1}' > `pwd`/fastqc.version
