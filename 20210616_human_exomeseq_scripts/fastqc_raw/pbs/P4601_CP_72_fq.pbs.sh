
cd '/scratch/weissvl/shengq2/20210616_human_exomeseq/fastqc_raw/result/P4601_CP_72'

set -o pipefail



rm -f P4601_CP_72.fastqc.failed

fastqc  --extract -t 2 -o `pwd` "/data/h_vivian_weiss/4601/4601-CP-72_S73_L005_R1_001.fastq.gz" "/data/h_vivian_weiss/4601/4601-CP-72_S73_L005_R2_001.fastq.gz" 2> >(tee P4601_CP_72.fastqc.stderr.log >&2)

status=$?
if [[ $status -ne 0 ]]; then
  touch P4601_CP_72.fastqc.failed
else
  touch P4601_CP_72.fastqc.succeed
fi

fastqc --version | cut -d ' ' -f2 | awk '{print "FastQC,"$1}' > `pwd`/fastqc.version
