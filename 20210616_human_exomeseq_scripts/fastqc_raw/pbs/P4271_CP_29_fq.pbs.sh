
cd '/scratch/weissvl/shengq2/20210616_human_exomeseq/fastqc_raw/result/P4271_CP_29'

set -o pipefail



rm -f P4271_CP_29.fastqc.failed

fastqc  --extract -t 2 -o `pwd` "/data/h_vivian_weiss/4271/4271-CP-29-CCGGAGTG-TCCGCCAA_S118_R1_001.fastq.gz" "/data/h_vivian_weiss/4271/4271-CP-29-CCGGAGTG-TCCGCCAA_S118_R2_001.fastq.gz" 2> >(tee P4271_CP_29.fastqc.stderr.log >&2)

status=$?
if [[ $status -ne 0 ]]; then
  touch P4271_CP_29.fastqc.failed
else
  touch P4271_CP_29.fastqc.succeed
fi

fastqc --version | cut -d ' ' -f2 | awk '{print "FastQC,"$1}' > `pwd`/fastqc.version
