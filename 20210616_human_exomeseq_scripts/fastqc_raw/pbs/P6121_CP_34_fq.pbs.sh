
cd '/scratch/weissvl/shengq2/20210616_human_exomeseq/fastqc_raw/result/P6121_CP_34'

set -o pipefail



rm -f P6121_CP_34.fastqc.failed

fastqc  --extract -t 2 -o `pwd` "/data/h_vivian_weiss/6121/6121-CP-34_S1_L005_R1_001.fastq.gz" "/data/h_vivian_weiss/6121/6121-CP-34_S1_L005_R2_001.fastq.gz" 2> >(tee P6121_CP_34.fastqc.stderr.log >&2)

status=$?
if [[ $status -ne 0 ]]; then
  touch P6121_CP_34.fastqc.failed
else
  touch P6121_CP_34.fastqc.succeed
fi

fastqc --version | cut -d ' ' -f2 | awk '{print "FastQC,"$1}' > `pwd`/fastqc.version
