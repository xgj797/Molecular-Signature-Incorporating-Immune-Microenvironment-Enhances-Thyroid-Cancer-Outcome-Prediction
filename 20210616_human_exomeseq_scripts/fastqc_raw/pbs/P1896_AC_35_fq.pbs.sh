
cd '/scratch/weissvl/shengq2/20210616_human_exomeseq/fastqc_raw/result/P1896_AC_35'

set -o pipefail



rm -f P1896_AC_35.fastqc.failed

fastqc  --extract -t 2 -o `pwd` "/data/h_vivian_weiss/DNAseq/1896/1896-AC-35-CCGGTTCCTA-TCCATTGCCG_S78_R1_001.fastq.gz" "/data/h_vivian_weiss/DNAseq/1896/1896-AC-35-CCGGTTCCTA-TCCATTGCCG_S78_R2_001.fastq.gz" 2> >(tee P1896_AC_35.fastqc.stderr.log >&2)

status=$?
if [[ $status -ne 0 ]]; then
  touch P1896_AC_35.fastqc.failed
else
  touch P1896_AC_35.fastqc.succeed
fi

fastqc --version | cut -d ' ' -f2 | awk '{print "FastQC,"$1}' > `pwd`/fastqc.version
