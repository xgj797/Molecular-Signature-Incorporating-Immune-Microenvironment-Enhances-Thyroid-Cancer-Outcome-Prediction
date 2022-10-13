
cd '/scratch/weissvl/shengq2/20210616_human_exomeseq/fastqc_raw/result/P1988_AC_16'

set -o pipefail



rm -f P1988_AC_16.fastqc.failed

fastqc  --extract -t 2 -o `pwd` "/data/h_vivian_weiss/DNAseq/1988/1988-AC-16-GTTGGATGAA-TCATAGATTG_S98_R1_001.fastq.gz" "/data/h_vivian_weiss/DNAseq/1988/1988-AC-16-GTTGGATGAA-TCATAGATTG_S98_R2_001.fastq.gz" 2> >(tee P1988_AC_16.fastqc.stderr.log >&2)

status=$?
if [[ $status -ne 0 ]]; then
  touch P1988_AC_16.fastqc.failed
else
  touch P1988_AC_16.fastqc.succeed
fi

fastqc --version | cut -d ' ' -f2 | awk '{print "FastQC,"$1}' > `pwd`/fastqc.version
