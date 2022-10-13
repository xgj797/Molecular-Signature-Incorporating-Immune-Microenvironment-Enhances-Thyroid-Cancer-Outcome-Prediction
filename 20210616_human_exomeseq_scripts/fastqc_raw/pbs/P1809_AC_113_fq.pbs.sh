
cd '/scratch/weissvl/shengq2/20210616_human_exomeseq/fastqc_raw/result/P1809_AC_113'

set -o pipefail



rm -f P1809_AC_113.fastqc.failed

fastqc  --extract -t 2 -o `pwd` "/data/h_vivian_weiss/DNAseq/1809/1809-AC-113-TTACCTGGAA-GATTCTGAAT_S111_R1_001.fastq.gz" "/data/h_vivian_weiss/DNAseq/1809/1809-AC-113-TTACCTGGAA-GATTCTGAAT_S111_R2_001.fastq.gz" 2> >(tee P1809_AC_113.fastqc.stderr.log >&2)

status=$?
if [[ $status -ne 0 ]]; then
  touch P1809_AC_113.fastqc.failed
else
  touch P1809_AC_113.fastqc.succeed
fi

fastqc --version | cut -d ' ' -f2 | awk '{print "FastQC,"$1}' > `pwd`/fastqc.version
