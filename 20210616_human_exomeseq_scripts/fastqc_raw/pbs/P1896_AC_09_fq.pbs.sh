
cd '/scratch/weissvl/shengq2/20210616_human_exomeseq/fastqc_raw/result/P1896_AC_09'

set -o pipefail



rm -f P1896_AC_09.fastqc.failed

fastqc  --extract -t 2 -o `pwd` "/data/h_vivian_weiss/DNAseq/1896/1896-AC-9-AAGCGCGCTT-CACCTGTTGC_S52_R1_001.fastq.gz" "/data/h_vivian_weiss/DNAseq/1896/1896-AC-9-AAGCGCGCTT-CACCTGTTGC_S52_R2_001.fastq.gz" 2> >(tee P1896_AC_09.fastqc.stderr.log >&2)

status=$?
if [[ $status -ne 0 ]]; then
  touch P1896_AC_09.fastqc.failed
else
  touch P1896_AC_09.fastqc.succeed
fi

fastqc --version | cut -d ' ' -f2 | awk '{print "FastQC,"$1}' > `pwd`/fastqc.version
