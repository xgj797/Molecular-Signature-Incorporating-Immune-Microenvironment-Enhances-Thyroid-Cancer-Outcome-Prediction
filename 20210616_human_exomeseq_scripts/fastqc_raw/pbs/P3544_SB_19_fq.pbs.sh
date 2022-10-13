
cd '/scratch/weissvl/shengq2/20210616_human_exomeseq/fastqc_raw/result/P3544_SB_19'

set -o pipefail



rm -f P3544_SB_19.fastqc.failed

fastqc  --extract -t 2 -o `pwd` "/data/h_vivian_weiss/3544/3544-SB-19-GTTGTAGT-TATCACTC_S41_R1_001.fastq.gz" "/data/h_vivian_weiss/3544/3544-SB-19-GTTGTAGT-TATCACTC_S41_R2_001.fastq.gz" 2> >(tee P3544_SB_19.fastqc.stderr.log >&2)

status=$?
if [[ $status -ne 0 ]]; then
  touch P3544_SB_19.fastqc.failed
else
  touch P3544_SB_19.fastqc.succeed
fi

fastqc --version | cut -d ' ' -f2 | awk '{print "FastQC,"$1}' > `pwd`/fastqc.version
