
cd '/scratch/weissvl/shengq2/20210616_human_exomeseq/fastqc_raw/result'



fastqc  --extract -t 2 -o /scratch/weissvl/shengq2/20210616_human_exomeseq/fastqc_raw/result/P4601_CP_69 "/data/h_vivian_weiss/4601/4601-CP-69_S5_L005_R1_001.fastq.gz" "/data/h_vivian_weiss/4601/4601-CP-69_S5_L005_R2_001.fastq.gz" 
fastqc --version | cut -d ' ' -f2 | awk '{print "FastQC,"$1}' > /scratch/weissvl/shengq2/20210616_human_exomeseq/fastqc_raw/result/P4601_CP_69/fastqc.version
