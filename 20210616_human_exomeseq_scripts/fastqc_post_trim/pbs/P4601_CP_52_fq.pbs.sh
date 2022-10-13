
cd '/scratch/weissvl/shengq2/20210616_human_exomeseq/fastqc_post_trim/result'



fastqc  --extract -t 2 -o /scratch/weissvl/shengq2/20210616_human_exomeseq/fastqc_post_trim/result/P4601_CP_52 "/scratch/weissvl/shengq2/20210616_human_exomeseq/cutadapt/result/P4601_CP_52_clipped.1.fastq.gz" "/scratch/weissvl/shengq2/20210616_human_exomeseq/cutadapt/result/P4601_CP_52_clipped.2.fastq.gz" 
fastqc --version | cut -d ' ' -f2 | awk '{print "FastQC,"$1}' > /scratch/weissvl/shengq2/20210616_human_exomeseq/fastqc_post_trim/result/P4601_CP_52/fastqc.version
