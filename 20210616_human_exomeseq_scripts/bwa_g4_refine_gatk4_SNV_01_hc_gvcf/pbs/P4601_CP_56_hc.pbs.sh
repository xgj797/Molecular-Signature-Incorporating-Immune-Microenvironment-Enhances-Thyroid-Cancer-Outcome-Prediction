
cd '/scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_01_hc_gvcf/result'

set -o pipefail



gatk --java-options "-Xmx40G" \
  HaplotypeCaller  -ERC GVCF -L /data/cqs/references/exomeseq/IDT/Exome-IDT_V1.hg38.slop50.bed -XL /data/cqs/references/blacklist_files/hg38-blacklist.v2.bed \
  --native-pair-hmm-threads 8 \
  -R /data/cqs/references/broad/hg38/v0/bwa_index_0.7.17/Homo_sapiens_assembly38.fasta \
  -I /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine/result/P4601_CP_56.rmdup.recal.bam \
  -O P4601_CP_56.g.vcf.gz

status=$?
if [[ $status -ne 0 ]]; then
  touch P4601_CP_56.hc.failed
  rm -f P4601_CP_56.g.vcf.gz P4601_CP_56.g.vcf.gz.tbi
else
  touch P4601_CP_56.hc.succeed
fi

