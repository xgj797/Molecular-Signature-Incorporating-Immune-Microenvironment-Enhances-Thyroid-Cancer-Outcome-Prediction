
cd '/scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_02_CollectReadCounts/result'

set -o pipefail


  

gatk --java-options "-Xmx40G" CollectReadCounts  \
  -L /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_01_PreprocessIntervals/result/human_exomeseq.preprocessed.interval_list \
  --input /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine/result/P1809_AC_083.rmdup.recal.bam \
  --reference /data/cqs/references/broad/hg38/v0/Homo_sapiens_assembly38.fasta \
  --format HDF5 \
  --interval-merging-rule OVERLAPPING_ONLY \
  --output P1809_AC_083.count.hdf5

status=$?
if [[ $status -ne 0 ]]; then
  touch P1809_AC_083.failed
  rm -f P1809_AC_083.count.hdf5
else
  touch P1809_AC_083.succeed
fi
            
rm -rf .cache .conda .config .theano
