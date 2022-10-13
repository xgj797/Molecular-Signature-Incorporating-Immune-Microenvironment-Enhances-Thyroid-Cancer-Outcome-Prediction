
cd '/scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_06_toMAF/result'

set -o pipefail



Rscript --vanilla /data/cqs/softwares/ngsperl/lib/Annotation/annovar2maf.r /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_05_filter/result/human_exomeseq.freq0.001.filtered.tsv human_exomeseq.freq0.001.filtered.tsv.tmp hg38
  
if [[ -s human_exomeseq.freq0.001.filtered.tsv.tmp ]]; then
  python3 /data/cqs/softwares/ngsperl/lib/Annotation/annovar2maf.py -i human_exomeseq.freq0.001.filtered.tsv.tmp -o human_exomeseq.freq0.001.filtered.tsv.maf
fi

if [[ -s human_exomeseq.freq0.001.filtered.tsv.maf ]]; then
  rm human_exomeseq.freq0.001.filtered.tsv.tmp
fi


Rscript --vanilla /data/cqs/softwares/ngsperl/lib/Annotation/annovar2maf.r /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_05_filter/result/human_exomeseq.freq0.01.filtered.tsv human_exomeseq.freq0.01.filtered.tsv.tmp hg38
  
if [[ -s human_exomeseq.freq0.01.filtered.tsv.tmp ]]; then
  python3 /data/cqs/softwares/ngsperl/lib/Annotation/annovar2maf.py -i human_exomeseq.freq0.01.filtered.tsv.tmp -o human_exomeseq.freq0.01.filtered.tsv.maf
fi

if [[ -s human_exomeseq.freq0.01.filtered.tsv.maf ]]; then
  rm human_exomeseq.freq0.01.filtered.tsv.tmp
fi


Rscript --vanilla /data/cqs/softwares/ngsperl/lib/Annotation/annovar2maf.r /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_05_filter/result/human_exomeseq.freq0.1.filtered.tsv human_exomeseq.freq0.1.filtered.tsv.tmp hg38
  
if [[ -s human_exomeseq.freq0.1.filtered.tsv.tmp ]]; then
  python3 /data/cqs/softwares/ngsperl/lib/Annotation/annovar2maf.py -i human_exomeseq.freq0.1.filtered.tsv.tmp -o human_exomeseq.freq0.1.filtered.tsv.maf
fi

if [[ -s human_exomeseq.freq0.1.filtered.tsv.maf ]]; then
  rm human_exomeseq.freq0.1.filtered.tsv.tmp
fi

