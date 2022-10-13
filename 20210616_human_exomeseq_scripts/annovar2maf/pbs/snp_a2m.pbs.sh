
cd '/workspace/shengq2/vivian_weiss/20210616_human_exomeseq/annovar2maf/result'

set -o pipefail



Rscript --vanilla /data/cqs/softwares/ngsperl/lib/Annotation/annovar2maf.r /workspace/shengq2/vivian_weiss/20210616_human_exomeseq/to_maf/WES_7633_JN.freq0.001.filtered.missense.tsv WES_7633_JN.freq0.001.filtered.missense.tsv.tmp hg38
  
if [[ -s WES_7633_JN.freq0.001.filtered.missense.tsv.tmp ]]; then
  python3 /data/cqs/softwares/ngsperl/lib/Annotation/annovar2maf.py -i WES_7633_JN.freq0.001.filtered.missense.tsv.tmp -o WES_7633_JN.freq0.001.filtered.missense.tsv.maf
fi

if [[ -s WES_7633_JN.freq0.001.filtered.missense.tsv.maf ]]; then
  rm WES_7633_JN.freq0.001.filtered.missense.tsv.tmp
fi

