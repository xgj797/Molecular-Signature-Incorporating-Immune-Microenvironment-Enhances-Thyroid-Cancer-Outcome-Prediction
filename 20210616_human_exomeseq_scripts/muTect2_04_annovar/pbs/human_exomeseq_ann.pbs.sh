
cd '/scratch/weissvl/shengq2/20210616_human_exomeseq/muTect2_04_annovar/result/human_exomeseq'

set -o pipefail



if [[ ! -s human_exomeseq.filtered.annovar.hg38_multianno.txt && ! -s human_exomeseq.filtered.annovar.final.tsv ]]; then 
  convert2annovar.pl -format vcf4old /scratch/weissvl/shengq2/20210616_human_exomeseq/muTect2_03_filterDepth/result/human_exomeseq.filtered.vcf | cut -f1-7 | awk '{gsub(",\\*", "", $0); print}'> human_exomeseq.filtered.avinput 
  if [ -s human_exomeseq.filtered.avinput ]; then
    table_annovar.pl human_exomeseq.filtered.avinput /data/cqs/references/annovar/humandb/ -buildver hg38 -protocol refGene,avsnp150,cosmic70,exac03,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_eas,1000g2015aug_eur,1000g2015aug_sas,gnomad211_genome,clinvar_20190305 -operation g,f,f,f,f,f,f,f,f,f,f,f --remove --outfile human_exomeseq.filtered.annovar --remove
  fi
fi

echo find_protein_position_for_splicing=`date`
if [[ -s human_exomeseq.filtered.annovar.hg38_multianno.txt && ! -s human_exomeseq.filtered.annovar.splicing.hg38_multianno.txt ]]; then 
  python3 /data/cqs/softwares/ngsperl/lib/Annotation/annovarSplicing.py -i human_exomeseq.filtered.annovar.hg38_multianno.txt -d /data/cqs/references/annovar/humandb/ -o human_exomeseq.filtered.annovar.splicing.hg38_multianno.txt -b hg38 
fi

if [[ -s human_exomeseq.filtered.annovar.splicing.hg38_multianno.txt && ! -s human_exomeseq.filtered.annovar.final.tsv ]]; then
  sed -n '/^[^#]/q;p' /scratch/weissvl/shengq2/20210616_human_exomeseq/muTect2_03_filterDepth/result/human_exomeseq.filtered.vcf|sed '$ d' > human_exomeseq.filtered.annovar.final.tsv.header
  cat /scratch/weissvl/shengq2/20210616_human_exomeseq/muTect2_03_filterDepth/result/human_exomeseq.filtered.vcf | grep -v "^##" | cut -f7- > /scratch/weissvl/shengq2/20210616_human_exomeseq/muTect2_03_filterDepth/result/human_exomeseq.filtered.vcf.clean
  grep -v "^##" human_exomeseq.filtered.annovar.splicing.hg38_multianno.txt > human_exomeseq.filtered.annovar.splicing.hg38_multianno.txt.clean
  paste human_exomeseq.filtered.annovar.splicing.hg38_multianno.txt.clean /scratch/weissvl/shengq2/20210616_human_exomeseq/muTect2_03_filterDepth/result/human_exomeseq.filtered.vcf.clean > human_exomeseq.filtered.annovar.final.tsv.data
  cat human_exomeseq.filtered.annovar.final.tsv.header human_exomeseq.filtered.annovar.final.tsv.data > human_exomeseq.filtered.annovar.final.tsv
  rm /scratch/weissvl/shengq2/20210616_human_exomeseq/muTect2_03_filterDepth/result/human_exomeseq.filtered.vcf.clean human_exomeseq.filtered.annovar.splicing.hg38_multianno.txt.clean human_exomeseq.filtered.annovar.final.tsv.header human_exomeseq.filtered.annovar.final.tsv.data
fi
