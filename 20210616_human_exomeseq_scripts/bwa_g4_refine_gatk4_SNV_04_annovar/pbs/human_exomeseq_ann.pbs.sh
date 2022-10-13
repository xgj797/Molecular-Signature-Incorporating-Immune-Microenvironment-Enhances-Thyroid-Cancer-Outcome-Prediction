
cd '/scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_04_annovar/result/human_exomeseq'

set -o pipefail



if [[ ! -s human_exomeseq.maf_filtered.annovar.hg38_multianno.txt && ! -s human_exomeseq.maf_filtered.annovar.final.tsv ]]; then 
  convert2annovar.pl -format vcf4old /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_03_filterMAF/result/human_exomeseq.maf_filtered.vcf.gz | cut -f1-7 | awk '{gsub(",\\*", "", $0); print}'> human_exomeseq.maf_filtered.avinput 
  if [ -s human_exomeseq.maf_filtered.avinput ]; then
    table_annovar.pl human_exomeseq.maf_filtered.avinput /data/cqs/references/annovar/humandb/ -buildver hg38 -protocol refGene,avsnp150,cosmic70,exac03,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_eas,1000g2015aug_eur,1000g2015aug_sas,gnomad211_genome,clinvar_20190305 -operation g,f,f,f,f,f,f,f,f,f,f,f --remove --outfile human_exomeseq.maf_filtered.annovar --remove
  fi
fi

echo find_protein_position_for_splicing=`date`
if [[ -s human_exomeseq.maf_filtered.annovar.hg38_multianno.txt && ! -s human_exomeseq.maf_filtered.annovar.splicing.hg38_multianno.txt ]]; then 
  python3 /data/cqs/softwares/ngsperl/lib/Annotation/annovarSplicing.py -i human_exomeseq.maf_filtered.annovar.hg38_multianno.txt -d /data/cqs/references/annovar/humandb/ -o human_exomeseq.maf_filtered.annovar.splicing.hg38_multianno.txt -b hg38 
fi

if [[ -s human_exomeseq.maf_filtered.annovar.splicing.hg38_multianno.txt && ! -s human_exomeseq.maf_filtered.annovar.final.tsv ]]; then
  sed -n '/^[^#]/q;p' /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_03_filterMAF/result/human_exomeseq.maf_filtered.vcf.gz|sed '$ d' > human_exomeseq.maf_filtered.annovar.final.tsv.header
  zcat /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_03_filterMAF/result/human_exomeseq.maf_filtered.vcf.gz | grep -v "^##" | cut -f7- > /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_03_filterMAF/result/human_exomeseq.maf_filtered.vcf.gz.clean
  grep -v "^##" human_exomeseq.maf_filtered.annovar.splicing.hg38_multianno.txt > human_exomeseq.maf_filtered.annovar.splicing.hg38_multianno.txt.clean
  paste human_exomeseq.maf_filtered.annovar.splicing.hg38_multianno.txt.clean /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_03_filterMAF/result/human_exomeseq.maf_filtered.vcf.gz.clean > human_exomeseq.maf_filtered.annovar.final.tsv.data
  cat human_exomeseq.maf_filtered.annovar.final.tsv.header human_exomeseq.maf_filtered.annovar.final.tsv.data > human_exomeseq.maf_filtered.annovar.final.tsv
  rm /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_03_filterMAF/result/human_exomeseq.maf_filtered.vcf.gz.clean human_exomeseq.maf_filtered.annovar.splicing.hg38_multianno.txt.clean human_exomeseq.maf_filtered.annovar.final.tsv.header human_exomeseq.maf_filtered.annovar.final.tsv.data
fi
