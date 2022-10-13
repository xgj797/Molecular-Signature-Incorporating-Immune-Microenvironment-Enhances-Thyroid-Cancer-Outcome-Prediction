
cd '/scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_05_filter/result'

set -o pipefail


 
python3 /data/cqs/softwares/ngsperl/lib/Annotation/filterAnnovar.py --exac_key ExAC_ALL --g1000_key 1000g2015aug_all --gnomad_key AF -i /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_04_annovar/result/human_exomeseq/human_exomeseq.maf_filtered.annovar.final.tsv -t 0.001 -o human_exomeseq.freq0.001      
 
python3 /data/cqs/softwares/ngsperl/lib/Annotation/filterAnnovar.py --exac_key ExAC_ALL --g1000_key 1000g2015aug_all --gnomad_key AF -i /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_04_annovar/result/human_exomeseq/human_exomeseq.maf_filtered.annovar.final.tsv -t 0.01 -o human_exomeseq.freq0.01      
 
python3 /data/cqs/softwares/ngsperl/lib/Annotation/filterAnnovar.py --exac_key ExAC_ALL --g1000_key 1000g2015aug_all --gnomad_key AF -i /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_04_annovar/result/human_exomeseq/human_exomeseq.maf_filtered.annovar.final.tsv -t 0.1 -o human_exomeseq.freq0.1      
 
python3 /data/cqs/softwares/ngsperl/lib/Annotation/filterAnnovar.py --exac_key ExAC_ALL --g1000_key 1000g2015aug_all --gnomad_key AF -i /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_SNV_04_annovar/result/human_exomeseq/human_exomeseq.maf_filtered.annovar.final.tsv -t 1.0 -o human_exomeseq.freq1.0      
