
cd '/scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_11_AnnotationGenesPlot/result'

set -o pipefail





python3 /data/cqs/softwares/ngsperl/lib/CQS/../Visualization/plotCNV.py   -i /scratch/weissvl/shengq2/20210616_human_exomeseq/annotation_genes_locus/result/human_exomeseq.bed -c /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_07_CombineGCNV/result/human_exomeseq.txt -s /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_08_SizeFactor/result/human_exomeseq.txt.sizefactor -b human_exomeseq__fileList1.list -o human_exomeseq.position.txt
