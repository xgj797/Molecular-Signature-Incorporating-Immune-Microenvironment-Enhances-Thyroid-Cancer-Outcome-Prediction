
cd '/scratch/weissvl/shengq2/20210616_human_exomeseq/muTect2_05_filter/result'

set -o pipefail


 
R --vanilla -f /data/cqs/softwares/ngsperl/lib/Annotation/oncoPrint.r --args /scratch/weissvl/shengq2/20210616_human_exomeseq/muTect2_05_filter/result/human_exomeseq.freq0.001.snv.missense.tsv human_exomeseq.freq0.001.snv.missense.oncoprint.tsv onco_options.txt BRAF RAS NTRK2 CDKN2A CDKN2B NF1 KMT2D RB1 MMR ARID2 ATM
 
R --vanilla -f /data/cqs/softwares/ngsperl/lib/Annotation/CBioPortal.r --args /scratch/weissvl/shengq2/20210616_human_exomeseq/muTect2_05_filter/result/human_exomeseq.freq0.001.snv.missense.tsv human_exomeseq.freq0.001.snv.missense.cBioPortal . BRAF RAS NTRK2 CDKN2A CDKN2B NF1 KMT2D RB1 MMR ARID2 ATM
 
R --vanilla -f /data/cqs/softwares/ngsperl/lib/Annotation/oncoPrint.r --args /scratch/weissvl/shengq2/20210616_human_exomeseq/muTect2_05_filter/result/human_exomeseq.freq0.01.snv.missense.tsv human_exomeseq.freq0.01.snv.missense.oncoprint.tsv onco_options.txt BRAF RAS NTRK2 CDKN2A CDKN2B NF1 KMT2D RB1 MMR ARID2 ATM
 
R --vanilla -f /data/cqs/softwares/ngsperl/lib/Annotation/CBioPortal.r --args /scratch/weissvl/shengq2/20210616_human_exomeseq/muTect2_05_filter/result/human_exomeseq.freq0.01.snv.missense.tsv human_exomeseq.freq0.01.snv.missense.cBioPortal . BRAF RAS NTRK2 CDKN2A CDKN2B NF1 KMT2D RB1 MMR ARID2 ATM
 
R --vanilla -f /data/cqs/softwares/ngsperl/lib/Annotation/oncoPrint.r --args /scratch/weissvl/shengq2/20210616_human_exomeseq/muTect2_05_filter/result/human_exomeseq.freq0.1.snv.missense.tsv human_exomeseq.freq0.1.snv.missense.oncoprint.tsv onco_options.txt BRAF RAS NTRK2 CDKN2A CDKN2B NF1 KMT2D RB1 MMR ARID2 ATM
 
R --vanilla -f /data/cqs/softwares/ngsperl/lib/Annotation/CBioPortal.r --args /scratch/weissvl/shengq2/20210616_human_exomeseq/muTect2_05_filter/result/human_exomeseq.freq0.1.snv.missense.tsv human_exomeseq.freq0.1.snv.missense.cBioPortal . BRAF RAS NTRK2 CDKN2A CDKN2B NF1 KMT2D RB1 MMR ARID2 ATM
 
R --vanilla -f /data/cqs/softwares/ngsperl/lib/Annotation/oncoPrint.r --args /scratch/weissvl/shengq2/20210616_human_exomeseq/muTect2_05_filter/result/human_exomeseq.freq1.0.snv.missense.tsv human_exomeseq.freq1.0.snv.missense.oncoprint.tsv onco_options.txt BRAF RAS NTRK2 CDKN2A CDKN2B NF1 KMT2D RB1 MMR ARID2 ATM
 
R --vanilla -f /data/cqs/softwares/ngsperl/lib/Annotation/CBioPortal.r --args /scratch/weissvl/shengq2/20210616_human_exomeseq/muTect2_05_filter/result/human_exomeseq.freq1.0.snv.missense.tsv human_exomeseq.freq1.0.snv.missense.cBioPortal . BRAF RAS NTRK2 CDKN2A CDKN2B NF1 KMT2D RB1 MMR ARID2 ATM
 
python3 /data/cqs/softwares/ngsperl/lib/Annotation/geneDetails.py -i /scratch/weissvl/shengq2/20210616_human_exomeseq/muTect2_05_filter/result/human_exomeseq.freq0.001.filtered.missense.tsv -o human_exomeseq.freq0.001.filtered.missense.geneDetails.txt -g BRAF,RAS,NTRK2,CDKN2A,CDKN2B,NF1,KMT2D,RB1,MMR,ARID2,ATM
 
python3 /data/cqs/softwares/ngsperl/lib/Annotation/geneDetails.py -i /scratch/weissvl/shengq2/20210616_human_exomeseq/muTect2_05_filter/result/human_exomeseq.freq0.01.filtered.missense.tsv -o human_exomeseq.freq0.01.filtered.missense.geneDetails.txt -g BRAF,RAS,NTRK2,CDKN2A,CDKN2B,NF1,KMT2D,RB1,MMR,ARID2,ATM
 
python3 /data/cqs/softwares/ngsperl/lib/Annotation/geneDetails.py -i /scratch/weissvl/shengq2/20210616_human_exomeseq/muTect2_05_filter/result/human_exomeseq.freq0.1.filtered.missense.tsv -o human_exomeseq.freq0.1.filtered.missense.geneDetails.txt -g BRAF,RAS,NTRK2,CDKN2A,CDKN2B,NF1,KMT2D,RB1,MMR,ARID2,ATM
 
python3 /data/cqs/softwares/ngsperl/lib/Annotation/geneDetails.py -i /scratch/weissvl/shengq2/20210616_human_exomeseq/muTect2_05_filter/result/human_exomeseq.freq1.0.filtered.missense.tsv -o human_exomeseq.freq1.0.filtered.missense.geneDetails.txt -g BRAF,RAS,NTRK2,CDKN2A,CDKN2B,NF1,KMT2D,RB1,MMR,ARID2,ATM
 
python3 /data/cqs/softwares/ngsperl/lib/Annotation/geneFilter.py -i /scratch/weissvl/shengq2/20210616_human_exomeseq/muTect2_05_filter/result/human_exomeseq.freq0.001.filtered.missense.tsv -o human_exomeseq.freq0.001.filtered.missense.geneFilter.txt -g BRAF,RAS,NTRK2,CDKN2A,CDKN2B,NF1,KMT2D,RB1,MMR,ARID2,ATM
 
python3 /data/cqs/softwares/ngsperl/lib/Annotation/geneFilter.py -i /scratch/weissvl/shengq2/20210616_human_exomeseq/muTect2_05_filter/result/human_exomeseq.freq0.01.filtered.missense.tsv -o human_exomeseq.freq0.01.filtered.missense.geneFilter.txt -g BRAF,RAS,NTRK2,CDKN2A,CDKN2B,NF1,KMT2D,RB1,MMR,ARID2,ATM
 
python3 /data/cqs/softwares/ngsperl/lib/Annotation/geneFilter.py -i /scratch/weissvl/shengq2/20210616_human_exomeseq/muTect2_05_filter/result/human_exomeseq.freq0.1.filtered.missense.tsv -o human_exomeseq.freq0.1.filtered.missense.geneFilter.txt -g BRAF,RAS,NTRK2,CDKN2A,CDKN2B,NF1,KMT2D,RB1,MMR,ARID2,ATM
 
python3 /data/cqs/softwares/ngsperl/lib/Annotation/geneFilter.py -i /scratch/weissvl/shengq2/20210616_human_exomeseq/muTect2_05_filter/result/human_exomeseq.freq1.0.filtered.missense.tsv -o human_exomeseq.freq1.0.filtered.missense.geneFilter.txt -g BRAF,RAS,NTRK2,CDKN2A,CDKN2B,NF1,KMT2D,RB1,MMR,ARID2,ATM
