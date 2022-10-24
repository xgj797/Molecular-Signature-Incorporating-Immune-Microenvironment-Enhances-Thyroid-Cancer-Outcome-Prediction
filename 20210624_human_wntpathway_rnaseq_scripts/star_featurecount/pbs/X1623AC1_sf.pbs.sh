
res_dir='/scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/star_featurecount/result'
tmp_dir=$(mktemp -d -t ci-$(date +%Y-%m-%d-%H-%M-%S)-XXXXXXXXXX)

tmp_cleaner()
{
rm -rf ${tmp_dir}
exit -1
}
trap 'tmp_cleaner' TERM

echo using tmp_dir=$tmp_dir
cd $tmp_dir


set -o pipefail




if [[ -e /scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/star_featurecount/result/X1623AC1.star.failed ]]; then
  rm /scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/star_featurecount/result/X1623AC1.star.failed
fi

if [[ -e /scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/star_featurecount/result/X1623AC1.featureCount.failed ]]; then
  rm /scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/star_featurecount/result/X1623AC1.featureCount.failed
fi

echo localize start at `date`

echo /scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/cutadapt/result/X1623AC1_clipped.1.fastq.gz      
if [[ ! -s /scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/cutadapt/result/X1623AC1_clipped.1.fastq.gz ]]; then
  echo file not exists: /scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/cutadapt/result/X1623AC1_clipped.1.fastq.gz
  touch $res_dir/X1623AC1_clipped.1.fastq.gz.not.exist
  rm X1623AC1_clipped.1.fastq.gz
  exit 1
fi

for i in {1..5}; do 
  echo iteration $i ...
  cp -f /scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/cutadapt/result/X1623AC1_clipped.1.fastq.gz X1623AC1_clipped.1.fastq.gz
  diff /scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/cutadapt/result/X1623AC1_clipped.1.fastq.gz X1623AC1_clipped.1.fastq.gz
  status=$?
  if [[ $status -eq 0 ]]; then
    break
  fi
done

if [[ $status -ne 0 ]]; then
  echo file copied is not idential to original file, quit.
  touch $res_dir/X1623AC1_clipped.1.fastq.gz.copy.failed
  rm X1623AC1_clipped.1.fastq.gz
  exit 1
fi

echo /scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/cutadapt/result/X1623AC1_clipped.2.fastq.gz      
if [[ ! -s /scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/cutadapt/result/X1623AC1_clipped.2.fastq.gz ]]; then
  echo file not exists: /scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/cutadapt/result/X1623AC1_clipped.2.fastq.gz
  touch $res_dir/X1623AC1_clipped.2.fastq.gz.not.exist
  rm X1623AC1_clipped.1.fastq.gz X1623AC1_clipped.2.fastq.gz
  exit 1
fi

for i in {1..5}; do 
  echo iteration $i ...
  cp -f /scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/cutadapt/result/X1623AC1_clipped.2.fastq.gz X1623AC1_clipped.2.fastq.gz
  diff /scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/cutadapt/result/X1623AC1_clipped.2.fastq.gz X1623AC1_clipped.2.fastq.gz
  status=$?
  if [[ $status -eq 0 ]]; then
    break
  fi
done

if [[ $status -ne 0 ]]; then
  echo file copied is not idential to original file, quit.
  touch $res_dir/X1623AC1_clipped.2.fastq.gz.copy.failed
  rm X1623AC1_clipped.1.fastq.gz X1623AC1_clipped.2.fastq.gz
  exit 1
fi

ls *
echo localize end at `date`


status=1

if [[ -s X1623AC1_clipped.1.fastq.gz ]]; then
  echo performing star ...
  STAR --twopassMode Basic --outSAMmapqUnique 60 --outSAMprimaryFlag AllBestScore --outSAMattrRGline ID:X1623AC1 SM:X1623AC1 LB:X1623AC1 PL:ILLUMINA PU:ILLUMINA --runThreadN 8 --genomeDir /data/cqs/references/gencode/GRCh38.p13/STAR_index_2.7.8a_v38_sjdb100 --readFilesIn X1623AC1_clipped.1.fastq.gz X1623AC1_clipped.2.fastq.gz  --readFilesCommand zcat --outFileNamePrefix X1623AC1_ --outSAMtype BAM Unsorted
  status=$?
  if [[ $status -eq 0 ]]; then
    touch X1623AC1.star.succeed
  else
    rm X1623AC1_Aligned.out.bam
    touch X1623AC1.star.failed
  fi

  STAR --version | awk '{print "STAR,v"$1}' > X1623AC1.count.star.version
  rm -rf X1623AC1__STARgenome X1623AC1__STARpass1 X1623AC1_Log.progress.out
fi

if [[ $status -eq 0 ]]; then
  echo bamStat=`date` 
  python3 /data/cqs/softwares/ngsperl/lib/Alignment/bamStat.py -i X1623AC1_Aligned.out.bam -o X1623AC1.bamstat
fi

if [[ $status -eq 0 ]]; then
  echo bamSort=`date` 
  samtools sort -m 5G -T X1623AC1 -t 8 -o X1623AC1_Aligned.sortedByCoord.out.bam X1623AC1_Aligned.out.bam && touch X1623AC1_Aligned.sortedByCoord.out.bam.succeed
  if [[ ! -e X1623AC1_Aligned.sortedByCoord.out.bam.succeed ]]; then
    rm X1623AC1_Aligned.sortedByCoord.out.bam
  else
    samtools index X1623AC1_Aligned.sortedByCoord.out.bam
    samtools idxstats X1623AC1_Aligned.sortedByCoord.out.bam > X1623AC1_Aligned.sortedByCoord.out.bam.chromosome.count
  fi
fi  


if [[ $status -eq 0 && -s X1623AC1_Aligned.sortedByCoord.out.bam ]]; then
  
  
  if [ ! -s X1623AC1.splicing.bed ]; then
    awk {'if($4=="2") print ""$1"\t"$2-$9-1"\t"$3+$9"\tJUNC000"NR"\t"$7+$8"\t-\t"$2-$9-1"\t"$3+$9"\t255,0,0\t2\t"$9","$9"\t","0,"$3-$2+$9+1; else if($4=="1") print ""$1"\t"$2-$9-1"\t"$3+$9"\tJUNC000"NR"\t"$7+$8"\t+\t"$2-$9-1"\t"$3+$9"\t0,0,255\t2\t"$9","$9"\t","0,"$3-$2+$9+1'} X1623AC1_SJ.out.tab > X1623AC1.splicing.bed
  fi
fi

if [[ $status -eq 0 ]]; then
  echo performing featureCounts ...
  featureCounts -g gene_id -t exon -p --countReadPairs -T 8 -a /data/cqs/references/gencode/GRCh38.p13/gencode.v38.annotation.gtf -o X1623AC1.count X1623AC1_Aligned.out.bam
  status=$?
  if [[ $status -eq 0 ]]; then
    touch X1623AC1.featureCount.succeed
  else
    touch X1623AC1.featureCount.failed
    rm X1623AC1.count
  fi

  featureCounts -v 2>&1 | grep featureCounts | cut -d ' ' -f2 | awk '{print "featureCounts,"$1}' > X1623AC1.count.featureCounts.version
fi 

if [[ -s X1623AC1.count && -s X1623AC1.bamstat ]]; then
  rm X1623AC1_Aligned.out.bam 
fi

rm X1623AC1_clipped.1.fastq.gz X1623AC1_clipped.2.fastq.gz

if [[ -d $tmp_dir && $tmp_dir != '/' ]]; then
  echo copy result from $tmp_dir to $res_dir
  #if the pbs was generated again during task is running, copy may be unpredictable. 
  #make sure to change to tmp_dir before copy result

  cd $tmp_dir
  cd $tmp_dir
  cd $tmp_dir
  cd $tmp_dir
  cd $tmp_dir
  cd $tmp_dir
  cd $tmp_dir
  cd $tmp_dir
  cd $tmp_dir
  cd $tmp_dir

  cp -p -r * $res_dir
  cd $res_dir
  echo delete tmp folder $tmp_dir
  rm -rf $tmp_dir
  echo move file and clean tmp folder done.
fi