
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




if [[ -e /scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/star_featurecount/result/X3607SB8.star.failed ]]; then
  rm /scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/star_featurecount/result/X3607SB8.star.failed
fi

if [[ -e /scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/star_featurecount/result/X3607SB8.featureCount.failed ]]; then
  rm /scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/star_featurecount/result/X3607SB8.featureCount.failed
fi

echo localize start at `date`

echo /scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/cutadapt/result/X3607SB8_clipped.1.fastq.gz      
if [[ ! -s /scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/cutadapt/result/X3607SB8_clipped.1.fastq.gz ]]; then
  echo file not exists: /scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/cutadapt/result/X3607SB8_clipped.1.fastq.gz
  touch $res_dir/X3607SB8_clipped.1.fastq.gz.not.exist
  rm X3607SB8_clipped.1.fastq.gz
  exit 1
fi

for i in {1..5}; do 
  echo iteration $i ...
  cp -f /scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/cutadapt/result/X3607SB8_clipped.1.fastq.gz X3607SB8_clipped.1.fastq.gz
  diff /scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/cutadapt/result/X3607SB8_clipped.1.fastq.gz X3607SB8_clipped.1.fastq.gz
  status=$?
  if [[ $status -eq 0 ]]; then
    break
  fi
done

if [[ $status -ne 0 ]]; then
  echo file copied is not idential to original file, quit.
  touch $res_dir/X3607SB8_clipped.1.fastq.gz.copy.failed
  rm X3607SB8_clipped.1.fastq.gz
  exit 1
fi

echo /scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/cutadapt/result/X3607SB8_clipped.2.fastq.gz      
if [[ ! -s /scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/cutadapt/result/X3607SB8_clipped.2.fastq.gz ]]; then
  echo file not exists: /scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/cutadapt/result/X3607SB8_clipped.2.fastq.gz
  touch $res_dir/X3607SB8_clipped.2.fastq.gz.not.exist
  rm X3607SB8_clipped.1.fastq.gz X3607SB8_clipped.2.fastq.gz
  exit 1
fi

for i in {1..5}; do 
  echo iteration $i ...
  cp -f /scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/cutadapt/result/X3607SB8_clipped.2.fastq.gz X3607SB8_clipped.2.fastq.gz
  diff /scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/cutadapt/result/X3607SB8_clipped.2.fastq.gz X3607SB8_clipped.2.fastq.gz
  status=$?
  if [[ $status -eq 0 ]]; then
    break
  fi
done

if [[ $status -ne 0 ]]; then
  echo file copied is not idential to original file, quit.
  touch $res_dir/X3607SB8_clipped.2.fastq.gz.copy.failed
  rm X3607SB8_clipped.1.fastq.gz X3607SB8_clipped.2.fastq.gz
  exit 1
fi

ls *
echo localize end at `date`


status=1

if [[ -s X3607SB8_clipped.1.fastq.gz ]]; then
  echo performing star ...
  STAR --twopassMode Basic --outSAMmapqUnique 60 --outSAMprimaryFlag AllBestScore --outSAMattrRGline ID:X3607SB8 SM:X3607SB8 LB:X3607SB8 PL:ILLUMINA PU:ILLUMINA --runThreadN 8 --genomeDir /data/cqs/references/gencode/GRCh38.p13/STAR_index_2.7.8a_v38_sjdb100 --readFilesIn X3607SB8_clipped.1.fastq.gz X3607SB8_clipped.2.fastq.gz  --readFilesCommand zcat --outFileNamePrefix X3607SB8_ --outSAMtype BAM Unsorted
  status=$?
  if [[ $status -eq 0 ]]; then
    touch X3607SB8.star.succeed
  else
    rm X3607SB8_Aligned.out.bam
    touch X3607SB8.star.failed
  fi

  STAR --version | awk '{print "STAR,v"$1}' > X3607SB8.count.star.version
  rm -rf X3607SB8__STARgenome X3607SB8__STARpass1 X3607SB8_Log.progress.out
fi

if [[ $status -eq 0 ]]; then
  echo bamStat=`date` 
  python3 /data/cqs/softwares/ngsperl/lib/Alignment/bamStat.py -i X3607SB8_Aligned.out.bam -o X3607SB8.bamstat
fi

if [[ $status -eq 0 ]]; then
  echo bamSort=`date` 
  samtools sort -m 5G -T X3607SB8 -t 8 -o X3607SB8_Aligned.sortedByCoord.out.bam X3607SB8_Aligned.out.bam && touch X3607SB8_Aligned.sortedByCoord.out.bam.succeed
  if [[ ! -e X3607SB8_Aligned.sortedByCoord.out.bam.succeed ]]; then
    rm X3607SB8_Aligned.sortedByCoord.out.bam
  else
    samtools index X3607SB8_Aligned.sortedByCoord.out.bam
    samtools idxstats X3607SB8_Aligned.sortedByCoord.out.bam > X3607SB8_Aligned.sortedByCoord.out.bam.chromosome.count
  fi
fi  


if [[ $status -eq 0 && -s X3607SB8_Aligned.sortedByCoord.out.bam ]]; then
  
  
  if [ ! -s X3607SB8.splicing.bed ]; then
    awk {'if($4=="2") print ""$1"\t"$2-$9-1"\t"$3+$9"\tJUNC000"NR"\t"$7+$8"\t-\t"$2-$9-1"\t"$3+$9"\t255,0,0\t2\t"$9","$9"\t","0,"$3-$2+$9+1; else if($4=="1") print ""$1"\t"$2-$9-1"\t"$3+$9"\tJUNC000"NR"\t"$7+$8"\t+\t"$2-$9-1"\t"$3+$9"\t0,0,255\t2\t"$9","$9"\t","0,"$3-$2+$9+1'} X3607SB8_SJ.out.tab > X3607SB8.splicing.bed
  fi
fi

if [[ $status -eq 0 ]]; then
  echo performing featureCounts ...
  featureCounts -g gene_id -t exon -p --countReadPairs -T 8 -a /data/cqs/references/gencode/GRCh38.p13/gencode.v38.annotation.gtf -o X3607SB8.count X3607SB8_Aligned.out.bam
  status=$?
  if [[ $status -eq 0 ]]; then
    touch X3607SB8.featureCount.succeed
  else
    touch X3607SB8.featureCount.failed
    rm X3607SB8.count
  fi

  featureCounts -v 2>&1 | grep featureCounts | cut -d ' ' -f2 | awk '{print "featureCounts,"$1}' > X3607SB8.count.featureCounts.version
fi 

if [[ -s X3607SB8.count && -s X3607SB8.bamstat ]]; then
  rm X3607SB8_Aligned.out.bam 
fi

rm X3607SB8_clipped.1.fastq.gz X3607SB8_clipped.2.fastq.gz

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
