#!/bin/bash
#SBATCH --mail-user=quanhu.sheng.1@vumc.org
#SBATCH --mail-type=FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mem=44G
#SBATCH --constraint=haswell
#SBATCH -o /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/log/human_exomeseq_ITER_034_gcc.log

if [ -n "${SLURM_JOB_ID}" ] ; then
  smemwatch -k 99 -d 50 $$ &
fi






set -o pipefail

cd '/scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_034'


if [[ !(1 -eq $1) ]]; then
  if [[ ( -s /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_034/gcc-calls/SAMPLE_0/sample_name.txt ) || ( -d /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_034/gcc-calls/SAMPLE_0/sample_name.txt ) ]]; then
    echo job has already been done. if you want to do again, delete /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_034/gcc-calls/SAMPLE_0/sample_name.txt and submit job again.
    exit 0
  fi
fi

echo GATK4::GermlineCNVCallerScatter_start=`date`
echo working in /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/result/human_exomeseq_ITER_034 ...
 

export R_LIBS=
export PYTHONPATH=
export JAVA_HOME=
 


singularity exec -c -B /gpfs51,/dors,/gpfs23,/scratch,/gpfs52,/data,/home,/tmp -H `pwd` -e /data/cqs/softwares/singularity/cqs-gatk4.simg  bash /scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_05_GermlineCNVCaller_2_scatterCall/pbs/human_exomeseq_ITER_034_gcc.pbs.sh 

echo GATK4::GermlineCNVCallerScatter_end=`date`

exit 0

