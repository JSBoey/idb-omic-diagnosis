#!/bin/bash -e
#SBATCH --job-name=train_classifier
#SBATCH --account=uoo03475
#SBATCH --time=3:00:00
#SBATCH --mem=80G
#SBATCH --cpus-per-task=2
#SBATCH --partition=milan
#SBATCH --array=0-2
#SBATCH --output=slurm_out/%x.%j.%A_%a.out
#SBATCH --error=slurm_err/%x.%j.%A_%a.err

module purge
module load QIIME2/2023.5

a=(SILVA_NR99_138.1_derep_segment \
   Greengenes2_2022.10_derep_segment \
   UNITE_dev99_10.0)
b=${a[$SLURM_ARRAY_TASK_ID]}
db=$(echo ${b%%_derep_segment} | sed 's/dev99_//g')

seq=db/${db}/${b}_refseq.qza
tax=db/${db}/${b}_taxonomy.qza
cls=db/${db}_classifier.qza

qiime feature-classifier fit-classifier-naive-bayes \
    --p-verbose \
    --i-reference-reads ${seq} \
    --i-reference-taxonomy ${tax} \
    --o-classifier ${cls}
