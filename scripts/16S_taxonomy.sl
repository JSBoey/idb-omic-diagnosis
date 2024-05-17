#!/bin/bash -e
#SBATCH --job-name=16S_taxonomy
#SBATCH --account=uoo03475
#SBATCH --time=3:00:00
#SBATCH --mem=80G
#SBATCH --cpus-per-task=12
#SBATCH --partition=milan
#SBATCH --array=0-1
#SBATCH --output=slurm_out/%x.%j.%A_%a.out
#SBATCH --error=slurm_err/%x.%j.%A_%a.err

module purge
module load QIIME2/2023.5

a=(Greengenes2_2022.10 SILVA_NR99_138.1)
cls=db/${a[$SLURM_ARRAY_TASK_ID]}_classifier.qza
db=${a[$SLURM_ARRAY_TASK_ID]%%_*}

reads=16S_rep_seq.qza

qiime feature-classifier classify-sklearn \
    --verbose \
    --p-n-jobs $SLURM_CPUS_PER_TASK \
    --i-reads 2.asv/${reads} \
    --i-classifier ${cls} \
    --o-classification 4.taxonomy/${reads%%_*}.${db}.qza