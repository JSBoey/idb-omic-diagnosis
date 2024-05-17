#!/bin/bash -e
#SBATCH --job-name=ITS_taxonomy
#SBATCH --account=uoo03475
#SBATCH --time=3:00:00
#SBATCH --mem=80G
#SBATCH --cpus-per-task=12
#SBATCH --partition=milan
#SBATCH --output=slurm_out/%x.%j.out
#SBATCH --error=slurm_err/%x.%j.err

module purge
module load QIIME2/2023.5

cls=db/UNITE_10.0_classifier.qza
db=$(basename ${cls%%_*})

reads=ITS_ITSxpress_rep_seq.se.qza

qiime feature-classifier classify-sklearn \
    --verbose \
    --p-reads-per-batch 10000 \
    --p-n-jobs $SLURM_CPUS_PER_TASK \
    --i-reads 2.asv/${reads} \
    --i-classifier ${cls} \
    --o-classification 4.taxonomy/${reads%%_*}.${db}.qza