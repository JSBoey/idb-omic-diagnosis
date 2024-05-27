#!/bin/bash -e
#SBATCH --job-name=dada2_single
#SBATCH --account=uoo03475
#SBATCH --time=1:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=24
#SBATCH --partition=milan
#SBATCH --array=0-1
#SBATCH --output=slurm_out/%x.%j.%A_%a.out
#SBATCH --error=slurm_err/%x.%j.%A_%a.err

module purge
module load QIIME2/2023.5

cd $PWD/2.asv
echo "Denoise single-end forward reads"

a=(16S ITS_ITSxpress)
b=${a[$SLURM_ARRAY_TASK_ID]}

qiime dada2 denoise-single \
    --verbose \
    --i-demultiplexed-seqs ${b}_demux.qza \
    --p-trunc-len 0 \
    --p-n-threads $SLURM_CPUS_PER_TASK \
    --o-table ${b}_feature_table.se.qza \
    --o-representative-sequences ${b}_rep_seq.se.qza \
    --o-denoising-stats ${b}_denoise_stats.se.qza

qiime feature-table tabulate-seqs \
    --i-data ${b}_rep_seq.se.qza \
    --o-visualization ${b}_rep_seq.se.qzv

qiime metadata tabulate \
    --m-input-file ${b}_denoise_stats.se.qza \
    --o-visualization ${b}_denoise_stats.se.qzv
