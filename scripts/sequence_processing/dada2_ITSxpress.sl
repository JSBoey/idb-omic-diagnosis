#!/bin/bash -e
#SBATCH --job-name=dada2_ITSxpress
#SBATCH --account=uoo03475
#SBATCH --time=1:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=24
#SBATCH --partition=milan
#SBATCH --output=slurm_out/%x.%j.out
#SBATCH --error=slurm_err/%x.%j.err

module purge
module load QIIME2/2023.5

cd $PWD/2.asv
echo "Denoise ITSxpress trimmed ITS reads"

b=ITS_ITSxpress

qiime dada2 denoise-paired \
    --verbose \
    --i-demultiplexed-seqs ${b}_demux.qza \
    --p-trunc-len-f 0 \
    --p-trunc-len-r 0 \
    --p-n-threads $SLURM_CPUS_PER_TASK \
    --o-table ${b}_feature_table.qza \
    --o-representative-sequences ${b}_rep_seq.qza \
    --o-denoising-stats ${b}_denoise_stats.qza

qiime feature-table tabulate-seqs \
    --i-data ${b}_rep_seq.qza \
    --o-visualization ${b}_rep_seq.qzv

qiime metadata tabulate \
    --m-input-file ${b}_denoise_stats.qza \
    --o-visualization ${b}_denoise_stats.qzv

