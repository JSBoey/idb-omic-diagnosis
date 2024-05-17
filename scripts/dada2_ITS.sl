#!/bin/bash -e
#SBATCH --job-name=dada2_paired_ITS
#SBATCH --account=uoo03475
#SBATCH --time=2:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=24
#SBATCH --partition=milan
#SBATCH --output=slurm_out/%x.%j.out
#SBATCH --error=slurm_err/%x.%j.err

module purge
module load QIIME2/2023.5

cd $PWD/2.asv
echo "Denoise ITS"

qiime dada2 denoise-paired \
    --verbose \
    --i-demultiplexed-seqs ITS_demux.qza \
    --p-trunc-len-f 212 \
    --p-trunc-len-r 0 \
    --p-n-threads $SLURM_CPUS_PER_TASK \
    --o-table ITS_feature_table.qza \
    --o-representative-sequences ITS_rep_seq.qza \
    --o-denoising-stats ITS_denoise_stats.qza

qiime feature-table tabulate-seqs \
    --i-data ITS_rep_seq.qza \
    --o-visualization ITS_rep_seq.qzv

qiime metadata tabulate \
    --m-input-file ITS_denoise_stats.qza \
    --o-visualization ITS_denoise_stats.qzv

