#!/bin/bash -e
#SBATCH --job-name=dada2_paired_16S
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
echo "Denoise 16S"

qiime dada2 denoise-paired \
    --verbose \
    --i-demultiplexed-seqs 16S_demux.qza \
    --p-trunc-len-f 0 \
    --p-trunc-len-r 0 \
    --p-n-threads 8 \
    --o-table 16S_feature_table.qza \
    --o-representative-sequences 16S_rep_seq.qza \
    --o-denoising-stats 16S_denoise_stats.qza

qiime feature-table tabulate-seqs \
    --i-data 16S_rep_seq.qza \
    --o-visualization 16S_rep_seq.qzv

qiime metadata tabulate \
    --m-input-file 16S_denoise_stats.qza \
    --o-visualization 16S_denoise_stats.qzv
