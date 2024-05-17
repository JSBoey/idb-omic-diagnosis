#!/bin/bash -e
#SBATCH --job-name=phylogeny
#SBATCH --account=uoo03475
#SBATCH --time=2:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=32
#SBATCH --partition=milan
#SBATCH --array=0-1
#SBATCH --output=slurm_out/%x.%j.%A_%a.out
#SBATCH --error=slurm_err/%x.%j.%A_%a.err

module purge
module load QIIME2/2023.5

a=(16S_rep_seq.qza \
   ITS_ITSxpress_rep_seq.se.qza)
b=${a[$SLURM_ARRAY_TASK_ID]}
gene=3.phylogeny/${b%%_*}

qiime phylogeny align-to-tree-mafft-fasttree \
    --verbose \
    --p-n-threads $SLURM_CPUS_PER_TASK \
    --i-sequences 2.asv/${b} \
    --o-alignment ${gene}.alignment.qza \
    --o-masked-alignment ${gene}.masked_alignment.qza \
    --o-tree ${gene}.unrooted_tree.qza \
    --o-rooted-tree ${gene}.rooted_tree.qza

