#!/bin/bash -e
#SBATCH --job-name=curate_silva
#SBATCH --account=uoo03475
#SBATCH --time=1:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=32
#SBATCH --partition=milan
#SBATCH --output=slurm_out/%x.%j.out
#SBATCH --error=slurm_err/%x.%j.err

module purge
module load QIIME2/2023.5

silva=db/SILVA_NR99_138.1/SILVA_NR99_138.1_refseq.qza
taxa=db/SILVA_NR99_138.1/SILVA_NR99_138.1_taxonomy.qza

# Remove low-quality sequences
cull=${silva/_refseq/_cull_refseq}

qiime rescript cull-seqs \
    --p-n-jobs $SLURM_CPUS_PER_TASK \
    --i-sequences ${silva} \
    --o-clean-sequences ${cull}

# Remove short sequences
filt=${cull/_cull/_filt}
discard=${cull/_cull/discard}

qiime rescript filter-seqs-length-by-taxon \
    --i-sequences ${cull} \
    --i-taxonomy ${taxa} \
    --p-labels Archaea Bacteria Eukaryota \
    --p-min-lens 900 1200 1400 \
    --o-filtered-seqs ${filt} \
    --o-discarded-seqs ${discard}

# Dereplicate
derep_seq=${filt/_filt/_derep}
derep_tax=${taxa/_taxonomy/_derep_taxonomy}

qiime rescript dereplicate \
    --verbose \
    --p-threads $SLURM_CPUS_PER_TASK \
    --i-sequences ${filt}  \
    --i-taxa ${taxa} \
    --p-mode 'uniq' \
    --o-dereplicated-sequences ${derep_seq} \
    --o-dereplicated-taxa ${derep_tax}

# In-silico PCR
segment=${derep_seq/_derep/_segment}
F16S='CCTACGGRRBGCASCAGKVRVGAAT'
R16S='GGACTACNVGGGTWTCTAATCC'

qiime feature-classifier extract-reads \
    --i-sequences ${derep_seq} \
    --p-f-primer ${F16S} \
    --p-r-primer ${R16S} \
    --p-n-jobs $SLURM_CPUS_PER_TASK \
    --p-read-orientation 'forward' \
    --o-reads ${segment}

# Dereplicate segments
derep_segment_seq=${segment/_segment/_derep_segment}
derep_segment_tax=${derep_tax/_derep/_derep_segment}

qiime rescript dereplicate \
    --verbose \
    --p-threads $SLURM_CPUS_PER_TASK \
    --i-sequences ${segment}  \
    --i-taxa ${derep_tax} \
    --p-mode 'uniq' \
    --o-dereplicated-sequences ${derep_segment_seq} \
    --o-dereplicated-taxa ${derep_segment_tax}
