#!/bin/bash -e
#SBATCH --job-name=curate_gg2
#SBATCH --account=uoo03475
#SBATCH --time=2:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=32
#SBATCH --partition=milan
#SBATCH --output=slurm_out/%x.%j.out
#SBATCH --error=slurm_err/%x.%j.err

module purge
module load QIIME2/2023.5

gg=db/Greengenes2_2022.10/Greengenes2_2022.10_refseq.qza
tax=db/Greengenes2_2022.10/Greengenes2_2022.10_taxonomy.qza

F16S='CCTACGGRRBGCASCAGKVRVGAAT'
R16S='GGACTACNVGGGTWTCTAATCC'

# In-silico PCR
segment=${gg/_refseq/_segment_refseq}

qiime feature-classifier extract-reads \
    --i-sequences ${gg} \
    --p-f-primer ${F16S} \
    --p-r-primer ${R16S} \
    --p-n-jobs $SLURM_CPUS_PER_TASK \
    --p-read-orientation 'forward' \
    --o-reads ${segment}

# Dereplicate segments
derep_segment_seq=${segment/_segment/_derep_segment}
derep_segment_tax=${tax/_taxonomy/_derep_segment_taxonomy}

qiime rescript dereplicate \
    --verbose \
    --p-threads $SLURM_CPUS_PER_TASK \
    --i-sequences ${segment}  \
    --i-taxa ${tax} \
    --p-mode 'uniq' \
    --o-dereplicated-sequences ${derep_segment_seq} \
    --o-dereplicated-taxa ${derep_segment_tax}

