#!/bin/bash -e
#SBATCH --job-name=ITSxpress
#SBATCH --account=uoo03475
#SBATCH --time=30:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=24
#SBATCH --partition=milan
#SBATCH --array=0-25
#SBATCH --output=slurm_out/%x.%j.%A_%a.out
#SBATCH --error=slurm_err/%x.%j.%A_%a.err

module purge
module load Mamba/23.1.0-1

a=(0.fq/1.adapter_trimmed/*ITS_at.1.fq.gz)
r1=${a[$SLURM_ARRAY_TASK_ID]}
r2=${r1/.1./.2.}
b=$(basename ${r1} _at.1.fq.gz)

o1=0.fq/2.ITSxpress_trimmed/${b}_it.1.fq.gz
o2=${o1/.1./.2.}

id=1.0
taxa=Fungi
region=ITS2
log=1.qc/2.ITSxpress_trimmed/${b}.itsxpress.log

mkdir -p tmp/

source ~/.bashrc
conda activate conda_env/itsxpress

echo "ITSxpress trimming sample ${b}"

itsxpress \
    --fastq ${r1} \
    --fastq2 ${r2} \
    --outfile ${o1} \
    --outfile2 ${o2} \
    --region ${region} \
    --taxa ${taxa} \
    --cluster_id ${id} \
    --log ${log} \
    --tempdir tmp/ \
    --threads $SLURM_CPUS_PER_TASK

mamba deactivate
