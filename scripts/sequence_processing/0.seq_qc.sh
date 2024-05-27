#!/bin/bash -e

module purge
module load \
  FastQC/0.12.1 \
  MultiQC/1.13-gimkl-2022a-Python-3.10.5

for i in 16S ITS; do
  mkdir -p 1.sequence_processing/${i}/0.qc
  in=1.sequence_processing/${i}/0.rawdata
  out=1.sequence_processing/${i}/0.qc
  fastqc --threads 8 --outdir ${out} ${in}/*.fq.gz
  multiqc --title ${i} --comment "Raw FASTQ" --outdir ${out} ${out}
done

module purge
