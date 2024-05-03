# Directories

Project directory: `uoo03475`
- `Data from Azenta/R1_16S/Rawdata` Raw data for 16S amplicons
- `Data from Azenta/R2_ITS/Rawdata` Raw data for 18S amplicons
- `1.sequence_processing` Processed sequences
- `db` SILVA and UNITE databases
- `nobackup_space` -> `nobackup` mirror of project

No backup directory

# 0. Raw sequence checks

Checking file integrity

```bash
cd /nesi/project/uoo03475/Data\ from\ Azenta/

for i in R*/Rawdata; do
  cd ${i}
  md5sum -c md5.md5
  cd ../../
done
```

File integrity OK.

Leverage low-latency`nobackup` space for all further analyses.

```bash
cd /nesi/nobackup/uoo03475

for i in R{1_16S,2_ITS}; do
  seqtype=$(echo ${i} | sed 's/R._//')
  mkdir -p 1.sequence_processing/${seqtype}/0.rawdata
  cp --verbose \
    project_space/Data\ from\ Azenta/${i}/Rawdata/* \
    1.sequence_processing/${seqtype}/0.rawdata/
done

for i in 1.sequence_processing/*/0.rawdata; do
  cd ${i}
  md5sum -c md5.md5
  cd ../../../
done
```

Successful copy.

Check sequence quality

`scripts/0.seq_qc.sh`

```bash
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
  multiqc --title ${i} -comment "Raw FASTQ" --outdir ${out} ${in}
done

```

