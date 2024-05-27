#!/bin/bash -e

module purge
module load cutadapt/4.4-gimkl-2022a-Python-3.11.3

# Adapter trimming
# Also trim poly-G tails
echo "Trimming adapters"

A1='AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
A2='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'

for r1 in 0.fq/0.raw/*_L1.1.fq.gz; do
    r2=${r1/\.1\./\.2\.}
    bn=$(basename ${r1} _L1.1.fq.gz)
    od=0.fq/1.adapter_trimmed

    cutadapt \
        -j 8 \
        --nextseq-trim=20 \
	--poly-a \
        --json 1.qc/1.adapter_trimmed/${bn}_at.cutadapt.json \
        -a ${A1} \
        -A ${A2} \
        -o ${od}/${bn}_at.1.fq.gz \
        -p ${od}/${bn}_at.2.fq.gz \
        ${r1} ${r2} \
        > 1.qc/1.adapter_trimmed/${bn}_at.cutadapt.log
done


# Primer trimming (16S)
echo "Trimming primers: 16S"

F16S='^CCTACGGRRBGCASCAGKVRVGAAT'
R16S='^GGACTACNVGGGTWTCTAATCC'

for r1 in 0.fq/1.adapter_trimmed/*16S_at.1.fq.gz; do
    r2=${r1/\.1\./\.2\.}
    bn=$(basename ${r1} _at.1.fq.gz)
    od=0.fq/2.primer_trimmed

    cutadapt \
        -j 8 --minimum-length 1 \
        --json 1.qc/2.primer_trimmed/${bn}_pt.cutadapt.json \
        -g ${F16S} \
        -G ${R16S} \
        -o ${od}/${bn}_pt.1.fq.gz \
        -p ${od}/${bn}_pt.2.fq.gz \
        ${r1} ${r2} \
        > 1.qc/2.primer_trimmed/${bn}_pt.cutadapt.log
done


# Primer trimming (ITS)
echo "Trimming primers: ITS"

FITS='^GTGAATCATCGARTC'
RITS='^TCCTCCGCTTATTGAT'

for r1 in 0.fq/1.adapter_trimmed/*ITS_at.1.fq.gz; do
    r2=${r1/\.1\./\.2\.}
    bn=$(basename ${r1} _at.1.fq.gz)
    od=0.fq/2.primer_trimmed

    cutadapt \
        -j 8 --minimum-length 1 \
        --json 1.qc/2.primer_trimmed/${bn}_pt.cutadapt.json \
        -g ${FITS} \
        -G ${RITS} \
        -o ${od}/${bn}_pt.1.fq.gz \
        -p ${od}/${bn}_pt.2.fq.gz \
        ${r1} ${r2} \
        > 1.qc/2.primer_trimmed/${bn}_pt.cutadapt.log
done

